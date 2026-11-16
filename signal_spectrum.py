# -*- coding: utf-8 -*-
"""
Python translation of R/eseis::signal_spectrum()

Calculate the spectrum of a time series using either:
- modified periodogram ('periodogram'), or
- autoregressive spectral estimation ('autoregressive')

Optional features:
- running-mean smoothing (n samples)
- spline interpolation to a fixed resolution (res), with linear or log frequency spacing

Returns: dict {'frequency': np.ndarray, 'power': np.ndarray}
Or, for 'eseis-like' inputs (dict with keys 'signal' and 'meta': {'dt': ...}), returns
the same dict updated with a 'spectrum' key containing the frequency/power table
and a minimal 'history' entry.

Dependencies:
- numpy, scipy (signal, interpolate)
- Optional: spectrum (for Burg AR PSD), statsmodels (fallback AR)

Notes:
- Frequencies returned are Hz, from 0 to Nyquist (fs/2).
- PSD ('power') is linear (x^2/Hz). Convert to dB for plotting via 10*log10(power).
- dt is the *sampling period* (seconds); fs=1/dt is sampling frequency (Hz).

R sources: eseis::signal_spectrum (periodogram, AR; smoothing; interpolation).
  See: https://rdrr.io/github/coffeemuggler/eseis/src/R/signal_spectrum.R
"""
import warnings
from typing import Any, Dict, List, Optional, Tuple, Union

import numpy as np
from scipy.signal import periodogram, windows
from scipy.interpolate import interp1d

# Optional AR dependencies
try:
    # Burg AR (closest match to R's spec.ar common usage)
    from spectrum import burg
    _HAVE_SPECTRUM = True
except Exception:
    _HAVE_SPECTRUM = False

try:
    # Yule-Walker fallback
    import statsmodels.api as sm
    _HAVE_STATSMODELS = True
except Exception:
    _HAVE_STATSMODELS = False


def _running_mean(x: np.ndarray, k: int) -> np.ndarray:
    """Running mean identical in spirit to caTools::runmean."""
    if k is None or k <= 1:
        return x
    # 'same' keeps vector length aligned with R behavior (border effects smoothed)
    kernel = np.ones(int(k), dtype=float) / float(k)
    return np.convolve(x, kernel, mode='same')


def _ar_psd_burg(x: np.ndarray, fs: float, order: Optional[int] = None) -> Tuple[np.ndarray, np.ndarray]:
    """
    Burg AR spectrum. If order not provided, use a modest default (e.g., 16)
    which is common for smooth ambient noise spectra.

    Returns f (Hz) and PSD (linear).
    """
    if order is None:
        order = 16
    if not _HAVE_SPECTRUM:
        raise ImportError("Burg AR requires the 'spectrum' package. Install `pip install spectrum`.")

    # Burg returns AR coefficients; we evaluate PSD on a frequency grid
    # Use a dense frequency grid up to Nyquist
    nfft = 4096
    f = np.linspace(0.0, fs / 2.0, num=nfft//2 + 1)

    # Estimate AR model
    ar = burg.Burg(x, order=order)
    ar()

    # Compute AR PSD: sigma2 / |1 - sum a_k e^{-j2πf k/fs}|^2
    # spectrum library has a method for PSD over normalized freq; we adapt to Hz.
    # Use ar.spectrum() which returns normalized freq [0, 0.5]; scale to Hz.
    psd_norm = ar.spectrum()  # default length nfft; normalized frequencies [0, 0.5]
    # Build matching Hz axis for psd_norm:
    f_norm = np.linspace(0.0, 0.5, len(psd_norm))
    f_hz = f_norm * fs
    # Clip to [0, fs/2]
    m = (f_hz >= 0) & (f_hz <= fs / 2.0)
    return f_hz[m], np.asarray(psd_norm)[m]


def _ar_psd_yule_walker(x: np.ndarray, fs: float, order: Optional[int] = None) -> Tuple[np.ndarray, np.ndarray]:
    """
    Yule-Walker AR PSD fallback using statsmodels.
    """
    if not _HAVE_STATSMODELS:
        raise ImportError("Yule-Walker fallback requires 'statsmodels'. Install `pip install statsmodels`.")
    if order is None:
        order = 16

    # Estimate AR params via Yule-Walker
    rho, sigma = sm.tsa.acf(x, fft=True, nlags=order, adjusted=False), None
    # Fit AR via Yule-Walker (statsmodels has yule_walker helper)
    ar_params, sigma = sm.tsa.stattools.yule_walker(x, order=order)
    # Build PSD on Hz grid
    nfft = 4096
    f = np.linspace(0.0, fs / 2.0, num=nfft//2 + 1)
    w = 2.0 * np.pi * f / fs
    # AR transfer function magnitude squared: sigma^2 / |1 - sum a_k e^{-j w k}|^2
    denom = np.ones_like(w, dtype=complex)
    for k, a_k in enumerate(ar_params, start=1):
        denom -= a_k * np.exp(-1j * w * k)
    psd = (sigma**2) / np.abs(denom)**2
    return f, psd.real


def signal_spectrum(
    data: Union[np.ndarray, Dict[str, Any], List[Any]],
    dt: Optional[float] = None,
    method: str = "periodogram",
    n: Optional[int] = None,
    res: Optional[int] = None,
    log: bool = False,
    **kwargs
) -> Union[Dict[str, np.ndarray], Dict[str, Any], List[Any]]:
    """
    Calculate the spectrum of a time series (R/eseis::signal_spectrum port).

    Parameters
    ----------
    data : np.ndarray | dict | list
        Numeric vector (1D numpy array), an 'eseis-like' dict {'signal': array, 'meta': {'dt': ...}},
        or a list of such elements.
    dt : float, optional
        Sampling period [s]. If omitted and 'data' is not eseis-like, defaults to 1/200 s, with a warning.
        If 'data' is eseis-like, dt is read from data['meta']['dt'].
    method : {'periodogram', 'autoregressive'}
        Spectrum estimator. Default 'periodogram' (modified periodogram).
    n : int, optional
        Running-mean smoothing window (in samples). If provided, applies a running mean to the PSD.
    res : int, optional
        Target resolution (#points) for the spectrum; if provided, cubic spline interpolation is performed.
    log : bool
        If True, interpolation uses log-spaced frequency values; otherwise linear spacing.
    **kwargs : dict
        Additional keyword args passed to the estimator. For periodogram, typical kwargs include 'window', 'scaling'.

    Returns
    -------
    dict or eseis-like dict or list
        If input is an array: {'frequency': f, 'power': Pxx}
        If input is eseis-like dict: returns that dict updated with 'spectrum'
            (DataFrame-like dict) and minimal 'history'.
        If input is list: list of outputs for each element.
    """
    # Handle list input
    if isinstance(data, list):
        return [signal_spectrum(d, dt=dt, method=method, n=n, res=res, log=log, **kwargs) for d in data]

    # Eseis-like object?
    eseis_class = False
    if isinstance(data, dict) and ('signal' in data) and ('meta' in data) and ('dt' in data['meta']):
        eseis_class = True
        eseis_data = data
        x = np.asarray(eseis_data['signal'], dtype=float)
        dt = eseis_data['meta']['dt']
    else:
        # Raw numeric vector
        x = np.asarray(data, dtype=float)
        if dt is None:
            warnings.warn("Sampling frequency missing! Set to 1/200", RuntimeWarning)
            dt = 1.0 / 200.0  # match R's default
    fs = 1.0 / float(dt)  # sampling frequency in Hz

    # Compute spectrum
    method = (method or "periodogram").lower()
    if method == "periodogram":
        # SciPy periodogram (linear PSD) – closest to R spec.pgram default behavior
        # You can pass window/scaling in kwargs; e.g., window='hann', scaling='density'

        if 'window' not in kwargs:
            kwargs.setdefault('scaling', 'density')
            kwargs.setdefault('detrend', 'constant')
            kwargs.setdefault('nfft', x.size)
            f, Pxx = periodogram(x,fs=fs,**kwargs)

        if res is not None: # replaces cubic interpolation with linear interpolation
            if log:
                fmin_pos = np.min(f[f>0.0]) if np.any(f>0.0) else f[1]
                f_int = np.exp(np.linspace(np.log(fmin_pos), np.log(np.max(f)), num=int(res)))
            else:
                f_int = np.linspace(np.min(f), np.max(f), num=int(res))
            
            Pxx = np.interp(f_int, f, Pxx)
            f = f_int
        
        
        # N = len(x)
        # ntaper = int(0.1 * N)
        # windowz = np.hanning(ntaper * 2)
        # win_full = np.ones(N)
        # win_full[:ntaper] = windowz[:ntaper]
        # win_full[-ntaper:] = windowz[-ntaper:]
        # x = x - np.mean(x)
        # f, Pxx = periodogram(x, fs=fs, window=win_full, scaling="density")

        #f, Pxx = periodogram(x, fs=fs, **kwargs)
    elif method == "autoregressive":
        # Try Burg AR first; fallback to Yule-Walker
        try:
            f, Pxx = _ar_psd_burg(x, fs, order=kwargs.get('order', None))
        except Exception:
            f, Pxx = _ar_psd_yule_walker(x, fs, order=kwargs.get('order', None))
    else:
        raise ValueError("Keyword for method not supported! Use 'periodogram' or 'autoregressive'.")

    # Optional running-mean smoothing (n samples)
    if n is not None:
        Pxx = _running_mean(Pxx, int(n))

    # Optional interpolation to 'res' points (linear or log frequency spacing)
    if res is not None:
        # Build interpolation target frequency vector
        if log:
            # Use log-spaced frequencies between the smallest positive freq and max
            # (avoid f=0 for log spacing)
            fmin_pos = np.min(f[f > 0.0]) if np.any(f > 0.0) else f[1]
            f_int = np.exp(np.linspace(np.log(fmin_pos), np.log(np.max(f)), num=int(res)))
        else:
            f_int = np.linspace(np.min(f), np.max(f), num=int(res))

        # Cubic spline interpolation akin to stats::spline
        # Interpolate both frequency 'mapping' and power like the R code does
        # (Frequency mapping step is redundant in Python but kept for parity.)
        s_power = interp1d(f, Pxx, kind='cubic', bounds_error=False, fill_value="extrapolate")
        s_freq  = interp1d(f, f,   kind='cubic', bounds_error=False, fill_value="extrapolate")

        Pxx = s_power(f_int)
        f   = s_freq(f_int)

    # Compose output
    out_table = {'frequency': f, 'power': Pxx}

    if eseis_class:
        # Update 'eseis-like' dict with spectrum and minimal history
        eseis_data['spectrum'] = out_table
        # Minimal history entry (mirrors R semantics at a high level)
        hist_entry = {
            'time': 'now',  # left as string; add datetime if desired
            'call': 'signal_spectrum()',
            'arguments': {'dt': dt, 'n': n, 'res': res, 'method': method}
        }
        # Append to history
        try:
            eseis_data.setdefault('history', []).append(hist_entry)
            eseis_data['meta']['type'] = 'spectrum'
        except Exception:
            pass
        return eseis_data


    out_table = {'frequency': f, 'power': Pxx}
    
    if np.any(m):
        band = Pxx[m]
        print("Freq band:", f[m][0], f[m][-1])
        print("dB range in band:",
            10*np.log10(band).min(),
            10*np.log10(band).max())
    else:
        print("No values in 1–150 Hz band!")

    return out_table # returns linear so need to plot as: psd_db = 10.0*np.log10(Pxx)

out_table = {'frequency': f, 'power': Pxx}

# === debug: confirm AFTER smoothing + interpolation ===
m = (f >= 1.0) & (f <= 150.0)
if np.any(m):
    band = Pxx[m]
    # Use a sane floor to avoid absurdly low dB from near-zero PSD
    band = np.clip(band, 1e-300, None)
    print("Freq band:", f[m][0], f[m][-1])
    print("dB range in band:", 10*np.log10(band).min(), 10*np.log10(band).max())
else:
    print("No values in 1–150 Hz band!")
# =================================================================

if eseis_class:
    eseis_data['spectrum'] = out_table
    # ... (history updates)
    return eseis_data

return out_table




# optional:
#m = (f >= 1.0) & (f <= 150.0)
#plt.plot(f[m], 10*np.log10(Pxx[m]))





# #### useage examples
# import numpy as np
# from signal_spectrum import signal_spectrum

# # Example 1: numeric vector, default periodogram
# dt = 1.0 / 300.0  # 300 Hz sampling
# t = np.arange(0, 180.0, dt)                      # 180 s window
# x = 0.5*np.sin(2*np.pi*20*t) + 0.2*np.sin(2*np.pi*70*t) + 0.05*np.random.randn(t.size)

# spec = signal_spectrum(x, dt=dt, method="periodogram", n=51, res=500, log=False)
# f, Pxx = spec['frequency'], spec['power']
# # plot in dB like your models:
# # plt.semilogx(f, 10*np.log10(Pxx))

# # Example 2: autoregressive spectrum (Burg), order 16
# spec_ar = signal_spectrum(x, dt=dt, method="autoregressive", res=500, log=True, order=16)
# f_ar, Pxx_ar = spec_ar['frequency'], spec_ar['power']

# # Example 3: eseis-like object (dict)
# eseis_obj = {'signal': x, 'meta': {'dt': dt, 'type': 'signal'}, 'history': []}
# spec_obj = signal_spectrum(eseis_obj, method="periodogram", n=31)
# # spec_obj['spectrum'] -> {'frequency': f, 'power': Pxx}






# ### then to integrate, after reading data and removing instrument response:
# from obspy import read
# st = read("ZE.2416..GPZ_20240726120000.mseed")
# tr = st[0].copy()
# # tr.remove_response(inventory=inv, pre_filt=(0.5,0.8,150,180), output="VEL")  # m/s
# dt = tr.stats.delta

# spec = signal_spectrum(x, dt=dt, method="periodogram", n=51, res=800, log=False)
# f, Pxx = spec['frequency'], spec['power']





# ### then put this in the pyseis function
# # Combined model power (linear), then convert to dB
# p_w_db = 10*np.log10(p_w["power"])
# p_s_db = 10*np.log10(p_s["power"])
# p_sum_db = 10*np.log10(p_w["power"] + p_s["power"])
# # Observed (in dB) clipped to [1,150] Hz:
# m = (f >= 1.0) & (f <= 150.0)
# plt.plot(f[m], 10*np.log10(Pxx[m]), color='orange', label='Observed (periodogram)')
# plt.plot(p_w["frequency"], p_w_db, label='Turbulence')
# plt.plot(p_s["frequency"], p_s_db, label='Bedload')
# plt.plot(p_w["frequency"], p_sum_db, label='Combined')
# plt.ylim(-160, -110); plt.xlim(1, 150)
# plt.legend(); plt.grid(True, which='both', ls='--', alpha=0.4)
# plt.show()