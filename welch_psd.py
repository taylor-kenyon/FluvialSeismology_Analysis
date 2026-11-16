## welch psd attempt

import numpy as np
from scipy.signal import welch

def welch_psd(S, inv, fmin=1.0, fmax=150.0,
                      nperseg=60 * 300,   # 300 Hz * 60 s = 18,000 samples (~60 s segments/subwindows)
                      noverlap=None,      # default: 50% if None
                      window='hann',
                      scaling='density'):
    """
    Compute Welch PSD (density) for the selected component over [t1, t2],
    remove instrument response, and return frequency/PSD in dB clipped to [fmin,fmax].
    """
    # 1) pick main component to compare (GPZ)
    tr = None
    for _tr in S:
        if _tr.stats.channel == channel:  # use the same 'channel' parsed earlier
            tr = _tr.copy()
            break
    if tr is None:
        # fallback to the first trace
        tr = S[0].copy()

    # 2) remove instrument response to physical units (velocity, m/s)
    #    use pre-filter to avoid amplifying noise near poles/zeros
    pre_filt = (0.5, 0.8, 150.0, 180.0)  # adjust for your band/sensor
    tr.remove_response(inventory=inv, pre_filt=pre_filt, output="VEL")  # m/s

    # 3) welch PSD
    fs = tr.stats.sampling_rate
    if noverlap is None:
        noverlap = nperseg // 2
    f, Pxx = welch(tr.data, fs=fs, window=window, nperseg=nperseg,
                   noverlap=noverlap, detrend='constant',
                   return_onesided=True, scaling=scaling)

    # 4) clip to [fmin, fmax] and convert to dB
    m = (f >= fmin) & (f <= fmax)
    f = f[m]
    Pxx = Pxx[m]
    psd_db = 10.0*np.log10(Pxx)  # dB

    return f, psd_db



## after calling the above, run:
# after existing npz export (or instead of), get welch PSD:
freq_welch, psd_welch_db = welch_psd(S.slice(t1, t2), inv,
                                             fmin=1.0, fmax=150.0,
                                             nperseg=60 * int(S[0].stats.sampling_rate),
                                             noverlap=None)

# then overlay in model script:
# plt.plot(freq_welch, psd_welch_db, color='orange', label='Observed Welch PSD')