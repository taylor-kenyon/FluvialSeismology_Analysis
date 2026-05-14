import os
import numpy as np
import matplotlib.pyplot as plt
from obspy import read, UTCDateTime, Stream
from scipy.signal import spectrogram
from matplotlib.ticker import MultipleLocator
from matplotlib.ticker import MaxNLocator, AutoMinorLocator
from glob import glob

# =========================================================
# CALIBRATION NODE
# =========================================================

def plot_calib_node(sta, calib_root, output_root):
    """
    Loads calibration node, plots waveform + spectrogram,
    and returns the experiment time window.
    """

    network, station, location, channel = sta.split('.')

    calib_sta = f"{network}.C{station}.{location}.{channel}"

    calib_pattern = os.path.join(
        calib_root,
        f"NOC{station}",
        channel,
        f"{calib_sta}_*.mseed"
    )

    print("\nCalibration search pattern:")
    print(calib_pattern)

    st = read(calib_pattern)
    st.merge(method=1)

    tr = st[0]

    print("\nCalibration trace:")
    print(tr)

    # experiment window
    test_start = tr.stats.starttime
    test_end = tr.stats.endtime

    # relative time axis
    time = tr.times(reftime=test_start)

    # output paths
    calib_dir = os.path.join(output_root, "calibration", station)
    os.makedirs(calib_dir, exist_ok=True)

    base_tag = f"{calib_sta}_{test_start.strftime('%Y%m%d%H%M%S')}"

    # =====================================================
    # WAVEFORM
    # =====================================================

    # plot waveform
    fig, ax = plt.subplots(figsize=(12, 4), dpi=600)

    ax.plot(time, tr.data, color='black', lw=0.5)

    ax.set_title(f"{calib_sta} Waveform")
    ax.set_xlabel("Time (s)")
    ax.set_ylabel("Amplitude")
    ax.xaxis.set_major_locator(MultipleLocator(120))
    ax.xaxis.set_minor_locator(MultipleLocator(30))
    ax.set_xlim(np.min(time), np.max(time))
    ax.set_ylim(np.min(tr.data)-10, np.max(tr.data)+10)
    fig.tight_layout()

    wave_path = os.path.join(calib_dir, f"{base_tag}_waveform.png")

    print("Saving:", wave_path)
    fig.savefig(wave_path, dpi=600, bbox_inches='tight')
    plt.close(fig)

    # =====================================================
    # SPECTROGRAM
    # =====================================================

    fs = tr.stats.sampling_rate
    nperseg = int(fs) # 1s FFT window
    #noverlap = int(fs * 0.5) # 50% overlap
    noverlap = int(0.80 * nperseg) # 80% overlap

    # plot spectrogram
    f, t, Sxx = spectrogram(
        tr.data,
        fs=fs,
        nperseg=nperseg,
        noverlap=noverlap
    )

    fig, ax = plt.subplots(figsize=(12, 6), dpi=600)

    pcm = ax.pcolormesh(
        t,
        f,
        10 * np.log10(Sxx + 1e-12),
        cmap='viridis',
        shading='auto'
    )

    ax.set_title(f"{calib_sta} Spectrogram")
    ax.set_xlabel("Time (s)")
    ax.set_ylabel("Frequency (Hz)")
    ax.set_ylim(0, 150)
    fig.colorbar(pcm, ax=ax, label="Power (dB)")

    ax.xaxis.set_major_locator(MultipleLocator(120)) # tick label every 10s
    ax.xaxis.set_minor_locator(MultipleLocator(30)) # tick unlabeled every 2s
    ax.tick_params(axis='x', which='minor', length=3.5)
    ax.tick_params(axis='x', which='major', length=6, width=1.1)
    ax.yaxis.set_major_locator(MultipleLocator(20)) # tick label every 20 Hz
    ax.yaxis.set_minor_locator(MultipleLocator(10)) # tick unlabeled every 10 Hz
    ax.tick_params(axis='y', which='minor', length=4)
    fig.tight_layout()

    spec_path = os.path.join(calib_dir, f"{base_tag}_spectrogram.png")

    print("Saving:", spec_path)
    fig.savefig(spec_path, dpi=600, bbox_inches='tight')
    plt.close(fig)

    return tr, calib_sta, test_start, test_end

# ***************************************************************************************************************************

# =========================================================
# ORIGINAL NODE
# =========================================================

def plot_og_node(sta, og_root, output_root, test_start, test_end):
    """
    Loads original node, trims to calibration window,
    plots waveform + spectrogram.
    """

    network, station, location, channel = sta.split('.')

    og_sta = f"{network}.{station}.{location}.{channel}"

    orig_dir = os.path.join(og_root, f"NO{station}", channel)

    # grab all candidate files, but not read
    files = glob(os.path.join(orig_dir, f"{og_sta}_*mseed"))
    print("\n")
    print("Total original node files found:", len(files))

    st = Stream()

    for f in files:
        st += read(f, starttime=test_start, endtime=test_end)

    if len(st) == 0:
        print("WARNING: No og node data found in calibration window")
        return

    st.merge(method=1)
    tr = st[0]

    # trim to experiment window
    #tr.trim(starttime=test_start, endtime=test_end)
    time = tr.times(reftime=test_start)

    # output paths
    og_dir = os.path.join(output_root, "original", station)
    os.makedirs(og_dir, exist_ok=True)

    base_tag = f"{og_sta}_{test_start.strftime('%Y%m%dT%H%M%S')}"

    # =====================================================
    # WAVEFORM
    # =====================================================

    # plot waveform
    fig, ax = plt.subplots(figsize=(12, 4), dpi=600)

    ax.plot(time, tr.data, color='black', lw=0.5)
    ax.set_title(f"{og_sta} Waveform")
    ax.set_xlabel("Time (s)")
    ax.set_ylabel("Amplitude")
    ax.xaxis.set_major_locator(MultipleLocator(120))
    ax.xaxis.set_minor_locator(MultipleLocator(30))
    ax.set_xlim(np.min(time), np.max(time))
    ax.set_ylim(np.min(tr.data)-10, np.max(tr.data)+10)
    fig.tight_layout()

    wave_path = os.path.join(og_dir, f"{base_tag}_waveform.png")

    print("Saving:", wave_path)
    fig.savefig(wave_path, dpi=600, bbox_inches='tight')
    plt.close(fig)

    # =====================================================
    # SPECTROGRAM
    # =====================================================

    fs = tr.stats.sampling_rate
    nperseg = int(fs)
    #noverlap = int(fs * 0.5) # 50% overlap
    noverlap = int(0.80 * nperseg) # 80% overlap

    # plot spectrogram
    f, t, Sxx = spectrogram(
        tr.data,
        fs=fs,
        nperseg=nperseg,
        noverlap=noverlap
    )

    fig, ax = plt.subplots(figsize=(12, 6), dpi=600)

    pcm = ax.pcolormesh(
        t,
        f,
        10 * np.log10(Sxx + 1e-12),
        cmap='viridis',
        shading='auto'
    )

    # ticks
    ax.xaxis.set_major_locator(MultipleLocator(120)) # tick label every 10s
    ax.xaxis.set_minor_locator(MultipleLocator(30)) # tick unlabeled every 2s
    ax.tick_params(axis='x', which='minor', length=3.5)
    ax.tick_params(axis='x', which='major', length=6, width=1.1)
    ax.yaxis.set_major_locator(MultipleLocator(20)) # tick label every 20 Hz
    ax.yaxis.set_minor_locator(MultipleLocator(10)) # tick unlabeled every 10 Hz
    ax.tick_params(axis='y', which='minor', length=4)
    
    ax.set_title(f"{og_sta} Spectrogram")
    ax.set_xlabel("Time (s)")
    ax.set_ylabel("Frequency (Hz)")
    ax.set_ylim(0, 150)
    fig.colorbar(pcm, ax=ax, label="Power (dB)")
    fig.tight_layout()

    spec_path = os.path.join(og_dir, f"{base_tag}_spectrogram.png")

    print("Saving:", spec_path)
    fig.savefig(spec_path, dpi=600, bbox_inches='tight')
    plt.close(fig)

    return tr, og_sta#, test_start

***************************************************************************************************************************
# =====================================================
# IMPACT WINDOWS
# =====================================================

def impact_window(tr, test_start, t1, t2, fs, out_dir, tag):
    """Function to plot waveform and spectrogram over a selected time window,
       - relative to global/calibration testing window."""

    # slice data
    sub = tr.slice(test_start + t1, test_start + t2)
    # global relative time axis
    time = sub.times(reftime=test_start)

    impact_dir = os.path.join(out_dir, "impact_windows", station)
    os.makedirs(impact_dir, exist_ok=True)

    base_tag = f"{tag}_{int(t1)}-{int(t2)}"

    window_start = test_start + t1
    window_end = test_start + t2

    base_tag = (f"{tag}_{window_start.strftime('%Y%m%d')}_{window_start.strftime('%H%M%S')}_{window_end.strftime('%H%M%S')}")

    # =====================================================
    # WAVEFORM
    # =====================================================

    # plot waveform
    fig, ax = plt.subplots(figsize=(12, 4), dpi=600)
    ax.plot(time, sub.data, color='black', lw=0.5)
    ax.set_xlim(t1, t2)
    #ax.set_ylim(np.min(sub.data)-30, np.max(sub.data)+30) # a bit too drastic/little when amps are very different
    ax.margins(y=0.05)
    ax.set_title(
        f"{tag} Waveform {t1}-{t2}s\n"
        f"Amp=({np.min(sub.data):.3g}, {np.max(sub.data):.3g})")
    ax.set_xlabel("Time (s)")
    ax.set_ylabel("Amplitude")

    # ticks
    ax.xaxis.set_major_locator(MultipleLocator(10)) # labeled every 10s
    ax.xaxis.set_minor_locator(MultipleLocator(2)) # unlabeled every 2s
    #ax.yaxis.set_major_locator(MultipleLocator(500)) # unlabeled every 500 amp
    #ax.yaxis.set_minor_locator(MultipleLocator(10)) # unlabeled every 10 amp
    ax.yaxis.set_major_locator(MaxNLocator(nbins=8)) # sets 8 major ticks labeled automatically
    ax.minorticks_on()
    # or: ax.yaxis.set_minor_locator(AutoMinorLocator())
    ax.tick_params(axis='x', which='minor', length=3.5)
    ax.tick_params(axis='x', which='major', length=6, width=1.1)
    ax.yaxis.set_major_locator(MultipleLocator(500)) # labeled every 500 magnitude
    fig.tight_layout()
    #plt.show()

    wave_path = os.path.join(impact_dir, f"{base_tag}_waveform.png")
    print("Saving waveform to:", wave_path)

    fig.savefig(wave_path, dpi=600, bbox_inches='tight')

    plt.close(fig)

    # =====================================================
    # SPECTROGRAM
    # =====================================================

    fs = tr.stats.sampling_rate
    # change back to shading = 'auto' ??
    nperseg = int(0.25 * fs) # 0.25*1000 = 250 samples = 0.25s FFT windows
    #noverlap = int(0.20 * fs) # 0.20*1000 = 200 samples of overlap (new FFT begins 50 samples later) = 0.05 ms
    noverlap = int(0.80 * nperseg) # this is same as commented out above but implies 80% overlap more clearly
    # i.e. 0.25(1 - 0.8) = 0.05s = 0

    f, spec_t, Sxx = spectrogram(
        sub.data,
        fs=fs,
        nperseg=nperseg,#int(fs),
        noverlap=noverlap#int(fs * 0.5)
    )

    # shift spectrogram time to global relative axis
    spec_t = spec_t + t1 # shift to global time

    # plot spectrogram
    fig, ax = plt.subplots(figsize=(12, 6), dpi=600)

    pcm = ax.pcolormesh(
        spec_t,
        f,
        10*np.log10(Sxx + 1e-12),
        cmap='viridis',
        shading="gouraud"  # or "nearest"
        )

    ax.set_xlim(t1, t2)
    ax.set_ylim(0, 150)
    ax.set_title(f"{tag} Spectrogram {t1}-{t2}s")
    ax.set_xlabel("Time (s)")
    ax.set_ylabel("Frequency (Hz)")
    fig.colorbar(pcm, ax=ax, label="Power (dB)")

    # ticks
    ax.xaxis.set_major_locator(MultipleLocator(10))
    ax.xaxis.set_minor_locator(MultipleLocator(2))
    ax.tick_params(axis='x', which='minor', length=3.5)
    ax.tick_params(axis='x', which='major', length=6, width=1.1)
    ax.yaxis.set_major_locator(MultipleLocator(20))
    ax.yaxis.set_minor_locator(MultipleLocator(10))
    ax.tick_params(axis='y', which='minor', length=4)

    fig.tight_layout()

    spec_path = os.path.join(impact_dir, f"{base_tag}_spectrogram.png")
    print("Saving spectrogram to:", spec_path)
    fig.savefig(spec_path, dpi=600, bbox_inches='tight')
    plt.close(fig)
    #plt.show()

# =========================================================
# MAIN
# =========================================================

if __name__ == "__main__":

    """To Use This Script:
        1) edit 'STA' below to select desired field location, calibration node and original node will plot automatically
            -- can change paths of where to save files as well
        2) run in terminal
        3) view raw waveform plots and spectrograms to decide impact windows of interest
        4) enter impact windows of interest (eg: "60-120") --> enter
        5) continue doing so with as many windows of interest as desired OR put q --> enter to finish
       """
    """Notes:
        -- full experiment window: 1s FFT window, 80% overlap
        -- impact set window: 0.25s FFT window, 80% overlap
        -- adjust as desired """

    # -----------------------------------------------------

    STA = "ZE.2516..GPZ"
    network, station, location, channel = STA.split('.')

    OG_ROOT = r"D:/2025_MSEED_From_Nodes_rename"
    CALIB_ROOT = r"D:/calibration_rename"

    OUTPUT_ROOT = r"C:/Users/zzawol/OneDrive - DOI/Documents/MORA25/FluvialSeismology_Analysis/calib_testing"

    # -----------------------------------------------------
    # CALIBRATION FIRST (defines experiment window)
    # -----------------------------------------------------

    trC, calib_sta, t0, t1 = plot_calib_node(STA, CALIB_ROOT, OUTPUT_ROOT)

    # -----------------------------------------------------
    # ORIGINAL NODE (trimmed to same window)
    # -----------------------------------------------------

    trO, og_sta = plot_og_node(STA, OG_ROOT, OUTPUT_ROOT, t0, t1) #, t0O

    # -----------------------------------------------------
    # SELECT IMPACT WINDOW
    # -----------------------------------------------------

    while True:
        print("Enter time windows in seconds (hyphen-separated) to visualize data more closely.")
        user_input = input("Example: 160-250, 1040-1075"
        "\nor type 'q' to quit\n")

        if user_input.lower() == 'q':
            break

        try:
            windows = []
            for item in user_input.split(','):
                start, end = item.split('-')
                windows.append((float(start), float(end)))

            for w_start, w_end in windows:

                # calibration impact window
                impact_window(tr=trC, test_start=t0, t1=w_start, t2=w_end, fs=trC.stats.sampling_rate, out_dir=OUTPUT_ROOT, tag=calib_sta)
                # original impact window
                impact_window(tr=trO, test_start=t0, t1=w_start, t2=w_end, fs=trO.stats.sampling_rate, out_dir=OUTPUT_ROOT, tag=og_sta)

        except Exception as e:
            print("Invalid input:", e)


    print("\nDONE.")