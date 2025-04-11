import obspy 
from obspy import read, UTCDateTime
import numpy as np
import matplotlib.pyplot as plt
import os
import waveformUtils

# function to plot spectrogram
def plot_spectrogram(S, t1, t2):
    """Function to plot Waveforms, RSAM, and Spectrogram"""

    # f, spectraF, spectraG, sttimes, tvec, RSAM = waveformUtils.plotSpectraTime(S, minfreq=1, maxfreq=100, winlength=7200, step=7200,
    #                                                                            spacing=2, ampScalar=3, normalize=True,
    #                                                                            specWin=60, specStep=30, rsamWin=1, saveFig=False)
    
    # plt.figure()
    # f.show()
    # remember to also return f when wanting this plot

    S.merge(method=1)  # merge traces to avoid duplicates

    f1 = waveformUtils.multiDaySpectrogram(S, averageLength=3600, fftLength=60, minFreq=0.5, maxFreq=150, cmap='turbo',
                                           dateLimits=None, vmin=0.4, vmax=0.9, plotAverage=True)

    plt.figure()
    #f1.show()  # comment this out when running file independently, uncomment for whole screening process
    plt.show(block=False)
    plt.close()

    return f1  #, f


# ---------------------------
# if running this script independently
# ---------------------------
if __name__ == "__main__":
    # setting up parameters for standalone use

    # define start and end time for the data range
    start_date = UTCDateTime(2023, 7, 23)  # Start time
    end_date = UTCDateTime(2023, 7, 24)    # End time
    delta = 86400  # 1 day in seconds

    stas = [2301, 2302, 2303, 2304, 2305, 2306, 2307, 2308, 2309, 2310, 2311, 2312, 2313, 2314, 2315, 2316, 
            'G2301', 'G2302', 'G2303', 'G2304', 'G2305', 'G2306', 'G2307', 'G2308', 'G2309', 'G2310', 'G2311',
            'G2312', 'G2313', 'G2314', 'G2315', 'G2316'] 
    stas = ['ZE.2304..GP1']  # station IDs

    # set station info
    nowsta = stas[0]
    network, station, location, channel = nowsta.split('.')
    print(f"Station: {nowsta}")

    # define path where .mseed files are stored
    path = f'C:/Users/zzawol/Documents/seismic-data-iris/seismic_data/NO{station}/{channel}'

    # check if the path exists, create if not
    if not os.path.exists(path):
        print(f"Path '{path}' does not exist. Creating directory...")
        os.makedirs(path, exist_ok=True)
    else:
        print(f"Path '{path}' already exists.")

    file_pattern = f'{path}/{network}.{station}.*.{channel}*.mseed'

    # read and run plot
    end = end_date + (6*60)
    S = obspy.read(file_pattern, starttime=start_date, endtime=end)
    plot_spectrogram(S, t1=start_date, t2=end)