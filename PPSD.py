import obspy
from obspy import read, UTCDateTime
from obspy.io.xseed import Parser
from obspy.signal import PPSD
import numpy as np
import matplotlib.pyplot as plt
from obspy.imaging.cm import pqlx
import os
from obspy.clients.fdsn import Client
from obspy.clients import fdsn
from fractions import Fraction

# function to plot PPSD and temporal plots
def plot_ppsd(S, t1, t2, stas, user_values):
    """Function to calculate and plot PPSD and Temporal Plot"""
    S.merge(method=1)  # merge traces to avoid duplicates
    # define path where .mseed files are stored

    # set station info
    nowsta = stas[0]
    network, station, location, channel = nowsta.split('.')
    if network == "ZE":
        S.resample(300) # resample
        client = Client("IRISPH5", timeout=600) # 3C data
    elif network == "ZD":
        S.resample(120) # resample
        client = Client("IRIS", timeout=600) # GEM data
    else:
        raise ValueError(f"Unknown network code: {network}")

    # metadata    
    inv = client.get_stations(network=network, station=station, location=location, channel=channel, level='response')

    # initialize PPSD object
    ppsd = None
    for tr in S:
        if ppsd is None:
            ppsd = PPSD(tr.stats, metadata=inv, db_bins=(-200,15,1.0),period_limits=(1/150, 10.0)) # if xaxis_frequency=False, set period_lim=(1/150,10)
        ppsd.add(tr)
        print('Trace added to PPSD')

    if ppsd is not None and len(ppsd.times_processed) > 0:
        print(f"Data accumulated for {(t2-t1)/(60*60*24)}-day periods starting on {t1}")
        ppsd.plot(period_lim=(0.1, 150), cmap=pqlx, xaxis_frequency=True)
        plt.close()
        #print(ppsd.times_processed)

        if user_values:
            print("Plotting temporal plot now...")
            ppsd.plot_temporal(user_values) # plot temporal plot with user values
            plt.show()
            plt.close()

    else:
        print(f"No PPSD data accumulated for {t1} to {t2}")

    return ppsd


# ---------------------------
# if running this script independently
# ---------------------------
if __name__ == '__main__':
    print("Running PPSD plotting script independently...")

    # define time range
    start_date = UTCDateTime(2023, 7, 23)
    end_date = UTCDateTime(2023, 7, 24)
    delta = 86400
    stas = [2301, 2302, 2303, 2304, 2305, 2306, 2307, 2308, 2309, 2310, 2311, 2312, 2313, 2314, 2315, 2316, 
        'G2301', 'G2302', 'G2303', 'G2304', 'G2305', 'G2306', 'G2307', 'G2308', 'G2309', 'G2310', 'G2311',
        'G2312', 'G2313', 'G2314', 'G2315', 'G2316'] 
    stas = ['ZE.2304..GP1'] # 3C station IDs
    #stas = ['ZD.G2411.01.HDF'] # GEM station IDs

    # set station info
    nowsta = stas[0]
    network, station, location, channel = nowsta.split('.')


    # can probably delete the whole client section here, it's in main function - same for spectraPlay?
    # can add path part to inside main function?

    # define path where .mseed files are stored
    if network == "ZE":
        path = f'C:/Users/zzawol/Documents/iris-data/seismic_data/NO{station}/{channel}'
        client = Client("IRISPH5", timeout=600) # 3C metadata
    elif network == "ZD":
        path = f'C:/Users/zzawol/Documents/iris-data/infrasound_data/{station}'
        client = Client("IRIS", timeout=600) # GEM metadata        ## not "IRISDMC" ?
    else:
        raise ValueError(f"Unknown network code: {network}")
    
    # create inventory for metadata  
    inv = client.get_stations(network=network, station=station, location=location, channel=channel, level='response')

    # check if the path exists, create if not
    if not os.path.exists(path):
        print(f"Path '{path}' does not exist. Creating directory...")
        os.makedirs(path, exist_ok=True)
    else:
        print(f"Path '{path}' already exists.")

    file_pattern = f'{path}/{network}.{station}.{location}.{channel}*.mseed'

    # read waveform
    S = read(file_pattern, starttime=start_date, endtime=end_date)

    # set frequency values for binning
    user_freqs = [120, 70, 20] # default frequencies in Hz
    user_periods = [1 / freq for freq in user_freqs] # convert to default periods for temporal plot

    # call function
    plot_ppsd(S, start_date, end_date, stas, user_values=user_periods)

# run in jupyter by: %run PPSD.py
# or run in terminal by: python PPSD.py






# call function independently - uncomment below for user input freqs

## OPTION 1: promopt for user input for temporal plot values
# likely wouldn't run this way without spectra since need to know freq bins for temporal plot

# user input frequency values
# user_inputs = input("Enter frequency values between 0.1 and 150 Hz (comma-separated, takes in fractions and decimals): ").split(',')

# # convert fractions and regular floats
# user_freqs = []
# user_periods = []
# for freq in user_inputs:
#     try:
#         freq = freq.strip()  # clean up extra spaces
#         if '/' in freq:
#             user_freqs.append(float(Fraction(freq)))  # handle fraction
#         else:
#             user_freqs.append(float(freq))  # handle float
#     except ValueError:
#         print(f"Invalid input: {freq}. Using default values.")

# # convert frequencies to periods
# user_periods = [1 / freq for freq in user_freqs]  # converts to periods for temporal plot
# #print(f'User-defined periods: {user_periods}')

# # if user didn't provide any valid values, set defaults
# if not user_freqs:
#     print("No valid input detected. Using default frequency values.")
#     user_freqs = [120, 70, 20]  # default frequencies in Hz
#     user_periods = [1 / freq for freq in user_freqs]  # convert to default periods for temporal plot

# # format periods for display
# round_periods = ['%.4f' % per for per in user_periods]
# print(f"Inputted frequencies: {user_freqs} --> periods: {round_periods}") # for checking


## OPTION 2: input temporal plot values within code
#user_freqs = [120, 70, 20] # default frequencies in Hz
#user_periods = [1 / freq for freq in user_freqs] # convert to default periods for temporal plot

# call function

#S = obspy.read(file_pattern, starttime=start_date, endtime=end_date)
#plot_ppsd(S, inv, t1=start_date, t2=end_date, user_values=user_periods) 
