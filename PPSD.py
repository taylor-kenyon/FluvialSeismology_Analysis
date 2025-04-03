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

# setting up parameters

# define start and end time for the data range
start_date = UTCDateTime(2023, 8, 20)  # Start time
end_date = UTCDateTime(2023, 8, 21)  # End time
delta = 86400  # 1 day in seconds
stas = [2301, 2302, 2303, 2304, 2305, 2306, 2307, 2308, 2309, 2310, 2311, 2312, 2313, 2314, 2315, 2316, 
        'G2301', 'G2302', 'G2303', 'G2304', 'G2305', 'G2306', 'G2307', 'G2308', 'G2309', 'G2310', 'G2311',
        'G2312', 'G2313', 'G2314', 'G2315', 'G2316'] 
stas = ['ZE.2303..GPZ']  # Station IDs

# set station info
nowsta = stas[0]
network, station, location, channel = nowsta.split('.')
print(f"Using station: {nowsta}")

# define path where .mseed files are stored
path = 'C:/Users/zzawol/Documents/seismic-data-iris/seismic_data/NO2303/GPZ'

# check if the path exists, create if not
if not os.path.exists(path):
    print(f"Path '{path}' does not exist. Creating directory...")
    os.makedirs(path, exist_ok=True)
else:
    print(f"Path '{path}' already exists.")

file_pattern = f'{path}/{network}.{station}.*.{channel}*.mseed'

# set up metadata
client = Client("IRISPH5", timeout=600)
#inv = client.get_stations(network='ZE', station=stas[0], channel='GPZ', level='response')  # adjust as needed
inv = client.get_stations(network=network, station=station, location=location, channel=channel, level='response') 

# function to plot ppsd and temporal plot
def plot_ppsd(S, inv, t1, t2, user_values):
    """Function to calculate and plot PPSD and Temporal Plot"""
    S.merge(method=1) # merge traces to avoid duplicates
    S.resample(300) # resample
    
    # initialize PPSD object
    ppsd = None
    for tr in S:
        if ppsd is None:
            ppsd = PPSD(tr.stats, metadata=inv, period_limits=(1/150, 10.0))
            # this format possibly avoids re-iterating over individual traces/re-creating PPSD object within loop
            #ppsd = PPSD(S[0].stats, metadata=inv, period_limits=(1/150, 10.0))  # adjust period limits

        ppsd.add(tr)
        print('Trace added to PPSD')

    # plot PPSD for the current day
    if ppsd is not None and len(ppsd.times_processed) > 0:
        print(f"Data accumulated for {(t2-t1)/(60*60*24)}-day periods starting on {t1}")
        ppsd.plot(period_lim=(0.1,150),cmap=pqlx,xaxis_frequency=True) # if xaxis_frequency=False, set period_lim=(1/150,10)
        plt.close() # prevents from displaying each iteration
        #print(ppsd.times_processed)
        
        # plot temporal plot if user values are provided
        if user_values:
            print(f"Plotting temporal plot for periods: {user_values}")
            ppsd.plot_temporal(user_values) # plot temporal plot with user values
            plt.show()
            plt.close()
            
        # save the plot as .png
        #plot_filename = f"ppsd_plot_{t1.date}.png" # change title if more than one date
        #plt.savefig(plot_filename)
        #print(f"Plot saved as {plot_filename}")

    else:
        print(f"No PPSD data accumulated for {t1} to {t2}")
        
    return ppsd

# call function independently
# likely wouldn't run this without spectra since need to know frequency bins for temporal plot

# user input frequency values
user_inputs = input("Enter frequency values between 0.1 and 150 Hz (comma-separated, takes in fractions and decimals): ").split(',')

# convert fractions and regular floats
user_freqs = []
user_periods = []
for freq in user_inputs:
    try:
        freq = freq.strip()  # clean up extra spaces
        if '/' in freq:
            user_freqs.append(float(Fraction(freq)))  # handle fraction
        else:
            user_freqs.append(float(freq))  # handle float
    except ValueError:
        print(f"Invalid input: {freq}. Using default values.")

# convert frequencies to periods
user_periods = [1 / freq for freq in user_freqs]  # converts to periods for temporal plot
#print(f'User-defined periods: {user_periods}')

# if user didn't provide any valid values, set defaults
if not user_freqs:
    print("No valid input detected. Using default frequency values.")
    user_freqs = [100, 10, 1]  # default frequencies in Hz
    user_periods = [1 / freq for freq in user_freqs]  # convert to default periods for temporal plot

# format periods for display
round_periods = ['%.4f' % per for per in user_periods]
print(f"Inputted frequencies: {user_freqs} --> periods: {round_periods}") # for checking

# call plot_ppsd function with user-provided values
S = obspy.read(file_pattern, starttime=start_date, endtime=end_date)
plot_ppsd(S, inv, t1=start_date, t2=end_date, user_values=user_periods)