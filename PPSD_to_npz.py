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
        S.resample(300) # resample from 350 to 300
        client = Client("IRISPH5", timeout=600) # 3C data
    elif network == "ZD":
        S.resample(100) # resample - not necessary but to make sure
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
        # save mean PSD in an .npz
        freq, mean_psd = ppsd.get_mean()
        npz_path = f'C:/Users/zzawol/OneDrive - DOI/Documents/MORA25/iris-data/seismic_data/NO{station}/npzs'
        if not os.path.exists(npz_path):
            print(f"Path '{npz_path}' does not exist. Creating directory...")
            os.makedirs(npz_path, exist_ok=True)
        else:
            print(f"Path '{npz_path}' already exists.")

        out_name = f"{npz_path}/{network}.{station}.{location}.{channel}_{start_date.date}-{end_date.date}.npz"
        # for file anming maybe use along the lines of before:
        #      f'{tr.id}_{tr.stats.starttime.strftime("%Y%m%d%H%M%S")}
        np.savez(out_name, freq=freq, mean_psd=mean_psd)
        print(f"Saved mean PSD to {out_name}")
        print(f"Data accumulated for {(t2-t1)/(60*60*24)}-day periods starting on {t1}")
        ppsd.plot(period_lim=(0.1,150), cmap=pqlx, xaxis_frequency=True)
        plt.close()

        
        # psd_shape = ppsd.psd_values
        # print("vals shape:", psd_shape.shape)
        # psd_len = ppsd.extract_psd_values
        # print("vals length:", len(psd_len))
        # print("ppsd type:", type(ppsd))
        # print("psd type:", type(ppsd.psd_values))
        # print("shape:", ppsd.shape)
        #print(ppsd.times_processed)

        if user_values:
            print("Plotting temporal plot now...")
            ppsd.plot_temporal(user_values) # plot temporal plot with user values
            plt.show()
            plt.close()

        #round_periods = ['%.4f' % per for per in user_periods]
        #print(f"Inputted frequencies: {user_freqs} --> periods: {round_periods}")

    else:
        print(f"No PPSD data accumulated for {t1} to {t2}")


### adding this for npz's ###########################################################
    # Save PPSD data to .npz for later plotting
#         out_dir = "C:/Users/iris-data/ppsd_npz"  # or wherever you want to store them
#         os.makedirs(out_dir, exist_ok=True)

#         out_file = f"{network}.{station}.{location}.{channel}_{t1.date}_{t2.date}_ppsd.npz"
#         out_path = os.path.join(out_dir, out_file)

#         np.savez(
#             out_path,
#             freqs=ppsd.period_bin_centers,
#             psd=ppsd.get_mean(),
#             times=np.array(ppsd.times_processed)
#         )
#         print(f"PPSD data saved to {out_path}")

# #################### then run this in model_bedload.py #########################
# import numpy as np
# import matplotlib.pyplot as plt

# data = np.load("path/to/your_ppsd_file.npz")
# freqs = data["freqs"]
# psd = data["psd"]

# plt.plot(1 / freqs, psd, color="orange")  # Convert periods to frequency
# plt.xscale("log")
# plt.xlabel("Frequency (Hz)")
# plt.ylabel("Power (dB)")
# plt.title("Mean PSD from PPSD")
# plt.grid(True)
# plt.show()


    return ppsd


# ---------------------------
# if running this script independently
# ---------------------------

if __name__ == '__main__':
    print("Running PPSD plotting script independently...")
    
    # setting up parameters for independent use

    # define time range
    start_date = UTCDateTime(2024, 7, 26, 12)
    end_date = UTCDateTime(2024, 7, 26, 13)
    delta = 12*60*60 #86400
    stas = [2301, 2302, 2303, 2304, 2305, 2306, 2307, 2308, 2309, 2310, 2311, 2312, 2313, 2314, 2315, 2316, 
        'G2301', 'G2302', 'G2303', 'G2304', 'G2305', 'G2306', 'G2307', 'G2308', 'G2309', 'G2310', 'G2311',
        'G2312', 'G2313', 'G2314', 'G2315', 'G2316'] 
    stas = ['ZE.2416..GPZ'] # 3C station IDs
    #stas = ['ZD.G2411.01.HDF'] # GEM station IDs

    # set station info
    nowsta = stas[0]
    network, station, location, channel = nowsta.split('.')

    # define path where .mseed files are stored and set frequencies for binning
    if network == "ZE":
        # define path to store 3C files
        path = f'C:/Users/zzawol/OneDrive - DOI/Documents/MORA25/iris-data/seismic_data/NO{station}/{channel}'
        # set frequency values for binning
        user_freqs = [120, 70, 20] # default seismic frequencies
        user_periods = [1 / freq for freq in user_freqs] # convert to periods
    elif network == "ZD":
        # define path where GEM files are stored
        path = f'C:/Users/zzawol/Documents/OneDrive - DOI/Documents/MORA25/iris-data/infrasound_data/{station}'
        # set frequency values for binning
        user_freqs = [1, 10, 20] 
        user_periods = [1 / freq for freq in user_freqs] # convert to periods
    else:
        raise ValueError(f"Unknown network code: {network}")
    
    # check if the path exists, create if not
    if not os.path.exists(path):
        print(f"Path '{path}' does not exist. Creating directory...")
        os.makedirs(path, exist_ok=True)
    else:
        print(f"Path '{path}' already exists.")

    file_pattern = f'{path}/{network}.{station}.{location}.{channel}*.mseed'

    # ## loop through time windows
    # current_time = start_date
    # while current_time < end_date:
    #     next_time = current_time + delta
    #     print("\nProcessing: {current_time.date} to {next_time.date}")

    #     match_files = []
    #     for fname in os.listdir(path): # or file_pattern?
    #         if fname.endswith(".mseed") and fname.startswith(f"{network}.{station}.{location}.{channel}_"):
    #             try:
    #                 timestamp_str = fname.split("_")[-1].replace(".mseed", "")
    #                 file_time = UTCDateTime(timestamp_str)
    #                 #print(f"File: {fname} --> parsed time: {file_time}")
    #                 if current_time <= file_time < next_time:
    #                     match_files.append(os.path.join(path, fname))
    #             except Exception as e:
    #                 print(f"Skipping file {fname}: {e}")

    #         if not match_files:
    #             print("No files found in this window.")
    #             current_time = next_time
    #             continue
            
    #         try:
    #             S = read(match_files[0])
    #             for f in match_files[1]:
    #                 S += read(f)
    #             print(f"Read {len(match_files)} files.")
    #         except Exception as e:
    #             print(f"Error reading files: {e}")
    #             current_time = next_time
    #             continue
            
    #         plot_ppsd(S, current_time, next_time, [nowsta], user_values=user_periods)
    #         current_time = next_time


    ### without time filtering/looping
    mseed_files = []
    for fname in os.listdir(path):
        if fname.endswith(".mseed") and fname.startswith(f"{nowsta}_"):
            try:
                timestamp_str = fname.split("_")[-1].replace(".mseed", "")
                file_time = UTCDateTime(timestamp_str)
                if start_date <= file_time < end_date:
                    mseed_files.append(os.path.join(path, fname))
            except Exception as e:
                print(f"Skipping file {fname}: {e}")

    if not mseed_files:
        print("No matching files found in the given time range.")
        exit()

    print(f"Found {len(mseed_files)} files:")
    for f in mseed_files:
        print(" -", os.path.basename(f))

    # read data
    S = read(mseed_files[0])
    for f in mseed_files[1:]:
        S += read(f)
        S.merge(method=1)














    # read waveform
    S = read(file_pattern, starttime=start_date, endtime=end_date)

    round_periods = ['%.4f' % per for per in user_periods]
    print(f"Inputted frequencies: {user_freqs} --> periods: {round_periods}")

    # call function
    plot_ppsd(S, start_date, end_date, stas, user_values=user_periods)

# run in jupyter by: %run PPSD.py
# or run in terminal by: python PPSD.py
