import obspy
from obspy import read, UTCDateTime
from obspy.clients import fdsn
import os

# function to get data from IRIS
def get_iris_data(t1, t2, stas, path):
    """Function to get data from IRIS if not already downloaded"""

    # set up IRIS client
    DATASELECT = 'http://service.iris.edu/ph5ws/dataselect/1'
    c = fdsn.client.Client(
        service_mappings={
            'dataselect': DATASELECT,
        },
    )

    # IRIS credentials
    username = 'zoe_zawol@partner.nps.gov'
    password = 'rJXKed4LZUHUE05g'
    c.set_credentials(username, password)

    # check if path exists, create if not
    if not os.path.exists(path):
        print(f"Path '{path}' does not exist. Creating directory...")
        os.makedirs(path, exist_ok=True)
    else:
        print(f"Path '{path}' already exists.")

    # set station info
    nowsta = stas[0]
    network, station, location, channel = nowsta.split('.')
    print(f"Station: {nowsta}")
    file_pattern = f'{path}/{network}.{station}.*.{channel}*.mseed'

    tNow = t1 # start time
    extract_delta = 1800 # 30 mins in seconds

    while tNow < t2:
        for nowsta in stas:
            print(f"Checking data for station {nowsta} from {tNow} to {tNow + extract_delta}")

            # create list of existing .mseed files that already exist in the path
            existing_files = [f for f in os.listdir(path) if f.endswith('.mseed')]

            # construct filename pattern to match for this time period
            filename_pattern = f'{network}.{station}..{channel}_{tNow.strftime("%Y%m%d%H%M%S")}.mseed'

            # check if any required .mseed files are already downloaded
            file_exists = False
            for filename in existing_files:
                if filename_pattern in filename:
                    file_exists = True
                    break

            if file_exists:
                print(f"Data for {nowsta} from {tNow} to {tNow + extract_delta} already stored locally. Skipping download.")
            else:
                print(f"Retrieving data for {nowsta} from {tNow} to {tNow + extract_delta}")
                try:
                    # retrieve data from IRIS
                    S = c.get_waveforms(network, station, location, channel, tNow, tNow + extract_delta)
                    print(f"Data for {nowsta} retrieved successfully.")

                    print('Resampling data...')
                    S.resample(350)  # resample to 350 Hz/s

                    # save data to .mseed files
                    print('Saving data to .mseed files...')
                    for tr in S:
                        # create filename using starttime of trace with .mseed file format
                        filename = os.path.join(path, f'{tr.id}_{tr.stats.starttime.strftime("%Y%m%d%H%M%S")}.mseed')
                        tr.write(filename)
                        print(f"Saved {filename}")
                except Exception as e:
                    print(f"Error retrieving data for {nowsta} from {tNow} to {tNow + extract_delta}: {e}")
        
        tNow += extract_delta  # move to next interval

    print('Data retrieval and saving completed.')


# format for running independently 
# these parameteres will be ignored when Screening_Tool calls getData since we will set the parameters in Screening_Tool
# so we dont have to comment and uncomment when using this file in Screening_Tool vs independently 

if __name__ == "__main__":
    # set parameters here for standalone testing
    print("using function direct method")
    start_date = UTCDateTime(2023, 8, 14, 14) # Start time
    end_date = UTCDateTime(2023, 8, 14, 18) # End time
    stas = [2301, 2302, 2303, 2304, 2305, 2306, 2307, 2308, 2309, 2310, 2311, 2312, 2313, 2314, 2315, 2316, 
            'G2301', 'G2302', 'G2303', 'G2304', 'G2305', 'G2306', 'G2307', 'G2308', 'G2309', 'G2310', 'G2311',
            'G2312', 'G2313', 'G2314', 'G2315', 'G2316'] 
    stas = ['ZE.2301..GPZ'] # station IDs
    nowsta = stas[0]
    network, station, location, channel = nowsta.split('.')
    path = f'C:/Users/zzawol/Documents/seismic-data-iris/seismic_data/NO{station}/{channel}'
    get_iris_data(t1=start_date, t2=end_date, stas=stas, path=path)

# run in jupyter by: %run getData.py
# or run in terminal by: python getData.py