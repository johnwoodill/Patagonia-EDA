import pandas as pd
import numpy as np 
import feather
import os as os
#os.environ["PROJ_LIB"] = "/Users/john/miniconda3/share/proj"; #fixr
#from mpl_toolkits.basemap import Basemap
import matplotlib.pyplot as plt
import glob
import cartopy.crs as ccrs
import cartopy.io.img_tiles as cimgt
import cartopy
import datetime

def spherical_dist_populate(data=None, lat_lis=None, lon_lis=None, r=3958.75):
    
    if data is not None:
        data = data.dropna()
        mmsi = data.mmsi
        data = data.sort_values('mmsi')
        lat_lis = data['lat']
        lon_lis = data['lon']
        timestamp = data['timestamp'].iat[0]
    
    
    lat_mtx = np.array([lat_lis]).T * np.pi / 180
    lon_mtx = np.array([lon_lis]).T * np.pi / 180

    cos_lat_i = np.cos(lat_mtx)
    cos_lat_j = np.cos(lat_mtx)
    cos_lat_J = np.repeat(cos_lat_j, len(lat_mtx), axis=1).T

    lat_Mtx = np.repeat(lat_mtx, len(lat_mtx), axis=1).T
    cos_lat_d = np.cos(lat_mtx - lat_Mtx)

    lon_Mtx = np.repeat(lon_mtx, len(lon_mtx), axis=1).T
    cos_lon_d = np.cos(lon_mtx - lon_Mtx)

    mtx = r * np.arccos(cos_lat_d - cos_lat_i*cos_lat_J*(1 - cos_lon_d))
    
    # Build data.frame
    matdat = pd.DataFrame(mtx)
    matdat.columns = mmsi[:]
    matdat = matdat.set_index(mmsi[:])
    
    # Stack and form three column data.frame 
    tmatdat = matdat.stack()
    lst = tmatdat.index.tolist()
    vessel_A = pd.Series([item[0] for item in lst])
    vessel_B = pd.Series([item[1] for item in lst])
    distance = tmatdat.values

    # Get lat/lon per mmsi
    posdat = data[['mmsi', 'lat', 'lon']]
    posdat = posdat.sort_values('mmsi')
    
    # Build data frame
    odat = pd.DataFrame({'timestamp': timestamp, 'vessel_A': vessel_A, 'vessel_B': vessel_B, 'distance': distance})
    odat = odat.sort_values(['vessel_A', 'distance'])
        
    # Get 10-NN
    odat = odat.sort_values('distance').groupby('vessel_A', as_index=False).nth([0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10])
    odat = odat.sort_values(['vessel_A', 'distance'])
    
    # Merge in vessel_B lat/lon
    posdat.columns = ['mmsi', 'vessel_B_lat', 'vessel_B_lon']
    odat = odat.merge(posdat, how='left', left_on='vessel_B', right_on='mmsi')
    
    # Merge in vessel_A lat/lon
    posdat.columns = ['mmsi', 'vessel_A_lat', 'vessel_A_lon']
    odat = odat.merge(posdat, how='left', left_on='vessel_A', right_on='mmsi')
    
    odat['rank'] = odat.groupby(['vessel_A'], as_index=False).cumcount()
    odat = odat.reset_index(drop=True)
    odat = odat[['timestamp', 'vessel_A', 'vessel_B', 'vessel_A_lat', 'vessel_A_lon', 'vessel_B_lat', 'vessel_B_lon', 'rank', 'distance']]
    odat = odat.sort_values(['vessel_A', 'rank'])
    
    return odat

def interp_hr(data):
    
    indat = data
    # Sort data by timestamp
    data = data.sort_values('timestamp')
    data['timestamp'] = data['timestamp'].dt.round('min')
    data = data.groupby('timestamp', as_index=False)[['lat', 'lon']].agg('mean')
    
    # Begin and end date data frame
    #start = pd.Timestamp(data['timestamp'].iat[0]).strftime('%m-%d-%Y 00:00')
    #end = pd.Timestamp(data['timestamp'].iat[-1]).strftime('%m-%d-%Y  23:59')
    
    start = pd.Timestamp("03-10-2016 00:00")
    end = pd.Timestamp("03-20-2016 23:59")
    
    #start = pd.Timestamp(f"{}")
    
    # Merge and interpolate between start and end
    pdat = pd.DataFrame({'timestamp': pd.date_range(start=start, end=end, freq='min')})
    
    pdat = pdat.merge(data, on='timestamp', how='left')
    pdat['lat'] = pd.Series(pdat['lat']).interpolate()
    pdat['lon'] = pd.Series(pdat['lon']).interpolate()
    
    # Keep on the hour
    pdat = pdat[pd.Series(pdat['timestamp']).dt.minute == 00].reset_index(drop=True)
    pdat = pdat.fillna(method='bfill')
    pdat = pdat.fillna(method='ffill')
    pdat['mmsi'] = indat['mmsi'].iat[0]
    
    pdat = pdat.reset_index(drop=True)
    #print(pdat)
    return pdat


def calc_dist(data):
    
    #indat = pd.read_feather(data_loc)
    indat = data
    #indat = indat.sort_values('mmsi')
    
#     # Puerto Mardryn Area
#     lon1 = -62.9
#     lon2 = -54.4
#     lat1 = -45
#     lat2 = -40
    
    lon1 = -62.925222
    lon2 = -54.418267
    lat1 = -45
    lat2 = -40
    indat = indat[(indat['lon'] >= lon1) & (indat['lon'] <= lon2)] 
    indat = indat[(indat['lat'] >= lat1) & (indat['lat'] <= lat2)]

    print(f"{datetime.datetime.now()}: Interpolating by MMSI [2/4]")
    # Group by mmis and interpolate to hour
    outdat = indat.groupby('mmsi', as_index=False).apply(interp_hr)

    print(f"{datetime.datetime.now()}: Calculating NN [3/4]")
    # Calc dist.
    odat = outdat.groupby('timestamp', as_index=False).apply(spherical_dist_populate)

    return odat

if __name__ == "__main__":
    
    files = glob.glob("/home/server/pi/homes/woodilla/Data/GFW_point/Patagonia_Shelf/feather" + "/*.feather")

    nfiles = sorted(files)
    nfiles = nfiles[68:79]
    
    print(f"{datetime.datetime.now()}: Binding data [1/4]")
    list_ = []
    for file in nfiles:
        df = pd.read_feather(file)
        list_.append(df)
        mdat = pd.concat(list_)
        

    df = calc_dist(mdat)    
    
    #print(mdat)
    print(f"{datetime.datetime.now()}: Saving [4/4]")
    df = df.reset_index(drop=True)
    df.to_feather('~/Data/GFW_point/Patagonia_Shelf/complete/Puerto_Madryn_2016-03-10_2016-03-20.feather')