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

def interp_hr(data, start, end):
    
    indat = data
    # Sort data by timestamp
    data = data.sort_values('timestamp')
    data['timestamp'] = data['timestamp'].dt.round('min')
    data = data.groupby('timestamp', as_index=False)[['lat', 'lon']].agg('mean')
    
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


def calc_dist(data, start, end, lon1, lon2, lat1, lat2):
    
    indat = data

    indat = indat[(indat['lon'] >= lon1) & (indat['lon'] <= lon2)] 
    indat = indat[(indat['lat'] >= lat1) & (indat['lat'] <= lat2)]

    print(f"{datetime.datetime.now()}: Interpolating by MMSI [2/4]")
    # Group by mmis and interpolate to hour
    outdat = indat.groupby('mmsi', as_index=False).apply(interp_hr, start, end)

    print(f"{datetime.datetime.now()}: Calculating NN [3/4]")
    # Calc dist.
    odat = outdat.groupby('timestamp', as_index=False).apply(spherical_dist_populate)

    return odat

if __name__ == "__main__":
    
    files = glob.glob("/home/server/pi/homes/woodilla/Data/GFW_point/Patagonia_Shelf/feather" + "/*.feather")

    nfiles = sorted(files)
    nfiles = nfiles[63:84]
    nfiles

    print(f"{datetime.datetime.now()}: Binding data [1/4]")
    list_ = []
    for file in nfiles:
        df = pd.read_feather(file)
        list_.append(df)
        mdat = pd.concat(list_)
    
      
    #start = pd.Timestamp("03-10-2016 00:00")
    #end = pd.Timestamp("03-20-2016 23:59")
    
    start_year = min(pd.DatetimeIndex(mdat['timestamp']).year)
    end_year = max(pd.DatetimeIndex(mdat['timestamp']).year)
    
    start_month = min(pd.DatetimeIndex(mdat['timestamp']).month)
    end_month = max(pd.DatetimeIndex(mdat['timestamp']).month)
    
    start_day = min(pd.DatetimeIndex(mdat['timestamp']).day)
    end_day = max(pd.DatetimeIndex(mdat['timestamp']).day)
    
    start = pd.Timestamp(f"{start_year} - {start_month} - {start_day} 00:00")
    end = pd.Timestamp(f"{end_year} - {end_month} - {end_day} 23:59")
    
    # Puerto Mardryn Area
    region = 1
    lon1 = -62.925222
    lon2 = -54.418267
    lat1 = -45
    lat2 = -40

    # Puerto Mardryn Area #2
#     region = 2
#     lon1 = -62.925222 - 1
#     lon2 = -54.418267 + 1
#     lat1 = -45 - 1
#     lat2 = -40 + 1

    # Puerto Mardryn Area #3
#     region = 3
#     lon1 = -62.925222 - 2
#     lon2 = -54.418267 + 2
#     lat1 = -45 - 2
#     lat2 = -40 + 2
    
    df = calc_dist(mdat, start=start, end=end, lon1=lon1, lon2=lon2, lat1=lat1, lat2=lat2)
    
    #print(mdat)
    print(f"{datetime.datetime.now()}: Saving: ~/Data/GFW_point/Patagonia_Shelf/complete/Puerto_Madryn_region{region}_{start_year}-{start_month}-{start_day}_{end_year}-{end_month}-{end_day}.feather [4/4]")
    df = df.reset_index(drop=True)
    df.to_feather(f"~/Data/GFW_point/Patagonia_Shelf/complete/Puerto_Madryn_region{region}_{start_year}-{start_month}-{start_day}_{end_year}-{end_month}-{end_day}.feather")
    
    
