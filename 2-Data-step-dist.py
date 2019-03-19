import pandas as pd
import numpy as np
import glob
from scipy import stats
import multiprocessing
from scipy.stats import kurtosis, skew


# Calculate spherical distances between lat/long and populate matrix
# Returns distance in miles
def spherical_dist_populate(lat_lis, lon_lis, r=3958.75):
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
    return mtx


def mp_gfw_dist(data_loc):
    
    indat = pd.read_feather(data_loc)
    indat = indat.sort_values('mmsi')
    outdat = indat
    
    # Average lat/long
    outdat['lat_avg'] = outdat.groupby('mmsi').lat.transform('mean')
    outdat['lon_avg'] = outdat.groupby('mmsi').lon.transform('mean')
    outdat = outdat.groupby('mmsi').first().reset_index()

    # (4) Determine if stationary where distance_traveled > 1
    outdat['stationary'] = np.where(outdat['dist'] > 1, 0, 1)

    # Organize data
    outdat = outdat[['timestamp', 'year', 'month', 'day', 'mmsi', 'lat_avg', 'lon_avg', \
                     'segment_id', 'message_id', 'type', 'speed', 'mph', 'dist', 'travel_time', 'stationary', \
                     'course', 'heading', 'shipname', 'callsign', \
                     'destination', 'elevation_m', 'distance_from_shore_m', 'distance_from_port_m', 'nnet_score', \
                     'logistic_score', 'flag', 'geartype', 'length', 'tonnage', 'engine_power', 'active_2012', \
                     'active_2013', 'active_2014', 'active_2015', 'active_2016']]
        
    # f-string for date
    date = f"{outdat['year'][1]}-" + f"{outdat['month'][1]}".zfill(2) + f"-" + f"{outdat['day'][1]}".zfill(2)

    # Get avg lat/lon per mmsi
    posdat = outdat[['mmsi', 'lat_avg', 'lon_avg']]
    posdat = posdat.groupby('mmsi', as_index=False).mean()
    posdat = posdat.sort_values('mmsi')

    # Build distance matrix
    matdat = pd.DataFrame(spherical_dist_populate(outdat['lat_avg'], outdat['lon_avg']))
    matdat = matdat.rename(index=outdat.mmsi, columns = outdat.mmsi)    
        
    # Stack and form three column data.frame 
    tmatdat = matdat.stack()
    lst = tmatdat.index.tolist()
    vessel_A = [item[0] for item in lst]
    vessel_B = [item[1] for item in lst]
    distance = tmatdat.values

    # Build data frame
    odat = pd.DataFrame({'date': date, 'vessel_A': vessel_A, 'vessel_B':vessel_B, 'distance': distance})
    odat = odat.sort_values(['vessel_A', 'distance'])
    #odat = odat.reset_index(drop=True)
    
    # Get 10-NN
    odat = odat.sort_values(['vessel_A', 'distance'])
    odat = odat.sort_values('distance').groupby('vessel_A', as_index=False).nth([1,2, 3, 4, 5, 6, 7, 8, 9, 10])
    odat = odat.sort_values(['vessel_A', 'distance'])
    
    # Merge in vessel_B lat/lon
    posdat.columns = ['mmsi', 'vessel_B_lat', 'vessel_B_lon']
    odat = odat.merge(posdat, how='left', left_on='vessel_B', right_on='mmsi')
    
    # Merge in vessel_A lat/lon
    posdat.columns = ['mmsi', 'vessel_A_lat', 'vessel_A_lon']
    odat = odat.merge(posdat, how='left', left_on='vessel_A', right_on='mmsi')
    

    odat['rank'] = odat.groupby(['vessel_A'], as_index=False).cumcount() + 1
    odat = odat.reset_index(drop=True)
    odat = odat[['date', 'vessel_A', 'vessel_B', 'vessel_A_lat', 'vessel_A_lon', 'vessel_B_lat', 'vessel_B_lon', 'rank', 'distance']]
    odat.to_feather('/home/server/pi/homes/woodilla/Data/GFW_point/Patagonia_Shelf/distance_data/NearestNeighbor/by_date/' + odat['date'][1] + '_NN.feather')
    print(odat['date'][0])
    return 0
            
if __name__ == '__main__':

    # Get list of files in folder
    
    files = glob.glob("/home/server/pi/homes/woodilla/Data/GFW_point/Patagonia_Shelf/feather" + "/*.feather")
    
    nfiles = sorted(files)
    
    # Multiprocess files
    pool = multiprocessing.Pool(10, maxtasksperchild=1)         
    pool.map(mp_gfw_dist, nfiles)
    pool.close()
    
    print('By date processing complete. Begin processing unique mmsi')
    
    subdir = '/home/server/pi/homes/woodilla/Data/GFW_point/Patagonia_Shelf/distance_data/NearestNeighbor/by_date/'
    allFiles = glob.glob(subdir + "*.feather")
    list_ = []

    # Append files in subdir
    for file_ in allFiles:
        df = pd.read_feather(file_)
        list_.append(df)
        dat = pd.concat(list_, axis = 0, ignore_index = True)  
        
    dat.to_feather('~/Data/GFW_point/Patagonia_Shelf/complete/NN_2016-2018.feather')

    grouped = dat.groupby('vessel_A')
        
    for name, group in grouped:
        mmsi = name
        days = len(group['date'].unique())
        group = group.sort_values('date')
        group = group.reset_index(drop=True)
        #print(group)
        
        if days > 5:
            filename = f"{mmsi}_{days}_days.feather" 
            print(filename)
            group.to_feather('~/Data/GFW_point/Patagonia_Shelf/distance_data/NearestNeighbor/by_mmsi/' + filename)

    

        