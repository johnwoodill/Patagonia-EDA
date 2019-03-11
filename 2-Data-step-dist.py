import pandas as pd
import numpy as np
import glob
from scipy import stats
import multiprocessing
from scipy.stats import kurtosis, skew

# Calculate spherical distances between lat/long and populate matrix
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
        
    # Organize data
    outdat = outdat[['timestamp', 'year', 'month', 'day', 'mmsi', 'lat_avg', 'lon_avg', \
                     'segment_id', 'message_id', 'type', 'speed', 'course', 'heading', 'shipname', 'callsign', \
                      'destination', 'elevation_m', 'distance_from_shore_m', 'distance_from_port_m', 'nnet_score', \
                      'logistic_score', 'flag', 'geartype', 'length', 'tonnage', 'engine_power', 'active_2012', \
                      'active_2013', 'active_2014', 'active_2015', 'active_2016']]
        
    # f-string for date
    date = f"{outdat['year'][1]}-" + f"{outdat['month'][1]}".zfill(2) + f"-" + f"{outdat['day'][1]}".zfill(2)
        
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
    odat = odat.reset_index(drop=True)
            
    # Save
    odat.to_feather('/home/server/pi/homes/woodilla/Data/GFW_point/Patagonia_Shelf/distance_data/circle_measure/' + date + '_dmatrix' + '.feather')
    return 0

if __name__ == '__main__':

    # Get list of files in folder
    files = glob.glob("/home/server/pi/homes/woodilla/Data/GFW_point/Patagonia_Shelf/feather" + "/*.feather")
    nfiles = sorted(files)
    
    # Multiprocess files
    pool = multiprocessing.Pool(10, maxtasksperchild=1)         
    pool.map(mp_gfw_dist, nfiles)
    pool.close()
    

        