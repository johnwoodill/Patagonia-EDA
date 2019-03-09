import pandas as pd
import numpy as np
import glob
from scipy import stats
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

# Get list of files in folder
files = []
for file in glob.glob("/home/server/pi/homes/woodilla/Data/GFW_point/Patagonia_Shelf/feather" + "/*.feather"):
    files.append(file)
    nfiles = sorted(files)

if __name__ == '__main__':

    # Loop through processed files
    for i in nfiles:
        indat = pd.read_feather(i)
        indat = indat.sort_values('mmsi')
        outdat = indat
        
        # Average lat/long
        outdat['lat_avg'] = outdat.groupby('mmsi').lat.transform('mean')
        outdat['lon_avg'] = outdat.groupby('mmsi').lon.transform('mean')
        outdat = outdat.groupby('mmsi').first().reset_index()
        
        # Organize data
        outdat = outdat[['timestamp', 'year', 'month', 'day', 'mmsi', 'lat', 'lon', \
                        'segment_id', 'message_id', 'type', 'speed', 'course', 'heading', 'shipname', 'callsign', \
                         'destination', 'elevation_m', 'distance_from_shore_m', 'distance_from_port_m', 'nnet_score', \
                         'logistic_score', 'flag', 'geartype', 'length', 'tonnage', 'engine_power', 'active_2012', \
                         'active_2013', 'active_2014', 'active_2015', 'active_2016']]
        
        # f-string for date
        date = f"{outdat['year'][1]}-" + f"{outdat['month'][1]}".zfill(2) + f"-" + f"{outdat['day'][1]}".zfill(2)
        
        # Build distance matrix
        matdat = pd.DataFrame(spherical_dist_populate(outdat['lat'], outdat['lon']))
        matdat = matdat.rename(index=outdat.mmsi, columns = outdat.mmsi)    
        
        # Collect, stack, and form three column data.frame 
        tmatdat = matdat.where(np.triu(np.ones(matdat.shape)).astype(np.bool))
        tmatdat = tmatdat.stack().reset_index()
        
        # Add date
        tmatdat['date'] = date
        
        # Column names: vessel_A = reference boat; vessel_B = boat to calc distance to; distance = distance
        tmatdat.columns = ['vessel_A', 'vessel_B', 'distance', 'date']
        tmatdat = tmatdat[['date', 'vessel_A', 'vessel_B', 'distance']]
        
        # Save
        tmatdat.to_feather('/home/server/pi/homes/woodilla/Data/GFW_point/Patagonia_Shelf/distance_data/circle_measure/' + date + '_dmatrix' + '.feather')
