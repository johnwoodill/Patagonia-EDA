import pandas as pd
import numpy as np 
import feather
import os as os
#os.environ["PROJ_LIB"] = "/Users/john/miniconda3/share/proj"; #fixr
#from mpl_toolkits.basemap import Basemap
import matplotlib.pyplot as plt
#import glob
%matplotlib inline
import cartopy.crs as ccrs
import cartopy.io.img_tiles as cimgt
import cartopy

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
    start = pd.Timestamp(data['timestamp'].iat[0]).strftime('%m-%d-%Y %H:%M')
    end = pd.Timestamp(data['timestamp'].iat[-1]).strftime('%m-%d-%Y %H:%M')
    
    # Merge and interpolate between start and end
    pdat = pd.DataFrame({'timestamp': pd.date_range(start=start, end=end, freq='min')})
    
    pdat = pdat.merge(data, on='timestamp', how='left')
    pdat['lat'] = pd.Series(pdat['lat']).interpolate()
    pdat['lon'] = pd.Series(pdat['lon']).interpolate()
    
    # Keep on the hour
    pdat = pdat[pd.Series(pdat['timestamp']).dt.minute == 00].reset_index(drop=True)
    pdat['mmsi'] = indat['mmsi'].iat[0]
    pdat = pdat.reset_index(drop=True)
    return pdat


def calc_dist(data_loc):
    indat = pd.read_feather(data_loc)
    indat = indat.sort_values('mmsi')
    
    # Puerto Mardryn Area
    lon1 = -65.925222
    lon2 = -54.418267
    lat1 = -45.145745
    lat2 = -38.147991
    rdat = indat[(indat['lon'] >= lon1) & (indat['lon'] <= lon2) & (indat['lat'] >= lat1) & (indat['lat'] <= lat2)]
    mmsi_keep = rdat['mmsi'].unique()
    outdat = indat[indat['mmsi'].isin(mmsi_keep)]

    # inter. data
    outdat = outdat.groupby('mmsi', as_index=False).apply(interp_hr)

    # Calc dist.
    odat = outdat.groupby('timestamp', as_index=False).apply(spherical_dist_populate)

    return odat



# Puerto Mardryn Area
lon1 = -65.925222
lon2 = -54.418267
lat1 = -45.145745
lat2 = -38.147991

files = glob.glob("/home/server/pi/homes/woodilla/Data/GFW_point/Patagonia_Shelf/feather" + "/*.feather")
    
nfiles = sorted(files)
nfiles = nfiles[68:79]
nfiles
test = calc_dist(nfiles[1])


list_ = []
for file in nfiles:
    
    df = calc_dist(file)
    list_.append(df)
    mdat = pd.concat(list_)

print(mdat)
mdat = mdat.reset_index(drop=True)
mdat.to_feather('~/Data/GFW_point/Patagonia_Shelf/complete/Puerto_Madryn_2016-03-10_2016-03-20.feather')

mdat = pd.read_feather('~/Data/GFW_point/Patagonia_Shelf/complete/Puerto_Madryn_2016-03-10_2016-03-20.feather')

# Patagonia shelf
lon1 = -77
lon2 = -22
lat1 = -58
lat2 = -23

#
mlon1 = -65.925222
mlon2 = -54.418267
mlat1 = -45.145745
mlat2 = -38.147991



my_dpi=96
fig = plt.figure(figsize=(2600/my_dpi, 1800/my_dpi), dpi=my_dpi)
ax = plt.axes(projection=ccrs.PlateCarree())
ax.set_extent([lon1, lon2, lat1, lat2], crs=ccrs.PlateCarree())
ax.coastlines()
#ax.stock_img()

ax.plot(mlon1, mlat1, 'o', markersize=10, color = 'red', transform=ccrs.PlateCarree())
ax.plot(mlon2, mlat2, 'o', markersize=10, color = 'red', transform=ccrs.PlateCarree())
ax.plot(mlon1, mlat2, 'o', markersize=10, color = 'red', transform=ccrs.PlateCarree())
ax.plot(mlon2, mlat1, 'o', markersize=10, color = 'red', transform=ccrs.PlateCarree())



# Specific time
my_dpi=96
fig = plt.figure(figsize=(2600/my_dpi, 1800/my_dpi), dpi=my_dpi)
ax = plt.axes(projection=ccrs.PlateCarree())
#ax.background_img(name='BM', resolution='low')
ax.set_extent([lon1, lon2, lat1, lat2], crs=ccrs.PlateCarree())
ax.coastlines()
#ax.stock_img()


#tdat = mdat[mdat['timestamp'] == "2016-03-10 00:00:00"]
for mmsi, track in test.groupby('vessel_A'):

    x = track.vessel_A_lon.values
    y = track.vessel_A_lat.values
    
    ax.plot(x, y, label=mmsi, transform=ccrs.PlateCarree())
    ax.plot(x, y, 'o', markersize=1, color = 'black', label=mmsi, transform=ccrs.PlateCarree())
ax.plot(-25.463418, -24.623938, 'o', markersize=10, color = 'red', transform=ccrs.PlateCarree())
ax.plot(pdat['lon'], pdat['lat'], label=mmsi, transform=ccrs.PlateCarree())
ax.plot(pdat['lon'], pdat['lat'], 'o', markersize=1, color = 'red', label=mmsi, transform=ccrs.PlateCarree())




# Animated Map

# Patagonia shelf
lon1 = -77
lon2 = -22
lat1 = -58
lat2 = -23

# Puerto Mardryn Area
lon1 = -65.925222
lon2 = -54.418267
lat1 = -45.145745
lat2 = -38.147991

# Zoom out Puerto Mardryn Area
lon1 = -68
lon2 = -51
lat1 = -48
lat2 = -35

#gdat = mdat[mdat['vessel_A'] == 900028621]
#vessel_b = list(gdat['vessel_B'].unique())
#vessel_b.append(gdat['vessel_A'].iat[0])

#gdat = mdat[mdat['vessel_B'].isin(vessel_b)]

gdat = mdat[mdat['timestamp'] == "2016-03-10 13:00:00"]
gdat = gdat[gdat['rank']  <= 5]
#gdat = mdat

for mmsi, group in gdat.groupby('timestamp'):
    #lon = group.vessel_A_lon.values[1]
    #lat = group.vessel_A_lat.values[1]
    #lon1 = lon - 1
    #lon2 = lon + 1
    #lat1 = lat - 1
    #lat2 = lat + 1
    
    # Zoom out Puerto Mardryn Area
    lon1 = -68
    lon2 = -51
    lat1 = -48
    lat2 = -35
    
    my_dpi=96
    
    #request = cimgt.GoogleTiles()
    #fig, ax = make_map(projection=request.crs)
    
    fig = plt.figure(figsize=(2600/my_dpi, 1800/my_dpi), dpi=my_dpi)
    ax = plt.axes(projection=cimgt.StamenTerrain().crs)
    #ax.stock_img()
    #ax.add_feature(cfeature.LAND)
    #ax.add_feature(cfeature.COASTLINE)
    # add some features to make the map a little more polished
    ax.set_extent([lon1, lon2, lat1, lat2], crs=cimgt.StamenTerrain().crs)
    ax.add_feature(cartopy.feature.LAND)
    ax.add_feature(cartopy.feature.OCEAN)
    ax.coastlines('50m')
    #ax.add_feature(states_provinces, edgecolor='gray')
    #cimgt.StamenTerrain().crs 
    #ax.coastlines()
    #ax.stock_img()

    VA_x = group.vessel_A_lon.values
    VA_y = group.vessel_A_lat.values
    
    VB_x = group.vessel_B_lon.values
    VB_y = group.vessel_B_lat.values
    
    # All ships
    #ax.plot(VB_x, VB_y, 'o', markersize=2, color = 'red', label=mmsi, transform=ccrs.PlateCarree())
    
    # Connected ships
    ax.plot([VB_x, VA_x], [VB_y,  VA_y], marker = 'o', linewidth=1, markersize=2, color = 'black', label=mmsi, transform=ccrs.PlateCarree())
    ax.plot(VA_x, VA_y, 'o', markersize=2, color = 'red', label=mmsi, transform=ccrs.PlateCarree())
    
    # Plot timestamp
    ax.annotate(f"Timestamp: {group['timestamp'].iat[0]}", xy=(0.22, .85), xycoords='figure fraction', fontsize=20, color='blue')
    #plt.annotate(f"MMSI: {group['vessel_A'].iat[0]}", xy=(0.7, .1), xycoords='figure fraction', fontsize=20, color='blue') 
    filepath = f"/home/server/pi/homes/woodilla/Projects/Patagonia-EDA/figures/{group['timestamp'].iat[0]}.jpg"
    
    #plt.savefig(filepath, dpi=300)

lon1 = -165
lon2 = 180
lat1 = -70
lat2 = 75
                 
import matplotlib.pyplot as plt
from matplotlib.transforms import offset_copy

import cartopy.crs as ccrs
import cartopy.io.img_tiles as cimgt

stamen_terrain = cimgt.StamenTerrain()

    # Create a GeoAxes in the tile's projection.
ax = plt.axes(projection=stamen_terrain.crs)

    # Limit the extent of the map to a small longitude/latitude range.
ax.set_extent([lon1, lon2, lat1, lat2])

    # Add the Stamen data at zoom level 8.
ax.add_image(stamen_terrain, 8)
                
ax.plot(gdat['vessel_A_lon'], gdat['vessel_A_lat'], 'o', markersize=2, color = 'red', label=mmsi, transform=ccrs.PlateCarree())                
#geodetic_transform = ccrs.Geodetic()._as_mpl_transform(ax)
#text_transform = offset_copy(geodetic_transform, units='dots', x=-25)

    # Add text 25 pixels to the left of the volcano.
#plt.text(-19.613333, 63.62, u'Eyjafjallajökull',
#             verticalalignment='center', horizontalalignment='right',
#             transform=text_transform,
#             bbox=dict(facecolor='sandybrown', alpha=0.5, boxstyle='round'))
plt.show()
