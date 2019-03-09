import pandas as pd
import glob
import multiprocessing

# Get list of files in folder
files = []
for file in glob.glob("/home/server/pi/homes/woodilla/Data/GFW_point/Patagonia_Shelf/distance_data/circle_measure" + "/*.feather"):
    files.append(file)
nfiles = sorted(files)

# Get 10 closests NN
def gfw_process_nn(data_loc):
    indat = pd.read_feather(data_loc)
    outdat = indat.sort_values('distance').groupby('vessel_A').nth([1,2, 3, 4, 5, 6, 7, 8, 9, 10])
    outdat = outdat.sort_values(['vessel_A', 'distance'])
    outdat['rank'] = outdat.groupby(['vessel_A']).cumcount() + 1
    outdat = outdat.reset_index()
    outdat.to_feather('/home/server/pi/homes/woodilla/Data/GFW_point/Patagonia_Shelf/distance_data/NearestNeighbor/' + outdat['date'][1] + '_NN.feather')
    return 0

pool = multiprocessing.Pool(5, maxtasksperchild=1)         
pool.map(gfw_process_nn, nfiles)
pool.close()

# n processed files
nprocess = glob.glob('/home/server/pi/homes/woodilla/Data/GFW_point/Patagonia_Shelf/distance_data/NearestNeighbor/' + "*.feather")

# If all files have been processed
if len(nprocess) == len(nfiles):
    # Get subdirectory list of files
    subdir = '/home/server/pi/homes/woodilla/Data/GFW_point/Patagonia_Shelf/distance_data/NearestNeighbor/'
    allFiles = glob.glob(subdir + "*.feather")
    list_ = []

    # Append files in subdir
    for file_ in allFiles:
        df = pd.read_feather(file_)
        list_.append(df)
        dat = pd.concat(list_, axis = 0, ignore_index = True)  

# Save file
dat = dat.sort_values(['date', 'vessel_A', 'distance'])
dat = dat[['date', 'vessel_A', 'vessel_B', 'distance', 'rank']]
dat = dat.reset_index(drop=True)
dat.to_feather('~/Data/GFW_point/Patagonia_Shelf/complete/NN_2016-2018.feather')
    
    