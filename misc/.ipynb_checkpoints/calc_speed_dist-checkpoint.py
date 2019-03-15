import pandas as pd
import os
import glob as glob
import multiprocessing

def GFW_directories():
    '''Get all GFW_point directions'''
    
    dirs = os.listdir(GFW_DIR)
    # Remove subfolders 'BK' and 'identities'
    if 'BK' in dirs:
        dirs.remove('BK')
    
    if 'identities' in dirs:
        dirs.remove('identities')
    
    return dirs

def process_gfw_speed(i):
    subdir = GFW_DIR + i
    allFiles = glob.glob(subdir + "/*.csv")
    list_ = []

    if len(allFiles) > 0:
    # Append files in subdir
        for file_ in allFiles:
            df = pd.read_csv(file_, index_col=None, header=0)
            list_.append(df)
            dat = pd.concat(list_, axis = 0, ignore_index = True)
    
    outdat = dat[['timestamp', 'speed']]
    outdat = outdat.dropna()
    outdat = outdat.reset_index(drop=True)
    outdat.to_feather('~/Data/GFW_point/speed_dist/' + i + '.feather')
    #outdat = outdat.append(speed, ignore_index=True)
    print(i)

if __name__ == '__main__':
    
    GFW_DIR = '/data2/GFW_point/'

    gfw_list_dirs = sorted(GFW_directories())

    gfw_list_dirs = gfw_list_dirs[0:366]
    #i = gfw_list_dirs[0]
    # Get subdirectory list of files
    
    pool = multiprocessing.Pool(5, maxtasksperchild=1)         
    #pool.start()
    pool.map(process_gfw_speed, gfw_list_dirs)
    pool.close()
    


allFiles = glob.glob('/home/server/pi/homes/woodilla/Data/GFW_point/speed_dist' + '/*.feather')
list_ = []

if len(allFiles) > 0:
    # Append files in subdir
    for file_ in allFiles:
    df = pd.read_csv(file_, index_col=None, header=0)
    list_.append(df)
    dat = pd.concat(list_, axis = 0, ignore_index = True)
    
    
    
    