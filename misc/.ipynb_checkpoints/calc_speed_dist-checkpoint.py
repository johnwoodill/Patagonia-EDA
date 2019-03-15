import pandas as pd
import os
import glob as glob
import numpy as np
import multiprocessing as mp

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
    
def reader(filename):
    print(filename)
    dat = pd.read_feather(filename)
    return dat['speed']

if __name__ == '__main__':
    
    #GFW_DIR = '/data2/GFW_point/'

    #gfw_list_dirs = sorted(GFW_directories())

    #gfw_list_dirs = gfw_list_dirs[0:366]
    #i = gfw_list_dirs[0]
    # Get subdirectory list of files
    
    #pool = multiprocessing.Pool(5, maxtasksperchild=1)         
    #pool.start()
    #pool.map(process_gfw_speed, gfw_list_dirs)
    #pool.close()
    
    allFiles = sorted(glob.glob('/home/server/pi/homes/woodilla/Data/GFW_point/speed_dist' + '/*.feather'))[0:5]
    
    with mp.Pool(processes=(mp.cpu_count() - 1)) as pool:
        chunks = pool.map(reader, allFiles)
        
    df = pd.concat(chunks)
    df2 = sorted(pd.Series(df).unique())
    print(df2[-1000:])
    #pool = Pool(5) # number of cores you want to use
    #df_list = pool.map(reader, allFiles) #creates a list of the loaded df's
    #df = pd.concat(df_list) # concatenates all the df's into a single df
    #df = np.concatenate(df_list, axis=0)
    #df = df.reset_index(drop=True)
    quant_ = np.quantile(df, [.1, .2, .3, .4, .5, .6, .7, .8, .9, .95, .96, .97, .98, .99])
    qdat = pd.DataFrame({'quant': [.1, .2, .3, .4, .5, .6, .7, .8, .9, .95, .96, .97, .98, .99], 'value': quant_})
    qdat = qdat.reset_index(drop=True)
    
    print(qdat)
    
    #qdat.to_feather('/home/server/pi/homes/woodilla/2016_vessel_speed.feather')

    