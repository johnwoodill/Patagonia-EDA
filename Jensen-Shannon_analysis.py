import pandas as pd
from scipy import stats
import numpy as np
from dit.divergences import jensen_shannon_divergence
from scipy.spatial import distance
import scipy.spatial.distance as ssd
import random

# Jensen-Shannon divergences calculation
#https://stats.stackexchange.com/questions/29578/jensen-shannon-divergence-calculation-for-3-prob-distributions-is-this-ok

# @author: jonathanfriedman

def jsd(x,y): #Jensen-shannon divergence
    import warnings
    warnings.filterwarnings("ignore", category = RuntimeWarning)
    x = np.array(x)
    y = np.array(y)
    d1 = x*np.log2(2*x/(x+y))
    d2 = y*np.log2(2*y/(x+y))
    d1[np.isnan(d1)] = 0
    d2[np.isnan(d2)] = 0
    d = 0.5*np.sum(d1+d2)    
    return d

jsd(np.array([0.5,0.5,0]),np.array([0,0.1,0.9]))

def calc_pdf(data, center=False):
    values = data['distance'].sample(20000)
    in_pdf = stats.kde.gaussian_kde(values.ravel())
    out_list = list(in_pdf(values))
    
    if center == True:
        out_list /= np.sum(out_list)
    return out_list

# Import data
dat = pd.read_feather('~/Data/GFW_point/Patagonia_Shelf/complete/Puerto_Madryn_2016-03-10_2016-03-20.feather')
dat['distance'] = np.log(1 + dat['distance'])

# Calc PDF for each day (don't need)
dat.loc[:, 'day'] = pd.DatetimeIndex(dat['timestamp']).day
#pdfd = dat.groupby('day').apply(calc_pdf, center=True)

# Apply JS-dist to each permutation of days
jds_dmat = pd.DataFrame()
for i in range(10, 20):
    for j in range(10, 20):
        x = dat[dat['day'] == i][0:20000].distance
        y = dat[dat['day'] == j][0:20000].distance
        jds = distance.jensenshannon(x, y)
        outdat = pd.DataFrame({'day_a': [i], 'day_b': [j], 'jds': [round(jds, 4)]})
        jds_dmat = jds_dmat.append(outdat, ignore_index=True)
        
jds_dmat    
jds_dmat.to_feather('~/js_dmat.feather')


# Cluster across x, y


# Generate distance matrix
import scipy.spatial.distance as ssd
distMatrix = jds_dmat.pivot(index='day_a', columns='day_b', values='jds')
distArray = ssd.squareform(distMatrix)
distArray

from scipy.cluster.hierarchy import dendrogram, linkage

Z = linkage(distArray)
dn = dendrogram(Z)

distMatrix.to_dict()





