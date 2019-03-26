import pandas as pd
from scipy import stats
import numpy as np
from dit.divergences import jensen_shannon_divergence
import matplotlib.pyplot as plt
from scipy.spatial import distance
import seaborn as sns
import scipy.spatial.distance as ssd
import random
from pyclustering.cluster.kmedoids import kmedoids
from pyclustering.cluster import cluster_visualizer, cluster_visualizer_multidim
from pyclustering.cluster.kmedoids import kmedoids
from pyclustering.utils import read_sample
from pyclustering.utils import timedcall
from scipy.cluster.hierarchy import dendrogram, linkage
import scipy.spatial.distance as ssd

def d_matrix(dat, NN=1):
    dat = dat[dat['rank'] <= NN]
    dat['distance'] = np.log(1 + dat['distance'])


    # Calc PDF for each day (don't need)
    dat.loc[:, 'day'] = pd.DatetimeIndex(dat['timestamp']).day
    dat = dat[dat['day'] != 21]

    nobs = dat.groupby('day').agg('count')
    nobs = min(nobs['timestamp'])
    # Apply JS-dist to each permutation of days
    jds_dmat = pd.DataFrame()
    for i in range(10, 20):
        for j in range(10, 20):
            x = dat[dat['day'] == i][0:nobs].distance
            y = dat[dat['day'] == j][0:nobs].distance
            jds = distance.jensenshannon(x, y)
            outdat = pd.DataFrame({'day_a': [i], 'day_b': [j], 'jds': [round(jds, 4)]})
            jds_dmat = jds_dmat.append(outdat, ignore_index=True)
    # Generate distance matrix
    distMatrix = jds_dmat.pivot(index='day_a', columns='day_b', values='jds')
    distArray = ssd.squareform(distMatrix)
    return (distMatrix, distArray)

#jds_dmat    
#jds_dmat.to_feather('~/js_dmat.feather')
#jds_dmat = pd.read_feather('~/js_dmat.feather')

# Cluster across x, y




# H-Clustering
def h_cluster(distArray):
    Z = linkage(distArray, "ward")
    dn = dendrogram(Z)
    return dn

def k_medoids(distMatrix):

    # K-Mediods Clustering
    distMatrix = np.array(distMatrix)

    # K-Medoids Clustering
    initial_medoids = [2, 5]

    # create K-Medoids algorithm for processing distance matrix instead of points
    kmedoids_instance = kmedoids(distMatrix, initial_medoids, data_type='distance_matrix', ccore=True)

    # run cluster analysis and obtain results
    kmedoids_instance.process()

    clusters = kmedoids_instance.get_clusters()
    medoids = kmedoids_instance.get_medoids()
    print(f"Clusters: {clusters}   Medoids: {medoids}")
    

    final_list = []
    for i, l in enumerate(clusters):
        for num in l:
            final_list.append({'value': num, 'group': i})

    df = pd.DataFrame(final_list)
    df['value'] += 10

    ax = sns.scatterplot(x='value', y='group', data=df)
    ax.set(xlabel='March', ylabel='Cluster')
    ax.set_yticks([0, 1])
    ax.set_xticks([10, 11, 12, 13, 14, 15, 16, 17, 18, 19, 20])
    plt.show()

    
# Subset out
plon1 = -62.925222
plon2 = -54.418267
plat1 = -45
plat2 = -40
    
# Import data
dat = pd.read_feather('~/Data/GFW_point/Patagonia_Shelf/complete/Puerto_Madryn_2016-03-10_2016-03-20.feather')

dat = dat[(dat['vessel_A_lon'] >= plon1) & (dat['vessel_A_lon'] <= plon2)]
dat = dat[(dat['vessel_A_lat'] >= plat1) & (dat['vessel_A_lat'] <= plat2)]
dat = dat.reset_index(drop=True)
dat.to_feather('~/test.feather')

distMatrix, distArray = d_matrix(dat, NN=5)            
distMatrix


h_cluster(distArray)
k_medoids(distMatrix)


