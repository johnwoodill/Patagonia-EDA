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


def f_js(x, y):
    return distance.jensenshannon(x, y)

def d_matrix(dat, interval, NN=1):
    dat = dat[dat['rank'] <= NN]
    dat['distance'] = np.log(1 + dat['distance'])
    
    dat.loc[:, 'day'] = pd.DatetimeIndex(dat['timestamp']).day
    dat.loc[:, 'hour'] = pd.DatetimeIndex(dat['timestamp']).hour


    
    if interval == 'day':
        dat = dat.groupby(['vessel_A', 'day'], as_index=False)['distance'].mean()
        x = []
        g = dat.groupby(['day'])['distance']
        for k1, g1 in g:
            for k2, g2 in g:
                x += [(k1, k2, f_js(g1, g2))]

        distMatrix = pd.DataFrame(x).pivot(index=0, columns=1, values=2)
        distArray = ssd.squareform(distMatrix)
    
    if interval == 'dayhour':
        dat = dat.groupby(['vessel_A', 'timestamp'], as_index=False)['distance'].mean()
        x = []
        g = dat.groupby(['timestamp'])['distance']
        for k1, g1 in g:
            for k2, g2 in g:
                x += [(k1, k2, f_js(g1, g2))]

        distMatrix = pd.DataFrame(x).pivot(index=0, columns=1, values=2)
        distArray = ssd.squareform(distMatrix)

    return (distMatrix, distArray)


# H-Clustering
def h_cluster(distArray):
    Z = linkage(distArray, "ward")
    dn = dendrogram(Z, labels=['10', '11', '12', '13', '14', '15', '16', '17', '18', '19', '20'])
    return dn

def k_medoids(distMatrix, interval, init_medoids):

    # K-Mediods Clustering
    distMatrix = np.array(distMatrix)

    # K-Medoids Clustering
    
    #initial_medoids = [30, 90, 140]

    initial_medoids = init_medoids
    
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
        
#     if interval == 'day':
#         df['value'] += 10

#         ax = sns.scatterplot(x='value', y='group', data=df)
#         ax.set(xlabel='March', ylabel='Cluster')
#         ax.set_yticks([0, 1])
#         ax.set_xticks([10, 11, 12, 13, 14, 15, 16, 17, 18, 19, 20])
#         plt.show()
#         print(medoids)
        
#     if interval == 'dayhour':
#         ax = sns.scatterplot(x='value', y='group', data=df)
#         ax.set(xlabel='March', ylabel='Cluster')
#         ax.set_yticks(range(len(init_medoids)))
#         plt.xticks([0, 24, 48, 72, 96, 120, 144, 168, 192, 216, 240], ['10', '11', '12', '13', '14', '15', '16', '17', '18', '19', '20'])
#         plt.axvline(120, color='red')
#         plt.show()
#         print(medoids)
    
    return(df)

#---------------------------------------------------
# Robustness Checks (Completed 3/28/2019)
# 1. Main region (1) time interval 5-25 and 1-31
#---------------------------------------------------

# Puerto Madryn March 5-25
dat = pd.read_feather('~/Data/GFW_point/Patagonia_Shelf/complete/Puerto_Madryn_region1_2016-3-5_2016-3-25.feather')
# Import data
print("Robustness Checks: Region 1 March 5 - 25")
# # Day
distMatrix, distArray = d_matrix(dat, interval='day', NN=1)
pdat1 = k_medoids(distMatrix, interval='day', init_medoids=[2, 5, 8])
# h_cluster(distArray)
distMatrix = distMatrix.reset_index(drop=False)
distMatrix.columns = distMatrix.columns.astype(str)
distMatrix.to_feather('~/Data/GFW_point/Patagonia_Shelf/Puerto_Madryn/dist_matrices/dmat_Puerto_Madryn_region1_NN1_day_2016_03-05_2016-03-25.feather')
pdat1.to_feather('~/Projects/Patagonia-EDA/data/Puerto_Madryn_region1_NN1_k_medoids_day_2016-03-05_2016-03-25.feather')

distMatrix, distArray = d_matrix(dat, interval='day', NN=5)
pdat2 = k_medoids(distMatrix, interval='day', init_medoids=[2, 5, 8])
distMatrix = distMatrix.reset_index(drop=False)
distMatrix.columns = distMatrix.columns.astype(str)
#h_cluster(distArray)
distMatrix.to_feather('~/Data/GFW_point/Patagonia_Shelf/Puerto_Madryn/dist_matrices/dmat_Puerto_Madryn_region1_NN5_day_2016_03-05_2016-03-25.feather')
pdat2.to_feather('~/Projects/Patagonia-EDA/data/Puerto_Madryn_region1_NN5_k_medoids_day_2016-03-05_2016-03-25.feather')

# # Day by hour
distMatrix_dh, distArray_dh = d_matrix(dat, interval='dayhour', NN=1) 
pdat3 = k_medoids(distMatrix_dh, interval='dayhour', init_medoids=[30, 90, 140])
distMatrix_dh = distMatrix_dh.reset_index(drop=False)
distMatrix_dh.columns = distMdistMatrix_dhatrix.columns.astype(str)
distMatrix_dh.to_feather('~/Data/GFW_point/Patagonia_Shelf/Puerto_Madryn/dist_matrices/dmat_Puerto_Madryn_region1_NN1_day-hour_2016-03-05_2016-03-25.feather')
pdat3.to_feather('~/Projects/Patagonia-EDA/data/Puerto_Madryn_region1_NN1_k_medoids_dayhour_2016-03-05_2016-03-25.feather')

distMatrix_dh, distArray_dh = d_matrix(dat, interval='dayhour', NN=5) 
pdat4 = k_medoids(distMatrix_dh, interval='dayhour', init_medoids=[30, 90, 140])
distMatrix_dh = distMatrix_dh.reset_index(drop=False)
distMatrix_dh.columns = distMatrix_dh.columns.astype(str)
distMatrix_dh.to_feather('~/Data/GFW_point/Patagonia_Shelf/Puerto_Madryn/dist_matrices/dmat_Puerto_Madryn_region1_NN5_day-hour_2016-03-05_2016-03-25.feather')
pdat4.to_feather('~/Projects/Patagonia-EDA/data/Puerto_Madryn_region1_NN5_k_medoids_dayhour_2016-03-05_2016-03-25.feather')


# Puerto Madryn March 1-31
print("Robustness Checks: Region 1 March 1 - 31")
# Import data
dat = pd.read_feather('~/Data/GFW_point/Patagonia_Shelf/complete/Puerto_Madryn_region1_2016-3-1_2016-3-31.feather')

# # Day
distMatrix, distArray = d_matrix(dat, interval='day', NN=1)
pdat1 = k_medoids(distMatrix, interval='day', init_medoids=[2, 5, 8])
# h_cluster(distArray)
distMatrix = distMatrix.reset_index(drop=False)
distMatrix.columns = distMatrix.columns.astype(str)
distMatrix.to_feather('~/Data/GFW_point/Patagonia_Shelf/Puerto_Madryn/dist_matrices/dmat_Puerto_Madryn_region1_NN1_day_2016-03-01_2016-03-31.feather')
pdat1.to_feather('~/Projects/Patagonia-EDA/data/Puerto_Madryn_region1_NN1_k_medoids_day_2016-03-01_2016-03-31.feather')

distMatrix, distArray = d_matrix(dat, interval='day', NN=5)
pdat2 = k_medoids(distMatrix, interval='day', init_medoids=[2, 5, 8])
distMatrix = distMatrix.reset_index(drop=False)
distMatrix.columns = distMatrix.columns.astype(str)
#h_cluster(distArray)
distMatrix.to_feather('~/Data/GFW_point/Patagonia_Shelf/Puerto_Madryn/dist_matrices/dmat_Puerto_Madryn_region1_NN5_day_2016-03-01_2016-03-31.feather')
pdat2.to_feather('~/Projects/Patagonia-EDA/data/Puerto_Madryn_region1_NN5_k_medoids_day_2016-03-01_2016-03-31.feather')

# # Day by hour
distMatrix_dh, distArray_dh = d_matrix(dat, interval='dayhour', NN=1) 
pdat3 = k_medoids(distMatrix_dh, interval='dayhour', init_medoids=[30, 90, 140])
***
distMatrix_dh = distMatrix_dh.reset_index(drop=False)
distMatrix_dh = distMatrix_dh.rename(columns={x:y for x,y in zip(distMatrix_dh.columns,range(0,len(distMatrix_dh.columns)))})
distMatrix_dh.columns = distMatrix_dh.columns.astype(str)
distMatrix_dh = distMatrix_dh.iloc[:, 1:]
distMatrix_dh.to_feather('~/Data/GFW_point/Patagonia_Shelf/Puerto_Madryn/dist_matrices/dmat_Puerto_Madryn_region1_NN1_day-hour_2016-03-01_2016-03-31.feather')
***
pdat3.to_feather('~/Projects/Patagonia-EDA/data/Puerto_Madryn_region1_NN1_k_medoids_dayhour_2016-03-01_2016-03-31.feather')

distMatrix_dh, distArray_dh = d_matrix(dat, interval='dayhour', NN=5) 
pdat4 = k_medoids(distMatrix_dh, interval='dayhour', init_medoids=[30, 90, 140])
distMatrix_dh = distMatrix_dh.reset_index(drop=False)
distMatrix_dh.columns = distMatrix_dh.columns.astype(str)
distMatrix_dh.to_feather('~/Data/GFW_point/Patagonia_Shelf/Puerto_Madryn/dist_matrices/dmat_Puerto_Madryn_region1_NN5_day-hour_2016-03-01_2016-03-31.feather')
pdat4.to_feather('~/Projects/Patagonia-EDA/data/Puerto_Madryn_region1_NN5_k_medoids_dayhour_2016-03-01_2016-03-31.feather')

#---------------------------------------------------
# Robustness Checks
# 1. Region 2 time interval 10-20, 5-25, and 1-31
#---------------------------------------------------

# # Puerto Madryn March 10-20
# # Import data
print("Robustness Checks: Region 2 March 10 - 20")
dat = pd.read_feather('~/Data/GFW_point/Patagonia_Shelf/complete/Puerto_Madryn_region2_2016-3-10_2016-3-20.feather')

# # Day
distMatrix, distArray = d_matrix(dat, interval='day', NN=1)
pdat1 = k_medoids(distMatrix, interval='day', init_medoids=[2, 5, 8])
# h_cluster(distArray)
distMatrix = distMatrix.reset_index(drop=False)
distMatrix.columns = distMatrix.columns.astype(str)
distMatrix.to_feather('~/Data/GFW_point/Patagonia_Shelf/Puerto_Madryn/dist_matrices/dmat_Puerto_Madryn_region2_NN1_day_2016-03-10_2016-03-20.feather')
pdat1.to_feather('~/Projects/Patagonia-EDA/data/Puerto_Madryn_region2_NN1_k_medoids_day_2016-03-10_2016-03-20.feather')

distMatrix, distArray = d_matrix(dat, interval='day', NN=5)
pdat2 = k_medoids(distMatrix, interval='day', init_medoids=[2, 5, 8])
distMatrix = distMatrix.reset_index(drop=False)
distMatrix.columns = distMatrix.columns.astype(str)
#h_cluster(distArray)
distMatrix.to_feather('~/Data/GFW_point/Patagonia_Shelf/Puerto_Madryn/dist_matrices/dmat_Puerto_Madryn_region2_NN5_day_2016-03-10_2016-03-20.feather')
pdat2.to_feather('~/Projects/Patagonia-EDA/data/Puerto_Madryn_region2_NN5_k_medoids_day_2016-03-10_2016-03-20.feather')

# # Day by hour
distMatrix_dh, distArray_dh = d_matrix(dat, interval='dayhour', NN=1) 
pdat3 = k_medoids(distMatrix_dh, interval='dayhour', init_medoids=[30, 90, 140])
distMatrix_dh = distMatrix_dh.reset_index(drop=False)
distMatrix_dh.columns = distMatrix_dh.columns.astype(str)
distMatrix_dh.to_feather('~/Data/GFW_point/Patagonia_Shelf/Puerto_Madryn/dist_matrices/dmat_Puerto_Madryn_region2_NN1_day-hour_2016-03-10_2016-03-20.feather')
pdat3.to_feather('~/Projects/Patagonia-EDA/data/Puerto_Madryn_region2_NN1_k_medoids_dayhour_2016-03-10_2016-03-20.feather')

distMatrix_dh, distArray_dh = d_matrix(dat, interval='dayhour', NN=5) 
pdat4 = k_medoids(distMatrix_dh, interval='dayhour', init_medoids=[30, 90, 140])
distMatrix_dh = distMatrix_dh.reset_index(drop=False)
distMatrix_dh.columns = distMatrix_dh.columns.astype(str)
distMatrix_dh.to_feather('~/Data/GFW_point/Patagonia_Shelf/Puerto_Madryn/dist_matrices/dmat_Puerto_Madryn_region2_NN5_day-hour_2016-03-10_2016-03-20.feather')
pdat4.to_feather('~/Projects/Patagonia-EDA/data/Puerto_Madryn_region2_NN5_k_medoids_dayhour_2016-03-10_2016-03-20.feather')

# Puerto Madryn March 5-25
print("Robustness Checks: Region 2 March 5 - 25")
# Import data
dat = pd.read_feather('~/Data/GFW_point/Patagonia_Shelf/complete/Puerto_Madryn_region2_2016-3-5_2016-3-25.feather')

# # Day
distMatrix, distArray = d_matrix(dat, interval='day', NN=1)
pdat1 = k_medoids(distMatrix, interval='day', init_medoids=[2, 5, 8])
distMatrix = distMatrix.reset_index(drop=False)
distMatrix.columns = distMatrix.columns.astype(str)
distMatrix.to_feather('~/Data/GFW_point/Patagonia_Shelf/Puerto_Madryn/dist_matrices/dmat_Puerto_Madryn_region2_NN1_day_2016-03-05_2016-03-25.feather')
pdat1.to_feather('~/Projects/Patagonia-EDA/data/Puerto_Madryn_region2_NN1_k_medoids_day_2016-03-05_2016-03-25.feather')

distMatrix, distArray = d_matrix(dat, interval='day', NN=5)
pdat2 = k_medoids(distMatrix, interval='day', init_medoids=[2, 5, 8])
distMatrix = distMatrix.reset_index(drop=False)
distMatrix.columns = distMatrix.columns.astype(str)
distMatrix.to_feather('~/Data/GFW_point/Patagonia_Shelf/Puerto_Madryn/dist_matrices/dmat_Puerto_Madryn_region2_NN5_day_2016-03-05_2016-03-25.feather')
pdat2.to_feather('~/Projects/Patagonia-EDA/data/Puerto_Madryn_region2_NN5_k_medoids_day_2016-03-05_2016-03-25.feather')

# # Day by hour
distMatrix_dh, distArray_dh = d_matrix(dat, interval='dayhour', NN=1) 
pdat3 = k_medoids(distMatrix_dh, interval='dayhour', init_medoids=[30, 90, 140])
distMatrix_dh = distMatrix_dh.reset_index(drop=False)
distMatrix_dh.columns = distMatrix_dh.columns.astype(str)
distMatrix_dh.to_feather('~/Data/GFW_point/Patagonia_Shelf/Puerto_Madryn/dist_matrices/dmat_Puerto_Madryn_region2_NN1_day-hour_2016-03-05_2016-03-25.feather')
pdat3.to_feather('~/Projects/Patagonia-EDA/data/Puerto_Madryn_region2_NN1_k_medoids_dayhour_2016-03-05_2016-03-25.feather')

distMatrix_dh, distArray_dh = d_matrix(dat, interval='dayhour', NN=5) 
pdat4 = k_medoids(distMatrix_dh, interval='dayhour', init_medoids=[30, 90, 140])
distMatrix_dh = distMatrix_dh.reset_index(drop=False)
distMatrix_dh.columns = distMatrix_dh.columns.astype(str)
distMatrix_dh.to_feather('~/Data/GFW_point/Patagonia_Shelf/Puerto_Madryn/dist_matrices/dmat_Puerto_Madryn_region2_NN5_day-hour_2016-03-05_2016-03-25.feather')
pdat4.to_feather('~/Projects/Patagonia-EDA/data/Puerto_Madryn_region2_NN5_k_medoids_dayhour_2016-03-05_2016-03-25.feather')



# Puerto Madryn March 1-31
print("Robustness Checks: Region 2 March 1 - 31")
# Import data
dat = pd.read_feather('~/Data/GFW_point/Patagonia_Shelf/complete/Puerto_Madryn_region2_2016-3-1_2016-3-31.feather')

# Day
distMatrix, distArray = d_matrix(dat, interval='day', NN=1)
pdat1 = k_medoids(distMatrix, interval='day', init_medoids=[2, 5, 8])
distMatrix = distMatrix.reset_index(drop=False)
distMatrix.columns = distMatrix.columns.astype(str)
distMatrix.to_feather('~/Data/GFW_point/Patagonia_Shelf/Puerto_Madryn/dist_matrices/dmat_Puerto_Madryn_region2_NN1_day_2016-03-01_2016-03-31.feather')
pdat1.to_feather('~/Projects/Patagonia-EDA/data/Puerto_Madryn_region2_NN1_k_medoids_day_2016-03-01_2016-03-31.feather')

distMatrix, distArray = d_matrix(dat, interval='day', NN=5)
pdat2 = k_medoids(distMatrix, interval='day', init_medoids=[2, 5, 8])
distMatrix = distMatrix.reset_index(drop=False)
distMatrix.columns = distMatrix.columns.astype(str)
distMatrix.to_feather('~/Data/GFW_point/Patagonia_Shelf/Puerto_Madryn/dist_matrices/dmat_Puerto_Madryn_region2_NN5_day_2016-03-01_2016-03-31.feather')
pdat2.to_feather('~/Projects/Patagonia-EDA/data/Puerto_Madryn_region2_NN5_k_medoids_day_2016-03-01_2016-03-31.feather')

# # Day by hour
distMatrix_dh, distArray_dh = d_matrix(dat, interval='dayhour', NN=1) 
pdat3 = k_medoids(distMatrix_dh, interval='dayhour', init_medoids=[30, 90, 140])
distMatrix_dh = distMatrix_dh.reset_index(drop=False)
distMatrix_dh.columns = distMatrix_dh.columns.astype(str)
distMatrix_dh.to_feather('~/Data/GFW_point/Patagonia_Shelf/Puerto_Madryn/dist_matrices/dmat_Puerto_Madryn_region2_NN1_day-hour_2016-03-01_2016-03-31.feather')
pdat3.to_feather('~/Projects/Patagonia-EDA/data/Puerto_Madryn_region2_NN1_k_medoids_dayhour_2016-03-01_2016-03-31.feather')

distMatrix_dh, distArray_dh = d_matrix(dat, interval='dayhour', NN=5) 
pdat4 = k_medoids(distMatrix_dh, interval='dayhour', init_medoids=[30, 90, 140])
distMatrix_dh = distMatrix_dh.reset_index(drop=False)
distMatrix_dh.columns = distMatrix_dh.columns.astype(str)
distMatrix_dh.to_feather('~/Data/GFW_point/Patagonia_Shelf/Puerto_Madryn/dist_matrices/dmat_Puerto_Madryn_region2_NN5_day-hour_2016-03-01_2016-03-31.feather')
pdat4.to_feather('~/Projects/Patagonia-EDA/data/Puerto_Madryn_region2_NN5_k_medoids_dayhour_2016-03-01_2016-03-31.feather')





#---------------------------------------------------
# Robustness Checks
# 1. Region 3 time interval 10-20, 5-25, and 1-31
#---------------------------------------------------

# Puerto Madryn March 10-20
print("Robustness Checks: Region 3 March 10 - 20")
# Import data
dat = pd.read_feather('~/Data/GFW_point/Patagonia_Shelf/complete/Puerto_Madryn_region3_2016-3-10_2016-3-20.feather')

# # Day
distMatrix, distArray = d_matrix(dat, interval='day', NN=1)
pdat1 = k_medoids(distMatrix, interval='day', init_medoids=[2, 5, 8])
distMatrix = distMatrix.reset_index(drop=False)
distMatrix.columns = distMatrix.columns.astype(str)
distMatrix.to_feather('~/Data/GFW_point/Patagonia_Shelf/Puerto_Madryn/dist_matrices/dmat_Puerto_Madryn_region3_NN1_day_2016-03-10_2016-03-20.feather')
pdat1.to_feather('~/Projects/Patagonia-EDA/data/Puerto_Madryn_region3_NN1_k_medoids_day_2016-03-10_2016-03-20.feather')

distMatrix, distArray = d_matrix(dat, interval='day', NN=5)
pdat2 = k_medoids(distMatrix, interval='day', init_medoids=[2, 5, 8])
distMatrix = distMatrix.reset_index(drop=False)
distMatrix.columns = distMatrix.columns.astype(str)
distMatrix.to_feather('~/Data/GFW_point/Patagonia_Shelf/Puerto_Madryn/dist_matrices/dmat_Puerto_Madryn_region3_NN5_day_2016-03-10_2016-03-20.feather')
pdat2.to_feather('~/Projects/Patagonia-EDA/data/Puerto_Madryn_region3_NN5_k_medoids_day_2016-03-10_2016-03-20.feather')

# # Day by hour
distMatrix_dh, distArray_dh = d_matrix(dat, interval='dayhour', NN=1) 
pdat3 = k_medoids(distMatrix_dh, interval='dayhour', init_medoids=[30, 90, 140])
distMatrix_dh = distMatrix_dh.reset_index(drop=False)
distMatrix_dh.columns = distMatrix_dh.columns.astype(str)
distMatrix_dh.to_feather('~/Data/GFW_point/Patagonia_Shelf/Puerto_Madryn/dist_matrices/dmat_Puerto_Madryn_region3_NN1_day-hour_2016-03-10_2016-03-20.feather')
pdat3.to_feather('~/Projects/Patagonia-EDA/data/Puerto_Madryn_region3_NN1_k_medoids_dayhour_2016-03-10_2016-03-20.feather')

distMatrix_dh, distArray_dh = d_matrix(dat, interval='dayhour', NN=5) 
pdat4 = k_medoids(distMatrix_dh, interval='dayhour', init_medoids=[30, 90, 140])
distMatrix_dh = distMatrix_dh.reset_index(drop=False)
distMatrix_dh.columns = distMatrix_dh.columns.astype(str)
distMatrix_dh.to_feather('~/Data/GFW_point/Patagonia_Shelf/Puerto_Madryn/dist_matrices/dmat_Puerto_Madryn_region3_NN5_day-hour_2016-03-10_2016-03-20.feather')
pdat4.to_feather('~/Projects/Patagonia-EDA/data/Puerto_Madryn_region3_NN5_k_medoids_dayhour_2016-03-10_2016-03-20.feather')


# Puerto Madryn March 5-25
print("Robustness Checks: Region 3 March 5 - 25")
# Import data
dat = pd.read_feather('~/Data/GFW_point/Patagonia_Shelf/complete/Puerto_Madryn_region3_2016-3-5_2016-3-25.feather')

# # Day
distMatrix, distArray = d_matrix(dat, interval='day', NN=1)
pdat1 = k_medoids(distMatrix, interval='day', init_medoids=[2, 5, 8])
distMatrix = distMatrix.reset_index(drop=False)
distMatrix.columns = distMatrix.columns.astype(str)
distMatrix.to_feather('~/Data/GFW_point/Patagonia_Shelf/Puerto_Madryn/dist_matrices/dmat_Puerto_Madryn_region3_NN1_day_2016-03-05_2016-03-25.feather')
pdat1.to_feather('~/Projects/Patagonia-EDA/data/Puerto_Madryn_region3_NN1_k_medoids_day_2016-03-05_2016-03-25.feather')

distMatrix, distArray = d_matrix(dat, interval='day', NN=5)
pdat2 = k_medoids(distMatrix, interval='day', init_medoids=[2, 5, 8])
distMatrix = distMatrix.reset_index(drop=False)
distMatrix.columns = distMatrix.columns.astype(str)
distMatrix.to_feather('~/Data/GFW_point/Patagonia_Shelf/Puerto_Madryn/dist_matrices/dmat_Puerto_Madryn_region3_NN5_day_2016-03-05_2016-03-25.feather')
pdat2.to_feather('~/Projects/Patagonia-EDA/data/Puerto_Madryn_region3_NN5_k_medoids_day_2016-03-05_2016-03-25.feather')

# # Day by hour
distMatrix_dh, distArray_dh = d_matrix(dat, interval='dayhour', NN=1) 
pdat3 = k_medoids(distMatrix_dh, interval='dayhour', init_medoids=[30, 90, 140])
distMatrix_dh = distMatrix_dh.reset_index(drop=False)
distMatrix_dh.columns = distMatrix_dh.columns.astype(str)
distMatrix_dh.to_feather('~/Data/GFW_point/Patagonia_Shelf/Puerto_Madryn/dist_matrices/dmat_Puerto_Madryn_region3_NN1_day-hour_2016-03-05_2016-03-25.feather')
pdat3.to_feather('~/Projects/Patagonia-EDA/data/Puerto_Madryn_region3_NN1_k_medoids_dayhour_2016-03-05_2016-03-25.feather')

distMatrix_dh, distArray_dh = d_matrix(dat, interval='dayhour', NN=5) 
pdat4 = k_medoids(distMatrix_dh, interval='dayhour', init_medoids=[30, 90, 140])
distMatrix_dh = distMatrix_dh.reset_index(drop=False)
distMatrix_dh.columns = distMatrix_dh.columns.astype(str)
distMatrix_dh.to_feather('~/Data/GFW_point/Patagonia_Shelf/Puerto_Madryn/dist_matrices/dmat_Puerto_Madryn_region3_NN5_day-hour_2016-03-05_2016-03-25.feather')
pdat4.to_feather('~/Projects/Patagonia-EDA/data/Puerto_Madryn_region3_NN5_k_medoids_dayhour_2016-03-05_2016-03-25.feather')




# Puerto Madryn March 1-31
print("Robustness Checks: Region 3 March 1 - 31")
# Import data
dat = pd.read_feather('~/Data/GFW_point/Patagonia_Shelf/complete/Puerto_Madryn_region3_2016-3-1_2016-3-31.feather')

# Day
distMatrix, distArray = d_matrix(dat, interval='day', NN=1)
pdat1 = k_medoids(distMatrix, interval='day', init_medoids=[2, 5, 8])
distMatrix = distMatrix.reset_index(drop=False)
distMatrix.columns = distMatrix.columns.astype(str)
distMatrix.to_feather('~/Data/GFW_point/Patagonia_Shelf/Puerto_Madryn/dist_matrices/dmat_Puerto_Madryn_region3_NN1_day_2016-03-01_2016-03-31.feather')
pdat1.to_feather('~/Projects/Patagonia-EDA/data/Puerto_Madryn_region3_NN1_k_medoids_day_2016-03-01_2016-03-31.feather')

distMatrix, distArray = d_matrix(dat, interval='day', NN=5)
pdat2 = k_medoids(distMatrix, interval='day', init_medoids=[2, 5, 8])
distMatrix = distMatrix.reset_index(drop=False)
distMatrix.columns = distMatrix.columns.astype(str)
distMatrix.to_feather('~/Data/GFW_point/Patagonia_Shelf/Puerto_Madryn/dist_matrices/dmat_Puerto_Madryn_region3_NN5_day_2016-03-01_2016-03-31.feather')
pdat2.to_feather('~/Projects/Patagonia-EDA/data/Puerto_Madryn_region3_NN5_k_medoids_day_2016-03-01_2016-03-31.feather')

# # Day by hour
distMatrix_dh, distArray_dh = d_matrix(dat, interval='dayhour', NN=1) 
pdat3 = k_medoids(distMatrix_dh, interval='dayhour', init_medoids=[30, 90, 140])
distMatrix_dh = distMatrix_dh.reset_index(drop=False)
distMatrix_dh.columns = distMatrix_dh.columns.astype(str)
distMatrix_dh.to_feather('~/Data/GFW_point/Patagonia_Shelf/Puerto_Madryn/dist_matrices/dmat_Puerto_Madryn_region3_NN1_day-hour_2016-03-01_2016-03-31.feather')
pdat3.to_feather('~/Projects/Patagonia-EDA/data/Puerto_Madryn_region3_NN1_k_medoids_dayhour_2016-03-01_2016-03-31.feather')

distMatrix_dh, distArray_dh = d_matrix(dat, interval='dayhour', NN=5) 
pdat4 = k_medoids(distMatrix_dh, interval='dayhour', init_medoids=[30, 90, 140])
distMatrix_dh = distMatrix_dh.reset_index(drop=False)
distMatrix_dh.columns = distMatrix_dh.columns.astype(str)
distMatrix_dh.to_feather('~/Data/GFW_point/Patagonia_Shelf/Puerto_Madryn/dist_matrices/dmat_Puerto_Madryn_region3_NN5_day-hour_2016-03-01_2016-03-31.feather')
pdat4.to_feather('~/Projects/Patagonia-EDA/data/Puerto_Madryn_region3_NN5_k_medoids_dayhour_2016-03-01_2016-03-31.feather')
