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

def d_matrix(dat, interval, NN=1):
    dat = dat[dat['rank'] <= NN]
    dat['distance'] = np.log(1 + dat['distance'])

    
    # Calc PDF for each day (don't need)
    dat.loc[:, 'day'] = pd.DatetimeIndex(dat['timestamp']).day
    dat.loc[:, 'hour'] = pd.DatetimeIndex(dat['timestamp']).hour

    
    min_day = min(dat['day'])
    max_day = max(dat['day'])
    min_hour = min(dat['hour'])
    max_hour = max(dat['hour'])
    

    #nobs_day = dat.groupby('day').agg('count')
    #nobs_day = min(nobs_day['timestamp'])
    
    #nobs_hour = dat.groupby(['day', 'hour']).agg('count')
    #nobs_hour = min(nobs_hour['timestamp'])
    
    if interval == 'day':
        # Apply JS-dist to each permutation of days
        jds_dmat = pd.DataFrame()
        for i in range(min_day, max_day + 1):
            for j in range(min_day, max_day + 1):
                x = dat[dat['day'] == i].distance
                y = dat[dat['day'] == j].distance
                jds = distance.jensenshannon(x, y)
                outdat = pd.DataFrame({'day_a': [i], 'day_b': [j], 'jds': [round(jds, 4)]})
                jds_dmat = jds_dmat.append(outdat, ignore_index=True)
        # Generate distance matrix
        distMatrix = jds_dmat.pivot(index='day_a', columns='day_b', values='jds')
        distArray = ssd.squareform(distMatrix)
    
    indat = dat
    if interval == 'dayhour':
        #vessels = dat['vessel_A'].unique()
        #vessels = np.random.choice(vessels, 45)
        #indat = dat[dat['vessel_A'].isin(vessels)]
        # Apply JS-dist to each permutation of days
        jds_dmat = pd.DataFrame()
        for k in range(min_day, max_day + 1):
            for l in range(min_day, max_day + 1):
                for i in range(min_hour, max_hour + 1):
                    for j in range(min_hour, max_hour + 1):
                        x = indat[(indat['hour'] == i) & (indat['day'] == k)].distance
                        y = indat[(indat['hour'] == j) & (indat['day'] == l)].distance
                        jds = distance.jensenshannon(x, y)
                        
# outmat = [distance.jensenshannon(x, y) 
#     for x in indat[(['hour'] == i) & (indat['day'] == k)].distance 
#     for i in range(min_hour, max_hour + 1) 
#     for k in range(min_day, max_day + 1) 
#     for y in indat[(indat['hour'] == j) &(indat['day'] == l)].distance 
#     for j in range(min_hour, max_hour + 1) 
#     for l in range(min_day, max_day + 1)
          
                        outdat = pd.DataFrame({'day_hour_a': f"{k}_{i}", 'day_hour_b': f"{l}_{j}", 'jds': [round(jds, 4)]})
                        jds_dmat = jds_dmat.append(outdat, ignore_index=True)
                    # Generate distance matrix
        distMatrix = jds_dmat.pivot(index='day_hour_a', columns='day_hour_b', values='jds')
        distArray = ssd.squareform(distMatrix)
    
    return (distMatrix, distArray)

#jds_dmat    
#jds_dmat.to_feather('~/js_dmat.feather')
#jds_dmat = pd.read_feather('~/js_dmat.feather')

# Cluster across x, y




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

# Puerto Madryn March 10-20
# Import data
#dat = pd.read_feather('~/Data/GFW_point/Patagonia_Shelf/complete/Puerto_Madryn_2016-03-10_2016-03-20.feather')

# # Day
#distMatrix, distArray = d_matrix(dat, interval='day', NN=1)
#pdat1 = k_medoids(distMatrix, interval='day', init_medoids=[2, 5, 8])
# h_cluster(distArray)
# pdat1.to_feather('~/Projects/Patagonia-EDA/data/k_medoids_day_NN1_3-10-2016-3-20-2016.feather')

# distMatrix, distArray = d_matrix(dat, interval='day', NN=5)
# pdat2 = k_medoids(distMatrix, interval='day', init_medoids=[2, 5, 8])
# h_cluster(distArray)
# pdat2.to_feather('~/Projects/Patagonia-EDA/data/k_medoids_day_NN5_3-10-2016-3-20-2016.feather')

# # Day by hour
# distMatrix_dh, distArray_dh = d_matrix(dat, interval='dayhour', NN=1) 
# pdat3 = k_medoids(distMatrix_dh, interval='dayhour', init_medoids=[30, 90, 140])
# pdat3.to_feather('~/Projects/Patagonia-EDA/data/k_medoids_day_hour_NN1_3-10-2016-3-20-2016.feather')

# distMatrix_dh, distArray_dh = d_matrix(dat, interval='dayhour', NN=5) 
# pdat4 = k_medoids(distMatrix_dh, interval='dayhour', init_medoids=[30, 90, 140])
# pdat4.to_feather('~/Projects/Patagonia-EDA/data/k_medoids_day_hour_NN5_3-10-2016-3-20-2016.feather')


# Puerto Madryn March 5-25
# Import data
print("Loading 3-5-2016 Data")
dat = pd.read_feather('~/Data/GFW_point/Patagonia_Shelf/complete/Puerto_Madryn_2016-3-5_2016-3-25.feather')

# # Day
# distMatrix, distArray = d_matrix(dat, interval='day', NN=1)
# pdat5 = k_medoids(distMatrix, interval='day', init_medoids=[5, 12, 17])
# pdat5.to_feather('~/Projects/Patagonia-EDA/data/k_medoids_day_NN1_3-05-2016-3-25-2016.feather')

# # Day
# distMatrix, distArray = d_matrix(dat, interval='day', NN=5)
# pdat6 = k_medoids(distMatrix, interval='day', init_medoids=[5, 12, 17])
# pdat6.to_feather('~/Projects/Patagonia-EDA/data/k_medoids_day_NN5_3-05-2016-3-25-2016.feather')

# # Day by hour
# distMatrix_dh, distArray_dh = d_matrix(dat, interval='dayhour', NN=1) 
# pdat7 = k_medoids(distMatrix_dh, interval='dayhour', init_medoids=[80, 264, 400])
# pdat7.to_feather('~/Projects/Patagonia-EDA/data/k_medoids_day_hour_NN1_3-05-2016-3-25-2016.feather')

print('Processing dayhour 3-5-2016 NN=5')
distMatrix_dh, distArray_dh = d_matrix(dat, interval='dayhour', NN=5) 
pdat8 = k_medoids(distMatrix_dh, interval='dayhour', init_medoids=[80, 264, 400])
pdat8.to_feather('~/Projects/Patagonia-EDA/data/k_medoids_day_hour_NN5_3-05-2016-3-25-2016.feather')



# Puerto Madryn March 1-31
# Import data
dat = pd.read_feather('~/Data/GFW_point/Patagonia_Shelf/complete/Puerto_Madryn_2016-3-1_2016-3-31.feather')

# Day
print('Processing day 3-1-2016 NN=1')
distMatrix, distArray = d_matrix(dat, interval='day', NN=1)
pdat5 = k_medoids(distMatrix, interval='day', init_medoids=[8, 15, 23])
pdat5.to_feather('~/Projects/Patagonia-EDA/data/k_medoids_day_NN1_3-01-2016-3-31-2016.feather')

# Day
print('Processing day 3-1-2016 NN=5')
distMatrix, distArray = d_matrix(dat, interval='day', NN=5)
pdat6 = k_medoids(distMatrix, interval='day', init_medoids=[8, 15, 23])
pdat6.to_feather('~/Projects/Patagonia-EDA/data/k_medoids_day_NN5_3-01-2016-3-31-2016.feather')

# Day by hour
print('Processing dayhour 3-1-2016 NN=1')
distMatrix_dh, distArray_dh = d_matrix(dat, interval='dayhour', NN=1) 
pdat7 = k_medoids(distMatrix_dh, interval='dayhour', init_medoids=[192, 360, 504])
pdat7.to_feather('~/Projects/Patagonia-EDA/data/k_medoids_day_hour_NN1_3-01-2016-3-31-2016.feather')

print('Processing dayhour 3-1-2016 NN=5')
distMatrix_dh, distArray_dh = d_matrix(dat, interval='dayhour', NN=5) 
pdat8 = k_medoids(distMatrix_dh, interval='dayhour', init_medoids=[192, 360, 504])
pdat8.to_feather('~/Projects/Patagonia-EDA/data/k_medoids_day_hour_NN5_3-01-2016-3-31-2016.feather')