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

    dat = dat[dat['day'] != 21]
    
    min_day = min(dat['day'])
    max_day = max(dat['day'])
    min_hour = min(dat['hour'])
    max_hour = max(dat['hour'])
    

    nobs_day = dat.groupby('day').agg('count')
    nobs_day = min(nobs_day['timestamp'])
    
    nobs_hour = dat.groupby(['day', 'hour']).agg('count')
    nobs_hour = min(nobs_hour['timestamp'])
    
    if interval == 'day':
        # Apply JS-dist to each permutation of days
        jds_dmat = pd.DataFrame()
        for i in range(min_day, max_day):
            for j in range(min_day, max_day):
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
        for k in range(min_day, max_day):
            for l in range(min_day, max_day):
                for i in range(min_hour, max_hour):
                    for j in range(min_hour, max_hour):
                        x = indat[(indat['hour'] == i) & (indat['day'] == k)].distance
                        y = indat[(indat['hour'] == j) & (indat['day'] == l)].distance
                        jds = distance.jensenshannon(x, y)
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
    
    if interval == 'day':
        df['value'] += 10

        ax = sns.scatterplot(x='value', y='group', data=df)
        ax.set(xlabel='March', ylabel='Cluster')
        ax.set_yticks([0, 1])
        ax.set_xticks([10, 11, 12, 13, 14, 15, 16, 17, 18, 19, 20])
        plt.show()
        print(medoids)
        
    if interval == 'dayhour':
        ax = sns.scatterplot(x='value', y='group', data=df)
        ax.set(xlabel='March', ylabel='Cluster')
        ax.set_yticks(range(len(init_medoids)))
        plt.xticks([0, 24, 48, 72, 96, 120, 144, 168, 192, 216, 240], ['10', '11', '12', '13', '14', '15', '16', '17', '18', '19', '20'])
        plt.axvline(120, color='red')
        plt.show()
        print(medoids)
        

    
# Import data
dat = pd.read_feather('~/Data/GFW_point/Patagonia_Shelf/complete/Puerto_Madryn_2016-03-10_2016-03-20.feather')

# Day
distMatrix, distArray = d_matrix(dat, interval='day', NN=5)
k_medoids(distMatrix, interval='day', init_medoids=[2, 5, 8])
h_cluster(distArray)

# Day by hour
distMatrix, distArray = d_matrix(dat, interval='dayhour', NN=3) 
k_medoids(distMatrix, interval='dayhour', init_medoids=[30, 90, 140])
h_cluster(distArray)


