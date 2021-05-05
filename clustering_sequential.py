# CS205
# Perform clustering on the single-cell data - Sequential

import numpy as np
import pandas as pd
import scanpy as sc
from anndata import AnnData
from sklearn.cluster import KMeans


# K-means clustering
def kmeans_clustering(exp_data, n_clusters):
    clusters = KMeans(n_clusters = n_clusters).fit_predict(exp_data)
    pd.DataFrame(clusters).to_csv("kmeans_clusters.csv")
    return clusters

# Louvain clustering
def louvain_clustering(exp_data):
    # Create a KNN graph from the data (community detection is performed on the graph)
    knn_graph = sc.pp.neighbors(AnnData(exp_data), n_neighbors = 20, n_pcs = 50, copy = True)
    
    # Perform Louvain clustering using scanpy
    louvain_data = sc.tl.louvain(knn_graph, copy = True)
    louvain_communities = pd.DataFrame(louvain_data.obs["louvain"])  # Pull the cluster assignments
    louvain_communities.to_csv("louvain_clusters.csv")
    
    # Convert the clusters to a list
    louvain_communites = louvain_communities["louvain"]
    louvain_communities = list(map(int, louvain_communities.to_numpy()))
    
    return louvain_communities
