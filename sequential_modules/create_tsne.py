# CS205
# Create a TSNE plot - Sequential

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import scprep
from sklearn.manifold import TSNE


def create_tsne(exp_data, clusters):
    # Perform PCA on the data
    sc_pca = scprep.reduce.pca(exp_data, n_components = 100, return_singular_values = False)
    
    # Perform TSNE on the data
    tsne_operator = TSNE(n_components = 2, perplexity = 30, random_state = 10)  # Perplexity usually set 5-50 (loosely number of neighbors)
    sc_tsne = pd.DataFrame(tsne_operator.fit_transform(sc_pca.iloc[:,0:100]))  # Perform tSNE on first 100 PC's
    sc_tsne.to_csv("TSNE_Data.csv")
    
    # Visualize the clusters on a TSNE plot
    fig, axes = plt.subplots()
    for clust in np.unique(clusters):
        # Pull the indices associated with the given cluster
        indices_to_keep = np.array(clusters) == clust
        
        sc_tsne_subset = sc_tsne.loc[indices_to_keep,:]
        
        axes.scatter(sc_tsne_subset.iloc[:,0], sc_tsne_subset.iloc[:,1], label=str(clust))
        
    axes.set_xlabel("TSNE-1", fontsize = 14)
    axes.set_ylabel("TSNE-2", fontsize = 14)
    axes.set_title("TSNE", fontsize = 16)
    axes.legend(bbox_to_anchor = (1.01, 1), loc = "upper left", title = "Cluster", markerscale = 1, ncol = 2)
    fig.tight_layout()
    plt.savefig("tsne_plot.png")

