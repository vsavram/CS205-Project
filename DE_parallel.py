# CS205
# Perform differential expression analysis - Sequential
from __future__ import division
import numpy as np
import pandas as pd
import math
from scipy.stats import ranksums
import os
from mpi4py import MPI

def DE():

    comm = MPI.COMM_WORLD
    size = comm.Get_size()
    rank = comm.Get_rank()

    if rank == 0:
        exp_data = pd.read_csv('exp_data.csv', index_col='Unnamed: 0')
        metadata = pd.read_csv('metadata.csv', index_col='Unnamed: 0')
        clusters = pd.read_csv('clusters.csv', index_col='Unnamed: 0')

        # Cluster array
        clusters = np.array(clusters).flatten()
        iters =int(math.ceil(len(np.unique(clusters)) / (size - 1)))
        queue = []
        for x in range(iters):
            queue.extend(range(1, size))
        # Making first contact with worker nodes
        for i in range(size-1):
            comm.send(iters, dest=i + 1, tag=1000)

        # Iterate over every cluster
        for i, clust in enumerate(np.unique(clusters)):
            # Filter the expression data and metadata to only include cells from the given cluster and send
            # it to node
            cluster_exp = exp_data[clusters == clust]
            cluster_metadata = metadata[clusters == clust]
            dest = queue.pop()
            comm.send(cluster_exp, dest=dest, tag=0)
            comm.send(cluster_metadata, dest=dest, tag=1)
            comm.send(exp_data.columns, dest=dest, tag=2)
            comm.send(clust, dest=dest, tag=3)

        while queue:
            dest = queue.pop()
            comm.send(None, dest=dest, tag=0)
            comm.send(None, dest=dest, tag=1)
            comm.send(None, dest=dest, tag=2)
            comm.send(None, dest=dest, tag=3)
    else:
        iters = comm.recv(source=0, tag=1000)
        count = 0
        while count < iters:
            # Worker node recieves cluster data
            cluster_exp = comm.recv(source=0, tag=0)
            cluster_metadata = comm.recv(source=0, tag=1)
            exp_data_columns = comm.recv(source=0, tag=2)
            clust = comm.recv(source=0, tag=3)
            if not isinstance(cluster_exp, pd.DataFrame):
                count += 1

            else:
                # Create the directory used to store the differential expression results if it does not exist
                if not os.path.exists('differential_expression'):
                    os.makedirs('differential_expression')

                # Determine the data associated with the COVID group and healthy group
                covid_group = cluster_exp[cluster_metadata['Sample Characteristic[disease]'].values == 'COVID-19']
                healthy_group = cluster_exp[cluster_metadata['Sample Characteristic[disease]'].values == 'normal']
                    
                    # Create lists used to store the log2 fold changes and p-values associated with each gene
                log_FC_list, p_val_list = [],[]
                    # Iterate over every gene
                for gene in cluster_exp.columns:
                        
                        # Compute the log2 fold change
                    log_FC = np.log2(np.mean(covid_group[gene].values)/np.mean(healthy_group[gene].values))
                    log_FC_list.append(log_FC)
                        # Perform a Wilcoxon rank sum test and pull the p-value
                    p_val = ranksums(covid_group[gene].values, healthy_group[gene].values)
                    p_val_list.append(p_val)
                        
                cluster_df = pd.DataFrame({'gene': exp_data_columns,
                                              'logFC': log_FC_list,
                                              'p_val': p_val_list})    
                cluster_df.to_csv("./differential_expression/cluster" + str(clust) + ".csv")
                count += 1

    return rank
    
# Read in data

ini = MPI.Wtime()
rank = DE()
end = MPI.Wtime()

if rank == 0:
    print(end - ini)

