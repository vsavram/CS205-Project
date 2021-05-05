# CS205
# Perform differential expression analysis - Sequential

import numpy as np
import pandas as pd
from scipy.stats import ranksums
import os

def DE(exp_data, clusters, metadata):
    
    clusters = np.array(clusters).flatten()
    
    # Create the directory used to store the differential expression results if it does not exist
    if not os.path.exists('differential_expression'):
        os.makedirs('differential_expression')
        
    # Iterate over every cluster
    for clust in np.unique(clusters):
        
        # Filter the expression data and metadata to only include cells from the given cluster
        cluster_exp = exp_data.loc[clusters == clust,:]
        cluster_metadata = metadata.loc[clusters == clust,:]
        
        # Determine the data associated with the COVID group and healthy group
        covid_group = cluster_exp.loc[(cluster_metadata['Sample Characteristic[disease]'] == 'COVID-19').to_numpy(),:]
        healthy_group = cluster_exp.loc[(cluster_metadata['Sample Characteristic[disease]'] == 'normal').to_numpy(),:]
        
        # Create lists used to store the log2 fold changes and p-values associated with each gene
        log_FC_list,p_val_list = [],[]
        # Iterate over every gene
        for gene in cluster_exp.columns:
            
            # Compute the log2 fold change
            log_FC = np.log2(np.mean(covid_group[gene])/np.mean(healthy_group[gene]))
            log_FC_list.append(log_FC)
            # Perform a Wilcoxon rank sum test and pull the p-value
            p_val = ranksums(covid_group[gene], healthy_group[gene])
            p_val_list.append(p_val)
            
        cluster_df = pd.DataFrame({'gene': exp_data.columns,
                                  'logFC': log_FC_list,
                                  'p_val': p_val_list})    
        cluster_df.to_csv("./differential_expression/cluster" + str(clust) + ".csv")
        