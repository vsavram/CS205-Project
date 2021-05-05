# Run file for sequential execution

import pandas as pd
import numpy as np
import argparse
import time

from preprocessing_sequential import *
from clustering_sequential import *
from create_tsne import *
from DE_sequential import DE


def get_options(args=None):
    parser = argparse.ArgumentParser(description="Sequential Execution")

    # Define the arguments for specifying the data location
    parser.add_argument('--raw_data_path', type=str, 
                        default='../data/raw_files/covid_filtered_counts.csv', 
                        help='Path to raw expression data')
    parser.add_argument('--metadata_path', type=str, 
                    default='../data/metadata.tsv', 
                    help='Path to metadata')
    parser.add_argument('--gene_length_path', type=str, 
                default='../data/gene_lengths.csv', 
                help='Path to the gene lengths')

    # Define an argument that specifies the threshold for removing lowly expressed genes
    parser.add_argument('--perc_zero', type=float, default=0.95, 
                        help='Threshold for removing lowly expressed genes')
    
    # Define an argument that specifies the type of clustering
    parser.add_argument('--clustering', type=str, default='kmeans', 
                        help='The type of clustering [kmeans, louvain]')

    opts = parser.parse_args(args)

    return opts


if __name__ == "__main__":

    # Pull the arguments
    opts = get_options()
    
    # Import the raw expression data and the metadata
    exp_data = pd.read_csv(opts.raw_data_path)
    metadata = pd.read_table(opts.metadata_path)
    metadata = metadata[['Assay', 'Sample Characteristic[disease]']]
    # Import the gene lengths
    gene_lengths = pd.read_csv(opts.gene_length_path)['gene length'].to_numpy()
    
    start_time = time.time()
    
    # Preprocess the raw data
    exp_data,gene_lengths = remove_lowly_expressed(exp_data, opts.perc_zero, gene_lengths)
    exp_data = gene_length_norm(exp_data, gene_lengths)
    exp_data = sequencing_depth_norm(exp_data)
    exp_data = log_transform(exp_data)
    
    end_time = time.time()
    preprocessing_time = end_time - start_time
    print(f"The execution time for sequential preprocessing is {preprocessing_time:0.4f} seconds.")
    
    start_time = time.time()

    # Peform clustering
    assert (opts.clustering in ['kmeans', 'louvain']), 'Did not specify a proper clustering method'
    if opts.clustering == 'kmeans':
        clusters = kmeans_clustering(exp_data, 10)
    else:
        clusters = louvain_clustering(exp_data)
       
    end_time = time.time()
    clustering_time = end_time - start_time
    print(f"The execution time for sequential clustering is {clustering_time:0.4f} seconds.")

    start_time = time.time()

    # Create a TSNE plot
    create_tsne(exp_data, clusters)
    
    end_time = time.time()
    tsne_time = end_time - start_time
    print(f"The execution time for TSNE visualization is {tsne_time:0.4f} seconds.")

    start_time = time.time()

    # Perform differential expression analysis
    DE(exp_data, clusters, metadata)
    
    end_time = time.time()
    DE_time = end_time - start_time
    print(f"The execution time for differential expression analysis is {DE_time:0.4f} seconds.")
    
    # Save the execution times for each component of the analysis
    execution_time_df = pd.DataFrame({'component': ['preprocessing', 'clustering', 'tsne', 'DE'], 
                                      'time (s)': [preprocessing_time, clustering_time, tsne_time, DE_time]})
    execution_time_df.to_csv("sequential_execution_times.csv")
