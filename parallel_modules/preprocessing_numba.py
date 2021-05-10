# CS205
# Preprocess the single-cell data - CUDA

import pandas as pd
import numpy as np
import math
import argparse
import time
import numba
from numba import njit, prange

# Import the data
def get_options(args=None):
    parser = argparse.ArgumentParser(description="Sequential Execution")

    # Define the arguments for specifying the data location
    parser.add_argument('--raw_data_path', type=str, 
                        default='../data/raw_files/covid_filtered_counts.csv', 
                        help='Path to raw expression data')
    parser.add_argument('--gene_length_path', type=str, 
                default='../data/gene_lengths.csv', 
                help='Path to the gene lengths')

    # Define an argument that specifies the threshold for removing lowly expressed genes
    parser.add_argument('--perc_zero', type=float, default=0.95, 
                        help='Threshold for removing lowly expressed genes')

    opts = parser.parse_args(args)

    return opts

    
@njit(parallel=True)
def define_low_genes(raw_exp):
    # Determine lowly expressed genes (not expressed in over a given percentage of samples)
    lowly_expressed = np.zeros(raw_exp.shape[1])
    # Iterate over every column in the expression data
    for i in prange(raw_exp.shape[1]):
        # Determine the number of samples that do not express the given gene
        num_no_expression = np.sum(raw_exp[:,i] == 0)
        lowly_expressed[i] = num_no_expression
        
    return lowly_expressed    

def remove_low_genes(raw_exp, gene_lengths, lowly_expressed, perc_zero=0.95):
    cutoff = perc_zero*raw_exp.shape[1]
    # Remove lowly expressed genes
    raw_exp = raw_exp[:,~(np.array(lowly_expressed) > cutoff)]
    gene_lengths = gene_lengths[~~(np.array(lowly_expressed) > cutoff)]
    
    return raw_exp,gene_lengths
    
    
@njit(parallel=True)
def preprocess_njit(raw_exp, gene_lengths, save=False):
    """
    Parameters
    ----------
    raw_exp : pandas dataframe
        a cells by genes dataframe containing raw expression values.
    gene_lengths : numpy array
        an array containing the gene length for each gene in raw_exp.

    Returns
    -------
    preprocessed_exp : pandas dataframe
        the preprocessed expression values.

    """
    preprocessed_exp = np.ones(raw_exp.shape)
    
    # Perform gene length normalization
    for i in prange(raw_exp.shape[0]):
        for j in range(raw_exp.shape[1]):
            preprocessed_exp[i,j] = raw_exp[i,j]/gene_lengths[j]     
    
    
    # Perform sequencing depth normalization
    for i in prange(preprocessed_exp.shape[0]):
        sequencing_depth = np.sum(preprocessed_exp[i,:])
        preprocessed_exp[i,:] = preprocessed_exp[i,:]/sequencing_depth     
    
    
    # Log-transform the data
    for i in prange(preprocessed_exp.shape[0]):
        for j in range(preprocessed_exp.shape[1]):
            preprocessed_exp[i,j] = np.log2(preprocessed_exp[i,j]*1e6 + 1)  
    
        
    return preprocessed_exp



if __name__ == "__main__":

    # Pull the arguments
    opts = get_options()
    
    # Import the raw expression data and the metadata
    exp_data = pd.read_csv(opts.raw_data_path, index_col=0).to_numpy()
    # Import the gene lengths
    gene_lengths = np.genfromtxt(opts.gene_length_path, delimiter=',')
    
    start_time = time.time()
    
    # Preprocess the raw data
    lowly_expressed = define_low_genes(exp_data)
    exp_data,gene_lengths = remove_low_genes(exp_data, gene_lengths, lowly_expressed, opts.perc_zero)
    preprocessed_exp = preprocess_njit(exp_data, gene_lengths)
    
    end_time = time.time()
    preprocessing_time = end_time - start_time
    print("The execution time for preprocessing using numba is {0}".format(preprocessing_time))
    
    if save==True:
        np.savetxt("preprocessed_counts.csv", preprocessed_exp, delimiter=",")
    