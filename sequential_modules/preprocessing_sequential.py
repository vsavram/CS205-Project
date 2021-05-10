# CS205
# Preprocess the single-cell data - Sequential

import numpy as np
import pandas as pd


### raw_exp is an n x p matrix
### gene_lengths is a vector of length p


#---------------------------------------------------------------------------------------------------
# Remove lowly expressed genes
#---------------------------------------------------------------------------------------------------

def remove_lowly_expressed(raw_exp, perc_zero, gene_lengths):
    # Determine lowly expressed genes (not expressed in over a given percentage of samples)
    lowly_expressed = []
    cutoff = perc_zero*raw_exp.shape[1]
    
    # Iterate over every column in the expression data
    for col in raw_exp.columns:
        data = raw_exp[col]
        # Determine the number of samples that do not express the given gene
        num_no_expression = np.sum(data == 0)
        lowly_expressed.append(num_no_expression >= cutoff)
        
    # Remove lowly expressed genes
    raw_exp = raw_exp.loc[:,~np.array(lowly_expressed)]
    gene_lengths = gene_lengths[~np.array(lowly_expressed)]
    
    return raw_exp,gene_lengths

#---------------------------------------------------------------------------------------------------
# Perform gene length normalization and sequencing depth normalization
#---------------------------------------------------------------------------------------------------

def gene_length_norm(raw_exp, gene_lengths):
    # Perform gene length normalization
    for gene,gene_len in zip(raw_exp.columns,gene_lengths):
        raw_exp[gene] = raw_exp[gene]/gene_len
    return raw_exp

def sequencing_depth_norm(raw_exp):
    # Perform sequencing depth normalization
    for i in range(raw_exp.shape[0]):
        sequencing_depth = sum(raw_exp.iloc[i,:])
        raw_exp.iloc[i,:] = raw_exp.iloc[i,:]/sequencing_depth
    return raw_exp

#---------------------------------------------------------------------------------------------------
# Log trasnform the data
#---------------------------------------------------------------------------------------------------

def log_transform(raw_exp):
    for i in range(raw_exp.shape[0]):
        for j in range(raw_exp.shape[1]):
            raw_exp.iloc[i,j] = np.log2(raw_exp.iloc[i,j]*1e6 + 1)  
    return raw_exp