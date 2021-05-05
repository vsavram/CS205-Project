# CS205
# Preprocess the single-cell data - Sequential

import numpy as np
import pandas as pd


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
    raw_exp = raw_exp.div(gene_lengths.to_list(), axis=1)
    return raw_exp

def sequencing_depth_norm(raw_exp):
    # Perform sequencing depth normalization
    sequencing_depth = raw_exp.apply(lambda x: np.sum(x), axis=1)
    raw_exp = raw_exp.div(sequencing_depth, axis=0)
    return raw_exp

#---------------------------------------------------------------------------------------------------
# Log trasnform the data
#---------------------------------------------------------------------------------------------------

def log_transform(raw_exp):
    raw_exp = raw_exp*1e6
    exp_data = np.log2(raw_exp + 1)    
    return exp_data
