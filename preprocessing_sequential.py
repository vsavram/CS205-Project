# CS205
# Preprocess the single-cell data - Sequential

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns
import os
from scipy.io import mmread


# Import the expression data and the metadata
raw_exp = pd.read_csv("./Data/covid_filtered_counts.csv")
metadata = pd.read_table("./Data/metadata.tsv")


#---------------------------------------------------------------------------------------------------
# Remove lowly expressed genes
#---------------------------------------------------------------------------------------------------

def remove_lowly_expressed(raw_exp, perc_zero):
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
    
    return raw_exp   

#---------------------------------------------------------------------------------------------------
# Perform gene length normalization and sequencing depth normalization
#---------------------------------------------------------------------------------------------------

def gene_length_norm(raw_exp, gene_lengths):
    # Perform gene length normalization
    raw_exp = raw_exp.div(gene_lengths, axis=1)
    return raw_exp

def sequencing_depth_norm(raw_exp):
    # Perform sequencing depth normalization
    sequencing_depth = raw_exp.apply(lambda x: np.sum(x), axis=0)
    raw_exp = raw_exp.div(sequencing_depth, axis=0)
    return raw_exp

#---------------------------------------------------------------------------------------------------
# Log trasnform the data
#---------------------------------------------------------------------------------------------------

def log_transform(raw_exp):
    exp_data = np.log2(raw_exp + 1)    
    return exp_data
