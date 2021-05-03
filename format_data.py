# CS205
# Convert the data to csv format

import numpy as np
import pandas as pd
import os
from scipy.io import mmread


os.chdir("/Users/VICTOR/Desktop/Harvard/AC205/project")


# Import the raw expression data along with the row and column names
raw_exp = mmread("./Data/raw_files/covid_filtered_counts.mtx")
raw_exp = pd.DataFrame(raw_exp.toarray())  # Convert the sparse matrix to a dataframe
raw_samples = pd.read_table("./Data/raw_files/covid_exp_cols.txt", header=None).to_numpy().flatten()
raw_genes = pd.read_table("./Data/raw_files/covid_exp_rows.txt", header=None).to_numpy().flatten()
# Format the raw expression dataframe
raw_exp.index = raw_genes
raw_exp.columns = raw_samples

# Transpose the expression data
raw_exp = raw_exp.transpose()
raw_exp.to_csv("./Data/raw_files/covid_filtered_counts.csv")