# CS205
# Preprocess the single-cell data (PySpark Version)

from pyspark import SparkConf, SparkContext
from pyspark.sql import SparkSession
from pyspark.sql.context import SQLContext
from pyspark.sql import functions as F
from pyspark.sql.types import IntegerType
import numpy as np
import sys
import re


# Define the arguments
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


# Pull the arguments
opts = get_options()

    
spark = SparkSession.builder.appName("Preprocessing").getOrCreate()

# Suppress INFO messages
spark.sparkContext.setLogLevel('WARN')
          
# Import the expression data as a spark dataframe
df = spark.read.option("header",True).csv(opts.raw_data_path)

# Remove the first column which contains the sample ID's
samples = df.select('Unnamed: 0').collect()
df = df.drop('Unnamed: 0')
    
# Import the gene lengths
gene_lengths = np.genfromtxt(opts.gene_length_path, delimiter=',')


#---------------------------------------------------------------------------------------------------
# Remove lowly expressed genes
#---------------------------------------------------------------------------------------------------

# Determine the threshold for removing lowly expressed genes
perc_zero = opts.perc_zero

# Determine the number of zeros in each column
num_zero = df.select([F.sum(F.when(F.col(c) == 0, 1)).alias(c) for c in df.columns])

# Determine the columns to keep
num_samples = df.count()
cols_to_remove = num_zero.select([F.when(F.col(c) > num_samples*perc_zero, True).otherwise(False).alias(c) for c in num_zero.columns])

# Pull the column names to remove
#cols_match = [key for (key, value) in cols_to_remove.collect()[0].asDict().items() if value == True]
cols_match = np.array(cols_to_remove.columns)[np.array(cols_to_remove.collect())[0]]

# Drop columns
df = df.drop(*cols_match)
# Drop gene lengths
gene_lengths = gene_lengths[np.array(cols_to_remove.collect())[0]]


#---------------------------------------------------------------------------------------------------
# Perform gene length normalization and sequencing depth normalization
#---------------------------------------------------------------------------------------------------

# Perform gene length normalization
df = df.select(*((df[x] / gene_len).alias(x) for x,gene_len in zip(df.columns,gene_lengths)))

# Compute the row sums
sequencing_depth = sum([df[x] for x in df.columns])
# Normalize the samples by sequencing depth
df = df.select(*((df[x] / sequencing_depth).alias(x) for x in df.columns))


#---------------------------------------------------------------------------------------------------
# Log transform the data
#---------------------------------------------------------------------------------------------------

# Scale the data by 1e6
df = df.select(*((df[x] * 1e6).alias(x) for x in df.columns))
# Log transform the data
df = df.select(*((F.log2(df[x] + 1)).alias(x) for x in df.columns))

