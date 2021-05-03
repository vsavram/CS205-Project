# CS205
# Preprocess the single-cell data (PySpark Version)

from pyspark import SparkConf, SparkContext
from pyspark.sql import SparkSession
from pyspark.sql.context import SQLContext
from pyspark.sql import functions as F
from pyspark.sql.types import IntegerType
import sys
import re

spark = SparkSession.builder.master("local").appName("Preprocessing").getOrCreate()

# Suppress INFO messages
spark.sparkContext.setLogLevel('WARN')
          
# Import the expression data as a spark dataframe
df = spark.read.option("header",True).csv("./Data/raw_files/covid_filtered_counts.csv")

# Remove the first column which contains the sample ID's
samples = df.select('_c0').collect()
df = df.drop('_c0')

# Set the Spark configurations
conf = SparkConf().setMaster('local').setAppName('Preprocessing')
sc = SparkContext(conf = conf)
# Set the SQL context and the Spark session
sqlContext = SQLContext(sc)
spark = SparkSession(sc)

# Import the expression data 
df = sqlContext.read.csv("./Data/raw_files/covid_filtered_counts.csv", header=True)

# Create a RDD based on the input file
RDDvar = sc.textFile("./Data/raw_files/covid_filtered_counts.csv")

    
    
#---------------------------------------------------------------------------------------------------
# Remove lowly expressed genes
#---------------------------------------------------------------------------------------------------

# Replace zero values with None
df = df.replace(0,None)
# 

# Define a function used to compute the number of zero elements in a given column
def num_zero(col):
  return np.sum(np.array(col) == 0)
# Convert the function to a PySpark UDF
spark.udf.register("num_zero", num_zero, IntegerType())

# Determine the number of samples
num_samples = df.count()

# Determine lowly expressed genes (not expressed in over 95% of samples) and drop them
df = df.select([count(when(isnull(c), c)).alias(c) for c in df.columns]).show()


# Determine lowly expressed genes (not expressed in over 95% of samples) and drop them
columns_to_drop = [col for col in df.schema.names if num_zero(df.select(col).collect()) > num_samples*0.95]
df = df.drop(*columns_to_drop)
                    
# Determine lowly expressed genes (not expressed in over 95% of samples)
rows_to_remove = df.select(['recclass', 'mass (g)'])

# Filter out rows where the mass is unknown
df = df.filter(df['mass (g)'].isNotNull())

# Compute the average mass per meteor type
df = df.groupby('recclass').agg(F.avg('mass (g)')).orderBy('recclass')

# Print the dataframe to the console
df.persist
df.show(df.count(), False)



# Import the raw expression data along with the row and column names
raw_exp = mmread("./Data/raw_files/covid_filtered_counts.mtx")
raw_exp = pd.DataFrame(raw_exp.toarray())  # Convert the sparse matrix to a dataframe
raw_samples = pd.read_table("./Data/raw_files/covid_exp_cols.txt", header=None).to_numpy().flatten()
raw_genes = pd.read_table("./Data/raw_files/covid_exp_rows.txt", header=None).to_numpy().flatten()
# Format the raw expression dataframe
raw_exp.index = raw_genes
raw_exp.columns = raw_samples

raw_exp.to_csv("./Data/raw_files/covid_filtered_counts.csv")

# Import the metadata
metadata = pd.read_table("./Data/metadata.tsv")


#---------------------------------------------------------------------------------------------------
# Remove lowly expressed genes
#---------------------------------------------------------------------------------------------------

# Determine lowly expressed genes (not expressed in over 95% of samples)
lowly_expressed = []
cutoff = 0.95*raw_exp.shape[0]
for row_index in range(raw_exp.shape[0]):
    data = raw_exp[gene]
    # Determine the number of samples that do not express the given gene
    num_no_expression = np.sum(raw_exp[row_index] == 0)
    lowly_expressed.append(num_no_expression >= cutoff)
    
# Remove lowly expressed genes
raw_exp = tpm_data.loc[~np.array(lowly_expressed),:]
raw_genes = raw_genes[~np.array(lowly_expressed)]

# Log transform the tpm data
log_tpm = np.log2(tpm + 1)


#---------------------------------------------------------------------------------------------------
# Perform gene length normalization and sequencing depth normalization
#---------------------------------------------------------------------------------------------------

# Perform gene length normalization
gene_lengths
raw_exp = raw_exp.div(gene_lengths, axis=0)

# Perform sequencing depth normalization
sequencing_depth = raw_exp.apply(lambda x: np.sum(x), axis=1)
raw_exp = raw_exp.div(sequencing_depth, axis=1)


#---------------------------------------------------------------------------------------------------
# Log trasnform the data
#---------------------------------------------------------------------------------------------------

raw_exp = np.log2(raw_exp + 1)
