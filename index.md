# Hybrid-Parallelism for Single-Cell Gene Expression Analyses
<p align="center">
  <img src="https://user-images.githubusercontent.com/29682604/117387165-4e8bd580-aeb6-11eb-96ec-fefc0f273e3e.png">
</p>

## Introduction and Motivation

Advances in next-generation sequencing (NGS) technology have allowed researchers to perform high-throughput sequencing at the single cell level. As opposed to bulk gene expression sequencing, single-cell sequencing generates gene expression profiles for every cell (to an extent) in a given sample. These data provide insight into the dynamics of and interactions between different cell populations as well as the discrepancies between cell populations/subpopulations across different experimental conditions. single-cell RNA sequencing (scRNA-seq) is for both rapidly generating hypotheses and for drawing biologically meaningful insights. Heavy reductions in the cost per cell across multiple sequencing platforms (e.g. 10X Genomics, SORT-seq) have made this type of data generation accesible to many labs/institutions. Furthermore, NGS data has become an intergral part of biomedical research and a complement to wet-lab experimentation. 

scRNA-seq data are large, typically spanning 10,000's or 100,000's of cells and 10,000's of genes. As well, these datasets are inherently complex and are computationally intensive to analyze. The problem can be framed in the context of big data and big compute, necessitating savvy ways of handling this data. We have built a hypbrid-parallel framework that integrates into the standard single-cell analysis pipeline from data preprocessing to performing differential gene expression analysis. We outline this standard pipeline below, before discussing the programming model used to optimize these analyses.

## The Standard Single-Cell Analysis Pipeline

A standard scRNA-seq workflow involves the following steps:
1. **Preprocessing** - Removing lowly expressed genes, performing gene length/sequencing depth normalization, and log-transforming the data.
2. **Clustering** - Clustering cells
3. **Visualization** - Visualizing the clusters on a dimensionality-reduced plot, typically using a non-linear dimensionality reduction technique such as tSNE or UMAP.
4. **Cluster Annotation** - Mapping each cluster to a specific cell type (e.g. cluster 1 corresponds to macrophages).
5. **Differential Expression Analysis** - Determining which genes in a given cluster are expressed at significantly different levels between two groups.

![single_cell_pipeline](https://user-images.githubusercontent.com/29682604/117387384-d245c200-aeb6-11eb-8812-def86d7d02ee.png)

Raw single-cell data is initially preprocessed, prepping it for downstream analyses. Initially, lowly expressed genes are removed. Lowly expressed genes are determined by the proportion of cells that do not express the given gene, using a pre-specified threshold (e.g. 95%). Preprocessing includes gene length normalization and sequencing depth normalization in order to make data comparable across genes and across cells. The last step of preprocessing log transforms the data with an initial small offset (typically an addition of 1). The expression data is subsequently clustered. Common approaches include $k$-means and Louvain/Leiden clustering. Ideally, cells roughly together based on cell type. Clusters are annotated using marker gene expression levels. Marker genes are gene that are exclusively (ideally) expressed in one or a few cell types. Cluster annotation assigns a cell type to each cluster. For example, cluster 1 might be identified as macrophages and cluset 2 might be identified as neutrophils. If samples are derived from different experimental conditions, the last step in the standard single-cell analysis pipeline compares the expression levels across groups. Differential expression analysis determines which genes are expressed at significantly different levels between two groups. 

## Data

The data used to assess the performance and scalability of our programming model was taken from Y. Zheng et al 2020 and can be found [here](https://www.ncbi.nlm.nih.gov/pmc/articles/PMC7417788/). Samples are taken from human blood. The dataset is advantageous to work with because human blood is comprised of a diverse set of immune cells, many of which are well documented. We show that our programming model handles a heterogenous dataset as is typically found when working with single-cell data. The dataset is in tabular form, spanning 30,000 cells (rows) and 20,896 genes (columns). Samples were taken from patients who were infected with COVID-19 and from healthy donors. Associated data files include the cell metadata, indicating whether a given cell came from a COVID-19 infected patient or from a healthy donor, and a file containing the gene lengths for the 20,896 genes. 

## Hybrid Parallelization Programming Model

A hybrid parallelization scheme is used to optimize the single-cell pipeline. Big data flow processing using the Spark framework is used to optimize the preprocessing functions. The message passing interface (MPI) framework is used to optimize clustering and differential expression. Visualization using tSNE is left as a sequential implementation. The details of the programming model are outlined in the following sections.

All implementations are built in Python. 

### Big Data Flow Preprocessing with Spark

The Spark programming model is used to optimize data preprocessing. Specifically, PySpark is used in order convert the expression matrix into a resilient distributed dataset (RDD) that can be partitioned across multiple nodes. The preprocessing stage is broken down into column-wise operations and row-wise operations. 
<img align="right" src="https://user-images.githubusercontent.com/29682604/117394283-e2b06980-aec3-11eb-8001-44dafa8edc51.png">
Given an expression matrix of cells by genes, the columns are initially filtered, removing genes that are not expressed in a large proportion of cells (the user specifies this proportion). Subsequently, every column is divided by the length of the corresponding gene (length in kilobases) in order to make the data comparable across genes. Every row is divided by the summation of the given row (i.e. the sequencing depth) making the data comparable across cells. Lastly, the expression data is multiplied by 1e6 and log-transformed with an initial addition of 1 in order to adjust for outliers. 

### Cell Clustering with MPI

### Differential Expression Analysis with MPI

## Runing the Model and Reproducibility Information

The following sections provide information regarding how to run the sequential and parallel implementations as well as information on how to reproduce the performance results. Implementations were developed in Python and run using an AWS infrastructure. EC2 instances were used to perform computations and S3 was used for storage of data files.

### Sequential Implementation

**Reproducibility Information:** The sequential implemenation was run on an m5.2xlarge AWS instance using the AMI \textbf{Ubuntu Server 20.04 LTS (HVM), SSD Volume Type}. The linux kernel version is 5.4.0-1038-aws. The instance has 8 vCPU's, 8 cores in total (i.e. 1 core per vCPU), 32 GiB of main memory, 32 K of L1d cache memory, 32 K of L1i cache memory, 256 K of L2 cache memory, and 46080 K of L3 cache memory. The CPU clock rate is 2.3 GHz. \\
By default, the m5.2xlarge instance has 8 G of disk space. Given that m5 intances are back by EBS, the disk space can by dynamically resized. Resizing was not needed for this execution because the 

The Python version used is 2.7.17. The following dependencies are required and can be installed by running the following command.
> $ pip install -r requirements.txt
Dependencies:
* matplotlib 
* pandas 
* numpy 
* anndata 
* sklearn 
* scprep 
* scipy 

The following files must be in the same directory in order to execute the sequential implementation.
* `run_sequential.py`
* `preprocessing_sequential.py`
* `clustering_sequential.py`
* `create_tsne.py`
* `DE_sequential.py`

Sequential execution can be run using the following command where --raw_data_path specifies the path to the raw cells by genes expression matrix, --metadata_path specifies the path to the metadata file containing cell-specific metadata, and --gene_length_path specifies the path to the file containing the gene lengths.
> $ ./run_sequential.py --raw_data_path 'Data/covid_filtered_counts_subset.csv' --metadata_path 'Data/metadata_subset.csv' --gene_length_path 'Data/gene_lengths.csv'


The Java version used is 1.8.0\_282. \\
The Scala version used is 2.11.12. \\
The Python version used is 2.7.17. \\
The Spark version used is 2.2.0 \\
