# Hybrid-Parallelism for Single-Cell Gene Expression Analyses

![gene](https://user-images.githubusercontent.com/29682604/117387165-4e8bd580-aeb6-11eb-96ec-fefc0f273e3e.png?style=centerme)

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

The data used to assess the performance and scalability of our programming model was taken from Y. Zheng et al 2020 and can be found [here](https://www.ncbi.nlm.nih.gov/pmc/articles/PMC7417788/). Samples are taken from human blood. The dataset is advantageous to work with because human blood is comprised of a diverse set of immune cells, many of which are well documented. We show that our programming model handles a heterogenous dataset as is typically found when working with single-cell data. The dataset is in tabular form, spanning 30,000 cells (rows) and 20,896 genes (columns). Samples were taken from patients who were infected with COVID-19 and from healthy donors. 

## Hybrid Parallelization Programming Model

A hybrid parallelization scheme is used to optimize the single-cell pipeline. Big data flow processing using the Spark framework is used to optimize the preprocessing functions. The message passing interface (MPI) framework is used to optimize clustering and differential expression. The details of the programming model are outlined in the following sections.

### Big Data Flow Preprocessing with Spark

![Screen Shot 2021-05-06 at 11 36 03 PM](https://user-images.githubusercontent.com/29682604/117394283-e2b06980-aec3-11eb-8001-44dafa8edc51.png)

### Cell Clustering with MPI

### Differential Expression Analysis with MPI

