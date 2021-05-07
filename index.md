# Hybrid-Parallelism for Single-Cell Gene Expression Analyses

![gene](https://user-images.githubusercontent.com/29682604/117387165-4e8bd580-aeb6-11eb-96ec-fefc0f273e3e.png)

## Introduction and Motivation

Advances in next-generation sequencing (NGS) technology have allowed researchers to perform high-throughput sequencing at the single cell level. As opposed to bulk gene expression sequencing, single-cell sequencing generates gene expression profiles for every cell (to an extent) in a given sample. These data provide insight into the dynamics of and interactions between different cell populations as well as the discrepancies between cell populations/subpopulations across different experimental conditions. single-cell RNA sequencing (scRNA-seq) is for both rapidly generating hypotheses and for drawing biologically meaningful insights. Heavy reductions in the cost per cell across multiple sequencing platforms (e.g. 10X Genomics, SORT-seq) have made this type of data generation accesible to many labs/institutions. Furthermore, NGS data has become an intergral part of biomedical research and a complement to wet-lab experimentation. 

scRNA-seq data are large, typically spanning 10,000's or 100,000's of cells and 10,000's of genes. As well, these datasets are inherently complex and are computationally intensive to analyze. The problem can be framed in the context of big data and big compute, necessitating savvy ways of handling this data. We have built a hypbrid-parallel framework that integrates into the standard single-cell analysis pipeline from data preprocessing to performing differential gene expression analysis. We outline this standard pipeline below, before discussing the programming model used to optimize these analyses.

##### The Standard Single-Cell Analysis Pipeline

A standard scRNA-seq workflow involves the following steps:
1. **Preprocessing** - Removing lowly expressed genes, performing gene length/sequencing depth normalization, and log-transforming the data.
2. **Clustering** - Clustering cells
3. **Visualization** - Visualizing the clusters on a dimensionality-reduced plot, typically using a non-linear dimensionality reduction technique such as tSNE or UMAP.
4. **Cluster Annotation** - Mapping each cluster to a specific cell type (e.g. cluster 1 corresponds to macrophages).
5. **Differential Expression Analysis** - Determining which genes in a given cluster are expressed at significantly different levels between two groups.

![single_cell_pipeline](https://user-images.githubusercontent.com/29682604/117387384-d245c200-aeb6-11eb-8812-def86d7d02ee.png)

```markdown
Syntax highlighted code block

# Header 1
## Header 2
### Header 3

- Bulleted
- List

1. Numbered
2. List

**Bold** and _Italic_ and `Code` text

[Link](url) and ![Image](src)
```

For more details see [GitHub Flavored Markdown](https://guides.github.com/features/mastering-markdown/).
