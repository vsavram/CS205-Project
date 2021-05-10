# Hybrid-Parallelism for Single-Cell Gene Expression Analyses
<h3 align="center">
Contributors: John Alling, Victor Avram, Diego Zertuche
</h3>
<p align="center">
  <img src="https://user-images.githubusercontent.com/29682604/117387165-4e8bd580-aeb6-11eb-96ec-fefc0f273e3e.png">
</p>

## Introduction and Motivation
---

Advances in next-generation sequencing (NGS) technology have allowed researchers to perform high-throughput sequencing at the single cell level. As opposed to bulk gene expression sequencing, single-cell sequencing generates gene expression profiles for every cell (to an extent) in a given sample. These data provide insight into the dynamics of and interactions between different cell populations as well as the discrepancies between cell populations/subpopulations across different experimental conditions. single-cell RNA sequencing (scRNA-seq) is for both rapidly generating hypotheses and for drawing biologically meaningful insights. Heavy reductions in the cost per cell across multiple sequencing platforms (e.g. 10X Genomics, SORT-seq) have made this type of data generation accesible to many labs/institutions. Furthermore, NGS data has become an intergral part of biomedical research and a complement to wet-lab experimentation. 

scRNA-seq data are large, typically spanning 10,000's or 100,000's of cells and 10,000's of genes. As well, these datasets are inherently complex and are computationally intensive to analyze. The problem can be framed in the context of big data and big compute, necessitating savvy ways of handling this data. We have built a hypbrid-parallel framework that integrates into the standard single-cell analysis pipeline from data preprocessing to performing differential gene expression analysis. We outline this standard pipeline below, before discussing the programming model used to optimize these analyses.

### The Standard Single-Cell Analysis Pipeline

A standard scRNA-seq workflow involves the following steps:
1. **Preprocessing** - Removing lowly expressed genes, performing gene length/sequencing depth normalization, and log-transforming the data.
2. **Clustering** - Clustering cells
3. **Visualization** - Visualizing the clusters on a dimensionality-reduced plot, typically using a non-linear dimensionality reduction technique such as tSNE or UMAP.
4. **Cluster Annotation** - Mapping each cluster to a specific cell type (e.g. cluster 1 corresponds to macrophages).
5. **Differential Expression Analysis** - Determining which genes in a given cluster are expressed at significantly different levels between two groups.

![single_cell_pipeline](https://user-images.githubusercontent.com/29682604/117387384-d245c200-aeb6-11eb-8812-def86d7d02ee.png)

Raw single-cell data is initially preprocessed, prepping it for downstream analyses. Initially, lowly expressed genes are removed. Lowly expressed genes are determined by the proportion of cells that do not express the given gene, using a pre-specified threshold (e.g. 95%). Preprocessing includes gene length normalization and sequencing depth normalization in order to make data comparable across genes and across cells. The last step of preprocessing log transforms the data with an initial small offset (typically an addition of 1). The expression data is subsequently clustered. Common approaches include $k$-means and Louvain/Leiden clustering. Ideally, cells roughly together based on cell type. Clusters are annotated using marker gene expression levels. Marker genes are gene that are exclusively (ideally) expressed in one or a few cell types. Cluster annotation assigns a cell type to each cluster. For example, cluster 1 might be identified as macrophages and cluset 2 might be identified as neutrophils. If samples are derived from different experimental conditions, the last step in the standard single-cell analysis pipeline compares the expression levels across groups. Differential expression analysis determines which genes are expressed at significantly different levels between two groups. 

## Data
---

The data used to assess the performance and scalability of our programming model was taken from Y. Zheng et al 2020 and can be found [here](https://www.ncbi.nlm.nih.gov/pmc/articles/PMC7417788/) [1]. Samples are taken from human blood. The dataset is advantageous to work with because human blood is comprised of a diverse set of immune cells, many of which are well documented. We show that our programming model handles a heterogenous dataset as is typically found when working with single-cell data. The dataset is in tabular form, spanning 30,000 cells (rows) and 20,896 genes (columns). Samples were taken from patients who were infected with COVID-19 and from healthy donors. Associated data files include the cell metadata, indicating whether a given cell came from a COVID-19 infected patient or from a healthy donor, and a file containing the gene lengths for the 20,896 genes. 

## Hybrid Parallelization Programming Model
---

A hybrid parallelization scheme is used to optimize the single-cell pipeline. Big data flow processing using the Spark framework is used to optimize the preprocessing functions. In order to alleviate communication overhead associated with a network of nodes, Numba (a compiler for python that performs just-in-time translation of numeric functions) is used to enhance the preprocessing optimization.The message passing interface (MPI) framework is used to optimize clustering and differential expression, specifically using mpi4py. Visualization using tSNE is left as a sequential implementation. The details of the programming model are outlined in the following sections.

All implementations are built in Python. 

### Big Data Flow Preprocessing with Spark

<img align="right" src="https://user-images.githubusercontent.com/29682604/117394283-e2b06980-aec3-11eb-8001-44dafa8edc51.png">
The Spark programming model is used to optimize data preprocessing. Specifically, PySpark is used in order convert the expression matrix into a resilient distributed dataset (RDD) that can be partitioned across multiple nodes. The preprocessing stage is broken down into column-wise operations and row-wise operations. Given an expression matrix of cells by genes, the columns are initially filtered, removing genes that are not expressed in a large proportion of cells (the user specifies this proportion). Subsequently, every column is divided by the length of the corresponding gene (length in kilobases) in order to make the data comparable across genes. Every row is divided by the summation of the given row (i.e. the sequencing depth) making the data comparable across cells. Lastly, the expression data is multiplied by 1e6 and log-transformed with an initial addition of 1 in order to adjust for outliers. 

### Preprocessing Optimization with Numba

Numba is a compiler for Python that specifically works with numerical functions to accelerate computing. Numba uses the LLVM compiler architecture to perform this optimization. THe just-in-time translation of a given function (the bytecode is translated to machine code just before execution), causes improved execution speeds. Numba is used to peform this preprocessing optimzation. Furthermore, Numba is also used as a shared memory multi-threading framework, able to provide local parallelism. With a fine-grained shared memory approach, we can circumvent the overhead associated with distsributed models that require communication and synchronization across nodes. 

### Cell Clustering with MPI

Cells are grouped into clusters based on like characteristics from the preprocessing step. K-Means clustering is a common algorithm for this approach, which is parallelized using the MPI framework. The MPI programming model allows different nodes in a cluster to communicate with each other through a 'message passing interface', where data can be shared, distributed and collected throughout the different computing nodes to allow computations to be carried out to different chunks of data at different nodes, and thus achieve parallelization. An important step of K-Means is finding the L2 norm of each data point with respect to a centroid. This makes for a good opportunity for parallelization, so the data is split such that each node receives an equal portion of the data, increasing the compute throughput overall.

### Differential Expression Analysis with MPI

Similarly to the cell clustering stage, the MPI framework can be applied to differential gene expression analysis. For the differential expression analysis, a statistical analysis, in this case a rank sum test, is performed independently for each cluster of genes to determine if there is a statistical difference in gene expression in the control and study group in that particular cluster. This leads to a natural integration to the MPI programming model, where a main node will send one cluster data to each node so the rank sum tests of each cluster can be performed parallelly.

## Runing the Model and Reproducibility Information
---

The following sections provide information regarding how to run the sequential and parallel implementations as well as information on how to reproduce the performance results. Implementations were developed in Python and run using an AWS infrastructure. EC2 instances were used to perform computations and S3 was used for storage of data files.

### Sequential Implementation

**Reproducibility Information:** The sequential implemenation was run on an m5.2xlarge AWS instance using the AMI \textbf{Ubuntu Server 20.04 LTS (HVM), SSD Volume Type}. The linux kernel version is 5.4.0-1038-aws. The instance has 8 vCPU's, 8 cores in total (i.e. 1 core per vCPU), 32 GiB of main memory, 32 K of L1d cache memory, 32 K of L1i cache memory, 256 K of L2 cache memory, and 46080 K of L3 cache memory. The CPU clock rate is 2.3 GHz. \\
By default, the m5.2xlarge instance has 8 G of disk space. Given that m5 intances are back by EBS, the disk space can by dynamically resized. Resizing was not needed for this execution.

The Python version used is 2.7.17. The following dependencies are required and can be installed by running the following command.
> $ pip install -r requirements.txt

Dependencies:
* matplotlib 
* pandas 
* numpy 
* argparse
* anndata 
* sklearn 
* scprep 
* scanpy
* scipy 

The following files must be in the same directory in order to execute the sequential implementation.
* `run_sequential.py`
* `preprocessing_sequential.py`
* `clustering_sequential.py`
* `create_tsne.py`
* `DE_sequential.py`

Sequential execution can be run using the following command where `--`raw_data_path specifies the path to the raw cells by genes expression matrix, `--`metadata_path specifies the path to the metadata file containing cell-specific metadata, and `--`gene_length_path specifies the path to the file containing the gene lengths.
> $ ./run_sequential.py  `--`raw_data_path  'Data/covid_filtered_counts_subset.csv'  `--`metadata_path  'Data/metadata_subset.csv'  `--`gene_length_path  'Data/gene_lengths.csv'

The execution times for each step of the preprocessing are printed to the console. The executions times are also saved to a file **sequential_execution_times.csv**.

### PySpark Preprocessing

**Reproducibility Information:** The PySpark implemenation of preprocessing was run on a cluster of m5.2xlarge AWS instances using the AMI **Ubuntu Server 20.04 LTS (HVM), SSD Volume Type**. The linux kernel version is 5.4.0-1038-aws. The instance has 8 vCPU's, 8 cores in total (i.e. 1 core per vCPU), 32 GiB of main memory, 32 K of L1d cache memory, 32 K of L1i cache memory, 256 K of L2 cache memory, and 46080 K of L3 cache memory. The CPU clock rate is 2.3 GHz. \\
By default, the m5.2xlarge instance has 8 G of disk space. Given that m5 intances are back by EBS, the disk space can by dynamically resized. Resizing was not needed for this execution. \\
Up to 8 nodes were used for this implementation. The network bandwidth is up to 10 Gbps. 

The Java version used is 1.8.0\_282. \\
The Scala version used is 2.11.12. \\
The Python version used is 2.7.17. \\
The Spark version used is 2.2.0 \\

The Python version used is 2.7.17. The following dependencies are required and can be installed by running the following command.
> $ pip install -r requirements_pyspark.txt

PySpark preprocessing can be run using the following command where `--`raw_data_path specifies the path to the raw cells by genes expression matrix, `--`metadata_path specifies the path to the metadata file containing cell-specific metadata, and `--`gene_length_path specifies the path to the file containing the gene lengths.
> $ spark-submit preprocessing_spark.py  `--`raw_data_path  'Data/covid_counts.csv'  `--`metadata_path  'Data/metadata.csv'  `--`gene_length_path  'Data/gene_lengths.csv'

Additionally, the `--`num-executors and `--`executor-cores arguments can be used to specify the number of worker nodes used and the number of cores per node used for execution. An example is given below.
> $ spark-submit preprocessing_spark.py  --num-executors 2 --executor-cores 4  `--`raw_data_path  'Data/covid_counts.csv'  `--`metadata_path  'Data/metadata.csv'  `--`gene_length_path  'Data/gene_lengths.csv'

The execution time is printed to the console.

### Numba Preprocessing

**Reproducibility Information:** The shared memory Numba implementation of preprocessing was run on an m5.2xlarge AWS instance using the AMI **Ubuntu Server 20.04 LTS (HVM), SSD Volume Type**. The linux kernel version is 5.4.0-1038-aws. The instance has 8 vCPU's, 8 cores in total (i.e. 1 core per vCPU), 32 GiB of main memory, 32 K of L1d cache memory, 32 K of L1i cache memory, 256 K of L2 cache memory, and 46080 K of L3 cache memory. The CPU clock rate is 2.3 GHz. \\
The shared memory Numba implementation of preprocessing was also run on an m5.4xlarge AWS instance using the AMI **Ubuntu Server 20.04 LTS (HVM), SSD Volume Type**. The linux kernel version is 5.4.0-1038-aws. The instance has 16 vCPU's, 16 cores in total (i.e. 1 core per vCPU), 32 GiB of main memory, 32 K of L1d cache memory, 32 K of L1i cache memory, 256 K of L2 cache memory, and 46080 K of L3 cache memory. The CPU clock rate is 2.3 GHz. \\
By default, the m5.2xlarge and m5.4xlarge instances has 8 G of disk space. Given that m5 intances are back by EBS, the disk space can by dynamically resized. Resizing was not needed for this execution.

The Python version used is 2.7.17. The following dependencies are required and can be installed by running the following command.
> $ pip install -r requirements_numba.txt

Dependencies:
* pandas 
* numpy 
* argparse
* numba 

Numba preprocessing can be run using the following command where `--`raw_data_path specifies the path to the raw cells by genes expression matrix, `--`metadata_path specifies the path to the metadata file containing cell-specific metadata, and `--`gene_length_path specifies the path to the file containing the gene lengths.
> $ ./preprocessing_numba.py  `--`raw_data_path  'Data/covid_counts.csv'  `--`metadata_path  'Data/metadata.csv'  `--`gene_length_path  'Data/gene_lengths.csv'

The execution time is printed to the console.

### Cell Clustering with MPI
**Reproducibility Information:** The MPI cluster was set up with AWS, using 16 nodes, each with an Ubuntu Server 18.04 LTS (HVM) image, EBS General Purpose (SSD) Volume Type and the instance type m4.xlarge. This instance has 4 Intel(R) Xeon(R) CPU E5-2686 v4 @ 2.30GHz, 2 cores per socket with 2 thread per core, and a Xen hipervisor, 16 GB of disk space and cache storages of: L1d cache 32K, L1i cache 32K, L2 cache 256K, L3 cache 46080K. We are using a VPC to create the cluster and the connections made are private.

The MPI version being used in all nodes is the HYDRA build 3.3a2, which is using gcc -Wl,-Bsymbolic-functions -Wl,-z,relro. The Python version used is in all nodes is 2.7.17.

Dependencies (all nodes):
* pandas 
* numpy 
* mpi4py

The Python version used is 2.7.17. The dependencies can be installed by running the following command.
> $ pip install -r requirements_clustering_mpi.txt

To execute the parallel implementation of K-Means based cell clustering, you can use the command
> $ mpiexec -n NUM_NODES python clustering_parallel.py

where NUM_NODES is the number of nodes in the cluster that you want to utilize. The file named `clusters.csv` should be in the same directory as the `clustering_parallel.py` file.

### Differential Expression Analysis with MPI
**Reproducibility Information:** The MPI cluster was set up with AWS, using 16 nodes, each with an Ubuntu Server 18.04 LTS (HVM) image, EBS General Purpose (SSD) Volume Type and the instance type m4.xlarge. This instance has 4 Intel(R) Xeon(R) CPU E5-2686 v4 @ 2.30GHz, 2 cores per socket with 2 thread per core, and a Xen hipervisor, 16 GB of disk space and cache storages of: L1d cache 32K, L1i cache 32K, L2 cache 256K, L3 cache 46080K. We are using a VPC to create the cluster and the connections made are private.

The MPI version being used in all nodes is the HYDRA build 3.3a2, which is using gcc -Wl,-Bsymbolic-functions -Wl,-z,relro. The Python version used is in all nodes is 2.7.17.

Dependencies (all nodes):
* pandas 
* numpy 
* scipy
* mpi4py

The Python version used is 2.7.17. The dependencies can be installed by running the following command.
> $ pip install -r requirements_DE_mpi.txt

To execute the parallel implementation of the Differential Expression Analysis, you can use the command
> $ mpiexec -n NUM_NODES python DE_parallel.py

where NUM_NODES is the number of nodes in the cluster that you want to utilize. The files named `exp_data.csv`, `metadata.csv`, `clusters.csv` should be in the same directory as the `DE_parallel.py` file. To get the `exp_data.csv` and `clusters.csv` please run the sequential script `run_sequential.py` before the parallel implementation.


## Performance
---

The execution times for each module in the single-cell analysis pipeline using the sequential implementation are provided below.
* Preprocessing execution time - 3088.9105 s
* Clustering with k-means execution time - 295.7333 s
* tSNE visualization execution time - 113.5847 s
* Differential expression analysis execution time - 7.1552 s
* Data I/O - 495.8741 s

The total execution time, not taking into account data I/O is 3505.3837 seconds (roughly 58 minutes). The run time taking into account I/O is 4001.2578 seconds (67 minutes).

We treat tSNE visualization as an inherently sequential portion of the code. The proportion of the code that can be parallelized is *c* and the expression for the the theoretical speedup *S<sub>T</sub>(1,p)* as a function of the number of processors *p* as dictated by Amdahl's law is given below. 
<p align="center">
  <img src="https://user-images.githubusercontent.com/29682604/117676275-63e55600-b17b-11eb-9dad-5435fdb88fa9.png">
</p>

If taking into account I/O, the proportion of the code that can be parallelized is *c* and the expression for the the theoretical speedup *S<sub>T</sub>(1,p)* as a function of the number of processors *p* as dictated by Amdahl's law is given below. 
<p align="center">
  <img src="https://user-images.githubusercontent.com/29682604/117676459-91320400-b17b-11eb-82f6-cc8addbc4d27.png">
</p>

We see that when ignoring data I/O, we can get sizeable return in speedup up to roughly 80 processors, after which we begin to realize diminishing returns in performance. However, if data I/O is taken into account, diminishing returns are realized much more rapidly. After roughly 15 to 20 processors, the performance increases are supposedly marginal.

![Theoretical_Speedup_with_IO](https://user-images.githubusercontent.com/29682604/117601460-78463600-b11c-11eb-9756-c81ec3b3ee6f.png)

### Preprocessing Performance

<u>PySpark:</u>

The PySpark implementation for the preprocessing steps was not amenable for the types of operations being performed on the raw expression data. 
**Challenges** Spark (and by extension PySpark) requires data manipulation of a specific format. The preprocessing steps involve column removal, column-wise and row-wise operations, and transformations on the entire matrix. At certain points in preprocessing, the collection of data represented by the method `collect()` (e.g. column names or index names) is unavoidable. Collecting elements from an RDD (or PySpark Dataframe which is a wrapped RDD), causes data to be sent back to the master node. Repeated data transfer inhibits PySpark's ability to be used as a viable means of big data flow preprocessing in this context.

<u>Numba:</u>

The preprocessing implemented using Numba produced large speedups. The Numba preprocessing implementation was run on an AWS m5.2xlarge instance with 8 vCPU's with an execution time of 8.6836 seconds. The Numba preprocessing implementation was also run on an AWS m5.4xlarge instance with 16 vCPU's with an execution time of 6.4418 seconds. These execution times correspond to speedups of 356 and 480 respectively.
**Benefits over the Spark framework** The implementation using Numba used a shared memory framework, avoiding overheads associated with data transfer and communication across nodes. Furthermore, Numba automatically detects the given machines CPU availability and allocates tasks to cores using an optimal scheme. However, Numba comes with some drawbacks. For example, functions must work on numerics and cannot work with strings or boolean expressions. While this was not an impediment for this use case, this drawback might pose a problem when being used in another context. 


### Cell Clustering Performance

While K-Means Clustering is capable of being run on relatively large datasets, the parallelized implementation ran into constant memory issues. This is likely a coding problem and not an algorithmic problem, as small-scale examples were able to be evaluated accurately. All attempts to run clustering on the full 30000x20000 data resulted in an `EXIT CODE: 9` error, which indicated to us that the code was using too much RAM and had to be exited. However, as a proof of concept, multiple nodes were used to evaluate a smaller version of the dataset to determine what the speedup looks like. A plot of that speedup vs. number of nodes is shown:

![cluster_speedup](https://user-images.githubusercontent.com/70713520/117595918-0a940d00-b110-11eb-93e5-e60195ed3e18.png)

Interestingly, the speedup gains as the number of nodes increases is relatively minor, maximizing at 8 nodes. The lowest speedup was actually at 2 nodes, where overall time was slower than a single node. This is likely due to the introduction of overhead communication costs between the two nodes. Aside from the 2 node case, it's interesting to see that there are "peaks" in speedup that occur when the number of nodes is a power of 2. This is possibly just a coincidence, but it's also possible splitting amongst *2<sub>n</sub>* nodes is inherently efficient in MPI.

### Differential Expression Analysis Performance
For the parallel execution of the Differential Expression Analysis, a varying number of nodes were tried out, from a range of two nodes (1 main and 1 worker node) to 16 nodes (1 main and 15 worker nodes). The execution time was the following:

| Nodes      | Execution Time | Speed Up |
| ----------- | ----------- | -------- |
| 1      | 7.1552       | 1 |
| 2      | 7.1788       | 0.997 |
| 3      | 4.294       | 1.666 |
| 4      | 3.567       | 2.006 |
| 5      | 2.977       | 2.404 |
| 6      | 2.442       | 2.931 |
| 7      | 2.376       | 3.011 |
| 8      | 2.287       | 3.129 |
| 9      | 2.186       | 3.329 |
| 10      | 2.033       | 3.473 |
| 11      | 1.858       | 3.851 |
| 12      | 1.827       | 3.915 |
| 13      | 1.859       | 3.849 |
| 14      | 1.859       | 3.848 |
| 15      | 1.873       | 3.820 |
| 16      | 1.876       | 3.813 |

![Differential Expression Speedup](https://user-images.githubusercontent.com/5700807/117595278-03680180-b106-11eb-926d-f10423613ce2.png)

Looking at the performance, this was exactly the speedup that we expected to see. Given the implementation of the application, were each cluster data is sent to a separate node, we expected to see an increase in the speedup when increasing the number of noes until the number of clusters was equal to the number of worker nodes used. In this case, there were 10 clusters in the data, so the maximum speed up was achieved when there were 10 worker nodes, i.e. 11 total nodes, as each node is assigned a single cluster. When we had 2 nodes, we could see that the performance was actually worse than that of the sequential process, and this is due to the overhead cost of setting up the MPI framework and the data transfer overhead. We could also see that when the number of nodes is bigger than the number of clusters, the speed up doesn't increase and actually slightly goes down, because there are nodes that actually end up doing no work and the increase in execution time is due to the increased overhead time of setting up the nodes.   


### Performance Comparison
Of the five steps in our data pipeline (preprocessing, clustering, visualization, differential expression analysis, data I/O), we worked on solutions to improve runtime in three of these steps. Successful performance gains were made in two of these steps: preprocessing and differential expression analysis. Although the clustering stage shows speedup on the small scale, it was not able to be applied to the full dataset, so it's runtime improvements will not be taken into consideration for the overall case. 

The overall speedup for each stage from sequential to parallel can be observed in the following table:

| Step in Data Pipeline      | Sequential Time (s) | Parallel Time (s) | Speed Up |
| ----------- | ----------- | -------- | ------ |
| Preprocessing  | 3088.9 | 5.44 | 567.81 |
| Clustering     | 295.7| 295.7 | 1.00 |
| Visualization  | 113.5| 113.5 | 1.00  |
| Differential Expression Analysis  |7.15 | 1.87 | 3.82 |
| Data I/O       | 495.8| 495.8 | 1.00 |

These individual speedups resulted in an overall speedup to the system of 4.39.

## Takeaways and Future Work
---

There were three major lessons learned throughout this experiment. Firstly, the lowest level operations need to be taken into consideration when parallelizing a task. It was not initially known how we were going to parallelize K-Means clustering, so in-depth analysis of the algorithm was required. Secondly, the transition from sequential code to parallel code is never as trivial as it may seem. Despite the many useful libraries and programming frameworks for moving sequential code to parallel code, it again requires knowing the granular details about the operations that are being made to the data. Where and how to apply these frameworks becomes complicated very quickly. Lastly, it has been critical to consider where and how data is stored, and how that data is accessed. Lots of time was spent moving the data to the compute instead of the other way around, which resulted in collective hours of waiting for data to upload.

Changes that would be made to this experiment in the future are focused mainly around the clustering implementation. The clustering stage is the last step in this pipeline that can be readily changed by parallelization, so fixing the implementation such that it is compatible with the large-scale dataset would be an important change. We would also like to increase the speedup achieved on the preprocessing step. Numba allows both CPU and GPU acceleration. While we focused on CPU acceleration for our implementation, utiliziing GPU's could better preprocessing performance. The overhead associated with data transfer between CPU and GPU would likely be mitigated by the performance increase on compute intensize operations across the expression matrix, but actual implementations across multiple problem sizes are needed in order to verify this claim.

## References
---
[1] Zheng Y, Liu X, Le W, et al. A human circulating immune cell landscape in aging and COVID-19. Protein & Cell. 2020 Oct;11(10):740-770. DOI: 10.1007/s13238-020-00762-2.
