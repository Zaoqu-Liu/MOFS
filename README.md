# MOFS: Multimodality Fusion Subtyping <img src="man/logo.jpg" alt="logo" align="right" height="200" width="180"/>

<!-- badges: start -->

[![Version](https://img.shields.io/badge/Version-1.0.0-orange)](https://github.com/Zaoqu-Liu/MOFS/)
[![Package](https://img.shields.io/badge/R_Package-MOFSR-blue)](https://github.com/Zaoqu-Liu/MOFS/MOFSR/MOFSR)
[![RCMD-Check](https://img.shields.io/badge/R%20CMD%20Check-yellow)](https://github.com/Zaoqu-Liu/MOFS/MOFSR)
[![Hits](https://hits.seeyoufarm.com/api/count/incr/badge.svg?url=https%3A%2F%2Fgithub.com%2FZaoqu-Liu%2FMOFS&count_bg=%2379C83D&title_bg=%23555555&icon=&icon_color=%23E7E7E7&title=hits&edge_flat=false)](https://hits.seeyoufarm.com)
[![License](https://img.shields.io/badge/License-GPL3-red)](https://github.com/Zaoqu-Liu/MOFS?tab=GPL-3.0-1-ov-file)
[![Manual](https://img.shields.io/badge/Manual-User-2a9d8f)](https://github.com/Zaoqu-Liu/MOFS/blob/master/MOFSR%20User%20Manual.pdf)
[![ShinyApp](https://img.shields.io/badge/Shiny-APP-f28482)](https://github.com/Zaoqu-Liu/MOFS/tree/master/MOFS-ShinyApp)
[![RCMD-Check](https://img.shields.io/badge/Feedback-c77dff)](liuzaoqu@163.com)
<!-- badges: end -->

The **MOFS (Multimodality Fusion Subtyping)** framework is a comprehensive approach for integrating and analyzing multi-layer biological data to achieve clinically relevant disease subtype classification. The core principle of MOFS lies in combining diverse biological modalities—such as genomic, transcriptomic, proteomic, histopathological, and radiological data—into a unified analytical framework. By integrating these heterogeneous data sources, MOFS offers a more holistic and nuanced view of disease biology, enabling the identification of underlying patterns that might be overlooked in single-modality analyses.

Integrating multimodal data reveals causal features that may be obscured in single-modality analyses, offering a more complete understanding of diseases. Multimodal data fusion can be categorized into early, intermediate, and late fusion based on the timing of integration. **Intermediate fusion** integrates data during clustering, allowing for the identification of multimodal joint clusters, capturing dependencies between different omics layers, and revealing underlying biological mechanisms. Intermediate fusion is generally considered more advanced than both early and late fusion, but it demands more sophisticated integration algorithms.

In the MOFS framework, multiple intermediate fusion strategies are employed based on diverse principles to maximize the interpretative power and reliability of the derived subtypes. Following this, late fusion is applied to combine the results obtained from the different algorithms, culminating in a final, consensus-driven clustering outcome. It is important to note that the framework does not rely on any specific intermediate fusion algorithm. Rather, it evaluates the overall consensus from multiple clustering results, ensuring robustness and reliability.

The MOFS framework represents a powerful tool for disease subtype classification by integrating multiple biological data modalities. Through the integration of intermediate fusion and late fusion strategies, MOFS uncovers novel disease subtypes and mechanisms that may be missed in single-modality approaches. This integrated approach not only enhances the accuracy and robustness of disease classification but also provides deeper insights into disease biology, potentially guiding more personalized diagnostic and therapeutic strategies.


<img src="man/abstract.jpg" width="100%" />


# MOFSR tool is designed for multimodal data fusion and analysis

### Feature Selection and Evaluation
1. **Statistical Metrics for Feature Evaluation**: The package provides several functions for calculating statistical metrics on data frames, which can be used to evaluate features. `MAD.df` calculates the median absolute deviation for each column in a data frame, `SD.df` calculates the standard deviation, `CV` calculates the coefficient of variation for a numeric vector, and `CV.df` calculates the coefficient of variation for each column in a data frame. These metrics can be used to assess the variability and distribution of features, which can be helpful in feature selection, especially when selecting hypervariable features. For example, in the `Select.Features` function, these statistical measures could potentially be used to identify features with high variance (e.g., based on standard deviation or coefficient of variation) or other characteristics relevant to the specific feature selection criteria. This allows for a more informed selection of features within a single modality or across multiple modalities, contributing to the overall goal of identifying the most relevant and informative features for further analysis, such as clustering or classification tasks within the MOFS framework.
2. **Optimal Feature Combination**: The `Find.OptClusterFeatures` plays a crucial role in this aspect. It aims to identify the most suitable combination of features for multi-modality clustering analysis. By integrating diverse modalities (such as mutation, CNV, RNA, protein, pathology, radiology), it explores the optimal number of clusters. This function utilizes Nonnegative Matrix Factorization (NMF) for cluster consistency analysis and Multi-block Principal Component Analysis (mbPCA) for assessing cluster separation. The combination of these methods not only helps in discerning distinct biological subgroups but also evaluates the stability and biological relevance of the clustering. In the process, it calculates metrics like the CPI (Cluster Performance Index) and GAP (Gap statistic) scores for all tested feature combinations and cluster numbers. The CPI and GAP scores provide quantitative measures to evaluate the quality of clustering solutions. A higher CPI or GAP score indicates a better clustering performance, and thus, the function returns a list containing the optimal feature combination along with the associated clustering score (including CPI and GAP scores) and all the results of the tested combinations, enabling users to make informed decisions about the feature selection.

### Multimodal Data Integration and Analysis
1. **Clustering Algorithm Integration**: It offers a variety of clustering algorithms, such as Consensus Iterative Multi-view Learning (CIMLR), Consensus Principal Component Analysis (CPCA), Integrative Non-negative Matrix Factorization (IntNMF), Low-Rank Approximation Clustering (LRAcluster), Multiple Co-Inertia Analysis (MCIA), Multimodality Fusion Subtyping (MOFS), NEMO, PINSPlus, Regularized Generalized Canonical Correlation Analysis (RGCCA), Sparse Generalized Canonical Correlation Analysis (SGCCA), Similarity Network Fusion (SNF), etc. These algorithms are used to integrate different modalities of data, like RNA, protein, genomic data, to identify shared patterns and structures in the data, thereby discovering potential disease subtypes or biological mechanisms.
2. **Comprehensive Evaluation of Multiple Algorithms**: The MOFS framework employs multiple intermediate fusion strategies and then uses late fusion to combine the results of different algorithms, ultimately obtaining a consensus-driven clustering result.

### Cluster Quality Evaluation
1. **Calinski-Harabasz Index with CalCHI**: The `CalCHI` function calculates the Calinski-Harabasz index for hierarchical clustering. It assesses cluster quality by comparing between-cluster and within-cluster dispersion. Higher values suggest better-defined clusters. It takes a clustering result and optional distance matrix, returning index values for different cluster numbers to help find the optimal number.
2. **PAC Calculation via CalPAC**: `CalPAC` calculates the Proportion of Ambiguous Clustering (PAC) for consensus clustering. A lower PAC indicates a more stable clustering solution. It uses consensus results and related parameters, returning PAC values for different cluster numbers.
3. **PCA Analysis by RunPCA**: `RunPCA` performs PCA on data. It helps visualize sample relationships in reduced dimensions, revealing data structure and patterns. Using the `FactoMineR` package, it provides insights and can assist in understanding clusters.
4. **Silhouette Coefficient**: The silhouette coefficient could be used to measure how well each data point lies within its cluster compared to other clusters. It provides a value between -1 and 1, where a higher value indicates better clustering. This metric can help evaluate the compactness and separation of clusters, complementing the other cluster quality evaluation methods.

### Classification Prediction Function
1. **Support for Multiple Classifiers**: It includes a variety of classifiers, such as AdaBoost, Decision Tree (DT), Elastic Net (Enet), Enrichment-based Neural Network (Enrichment), Gradient Boosted Decision Trees (GBDT), k-Nearest Neighbors (kNN), LASSO, Linear Discriminant Analysis (LDA), Naive Bayes (NBayes), Neural Network (NNet), PCA-based Neural Network (PCA), Random Forest (RF), Ridge Regression (Ridge), Stepwise Logistic Regression (StepLR), Support Vector Machine (SVM), XGBoost, etc. These classifiers can be used to predict cluster assignments for test data (single-modality data) based on trained models (single-modality data) and cluster markers.
2. **Ensemble Learning**: The `RunEnsemble` function runs an ensemble of different classification models, considering the consensus among the models, and can also perform survival analysis to filter models based on trends in clinical outcomes, providing more comprehensive and accurate prediction results.

### Functional Enrichment
1. **Gene Set Variation Analysis with RunGSVA**: The `RunGSVA` function estimates gene-set enrichment scores across samples. It offers methods like "gsva" for detecting subtle pathway activity changes, "ssgsea" for individual sample analysis, "zscore", and "plage". With parameters for expression data, gene sets, and method-specific options, it returns a gene-set by sample matrix of enrichment scores. This enables analysis of differentially enriched gene sets between groups, providing insights into perturbed biological functions and pathways.
2. **Single-Sample Pathway Activity Analysis using ssMwwGST**: The `ssMwwGST` function performs single-sample pathway activity analysis via the MWW-GST method. It calculates gene means and standard deviations, normalizes expression, assesses pathway enrichment, and corrects p-values. Taking gene expression data and gene sets as inputs, it returns matrices of NES, p-values, and FDR-adjusted p-values. These results help prioritize biologically relevant pathways, aiding in understanding underlying mechanisms and functions related to phenotypes or conditions.


## Installation
```R
if (!requireNamespace("devtools", quietly = TRUE)) {
  install.packages("devtools")
}
devtools::install_github('Zaoqu-Liu/MOFS', subdir = 'MOFSR')
```

## Multi-Modality Feature Selection

This R script performs multi-modality clustering by integrating diverse biological data types, such as mutation, CNV, RNA and protein expression, pathology, and radiology. Its main goal is to determine the optimal number of clusters in a multi-modality dataset using Nonnegative Matrix Factorization (NMF) and Principal Component Analysis (PCA). These approaches are effective in identifying distinct biological subgroups within complex datasets, aiding in the understanding of disease mechanisms and discovery of new therapeutic targets.

The script combines different feature subset sizes for each modality, iteratively evaluating the clustering performance. Two metrics are used to assess clustering quality: Cluster Prediction Index (CPI) and GAP statistic. CPI, calculated using the IntNMF package, measures the stability of clustering by quantifying how consistently features and samples group across multiple iterations. Higher CPI values indicate more stable and robust clustering solutions. GAP statistic, calculated using the mogsa package, helps identify the optimal number of clusters by comparing the dispersion of the actual dataset with that of a random reference dataset. The peak value of the GAP statistic indicates the most meaningful clustering structure.

Feature selection is done using Median Absolute Deviation (MAD) to rank features by variability, with higher MAD values suggesting greater potential for cluster differentiation. The script employs parallel computing with the future.apply package to accelerate the analysis, considering numerous combinations of feature subset sizes and cluster numbers.

```R
# A list of multi-modality datasets
data_list <- list(
  mutation,  # Mutation data (placeholder)
  cnv,       # Copy Number Variation data (placeholder)
  rna,       # RNA expression data
  protein,   # Protein expression data
  pathology, # Pathology data
  radiology  # Radiology imaging data
)

# Create a dataframe of all possible combinations of feature subset sizes for modality data
mrna_size <- seq(3000, 8000, 500)  # Sequence of mRNA feature subset sizes to be tested
protein_size <- seq(500, 2000, 500)  # Sequence of protein feature subset sizes to be tested
pathology_size <- seq(200, 600, 200)  # Sequence of pathology feature subset sizes to be tested
radiology_size <- seq(500, 2000, 500)  # Sequence of radiology feature subset sizes to be tested
feature_combinations <- expand.grid(mrna_size, protein_size, pathology_size, radiology_size)  # Create all possible combinations

# Calculate variance of each feature to determine their significance
# Use Median Absolute Deviation (MAD) to rank features for each data type
rna_data <- data_list$rna
mRNA_ID <- names(sort(apply(rna_data, 1, mad), decreasing = TRUE))  # Rank RNA features by MAD

protein_data <- data_list$protein
Protein_ID <- names(sort(apply(protein_data, 1, mad), decreasing = TRUE))  # Rank protein features by MAD

pathology_data <- data_list$pathology
Pathology_ID <- names(sort(apply(pathology_data, 1, mad), decreasing = TRUE))  # Rank pathology features by MAD

radiology_data <- data_list$radiology
Radiology_ID <- names(sort(apply(radiology_data, 1, mad), decreasing = TRUE))  # Rank radiology features by MAD

# Set up parallel computing
library(future.apply)
plan('multisession')  # Use multiple sessions for parallel computing

# Define the range of cluster numbers to be tested
try_num_clusters <- 2:6  # Range of clusters to try during clustering

# Loop over each row of the feature combination dataframe using parallel execution
results <- future_lapply(1:nrow(feature_combinations), function(row_n) {
  Sys.sleep(0.01)  # Short pause to simulate progress
  progress(row_n, nrow(feature_combinations))  # Show progress
  
  # Extract feature subset sizes for the current row
  subset_sizes <- as.numeric(feature_combinations[row_n, ])
  
  # Create a working copy of the data
  working_data <- data_list
  working_data$rna <- working_data$rna[mRNA_ID[1:subset_sizes[1]], ]  # Subset RNA data based on top-ranked features
  working_data$protein <- working_data$protein[Protein_ID[1:subset_sizes[2]], ]  # Subset protein data based on top-ranked features
  working_data$pathology <- working_data$pathology[Pathology_ID[1:subset_sizes[3]], ]  # Subset pathology data based on top-ranked features
  working_data$radiology <- working_data$radiology[Radiology_ID[1:subset_sizes[4]], ]  # Subset radiology data based on top-ranked features
  
  # Normalize each data type
  normalized_data <- lapply(working_data, function(dataset) {
    if (!all(dataset >= 0)) {
      dataset <- pmax(dataset + abs(min(dataset)), 0) + .Machine$double.eps  # Ensure non-negativity by shifting values
    }
    dataset <- dataset / max(dataset)  # Scale to max value of 1
    return(as.matrix(dataset))
  })
  
  # Transpose data to have features in columns (necessary for downstream analysis)
  normalized_data <- lapply(normalized_data, function(x) t(x) + .Machine$double.eps)
  
  # Perform optimal cluster number selection using NMF
  # CPI (Cluster Prediction Index) quantifies the consistency of clustering results across multiple iterations
  # CPI is calculated by evaluating the stability of clustering memberships across resampled data subsets
  opt_k_nmf <- IntNMF::nmf.opt.k(
    dat = normalized_data,        # Input data
    n.runs = 5,                   # Number of times NMF should run for each cluster number
    n.fold = 5,                   # Number of folds for cross-validation
    k.range = try_num_clusters,   # Range of clusters to try
    result = TRUE,                # Show results
    make.plot = FALSE,            # Do not plot intermediate results
    maxiter = 1000,               # Maximum number of iterations
    st.count = 10,                # Stop count threshold
    progress = FALSE              # Disable internal progress
  )
  cpi_df <- as.data.frame(opt_k_nmf)  # Convert the result to a dataframe
  cpi_df$mean <- rowMeans(cpi_df)     # Calculate mean CPI across iterations (higher CPI indicates better clustering consistency)
  
  # Perform multi-block PCA (mbPCA) for integration analysis
  # mbPCA is used to extract global components that summarize the shared variation across different modality datasets
  mbpca_result <- mogsa::mbpca(
    x = working_data,              # Input multi-modality data
    ncomp = 3,                     # Number of principal components to compute
    k = 0.5,                       # Proportion of variables to include
    method = "globalScore",       # Integration scoring method
    option = "uniform",           # Standardize across datasets
    center = TRUE,                 # Center data
    scale = TRUE,                  # Scale data
    moa = TRUE,                    # Return as moa object
    svd.solver = "fast",          # Use fast SVD solver
    maxiter = 1000,                # Maximum number of iterations
    verbose = FALSE                # Suppress output
  )
  # Perform clustering gap analysis using hierarchical clustering
  # GAP statistic measures the robustness of clustering by comparing intra-cluster dispersion to that of a reference distribution
  # The GAP statistic helps to determine the optimal number of clusters by identifying the point where the GAP value reaches its maximum
  gap_analysis <- mogsa::moGap(mbpca_result, K.max = max(try_num_clusters), cluster = "hclust", plot = FALSE)
  gap_df <- as.data.frame(gap_analysis$Tab)[-1, ]
  
  # Collect results for the current combination
  tmp <- data.frame(
    mRNA_num = subset_sizes[1],  # Number of RNA features used
    Protein_num = subset_sizes[2],  # Number of protein features used
    Pathology_num = subset_sizes[3],  # Number of pathology features used
    Radiology_num = subset_sizes[4],  # Number of radiology features used
    K = try_num_clusters,  # Number of clusters tested
    CPI = cpi_df$mean,  # Mean CPI value
    GAP = gap_df$gap  # GAP statistic value for cluster separation
  )
  return(tmp)
})

# Combine results from all iterations into one dataframe
combined_results <- Reduce(rbind, results)
# Calculate a combined score (sum of CPI and GAP)
# Higher score indicates better clustering consistency and separation
combined_results$score <- combined_results$CPI + combined_results$GAP
```

### Multi-Modality Clustering

This R script performs clustering on a multi-omics dataset using various clustering algorithms and methods, allowing for the exploration of the optimal number of clusters. The analysis employs both intermediate fusion and late-stage fusion approaches to achieve robust clustering results. The intermediate fusion strategy involves combining information across different omics types before running clustering algorithms, whereas the late-stage fusion approach clusters individual data types first and then combines their results. Each approach offers distinct advantages, such as improved interpretability and robustness in detecting cluster-specific features.

The initial steps involve using a variety of clustering algorithms on the multi-omics data. Then, an intermediate fusion approach is applied, allowing the fusion of different data types before applying clustering methods. This fusion approach allows for better integration and consistent interpretation of complex biological data. Afterward, the Jaccard distance is calculated for clustering consensus analysis, followed by consensus clustering (COCA) to derive the final clustering solution.

The script further evaluates the clustering performance using silhouette scores, PCA visualization, and clustering metrics such as the Proportion of Ambiguous Clustering (PAC) score and Calinski-Harabasz index (CHI). The final result is stored in a list, which contains all intermediate and final clustering outcomes, providing a comprehensive summary of the multi-omics clustering process.
```R
# Initialize an empty list to store the clustering results
# This part runs multiple clustering algorithms separately, storing results for each
res <- list()
res[['CPCA']] <- RunCPCA(data_list, cluster_num)  # Run CPCA clustering
res[['CIMLR']] <- RunCIMLR(data_list, cluster_num)  # Run CIMLR clustering
res[['iClusterBayes']] <- RuniClusterBayes(data_list, cluster_num)  # Run iClusterBayes clustering
res[['IntNMF']] <- RunIntNMF(data_list, cluster_num)  # Run IntNMF clustering
res[['LRAcluster']] <- RunLRAcluster(data_list, cluster_num)  # Run LRAcluster clustering
res[['MCIA']] <- RunMCIA(data_list, cluster_num)  # Run MCIA clustering
res[['NEMO']] <- RunNEMO(data_list, cluster_num)  # Run NEMO clustering
res[['PINSPlus']] <- RunPINSPlus(data_list, cluster_num)  # Run PINSPlus clustering
res[['RGCCA']] <- RunRGCCA(data_list, cluster_num)  # Run RGCCA clustering
res[['SGCCA']] <- RunSGCCA(data_list, cluster_num)  # Run SGCCA clustering
res[['SNF']] <- RunSNF(data_list, cluster_num)  # Run SNF clustering
```

#### Alternatively, you can use an intermediate fusion step by using RunIF

```R
# This allows specifying a method in the 'RunIF' function
res <- list()
res[['CPCA']] <- RunIF(data = data_list, method = 'CPCA', cluster_num)  # Run CPCA with intermediate fusion
res[['CIMLR']] <- RunIF(data = data_list, method = 'CIMLR', cluster_num)  # Run CIMLR with intermediate fusion
res[['iClusterBayes']] <- RunIF(data = data_list, method = 'iClusterBayes', cluster_num)  # Run iClusterBayes with intermediate fusion
res[['IntNMF']] <- RunIF(data = data_list, method = 'IntNMF', cluster_num)  # Run IntNMF with intermediate fusion
res[['LRAcluster']] <- RunIF(data = data_list, method = 'LRAcluster', cluster_num)  # Run LRAcluster with intermediate fusion
res[['MCIA']] <- RunIF(data = data_list, method = 'MCIA', cluster_num)  # Run MCIA with intermediate fusion
res[['NEMO']] <- RunIF(data = data_list, method = 'NEMO', cluster_num)  # Run NEMO with intermediate fusion
res[['PINSPlus']] <- RunIF(data = data_list, method = 'PINSPlus', cluster_num)  # Run PINSPlus with intermediate fusion
res[['RGCCA']] <- RunIF(data = data_list, method = 'RGCCA', cluster_num)  # Run RGCCA with intermediate fusion
res[['SGCCA']] <- RunIF(data = data_list, method = 'SGCCA', cluster_num)  # Run SGCCA with intermediate fusion
res[['SNF']] <- RunIF(data = data_list, method = 'SNF', cluster_num)  # Run SNF with intermediate fusion
```

```R
# Compute binary cluster membership matrix and Jaccard distance
bm <- get.binary.clusters(res)  # Obtain binary cluster membership matrix from results
sm <- get.Jaccard.Distance(bm)  # Calculate Jaccard distance based on the membership matrix

# Run consensus clustering (COCA) on the Jaccard distance matrix
coca_res <- RunCOCA(jaccard.matrix = sm,  # Input Jaccard distance matrix
                    max.clusters = 6,  # Maximum number of clusters to consider
                    linkage.method = "ward.D2",  # Linkage method for hierarchical clustering
                    clustering.algorithm = 'pam',  # Clustering algorithm to use
                    distance.metric = "euclidean",  # Distance metric to use for clustering
                    resampling.iterations = 10000,  # Number of iterations for resampling
                    resample.proportion = 0.7)  # Proportion of samples to resample in each iteration

# Calculate cluster metrics for evaluating clustering performance
PAC <- CalPAC(coca_res$fit)  # Calculate Proportion of Ambiguous Clustering (PAC) score
CHI <- CalCHI(hclust(as.dist(sm), method = "average"), max_clusters = 6)  # Calculate Calinski-Harabasz index for cluster evaluation

# Calculate consensus matrix and silhouette score
consMatrix <- 1 - coca_res$optimal$consensusMatrix  # Calculate consensus matrix (1 - optimal consensus matrix from COCA)
rownames(consMatrix) <- colnames(consMatrix) <- colnames(sm)  # Set row and column names for the consensus matrix
consMatrix <- as.dist(consMatrix)  # Convert consensus matrix to distance matrix format
tmp <- coca_res$Cluster$Cluster  # Extract cluster assignments from COCA results
names(tmp) <- colnames(sm)  # Set cluster assignment names to sample names
aSil <- cluster::silhouette(x = tmp, dist = consMatrix)  # Calculate silhouette scores for each sample

# Plot silhouette scores
plot(aSil, col = c('red','blue','yellow'), main = '')  # Plot silhouette scores with specific colors
abline(v = 0.5)  # Add a vertical line at silhouette score of 0.5 to indicate threshold

# Identify core samples based on silhouette scores
colnames(sm)[aSil[, 3] < 0.4]  # Identify samples with silhouette scores less than 0.4 (not well-clustered)
core_set <- colnames(sm)[aSil[, 3] >= 0.4]  # Define core set as samples with silhouette scores greater or equal to 0.4
sm_core <- sm[core_set, core_set]  # Subset Jaccard distance matrix for core samples

# Extract cluster assignment for the core set
cluster <- get.class(coca_res$fit, 3)  # Get class assignments from COCA results for 3 clusters
cluster <- cluster[match(core_set, cluster$ID), ]  # Match cluster assignment to core set

# Run PCA on the core similarity matrix
ddb.pca <- RunPCA(sm_core)  # Perform Principal Component Analysis (PCA) on core Jaccard similarity matrix

# Visualize PCA results
factoextra::fviz_pca_ind(
  X = ddb.pca,  # Input PCA result
  geom.ind = "point",  # Use point geometry for individuals
  pointshape = 21,  # Set point shape
  fill.ind = cluster$Cluster,  # Fill points based on cluster assignment
  palette = "npg",  # Use a color palette
  alpha.ind = 0.7,  # Set point transparency
  addEllipses = TRUE  # Add ellipses for each cluster
) +
  theme_classic() +
  theme(plot.title = element_text(hjust = 0.5))  # Set classic theme and center the title

# Store all the results in a list and save to a file
allres <- list(All_Cluster = resm,  # All cluster results
               JaccardDis = sm,  # Jaccard distance matrix
               coca_res = coca_res,  # Consensus clustering results
               Silhouette = aSil,  # Silhouette scores
               Core_set = core_set,  # Core set of well-clustered samples
               PCA = ddb.pca,  # PCA results
               PAC = PAC,  # PAC score
               Calinsky = aCalinsky)  # Calinski-Harabasz index
```

#### Alternatively, you can use the RunMOFS function to get all the above results directly

```R
allres <- RunMOFS(data = data_list,  # Input multi-modality dataset
                  methods = c("CPCA", "iClusterBayes", "IntNMF", "LRAcluster", "MCIA", "NEMO", "PINSPlus", "RGCCA", "SGCCA", "SNF", "CIMLR"),  # Methods to run
                  max.clusters = 6,  # Maximum number of clusters to consider
                  linkage.method = "ward.D2",  # Linkage method for clustering
                  clustering.algorithm = "pam",  # Clustering algorithm to use
                  distance.metric = "euclidean",  # Distance metric to use
                  resampling.iterations = 10000,  # Number of resampling iterations
                  resample.proportion = 0.7,  # Proportion of samples to resample
                  silhouette.cutoff = 0.4)  # Silhouette cutoff value for core sample selection
```

<p align="center">
<a href="https://clustrmaps.com/site/1c32n"  title="Visit tracker"><img src="//www.clustrmaps.com/map_v2.png?d=_vogv3DGZKPfHE6Dm7QUcVRJ6R6O-VuMWV1JVQGCmAM&cl=ffffff" /></a>
</p>


|Package                   |Title                                                                                                                                                             |Version     |License                                     |
|:-------------------------|:-----------------------------------------------------------------------------------------------------------------------------------------------------------------|:-----------|:-------------------------------------------|
|MOFSR                     |Multimodality Fusion Subtyping                                                                                                                                    |1.0.0       |GPL (>= 3)                                  |
|ComplexHeatmap            |Make Complex Heatmaps                                                                                                                                             |2.20.0      |MIT + file LICENSE                          |
|circlize                  |Circular Visualization                                                                                                                                            |0.4.16      |MIT + file LICENSE                          |
|furrr                     |Apply Mapping Functions in Parallel using Futures                                                                                                                 |0.3.1       |MIT + file LICENSE                          |
|future                    |Unified Parallel and Distributed Processing in R for Everyone                                                                                                     |1.33.2      |LGPL (>= 2.1)                               |
|neuralnet                 |Training of Neural Networks                                                                                                                                       |1.44.2      |GPL (>= 2)                                  |
|ggpol                     |Visualizing Social Science Data with 'ggplot2'                                                                                                                    |0.0.7       |MIT + file LICENSE                          |
|NeuralNetTools            |Visualization and Analysis Tools for Neural Networks                                                                                                              |1.5.3       |CC0                                         |
|ggsankey                  |Sankey, Alluvial and Sankey Bump Plots                                                                                                                            |0.0.99999   |MIT + file LICENSE                          |
|mlr3verse                 |Easily Install and Load the 'mlr3' Package Family                                                                                                                 |0.3.0       |LGPL-3                                      |
|mlr3                      |Machine Learning in R - Next Generation                                                                                                                           |0.21.1      |LGPL-3                                      |
|randomForest              |Breiman and Cutlers Random Forests for Classification and
Regression                                                                                               |4.7-1.2     |GPL (>= 2)                                  |
|psych                     |Procedures for Psychological, Psychometric, and Personality
Research                                                                                               |2.4.3       |GPL (>= 2)                                  |
|e1071                     |Misc Functions of the Department of Statistics, Probability
Theory Group (Formerly: E1071), TU Wien                                                                |1.7-14      |GPL-2 &#124; GPL-3                          |
|gbm                       |Generalized Boosted Regression Models                                                                                                                             |2.1.9       |GPL (>= 2) &#124; file LICENSE              |
|doMC                      |Foreach Parallel Adaptor for 'parallel'                                                                                                                           |1.3.8       |GPL-2                                       |
|yaGST                     |Competitive gene set and regulon tests.                                                                                                                           |2017.08.25  |GPL (>= 3)                                  |
|nnet                      |Feed-Forward Neural Networks and Multinomial Log-Linear Models                                                                                                    |7.3-19      |GPL-2 &#124; GPL-3                          |
|glmnet                    |Lasso and Elastic-Net Regularized Generalized Linear Models                                                                                                       |4.1-8       |GPL-2                                       |
|Matrix                    |Sparse and Dense Matrix Classes and Methods                                                                                                                       |1.7-0       |GPL (>= 2) &#124; file LICENCE              |
|adabag                    |Applies Multiclass AdaBoost.M1, SAMME and Bagging                                                                                                                 |5.0         |GPL (>= 2)                                  |
|doParallel                |Foreach Parallel Adaptor for the 'parallel' Package                                                                                                               |1.0.17      |GPL-2                                       |
|iterators                 |Provides Iterator Construct                                                                                                                                       |1.0.14      |Apache License (== 2.0)                     |
|foreach                   |Provides Foreach Looping Construct                                                                                                                                |1.5.2       |Apache License (== 2.0)                     |
|caret                     |Classification and Regression Training                                                                                                                            |6.0-94      |GPL (>= 2)                                  |
|lattice                   |Trellis Graphics for R                                                                                                                                            |0.22-6      |GPL (>= 2)                                  |
|rpart                     |Recursive Partitioning and Regression Trees                                                                                                                       |4.1.23      |GPL-2 &#124; GPL-3                          |
|maftools                  |Summarize, Analyze and Visualize MAF Files                                                                                                                        |2.20.0      |MIT + file LICENSE                          |
|patchwork                 |The Composer of Plots                                                                                                                                             |1.2.0       |MIT + file LICENSE                          |
|survminer                 |Drawing Survival Curves using 'ggplot2'                                                                                                                           |0.4.9       |GPL-2                                       |
|ggpubr                    |'ggplot2' Based Publication Ready Plots                                                                                                                           |0.6.0       |GPL (>= 2)                                  |
|survival                  |Survival Analysis                                                                                                                                                 |3.6-4       |LGPL (>= 2)                                 |
|fstcore                   |R Bindings to the 'Fstlib' Library                                                                                                                                |0.9.18      |MPL-2.0 &#124; file LICENSE                 |
|Biobase                   |Biobase: Base functions for Bioconductor                                                                                                                          |2.64.0      |Artistic-2.0                                |
|BiocGenerics              |S4 generic functions used in Bioconductor                                                                                                                         |0.50.0      |Artistic-2.0                                |
|DynDoc                    |Dynamic document tools                                                                                                                                            |1.80.0      |Artistic-2.0                                |
|widgetTools               |Creates an interactive tcltk widget                                                                                                                               |1.82.0      |LGPL                                        |
|MCMCpack                  |Markov Chain Monte Carlo (MCMC) Package                                                                                                                           |1.7-1       |GPL-3                                       |
|coda                      |Output Analysis and Diagnostics for MCMC                                                                                                                          |0.19-4.1    |GPL (>= 2)                                  |
|RColorBrewer              |ColorBrewer Palettes                                                                                                                                              |1.1-3       |Apache License 2.0                          |
|ROCR                      |Visualizing the Performance of Scoring Classifiers                                                                                                                |1.0-11      |GPL (>= 2)                                  |
|verification              |Weather Forecast Verification Utilities                                                                                                                           |1.42        |GPL (>= 2)                                  |
|dtw                       |Dynamic Time Warping Algorithms                                                                                                                                   |1.23-1      |GPL (>= 2)                                  |
|proxy                     |Distance and Similarity Measures                                                                                                                                  |0.4-27      |GPL-2                                       |
|CircStats                 |Circular Statistics, from "Topics in Circular Statistics" (2001)                                                                                                  |0.2-6       |GPL-2                                       |
|MASS                      |Support Functions and Datasets for Venables and Ripley's MASS                                                                                                     |7.3-60.2    |GPL-2 &#124; GPL-3                          |
|boot                      |Bootstrap Functions (Originally by Angelo Canty for S)                                                                                                            |1.3-30      |Unlimited                                   |
|fields                    |Tools for Spatial Data                                                                                                                                            |15.2        |GPL (>= 2)                                  |
|viridisLite               |Colorblind-Friendly Color Maps (Lite Version)                                                                                                                     |0.4.2       |MIT + file LICENSE                          |
|spam                      |SPArse Matrix                                                                                                                                                     |2.10-0      |LGPL-2 &#124; BSD_3_clause + file LICENSE   |
|gtools                    |Various R Programming Tools                                                                                                                                       |3.9.5       |GPL-2                                       |
|lubridate                 |Make Dealing with Dates a Little Easier                                                                                                                           |1.9.3       |GPL (>= 2)                                  |
|forcats                   |Tools for Working with Categorical Variables (Factors)                                                                                                            |1.0.0       |MIT + file LICENSE                          |
|stringr                   |Simple, Consistent Wrappers for Common String Operations                                                                                                          |1.5.1       |MIT + file LICENSE                          |
|dplyr                     |A Grammar of Data Manipulation                                                                                                                                    |1.1.4       |MIT + file LICENSE                          |
|purrr                     |Functional Programming Tools                                                                                                                                      |1.0.2       |MIT + file LICENSE                          |
|readr                     |Read Rectangular Text Data                                                                                                                                        |2.1.5       |MIT + file LICENSE                          |
|tidyr                     |Tidy Messy Data                                                                                                                                                   |1.3.1       |MIT + file LICENSE                          |
|tibble                    |Simple Data Frames                                                                                                                                                |3.2.1       |MIT + file LICENSE                          |
|ggplot2                   |Create Elegant Data Visualisations Using the Grammar of Graphics                                                                                                  |3.5.1       |MIT + file LICENSE                          |
|tidyverse                 |Easily Install and Load the 'Tidyverse'                                                                                                                           |2.0.0       |MIT + file LICENSE                          |
|ssgsea.GBM.classification |What the package does (short line)                                                                                                                                |1.0         |What license is it under?                   |
|pheatmap                  |Pretty Heatmaps                                                                                                                                                   |1.0.12      |GPL-2                                       |
|DBI                       |R Database Interface                                                                                                                                              |1.2.2       |LGPL (>= 2.1)                               |
|bslib                     |Custom 'Bootstrap' 'Sass' Themes for 'shiny' and 'rmarkdown'                                                                                                      |0.7.0       |MIT + file LICENSE                          |
|lpSolve                   |Interface to 'Lp_solve' v. 5.5 to Solve Linear/Integer Programs                                                                                                   |5.6.20      |LGPL-2                                      |
|httr                      |Tools for Working with URLs and HTTP                                                                                                                              |1.4.7       |MIT + file LICENSE                          |
|registry                  |Infrastructure for R Package Registries                                                                                                                           |0.5-1       |GPL-2                                       |
|CoGAPS                    |Coordinated Gene Activity in Pattern Sets                                                                                                                         |3.25.2      |BSD_3_clause + file LICENSE                 |
|BiocParallel              |Bioconductor facilities for parallel evaluation                                                                                                                   |1.38.0      |GPL-2 &#124; GPL-3                          |
|prettyunits               |Pretty, Human Readable Formatting of Quantities                                                                                                                   |1.2.0       |MIT + file LICENSE                          |
|mlr3cluster               |Cluster Extension for 'mlr3'                                                                                                                                      |0.1.10      |LGPL-3                                      |
|yulab.utils               |Supporting Functions for Packages Maintained by 'YuLab-SMU'                                                                                                       |0.1.4       |Artistic-2.0                                |
|ggplotify                 |Convert Plot to 'grob' or 'ggplot' Object                                                                                                                         |0.1.2       |Artistic-2.0                                |
|sparseMatrixStats         |Summary Statistics for Rows and Columns of Sparse Matrices                                                                                                        |1.16.0      |MIT + file LICENSE                          |
|spatstat.geom             |Geometrical Functionality of the 'spatstat' Family                                                                                                                |3.3-2       |GPL (>= 2)                                  |
|survMisc                  |Miscellaneous Functions for Survival Data                                                                                                                         |0.5.6       |GPL-2                                       |
|pillar                    |Coloured Formatting for Columns                                                                                                                                   |1.9.0       |MIT + file LICENSE                          |
|Rgraphviz                 |Provides plotting capabilities for R graph objects                                                                                                                |2.48.0      |EPL                                         |
|R6                        |Encapsulated Classes with Reference Semantics                                                                                                                     |2.5.1       |MIT + file LICENSE                          |
|mime                      |Map Filenames to MIME Types                                                                                                                                       |0.12        |GPL                                         |
|reticulate                |Interface to 'Python'                                                                                                                                             |1.37.0      |Apache License 2.0                          |
|uwot                      |The Uniform Manifold Approximation and Projection (UMAP) Method
for Dimensionality Reduction                                                                       |0.2.2       |GPL (>= 3)                                  |
|gridtext                  |Improved Text Rendering Support for 'Grid' Graphics                                                                                                               |0.1.5       |MIT + file LICENSE                          |
|Rttf2pt1                  |'ttf2pt1' Program                                                                                                                                                 |1.3.12      |file LICENSE                                |
|viridis                   |Colorblind-Friendly Color Maps for R                                                                                                                              |0.6.5       |MIT + file LICENSE                          |
|genekitr                  |Gene Analysis Toolkit                                                                                                                                             |1.2.5       |GPL-3                                       |
|Rhdf5lib                  |hdf5 library as an R package                                                                                                                                      |1.26.0      |Artistic-2.0                                |
|dorothea                  |Collection Of Human And Mouse TF Regulons                                                                                                                         |1.16.0      |GPL-3 + file LICENSE                        |
|Hmisc                     |Harrell Miscellaneous                                                                                                                                             |5.1-2       |GPL (>= 2)                                  |
|fpc                       |Flexible Procedures for Clustering                                                                                                                                |2.2-13      |GPL                                         |
|KMsurv                    |Data sets from Klein and Moeschberger (1997), Survival Analysis                                                                                                   |0.1-5       |GPL (>= 3)                                  |
|parallelly                |Enhancing the 'parallel' Package                                                                                                                                  |1.37.1      |LGPL (>= 2.1)                               |
|SpatialExperiment         |S4 Class for Spatially Resolved -omics Data                                                                                                                       |1.14.0      |GPL-3                                       |
|GlobalOptions             |Generate Functions to Get or Set Global Options                                                                                                                   |0.1.2       |MIT + file LICENSE                          |
|caTools                   |Tools: Moving Window Statistics, GIF, Base64, ROC AUC, etc                                                                                                        |1.18.2      |GPL-3                                       |
|FNN                       |Fast Nearest Neighbor Search Algorithms and Applications                                                                                                          |1.1.4       |GPL (>= 2)                                  |
|polyclip                  |Polygon Clipping                                                                                                                                                  |1.10-6      |BSL                                         |
|NMF                       |Algorithms and Framework for Nonnegative Matrix Factorization
(NMF)                                                                                                |0.30.4.900  |GPL (>=2)                                   |
|beachmat                  |Compiling Bioconductor to Handle Each Matrix Type                                                                                                                 |2.20.0      |GPL-3                                       |
|htmltools                 |Tools for HTML                                                                                                                                                    |0.5.8.1     |GPL (>= 2)                                  |
|fansi                     |ANSI Control Sequence Aware String Functions                                                                                                                      |1.0.6       |GPL-2 &#124; GPL-3                          |
|remotes                   |R Package Installation from Remote Repositories, Including
'GitHub'                                                                                                |2.5.0       |MIT + file LICENSE                          |
|commonmark                |High Performance CommonMark and Github Markdown Rendering in R                                                                                                    |1.9.1       |BSD_2_clause + file LICENSE                 |
|ggrepel                   |Automatically Position Non-Overlapping Text Labels with
'ggplot2'                                                                                                  |0.9.5       |GPL-3 &#124; file LICENSE                   |
|car                       |Companion to Applied Regression                                                                                                                                   |3.1-2       |GPL (>= 2)                                  |
|fgsea                     |Fast Gene Set Enrichment Analysis                                                                                                                                 |1.30.0      |MIT + file LICENCE                          |
|spatstat.utils            |Utility Functions for 'spatstat'                                                                                                                                  |3.1-0       |GPL (>= 2)                                  |
|HDO.db                    |A set of annotation maps describing the entire Human Disease
Ontology                                                                                              |0.99.1      |Artistic-2.0                                |
|clusterProfiler           |A universal enrichment tool for interpreting omics data                                                                                                           |4.13.0      |Artistic-2.0                                |
|maps                      |Draw Geographical Maps                                                                                                                                            |3.4.2       |GPL-2                                       |
|scatterplot3d             |3D Scatter Plot                                                                                                                                                   |0.3-44      |GPL-2                                       |
|mlr3learners              |Recommended Learners for 'mlr3'                                                                                                                                   |0.8.0       |LGPL-3                                      |
|clue                      |Cluster Ensembles                                                                                                                                                 |0.3-65      |GPL-2                                       |
|goftest                   |Classical Goodness-of-Fit Tests for Univariate Distributions                                                                                                      |1.2-3       |GPL (>= 2)                                  |
|fitdistrplus              |Help to Fit of a Parametric Distribution to Non-Censored or
Censored Data                                                                                          |1.2-1       |GPL (>= 2)                                  |
|scatterpie                |Scatter Pie Plot                                                                                                                                                  |0.2.2       |Artistic-2.0                                |
|tidyselect                |Select from a Set of Strings                                                                                                                                      |1.2.1       |MIT + file LICENSE                          |
|scattermore               |Scatterplots with More Points                                                                                                                                     |1.2         |GPL (>= 3)                                  |
|RSQLite                   |SQLite Interface for R                                                                                                                                            |2.3.6       |LGPL (>= 2.1)                               |
|cowplot                   |Streamlined Plot Theme and Plot Annotations for 'ggplot2'                                                                                                         |1.1.3       |GPL-2                                       |
|GenomeInfoDbData          |Species and taxonomy ID look up tables used by GenomeInfoDb                                                                                                       |1.2.12      |Artistic-2.0                                |
|utf8                      |Unicode Text Processing                                                                                                                                           |1.2.4       |Apache License (== 2.0) &#124; file LICENSE |
|ScaledMatrix              |Creating a DelayedMatrix of Scaled and Centered Values                                                                                                            |1.12.0      |GPL-3                                       |
|sessioninfo               |R Session Information                                                                                                                                             |1.2.2       |GPL-2                                       |
|spatstat.data             |Datasets for 'spatstat' Family                                                                                                                                    |3.1-2       |GPL (>= 2)                                  |
|gridExtra                 |Miscellaneous Functions for "Grid" Graphics                                                                                                                       |2.3         |GPL (>= 2)                                  |
|fs                        |Cross-Platform File System Operations Based on 'libuv'                                                                                                            |1.6.4       |MIT + file LICENSE                          |
|xgboost                   |Extreme Gradient Boosting                                                                                                                                         |1.7.8.1     |Apache License (== 2.0) &#124; file LICENSE |
|sctransform               |Variance Stabilizing Transformations for Single Cell UMI Data                                                                                                     |0.4.1       |GPL-3 &#124; file LICENSE                   |
|future.apply              |Apply Function to Elements in Parallel using Futures                                                                                                              |1.11.2      |GPL (>= 2)                                  |
|ggVennDiagram             |A 'ggplot2' Implement of Venn Diagram                                                                                                                             |1.5.2       |GPL-3                                       |
|graph                     |graph: A package to handle graph data structures                                                                                                                  |1.82.0      |Artistic-2.0                                |
|ipred                     |Improved Predictors                                                                                                                                               |0.9-14      |GPL (>= 2)                                  |
|mlr3hyperband             |Hyperband for 'mlr3'                                                                                                                                              |0.6.0       |LGPL-3                                      |
|uuid                      |Tools for Generating and Handling of UUIDs                                                                                                                        |1.2-0       |MIT + file LICENSE                          |
|bbotk                     |Black-Box Optimization Toolkit                                                                                                                                    |1.3.0       |LGPL-3                                      |
|BioEnricher               |Bioinformatics analysis and visualization pipeline                                                                                                                |4.0.3       |MIT + file LICENSE                          |
|Rtsne                     |T-Distributed Stochastic Neighbor Embedding using a Barnes-Hut
Implementation                                                                                      |0.17        |file LICENSE                                |
|lazyeval                  |Lazy (Non-Standard) Evaluation                                                                                                                                    |0.2.2       |GPL-3                                       |
|sass                      |Syntactically Awesome Style Sheets ('Sass')                                                                                                                       |0.4.9       |MIT + file LICENSE                          |
|scales                    |Scale Functions for Visualization                                                                                                                                 |1.3.0       |MIT + file LICENSE                          |
|carData                   |Companion to Applied Regression Data Sets                                                                                                                         |3.0-5       |GPL (>= 2)                                  |
|munsell                   |Utilities for Using Munsell Colours                                                                                                                               |0.5.1       |MIT + file LICENSE                          |
|lava                      |Latent Variable Models                                                                                                                                            |1.8.0       |GPL-3                                       |
|treeio                    |Base Classes and Functions for Phylogenetic Tree Input and
Output                                                                                                  |1.28.0      |Artistic-2.0                                |
|profvis                   |Interactive Visualizations for Profiling R Code                                                                                                                   |0.3.8       |GPL-3 &#124; file LICENSE                   |
|KEGGgraph                 |KEGGgraph: A graph approach to KEGG PATHWAY in R and
Bioconductor                                                                                                  |1.64.0      |GPL (>= 2)                                  |
|bitops                    |Bitwise Operations                                                                                                                                                |1.0-7       |GPL (>= 2)                                  |
|labeling                  |Axis Labeling                                                                                                                                                     |0.4.3       |MIT + file LICENSE &#124; Unlimited         |
|KEGGREST                  |Client-side REST access to the Kyoto Encyclopedia of Genes and
Genomes (KEGG)                                                                                      |1.44.0      |Artistic-2.0                                |
|promises                  |Abstractions for Promise-Based Asynchronous Programming                                                                                                           |1.3.0       |MIT + file LICENSE                          |
|SCOPE                     |Single Cell Omics Pipeline and Environment                                                                                                                        |1.0.0       |GPL (>= 3)                                  |
|shape                     |Functions for Plotting Graphical Shapes, Colors                                                                                                                   |1.4.6.1     |GPL (>= 3)                                  |
|rhdf5filters              |HDF5 Compression Filters                                                                                                                                          |1.16.0      |BSD_2_clause + file LICENSE                 |
|lgr                       |A Fully Featured Logging Framework                                                                                                                                |0.4.4       |MIT + file LICENSE                          |
|zoo                       |S3 Infrastructure for Regular and Irregular Time Series (Z's
Ordered Observations)                                                                                 |1.8-12      |GPL-2 &#124; GPL-3                          |
|GeneNMF                   |Non-Negative Matrix Factorization for Single-Cell Omics                                                                                                           |0.6.2       |GPL-3                                       |
|GSVA                      |Gene Set Variation Analysis for Microarray and RNA-Seq Data                                                                                                       |1.52.0      |GPL (>= 2)                                  |
|fst                       |Lightning Fast Serialization of Data Frames                                                                                                                       |0.9.8       |AGPL-3 &#124; file LICENSE                  |
|DelayedArray              |A unified framework for working transparently with on-disk and
in-memory array-like datasets                                                                       |0.30.0      |Artistic-2.0                                |
|RSpectra                  |Solvers for Large-Scale Eigenvalue and SVD Problems                                                                                                               |0.16-1      |MPL (>= 2)                                  |
|multcomp                  |Simultaneous Inference in General Parametric Models                                                                                                               |1.4-25      |GPL-2                                       |
|assertthat                |Easy Pre and Post Assertions                                                                                                                                      |0.2.1       |GPL-3                                       |
|tools                     |Tools for Package Development                                                                                                                                     |4.4.0       |Part of R 4.4.0                             |
|ape                       |Analyses of Phylogenetics and Evolution                                                                                                                           |5.8         |GPL-2 &#124; GPL-3                          |
|shiny                     |Web Application Framework for R                                                                                                                                   |1.8.1.1     |GPL-3 &#124; file LICENSE                   |
|SingleCellExperiment      |S4 Classes for Single Cell Data                                                                                                                                   |1.26.0      |GPL-3                                       |
|mlr3filters               |Filter Based Feature Selection for 'mlr3'                                                                                                                         |0.8.1       |LGPL-3                                      |
|rlang                     |Functions for Base Types and Core R and 'Tidyverse' Features                                                                                                      |1.1.3       |MIT + file LICENSE                          |
|SCEVAN                    |SCEVAN                                                                                                                                                            |1.0.1       |GPL-2                                       |
|generics                  |Common S3 Generics not Provided by Base R Methods Related to
Model Fitting                                                                                         |0.1.3       |MIT + file LICENSE                          |
|diceR                     |Diverse Cluster Ensemble in R                                                                                                                                     |2.2.0       |MIT + file LICENSE                          |
|ggridges                  |Ridgeline Plots in 'ggplot2'                                                                                                                                      |0.5.6       |GPL-2 &#124; file LICENSE                   |
|BiocSingular              |Singular Value Decomposition for Bioconductor Packages                                                                                                            |1.20.0      |GPL-3                                       |
|markdown                  |Render Markdown with 'commonmark'                                                                                                                                 |1.12        |MIT + file LICENSE                          |
|extrafont                 |Tools for Using Fonts                                                                                                                                             |0.19        |GPL-2                                       |
|evaluate                  |Parsing and Evaluation Tools that Provide More Details than the
Default                                                                                            |0.23        |MIT + file LICENSE                          |
|GenomeInfoDb              |Utilities for manipulating chromosome names, including modifying
them to follow a particular naming style                                                          |1.40.0      |Artistic-2.0                                |
|fastcluster               |Fast Hierarchical Clustering Routines for R and 'Python'                                                                                                          |1.2.6       |FreeBSD &#124; GPL-2 &#124; file LICENSE    |
|spacesXYZ                 |CIE XYZ and some of Its Derived Color Spaces                                                                                                                      |1.3-0       |GPL (>= 3)                                  |
|reshape2                  |Flexibly Reshape Data: A Reboot of the Reshape Package                                                                                                            |1.4.4       |MIT + file LICENSE                          |
|devtools                  |Tools to Make Developing R Packages Easier                                                                                                                        |2.4.5       |MIT + file LICENSE                          |
|bcellViper                |Human B-cell transcriptional interactome and normal human B-cell
expression data                                                                                   |1.40.0      |GPL (>=2)                                   |
|gelnet                    |Generalized Elastic Nets                                                                                                                                          |1.2.1       |GPL (>= 3)                                  |
|ConsRank                  |Compute the Median Ranking(s) According to the Kemeny's
Axiomatic Approach                                                                                         |2.1.4       |GPL-3                                       |
|geneset                   |Get Gene Sets for Gene Enrichment Analysis                                                                                                                        |0.2.7       |GPL-3                                       |
|colorspace                |A Toolbox for Manipulating and Assessing Colors and Palettes                                                                                                      |2.1-0       |BSD_3_clause + file LICENSE                 |
|limSolve                  |Solving Linear Inverse Models                                                                                                                                     |1.5.7.1     |GPL                                         |
|ellipsis                  |Tools for Working with ...                                                                                                                                        |0.3.2       |MIT + file LICENSE                          |
|data.table                |Extension of `data.frame`                                                                                                                                         |1.15.4      |MPL-2.0 &#124; file LICENSE                 |
|withr                     |Run Code 'With' Temporarily Modified Global State                                                                                                                 |3.0.0       |MIT + file LICENSE                          |
|presto                    |Fast Functions for Differential Expression using Wilcox and AUC                                                                                                   |1.0.0       |GPL-3                                       |
|RCurl                     |General Network (HTTP/FTP/...) Client Interface for R                                                                                                             |1.98-1.14   |BSD_3_clause + file LICENSE                 |
|mlr3data                  |Collection of Machine Learning Data Sets for 'mlr3'                                                                                                               |0.9.0       |LGPL-3                                      |
|xtable                    |Export Tables to LaTeX or HTML                                                                                                                                    |1.8-4       |GPL (>= 2)                                  |
|plyr                      |Tools for Splitting, Applying and Combining Data                                                                                                                  |1.8.9       |MIT + file LICENSE                          |
|paradox                   |Define and Work with Parameter Spaces for Complex Algorithms                                                                                                      |1.0.1       |LGPL-3                                      |
|aplot                     |Decorate a 'ggplot' with Associated Information                                                                                                                   |0.2.2       |Artistic-2.0                                |
|systemfonts               |System Native Font Finding                                                                                                                                        |1.1.0       |MIT + file LICENSE                          |
|MatrixModels              |Modelling with Sparse and Dense Matrices                                                                                                                          |0.5-3       |GPL (>= 2)                                  |
|ggvenn                    |Draw Venn Diagram by 'ggplot2'                                                                                                                                    |0.1.10      |MIT + file LICENSE                          |
|mclust                    |Gaussian Mixture Modelling for Model-Based Clustering,
Classification, and Density Estimation                                                                      |6.1.1       |GPL (>= 2)                                  |
|httpuv                    |HTTP and WebSocket Server Library                                                                                                                                 |1.6.15      |GPL (>= 2) &#124; file LICENSE              |
|rmarkdown                 |Dynamic Documents for R                                                                                                                                           |2.27        |GPL-3                                       |
|robustbase                |Basic Robust Statistics                                                                                                                                           |0.99-4      |GPL (>= 2)                                  |
|openxlsx                  |Read, Write and Edit xlsx Files                                                                                                                                   |4.2.5.2     |MIT + file LICENSE                          |
|broom                     |Convert Statistical Objects into Tidy Tibbles                                                                                                                     |1.0.6       |MIT + file LICENSE                          |
|FactoMineR                |Multivariate Exploratory Data Analysis and Data Mining                                                                                                            |2.11        |GPL (>= 2)                                  |
|deldir                    |Delaunay Triangulation and Dirichlet (Voronoi) Tessellation                                                                                                       |2.0-4       |GPL (>= 2)                                  |
|GO.db                     |A set of annotation maps describing the entire Gene Ontology                                                                                                      |3.19.1      |Artistic-2.0                                |
|sandwich                  |Robust Covariance Matrix Estimators                                                                                                                               |3.1-0       |GPL-2 &#124; GPL-3                          |
|rhdf5                     |R Interface to HDF5                                                                                                                                               |2.48.0      |Artistic-2.0                                |
|tensor                    |Tensor product of arrays                                                                                                                                          |1.5         |GPL (>= 2)                                  |
|ragg                      |Graphic Devices Based on AGG                                                                                                                                      |1.3.2       |MIT + file LICENSE                          |
|vctrs                     |Vector Helpers                                                                                                                                                    |0.6.5       |MIT + file LICENSE                          |
|lifecycle                 |Manage the Life Cycle of your Package Functions                                                                                                                   |1.0.4       |MIT + file LICENSE                          |
|codetools                 |Code Analysis Tools for R                                                                                                                                         |0.2-20      |GPL                                         |
|DT                        |A Wrapper of the JavaScript Library 'DataTables'                                                                                                                  |0.33        |GPL-3 &#124; file LICENSE                   |
|mnormt                    |The Multivariate Normal and t Distributions, and Their Truncated
Versions                                                                                          |2.1.1       |GPL-2 &#124; GPL-3                          |
|recipes                   |Preprocessing and Feature Engineering Steps for Modeling                                                                                                          |1.0.10      |MIT + file LICENSE                          |
|nlme                      |Linear and Nonlinear Mixed Effects Models                                                                                                                         |3.1-164     |GPL (>= 2)                                  |
|mcmc                      |Markov Chain Monte Carlo                                                                                                                                          |0.9-8       |MIT + file LICENSE                          |
|mlr3tuningspaces          |Search Spaces for 'mlr3'                                                                                                                                          |0.5.1       |LGPL-3                                      |
|progress                  |Terminal Progress Bars                                                                                                                                            |1.2.3       |MIT + file LICENSE                          |
|pkgload                   |Simulate Package Installation and Attach                                                                                                                          |1.3.4       |GPL-3                                       |
|jquerylib                 |Obtain 'jQuery' as an HTML Dependency Object                                                                                                                      |0.1.4       |MIT + file LICENSE                          |
|Rcpp                      |Seamless R and C++ Integration                                                                                                                                    |1.0.13      |GPL (>= 2)                                  |
|rstudioapi                |Safely Access the RStudio API                                                                                                                                     |0.16.0      |MIT + file LICENSE                          |
|stringi                   |Fast and Portable Character String Processing Facilities                                                                                                          |1.8.4       |file LICENSE                                |
|pbapply                   |Adding Progress Bar to '*apply' Functions                                                                                                                         |1.7-2       |GPL (>= 2)                                  |
|hms                       |Pretty Time of Day                                                                                                                                                |1.1.3       |MIT + file LICENSE                          |
|cachem                    |Cache R Objects with Automatic Pruning                                                                                                                            |1.1.0       |MIT + file LICENSE                          |
|mlr3mbo                   |Flexible Bayesian Optimization                                                                                                                                    |0.2.7       |LGPL-3                                      |
|pROC                      |Display and Analyze ROC Curves                                                                                                                                    |1.18.5      |GPL (>= 3)                                  |
|BiocManager               |Access the Bioconductor Project Package Repository                                                                                                                |1.30.23     |Artistic-2.0                                |
|tidytree                  |A Tidy Tool for Phylogenetic Tree Data Manipulation                                                                                                               |0.4.6       |Artistic-2.0                                |
|listenv                   |Environments Behaving (Almost) as Lists                                                                                                                           |0.9.1       |LGPL (>= 2.1)                               |
|XVector                   |Foundation of external vector representation and manipulation in
Bioconductor                                                                                      |0.44.0      |Artistic-2.0                                |
|plotly                    |Create Interactive Web Graphics via 'plotly.js'                                                                                                                   |4.10.4      |MIT + file LICENSE                          |
|urlchecker                |Run CRAN URL Checks from Older R Versions                                                                                                                         |1.0.1       |GPL-3                                       |
|WGCNA                     |Weighted Correlation Network Analysis                                                                                                                             |1.72-5      |GPL (>= 2)                                  |
|ggtree                    |an R package for visualization of tree and annotation data                                                                                                        |3.12.0      |Artistic-2.0                                |
|enrichplot                |Visualization of Functional Enrichment Result                                                                                                                     |1.24.0      |Artistic-2.0                                |
|GetoptLong                |Parsing Command-Line Arguments and Simple Variable Interpolation                                                                                                  |1.0.5       |MIT + file LICENSE                          |
|pkgbuild                  |Find Tools Needed to Build R Packages                                                                                                                             |1.4.4       |MIT + file LICENSE                          |
|palmerpenguins            |Palmer Archipelago (Antarctica) Penguin Data                                                                                                                      |0.1.1       |CC0                                         |
|ggfun                     |Miscellaneous Functions for 'ggplot2'                                                                                                                             |0.1.4       |Artistic-2.0                                |
|HDF5Array                 |HDF5 backend for DelayedArray objects                                                                                                                             |1.32.0      |Artistic-2.0                                |
|ggsci                     |Scientific Journal and Sci-Fi Themed Color Palettes for
'ggplot2'                                                                                                  |3.1.0       |GPL (>= 3)                                  |
|SparseArray               |High-performance sparse data representation and manipulation in
R                                                                                                  |1.4.1       |Artistic-2.0                                |
|htmlwidgets               |HTML Widgets for R                                                                                                                                                |1.6.4       |MIT + file LICENSE                          |
|kernlab                   |Kernel-Based Machine Learning Lab                                                                                                                                 |0.9-32      |GPL-2                                       |
|Formula                   |Extended Model Formulas                                                                                                                                           |1.2-5       |GPL-2 &#124; GPL-3                          |
|prabclus                  |Functions for Clustering and Testing of Presence-Absence,
Abundance and Multilocus Genetic Data                                                                    |2.3-4       |GPL                                         |
|leaps                     |Regression Subset Selection                                                                                                                                       |3.1         |GPL (>= 2)                                  |
|dendextend                |Extending 'dendrogram' Functionality in R                                                                                                                         |1.17.1      |GPL-2 &#124; GPL-3                          |
|class                     |Functions for Classification                                                                                                                                      |7.3-22      |GPL-2 &#124; GPL-3                          |
|memoise                   |'Memoisation' of Functions                                                                                                                                        |2.0.1       |MIT + file LICENSE                          |
|crayon                    |Colored Terminal Output                                                                                                                                           |1.5.2       |MIT + file LICENSE                          |
|gridGraphics              |Redraw Base Graphics Using 'grid' Graphics                                                                                                                        |0.5-1       |GPL (>= 2)                                  |
|Seurat                    |Tools for Single Cell Genomics                                                                                                                                    |4.4.0       |MIT + file LICENSE                          |
|mlr3pipelines             |Preprocessing Operators and Pipelines for 'mlr3'                                                                                                                  |0.7.1       |LGPL-3                                      |
|S4Arrays                  |Foundation of array-like containers in Bioconductor                                                                                                               |1.4.0       |Artistic-2.0                                |
|xml2                      |Parse XML                                                                                                                                                         |1.3.6       |MIT + file LICENSE                          |
|preprocessCore            |A collection of pre-processing functions                                                                                                                          |1.66.0      |LGPL (>= 2)                                 |
|GOSemSim                  |GO-terms Semantic Similarity Measures                                                                                                                             |2.30.0      |Artistic-2.0                                |
|flexmix                   |Flexible Mixture Modeling                                                                                                                                         |2.3-19      |GPL (>= 2)                                  |
|ggtext                    |Improved Text Rendering Support for 'ggplot2'                                                                                                                     |0.1.2       |GPL-2                                       |
|UCSC.utils                |Low-level utilities to retrieve data from the UCSC Genome
Browser                                                                                                  |1.0.0       |Artistic-2.0                                |
|png                       |Read and write PNG images                                                                                                                                         |0.1-8       |GPL-2 &#124; GPL-3                          |
|progressr                 |An Inclusive, Unifying API for Progress Updates                                                                                                                   |0.14.0      |GPL (>= 3)                                  |
|tzdb                      |Time Zone Database Information                                                                                                                                    |0.4.0       |MIT + file LICENSE                          |
|emmeans                   |Estimated Marginal Means, aka Least-Squares Means                                                                                                                 |1.10.2      |GPL-2 &#124; GPL-3                          |
|fastmap                   |Fast Data Structures                                                                                                                                              |1.2.0       |MIT + file LICENSE                          |
|GSEABase                  |Gene set enrichment data structures and methods                                                                                                                   |1.66.0      |Artistic-2.0                                |
|tidygraph                 |A Tidy API for Graph Manipulation                                                                                                                                 |1.3.1       |MIT + file LICENSE                          |
|flashClust                |Implementation of optimal hierarchical clustering                                                                                                                 |1.01-2      |GPL (>= 2)                                  |
|pkgconfig                 |Private Configuration for 'R' Packages                                                                                                                            |2.0.3       |MIT + file LICENSE                          |
|cli                       |Helpers for Developing Command Line Interfaces                                                                                                                    |3.6.2       |MIT + file LICENSE                          |
|DOSE                      |Disease Ontology Semantic and Enrichment analysis                                                                                                                 |3.30.1      |Artistic-2.0                                |
|ggforce                   |Accelerating 'ggplot2'                                                                                                                                            |0.4.2       |MIT + file LICENSE                          |
|europepmc                 |R Interface to the Europe PubMed Central RESTful Web Service                                                                                                      |0.4.3       |GPL-3                                       |
|prodlim                   |Product-Limit Estimation for Censored Event History Analysis                                                                                                      |2023.08.28  |GPL (>= 2)                                  |
|pathview                  |a tool set for pathway based data integration and visualization                                                                                                   |1.44.0      |GPL (>=3.0)                                 |
|ggsignif                  |Significance Brackets for 'ggplot2'                                                                                                                               |0.6.4       |GPL-3 &#124; file LICENSE                   |
|gridBase                  |Integration of base and grid graphics                                                                                                                             |0.4-7       |GPL                                         |
|rlist                     |A Toolbox for Non-Tabular Data Manipulation                                                                                                                       |0.4.6.2     |MIT + file LICENSE                          |
|SummarizedExperiment      |SummarizedExperiment container                                                                                                                                    |1.34.0      |Artistic-2.0                                |
|ggalluvial                |Alluvial Plots in 'ggplot2'                                                                                                                                       |0.12.5      |GPL-3                                       |
|lmtest                    |Testing Linear Regression Models                                                                                                                                  |0.9-40      |GPL-2 &#124; GPL-3                          |
|textshaping               |Bindings to the 'HarfBuzz' and 'Fribidi' Libraries for Text
Shaping                                                                                                |0.3.7       |MIT + file LICENSE                          |
|RcppAnnoy                 |'Rcpp' Bindings for 'Annoy', a Library for Approximate Nearest
Neighbors                                                                                           |0.0.22      |GPL (>= 2)                                  |
|usethis                   |Automate Package and Project Setup                                                                                                                                |2.2.3       |MIT + file LICENSE                          |
|multcompView              |Visualizations of Paired Comparisons                                                                                                                              |0.1-10      |GPL                                         |
|DEoptimR                  |Differential Evolution Optimization in Pure R                                                                                                                     |1.1-3       |GPL (>= 2)                                  |
|timechange                |Efficient Manipulation of Date-Times                                                                                                                              |0.3.0       |GPL (>= 3)                                  |
|foreign                   |Read Data Stored by 'Minitab', 'S', 'SAS', 'SPSS', 'Stata',
'Systat', 'Weka', 'dBase', ...                                                                         |0.8-86      |GPL (>= 2)                                  |
|timeDate                  |Rmetrics - Chronological and Calendar Objects                                                                                                                     |4032.109    |GPL (>= 2)                                  |
|splines                   |Regression Spline Functions and Classes                                                                                                                           |4.4.0       |Part of R 4.4.0                             |
|mlr3viz                   |Visualizations for 'mlr3'                                                                                                                                         |0.10.0      |LGPL-3                                      |
|askpass                   |Password Entry Utilities for R, Git, and SSH                                                                                                                      |1.2.0       |MIT + file LICENSE                          |
|blob                      |A Simple S3 Class for Representing Vectors of Binary Data
('BLOBS')                                                                                                |1.2.4       |MIT + file LICENSE                          |
|impute                    |impute: Imputation for microarray data                                                                                                                            |1.78.0      |GPL-2                                       |
|annotate                  |Annotation for microarrays                                                                                                                                        |1.82.0      |Artistic-2.0                                |
|XML                       |Tools for Parsing and Generating XML Within R and S-Plus                                                                                                          |3.99-0.16.1 |BSD_3_clause + file LICENSE                 |
|network                   |Classes for Relational Data                                                                                                                                       |1.18.2      |GPL (>= 2)                                  |
|globals                   |Identify Global Objects in R Expressions                                                                                                                          |0.16.3      |LGPL (>= 2.1)                               |
|knitr                     |A General-Purpose Package for Dynamic Report Generation in R                                                                                                      |1.46        |GPL                                         |
|ica                       |Independent Component Analysis                                                                                                                                    |1.0-3       |GPL (>= 2)                                  |
|stats4                    |Statistical Functions using S4 Classes                                                                                                                            |4.4.0       |Part of R 4.4.0                             |
|compiler                  |The R Compiler Package                                                                                                                                            |4.4.0       |Part of R 4.4.0                             |
|rjson                     |JSON for R                                                                                                                                                        |0.2.21      |GPL-2                                       |
|triebeard                 |'Radix' Trees in 'Rcpp'                                                                                                                                           |0.4.1       |MIT + file LICENSE                          |
|DNAcopy                   |DNA Copy Number Data Analysis                                                                                                                                     |1.78.0      |GPL (>= 2)                                  |
|pkgmaker                  |Development Utilities for R Packages                                                                                                                              |0.32.8      |GPL (>=2)                                   |
|mlr3fselect               |Feature Selection for 'mlr3'                                                                                                                                      |1.2.1       |LGPL-3                                      |
|extrafontdb               |Package for holding the database for the extrafont package                                                                                                        |1.0         |GPL-2                                       |
|bit                       |Classes and Methods for Fast Memory-Efficient Boolean Selections                                                                                                  |4.0.5       |GPL-2 &#124; GPL-3                          |
|RcppML                    |Rcpp Machine Learning Library                                                                                                                                     |0.3.7       |GPL (>= 2)                                  |
|diptest                   |Hartigan's Dip Test Statistic for Unimodality - Corrected                                                                                                         |0.77-1      |GPL (>= 2)                                  |
|BiocNeighbors             |Nearest Neighbor Detection for Bioconductor Packages                                                                                                              |1.20.2      |GPL-3                                       |
|lsa                       |Latent Semantic Analysis                                                                                                                                          |0.73.3      |GPL (>= 2)                                  |
|glue                      |Interpreted String Literals                                                                                                                                       |1.7.0       |MIT + file LICENSE                          |
|sp                        |Classes and Methods for Spatial Data                                                                                                                              |2.1-4       |GPL (>= 2)                                  |
|estimability              |Tools for Assessing Estimability of Linear Predictions                                                                                                            |1.5.1       |GPL (>= 3)                                  |
|mlr3tuning                |Hyperparameter Optimization for 'mlr3'                                                                                                                            |1.2.0       |LGPL-3                                      |
|ggnetwork                 |Geometries to Plot Networks with 'ggplot2'                                                                                                                        |0.5.13      |GPL-3                                       |
|digest                    |Create Compact Hash Digests of R Objects                                                                                                                          |0.6.35      |GPL (>= 2)                                  |
|leiden                    |R Implementation of Leiden Clustering Algorithm                                                                                                                   |0.4.3.1     |GPL-3 &#124; file LICENSE                   |
|quadprog                  |Functions to Solve Quadratic Programming Problems                                                                                                                 |1.5-8       |GPL (>= 2)                                  |
|urltools                  |Vectorised Tools for URL Handling and Parsing                                                                                                                     |1.7.3       |MIT + file LICENSE                          |
|irlba                     |Fast Truncated Singular Value Decomposition and Principal
Components Analysis for Large Dense and Sparse Matrices                                                  |2.3.5.1     |GPL-3                                       |
|spacefillr                |Space-Filling Random and Quasi-Random Sequences                                                                                                                   |0.3.3       |MIT + file LICENSE                          |
|graphlayouts              |Additional Layout Algorithms for Network Visualizations                                                                                                           |1.1.1       |MIT + file LICENSE                          |
|rgl                       |3D Visualization Using OpenGL                                                                                                                                     |1.3.12      |GPL                                         |
|magick                    |Advanced Graphics and Image-Processing in R                                                                                                                       |2.8.3       |MIT + file LICENSE                          |
|GenomicRanges             |Representation and manipulation of genomic intervals                                                                                                              |1.56.0      |Artistic-2.0                                |
|mlr3misc                  |Helper Functions for 'mlr3'                                                                                                                                       |0.15.1      |LGPL-3                                      |
|spatstat.random           |Random Generation Functionality for the 'spatstat' Family                                                                                                         |3.3-1       |GPL (>= 2)                                  |
|SparseM                   |Sparse Linear Algebra                                                                                                                                             |1.81        |GPL (>= 2)                                  |
|zlibbioc                  |An R packaged zlib-1.2.5                                                                                                                                          |1.50.0      |Artistic-2.0 + file LICENSE                 |
|tidydr                    |Unify Dimensionality Reduction Results                                                                                                                            |0.0.5.001   |Artistic-2.0                                |
|dotCall64                 |Enhanced Foreign Function Interface Supporting Long Vectors                                                                                                       |1.1-1       |GPL (>= 2)                                  |
|tweenr                    |Interpolate Data for Smooth Animations                                                                                                                            |2.0.3       |MIT + file LICENSE                          |
|CellChat                  |Inference and analysis of cell-cell communication from
single-cell and spatially resolved transcriptomics data                                                     |2.1.2       |GPL-3                                       |
|ModelMetrics              |Rapid Calculation of Model Metrics                                                                                                                                |1.2.2.2     |GPL (>= 2)                                  |
|microbenchmark            |Accurate Timing Functions                                                                                                                                         |1.4.10      |BSD_2_clause + file LICENSE                 |
|ggraph                    |An Implementation of Grammar of Graphics for Graphs and Networks                                                                                                  |2.2.1       |MIT + file LICENSE                          |
|openssl                   |Toolkit for Encryption, Signatures and Certificates Based on
OpenSSL                                                                                               |2.2.0       |MIT + file LICENSE                          |
|rsvd                      |Randomized Singular Value Decomposition                                                                                                                           |1.0.5       |GPL (>= 3)                                  |
|gson                      |Base Class and Methods for 'gson' Format                                                                                                                          |0.1.0       |Artistic-2.0                                |
|igraph                    |Network Analysis and Visualization                                                                                                                                |2.0.3       |GPL (>= 2)                                  |
|mvtnorm                   |Multivariate Normal and t Distributions                                                                                                                           |1.2-5       |GPL-2                                       |
|qvalue                    |Q-value estimation for false discovery rate control                                                                                                               |2.36.0      |LGPL                                        |
|later                     |Utilities for Scheduling Functions to Execute Later with Event
Loops                                                                                               |1.3.2       |MIT + file LICENSE                          |
|modeltools                |Tools and Classes for Statistical Models                                                                                                                          |0.2-23      |GPL-2                                       |
|statnet.common            |Common R Scripts and Utilities Used by the Statnet Project
Software                                                                                                |4.9.0       |GPL-3 + file LICENSE                        |
|backports                 |Reimplementations of Functions Introduced Since R-3.0.0                                                                                                           |1.5.0       |GPL-2 &#124; GPL-3                          |
|rstatix                   |Pipe-Friendly Framework for Basic Statistical Tests                                                                                                               |0.7.2       |GPL-2                                       |
|shadowtext                |Shadow Text Grob and Layer                                                                                                                                        |0.1.3       |Artistic-2.0                                |
|AnnotationDbi             |Manipulation of SQLite-based annotations in Bioconductor                                                                                                          |1.66.0      |Artistic-2.0                                |
|sna                       |Tools for Social Network Analysis                                                                                                                                 |2.7-2       |GPL (>= 2)                                  |
|Mfuzz                     |Soft clustering of time series gene expression data                                                                                                               |2.64.0      |GPL-2                                       |
|quantreg                  |Quantile Regression                                                                                                                                               |5.97        |GPL (>= 2)                                  |
|miniUI                    |Shiny UI Widgets for Small Screens                                                                                                                                |0.1.1.1     |GPL-3                                       |
|gtable                    |Arrange 'Grobs' in Tables                                                                                                                                         |0.3.5       |MIT + file LICENSE                          |
|abind                     |Combine Multidimensional Arrays                                                                                                                                   |1.4-5       |LGPL (>= 2)                                 |
|xfun                      |Supporting Functions for Packages Maintained by 'Yihui Xie'                                                                                                       |0.44        |MIT + file LICENSE                          |
|Cairo                     |R Graphics Device using Cairo Graphics Library for Creating
High-Quality Bitmap (PNG, JPEG, TIFF), Vector (PDF, SVG,
PostScript) and Display (X11 and Win32) Output |1.6-2       |GPL-2 &#124; GPL-3                          |
|Biostrings                |Efficient manipulation of biological strings                                                                                                                      |2.72.0      |Artistic-2.0                                |
|dynamicTreeCut            |Methods for Detection of Clusters in Hierarchical Clustering
Dendrograms                                                                                           |1.63-1      |GPL (>= 2)                                  |
|org.Hs.eg.db              |Genome wide annotation for Human                                                                                                                                  |3.19.1      |Artistic-2.0                                |
|KernSmooth                |Functions for Kernel Smoothing Supporting Wand & Jones (1995)                                                                                                     |2.23-24     |Unlimited                                   |
|jsonlite                  |A Simple and Robust JSON Parser and Generator for R                                                                                                               |1.8.8       |MIT + file LICENSE                          |
|magrittr                  |A Forward-Pipe Operator for R                                                                                                                                     |2.0.3       |MIT + file LICENSE                          |
|tkWidgets                 |R based tk widgets                                                                                                                                                |1.82.0      |Artistic-2.0                                |
|svglite                   |An 'SVG' Graphics Device                                                                                                                                          |2.1.3       |GPL (>= 2)                                  |
|base64enc                 |Tools for base64 encoding                                                                                                                                         |0.1-3       |GPL-2 &#124; GPL-3                          |
|spatstat.univar           |One-Dimensional Probability Distribution Support for the
'spatstat' Family                                                                                         |3.0-0       |GPL (>= 2)                                  |
|TH.data                   |TH's Data Archive                                                                                                                                                 |1.1-2       |GPL-3                                       |
|matrixStats               |Functions that Apply to Rows and Columns of Matrices (and to
Vectors)                                                                                              |1.3.0       |Artistic-2.0                                |
|SeuratObject              |Data Structures for Single Cell Data                                                                                                                              |5.0.2       |MIT + file LICENSE                          |
|cols4all                  |Colors for all                                                                                                                                                    |0.7-1       |GPL-3                                       |
|km.ci                     |Confidence Intervals for the Kaplan-Meier Estimator                                                                                                               |0.5-6       |GPL (>= 2)                                  |
|fastmatch                 |Fast 'match()' Function                                                                                                                                           |1.1-4       |GPL-2                                       |
|gower                     |Gower's Distance                                                                                                                                                  |1.0.1       |GPL-3                                       |
|checkmate                 |Fast and Versatile Argument Checks                                                                                                                                |2.3.1       |BSD_3_clause + file LICENSE                 |
|hardhat                   |Construct Modeling Packages                                                                                                                                       |1.3.1       |MIT + file LICENSE                          |
|MatrixGenerics            |S4 Generic Summary Statistic Functions that Operate on
Matrix-Like Objects                                                                                         |1.16.0      |Artistic-2.0                                |
|SnowballC                 |Snowball Stemmers Based on the C 'libstemmer' UTF-8 Library                                                                                                       |0.7.1       |BSD_3_clause + file LICENSE                 |
|spatstat.sparse           |Sparse Three-Dimensional Arrays and Linear Algebra Utilities                                                                                                      |3.1-0       |GPL (>= 2)                                  |
|htmlTable                 |Advanced Tables for Markdown/HTML                                                                                                                                 |2.4.2       |GPL (>= 3)                                  |
|rngtools                  |Utility Functions for Working with Random Number Generators                                                                                                       |1.5.2       |GPL-3                                       |
|RANN                      |Fast Nearest Neighbour Search (Wraps ANN Library) Using L2
Metric                                                                                                  |2.6.2       |GPL (>= 3)                                  |
|S4Vectors                 |Foundation of vector-like and list-like containers in
Bioconductor                                                                                                 |0.42.0      |Artistic-2.0                                |
|spatstat.explore          |Exploratory Data Analysis for the 'spatstat' Family                                                                                                               |3.3-2       |GPL (>= 2)                                  |
|IRanges                   |Foundation of integer range manipulation in Bioconductor                                                                                                          |2.38.0      |Artistic-2.0                                |
|bit64                     |A S3 Class for Vectors of 64bit Integers                                                                                                                          |4.0.5       |GPL-2 &#124; GPL-3                          |
|cluster                   |"Finding Groups in Data": Cluster Analysis Extended Rousseeuw et
al.                                                                                               |2.1.6       |GPL (>= 2)                                  |
|farver                    |High Performance Colour Space Manipulation                                                                                                                        |2.1.2       |MIT + file LICENSE                          |
|zip                       |Cross-Platform 'zip' Compression                                                                                                                                  |2.3.1       |MIT + file LICENSE                          |
|gplots                    |Various R Programming Tools for Plotting Data                                                                                                                     |3.1.3.1     |GPL-2                                       |

