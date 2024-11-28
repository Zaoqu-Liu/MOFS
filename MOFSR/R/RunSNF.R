#' @title Run Similarity Network Fusion (SNF) for Multi-Modality Data Integration
#' @description This function performs clustering analysis using Similarity Network Fusion (SNF), which integrates multiple data types (e.g., different modalities such as RNA, protein, methylation) to provide a unified clustering result. SNF is an effective technique for capturing complementary information from multiple modalities and determining common subtypes.
#' @author Zaoqu Liu; Email: liuzaoqu@163.com
#' @param data A list of matrices where each element represents a different modality (e.g., RNA, protein, methylation). Each matrix should have rows as features and columns as samples.
#' @param N.clust Integer. Number of clusters for spectral clustering. This is the desired number of groups to partition the samples into (optional but recommended).
#' @param num.neighbors Integer. Number of nearest neighbors to consider in building the affinity matrix (default: 20). This parameter controls the local neighborhood size. It typically takes values between 10 and 30.
#' @param variance Numeric. The variance for the local model in building the affinity matrix (default: 0.5). This value affects the distance weighting between samples. It typically ranges from 0.3 to 0.8.
#' @param num.iterations Integer. Number of iterations for the similarity network fusion process (default: 20). This parameter controls how the affinity matrices from different data types are fused. Typically, 10 to 20 iterations are used.
#' @return A data frame with the following columns:
#'   - Sample: The sample identifier.
#'   - Cluster: The assigned cluster number for each sample.
#'   - Cluster2: The assigned cluster label, prefixed by 'SNF' to indicate that the clustering was performed using SNF.
#' @details This function uses the Similarity Network Fusion (SNF) approach to integrate multiple data matrices and then performs spectral clustering on the fused network to identify clusters. SNF is particularly useful when dealing with multi-modality datasets, as it takes into account the complementary nature of different data types to improve clustering robustness and biological interpretability.
#'
#' The function operates as follows:
#' 1. Each matrix in the input list is transposed so that rows represent samples and columns represent features.
#' 2. For each data type, an affinity matrix is computed using the distance between samples and then transformed using an exponential kernel with parameters `num.neighbors` and `variance`.
#' 3. The affinity matrices are fused using the SNF approach over `num.iterations` iterations to create a consensus similarity network.
#' 4. Spectral clustering is then applied to the fused affinity matrix to assign each sample to a cluster.
#' 5. The function returns a data frame containing the cluster assignment for each sample, along with additional information about the clustering process.
#' @references Wang B, Mezlini AM, Demir F, Fiume M, Tu Z, Brudno M, Haibe-Kains B, Goldenberg A. Similarity Network Fusion for Aggregating Data Types on a Genomic Scale. Nat Methods. 2014;11(3):333-337. doi:10.1038/nmeth.2810
#' @examples
#' # Example usage:
#' data1 <- matrix(rnorm(10000), nrow = 100, ncol = 100)
#' data2 <- matrix(rnorm(10000), nrow = 100, ncol = 100)
#' colnames(data1) <- colnames(data2) <- paste0("Sample", 1:100)
#' data_list <- list(data1, data2)
#'
#' # Run SNF clustering
#' result <- RunSNF(data = data_list, N.clust = 3)
#'
#' @export
RunSNF <- function(data, N.clust = NULL, num.neighbors = 20, variance = 0.5, num.iterations = 20) {
  # Check if required packages are installed
  if (!requireNamespace("SNFtool", quietly = TRUE)) {
    message("SNFtool package is not installed. Installing from CRAN...")
    install.packages("SNFtool")
  }

  # Input validation checks
  if (!is.list(data)) {
    stop("The 'data' parameter must be a list of matrices.")
  }

  # Transpose each data matrix to have rows as samples and columns as features
  data <- lapply(data, function(x) {
    t(x)
  })

  # Generate affinity matrices for each data type
  affinity.matrices <- lapply(data, function(modality) {
    modality <- as.matrix(modality) # Ensure data is a matrix
    affinity <- SNFtool::affinityMatrix(SNFtool::dist2(modality, modality), K = num.neighbors, sigma = variance) # Compute distance and build affinity matrix
    return(affinity)
  })

  # Fuse affinity matrices using SNF
  fused.network <- SNFtool::SNF(Wall = affinity.matrices, K = num.neighbors, t = num.iterations)

  # Perform spectral clustering on the fused affinity matrix
  cluster.assignments <- SNFtool::spectralClustering(fused.network, N.clust)

  # Create a data frame with cluster assignments
  cluster.results <- data.frame(
    Sample = rownames(data[[1]]),
    Cluster = cluster.assignments,
    Cluster2 = paste0("SNF", cluster.assignments),
    row.names = rownames(data[[1]]),
    stringsAsFactors = FALSE
  )

  # Print summary of clustering results
  cat(paste0("Cluster algorithm: SNF\n"))
  if (!is.null(N.clust)) {
    cat(paste0("Selection of cluster number: ", N.clust, "\n"))
    cluster_summary <- paste(
      paste0("C", 1:N.clust),
      round(table(cluster.results$Cluster) / nrow(cluster.results), 2) * 100,
      sep = " = "
    )
    cat(paste(cluster_summary, collapse = "% "), "%\n")
  } else {
    cat("Number of clusters was determined automatically by SNF.\n")
  }

  return(cluster.results)
}
