#' @title Run Regularized Generalized Canonical Correlation Analysis (RGCCA) for Multi-Modality Data Integration
#' @description This function performs clustering analysis using Regularized Generalized Canonical Correlation Analysis (RGCCA), which integrates multiple data modalities (e.g., RNA, protein, methylation) to provide a unified clustering result. RGCCA helps in capturing shared information between multiple data types and provides insights into the relationships between different modalities.
#' @author Zaoqu Liu; Email: liuzaoqu@163.com
#' @param data A list of matrices where each element represents a different modality (e.g., RNA, protein, methylation). Each matrix should have rows as features and columns as samples.
#' @param N.clust Integer. Number of clusters to create from the hierarchical clustering of the RGCCA components (optional but recommended).
#' @param connection.matrix Matrix. Connection matrix specifying the relationships between blocks (default: 1 - diagonal matrix of block length).
#' @param num.components Integer vector. Number of components to compute for each block (default: 1 component per block).
#' @param scheme Character. The RGCCA scheme to use, can be one of "centroid", "factorial", or "horst" (default: "centroid").
#' @param regularization Character or numeric vector. Regularization parameter for each block, can be "optimal" or a numeric value (default: "optimal").
#' @param scale Logical. Whether to scale each block to zero mean and unit variance (default: TRUE).
#' @param initialization Character. Initialization method for the RGCCA algorithm, either "svd" or "random" (default: "svd").
#' @param bias Logical. Whether to use biased or unbiased estimator of the variance/covariance (default: TRUE).
#' @param tolerance Numeric. Convergence tolerance (default: 1e-08).
#' @param verbose Logical. Whether to show progress messages (default: FALSE).
#' @param max.iterations Integer. Maximum number of iterations for the RGCCA algorithm (default: 1000).
#' @param clustering.algorithm Character. The clustering algorithm to use for hierarchical clustering (default: "ward.D2").
#' @return A data frame with the following columns:
#'   - Sample: The sample identifier.
#'   - Cluster: The assigned cluster number for each sample.
#'   - Cluster2: The assigned cluster label, prefixed by 'RGCCA' to indicate that the clustering was performed using RGCCA.
#' @details This function uses RGCCA to integrate multiple data matrices and then performs hierarchical clustering on the resulting components to identify clusters. RGCCA is particularly useful for identifying shared structures across multiple modalities, providing a comprehensive view of the relationships between different data types.
#'
#' The function operates as follows:
#' 1. Each matrix in the input list is transposed so that rows represent samples and columns represent features.
#' 2. RGCCA is performed to extract components that summarize the shared variation across different modalities.
#' 3. Hierarchical clustering is applied to the concatenated components to assign each sample to a cluster.
#' 4. The function returns a data frame containing the cluster assignment for each sample, along with additional information about the clustering process.
#' @references
#' Tenenhaus, M., Tenenhaus, A., & Groenen, P. J. (2017). Regularized generalized canonical correlation analysis: a framework for sequential multiblock component methods. Psychometrika, 82(3), 737-777.
#' Tenenhaus, A., Philippe, C., & Frouin, V. (2015). Kernel generalized canonical correlation analysis. Computational Statistics & Data Analysis, 90, 114-131.
#' @examples
#' # Example usage:
#' data1 <- matrix(rnorm(10000), nrow = 100, ncol = 100)
#' data2 <- matrix(rnorm(10000), nrow = 100, ncol = 100)
#' colnames(data1) <- colnames(data2) <- paste0("Sample", 1:100)
#' data_list <- list(data1, data2)
#'
#' # Run RGCCA clustering
#' result <- RunRGCCA(data = data.list, N.clust = 3)
#' @export
RunRGCCA <- function(data, N.clust = NULL, connection.matrix = 1 - diag(length(data)),
                     num.components = rep(1, length(data)), scheme = "centroid",
                     regularization = "optimal", scale = TRUE, initialization = "svd",
                     bias = TRUE, tolerance = 1e-08, verbose = FALSE,
                     clustering.algorithm = "ward.D2") {
  # Check if required packages are installed
  if (!requireNamespace("devtools", quietly = TRUE)) {
    message("Installing devtools package from CRAN...")
    install.packages("devtools")
  }
  if (!requireNamespace("RGCCA", quietly = TRUE)) {
    message("Installing RGCCA package from GitHub...")
    devtools::install_github("Zaoqu-Liu/RGCCA")
  }

  # Input validation checks
  if (!is.list(data)) {
    stop("The 'data' parameter must be a list of matrices.")
  }

  # Transpose each data matrix to have rows as samples and columns as features
  data <- lapply(data, function(x) t(x))

  # Perform RGCCA to extract components summarizing shared variation across modalities
  result.rgcca <- RGCCA::rgcca(
    A = data, C = connection.matrix, tau = regularization,
    ncomp = num.components, scheme = scheme, scale = scale,
    init = initialization, bias = bias, tol = tolerance,
    verbose = verbose
  )
  resDat <- do.call(cbind, result.rgcca$Y)

  # Perform hierarchical clustering on the concatenated components
  dist.matrix <- dist(resDat)
  hc <- hclust(dist.matrix, method = clustering.algorithm)
  clust.rgcca <- cutree(hc, k = N.clust)

  # Create a data frame with cluster assignments
  clustres <- data.frame(
    Sample = rownames(data[[1]]), Cluster = clust.rgcca,
    Cluster2 = paste0("RGCCA", clust.rgcca),
    row.names = rownames(data[[1]]), stringsAsFactors = FALSE
  )

  # Print summary of clustering results
  cat("Cluster algorithm: RGCCA\n")
  if (!is.null(N.clust)) {
    cat(paste0("Selection of cluster number: ", N.clust, "\n"))
    cluster_summary <- paste(
      paste0("C", 1:N.clust),
      round(table(clustres$Cluster) / nrow(clustres), 2) * 100,
      sep = " = "
    )
    cat(paste(cluster_summary, collapse = "% "), "%\n")
  } else {
    cat("Number of clusters was determined automatically by RGCCA.\n")
  }

  return(clustres)
}
