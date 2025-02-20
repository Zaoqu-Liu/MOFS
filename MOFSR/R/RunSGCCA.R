#' @title Run Sparse Generalized Canonical Correlation Analysis (SGCCA) for Multi-Modality Data Integration
#' @description This function performs clustering analysis using Sparse Generalized Canonical Correlation Analysis (SGCCA). SGCCA is an extension of Generalized Canonical Correlation Analysis that allows for sparse estimation, which is beneficial when working with high-dimensional datasets, making it effective for integrating multiple modalities (e.g., RNA, protein).
#' @author Zaoqu Liu; Email: liuzaoqu@163.com
#' @param data A list of matrices where each element represents a different modality (e.g., RNA, protein, methylation). Each matrix should have rows as features and columns as samples.
#' @param N.clust Integer. Number of clusters for hierarchical clustering (optional but recommended).
#' @param connection.matrix Matrix. A matrix describing the relationships between different modalities (default: complete design matrix).
#' @param num.components.per.modality Integer vector specifying the number of components to compute for each modality (default: 1 component per modality).
#' @param integration.scheme Character. The method used for integrating different data modalities. Options are "centroid", "horst", or "factorial" (default: "centroid").
#' @param sparsity.level Numeric vector specifying the regularization (sparsity) parameters for each modality, controlling the sparsity level (default: 0.5 for each modality).
#' @param scale.data Logical. Whether to scale each block to zero mean and unit variance (default: FALSE).
#' @param initialization.method Character. Method for initializing the SGCCA algorithm, either "svd" or "random" (default: "svd").
#' @param use.biased.variance Logical. Whether to use a biased estimator for the variance/covariance (default: TRUE).
#' @param convergence.tolerance Numeric. The convergence tolerance for the algorithm (default: .Machine$double.eps).
#' @param show.progress Logical. Whether to display progress messages during the execution (default: FALSE).
#' @param cluster.algorithm Character. The method used for hierarchical clustering (default: "ward.D2").
#' @return A data frame with the following columns:
#'   - Sample: The sample identifier.
#'   - Cluster: The assigned cluster number for each sample.
#'   - Cluster2: The assigned cluster label, prefixed by 'SGCCA' to indicate that the clustering was performed using SGCCA.
#' @details The function proceeds as follows:
#' 1. Each matrix in the input list is transposed so that rows represent samples and columns represent features.
#' 2. SGCCA is applied to find the canonical components for each modality.
#' 3. Combines the canonical components and applies hierarchical clustering to identify clusters.
#' 4. Returns a data frame with the assigned clusters for each sample.
#' @references
#' Tenenhaus, M., Tenenhaus, A., & Groenen, P. J. (2017). Regularized generalized canonical correlation analysis: a framework for sequential multiblock component methods. Psychometrika, 82(3), 737-777.
#' Tenenhaus, A., Philippe, C., & Frouin, V. (2015). Kernel generalized canonical correlation analysis. Computational Statistics & Data Analysis, 90, 114-131.
#' Tenenhaus, A., Philippe, C., Guillemot, V., Le Cao, K. A., Grill, J., & Frouin, V. (2014). Variable selection for generalized canonical correlation analysis. Biostatistics, 15(3), 569-583.
#' Tenenhaus, A., & Tenenhaus, M. (2011). Regularized generalized canonical correlation analysis. Psychometrika, 76(2), 257.
#' Van de Geer, J. P. (1984). Linear relations among K sets of variables. Psychometrika, 49(1), 79-94.
#' Schafer J. and Strimmer K. (2005). A shrinkage approach to large-scale covariance matrix estimation and implications for functional genomics. Statistical Applications in Genetics and Molecular Biology 4:32.
#' Tenenhaus et al. Variable Selection For Generalized Canonical Correlation Analysis. 2013. Submitted to Biostatistics.
#' @examples
#' # Example usage:
#' data1 <- matrix(rnorm(10000), nrow = 100, ncol = 100)
#' data2 <- matrix(rnorm(10000), nrow = 100, ncol = 100)
#' colnames(data1) <- colnames(data2) <- paste0("Sample", 1:100)
#' data_list <- list(data1, data2)
#'
#' # Run SGCCA clustering
#' result <- RunSGCCA(data = data_list, N.clust = 3)
#'
#' @export
RunSGCCA <- function(data, N.clust = NULL, connection.matrix = 1 - diag(length(data)),
                     num.components.per.modality = rep(1, length(data)), integration.scheme = "centroid",
                     sparsity.level = rep(0.5, length(data)), scale.data = FALSE, initialization.method = "svd",
                     use.biased.variance = TRUE, convergence.tolerance = .Machine$double.eps, show.progress = FALSE,
                     cluster.algorithm = "ward.D2") {
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
  data <- lapply(data, function(x) {
    t(x)
  })

  # Perform SGCCA
  result.sgcca <- RGCCA::sgcca(
    A = data, C = connection.matrix, c1 = sparsity.level,
    ncomp = num.components.per.modality, scheme = integration.scheme, scale = scale.data,
    init = initialization.method, bias = use.biased.variance, tol = convergence.tolerance, verbose = show.progress
  )
  resDat <- do.call(cbind, result.sgcca$Y)

  # Apply hierarchical clustering on the canonical components
  clust.sgcca <- resDat %>%
    dist() %>%
    hclust(method = cluster.algorithm) %>%
    cutree(k = N.clust)

  # Create a data frame with cluster assignments
  clustres <- data.frame(
    Sample = rownames(data[[1]]),
    Cluster = clust.sgcca,
    Cluster2 = paste0("SGCCA", clust.sgcca),
    row.names = rownames(data[[1]]),
    stringsAsFactors = FALSE
  )

  # Print summary of clustering results
  cat(paste0("Cluster algorithm: SGCCA\n"))
  if (!is.null(N.clust)) {
    cat(paste0("Selection of cluster number: ", N.clust, "\n"))
    cluster_summary <- paste(
      paste0("C", 1:N.clust),
      round(table(clustres$Cluster) / nrow(clustres), 2) * 100,
      sep = " = "
    )
    cat(paste(cluster_summary, collapse = "% "), "%\n")
  } else {
    cat("Number of clusters was determined automatically by SGCCA.\n")
  }

  return(clustres)
}
