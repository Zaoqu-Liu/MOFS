#' @title Run Bayesian Consensus Clustering (BCC) for Multi-Modality Data Integration
#' @description This function performs clustering analysis using Bayesian Consensus Clustering (BCC) to integrate multiple omics datasets. BCC is a Bayesian method that helps capture shared patterns across different datasets by finding a consensus clustering solution.
#' @author Zaoqu Liu; Email: liuzaoqu@163.com
#' @param data A list of matrices where each element represents a different modality (e.g., RNA, protein, methylation). Each matrix should have rows as features and columns as samples.
#' @param N.clust Integer. Number of clusters to create from the BCC components (optional but recommended).
#' @param max.iterations Integer. Maximum number of iterations for the Bayesian algorithm (default: 10).
#' @return A data frame with the following columns:
#'   - Sample: The sample identifier.
#'   - Cluster: The assigned cluster number for each sample.
#'   - Cluster2: The assigned cluster label, prefixed by 'BCC' to indicate that the clustering was performed using BCC.
#' @details This function uses BCC to integrate multiple data matrices and assign clusters to the samples. BCC uses a Bayesian approach to identify shared patterns across multiple datasets, providing a robust clustering solution.
#'
#' The function operates as follows:
#' 1. Each matrix in the input list is converted to a matrix to ensure compatibility.
#' 2. BCC is used to identify shared patterns across different modalities.
#' 3. The function returns a data frame containing the cluster assignment for each sample, along with additional information about the clustering process.
#' @references Lock EF, Dunson DB. Bayesian consensus clustering. Bioinformatics. 2013;29(20):2610-2616. doi:10.1093/bioinformatics/btt425.
#' @examples
#' # Example usage:
#' data1 <- matrix(rnorm(10000), nrow = 100, ncol = 100)
#' data2 <- matrix(rnorm(10000), nrow = 100, ncol = 100)
#' colnames(data1) <- colnames(data2) <- paste0("Sample", 1:100)
#' data_list <- list(data1, data2)
#'
#' # Run BCC clustering
#' result <- RunBCC(data = data_list, N.clust = 3)
#'
#' @export
RunBCC <- function(data = NULL, N.clust = NULL, max.iterations = 10) {
  # Set seed for reproducibility
  set.seed(1)

  # Check if required packages are installed
  if (!requireNamespace("bayesCC", quietly = TRUE)) {
    if (!requireNamespace("devtools", quietly = TRUE)) {
      install.packages("devtools")
    }
    devtools::install_github("ttriche/bayesCC")
  }

  if (!requireNamespace("dirmult", quietly = TRUE)) {
    install.packages("dirmult")
  }
  require(dirmult)

  # Input validation checks
  if (!is.list(data)) {
    stop("The 'data' parameter must be a list of matrices.")
  }

  # Convert each element to a matrix to ensure compatibility
  data <- lapply(data, as.matrix)

  # Perform Bayesian Consensus Clustering
  fit <- bayesCC::bayesCC(data, K = N.clust, IndivAlpha = TRUE, maxiter = max.iterations)

  # Create a data frame with cluster assignments
  cluster.assignments <- as.data.frame(fit$Cbest)
  cluster.assignments$Cluster <- apply(cluster.assignments, 1, which.max)

  cluster.results <- data.frame(
    Sample = colnames(data[[1]]),
    Cluster = cluster.assignments$Cluster,
    Cluster2 = paste0("BCC", cluster.assignments$Cluster),
    row.names = colnames(data[[1]]),
    stringsAsFactors = FALSE
  )

  # Print summary of clustering results
  cat(paste0("Cluster algorithm: BCC\n"))
  if (!is.null(N.clust)) {
    cat(paste0("Selection of cluster number: ", N.clust, "\n"))
    cluster_summary <- paste(
      paste0("C", 1:N.clust),
      round(table(cluster.results$Cluster) / nrow(cluster.results), 2) * 100,
      sep = " = "
    )
    cat(paste(cluster_summary, collapse = "% "), "%\n")
  } else {
    cat("Number of clusters was determined automatically by BCC.\n")
  }

  return(cluster.results)
}
