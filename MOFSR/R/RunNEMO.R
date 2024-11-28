#' @title Run NEMO for Multi-Modality Data Clustering
#' @description This function performs clustering analysis using the NEMO method, which is suitable for identifying subtypes in multi-omics data. NEMO allows integration of multiple modalities to discover coherent subtypes by partial data integration.
#' @author Zaoqu Liu; Email: liuzaoqu@163.com
#' @param data A list of matrices where each element represents a different modality (e.g., RNA, protein, methylation). Each matrix should have rows as features and columns as samples.
#' @param N.clust Integer. Number of clusters to evaluate during subtyping. If NULL, NEMO will determine the optimal number of clusters automatically. Default is NULL.
#' @return A data frame with the following columns:
#'   - Sample: The sample identifier, taken from the column names of the input data matrices.
#'   - Cluster: The assigned cluster number for each sample, indicating the subtype classification.
#'   - Cluster2: The assigned cluster label, prefixed by 'NEMO' to indicate that the clustering was performed using NEMO.
#' @details The function operates as follows:
#' 1. Uses the NEMO package to identify subtypes across the provided omics modalities.
#' 2. Returns a data frame containing the cluster assignment for each sample.
#' @references Rappoport N, Shamir R. NEMO: cancer subtyping by integration of partial multi-omic data. Bioinformatics. 2019;35(18):3348-3356. doi:10.1093/bioinformatics/btz058
#' @examples
#' # Example usage:
#' data1 <- matrix(rnorm(10000), nrow = 100, ncol = 100)
#' data2 <- matrix(rnorm(10000), nrow = 100, ncol = 100)
#' colnames(data1) <- colnames(data2) <- paste0("Sample", 1:100)
#' data_list <- list(data1, data2)
#'
#' # Run NEMO clustering
#' result <- RunNEMO(data = data_list, N.clust = 3)
#' @export
RunNEMO <- function(data, N.clust = NULL) {
  # Set seed for reproducibility
  set.seed(1)

  # Check if required packages are installed
  if (!requireNamespace("NEMO", quietly = TRUE)) {
    if (!requireNamespace("devtools", quietly = TRUE)) {
      install.packages("devtools")
    }
    devtools::install_github("Shamir-Lab/NEMO/NEMO")
  }
  if (!requireNamespace("SNFtool", quietly = TRUE)) {
    install.packages("SNFtool")
  }

  library(SNFtool)
  # Input validation checks
  if (!is.list(data)) {
    stop("The 'data' parameter must be a list of matrices.")
  }

  # Ensure that all data matrices have the same number of columns (samples)
  num_samples <- ncol(data[[1]])
  if (!all(sapply(data, ncol) == num_samples)) {
    stop("All matrices must have the same number of samples (columns).")
  }

  # Run NEMO clustering analysis
  fit <- NEMO::nemo.clustering(omics.list = data, num.clusters = N.clust, num.neighbors = NA)

  # Create a data frame with cluster assignments
  clustres <- data.frame(
    Sample = colnames(data[[1]]), # Sample identifiers
    Cluster = fit, # Cluster number assigned by NEMO
    Cluster2 = paste0("NEMO", fit), # Prefixed cluster labels
    row.names = colnames(data[[1]]),
    stringsAsFactors = FALSE
  )

  # Print summary of clustering results
  cat("Cluster algorithm: NEMO\n")
  if (!is.null(N.clust)) {
    cat(paste0("Selection of cluster number: ", N.clust, "\n"))
    cluster_summary <- paste(
      paste0("C", 1:N.clust),
      round(table(clustres$Cluster) / nrow(clustres), 2) * 100,
      sep = " = "
    )
    cat(paste(cluster_summary, collapse = "% "), "%\n")
  } else {
    cat("Number of clusters was determined automatically by NEMO.\n")
  }

  return(clustres)
}
