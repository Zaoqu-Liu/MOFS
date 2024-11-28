#' @title Run Consensus Iterative Multi-view Learning (CIMLR) for Multi-Modality Data Integration
#' @description This function performs clustering analysis using Consensus Iterative Multi-view Learning (CIMLR) to integrate multiple omics datasets. CIMLR is useful for discovering shared patterns across multiple datasets and identifying distinct subtypes.
#' @author Zaoqu Liu; Email: liuzaoqu@163.com
#' @param data A list of matrices where each element represents a different modality (e.g., RNA, protein, methylation). Each matrix should have rows as features and columns as samples.
#' @param N.clust Integer. Number of clusters to create (optional but recommended).
#' @param num.dimensions Integer. Number of dimensions for CIMLR (default: NA).
#' @param tuning.parameter Integer. Tuning parameter for CIMLR (default: 10).
#' @param cores.ratio Numeric. Ratio of the number of cores to be used when computing the multi-kernel (default: 1).
#' @return A data frame with the following columns:
#'   - Sample: The sample identifier.
#'   - Cluster: The assigned cluster number for each sample.
#'   - Cluster2: The assigned cluster label, prefixed by 'CIMLR' to indicate that the clustering was performed using CIMLR.
#' @details This function uses CIMLR to integrate multiple data matrices and assign clusters to the samples. CIMLR is particularly effective for multi-view learning and discovering common patterns among different data types.
#'
#' The function operates as follows:
#' 1. CIMLR is performed using the CIMLR package to extract components that summarize the shared variation across different modalities.
#' 2. Feature ranking is applied to identify the most important features across all data.
#' 3. The function returns a data frame containing the cluster assignment for each sample, along with additional information about the clustering process.
#' @references Ramazzotti D, Lal A, Wang B, Batzoglou S, Sidow A. Multi-omic tumor data reveal diversity of molecular mechanisms that correlate with survival. Nat Commun. 2018;9(1):4453. doi:10.1038/s41467-018-06921-8.
#' @examples
#' # Example usage:
#' data1 <- matrix(rnorm(10000), nrow = 100, ncol = 100)
#' data2 <- matrix(rnorm(10000), nrow = 100, ncol = 100)
#' colnames(data1) <- colnames(data2) <- paste0("Sample", 1:100)
#' data_list <- list(data1, data2)
#'
#' # Run CIMLR clustering
#' result <- RunCIMLR(data = data_list, N.clust = 3)
#'
#' @export
RunCIMLR <- function(data, N.clust = NULL, num.dimensions = NA, tuning.parameter = 10, cores.ratio = 1) {
  # Set seed for reproducibility
  set.seed(123)

  # Check if required packages are installed
  if (!requireNamespace("CIMLR", quietly = TRUE)) {
    if (!requireNamespace("devtools", quietly = TRUE)) {
      install.packages("devtools")
    }
    devtools::install_github("danro9685/CIMLR", ref = "R")
  }

  # Perform CIMLR clustering
  fit <- zquiet(CIMLR::CIMLR(
    X = data, c = N.clust, no.dim = num.dimensions,
    k = tuning.parameter, cores.ratio = cores.ratio
  ))

  # Prepare input data for feature ranking
  input_data <- do.call(rbind, lapply(seq_along(data), function(x) {
    data_matrix <- data[[x]]
    rownames(data_matrix) <- paste(rownames(data_matrix), names(data)[x], sep = "+")
    data_matrix
  }))

  # Perform feature ranking
  ranks <- zquiet(CIMLR::CIMLR_Feature_Ranking(A = fit$S, X = input_data))
  ranks$names <- rownames(input_data)[ranks$aggR]
  fit$selected_features <- ranks

  # Create a data frame with cluster assignments
  cluster_results <- data.frame(
    Sample = colnames(data[[1]]),
    Cluster = fit$y$cluster,
    Cluster2 = paste0("CIMLR", fit$y$cluster),
    row.names = colnames(data[[1]]),
    stringsAsFactors = FALSE
  )

  # Print summary of clustering results
  cat(paste0("Cluster algorithm: CIMLR\n"))
  if (!is.null(N.clust)) {
    cat(paste0("Selection of cluster number: ", N.clust, "\n"))
    cluster_summary <- paste(
      paste0("C", 1:N.clust),
      round(table(cluster_results$Cluster) / nrow(cluster_results), 2) * 100,
      sep = " = "
    )
    cat(paste(cluster_summary, collapse = "% "), "%\n")
  } else {
    cat("Number of clusters was determined automatically by CIMLR.\n")
  }

  return(cluster_results)
}
