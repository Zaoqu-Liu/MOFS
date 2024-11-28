#' @title Run PINSPlus for Multi-Modality Data Clustering
#' @description This function performs clustering analysis using the PINSPlus method, which is suitable for identifying subtypes in multi-omics data. PINSPlus allows integration of multiple modalities to discover coherent subtypes.
#' @author Zaoqu Liu; Email: liuzaoqu@163.com
#' @param data A list of matrices where each element represents a different modality (e.g., RNA, protein, methylation). Each matrix should have rows as features and columns as samples. All matrices should have the same number of samples.
#' @param N.clust Integer. Maximum number of clusters to evaluate during subtyping. If NULL, PINSPlus will determine the optimal number of clusters automatically. Default is NULL.
#' @param agreement.cutoff Numeric. Agreement threshold to be considered consistent (default: 0.5).
#' @param num.cores Integer. Number of cores for parallel processing (default: 10). Setting a higher value can speed up the computation.
#' @param sampled.set.size Integer. The number of sample size used for the sampling process when the dataset is large (default: 2000).
#' @param knn.k Integer. The value of k for the k-nearest neighbors algorithm. If not set, elbow method will be used to calculate k (default: NULL).
#' @return A data frame with the following columns:
#'   - Sample: The sample identifier, taken from the column names of the input data matrices.
#'   - Cluster: The assigned cluster number for each sample, indicating the subtype classification.
#'   - Cluster2: The assigned cluster label, prefixed by 'PINSPlus' to indicate that the clustering was performed using PINSPlus.
#' @details The function operates as follows:
#' 1. Transposes each data matrix so that rows represent samples and columns represent features.
#' 2. Uses the PINSPlus package to identify subtypes across the modalities.
#' 3. Returns a data frame containing the cluster assignment for each sample.
#' @references
#' Nguyen H, Shrestha S, Draghici S, Nguyen T. PINSPlus: a tool for tumor subtype discovery in integrated genomic data. Bioinformatics. 2019;35(16):2843-2846. doi:10.1093/bioinformatics/bty1049
#' Nguyen T, Tagett R, Diaz D, Draghici S. A novel method for data integration and disease subtyping. Genome Research. 2017;27(12):2025-2039.
#' Nguyen T. Horizontal and vertical integration of bio-molecular data. PhD thesis, Wayne State University. 2017.
#' Nguyen H, Tran D, Tran B, Roy M, Cassell A, Dascalu S, Draghici S, Nguyen T. SMRT: Randomized Data Transformation for Cancer Subtyping and Big Data Analysis. Frontiers in Oncology. 2021.
#' @examples
#' # Example usage:
#' data1 <- matrix(rnorm(10000), nrow = 100, ncol = 100)
#' data2 <- matrix(rnorm(10000), nrow = 100, ncol = 100)
#' colnames(data1) <- colnames(data2) <- paste0("Sample", 1:100)
#' data_list <- list(data1, data2)
#'
#' # Run PINSPlus clustering
#' result <- RunPINSPlus(data = data_list, N.clust = 3)
#' @export
RunPINSPlus <- function(data, N.clust = NULL, agreement.cutoff = 0.5,
                        num.cores = 10, sampled.set.size = 2000, knn.k = NULL) {
  # Check if required packages are installed
  if (!requireNamespace("PINSPlus", quietly = TRUE)) {
    message("PINSPlus package is not installed. Installing from CRAN...")
    install.packages("PINSPlus")
  }

  # Input validation checks
  if (!is.list(data)) {
    stop("The 'data' parameter must be a list of matrices.")
  }

  # Transpose each data matrix to have rows as samples and columns as features
  data <- lapply(data, t)

  # Ensure that all data matrices have the same number of rows (samples)
  num_samples <- nrow(data[[1]])
  if (!all(sapply(data, nrow) == num_samples)) {
    stop("All matrices must have the same number of samples (rows).")
  }

  # Set seed for reproducibility
  set.seed(123)

  # Run subtyping analysis using PINSPlus
  fit <- PINSPlus::SubtypingOmicsData(
    dataList = data, kMax = N.clust,
    agreementCutoff = agreement.cutoff, ncore = num.cores,
    sampledSetSize = sampled.set.size, knn.k = knn.k, verbose = FALSE
  )

  # Create a data frame with cluster assignments
  clustres <- data.frame(
    Sample = rownames(data[[1]]), # Sample identifiers
    Cluster = fit$cluster1, # Cluster number assigned by PINSPlus
    Cluster2 = paste0("PINSPlus", fit$cluster1), # Prefixed cluster labels
    row.names = rownames(data[[1]]),
    stringsAsFactors = FALSE
  )

  # Print summary of clustering results
  cat(paste0("Cluster algorithm: PINSPlus\n"))
  if (!is.null(N.clust)) {
    cat(paste0("Selection of cluster number: ", N.clust, "\n"))
    cluster_summary <- paste(
      paste0("C", 1:N.clust),
      round(table(clustres$Cluster) / nrow(clustres), 2) * 100,
      sep = " = "
    )
    cat(paste(cluster_summary, collapse = "% "), "%\n")
  } else {
    cat("Number of clusters was determined automatically by PINSPlus.\n")
  }

  return(clustres)
}
