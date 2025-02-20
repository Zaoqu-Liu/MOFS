#' @title Run Integrative Non-negative Matrix Factorization (IntNMF) for Multi-Modality Data Integration
#' @description This function performs clustering analysis using Integrative Non-negative Matrix Factorization (IntNMF) to integrate multiple omics datasets. IntNMF is a powerful tool for capturing shared patterns across multiple datasets by decomposing them into components that reflect shared and individual structures.
#' @author Zaoqu Liu; Email: liuzaoqu@163.com
#' @param data A list of matrices where each element represents a different modality (e.g., RNA, protein, methylation). Each matrix should have rows as features and columns as samples.
#' @param N.clust Integer. Number of clusters to create from the IntNMF components (optional but recommended).
#' @param max.iterations Integer. Maximum number of iterations for the NMF algorithm (default: 200).
#' @param stability.count Integer. Count for stability in connectivity matrix (default: 20).
#' @param num.initializations Integer. Number of initializations of the random matrices (default: 30).
#' @param use.nndsvd Logical. Whether to use non-negative double singular value decomposition (NNDSVD) for initialization (default: TRUE).
#' @param random.seed Logical. Whether to use a random seed for initialization of the algorithm (default: TRUE).
#' @param weight Numeric vector. Weight for each data matrix in the list (default: equal weights).
#' @return A data frame with the following columns:
#'   - Sample: The sample identifier.
#'   - Cluster: The assigned cluster number for each sample.
#'   - Cluster2: The assigned cluster label, prefixed by 'IntNMF' to indicate that the clustering was performed using IntNMF.
#' @details This function uses IntNMF to integrate multiple data matrices and then assigns clusters to the samples based on the resulting components. IntNMF is particularly useful for capturing shared and individual variation in multi-omics data.
#'
#' The function operates as follows:
#' 1. Each data matrix is normalized to ensure non-negative values and scaled to unit variance.
#' 2. IntNMF is performed using the IntNMF package to extract components that summarize the shared variation across different modalities.
#' 3. The function returns a data frame containing the cluster assignment for each sample, along with additional information about the clustering process.
#' @references Chalise P, Fridley BL. Integrative clustering of multi-level 'omic data based on non-negative matrix factorization algorithm. PLoS One. 2017;12(5):e0176278. doi:10.1371/journal.pone.0176278.
#'             Chalise P, Raghavan R and Fridley B (2016). InterSIM: Simulation tool for multiple integrative 'omic datasets. Computer Methods and Programs in Biomedicine, 128:69-74.
#' @examples
#' # Example usage:
#' data1 <- matrix(rnorm(10000), nrow = 100, ncol = 100)
#' data2 <- matrix(rnorm(10000), nrow = 100, ncol = 100)
#' colnames(data1) <- colnames(data2) <- paste0("Sample", 1:100)
#' data_list <- list(data1, data2)
#'
#' # Run IntNMF clustering
#' result <- RunIntNMF(data = data_list, N.clust = 3)
#'
#' @export
RunIntNMF <- function(data = NULL, N.clust = NULL, max.iterations = 5, stability.count = 20,
                      num.initializations = 30, use.nndsvd = TRUE, random.seed = TRUE, weight = NULL) {
  # Set seed for reproducibility
  if (random.seed) {
    set.seed(1)
  }

  # Check if required packages are installed
  if (!requireNamespace("IntNMF", quietly = TRUE)) {
    install.packages("IntNMF")
  }

  # Normalize each data matrix to ensure non-negative values
  data_normalized <- lapply(data, function(dd) {
    if (!all(dd >= 0)) {
      dd <- pmax(dd + abs(min(dd)), 0) + .Machine$double.eps
    }
    dd <- dd / max(dd)
    return(as.matrix(dd))
  })

  # Transpose each matrix to have rows as samples and columns as features
  data_normalized <- lapply(data_normalized, function(x) t(x) + .Machine$double.eps)

  # Set default weight if not provided
  if (is.null(weight)) {
    weight <- if (is.list(data_normalized)) rep(1, length(data_normalized)) else 1
  }

  # Perform IntNMF clustering
  result.intNMF <- IntNMF::nmf.mnnals(
    dat = data_normalized, k = N.clust, maxiter = max.iterations,
    st.count = stability.count, n.ini = num.initializations,
    ini.nndsvd = use.nndsvd, seed = random.seed, wt = weight
  )
  clust.intNMF <- result.intNMF$clusters

  # Create a data frame with cluster assignments
  clustres <- data.frame(
    Sample = colnames(data[[1]]),
    Cluster = as.numeric(clust.intNMF),
    Cluster2 = paste0("IntNMF", as.numeric(clust.intNMF)),
    row.names = colnames(data[[1]]),
    stringsAsFactors = FALSE
  )

  # Print summary of clustering results
  cat(paste0("Cluster algorithm: IntNMF\n"))
  if (!is.null(N.clust)) {
    cat(paste0("Selection of cluster number: ", N.clust, "\n"))
    cluster_summary <- paste(
      paste0("C", 1:N.clust),
      round(table(clustres$Cluster) / nrow(clustres), 2) * 100,
      sep = " = "
    )
    cat(paste(cluster_summary, collapse = "% "), "%\n")
  } else {
    cat("Number of clusters was determined automatically by IntNMF.\n")
  }

  return(clustres)
}
