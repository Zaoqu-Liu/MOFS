#' @title Run Multiple Co-Inertia Analysis (MCIA) for Multi-Modality Data Integration
#' @description This function performs clustering analysis using Multiple Co-Inertia Analysis (MCIA) to integrate multiple omics datasets. MCIA helps capture shared variation across different data modalities and provides insights into the common structures of the data.
#' @author Zaoqu Liu; Email: liuzaoqu@163.com
#' @param data A list of matrices where each element represents a different modality (e.g., RNA, protein, methylation). Each matrix should have rows as features and columns as samples.
#' @param N.clust Integer. Number of clusters to create from the hierarchical clustering of the MCIA components (optional but recommended).
#' @param n.components Integer. Number of components to retain for each modality (default: 10).
#' @param clustering.algorithm Character. The clustering algorithm to use for hierarchical clustering (default: "ward.D2").
#' @param scan.eigenvalues Logical. Whether to show the co-inertia analysis eigenvalue plot to help select the number of axes (default: FALSE).
#' @param use.nsc Logical. Whether to perform multiple non-symmetric correspondence analyses. Recommended to keep TRUE (default: TRUE).
#' @param use.svd Logical. Whether to use singular value decomposition to perform the analysis (default: TRUE).
#' @return A data frame with the following columns:
#'   - Sample: The sample identifier.
#'   - Cluster: The assigned cluster number for each sample.
#'   - Cluster2: The assigned cluster label, prefixed by 'MCIA' to indicate that the clustering was performed using MCIA.
#' @details This function uses MCIA to integrate multiple data matrices and then performs hierarchical clustering on the resulting components to identify clusters. MCIA is useful for identifying shared structures across multiple data modalities.
#'
#' The function operates as follows:
#' 1. Each matrix in the input list is transposed so that rows represent samples and columns represent features.
#' 2. MCIA is performed using the omicade4 package to extract components that summarize the shared variation across different modalities.
#' 3. Hierarchical clustering is applied to the concatenated components to assign each sample to a cluster.
#' 4. The function returns a data frame containing the cluster assignment for each sample, along with additional information about the clustering process.
#' @references Meng C, Kuster B, Culhane A, Gholami AM. A multivariate approach to the integration of multi-omics datasets. BMC Bioinformatics. 2013.
#' @examples
#' # Example usage:
#' data1 <- matrix(rnorm(10000), nrow = 100, ncol = 100)
#' data2 <- matrix(rnorm(10000), nrow = 100, ncol = 100)
#' colnames(data1) <- colnames(data2) <- paste0("Sample", 1:100)
#' data_list <- list(data1, data2)
#'
#' # Run MCIA clustering
#' result <- RunMCIA(data = data_list, N.clust = 3)
#'
#' @export
RunMCIA <- function(data = NULL, N.clust = NULL, n.components = 10, clustering.algorithm = "ward.D2",
                    scan.eigenvalues = FALSE, use.nsc = TRUE, use.svd = TRUE) {
  # Set seed for reproducibility
  set.seed(1)

  # Check if required packages are installed
  if (!requireNamespace("omicade4", quietly = TRUE)) {
    if (!requireNamespace("BiocManager", quietly = TRUE)) {
      install.packages("BiocManager")
    }
    BiocManager::install("omicade4")
  }

  # Input validation checks
  if (!is.list(data)) {
    stop("The 'data' parameter must be a list of matrices.")
  }

  # Perform MCIA using omicade4
  mcoin <- omicade4::mcia(
    df.list = data, cia.nf = n.components, cia.scan = scan.eigenvalues,
    nsc = use.nsc, svd = use.svd
  )

  # Perform hierarchical clustering on the MCIA components
  clust <- stats::hclust(stats::dist(mcoin$mcoa$SynVar), method = clustering.algorithm) %>%
    stats::cutree(k = N.clust)

  # Create a data frame with cluster assignments
  clustres <- data.frame(
    Sample = colnames(data[[1]]),
    Cluster = clust,
    Cluster2 = paste0("MCIA", clust),
    row.names = colnames(data[[1]]),
    stringsAsFactors = FALSE
  )

  # Print summary of clustering results
  cat(paste0("Cluster algorithm: MCIA\n"))
  if (!is.null(N.clust)) {
    cat(paste0("Selection of cluster number: ", N.clust, "\n"))
    cluster_summary <- paste(
      paste0("C", 1:N.clust),
      round(table(clustres$Cluster) / nrow(clustres), 2) * 100,
      sep = " = "
    )
    cat(paste(cluster_summary, collapse = "% "), "%\n")
  } else {
    cat("Number of clusters was determined automatically by MCIA.\n")
  }

  return(clustres)
}
