#' @title Run Consensus Principal Component Analysis (CPCA) for Multi-Modality Data Integration
#' @description This function performs clustering analysis using Consensus Principal Component Analysis (CPCA) to integrate multiple omics data. CPCA helps capture the shared variation across different data modalities and provides insights into the common structure of the data.
#' @author Zaoqu Liu; Email: liuzaoqu@163.com
#' @param data A list of matrices where each element represents a different modality (e.g., RNA, protein, methylation). Each matrix should have rows as features and columns as samples.
#' @param N.clust Integer. Number of clusters to create from the hierarchical clustering of the CPCA components (optional but recommended).
#' @param num.components Integer. Number of components to retain for each modality (default: 2).
#' @param integration.algorithm Character. Algorithm to use for CPCA, options include "globalScore" (default), "blockScore", or "blockLoading".
#' @param nonzero.coeff.k Numeric or Character. The number (if >= 1) or proportion (if 0 < k < 1) of non-zero coefficients for variable loadings (default: "all").
#' @param center.data Logical. Whether to center each block to zero mean (default: FALSE).
#' @param scale.data Logical. Whether to scale each block to unit variance (default: FALSE).
#' @param normalization.option Character. Normalization option, one of "lambda1", "inertia", or "uniform" (default: "uniform").
#' @param max.iterations Integer. Maximum number of iterations (default: 1000).
#' @param return.moa.object Logical. Whether to return an object of class `moa-class` (default: TRUE).
#' @param show.verbose Logical. Whether to print process information (default: FALSE).
#' @param svd.solver.method Character. SVD solver to use, one of "svd", "fast.svd", or "propack" (default: "fast.svd").
#' @param nonzero.coeff.obs Numeric or Character. Number or proportion of non-zero coefficients for observation scores (default: "all").
#' @param weight.variables Numeric, Vector, or List. Weights for variables (default: NA).
#' @param weight.observations Numeric, Vector, or List. Weights for observations (default: NA).
#' @param unit.length.variables Logical. Whether the loading vectors for each block should have unit length (default: FALSE).
#' @param unit.length.observations Logical. Whether the score vectors for each block should have unit length (default: FALSE).
#' @param retain.nonnegative Logical. Whether to retain only non-negative coefficients in loadings and scores (default: FALSE).
#' @param clustering.algorithm Character. The clustering algorithm to use for hierarchical clustering (default: "ward.D2").
#' @return A data frame with the following columns:
#'   - Sample: The sample identifier.
#'   - Cluster: The assigned cluster number for each sample.
#'   - Cluster2: The assigned cluster label, prefixed by 'CPCA' to indicate that the clustering was performed using CPCA.
#' @details This function uses CPCA to integrate multiple data matrices and then performs hierarchical clustering on the resulting components to identify clusters.
#' @references Meng C, Basunia A, Peters B, Gholami AM, Kuster B, Culhane AC. MOGSA: Integrative Single Sample Gene-set Analysis of Multiple Omics Data. Mol Cell Proteomics. 2019;18(8 suppl 1):S153-S168. doi:10.1074/mcp.TIR118.001251.
#' @examples
#' # Example usage:
#' data1 <- matrix(rnorm(10000), nrow = 100, ncol = 100)
#' data2 <- matrix(rnorm(10000), nrow = 100, ncol = 100)
#' colnames(data1) <- colnames(data2) <- paste0("Sample", 1:100)
#' data_list <- list(data1, data2)
#'
#' # Run CPCA clustering
#' result <- RunCPCA(data = data_list, N.clust = 3)
#'
#' @export
RunCPCA <- function(data = NULL, N.clust = NULL, num.components = 2,
                    integration.algorithm = "globalScore", nonzero.coeff.k = "all",
                    center.data = FALSE, scale.data = FALSE, normalization.option = "uniform",
                    max.iterations = 1000, return.moa.object = TRUE, show.verbose = FALSE,
                    svd.solver.method = "fast.svd", nonzero.coeff.obs = "all",
                    weight.variables = NA, weight.observations = NA,
                    unit.length.variables = FALSE, unit.length.observations = FALSE,
                    retain.nonnegative = FALSE, clustering.algorithm = "ward.D2") {
  # Set seed for reproducibility
  set.seed(1)

  # Check if required packages are installed
  if (!requireNamespace("mogsa", quietly = TRUE)) {
    if (!requireNamespace("BiocManager", quietly = TRUE)) {
      install.packages("BiocManager")
    }
    BiocManager::install("mogsa")
  }

  # Input validation checks
  if (!is.list(data)) {
    stop("The 'data' parameter must be a list of matrices.")
  }

  # Perform CPCA using mogsa
  moas <- mogsa::mbpca(
    x = data, ncomp = num.components, method = integration.algorithm, k = nonzero.coeff.k,
    center = center.data, scale = scale.data, option = normalization.option,
    maxiter = max.iterations, moa = return.moa.object, verbose = show.verbose,
    svd.solver = svd.solver.method, k.obs = nonzero.coeff.obs, w = weight.variables,
    w.obs = weight.observations, unit.p = unit.length.variables,
    unit.obs = unit.length.observations, pos = retain.nonnegative
  )

  # Extract the scores for each sample
  scrs <- mogsa::moaScore(moas)

  # Perform hierarchical clustering on the extracted scores
  dist <- stats::dist(scrs)
  clust.dend <- stats::hclust(dist, method = clustering.algorithm)
  clustres <- data.frame(
    Sample = colnames(data[[1]]),
    Cluster = stats::cutree(clust.dend, k = N.clust),
    Cluster2 = paste0("CPCA", stats::cutree(clust.dend, k = N.clust)),
    row.names = colnames(data[[1]]), stringsAsFactors = FALSE
  )

  # Print summary of clustering results
  cat(paste0("Cluster algorithm: CPCA\n"))
  if (!is.null(N.clust)) {
    cat(paste0("Selection of cluster number: ", N.clust, "\n"))
    cluster_summary <- paste(
      paste0("C", 1:N.clust),
      round(table(clustres$Cluster) / nrow(clustres), 2) * 100,
      sep = " = "
    )
    cat(paste(cluster_summary, collapse = "% "), "%\n")
  } else {
    cat("Number of clusters was determined automatically by CPCA.\n")
  }

  return(clustres)
}
