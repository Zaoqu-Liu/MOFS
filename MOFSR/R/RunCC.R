#' @title Run Consensus Clustering using ConsensusClusterPlus
#' @description This function performs consensus clustering on the given data using the ConsensusClusterPlus package.
#' @author Zaoqu Liu; Email: liuzaoqu@163.com
#' @param data A numeric matrix or data frame where rows represent features and columns represent samples.
#' @param maxK The maximum number of clusters to evaluate (default: 6).
#' @param reps Number of subsamples (default: 1000).
#' @param pItem Proportion of items to sample (default: 0.8).
#' @param pFeature Proportion of features to sample (default: 1).
#' @param clusterAlg The clustering algorithm to use, either "hc" for hierarchical or "km" for k-means (default: "hc").
#' @param distance The distance metric to use, "pearson", "spearman", "euclidean", etc. (default: "euclidean").
#' @param title Optional title for the results (default: "Consensus Clustering").
#' @param plot Whether to plot the consensus matrix and dendrogram (default: TRUE).
#' @return A list containing the consensus clustering results.
#' @details This function leverages the ConsensusClusterPlus package to perform consensus clustering, which is useful for identifying robust clusters in genomic data or other high-dimensional data.
#' @examples
#' # Example usage:
#' data <- matrix(rnorm(1000), nrow = 100, ncol = 10)
#' result <- RunCC(data, maxK = 4)
#' @export
RunCC <- function(data, maxK = 6, reps = 1000, pItem = 0.8, pFeature = 1,
                  clusterAlg = "hc", distance = "euclidean", title = "Consensus Clustering", plot = TRUE) {
  # Load required package
  if (!requireNamespace("ConsensusClusterPlus", quietly = TRUE)) {
    install.packages("ConsensusClusterPlus")
  }
  library(ConsensusClusterPlus)

  # Perform consensus clustering
  set.seed(1234)
  cc_result <- ConsensusClusterPlus(t(data),
    maxK = maxK, reps = reps, pItem = pItem,
    pFeature = pFeature, clusterAlg = clusterAlg,
    distance = distance, title = title, plot = plot
  )

  return(cc_result)
}
