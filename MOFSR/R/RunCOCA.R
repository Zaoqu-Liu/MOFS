#' @title Run Consensus Clustering Analysis (COCA)
#' @description This function performs Consensus Clustering Analysis (COCA) using the ConsensusClusterPlus package to identify stable clusters in the input data.
#' @author Zaoqu Liu; Email: liuzaoqu@163.com
#' @param jaccard.matrix A Jaccard distance matrix, typically obtained from binary data.
#' @param max.clusters Integer. The maximum number of clusters to evaluate (default: 6).
#' @param optimal.clusters Integer. The optimal number of clusters to select (default: 3).
#' @param linkage.method Character. The linkage method for hierarchical clustering (default: "ward.D2").
#' @param clustering.algorithm Character. The clustering algorithm to use (default: 'hc').
#' @param distance.metric Character. The distance metric to use (default: "euclidean").
#' @param resampling.iterations Integer. The number of resampling iterations (default: 10000).
#' @param resample.proportion Numeric. Proportion of items to resample in each iteration (default: 0.7).
#' @return A list containing the consensus clustering results, optimal cluster solution, PAC values, and final cluster assignments.
#' @details This function uses ConsensusClusterPlus to perform consensus clustering on the input Jaccard distance matrix, evaluates the stability of different clustering solutions using PAC, and returns the clustering assignments.
#' @examples
#' # Example usage:
#' jaccard_matrix <- CalJaccardDistance(data)
#' coca_result <- RunCOCA(jaccard.matrix = jaccard_matrix, max.clusters = 6, optimal.clusters = 3)
#'
#' @export
RunCOCA <- function(jaccard.matrix, max.clusters = 6, optimal.clusters = 3,
                    linkage.method = "ward.D2", clustering.algorithm = "hc", distance.metric = "euclidean",
                    resampling.iterations = 10000, resample.proportion = 0.7) {
  # Set seed for reproducibility
  set.seed(1)

  # Check if required packages are installed
  if (!requireNamespace("ConsensusClusterPlus", quietly = TRUE)) {
    if (!requireNamespace("BiocManager", quietly = TRUE)) {
      install.packages("BiocManager")
    }
    BiocManager::install("ConsensusClusterPlus")
  }

  # Perform Consensus Clustering
  consensus.results <- zquiet(ConsensusClusterPlus::ConsensusClusterPlus(
    d = as.dist(jaccard.matrix),
    maxK = max.clusters,
    reps = resampling.iterations, pItem = resample.proportion,
    innerLinkage = linkage.method, finalLinkage = linkage.method,
    clusterAlg = clustering.algorithm, distance = distance.metric,
    plot = NULL,
    seed = 123, verbose = FALSE
  ))

  # Extract optimal clustering solution
  optimal.solution <- consensus.results[[optimal.clusters]]

  # Calculate PAC values
  PAC <- CalPAC(consensus.results, range_clusters = 2:max.clusters)

  # Extract final cluster assignments
  cluster.assignments <- consensus.results[[optimal.clusters]][["consensusClass"]]
  cluster.results <- data.frame(Sample = names(cluster.assignments), Cluster = cluster.assignments)

  # Print summary of clustering results
  cat(paste0("Cluster algorithm: COCA\nSelection of cluster number: ", optimal.clusters, "\n"))
  cluster_summary <- paste(
    paste0("C", 1:optimal.clusters),
    round(table(cluster.results$Cluster) / nrow(cluster.results), 2) * 100,
    sep = " = "
  )
  cat(paste(cluster_summary, collapse = "% "), "%\n")

  return(list(fit = consensus.results, optimal = optimal.solution, PAC = PAC, Cluster = cluster.results))
}
