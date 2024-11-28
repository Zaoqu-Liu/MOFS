#' @title Calculate Proportion of Ambiguous Clustering (PAC)
#' @description This function calculates the Proportion of Ambiguous Clustering (PAC) to help evaluate the optimal number of clusters in a consensus clustering analysis.
#' @author Zaoqu Liu; Email: liuzaoqu@163.com
#' @param consensus_result A list containing consensus clustering results from ConsensusClusterPlus.
#' @param range_clusters Integer vector. The range of cluster numbers evaluated during consensus clustering.
#' @param x1 Numeric. Lower bound for defining the PAC (default: 0.1).
#' @param x2 Numeric. Upper bound for defining the PAC (default: 0.9).
#' @return A data frame containing the PAC values for each number of clusters.
#' @details The PAC is calculated as the difference between the cumulative distribution function values at two thresholds, x2 and x1, on the consensus matrix values. A lower PAC indicates a more stable clustering solution.
#' @examples
#' data <- mtcars
#' cc_res <- RunCC(data)
#' pac_values <- CalPAC(cc_res)
#' pac_values
#' @export
CalPAC <- function(consensus_result, range_clusters = 2:6, x1 = 0.1, x2 = 0.9) {
  PAC <- rep(NA, length(range_clusters))
  names(PAC) <- paste("K=", range_clusters, sep = "")
  for (i in range_clusters) {
    consensus_matrix <- consensus_result[[i]]$consensusMatrix
    Fn <- ecdf(consensus_matrix[lower.tri(consensus_matrix)])
    PAC[i - 1] <- Fn(x2) - Fn(x1)
  }
  PAC <- as.data.frame(PAC)
  PAC$K <- range_clusters
  return(PAC)
}
