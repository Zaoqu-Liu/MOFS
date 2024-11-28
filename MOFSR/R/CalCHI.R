#' @title Calinski-Harabasz Index Calculation
#' @description This function calculates the Calinski-Harabasz index for evaluating the quality of clustering solutions. It is used to determine the optimal number of clusters in a hierarchical clustering scenario.
#' @author Zaoqu Liu; Email: liuzaoqu@163.com
#' @param hclust_result A hierarchical clustering object (result of hclust function).
#' @param dist_matrix An optional distance matrix. If not provided, the cophenetic matrix is used. The matrix should be of type "euclidean".
#' @param max_clusters Integer. The maximum number of clusters to evaluate (default: round(1 + 3.3 * log(length(hclust_result$order), 10))).
#' @return A vector containing the Calinski-Harabasz index values for each number of clusters.
#' @details The function calculates the Calinski-Harabasz index for different numbers of clusters to help identify the optimal number of clusters. The index measures the ratio of between-cluster dispersion to within-cluster dispersion.
#' @examples
#' # Example usage:
#' data <- mtcars
#' dist_matrix <- dist(data)
#' hclust_result <- hclust(dist_matrix)
#' calinski_values <- CalCHI(hclust_result, dist_matrix)
#'
#' @export
CalCHI <- function(hclust_result, dist_matrix = NULL, max_clusters = round(1 + 3.3 * log(length(hclust_result$order), 10))) {
  msg <- ""

  # If clue package is not installed, install it from CRAN
  if (!requireNamespace("clue", quietly = TRUE)) {
    install.packages("clue")
  }

  # If dist_matrix is not provided or not "euclidean", use cophenetic matrix
  if (is.null(dist_matrix)) {
    dist_matrix <- sqrt(clue::as.cl_ultrametric(hclust_result))
  } else if (attr(dist_matrix, "method") != "euclidean") {
    dist_matrix <- sqrt(clue::as.cl_ultrametric(hclust_result))
  }

  # Convert distance matrix to squared matrix
  dist_matrix <- as.matrix(dist_matrix)^2

  # Calculate A matrix
  A <- -dist_matrix / 2
  A_bar <- apply(A, 1, mean)
  total_sum <- sum(diag(A) - 2 * A_bar + mean(A))
  n_samples <- length(hclust_result$order)

  # Initialize vector to store Calinski-Harabasz values
  calinski_values <- rep(0, max_clusters)

  # Loop over possible number of clusters
  for (g in 2:max_clusters) {
    cluster_assignment <- cutree(hclust_result, k = g)
    within_sum <- 0

    # Calculate within-cluster sum of squares
    for (k in 1:g) {
      if (sum(cluster_assignment == k) == 1) {
        next
      }
      A_cluster <- as.matrix(-dist_matrix / 2)[cluster_assignment == k, cluster_assignment == k]
      A_bar_cluster <- apply(A_cluster, 1, mean)
      within_sum <- within_sum + sum(diag(A_cluster) - 2 * A_bar_cluster + mean(A_cluster))
    }

    # Calculate between-cluster sum of squares
    between_sum <- total_sum - within_sum
    between_sum <- between_sum / (g - 1)
    within_sum <- within_sum / (n_samples - g)

    # Store Calinski-Harabasz index value
    calinski_values[g] <- between_sum / within_sum
  }

  # Assign class and message attribute to result
  class(calinski_values) <- "calinski_harabasz"
  attr(calinski_values, "message") <- msg

  return(calinski_values)
}
