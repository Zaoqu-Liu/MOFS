#' @title Jaccard Distance Calculation for Binary Matrix
#' @description This function calculates the Jaccard distance or similarity for a binary matrix. It is typically used to evaluate the similarity or dissimilarity between columns of a binary matrix.
#' @author Zaoqu Liu; Email: liuzaoqu@163.com
#' @param data A binary matrix where rows represent features and columns represent samples.
#' @param dissimilarity Logical. If TRUE, returns the Jaccard distance; if FALSE, returns the Jaccard similarity (default: TRUE).
#' @return A matrix containing the Jaccard distance or similarity between each pair of columns in the input matrix.
#' @details The function computes the Jaccard distance (or similarity) between each pair of columns in the input binary matrix. The Jaccard distance is calculated as 1 minus the Jaccard similarity.
#' @examples
#' data <- matrix(sample(0:1, 1500, replace = TRUE), nrow = 30, ncol = 50)
#' jaccard_dist <- get.Jaccard.Distance(as.data.frame(data), dissimilarity = TRUE)
#' jaccard_dist
#'
#' @export
get.Jaccard.Distance <- function(data, dissimilarity = TRUE) {
  # Initialize the Jaccard matrix
  data <- as.data.frame(data)
  jaccard_matrix <- matrix(0, nrow = ncol(data), ncol = ncol(data))
  rownames(jaccard_matrix) <- colnames(jaccard_matrix) <- colnames(data)

  # Calculate Jaccard similarity for each pair of columns
  for (i in 1:ncol(data)) {
    for (j in 1:ncol(data)) {
      a <- rownames(data[data[, i] != 0, i, drop = FALSE])
      b <- rownames(data[data[, j] != 0, j, drop = FALSE])
      jaccard_matrix[i, j] <- length(intersect(a, b)) / length(union(a, b))
    }
  }

  # Convert to Jaccard distance if specified
  if (dissimilarity) {
    jaccard_matrix <- 1 - jaccard_matrix
  }

  return(jaccard_matrix)
}
