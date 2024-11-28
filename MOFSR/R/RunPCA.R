#' @title Run PCA on Data
#' @description This function performs Principal Component Analysis (PCA) on the given data using the FactoMineR package.
#' @author Zaoqu Liu; Email: liuzaoqu@163.com
#' @param data A distance matrix or data frame on which to perform PCA.
#' @return The PCA result as produced by the FactoMineR package.
#' @details This function uses the FactoMineR package to perform PCA, which can be used to visualize the relationships between samples in reduced dimensional space.
#' @examples
#' # Example usage:
#' data <- matrix(rnorm(100), nrow = 10, ncol = 10)
#' pca_result <- RunPCA(data)
#'
#' @export
RunPCA <- function(data) {
  if (!requireNamespace("FactoMineR", quietly = TRUE)) {
    install.packages("FactoMineR")
  }
  ddb.pca <- FactoMineR::PCA(data, graph = FALSE)
  return(ddb.pca)
}
