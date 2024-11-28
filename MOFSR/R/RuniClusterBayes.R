#' @title Run Bayesian iCluster (iClusterBayes) for Multi-Modality Data Integration
#' @description This function performs clustering analysis using Bayesian iCluster (iClusterBayes) to integrate multiple omics datasets. iClusterBayes is an extension of iCluster that allows for Bayesian inference, which is useful for identifying shared patterns across multiple datasets with different distributions.
#' @author Zaoqu Liu; Email: liuzaoqu@163.com
#' @param data A list of matrices where each element represents a different modality (e.g., RNA, protein, methylation). Each matrix should have rows as features and columns as samples.
#' @param N.clust Integer. Number of clusters to create (optional but recommended).
#' @param data.type A character vector specifying the type of data for each modality. Options include "binomial", "gaussian", etc. Default is c("binomial", "gaussian", "gaussian", "gaussian", "gaussian", "gaussian").
#' @param num.burnin Integer. Number of burn-in iterations for the Bayesian algorithm (default: 18000).
#' @param num.draws Integer. Number of MCMC draws after burn-in (default: 12000).
#' @param prior.probabilities Numeric vector. Prior values for the inclusion probabilities for each modality (default: 0.5 for each modality).
#' @param proposal.sdev Numeric. Standard deviation of the proposal distribution (default: 0.05).
#' @param thinning.interval Integer. Thinning interval for MCMC sampling (default: 3).
#' @return A data frame with the following columns:
#'   - Sample: The sample identifier.
#'   - Cluster: The assigned cluster number for each sample.
#'   - Cluster2: The assigned cluster label, prefixed by 'iClusterBayes' to indicate that the clustering was performed using iClusterBayes.
#' @details This function uses iClusterBayes to integrate multiple data matrices and assign clusters to the samples. iClusterBayes performs Bayesian inference, which is useful for modeling different data distributions across multiple datasets.
#'
#' The function operates as follows:
#' 1. Each matrix in the input list is transposed so that rows represent samples and columns represent features.
#' 2. iClusterBayes is used to identify shared patterns across the different modalities.
#' 3. The function returns a data frame containing the cluster assignment for each sample, along with additional information about the clustering process.
#' @references
#' Mo Q, Shen R, et al. (2013) Pattern discovery and cancer gene identification in integrated cancer genomic data. PNAS, 110 (11) 4245-4250.
#' Shen R, Olshen AB, et al. (2012) Integrative Subtype Discovery in Glioblastoma Using iCluster, PLOS ONE 7(4):e35236.
#' Shen R, Olshen AB, et al. (2009) Integrative clustering of multiple genomic data types using a joint latent variable model with application to breast and lung cancer subtype analysis. Bioinformatics, 25(22):2906-12.
#' @examples
#' # Example usage:
#' data1 <- matrix(rnorm(10000), nrow = 100, ncol = 100)
#' data2 <- matrix(rnorm(10000), nrow = 100, ncol = 100)
#' colnames(data1) <- colnames(data2) <- paste0("Sample", 1:100)
#' data_list <- list(data1, data2)
#'
#' # Run iClusterBayes clustering
#' result <- RuniClusterBayes(
#'   data = data_list, N.clust = 3,
#'   data.type = c("gaussian", "gaussian")
#' )
#'
#' @export
RuniClusterBayes <- function(data = NULL, N.clust = NULL,
                             data.type = c("binomial", "gaussian", "gaussian", "gaussian", "gaussian", "gaussian"),
                             num.burnin = 20, num.draws = 10, prior.probabilities = rep(0.5, length(data)),
                             proposal.sdev = 0.05, thinning.interval = 3) {
  # Set seed for reproducibility
  set.seed(1)

  # Check if required packages are installed
  if (!requireNamespace("iClusterPlus", quietly = TRUE)) {
    if (!requireNamespace("BiocManager", quietly = TRUE)) {
      install.packages("BiocManager")
    }
    BiocManager::install("iClusterPlus")
  }

  # Input validation checks
  if (!is.list(data)) {
    stop("The 'data' parameter must be a list of matrices.")
  }

  # Transpose each matrix to have rows as samples and columns as features
  data <- lapply(data, t)

  # Perform iClusterBayes clustering based on the number of data types
  n_dat <- length(data)
  res <- switch(n_dat,
    `1` = iClusterPlus::iClusterBayes(
      dt1 = data[[1]], type = data.type, K = N.clust - 1,
      n.burnin = num.burnin, n.draw = num.draws, prior.gamma = prior.probabilities,
      sdev = proposal.sdev, thin = thinning.interval
    ),
    `2` = iClusterPlus::iClusterBayes(
      dt1 = data[[1]], dt2 = data[[2]],
      type = data.type, K = N.clust - 1, n.burnin = num.burnin,
      n.draw = num.draws, prior.gamma = prior.probabilities, sdev = proposal.sdev,
      thin = thinning.interval
    ),
    `3` = iClusterPlus::iClusterBayes(
      dt1 = data[[1]], dt2 = data[[2]],
      dt3 = data[[3]], type = data.type, K = N.clust - 1, n.burnin = num.burnin,
      n.draw = num.draws, prior.gamma = prior.probabilities, sdev = proposal.sdev,
      thin = thinning.interval
    ),
    `4` = iClusterPlus::iClusterBayes(
      dt1 = data[[1]], dt2 = data[[2]],
      dt3 = data[[3]], dt4 = data[[4]], type = data.type, K = N.clust - 1,
      n.burnin = num.burnin, n.draw = num.draws, prior.gamma = prior.probabilities,
      sdev = proposal.sdev, thin = thinning.interval
    ),
    `5` = iClusterPlus::iClusterBayes(
      dt1 = data[[1]], dt2 = data[[2]],
      dt3 = data[[3]], dt4 = data[[4]], dt5 = data[[5]],
      type = data.type, K = N.clust - 1, n.burnin = num.burnin,
      n.draw = num.draws, prior.gamma = prior.probabilities, sdev = proposal.sdev,
      thin = thinning.interval
    ),
    `6` = iClusterPlus::iClusterBayes(
      dt1 = data[[1]], dt2 = data[[2]],
      dt3 = data[[3]], dt4 = data[[4]], dt5 = data[[5]],
      dt6 = data[[6]], type = data.type, K = N.clust - 1, n.burnin = num.burnin,
      n.draw = num.draws, prior.gamma = prior.probabilities, sdev = proposal.sdev,
      thin = thinning.interval
    ),
    stop("Unsupported number of data types. Maximum is 6.")
  )

  # Create a data frame with cluster assignments
  clustres <- data.frame(
    Sample = rownames(data[[1]]),
    Cluster = res$clusters,
    Cluster2 = paste0("iClusterBayes", res$clusters),
    row.names = rownames(data[[1]]),
    stringsAsFactors = FALSE
  )

  # Print summary of clustering results
  cat(paste0("Cluster algorithm: iClusterBayes\n"))
  if (!is.null(N.clust)) {
    cat(paste0("Selection of cluster number: ", N.clust, "\n"))
    cluster_summary <- paste(
      paste0("C", 1:N.clust),
      round(table(clustres$Cluster) / nrow(clustres), 2) * 100,
      sep = " = "
    )
    cat(paste(cluster_summary, collapse = "% "), "%\n")
  } else {
    cat("Number of clusters was determined automatically by iClusterBayes.\n")
  }

  return(clustres)
}
