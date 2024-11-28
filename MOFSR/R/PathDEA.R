#' @title Pathway Differential Expression Analysis (PathDEA)
#' @description This function performs pathway differential expression analysis based on clustering results and pathway activity scores derived from ssMwwGST.
#' @author Zaoqu Liu; Email: liuzaoqu@163.com
#' @param Cluster_data A data frame where the first column must be the sample IDs and the second column must be the cluster assignments.
#' @param ssMwwGST_results A list of results from ssMwwGST, including NES (Normalized Enrichment Scores).
#' @param dea_FDR_threshold Numeric. The FDR threshold to use for filtering significant pathways. Default is 0.001.
#' @param dea_gap_threshold Numeric. The median gap threshold to use for filtering significant pathways. Default is 1.5.
#' @return A list containing:
#'   - dea_path: A list of data frames, each containing the differential expression analysis results for each cluster.
#'   - dea_path2: A list of data frames containing the filtered differential expression analysis results for each cluster based on FDR and median gap thresholds.
#'   - NES: A data frame of Normalized Enrichment Scores for each gene set and each sample.
#'   - Cluster: A data frame of sample IDs and their corresponding cluster assignments.
#' @details The function operates as follows:
#' 1. Extracts Normalized Enrichment Scores (NES) from the ssMwwGST results.
#' 2. Performs Wilcoxon rank-sum tests to compare pathway activity between clusters for each pathway.
#' 3. Calculates median and mean differences in pathway activity between clusters.
#' 4. Adjusts p-values using the Benjamini-Hochberg method to control the false discovery rate (FDR).
#' @examples
#' # Example usage:
#' Cluster_data <- data.frame(Sample = paste0("Sample", 1:10), Cluster = rep(1:2, each = 5))
#' ssMwwGST_results <- list(NES = matrix(rnorm(100),
#'   nrow = 10, ncol = 10,
#'   dimnames = list(paste0("Pathway", 1:10), paste0("Sample", 1:10))
#' ))
#' result <- PathDEA(Cluster_data, ssMwwGST_results)
#'
#' @export
PathDEA <- function(Cluster_data, ssMwwGST_results, dea_FDR_threshold = 0.001, dea_gap_threshold = 1.5) {
  # Set column names for Cluster_data
  colnames(Cluster_data) <- c("Sample", "Cluster")

  # Extract unique cluster labels
  Cluster_vec <- unique(Cluster_data$Cluster)

  # Extract Normalized Enrichment Scores (NES) from ssMwwGST results
  NES <- as.data.frame(ssMwwGST_results$NES)
  NES <- NES[, Cluster_data$Sample]

  # Perform differential expression analysis for each cluster
  dea_path <- purrr::map(Cluster_vec, function(x) {
    tt <- Cluster_data$Sample[Cluster_data$Cluster == x]
    cc <- Cluster_data$Sample[Cluster_data$Cluster != x]
    res <- apply(NES, 1, function(y) {
      names(y) <- colnames(NES)
      fit <- wilcox.test(y[tt], y[cc])
      median_gap <- median(y[tt]) - median(y[cc])
      mean_gap <- mean(y[tt]) - mean(y[cc])
      P <- fit$p.value
      Statistic <- as.numeric(fit$statistic)
      return(data.frame(
        Cluster = x, P = P, Statistic = Statistic,
        Median_gap = median_gap, Mean_gap = mean_gap
      ))
    }) %>% Reduce(rbind, .)
    res$FDR <- p.adjust(res$P, "BH")
    res$Pathway <- rownames(NES)
    regexp <- paste0(paste0(unique(stringr::str_match(res$Pathway, "([A-Z]{1,})_.*")[, 2]), "_"), collapse = "|")
    res$Pathway2 <- gsub(regexp, "", res$Pathway)
    res$Pathway2 <- gsub("_", " ", res$Pathway2)
    res$Pathway2 <- tolower(res$Pathway2)
    res$Pathway2 <- Hmisc::capitalize(res$Pathway2)
    res <- dplyr::select(res, Cluster, Pathway2, FDR, Median_gap, Mean_gap, everything())
    return(res)
  })

  # Name the list of differential expression analysis results by cluster
  names(dea_path) <- Cluster_vec

  # Filter significant pathways based on FDR and median gap thresholds
  dea_path2 <- purrr::map(dea_path, function(x) {
    x <- x[x$Median_gap >= dea_gap_threshold & x$FDR < dea_FDR_threshold, ]
    return(x)
  })

  # Return results
  return(list(dea_path = dea_path, dea_path2 = dea_path2, NES = NES, Cluster = Cluster_data))
}
