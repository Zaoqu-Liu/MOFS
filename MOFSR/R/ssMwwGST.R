#' @title Single-Sample Pathway Activity Analysis Using ssMwwGST
#' @description This function performs single-sample pathway activity analysis using the MWW-GST (Mann-Whitney-Wilcoxon Gene Set Test) method. It assesses pathway activity for each individual sample, leveraging the yaGST package to estimate pathway scores for each gene set.
#' @author Zaoqu Liu; Email: liuzaoqu@163.com
#' @param geData A numeric matrix of gene expression data. Rows represent genes, and columns represent samples.
#' @param geneSet A list of gene sets where each element is a vector of gene names included in the set. Each gene set should have at least 15 genes.
#' @param nCores Integer. The number of cores to use for parallel processing. Default is 8.
#' @return A list containing:
#'   - NES: A numeric matrix of Normalized Enrichment Scores (NES) for each gene set and each sample.
#'   - pValue: A numeric matrix of p-values for the enrichment of each gene set in each sample.
#'   - FDR: A numeric matrix of False Discovery Rate (FDR) adjusted p-values for each gene set in each sample.
#' @details The function operates as follows:
#' 1. Calculates the mean and standard deviation for each gene across samples.
#' 2. For each sample, the gene expression is standardized (z-score normalization).
#' 3. Uses the Mann-Whitney-Wilcoxon Gene Set Test (MWW-GST) to assess pathway enrichment for each gene set in each sample.
#' 4. Performs multiple testing correction (FDR) for p-values.
#' @references Frattini V, Pagnotta SM, Tala, et al. A metabolic function of FGFR3-TACC3 gene fusions in cancer. Nature. 2018;553(7687):222-227. doi:10.1038/nature25171
#' @examples
#' # Example usage:
#' geData <- matrix(rnorm(1000), nrow = 100, ncol = 10) # Generate random gene expression data
#' geneSet <- list(set1 = sample(rownames(geData), 20), set2 = sample(rownames(geData), 15))
#' # Run ssMwwGST pathway activity analysis
#' result <- ssMwwGST(geData, geneSet, nCores = 4)
#'
#' @export
ssMwwGST <- function(geData, geneSet, nCores = 8) {
  # Load required packages
  if (!requireNamespace("yaGST", quietly = TRUE)) {
    message("Installing yaGST package from GitHub...")
    devtools::install_github("miccec/yaGST")
  }
  if (!requireNamespace("doMC", quietly = TRUE)) {
    message("Installing doMC package from CRAN...")
    install.packages("doMC")
  }
  if (!requireNamespace("foreach", quietly = TRUE)) {
    message("Installing foreach package from CRAN...")
    install.packages("foreach")
  }

  # Load necessary libraries
  library(yaGST)
  library(doMC)
  library(foreach)

  # Calculate mean and standard deviation for each gene across samples
  means <- rowMeans(geData)
  sds <- apply(geData, 1, sd)

  # Register parallel backend
  doMC::registerDoMC(nCores)

  # Initialize list to store results for each sample
  ans <- foreach::foreach(ss = 1:ncol(geData), .packages = c("yaGST")) %dopar% {
    # Standardize gene expression for the current sample (z-score normalization)
    currentSample <- (geData[, ss] - means) / sds
    rankedList <- sort(currentSample, decreasing = TRUE)

    # Apply MWW-GST for each gene set
    aMwwGST <- lapply(geneSet, function(x) {
      if (length(x) >= 15) {
        yaGST::mwwGST(rankedList, geneSet = x, minLenGeneSet = 15, alternative = "two.sided", verbose = FALSE)
      } else {
        NULL
      }
    })
    aMwwGST <- aMwwGST[sapply(aMwwGST, length) != 0] # Remove empty results

    # Extract NES and p-values
    tmp_NES <- sapply(aMwwGST, function(x) x$log.pu)
    tmp_pValue <- sapply(aMwwGST, function(x) x$p.value)

    # Return results for the current sample
    list(tmp_NES = tmp_NES, tmp_pValue = tmp_pValue)
  }

  # Extract results for all samples
  NES <- sapply(ans, function(x) x$tmp_NES)
  pValue <- sapply(ans, function(x) x$tmp_pValue)
  colnames(NES) <- colnames(pValue) <- colnames(geData)

  # Calculate FDR for each gene set across all samples
  FDR <- t(apply(pValue, 1, function(x) p.adjust(x, method = "fdr")))

  # Create a result list containing NES, p-values, and FDR-adjusted p-values
  res <- list(NES = NES, pValue = pValue, FDR = FDR)

  return(res)
}
