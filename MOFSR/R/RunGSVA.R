#' @title Generate Single-Sample Gene-Set Enrichment Score
#' @description This function estimates gene-set enrichment scores across all samples using various methods.
#' @details This function supports multiple methods for estimating gene-set enrichment scores, including ssGSEA, GSVA, zscore, and plage. The scores are calculated for each gene set across all samples.
#' The `ssGSES` function is flexible and allows for customization of the minimum and maximum size of gene sets considered in the analysis. By providing different methods, the function can adapt to various types of gene-set enrichment analysis, each having its own strengths and suitable applications.
#' \itemize{
#'   \item \code{"gsva"}: Gene Set Variation Analysis, suitable for detecting subtle changes in pathway activity.
#'   \item \code{"ssgsea"}: Single-Sample Gene Set Enrichment Analysis, useful for individual sample analysis.
#'   \item \code{"zscore"}: Z-score transformation, a simpler approach to standardize expression values.
#'   \item \code{"plage"}: Pathway Level Analysis of Gene Expression, which focuses on correlating pathway components.
#' }
#' @author
#' Zaoqu Liu; Email: liuzaoqu@163.com
#' @param exp Numeric matrix containing the expression data or gene expression signatures, with samples in columns and genes in rows.
#' @param gene.list Gene sets provided either as a list object or as a GeneSetCollection object.
#' @param min.size Minimum size of the gene sets to be considered in the analysis. Default is 3.
#' @param max.size Maximum size of the gene sets to be considered in the analysis. Default is 1000.
#' @param method Method to employ in the estimation of gene-set enrichment scores per sample. Options are "gsva" (default), "ssgsea", "zscore", or "plage".
#' @param ssgsea.normalize Logical vector of length 1; if TRUE runs the ssGSEA method from Barbie et al. (2009) normalizing the scores by the absolute difference between the minimum and the maximum, as described in their paper. Otherwise this last normalization step is skipped.
#' @param ssgsea.alpha Numeric vector of length 1. The exponent defining the weight of the tail in the random walk performed by the ssGSEA (Barbie et al., 2009) method. The default value is 0.25 as described in the paper.
#' @param gsva.kcdf Character vector of length 1 denoting the kernel to use during the non-parametric estimation of the cumulative distribution function of expression levels across samples. By default, kcdf="Gaussian" which is suitable when input expression values are continuous, such as microarray fluorescent units in logarithmic scale, RNA-seq log-CPMs, log-RPKMs or log-TPMs. When input expression values are integer counts, such as those derived from RNA-seq experiments, then this argument should be set to kcdf="Poisson".
#' @param gsva.tau Numeric vector of length 1. The exponent defining the weight of the tail in the random walk performed by the GSVA (Hänzelmann et al., 2013) method. The default value is 1 as described in the paper.
#' @param gsva.maxDiff Logical vector of length 1 which offers two approaches to calculate the enrichment statistic (ES) from the KS random walk statistic. FALSE: ES is calculated as the maximum distance of the random walk from 0. TRUE (the default): ES is calculated as the magnitude difference between the largest positive and negative random walk deviations.
#' @param gsva.absRanking Logical vector of length 1 used only when maxDiff=TRUE. When absRanking=FALSE (default) a modified Kuiper statistic is used to calculate enrichment scores, taking the magnitude difference between the largest positive and negative random walk deviations. When absRanking=TRUE the original Kuiper statistic that sums the largest positive and negative random walk deviations, is used. In this latter case, gene sets with genes enriched on either extreme (high or low) will be regarded as ’highly’ activated.
#' @param verbose Logical indicating whether to print progress messages. Default is TRUE.
#' @param nCores The number of cores to use for parallel computation. Default is `parallel::detectCores() - 2`, which detects the number of cores available on the system and reserves 2 cores for other tasks.
#' @return A gene-set by sample matrix of gene-set enrichment scores.
#' @export
RunGSVA <- function(
    exp,
    gene.list,
    min.size = 3,
    max.size = 1000,
    method = "ssgsea",
    ssgsea.normalize = TRUE,
    ssgsea.alpha = 0.25,
    gsva.kcdf = "Gaussian",
    gsva.tau = 1,
    gsva.maxDiff = TRUE,
    gsva.absRanking = FALSE,
    verbose = TRUE,
    nCores = parallel::detectCores() - 3) {
  # Ensure the method is correctly specified
  if (!method %in% c("ssgsea", "gsva", "zscore", "plage")) {
    stop("Invalid method. Choose from 'ssgsea', 'gsva', 'zscore', or 'plage'.")
  }

  # Check if required packages are installed
  if (!requireNamespace("GSVA", quietly = TRUE)) {
    if (!requireNamespace("BiocManager", quietly = TRUE)) {
      install.packages("BiocManager")
    }
    BiocManager::install("GSVA")
  }

  # Select the appropriate parameters and method for GSVA
  if (method == "ssgsea") {
    params <- GSVA::ssgseaParam(as.matrix(exp), gene.list,
      minSize = min.size, maxSize = max.size,
      normalize = ssgsea.normalize, alpha = ssgsea.alpha
    )
  } else if (method == "gsva") {
    params <- GSVA::gsvaParam(as.matrix(exp), gene.list,
      minSize = min.size, maxSize = max.size,
      kcdf = gsva.kcdf, tau = gsva.tau,
      maxDiff = gsva.maxDiff, absRanking = gsva.absRanking
    )
  } else if (method == "zscore") {
    params <- GSVA::zscoreParam(as.matrix(exp), gene.list,
      minSize = min.size, maxSize = max.size
    )
  } else if (method == "plage") {
    params <- GSVA::plageParam(as.matrix(exp), gene.list,
      minSize = min.size, maxSize = max.size
    )
  }

  # Prepare parallel backend
  BPPARAM <- if (.Platform$OS.type == "windows") {
    BiocParallel::SnowParam(workers = nCores)
  } else {
    BiocParallel::MulticoreParam(workers = nCores)
  }

  # Calculate the gene-set enrichment scores
  score <- as.data.frame(t(GSVA::gsva(params,
    verbose = verbose,
    BPPARAM = BPPARAM
  )))

  return(score)
}
