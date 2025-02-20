#' @title Perform ssGSEA-based Subtyping Using Marker Gene Sets
#' @description This function performs single-sample Gene Set Enrichment Analysis (ssGSEA) to assign subtypes to samples based on marker gene sets. It generates necessary input files, runs ssGSEA with permutation testing, and predicts sample subtypes based on the marker enrichment results.
#' @author Zaoqu Liu; Email: liuzaoqu@163.com
#' @param data.test A matrix or data frame representing the input expression data, where rows are genes and columns are samples.
#' @param marker.list A named list of marker gene sets, where each list element corresponds to a specific subtype or category of interest.
#' @param dir.file Character. Directory for saving the output files (default: '.').
#' @param gct.filename Character. The filename for the generated GCT file (default: 'data.gct').
#' @param number.perms Integer. Number of permutations for ssGSEA analysis (default: 10).
#' @param tolerate.mixed Logical. Whether to allow "Mixed" predictions when multiple gene sets have the same minimum p-value (default: FALSE).
#' @return A data frame with the following columns:
#'   - `ID`: Sample identifiers.
#'   - `Predict`: Predicted subtype for each sample.
#'   - Columns with `_pval`: P-values for each marker gene set or subtype.
#' @details The function:
#'   1. Uses the provided marker gene sets to create MOD files for ssGSEA analysis.
#'   2. Generates a GCT file based on the input expression data.
#'   3. Runs ssGSEA with permutation testing to calculate enrichment scores for each marker gene set in every sample.
#'   4. Predicts subtypes for samples by identifying the marker gene set with the most significant enrichment (smallest p-value).
#'   5. If `tolerate.mixed` is TRUE and multiple gene sets share the same minimum p-value, the sample is labeled as "Mixed".
#' @references Wang Q, Hu B, Hu X, Kim H, Squatrito M, Scarpace L, et al. Tumor Evolution of Glioma-Intrinsic Gene Expression Subtypes Associates with Immunological Changes in the Microenvironment. Cancer Cell. July 2017;32(1):42-56.e6.
#' @examples
#' # Simulated expression data
#' data.test <- matrix(rnorm(10000), nrow = 100, ncol = 100)
#' rownames(data.test) <- paste0("Gene", 1:100)
#' colnames(data.test) <- paste0("Sample", 1:100)
#'
#' # Example marker list
#' marker.list <- list(
#'   Subtype1 = c("Gene1", "Gene2", "Gene3"),
#'   Subtype2 = c("Gene4", "Gene5", "Gene6")
#' )
#'
#' # Run ssGSEA-based subtyping
#' result <- Classification.ssGSEA(
#'   data.test = data.test,
#'   marker.list = marker.list,
#'   dir.file = "./results",
#'   gct.filename = "test_data.gct",
#'   number.perms = 10,
#'   tolerate.mixed = TRUE
#' )
#' print(result)
#' @export
Classifier.ssGSEA <- function(data.test, marker.list, dir.file = ".",
                              gct.filename = "data.gct",
                              number.perms = 10,
                              tolerate.mixed = FALSE) {
  # Ensure dir.file ends with '/'
  if (substr(dir.file, nchar(dir.file), nchar(dir.file)) != "/") {
    dir.file <- paste0(dir.file, "/")
  }

  # Create Signature directory if it doesn't exist
  signature_dir <- paste0(dir.file, "Signature/")
  if (!dir.exists(signature_dir)) {
    dir.create(signature_dir, recursive = TRUE)
  }

  # Generate MOD files for each marker gene set
  for (model.name in names(marker.list)) {
    genes <- marker.list[[model.name]]
    generateMODfile(
      model.name = model.name,
      msig.up.genes = genes,
      cl = c("black", "lightgrey"),
      mod_file = paste0(signature_dir, model.name, ".mod")
    )
  }

  # Generate GCT file for input data
  generateGCTfile(
    as.data.frame(data.test),
    gct_file = paste0(dir.file, gct.filename)
  )

  # Run ssGSEA with permutation testing
  runSsGSEAwithPermutation(
    profile_data_file = paste0(dir.file, gct.filename),
    MOD_file = signature_dir,
    select.mods = gsub(".mod", "", list.files(signature_dir, pattern = "\\.mod$")),
    number_perms = number.perms
  )

  unlink(signature_dir, recursive = TRUE)

  # Read p-value results
  pvalue_path <- paste0(dir.file, "p_result_", gct.filename, ".txt")
  if (!file.exists(pvalue_path)) {
    stop("The expected p-value result file does not exist: ", pvalue_path)
  }

  pvalue <- read.table(pvalue_path, row.names = 1, header = TRUE)
  pvalue2 <- pvalue[, grepl("_pval", colnames(pvalue))]

  # Assign predicted subtypes based on minimum p-value
  pvalue2$Predict <- apply(pvalue2, 1, function(x) {
    min_pval <- min(x)
    min_idx <- which(x == min_pval)
    if (tolerate.mixed && length(min_idx) >= 2) {
      return("Mixed")
    } else {
      return(gsub("_pval", "", colnames(pvalue2)[min_idx[1]]))
    }
  })

  # Prepare final output
  pvalue2 <- cbind(ID = rownames(pvalue2), pvalue2)
  rownames(pvalue2) <- NULL

  return(pvalue2)
}
