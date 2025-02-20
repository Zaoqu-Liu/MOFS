#' @title MOFSR Environment Initialization
#' @description This function sets up the required environment for the MOFSR package by checking and installing the necessary R packages. It also ensures that specific GitHub packages are installed with the correct versions using the `force = TRUE` option.
#'
#' @author Zaoqu Liu; E-mail: liuzaoqu@163.com
#'
#' @details This function is designed to configure the environment for MOFSR by ensuring the installation of a wide range of R packages that are dependencies for the package. In addition to standard CRAN packages, specific packages hosted on GitHub (`genekitr2`, `pathview`, and `GSReg`) are reinstalled to guarantee that the correct versions are used.
#'
#' @return This function does not return any value. It performs the installation of required packages and outputs messages to indicate the status of the installation process.
#'
#' @param force A logical value. If `TRUE`, the function will force reinstallation of the specified GitHub packages (`genekitr2`, `pathview`, and `GSReg`) regardless of whether they are already installed or not. Default is `TRUE`.
#'
#' @examples
#' # Run the init function to set up the MOFSR environment
#' init()
#'
#' @export
init <- function() {
  if (!requireNamespace("BiocManager", quietly = TRUE)) {
    utils::install.packages("BiocManager")
  }
  if (!requireNamespace("remotes", quietly = TRUE)) {
    utils::install.packages("remotes")
  }

  # List of required packages
  required_packages <- c(
    "clue", "glmnet", "adabag", "dplyr", "rpart", "nnet", "tibble", "yaGST",
    "gbm", "caret", "MASS", "e1071", "psych", "randomForest", "mlr3", "mlr3verse",
    "furrr", "purrr", "future", "IntNMF", "dirmult", "ConsensusClusterPlus",
    "survival", "survminer", "progress", "SNFtool", "FactoMineR", "PINSPlus",
    "doMC"
  )
  for (pkg in required_packages) {
    utils::install.packages(pkg)
  }

  # Force reinstall of specific GitHub packages to ensure correct versions
  devtools::install_github("ttriche/bayesCC")
  devtools::install_github("danro9685/CIMLR", ref = "R")
  BiocManager::install("mogsa")
  BiocManager::install("GSVA")
  BiocManager::install("iClusterPlus")
  devtools::install_github("Zaoqu-Liu/LRAcluster")
  BiocManager::install("omicade4")
  devtools::install_github("Shamir-Lab/NEMO/NEMO")
  devtools::install_github("Zaoqu-Liu/RGCCA")
  devtools::install_github("Zaoqu-Liu/ssgsea.GBM.classification")


  # Additional initialization tasks can be added here
  message("Initialization complete!")
}
