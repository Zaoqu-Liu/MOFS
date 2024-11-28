#' @title Run Classifiers for Cluster Prediction
#' @description This function runs different classification models based on user input to predict cluster assignments for test data.
#' @author Zaoqu Liu; Email: liuzaoqu@163.com
#' @param algorithm A character string indicating the classifier to use. Supported algorithms include: "Adaboost", "DT", "Enet", "Enrichment", "GBDT", "LASSO", "LDA", "NBayes", "NNet", "PCA", "Ridge", "StepLR", "SVD", "SVM", "XGBoost", "kNN", "RF".
#' @param data.test A numeric matrix or data frame of test data. Rows represent genes, and columns represent samples.
#' @param data.train A numeric matrix or data frame of training data. Rows represent genes, and columns represent samples.
#' @param cluster.data A data frame where the first column must be the sample IDs and the second column must be the cluster assignments. The sample IDs must match the column names of the training data.
#' @param cluster.markers A list of data frames, each containing markers for a specific cluster, with columns 'Gene' indicating gene names.
#' @param scale A logical value indicating whether to scale the test data. Default is TRUE.
#' @return A data frame containing the prediction results based on the selected algorithm.
#' @details The function dynamically selects and runs a classification model based on user input. The supported classifiers include a range of machine learning models such as Random Forest, kNN, PCA, SVM, LASSO, Ridge, and others.
#' @examples
#' # Example usage:
#' data.test <- matrix(rnorm(1000), nrow = 100, ncol = 10)
#' data.train <- matrix(rnorm(1000), nrow = 100, ncol = 10)
#' cluster.data <- data.frame(Sample = paste0("Sample", 1:10), Cluster = rep(1:2, each = 5))
#' cluster.markers <- setNames(lapply(unique(cluster.data$Cluster), function(c) data.frame(Gene = paste0("Gene", sample(1:30, 10)), OR = runif(10, 0.5, 2), AUC = runif(10))), unique(cluster.data$Cluster))
#' result <- RunClassifier(algorithm = "RF", data.test, data.train, cluster.data, cluster.markers, scale = TRUE)
#'
#' @export
RunClassifier <- function(algorithm, data.test, data.train, cluster.data, cluster.markers, scale = TRUE) {
  # Convert data.test and data.train to matrices if they are data frames
  data.test <- as.matrix(data.test)
  data.train <- as.matrix(data.train)

  # Input parameter checks
  stopifnot(is.matrix(data.test), is.matrix(data.train), is.data.frame(cluster.data), is.list(cluster.markers))

  # Convert algorithm to lowercase for consistency
  algorithm <- tolower(algorithm)
  available_algorithms <- c("adaboost", "dt", "enet", "enrichment", "gbdt", "lasso", "lda", "nbayes", "nnet", "pca", "ridge", "steplr", "svd", "svm", "xgboost", "knn", "rf")

  # Check if the selected algorithm is available
  if (!algorithm %in% available_algorithms) {
    stop("Invalid algorithm. Available algorithms are: ", paste(available_algorithms, collapse = ", "))
  }

  # Create a switch-case to determine which classifier to run
  classifier_function <- switch(algorithm,
    "adaboost" = Classifier.Adaboost,
    "dt" = Classifier.DT,
    "enet" = Classifier.Enet,
    "enrichment" = Classifier.Enrichment,
    "gbdt" = Classifier.GBDT,
    "lasso" = Classifier.LASSO,
    "lda" = Classifier.LDA,
    "nbayes" = Classifier.NBayes,
    "nnet" = Classifier.NNet,
    "pca" = Classifier.PCA,
    "ridge" = Classifier.Ridge,
    "steplr" = Classifier.StepLR,
    "svd" = Classifier.SVD,
    "svm" = Classifier.SVM,
    "xgboost" = Classifier.XGBoost,
    "knn" = Classifier.kNN,
    "rf" = Classifier.RF
  )

  # Run the selected classifier with the scale parameter
  result <- classifier_function(data.test, data.train, cluster.data, cluster.markers, scale = scale)
  return(result)
}
