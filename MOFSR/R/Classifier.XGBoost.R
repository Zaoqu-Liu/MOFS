#' @title XGBoost Classifier for Cluster Prediction
#' @description This function performs classification using XGBoost to predict cluster assignments for test data based on trained models from training data and cluster markers.
#' @author Zaoqu Liu; Email: liuzaoqu@163.com
#' @param data.test A numeric matrix or data frame of test data. Rows represent genes, and columns represent samples.
#' @param data.train A numeric matrix or data frame of training data. Rows represent genes, and columns represent samples.
#' @param cluster.data A data frame where the first column must be the sample IDs and the second column must be the cluster assignments. The sample IDs must match the column names of the training data.
#' @param cluster.markers A list of data frames, each containing markers for a specific cluster, with columns 'Gene' indicating gene names.
#' @param scale A logical value indicating whether to scale the test data. Default is TRUE.
#' @return A data frame with:
#'   - ID: The sample identifier.
#'   - Probabilities: The probabilities for each cluster assignment.
#'   - Predict: The predicted cluster label for each sample.
#' @details The function operates as follows:
#' 1. Ensures that the `cluster.data` has the correct column names.
#' 2. Adds a one-hot encoded matrix for cluster assignments.
#' 3. Scales the test data for prediction.
#' 4. Selects genes that are common between the test and training datasets.
#' 5. Uses glmnet to identify the important markers for each cluster and trains an XGBoost model for classification.
#' 6. Predicts the cluster for test samples and provides probabilities for each cluster.
#' @examples
#' cluster.data <- data.frame(
#'   Sample = paste0("Sample", 1:60),
#'   Cluster = rep(paste0("C", 1:3), each = 20)
#' )
#' data.train <- matrix(rnorm(6000),
#'   nrow = 100,
#'   dimnames = list(paste0("Gene", 1:100), cluster.data$Sample)
#' )
#' data.test <- matrix(rnorm(5000),
#'   nrow = 100,
#'   dimnames = list(paste0("Gene", 1:100), paste0("P", 1:50))
#' )
#' cluster.markers <- setNames(
#'   lapply(
#'     unique(cluster.data$Cluster),
#'     function(cluster) {
#'       data.frame(Gene = sample(rownames(data.train), 10))
#'     }
#'   ),
#'   unique(cluster.data$Cluster)
#' )
#' result <- Classifier.XGBoost(
#'   data.test = data.test, data.train = data.train,
#'   cluster.data, cluster.markers
#' )
#' head(result)
#'
#' @export
Classifier.XGBoost <- function(data.test, data.train, cluster.data, cluster.markers, scale = TRUE) {
  # Install required packages if not already installed
  required_packages <- c("mlr3", "mlr3verse", "glmnet", "dplyr")
  for (pkg in required_packages) {
    if (!requireNamespace(pkg, quietly = TRUE)) {
      install.packages(pkg, dependencies = TRUE)
    }
  }

  # Convert data.test and data.train to matrices if they are data frames
  data.test <- as.matrix(data.test)
  data.train <- as.matrix(data.train)

  # Input parameter checks
  stopifnot(is.matrix(data.test), is.matrix(data.train), is.data.frame(cluster.data), is.list(cluster.markers))

  # Set column names for cluster.data
  cluster.data <- cluster.data[, 1:2]
  colnames(cluster.data) <- c("Sample", "Cluster")
  cluster.data$Cluster <- factor(cluster.data$Cluster, levels = sort(unique(cluster.data$Cluster)))

  # Ensure sample IDs in cluster.data match column names of data.train
  if (!all(cluster.data$Sample %in% colnames(data.train))) {
    stop("Sample IDs in cluster.data must match the column names of data.train!")
  }

  # Add one-hot encoded matrix for cluster assignments
  one_hot_matrix <- stats::model.matrix(~ 0 + Cluster, data = cluster.data)
  cluster.data <- cbind(cluster.data, one_hot_matrix)

  # Scale the test data if scale is TRUE
  if (isTRUE(scale)) {
    data.test <- t(base::scale(t(data.test)))
  }

  # Identify common genes between test and training datasets
  same.genes <- intersect(rownames(data.test), rownames(data.train))
  if (length(same.genes) == 0) {
    stop("No common genes found between test and training datasets!")
  }

  # Subset test and training data to only the intersecting genes
  data.test <- data.test[same.genes, , drop = FALSE]
  data.train <- data.train[same.genes, , drop = FALSE]

  # Identify important markers for each cluster
  cluster.markers2 <- lapply(levels(cluster.data$Cluster), function(cluster) {
    cluster_genes <- cluster.markers[[as.character(cluster)]]$Gene
    # Filter cluster_genes to only those present in data.train
    cluster_genes <- cluster_genes[cluster_genes %in% rownames(data.train)]
    if (length(cluster_genes) == 0) {
      stop(paste("No valid genes found for cluster", cluster))
    }
    data.train_subset <- t(stats::na.omit(data.train[cluster_genes, , drop = FALSE]))
    if (ncol(data.train_subset) == 0) {
      stop(paste("Training data subset is empty for cluster", cluster))
    }
    set.seed(1234)
    cv.fit <- glmnet::cv.glmnet(data.train_subset, as.numeric(cluster.data$Cluster == cluster),
      family = "binomial", alpha = 1, nfolds = 3, type.measure = "auc"
    )
    coefficients <- stats::coef(cv.fit, s = cv.fit$lambda.min)
    active.index <- which(as.numeric(coefficients) != 0)
    if (length(active.index) <= 1) {
      stop(paste("No significant markers found for cluster", cluster))
    }
    data.frame(Gene = coefficients@Dimnames[[1]][active.index][-1], Coef = coefficients[active.index][-1])
  })
  names(cluster.markers2) <- levels(cluster.data$Cluster)

  # Check if there are enough genes for classification
  if (any(sapply(cluster.markers2, function(x) nrow(x) < 1))) {
    stop("Gene intersections between two datasets were insufficient!")
  }

  # Create a list of markers
  markersfor <- unique(unlist(lapply(cluster.markers2, function(x) x$Gene)))
  markersfor <- unique(sapply(markersfor, function(x) strsplit(x, "\\.")[[1]][1]))

  # Prepare training data for XGBoost
  d <- merge(cluster.data[, 1:2], t(data.train[markersfor, , drop = FALSE]), by.x = "Sample", by.y = 0)
  if (ncol(d) <= 2) {
    stop("Training data for XGBoost is insufficient.")
  }
  d <- d[, -1]
  colnames(d) <- gsub("-", "_", colnames(d))
  rownames(data.test) <- gsub("-", "_", rownames(data.test))

  # Prepare training task and train XGBoost model
  task <- mlr3::as_task_classif(d, target = "Cluster", id = "LZQ")
  learner <- mlr3::lrn("classif.xgboost", predict_type = "prob")
  learner$train(task)

  # Prepare test data for prediction
  test_data_for_prediction <- t(data.test[rownames(data.test) %in% colnames(d)[-1], , drop = FALSE])
  if (ncol(test_data_for_prediction) == 0) {
    stop("Test data for prediction is empty.")
  }
  test_data_df <- as.data.frame(test_data_for_prediction)
  test_data_df$Cluster <- factor(rep(levels(cluster.data$Cluster)[1], nrow(test_data_df)), levels = levels(cluster.data$Cluster)) # Adding dummy cluster label for task creation
  test_task <- mlr3::as_task_classif(test_data_df, target = "Cluster", id = "LZQ")

  # Predict cluster assignments for test data
  prediction <- learner$predict(test_task)
  pred_df <- as.data.frame(prediction$prob)
  pred_df$ID <- rownames(test_data_df)
  colnames(pred_df)[-ncol(pred_df)] <- paste0(levels(cluster.data$Cluster), "_prob")
  pred_df$Predict <- prediction$response
  pred_df <- dplyr::select(pred_df, ID, dplyr::everything())

  return(pred_df)
}
