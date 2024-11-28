#' @title Enrichment-Based Neural Network Classifier for Cluster Prediction
#' @description This function performs classification using pathway enrichment analysis and a neural network to predict cluster assignments for test data based on trained models from training data and cluster markers.
#' @author Zaoqu Liu; Email: liuzaoqu@163.com
#' @param data.test A numeric matrix or data frame of test data. Rows represent genes, and columns represent samples.
#' @param data.train A numeric matrix or data frame of training data. Rows represent genes, and columns represent samples.
#' @param cluster.data A data frame where the first column must be the sample IDs and the second column must be the cluster assignments. The sample IDs must match the column names of the training data.
#' @param cluster.markers A list of data frames, each containing markers for a specific cluster, with columns 'Gene' indicating gene names, 'OR' as odds ratio, and 'AUC' as area under curve.
#' @param scale A logical value indicating whether to scale the test data. Default is TRUE.
#' @param nCores An integer indicating the number of cores to use for pathway enrichment analysis. Default is 5.
#' @return A data frame with:
#'   - ID: The sample identifier.
#'   - Cluster: The predicted cluster label for each sample.
#'   - NES: Normalized Enrichment Score (NES) for each cluster assignment.
#' @details The function operates as follows:
#' 1. Ensures that the `cluster.data` has the correct column names.
#' 2. Adds a one-hot encoded matrix for cluster assignments.
#' 3. Scales the test data for prediction.
#' 4. Selects genes that are common between the test and training datasets.
#' 5. Filters out markers with low odds ratio.
#' 6. Performs pathway enrichment analysis using the ssMwwGST method.
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
#' result <- Classifier.Enrichment(
#'   data.test = data.test, data.train = data.train,
#'   cluster.data, cluster.markers
#' )
#' head(result)
#'
#' @export
Classifier.Enrichment <- function(data.test, data.train, cluster.data, cluster.markers, scale = TRUE, nCores = 5) {
  # Install required packages if not already installed
  required_packages <- c("glmnet", "nnet", "dplyr", "tibble", "yaGST")
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

  # Filter and identify important markers for each cluster
  cluster.markers2 <- lapply(cluster.markers, function(x) {
    x <- x[x$OR > 1, ] %>% dplyr::arrange(desc(AUC))
    l <- ifelse(nrow(x) < 50, nrow(x), 50)
    return(x$Gene[1:l])
  })

  # Check if there are enough genes for classification
  if (any(sapply(cluster.markers2, function(x) length(x) < 3))) {
    stop("Gene intersections between two datasets were insufficient!")
  }

  # Perform pathway enrichment analysis
  res <- ssMwwGST(data.test, cluster.markers2, nCores = nCores)
  NES <- res$NES %>%
    t() %>%
    as.data.frame()
  NES <- as.data.frame(t(apply(NES, 1, exp)))
  NES$Cluster <- apply(NES, 1, which.max)
  NES$Predict <- levels(cluster.data$Cluster)[NES$Cluster]
  NES <- NES %>%
    tibble::rownames_to_column("ID") %>%
    dplyr::select(ID, everything(), Predict, -Cluster)
  colnames(NES)[-c(1, ncol(NES))] <- paste0(colnames(NES)[-c(1, ncol(NES))], "_prob")

  return(NES)
}
