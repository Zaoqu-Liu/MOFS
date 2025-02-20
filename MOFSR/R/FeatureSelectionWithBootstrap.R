#' @title Feature Selection with Bootstrap for Each Cluster
#' @description This function performs feature selection using a logistic regression model and bootstrapping for each cluster in the dataset.
#' For each cluster, the function evaluates which features are significantly associated with the cluster based on a logistic regression model with bootstrap sampling.
#' The output includes a count of significant features for each cluster based on the p-value threshold.
#'
#' @param data A data frame where the first column is the cluster variable (categorical), and the other columns are feature values.
#' @param p.no_bootstrap Numeric. The p-value threshold for the features to be selected based on the original data (default: 0.01).
#' @param p.bootstrap Numeric. The p-value threshold for the features to be selected based on bootstrapped samples (default: 0.05).
#' @param num.iteration Integer. The number of bootstrap iterations (default: 1000).
#' @param nCores Integer. The number of CPU cores to use for parallel processing (default: automatically set to all cores minus 3).
#'
#' @return A list of data frames, each containing the features selected for each cluster and the count of significant features based on bootstrapping.
#' @author Zaoqu Liu; Email: liuzaoqu@163.com
#' @export
FeatureSelectionWithBootstrap <- function(data, p.no_bootstrap = 0.01, p.bootstrap = 0.05, num.iteration = 1000, nCores = parallel::detectCores() - 3) {
  # Prepare the dataset
  colnames(data)[1] <- "Cluster" # Ensure the first column is named 'Cluster'

  tmp <- list() # Initialize list to store results

  # Check if the Cluster column has more than one unique value
  unique_clusters <- unique(data$Cluster)
  if (length(unique_clusters) == 2) {
    # If there are only two clusters, directly proceed without looping
    dd <- data
    dd$Cluster <- ifelse(dd$Cluster == unique_clusters[1], 1, 0) # Create binary cluster variable (1 for first cluster, 0 for others)

    # Initial logistic regression to identify significant features
    rg <- apply(dd[, -1], 2, function(x) {
      fit <- summary(stats::glm(dd$Cluster ~ x, family = "binomial"))
      return(fit$coefficients[2, 4]) # p-value for the feature
    })

    # Filter columns based on p-value threshold from the initial logistic regression
    dd <- dd[, c("Cluster", names(rg)[rg < p.no_bootstrap])]

    # Parallel processing setup
    future::plan(future::multisession, workers = nCores) # Set up parallel processing

    # Perform bootstrapping
    res <- furrr::future_map(colnames(dd[, 2:ncol(dd)]), function(x) {
      Mboot <- purrr::map_vec(1:num.iteration, function(y) {
        indices <- caret::createDataPartition(dd$Cluster, p = 0.7, list = FALSE) %>% as.numeric() # Create bootstrap samples
        data <- dd[indices, ]
        fit <- summary(stats::glm(data$Cluster ~ data[, x], family = "binomial", data = data))
        P <- fit$coefficients[2, 4] # p-value for the feature
      })
      return(Mboot) # Return the p-values for the bootstrapped samples
    }, .progress = TRUE)

    # Combine results from all features
    res2 <- Reduce(cbind, res) %>% as.data.frame()

    # Store the count of significant features based on bootstrap results
    tmp[[unique_clusters[1]]] <- data.frame(ID = colnames(dd)[2:ncol(dd)], Eli = apply(res2, 2, function(x) sum(x < p.bootstrap)))
  } else {
    # Loop over each unique cluster if there are more than two
    for (i in sort(unique_clusters)) {
      dd <- data
      dd$Cluster <- ifelse(dd$Cluster == i, 1, 0) # Create binary cluster variable (1 for current cluster, 0 for others)

      # Initial logistic regression to identify significant features
      rg <- apply(dd[, -1], 2, function(x) {
        fit <- summary(stats::glm(dd$Cluster ~ x, family = "binomial"))
        return(fit$coefficients[2, 4]) # p-value for the feature
      })

      # Filter columns based on p-value threshold from the initial logistic regression
      dd <- dd[, c("Cluster", names(rg)[rg < p.no_bootstrap])]

      # Parallel processing setup
      future::plan(works = nCores) # Set up parallel processing

      # Perform bootstrapping
      res <- furrr::future_map(colnames(dd[, 2:ncol(dd)]), function(x) {
        Mboot <- purrr::map_vec(1:num.iteration, function(y) {
          indices <- caret::createDataPartition(dd$Cluster, p = 0.7, list = FALSE) %>% as.numeric() # Create bootstrap samples
          data <- dd[indices, ]
          fit <- summary(stats::glm(data$Cluster ~ data[, x], family = "binomial", data = data))
          P <- fit$coefficients[2, 4] # p-value for the feature
        })
        return(Mboot) # Return the p-values for the bootstrapped samples
      }, .progress = TRUE)

      # Combine results from all features
      res2 <- Reduce(cbind, res) %>% as.data.frame()

      # Store the count of significant features based on bootstrap results
      tmp[[i]] <- data.frame(ID = colnames(dd)[2:ncol(dd)], Eli = apply(res2, 2, function(x) sum(x < p.bootstrap)))
    }
  }

  return(tmp) # Return the list of results for each cluster
}
