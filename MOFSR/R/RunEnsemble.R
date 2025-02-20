#' @title Run Ensemble of Multiple Classifiers for Cluster Prediction
#' @description This function runs an ensemble of different classification models to predict cluster assignments for test data, considering consensus among the models.
#' @author Zaoqu Liu; Email: liuzaoqu@163.com
#' @param data.test A numeric matrix or data frame of test data. Rows represent genes, and columns represent samples.
#' @param data.train A numeric matrix or data frame of training data. Rows represent genes, and columns represent samples.
#' @param cluster.data A data frame where the first column must be the sample IDs and the second column must be the cluster assignments. The sample IDs must match the column names of the training data.
#' @param cluster.markers A list of data frames, each containing markers for a specific cluster, with columns 'Gene' indicating gene names.
#' @param surdata A data frame containing survival information for the samples. The first column must be sample IDs.
#' @param time A character string specifying the column name in `surdata` representing survival time (default: "time").
#' @param event A character string specifying the column name in `surdata` representing the event status (default: "event").
#' @param methods A character vector specifying which classifiers to use in the ensemble. If NULL, all available methods will be used (default: NULL).
#' @param sur.trend.rank A character vector specifying the desired order of survival trends (e.g., c("C2", "C3", "C1")) to filter the models (default: NULL). Note: the order should be from risk to protective.
#' @param cutoff.P Numeric value for the p-value cutoff for survival analysis (default: 0.05).
#' @return A list containing the ensemble prediction results and optional survival analysis.
#' @details This function runs an ensemble of classifiers, checks the consistency among classifiers, and optionally performs survival analysis to filter models based on trends in clinical outcomes.
#' @examples
#' # Example usage:
#' data.test <- matrix(rnorm(1000), nrow = 100, ncol = 10)
#' data.train <- matrix(rnorm(1000), nrow = 100, ncol = 10)
#' cluster.data <- data.frame(Sample = paste0("Sample", 1:10), Cluster = rep(1:2, each = 5))
#' cluster.markers <- setNames(lapply(unique(cluster.data$Cluster), function(c) data.frame(Gene = paste0("Gene", sample(1:30, 10)), OR = runif(10, 0.5, 2), AUC = runif(10))), unique(cluster.data$Cluster))
#' surdata <- data.frame(ID = paste0("Sample", 1:10), time = runif(10, 1, 1000), event = sample(0:1, 10, replace = TRUE))
#' result <- RunEnsemble(data.test, data.train, cluster.data, cluster.markers, surdata, time = "time", event = "event", methods = c("rf", "xgboost"), sur.trend.rank = c("C1", "C3", "C2"))
#'
#' @export
RunEnsemble <- function(
    data.test, data.train,
    cluster.data, cluster.markers,
    surdata = NULL, time = "time", event = "event",
    methods = NULL, sur.trend.rank = NULL, cutoff.P = 0.05) {
  # Install required packages if not already installed
  required_packages <- c("survival", "survminer", "dplyr", "progress")
  for (pkg in required_packages) {
    if (!requireNamespace(pkg, quietly = TRUE)) {
      install.packages(pkg, dependencies = TRUE)
    }
  }

  # Set available methods
  available_methods <- c(
    "adaboost", "dt", "enet", "enrichment", "gbdt", "knn",
    "lasso", "lda", "nbayes", "nnet", "pca", "rf", "ridge",
    "steplr", "svd", "svm", "xgboost"
  )

  # Use user-specified methods or all available methods
  if (!is.null(methods)) {
    methods <- tolower(methods)
    if (!all(methods %in% available_methods)) {
      stop("Invalid methods specified. Available methods are: ", paste(available_methods, collapse = ", "))
    }
  } else {
    methods <- available_methods
  }

  # Create progress bar for ensemble classifiers
  pb <- progress::progress_bar$new(
    format = "Running classifiers [:bar] :current/:total (:percent) in :elapsed",
    total = length(methods),
    clear = FALSE,
    width = 60
  )

  # Run consensus classifiers
  res <- list()
  for (method in methods) {
    pb$tick()
    res[[toupper(method)]] <- tryCatch(
      RunClassifier(method, data.test,
        data.train, cluster.data,
        cluster.markers,
        scale = TRUE
      ),
      error = function(e) NA
    )
  }

  # Filter out failed models
  res <- res[sapply(res, function(x) {
    all(!is.na(x)) &&
      is.data.frame(x) && nrow(x) > 0
  })]

  # Normalize probabilities and prepare results
  res <- lapply(names(res), function(a) {
    b <- res[[a]]
    b[, -c(1, ncol(b))] <- b[, -c(1, ncol(b))] / rowSums(b[, -c(1, ncol(b))], na.rm = TRUE)
    b$Classifiers <- a
    return(b)
  }) %>% Reduce(rbind, .)

  if (!is.null(surdata)) {
    colnames(surdata)[1] <- "ID"
    res <- merge(res, surdata, by = "ID", all.x = TRUE)
  }

  # Perform survival analysis if needed
  if (!is.null(surdata) && !is.null(sur.trend.rank)) {
    message("Performing survival trend analysis...")
    k <- list()
    for (i in unique(res$Classifiers)) {
      d <- res[res$Classifiers == i, ]
      colnames(d)[colnames(d) == time] <- "time"
      colnames(d)[colnames(d) == event] <- "event"
      if (length(unique(d$Predict)) == length(sur.trend.rank)) {
        fit <- survival::survfit(survival::Surv(time, event) ~ Predict, as.data.frame(d))
        sum_fit <- suppressWarnings(as.data.frame(survminer::surv_summary(fit)))
        sum_fit$strata <- as.character(sum_fit$strata)
        sum_fit$strata <- sapply(sum_fit$strata, \(x){
          strsplit(x, "=")[[1]][2]
        })
        sum_fit2 <- dplyr::group_by(sum_fit, strata) %>%
          dplyr::summarise(vv = mean(surv)) %>%
          dplyr::arrange(dplyr::desc(vv))
        cluster_order <- sum_fit2$strata
        p <- survminer::surv_pvalue(fit, d)$pval
        if (all(sur.trend.rank == cluster_order)) {
          k[[i]] <- data.frame(Classifiers = i, P = p, Risk2Protective = paste(cluster_order, collapse = " -> "))
        } else {
          k[[i]] <- NULL
        }
      } else {
        k[[i]] <- NULL
      }
    }
    k <- dplyr::bind_rows(k)

    # Filter classifiers based on survival trend and p-value cutoff
    k2 <- k[k$P < cutoff.P & !is.na(k$P), ]
    if (nrow(k2) == 0) {
      stop("This dataset is too different from our training dataset for validation work!\n>>>>>>>Please use a new dataset...")
    }
    res2 <- res %>%
      dplyr::filter(Classifiers %in% k2$Classifiers) %>%
      dplyr::group_by(ID) %>%
      dplyr::summarise(dplyr::across(dplyr::contains("_prob"), mean, .names = "mean_{.col}"))
    res2$Ensemble_Predict <- apply(res2[, -1], 1, function(x) names(x)[which.max(x)])
    res2$Ensemble_Predict <- gsub("_prob|mean_", "", res2$Ensemble_Predict)
    res2 <- merge(res2, surdata, by = "ID", all.x = TRUE)
    data <- res2
    colnames(data)[colnames(data) == time] <- "time"
    colnames(data)[colnames(data) == event] <- "event"
    fit <- survival::survfit(survival::Surv(time, event) ~ Ensemble_Predict, as.data.frame(data))
    return(list(res = as.data.frame(res2), km = list(fit = fit, data = data), cutoff.P = cutoff.P))
  } else {
    res2 <- res[, c(1, grep("_prob", colnames(res)))] %>%
      dplyr::group_by(ID) %>%
      dplyr::summarise(dplyr::across(dplyr::contains("_prob"), mean, .names = "mean_{.col}"))
    res2$Ensemble_Predict <- apply(res2[, -1], 1, function(x) names(x)[which.max(x)])
    res2$Ensemble_Predict <- gsub("_prob|mean_", "", res2$Ensemble_Predict)
    if (!is.null(surdata)) {
      res2 <- merge(res2, surdata, by = "ID", all.x = TRUE)
    }
    return(list(res = as.data.frame(res2), Readme = "Did not require survival P cutoff value!"))
  }
}
