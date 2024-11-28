#' @title Optimal Feature Combination for Multi-Modality Clustering
#' @description This function selects the optimal combination of features for multi-modality clustering analysis by integrating various modalities (e.g., mutation, CNV, RNA, protein, pathology, radiology) to explore the optimal number of clusters. It uses Nonnegative Matrix Factorization (NMF) for cluster consistency analysis and Multi-block Principal Component Analysis (mbPCA) for assessing cluster separation. The combination of these methods helps identify distinct biological subgroups and evaluate the stability and biological relevance of clustering.
#' @author Zaoqu Liu
#' @param data_layers A named list of matrices representing different modalities. Each matrix should have features as rows and samples as columns. A minimum of two datasets is required.
#' @param feature_subset_sizes A list of sequences representing possible feature subset sizes for each modality. The names of this list must match the names in data_layers.
#' @param try_num_clusters Integer vector. The range of cluster numbers to be tested (default: 2:6).
#' @param n_runs Integer. Number of iterations for NMF for each cluster number (default: 5).
#' @param n_fold Integer. Number of folds for cross-validation in NMF (default: 5).
#' @return A list containing:
#'   - optimal_combination: A dataframe with the optimal feature combination and associated clustering score.
#'   - all_results: A dataframe containing CPI and GAP scores for all tested feature combinations and cluster numbers.
#' @export
Find.OptClusterFeatures <- function(
    data_layers,
    feature_subset_sizes,
    try_num_clusters = 2:6,
    n_runs = 5,
    n_fold = 5) {
  # Load necessary libraries
  if (!requireNamespace("IntNMF", quietly = TRUE)) {
    install.packages("IntNMF")
  }
  if (!requireNamespace("mogsa", quietly = TRUE)) {
    install.packages("mogsa")
  }
  if (!requireNamespace("future.apply", quietly = TRUE)) {
    install.packages("future.apply")
  }

  # Input validation
  if (length(data_layers) < 2) {
    stop("At least two data layers are required for clustering analysis.")
  }
  if (!all(sapply(data_layers, function(x) all(colnames(x) == colnames(data_layers[[1]]))))) {
    stop("All data layers must have matching column names (same samples).")
  }

  # Set up parallel computing
  future::plan("multisession")

  # Generate all possible combinations of feature subset sizes
  feature_combinations <- expand.grid(feature_subset_sizes)

  # Rank features by Median Absolute Deviation (MAD) for each data layer
  ranked_features <- lapply(data_layers, function(layer) {
    names(sort(apply(layer, 1, mad), decreasing = TRUE))
  })

  # Perform clustering analysis for each feature combination using parallel execution
  results <- future.apply::future_lapply(1:nrow(feature_combinations), function(row_n) {
    subset_sizes <- as.numeric(feature_combinations[row_n, ])

    # Subset data based on top-ranked features
    working_data <- mapply(function(layer, feature_ids, size) {
      layer[feature_ids[1:size], ]
    }, data_layers, ranked_features, subset_sizes, SIMPLIFY = FALSE)

    # Normalize data to ensure non-negativity and scale to a maximum value of 1
    normalized_data <- lapply(working_data, function(dataset) {
      dataset <- pmax(dataset + abs(min(dataset)), 0) + .Machine$double.eps
      dataset / max(dataset)
    })

    # Transpose data for downstream analysis
    normalized_data <- lapply(normalized_data, function(x) t(x) + .Machine$double.eps)

    # Optimal cluster number selection using NMF and calculating CPI
    opt_k_nmf <- IntNMF::nmf.opt.k(
      dat = normalized_data,
      n.runs = n_runs,
      n.fold = n_fold,
      k.range = try_num_clusters,
      result = TRUE,
      make.plot = FALSE,
      maxiter = 1000,
      st.count = 10,
      progress = FALSE
    )
    cpi_df <- as.data.frame(opt_k_nmf)
    cpi_df$mean <- rowMeans(cpi_df)

    # Perform multi-block PCA (mbPCA) for integration analysis
    mbpca_result <- mogsa::mbpca(
      x = working_data,
      ncomp = 3,
      k = 0.5,
      method = "globalScore",
      option = "uniform",
      center = TRUE,
      scale = TRUE,
      moa = TRUE,
      svd.solver = "fast",
      maxiter = 1000,
      verbose = FALSE
    )

    # Perform GAP analysis using hierarchical clustering
    gap_analysis <- mogsa::moGap(mbpca_result, K.max = max(try_num_clusters), cluster = "hclust", plot = FALSE)
    gap_df <- as.data.frame(gap_analysis$Tab)[-1, ]

    # Collect results for the current combination
    data.frame(
      Combination_ID = row_n,
      Feature_Subset_Sizes = paste(subset_sizes, collapse = ","),
      K = try_num_clusters,
      CPI = cpi_df$mean,
      GAP = gap_df$gap
    )
  })

  # Combine results from all iterations into one dataframe
  combined_results <- Reduce(rbind, results)
  combined_results$score <- combined_results$CPI + combined_results$GAP

  # Find the optimal feature combination based on the highest combined score
  optimal_combination <- combined_results[which.max(combined_results$score), ]

  return(list(optimal_combination = optimal_combination, all_results = combined_results))
}
