#' @title Run MultiModality Fusion Subtyping (MOFS) for Multi-Modality Data Integration
#' @description This function performs MultiModality Fusion Subtyping (MOFS) analysis by utilizing multiple clustering algorithms for multi-modality data integration. The user can flexibly select the desired clustering algorithms and adjust relevant parameters.
#' @author Zaoqu Liu; Email: liuzaoqu@163.com
#' @param data A list of matrices where each element represents a different modality (e.g., RNA, protein, methylation). Each matrix should have rows as features and columns as samples.
#' @param methods Character vector. The clustering algorithms to use. Options are: "CPCA", "iClusterBayes", "IntNMF", "LRAcluster", "MCIA", "NEMO", "PINSPlus", "RGCCA", "SGCCA", "SNF", "CIMLR", "BCC". At least two methods must be specified.
#' @param max.clusters Integer. The maximum number of clusters to evaluate during consensus clustering analysis (default: 6).
#' @param optimal.clusters Integer. The optimal number of clusters to select from the consensus clustering analysis (default: 3).
#' @param linkage.method Character. The linkage method to use for hierarchical clustering (default: "ward.D2").
#' @param clustering.algorithm Character. The clustering algorithm to use during consensus clustering (default: 'hc').
#' @param distance.metric Character. The distance metric to use for clustering (default: "euclidean").
#' @param resampling.iterations Integer. The number of resampling iterations for consensus clustering (default: 10000).
#' @param resample.proportion Numeric. The proportion of items to resample in each iteration for consensus clustering (default: 0.7).
#' @param silhouette.cutoff Numeric. Silhouette coefficient cutoff value for selecting core set samples (default: 0.4).
#' @param ... Additional parameters specific to the chosen clustering algorithms.
#' @return A list containing the results for each specified clustering algorithm, as well as the results of further analysis including consensus clustering, silhouette scores, core set identification, and PCA.
#' @details The function performs MultiModality Fusion Subtyping (MOFS) by running multiple clustering algorithms on the input multi-modality data. The results of each clustering algorithm are stored for further analysis, including binary cluster assignments, Jaccard distance calculation, consensus clustering analysis, Calinski-Harabasz index calculation, silhouette analysis, and PCA.
#'
#' The steps involved are:
#' 1. Running the specified clustering algorithms on the input data.
#' 2. Extracting binary cluster assignments from the clustering results.
#' 3. Calculating Jaccard distance between clusters.
#' 4. Performing consensus clustering analysis to identify stable clusters.
#' 5. Calculating the Calinski-Harabasz index to assess clustering quality.
#' 6. Performing silhouette analysis to evaluate cluster cohesion and separation.
#' 7. Identifying a core set of samples based on the silhouette coefficient cutoff.
#' 8. Performing PCA on the core set for dimensionality reduction and visualization.
#' @examples
#' # Example usage:
#' data1 <- matrix(rnorm(10000), nrow = 100, ncol = 100)
#' data2 <- matrix(rnorm(10000), nrow = 100, ncol = 100)
#' colnames(data1) <- colnames(data2) <- paste0("Sample", 1:100)
#' data_list <- list(data1, data2)
#'
#' # Run MultiModality Fusion Subtyping using CPCA and CIMLR
#' result <- RunMOFS(data = data_list, methods = c("CPCA", "CIMLR"), max.clusters = 6, optimal.clusters = 3, linkage.method = "ward.D2", clustering.algorithm = "hc", distance.metric = "euclidean", resampling.iterations = 10000, resample.proportion = 0.7, silhouette.cutoff = 0.4)
#'
#' @export
RunMOFS <- function(data, methods, max.clusters = 6, optimal.clusters = 3,
                    linkage.method = "ward.D2", clustering.algorithm = "hc",
                    distance.metric = "euclidean", resampling.iterations = 10000,
                    resample.proportion = 0.7, silhouette.cutoff = 0.4, ...) {
  set.seed(1)

  # Ensure at least two methods are provided
  if (length(methods) < 2) {
    stop("At least two clustering methods must be specified.")
  }

  # Convert methods to lowercase to handle case sensitivity
  methods <- tolower(methods)

  # Validate that the methods are among the allowed options
  valid.methods <- c("cpca", "iclusterbayes", "intnmf", "lracluster", "mcia", "nemo", "pinsplus", "rgcca", "sgcca", "snf", "cimlr", "bcc")
  if (!all(methods %in% valid.methods)) {
    stop("Invalid methods specified. Please choose from: CPCA, iClusterBayes, IntNMF, LRAcluster, MCIA, NEMO, PINSPlus, RGCCA, SGCCA, SNF, CIMLR, BCC.")
  }

  # Use lapply to run each clustering function via RunIF and store the results in a list
  res <- lapply(methods, function(m) {
    RunIF(data = data, method = m, N.clust = optimal.clusters, ...)
  })
  names(res) <- methods

  # Get binary cluster assignments for each method
  resm <- get.binary.clusters(res)

  # Calculate Jaccard distance between binary cluster assignments
  jaccard.distance <- get.Jaccard.Distance(resm, dissimilarity = TRUE)

  # Run Consensus Clustering Analysis (COCA)
  coca.results <- RunCOCA(
    jaccard.matrix = jaccard.distance, max.clusters = max.clusters, optimal.clusters = optimal.clusters,
    linkage.method = linkage.method, clustering.algorithm = clustering.algorithm, distance.metric = distance.metric,
    resampling.iterations = resampling.iterations, resample.proportion = resample.proportion
  )

  # Calculate PAC values
  PAC <- CalPAC(coca.results$fit)

  # Perform hierarchical clustering and calculate Calinski-Harabasz index
  hierarchical.clustering <- stats::hclust(stats::as.dist(jaccard.distance),
    method = "average"
  )
  calinski.harabasz.index <- CalCHI(hierarchical.clustering, max_clusters = max.clusters)

  # Silhouette analysis on the consensus matrix
  consensus.matrix <- 1 - coca.results$optimal$consensusMatrix
  rownames(consensus.matrix) <- colnames(consensus.matrix) <- colnames(resm)
  consensus.matrix <- as.dist(consensus.matrix)
  cluster.assignments <- coca.results$Cluster$Cluster
  names(cluster.assignments) <- colnames(resm)
  silhouette.analysis <- cluster::silhouette(x = cluster.assignments, dist = consensus.matrix)

  # Plot silhouette analysis
  plot(silhouette.analysis, col = c("red", "blue", "yellow"), main = "")

  # Select core set based on silhouette coefficient cutoff
  core.samples <- colnames(resm)[silhouette.analysis[, 3] >= silhouette.cutoff]

  # Perform PCA on the core set
  jaccard.core <- jaccard.distance[core.samples, core.samples]
  pca.results <- RunPCA(jaccard.core)

  # Get final cluster assignments for the core set
  core.cluster <- get.class(coca.results$fit, optimal.clusters)
  core.cluster <- core.cluster[match(core.samples, core.cluster$ID), ]

  # Visualize PCA with clusters
  print(factoextra::fviz_pca_ind(
    X = pca.results,
    geom.ind = "point",
    pointshape = 21,
    fill.ind = core.cluster$Cluster,
    palette = "npg",
    alpha.ind = 0.7,
    addEllipses = TRUE
  ) +
    ggplot2::theme_classic() +
    ggplot2::theme(plot.title = ggplot2::element_text(hjust = 0.5)))

  # Collect all results in a list
  all.results <- list(
    All_Cluster = resm,
    JaccardDistance = jaccard.distance,
    ConsensusClustering = coca.results,
    Silhouette = silhouette.analysis,
    CoreSet = core.samples,
    PCA = pca.results,
    PAC = PAC,
    CalinskiHarabasz = calinski.harabasz.index
  )

  return(all.results)
}
