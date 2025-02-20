#' @title Run Intermediate Fusion (IF) for Multi-Modality Data Integration
#' @description This function runs intermediate fusion (IF) analysis using a specified multi-modality clustering algorithm. Users can choose from a variety of clustering algorithms and adjust their respective parameters to perform data integration on multiple modalities, such as RNA, protein, and methylation.
#' @author Zaoqu Liu; Email: liuzaoqu@163.com
#' @param data A list of matrices where each element represents a different modality (e.g., RNA, protein, methylation). Each matrix should have rows as features and columns as samples.
#' @param algorithm Character. The integration algorithm to use. Options include "cpca", "iclusterbayes", "intnmf", "lracluster", "mcia", "nemo", "pinsplus", "rgcca", "sgcca", "snf", "cimlr", "bcc".
#' @param N.clust Integer. Number of clusters to create from the hierarchical clustering of the integrated components (optional but recommended).
#' @param data.types Character vector. Specifies the type of data for each modality (e.g., "binomial", "gaussian"). Default is a mixture of "binomial" and "gaussian".
#' @param BCC.max.iterations Integer. Maximum number of iterations for the BCC algorithm. Default is 10.
#' @param CIMLR.num.dimensions Integer. Number of dimensions for CIMLR. Default is NA.
#' @param CIMLR.tuning.parameter Numeric. Tuning parameter for CIMLR. Default is 10.
#' @param CIMLR.cores.ratio Numeric. Ratio of cores to use for CIMLR. Default is 1.
#' @param CPCA.num.components Integer. Number of components for CPCA. Default is 2.
#' @param CPCA.integration.algorithm Character. Integration algorithm for CPCA. Default is "globalScore".
#' @param CPCA.nonzero.coeff.k Character or numeric. Nonzero coefficient for CPCA. Default is "all".
#' @param CPCA.center.data Logical. Whether to center the data for CPCA. Default is FALSE.
#' @param CPCA.scale.data Logical. Whether to scale the data for CPCA. Default is FALSE.
#' @param CPCA.normalization.option Character. Normalization option for CPCA. Default is "uniform".
#' @param CPCA.max.iterations Integer. Maximum number of iterations for CPCA. Default is 1000.
#' @param CPCA.return.moa.object Logical. Whether to return MOA object for CPCA. Default is TRUE.
#' @param CPCA.show.verbose Logical. Whether to show verbose output for CPCA. Default is FALSE.
#' @param CPCA.svd.solver.method Character. SVD solver method for CPCA. Default is "fast.svd".
#' @param CPCA.nonzero.coeff.obs Character or numeric. Nonzero coefficient for observations in CPCA. Default is "all".
#' @param CPCA.weight.variables Numeric. Weight for variables in CPCA. Default is NA.
#' @param CPCA.weight.observations Numeric. Weight for observations in CPCA. Default is NA.
#' @param CPCA.unit.length.variables Logical. Whether to set unit length for variables in CPCA. Default is FALSE.
#' @param CPCA.unit.length.observations Logical. Whether to set unit length for observations in CPCA. Default is FALSE.
#' @param CPCA.retain.nonnegative Logical. Whether to retain nonnegative values in CPCA. Default is FALSE.
#' @param CPCA.clustering.algorithm Character. Clustering algorithm for CPCA. Default is "ward.D2".
#' @param iClusterBayes.num.burnin Integer. Number of burn-in iterations for iClusterBayes. Default is 20.
#' @param iClusterBayes.num.draws Integer. Number of draws for iClusterBayes. Default is 10.
#' @param iClusterBayes.prior.probabilities Numeric vector. Prior probabilities for iClusterBayes. Default is rep(0.5, length(data)).
#' @param iClusterBayes.proposal.sdev Numeric. Proposal standard deviation for iClusterBayes. Default is 0.05.
#' @param iClusterBayes.thinning.interval Integer. Thinning interval for iClusterBayes. Default is 3.
#' @param IntNMF.max.iterations Integer. Maximum number of iterations for IntNMF. Default is 5.
#' @param IntNMF.stability.count Integer. Stability count for IntNMF. Default is 20.
#' @param IntNMF.num.initializations Integer. Number of initializations for IntNMF. Default is 30.
#' @param IntNMF.use.nndsvd Logical. Whether to use NNDSVD for IntNMF. Default is TRUE.
#' @param IntNMF.random.seed Logical. Whether to use a random seed for IntNMF. Default is TRUE.
#' @param IntNMF.weight Numeric. Weight for IntNMF. Default is NULL.
#' @param LRAcluster.cluster.algorithm Character. Clustering algorithm for LRAcluster. Default is "ward.D2".
#' @param MCIA.n.components Integer. Number of components for MCIA. Default is 10.
#' @param MCIA.clustering.algorithm Character. Clustering algorithm for MCIA. Default is "ward.D2".
#' @param MCIA.scan.eigenvalues Logical. Whether to scan eigenvalues for MCIA. Default is FALSE.
#' @param MCIA.use.nsc Logical. Whether to use NSC for MCIA. Default is TRUE.
#' @param MCIA.use.svd Logical. Whether to use SVD for MCIA. Default is TRUE.
#' @param PINSPlus.agreement.cutoff Numeric. Agreement cutoff for PINSPlus. Default is 0.5.
#' @param PINSPlus.num.cores Integer. Number of cores to use for PINSPlus. Default is 10.
#' @param PINSPlus.sampled.set.size Integer. Sampled set size for PINSPlus. Default is 2000.
#' @param PINSPlus.knn.k Integer. K for k-NN in PINSPlus. Default is NULL.
#' @param RGCCA.connection.matrix Matrix. Connection matrix specifying the relationships between blocks for RGCCA. Default is 1 - diag(length(data)).
#' @param RGCCA.num.components Integer vector. Number of components for each block in RGCCA. Default is rep(1, length(data)).
#' @param RGCCA.scheme Character. Scheme for RGCCA. Default is "centroid".
#' @param RGCCA.regularization Character or numeric vector. Regularization parameter for RGCCA. Default is "optimal".
#' @param RGCCA.scale Logical. Whether to scale data for RGCCA. Default is TRUE.
#' @param RGCCA.initialization Character. Initialization method for RGCCA. Default is "svd".
#' @param RGCCA.bias Logical. Whether to use a biased estimator for RGCCA. Default is TRUE.
#' @param RGCCA.tolerance Numeric. Convergence tolerance for RGCCA. Default is 1e-08.
#' @param RGCCA.verbose Logical. Whether to show progress messages for RGCCA. Default is FALSE.
#' @param RGCCA.clustering.algorithm Character. Clustering algorithm for RGCCA. Default is "ward.D2".
#' @param SGCCA.connection.matrix Matrix. Connection matrix specifying the relationships between blocks for SGCCA. Default is 1 - diag(length(data)).
#' @param SGCCA.num.components.per.modality Integer vector. Number of components per modality for SGCCA. Default is rep(1, length(data)).
#' @param SGCCA.integration.scheme Character. Integration scheme for SGCCA. Default is "centroid".
#' @param SGCCA.sparsity.level Numeric vector. Sparsity level for each block in SGCCA. Default is rep(0.5, length(data)).
#' @param SGCCA.scale.data Logical. Whether to scale data for SGCCA. Default is FALSE.
#' @param SGCCA.initialization.method Character. Initialization method for SGCCA. Default is "svd".
#' @param SGCCA.use.biased.variance Logical. Whether to use a biased estimator for SGCCA. Default is TRUE.
#' @param SGCCA.convergence.tolerance Numeric. Convergence tolerance for SGCCA. Default is .Machine$double.eps.
#' @param SGCCA.show.progress Logical. Whether to show progress messages for SGCCA. Default is FALSE.
#' @param SGCCA.cluster.algorithm Character. Clustering algorithm for SGCCA. Default is "ward.D2".
#' @param SNF.num.neighbors Integer. Number of neighbors for SNF. Default is 20.
#' @param SNF.variance Numeric. Variance for SNF. Default is 0.5.
#' @param SNF.num.iterations Integer. Number of iterations for SNF. Default is 20.
#' @return A list containing the clustering results based on the selected algorithm.
#' @details This function allows the user to integrate multiple data modalities using a variety of different algorithms. Each algorithm has its own parameters that can be adjusted to fit the data and the research question. The function returns clustering results for each sample based on the selected algorithm.
#' @examples
#' # Example usage:
#' data1 <- matrix(rnorm(10000), nrow = 100, ncol = 100)
#' data2 <- matrix(rnorm(10000), nrow = 100, ncol = 100)
#' colnames(data1) <- colnames(data2) <- paste0("Sample", 1:100)
#' data_list <- list(data1, data2)
#'
#' # Run integration clustering using CPCA
#' result <- RunIF(data = data_list, algorithm = "cpca", N.clust = 3)
#' @export
RunIF <- function(data, algorithm, N.clust,
                  data.types = c("binomial", "gaussian", "gaussian", "gaussian", "gaussian", "gaussian"),
                  BCC.max.iterations = 10,
                  CIMLR.num.dimensions = NA,
                  CIMLR.tuning.parameter = 10,
                  CIMLR.cores.ratio = 1,
                  CPCA.num.components = 2,
                  CPCA.integration.algorithm = "globalScore",
                  CPCA.nonzero.coeff.k = "all",
                  CPCA.center.data = FALSE,
                  CPCA.scale.data = FALSE,
                  CPCA.normalization.option = "uniform",
                  CPCA.max.iterations = 1000,
                  CPCA.return.moa.object = TRUE,
                  CPCA.show.verbose = FALSE,
                  CPCA.svd.solver.method = "fast.svd",
                  CPCA.nonzero.coeff.obs = "all",
                  CPCA.weight.variables = NA,
                  CPCA.weight.observations = NA,
                  CPCA.unit.length.variables = FALSE,
                  CPCA.unit.length.observations = FALSE,
                  CPCA.retain.nonnegative = FALSE,
                  CPCA.clustering.algorithm = "ward.D2",
                  iClusterBayes.num.burnin = 20,
                  iClusterBayes.num.draws = 10,
                  iClusterBayes.prior.probabilities = rep(0.5, length(data)),
                  iClusterBayes.proposal.sdev = 0.05,
                  iClusterBayes.thinning.interval = 3,
                  IntNMF.max.iterations = 5,
                  IntNMF.stability.count = 20,
                  IntNMF.num.initializations = 30,
                  IntNMF.use.nndsvd = TRUE,
                  IntNMF.random.seed = TRUE,
                  IntNMF.weight = NULL,
                  LRAcluster.cluster.algorithm = "ward.D2",
                  MCIA.n.components = 10,
                  MCIA.clustering.algorithm = "ward.D2",
                  MCIA.scan.eigenvalues = FALSE,
                  MCIA.use.nsc = TRUE,
                  MCIA.use.svd = TRUE,
                  PINSPlus.agreement.cutoff = 0.5,
                  PINSPlus.num.cores = 10,
                  PINSPlus.sampled.set.size = 2000,
                  PINSPlus.knn.k = NULL,
                  RGCCA.connection.matrix = 1 - diag(length(data)),
                  RGCCA.num.components = rep(1, length(data)),
                  RGCCA.scheme = "centroid",
                  RGCCA.regularization = "optimal",
                  RGCCA.scale = TRUE,
                  RGCCA.initialization = "svd",
                  RGCCA.bias = TRUE,
                  RGCCA.tolerance = 1e-08,
                  RGCCA.verbose = FALSE,
                  RGCCA.clustering.algorithm = "ward.D2",
                  SGCCA.connection.matrix = 1 - diag(length(data)),
                  SGCCA.num.components.per.modality = rep(1, length(data)),
                  SGCCA.integration.scheme = "centroid",
                  SGCCA.sparsity.level = rep(0.5, length(data)),
                  SGCCA.scale.data = FALSE,
                  SGCCA.initialization.method = "svd",
                  SGCCA.use.biased.variance = TRUE,
                  SGCCA.convergence.tolerance = .Machine$double.eps,
                  SGCCA.show.progress = FALSE,
                  SGCCA.cluster.algorithm = "ward.D2",
                  SNF.num.neighbors = 20,
                  SNF.variance = 0.5,
                  SNF.num.iterations = 20,
                  ...) {
  # Set seed for reproducibility
  set.seed(1)

  # Validate input
  if (!is.list(data)) {
    stop("The 'data' parameter must be a list of matrices.")
  }

  # Convert method to lowercase to handle case sensitivity
  algorithm <- tolower(algorithm)

  # Ensure method is one of the allowed options
  valid.algorithms <- c("cpca", "iclusterbayes", "intnmf", "lracluster", "mcia", "nemo", "pinsplus", "rgcca", "sgcca", "snf", "cimlr", "bcc")
  if (!(algorithm %in% valid.algorithms)) {
    stop("Invalid algorithm specified. Please choose from: CPCA, iClusterBayes, IntNMF, LRAcluster, MCIA, NEMO, PINSPlus, RGCCA, SGCCA, SNF, CIMLR, BCC.")
  }

  # Use switch to call the appropriate function based on the selected method
  result <- switch(algorithm,
    "cpca" = RunCPCA(data,
      N.clust = N.clust,
      num.components = CPCA.num.components,
      integration.algorithm = CPCA.integration.algorithm,
      nonzero.coeff.k = CPCA.nonzero.coeff.k,
      center.data = CPCA.center.data,
      scale.data = CPCA.scale.data,
      normalization.option = CPCA.normalization.option,
      max.iterations = CPCA.max.iterations,
      return.moa.object = CPCA.return.moa.object,
      show.verbose = CPCA.show.verbose,
      svd.solver.method = CPCA.svd.solver.method,
      nonzero.coeff.obs = CPCA.nonzero.coeff.obs,
      weight.variables = CPCA.weight.variables,
      weight.observations = CPCA.weight.observations,
      unit.length.variables = CPCA.unit.length.variables,
      unit.length.observations = CPCA.unit.length.observations,
      retain.nonnegative = CPCA.retain.nonnegative,
      clustering.algorithm = CPCA.clustering.algorithm
    ),
    "iclusterbayes" = RuniClusterBayes(data,
      N.clust = N.clust,
      data.type = ifelse(data.types == "binary", "binomial", data.types),
      num.burnin = iClusterBayes.num.burnin,
      num.draws = iClusterBayes.num.draws,
      prior.probabilities = iClusterBayes.prior.probabilities,
      proposal.sdev = iClusterBayes.proposal.sdev,
      thinning.interval = iClusterBayes.thinning.interval
    ),
    "intnmf" = RunIntNMF(data,
      N.clust = N.clust,
      max.iterations = IntNMF.max.iterations,
      stability.count = IntNMF.stability.count,
      num.initializations = IntNMF.num.initializations,
      use.nndsvd = IntNMF.use.nndsvd,
      random.seed = IntNMF.random.seed,
      weight = IntNMF.weight
    ),
    "lracluster" = RunLRAcluster(data,
      N.clust = N.clust,
      data.types = ifelse(data.types == "binomial", "binary", data.types),
      data.names = NULL,
      cluster.algorithm = LRAcluster.cluster.algorithm
    ),
    "mcia" = RunMCIA(data,
      N.clust = N.clust,
      n.components = MCIA.n.components,
      clustering.algorithm = MCIA.clustering.algorithm,
      scan.eigenvalues = MCIA.scan.eigenvalues,
      use.nsc = MCIA.use.nsc,
      use.svd = MCIA.use.svd
    ),
    "nemo" = RunNEMO(data, N.clust = N.clust),
    "pinsplus" = RunPINSPlus(data,
      N.clust = N.clust,
      agreement.cutoff = PINSPlus.agreement.cutoff,
      num.cores = PINSPlus.num.cores,
      sampled.set.size = PINSPlus.sampled.set.size,
      knn.k = PINSPlus.knn.k
    ),
    "rgcca" = RunRGCCA(data,
      N.clust = N.clust,
      connection.matrix = RGCCA.connection.matrix,
      num.components = RGCCA.num.components,
      scheme = RGCCA.scheme,
      regularization = RGCCA.regularization,
      scale = RGCCA.scale,
      initialization = RGCCA.initialization,
      bias = RGCCA.bias,
      tolerance = RGCCA.tolerance,
      verbose = RGCCA.verbose,
      clustering.algorithm = RGCCA.clustering.algorithm
    ),
    "sgcca" = RunSGCCA(data,
      N.clust = N.clust,
      connection.matrix = SGCCA.connection.matrix,
      num.components.per.modality = SGCCA.num.components.per.modality,
      integration.scheme = SGCCA.integration.scheme,
      sparsity.level = SGCCA.sparsity.level,
      scale.data = SGCCA.scale.data,
      initialization.method = SGCCA.initialization.method,
      use.biased.variance = SGCCA.use.biased.variance,
      convergence.tolerance = SGCCA.convergence.tolerance,
      show.progress = SGCCA.show.progress,
      cluster.algorithm = SGCCA.cluster.algorithm
    ),
    "snf" = RunSNF(data,
      N.clust = N.clust,
      num.neighbors = SNF.num.neighbors,
      variance = SNF.num.neighbors,
      num.iterations = SNF.num.neighbors
    ),
    "cimlr" = RunCIMLR(data,
      N.clust = N.clust, num.dimensions = CIMLR.num.dimensions,
      tuning.parameter = CIMLR.tuning.parameter,
      cores.ratio = CIMLR.cores.ratio
    ),
    "bcc" = RunBCC(data, N.clust = N.clust, max.iterations = BCC.max.iterations)
  )

  return(result)
}
