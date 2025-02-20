generateMODfile <- function(model.name, msig.up.genes, cl, mod_file) {
  cat(paste0("model.creation.date", "\t", "'", Sys.time(), "'\n"), file = mod_file)
  cat(paste0("model.name", "\t", "'", model.name, "'\n"), file = mod_file, append = T)
  cat(paste0("model.description", "\t", "'", "Paper In Submit_TBD", "'\n"), file = mod_file, append = T)
  cat(paste0("sample.norm.type", "\t", "'", "rank", "'\n"), file = mod_file, append = T)
  cat(paste0("statistic", "\t", "'", "area.under.RES", "'\n"), file = mod_file, append = T)
  cat(paste0("weight", "\t", 0.25, "\n"), file = mod_file, append = T)
  cat(paste0("random.seed", "\t", 12345, "\n"), file = mod_file, append = T)
  cat(paste0("nperm", "\t", 0, "\n"), file = mod_file, append = T)
  cat(paste0("link.function", "\t", "'", "logit", "'\n"), file = mod_file, append = T)
  cat(paste0("cl", "\t", "c(", paste0(sQuote(cl, q = F), collapse = ","), ")\n"), file = mod_file, append = T)
  cat(paste0("msig.up.genes", "\t", "c(", paste0(sQuote(msig.up.genes, q = F), collapse = ","), ")"), file = mod_file, append = T)
}

generateGCTfile <- function(in_gct, gct_file) {
  in_gct <- data.frame(
    GeneID = rownames(in_gct),
    description = "na",
    in_gct,
    stringsAsFactors = F,
    check.names = F
  )
  cat("#1.2\n", file = gct_file)
  cat(nrow(in_gct), "\t", ncol(in_gct) - 2, "\n", file = gct_file, append = T)
  cat(paste(colnames(in_gct), collapse = "\t"), "\n", file = gct_file, append = T)
  for (i in 1:nrow(in_gct)) cat(paste(in_gct[i, ], collapse = "\t"), "\n", file = gct_file, append = T)
}


# The Broad Institute
# SOFTWARE COPYRIGHT NOTICE AGREEMENT
# This software and its documentation are copyright 2006 by the
# Broad Institute/Massachusetts Institute of Technology.
# All rights are reserved.
#
# This software is supplied without any warranty or guaranteed support
# whatsoever. Neither the Broad Institute nor MIT can be responsible for
# its use, misuse, or functionality.


# Molecular signatures tools

# Auxiliary functions and definitions

GSEA.GeneRanking <- function(A, class.labels, gene.labels, nperm, permutation.type = 0, sigma.correction = "GeneCluster", fraction = 1.0, replace = F, reverse.sign = F) {
  # This function ranks the genes according to the signal to noise ratio for the actual phenotype and also random permutations and bootstrap
  # subsamples of both the observed and random phenotypes. It uses matrix operations to implement the signal to noise calculation
  # in stages and achieves fast execution speed. It supports two types of permutations: random (unbalanced) and balanced.
  # It also supports subsampling and bootstrap by using masking and multiple-count variables.  When "fraction" is set to 1 (default)
  # the there is no subsampling or boostrapping and the matrix of observed signal to noise ratios will have the same value for
  # all permutations. This is wasteful but allows to support all the multiple options with the same code. Notice that the second
  # matrix for the null distribution will still have the values for the random permutations g
  # (null distribution). This mode (fraction = 1.0) is the defaults, the recommended one and the one used in the examples.
  # It is also the one that has be tested more thoroughly. The resampling and boostrapping options are intersting to obtain
  # smooth estimates of the observed distribution but its is left for the expert user who may want to perform some sanity
  # checks before trusting the code.
  #
  # Inputs:
  #   A: Matrix of gene expression values (rows are genes, columns are samples)
  #   class.labels: Phenotype of class disticntion of interest. A vector of binary labels having first the 1's and then the 0's
  #   gene.labels: gene labels. Vector of probe ids or accession numbers for the rows of the expression matrix
  #   nperm: Number of random permutations/bootstraps to perform
  #   permutation.type: Permutation type: 0 = unbalanced, 1 = balanced. For experts only (default: 0)
  #   sigma.correction: Correction to the signal to noise ratio (Default = GeneCluster, a choice to support the way it was handled in a previous package)
  #   fraction: Subsampling fraction. Set to 1.0 (no resampling). For experts only (default: 1.0)
  #   replace: Resampling mode (replacement or not replacement). For experts only (default: F)
  #   reverse.sign: Reverse direction of gene list (default = F)
  #
  # Outputs:
  #   s2n.matrix: Matrix with random permuted or bootstraps singal to noise rations (rows are genes, columns are permutations or bootstrap subsamplings
  #   obs.s2n.matrix: Matrix with observed signal to noise ratios (rows are genes, columns are boostraps subsamplings. If fraction is set to 1.0 then all the columns have the same values
  #   order.matrix: Matrix with the orderings that will sort the columns of the obs.s2n.matrix in decreasing s2n order
  #   obs.order.matrix: Matrix with the orderings that will sort the columns of the s2n.matrix in decreasing s2n order
  #
  # The Broad Institute
  # SOFTWARE COPYRIGHT NOTICE AGREEMENT
  # This software and its documentation are copyright 2003 by the
  # Broad Institute/Massachusetts Institute of Technology.
  # All rights are reserved.
  #
  # This software is supplied without any warranty or guaranteed support
  # whatsoever. Neither the Broad Institute nor MIT can be responsible for
  # its use, misuse, or functionality.

  A <- A + 0.00000001

  N <- length(A[, 1])
  Ns <- length(A[1, ])

  subset.mask <- matrix(0, nrow = Ns, ncol = nperm)
  reshuffled.class.labels1 <- matrix(0, nrow = Ns, ncol = nperm)
  reshuffled.class.labels2 <- matrix(0, nrow = Ns, ncol = nperm)
  class.labels1 <- matrix(0, nrow = Ns, ncol = nperm)
  class.labels2 <- matrix(0, nrow = Ns, ncol = nperm)

  order.matrix <- matrix(0, nrow = N, ncol = nperm)
  obs.order.matrix <- matrix(0, nrow = N, ncol = nperm)
  s2n.matrix <- matrix(0, nrow = N, ncol = nperm)
  obs.s2n.matrix <- matrix(0, nrow = N, ncol = nperm)

  obs.gene.labels <- vector(length = N, mode = "character")
  obs.gene.descs <- vector(length = N, mode = "character")
  obs.gene.symbols <- vector(length = N, mode = "character")

  M1 <- matrix(0, nrow = N, ncol = nperm)
  M2 <- matrix(0, nrow = N, ncol = nperm)
  S1 <- matrix(0, nrow = N, ncol = nperm)
  S2 <- matrix(0, nrow = N, ncol = nperm)

  gc()

  C <- split(class.labels, class.labels)
  class1.size <- length(C[[1]])
  class2.size <- length(C[[2]])
  class1.index <- seq(1, class1.size, 1)
  class2.index <- seq(class1.size + 1, class1.size + class2.size, 1)

  for (r in 1:nperm) {
    class1.subset <- sample(class1.index, size = ceiling(class1.size * fraction), replace = replace)
    class2.subset <- sample(class2.index, size = ceiling(class2.size * fraction), replace = replace)
    class1.subset.size <- length(class1.subset)
    class2.subset.size <- length(class2.subset)
    subset.class1 <- rep(0, class1.size)
    for (i in 1:class1.size) {
      if (is.element(class1.index[i], class1.subset)) {
        subset.class1[i] <- 1
      }
    }
    subset.class2 <- rep(0, class2.size)
    for (i in 1:class2.size) {
      if (is.element(class2.index[i], class2.subset)) {
        subset.class2[i] <- 1
      }
    }
    subset.mask[, r] <- as.numeric(c(subset.class1, subset.class2))
    fraction.class1 <- class1.size / Ns
    fraction.class2 <- class2.size / Ns

    if (permutation.type == 0) { # random (unbalanced) permutation
      full.subset <- c(class1.subset, class2.subset)
      label1.subset <- sample(full.subset, size = Ns * fraction.class1)
      reshuffled.class.labels1[, r] <- rep(0, Ns)
      reshuffled.class.labels2[, r] <- rep(0, Ns)
      class.labels1[, r] <- rep(0, Ns)
      class.labels2[, r] <- rep(0, Ns)
      for (i in 1:Ns) {
        m1 <- sum(!is.na(match(label1.subset, i)))
        m2 <- sum(!is.na(match(full.subset, i)))
        reshuffled.class.labels1[i, r] <- m1
        reshuffled.class.labels2[i, r] <- m2 - m1
        if (i <= class1.size) {
          class.labels1[i, r] <- m2
          class.labels2[i, r] <- 0
        } else {
          class.labels1[i, r] <- 0
          class.labels2[i, r] <- m2
        }
      }
    } else if (permutation.type == 1) { # proportional (balanced) permutation

      class1.label1.subset <- sample(class1.subset, size = ceiling(class1.subset.size * fraction.class1))
      class2.label1.subset <- sample(class2.subset, size = floor(class2.subset.size * fraction.class1))
      reshuffled.class.labels1[, r] <- rep(0, Ns)
      reshuffled.class.labels2[, r] <- rep(0, Ns)
      class.labels1[, r] <- rep(0, Ns)
      class.labels2[, r] <- rep(0, Ns)
      for (i in 1:Ns) {
        if (i <= class1.size) {
          m1 <- sum(!is.na(match(class1.label1.subset, i)))
          m2 <- sum(!is.na(match(class1.subset, i)))
          reshuffled.class.labels1[i, r] <- m1
          reshuffled.class.labels2[i, r] <- m2 - m1
          class.labels1[i, r] <- m2
          class.labels2[i, r] <- 0
        } else {
          m1 <- sum(!is.na(match(class2.label1.subset, i)))
          m2 <- sum(!is.na(match(class2.subset, i)))
          reshuffled.class.labels1[i, r] <- m1
          reshuffled.class.labels2[i, r] <- m2 - m1
          class.labels1[i, r] <- 0
          class.labels2[i, r] <- m2
        }
      }
    }
  }

  # compute S2N for the random permutation matrix

  P <- reshuffled.class.labels1 * subset.mask
  n1 <- sum(P[, 1])
  M1 <- A %*% P
  M1 <- M1 / n1
  gc()
  A2 <- A * A
  S1 <- A2 %*% P
  S1 <- S1 / n1 - M1 * M1
  S1 <- sqrt(abs((n1 / (n1 - 1)) * S1))
  gc()
  P <- reshuffled.class.labels2 * subset.mask
  n2 <- sum(P[, 1])
  M2 <- A %*% P
  M2 <- M2 / n2
  gc()
  A2 <- A * A
  S2 <- A2 %*% P
  S2 <- S2 / n2 - M2 * M2
  S2 <- sqrt(abs((n2 / (n2 - 1)) * S2))
  rm(P)
  rm(A2)
  gc()

  if (sigma.correction == "GeneCluster") { # small sigma "fix" as used in GeneCluster
    S2 <- ifelse(0.2 * abs(M2) < S2, S2, 0.2 * abs(M2))
    S2 <- ifelse(S2 == 0, 0.2, S2)
    S1 <- ifelse(0.2 * abs(M1) < S1, S1, 0.2 * abs(M1))
    S1 <- ifelse(S1 == 0, 0.2, S1)
    gc()
  }

  M1 <- M1 - M2
  rm(M2)
  gc()
  S1 <- S1 + S2
  rm(S2)
  gc()

  s2n.matrix <- M1 / S1

  if (reverse.sign == T) {
    s2n.matrix <- -s2n.matrix
  }
  gc()

  # compute S2N for the "observed" permutation matrix

  P <- class.labels1 * subset.mask
  n1 <- sum(P[, 1])
  M1 <- A %*% P
  M1 <- M1 / n1
  gc()
  A2 <- A * A
  S1 <- A2 %*% P
  S1 <- S1 / n1 - M1 * M1
  S1 <- sqrt(abs((n1 / (n1 - 1)) * S1))
  gc()
  P <- class.labels2 * subset.mask
  n2 <- sum(P[, 1])
  M2 <- A %*% P
  M2 <- M2 / n2
  gc()
  A2 <- A * A
  S2 <- A2 %*% P
  S2 <- S2 / n2 - M2 * M2
  S2 <- sqrt(abs((n2 / (n2 - 1)) * S2))
  rm(P)
  rm(A2)
  gc()

  if (sigma.correction == "GeneCluster") { # small sigma "fix" as used in GeneCluster
    S2 <- ifelse(0.2 * abs(M2) < S2, S2, 0.2 * abs(M2))
    S2 <- ifelse(S2 == 0, 0.2, S2)
    S1 <- ifelse(0.2 * abs(M1) < S1, S1, 0.2 * abs(M1))
    S1 <- ifelse(S1 == 0, 0.2, S1)
    gc()
  }

  M1 <- M1 - M2
  rm(M2)
  gc()
  S1 <- S1 + S2
  rm(S2)
  gc()

  obs.s2n.matrix <- M1 / S1
  gc()

  if (reverse.sign == T) {
    obs.s2n.matrix <- -obs.s2n.matrix
  }

  for (r in 1:nperm) {
    obs.order.matrix[, r] <- order(obs.s2n.matrix[, r], decreasing = T)
  }

  return(list(
    s2n.matrix = s2n.matrix,
    obs.s2n.matrix = obs.s2n.matrix,
    order.matrix = order.matrix,
    obs.order.matrix = obs.order.matrix
  ))
}

Wilcox.Score <- function(gene.list, gene.set) {
  library(exactRankTests)
  N <- length(gene.list)
  tag.indicator <- sign(match(gene.list, gene.set, nomatch = 0)) # notice that the sign is 0 (no tag) or 1 (tag)
  seq.index <- seq(1, N)
  gene.set.ranks <- seq.index[tag.indicator == 1]
  gene.set.comp.ranks <- seq.index[tag.indicator == 0]

  W <- wilcox.exact(x = gene.set.ranks, y = gene.set.comp.ranks, alternative = "two.sided", mu = 0, paired = FALSE, exact = F, conf.int = T, conf.level = 0.95)

  #      ES <- 1/W$p.value
  #      ES <- 1 - W$p.value
  #      ES <- W$statistic/N
  ES <- log(1 / W$p.value)

  return(list(ES = ES))
}





GSEA.EnrichmentScore <- function(gene.list, gene.set, weighted.score.type = 1, correl.vector = NULL) {
  #
  # Computes the weighted GSEA score of gene.set in gene.list.
  # The weighted score type is the exponent of the correlation
  # weight: 0 (unweighted = Kolmogorov-Smirnov), 1 (weighted), and 2 (over-weighted). When the score type is 1 or 2 it is
  # necessary to input the correlation vector with the values in the same order as in the gene list.
  #
  # Inputs:
  #   gene.list: The ordered gene list (e.g. integers indicating the original position in the input dataset)
  #   gene.set: A gene set (e.g. integers indicating the location of those genes in the input dataset)
  #   weighted.score.type: Type of score: weight: 0 (unweighted = Kolmogorov-Smirnov), 1 (weighted), and 2 (over-weighted)
  #  correl.vector: A vector with the coorelations (e.g. signal to noise scores) corresponding to the genes in the gene list
  #
  # Outputs:
  #   ES: Enrichment score (real number between -1 and +1)
  #   arg.ES: Location in gene.list where the peak running enrichment occurs (peak of the "mountain")
  #   RES: Numerical vector containing the running enrichment score for all locations in the gene list
  #   tag.indicator: Binary vector indicating the location of the gene sets (1's) in the gene list
  #
  # The Broad Institute
  # SOFTWARE COPYRIGHT NOTICE AGREEMENT
  # This software and its documentation are copyright 2003 by the
  # Broad Institute/Massachusetts Institute of Technology.
  # All rights are reserved.
  #
  # This software is supplied without any warranty or guaranteed support
  # whatsoever. Neither the Broad Institute nor MIT can be responsible for
  # its use, misuse, or functionality.

  tag.indicator <- sign(match(gene.list, gene.set, nomatch = 0)) # notice that the sign is 0 (no tag) or 1 (tag)
  no.tag.indicator <- 1 - tag.indicator
  N <- length(gene.list)
  Nh <- length(gene.set)
  Nm <- N - Nh
  if (weighted.score.type == 0) {
    correl.vector <- rep(1, N)
  }
  alpha <- weighted.score.type
  correl.vector <- abs(correl.vector**alpha)
  sum.correl.tag <- sum(correl.vector[tag.indicator == 1])
  norm.tag <- 1.0 / sum.correl.tag
  norm.no.tag <- 1.0 / Nm
  RES <- cumsum(tag.indicator * correl.vector * norm.tag - no.tag.indicator * norm.no.tag)
  max.ES <- max(RES)
  min.ES <- min(RES)
  if (max.ES > -min.ES) {
    #      ES <- max.ES
    ES <- signif(max.ES, digits = 5)
    arg.ES <- which.max(RES)
  } else {
    #      ES <- min.ES
    ES <- signif(min.ES, digits = 5)
    arg.ES <- which.min(RES)
  }
  return(list(ES = ES, arg.ES = arg.ES, RES = RES, indicator = tag.indicator))
}

GSEA.EnrichmentScore2 <- function(gene.list, gene.set, weighted.score.type = 1, correl.vector = NULL) {
  #
  # Computes the weighted GSEA score of gene.set in gene.list. It is the same calculation as in
  # GSEA.EnrichmentScore but faster (x8) without producing the RES, arg.RES and tag.indicator outputs.
  # This call is intended to be used to asses the enrichment of random permutations rather than the
  # observed one.
  # The weighted score type is the exponent of the correlation
  # weight: 0 (unweighted = Kolmogorov-Smirnov), 1 (weighted), and 2 (over-weighted). When the score type is 1 or 2 it is
  # necessary to input the correlation vector with the values in the same order as in the gene list.
  #
  # Inputs:
  #   gene.list: The ordered gene list (e.g. integers indicating the original position in the input dataset)
  #   gene.set: A gene set (e.g. integers indicating the location of those genes in the input dataset)
  #   weighted.score.type: Type of score: weight: 0 (unweighted = Kolmogorov-Smirnov), 1 (weighted), and 2 (over-weighted)
  #  correl.vector: A vector with the coorelations (e.g. signal to noise scores) corresponding to the genes in the gene list
  #
  # Outputs:
  #   ES: Enrichment score (real number between -1 and +1)
  #
  # The Broad Institute
  # SOFTWARE COPYRIGHT NOTICE AGREEMENT
  # This software and its documentation are copyright 2003 by the
  # Broad Institute/Massachusetts Institute of Technology.
  # All rights are reserved.
  #
  # This software is supplied without any warranty or guaranteed support
  # whatsoever. Neither the Broad Institute nor MIT can be responsible for
  # its use, misuse, or functionality.

  N <- length(gene.list)
  Nh <- length(gene.set)
  Nm <- N - Nh

  loc.vector <- vector(length = N, mode = "numeric")
  peak.res.vector <- vector(length = Nh, mode = "numeric")
  valley.res.vector <- vector(length = Nh, mode = "numeric")
  tag.correl.vector <- vector(length = Nh, mode = "numeric")
  tag.diff.vector <- vector(length = Nh, mode = "numeric")
  tag.loc.vector <- vector(length = Nh, mode = "numeric")

  loc.vector[gene.list] <- seq(1, N)
  tag.loc.vector <- loc.vector[gene.set]

  tag.loc.vector <- sort(tag.loc.vector, decreasing = F)

  if (weighted.score.type == 0) {
    tag.correl.vector <- rep(1, Nh)
  } else if (weighted.score.type == 1) {
    tag.correl.vector <- correl.vector[tag.loc.vector]
    tag.correl.vector <- abs(tag.correl.vector)
  } else if (weighted.score.type == 2) {
    tag.correl.vector <- correl.vector[tag.loc.vector] * correl.vector[tag.loc.vector]
    tag.correl.vector <- abs(tag.correl.vector)
  } else {
    tag.correl.vector <- correl.vector[tag.loc.vector]**weighted.score.type
    tag.correl.vector <- abs(tag.correl.vector)
  }

  norm.tag <- 1.0 / sum(tag.correl.vector)
  tag.correl.vector <- tag.correl.vector * norm.tag
  norm.no.tag <- 1.0 / Nm
  tag.diff.vector[1] <- (tag.loc.vector[1] - 1)
  tag.diff.vector[2:Nh] <- tag.loc.vector[2:Nh] - tag.loc.vector[1:(Nh - 1)] - 1
  tag.diff.vector <- tag.diff.vector * norm.no.tag
  peak.res.vector <- cumsum(tag.correl.vector - tag.diff.vector)
  valley.res.vector <- peak.res.vector - tag.correl.vector
  max.ES <- max(peak.res.vector)
  min.ES <- min(valley.res.vector)
  ES <- ifelse(max.ES > -min.ES, max.ES, min.ES)

  return(list(ES = ES))
}

GSEA.EnrichmentScore3 <- function(gene.list, gene.set, weighted.score.type = 1, correl.vector = NULL) {
  #
  # Computes the weighted GSEA score of gene.set in gene.list.
  # The weighted score type is the exponent of the correlation
  # weight: 0 (unweighted = Kolmogorov-Smirnov), 1 (weighted), and 2 (over-weighted). When the score type is 1 or 2 it is
  # necessary to input the correlation vector with the values in the same order as in the gene list.
  #
  # Inputs:
  #   gene.list: The ordered gene list (e.g. integers indicating the original position in the input dataset)
  #   gene.set: A gene set (e.g. integers indicating the location of those genes in the input dataset)
  #   weighted.score.type: Type of score: weight: 0 (unweighted = Kolmogorov-Smirnov), 1 (weighted), and 2 (over-weighted)
  #  correl.vector: A vector with the coorelations (e.g. signal to noise scores) corresponding to the genes in the gene list
  #
  # Outputs:
  #   ES: Enrichment score (real number between -1 and +1)
  #   arg.ES: Location in gene.list where the peak running enrichment occurs (peak of the "mountain")
  #   RES: Numerical vector containing the running enrichment score for all locations in the gene list
  #   tag.indicator: Binary vector indicating the location of the gene sets (1's) in the gene list
  #
  # The Broad Institute
  # SOFTWARE COPYRIGHT NOTICE AGREEMENT
  # This software and its documentation are copyright 2003 by the
  # Broad Institute/Massachusetts Institute of Technology.
  # All rights are reserved.
  #
  # This software is supplied without any warranty or guaranteed support
  # whatsoever. Neither the Broad Institute nor MIT can be responsible for
  # its use, misuse, or functionality.

  tag.indicator <- sign(match(gene.list, gene.set, nomatch = 0)) # notice that the sign is 0 (no tag) or 1 (tag)
  no.tag.indicator <- 1 - tag.indicator
  N <- length(gene.list)
  Nh <- length(gene.set)
  Nm <- N - Nh
  if (weighted.score.type == 0) {
    correl.vector <- rep(1, N)
  }
  alpha <- weighted.score.type
  correl.vector <- abs(correl.vector**alpha)
  sum.correl.tag <- sum(correl.vector[tag.indicator == 1])
  norm.tag <- 1.0 / sum.correl.tag
  norm.no.tag <- 1.0 / Nm
  RES <- cumsum(tag.indicator * correl.vector * norm.tag - no.tag.indicator * norm.no.tag)

  ES <- sum(RES)

  return(list(ES = ES))
}


GSEA.HeatMapPlot <- function(V, row.names = F, col.labels, col.classes, col.names = F, main = " ", xlab = " ", ylab = " ") {
  #
  # Plots a heatmap "pinkogram" of a gene expression matrix including phenotype vector and gene, sample and phenotype labels
  #
  # The Broad Institute
  # SOFTWARE COPYRIGHT NOTICE AGREEMENT
  # This software and its documentation are copyright 2003 by the
  # Broad Institute/Massachusetts Institute of Technology.
  # All rights are reserved.
  #
  # This software is supplied without any warranty or guaranteed support
  # whatsoever. Neither the Broad Institute nor MIT can be responsible for
  # its use, misuse, or functionality.

  n.rows <- length(V[, 1])
  n.cols <- length(V[1, ])
  row.mean <- apply(V, MARGIN = 1, FUN = mean)
  row.sd <- apply(V, MARGIN = 1, FUN = sd)
  row.n <- length(V[, 1])
  for (i in 1:n.rows) {
    if (row.sd[i] == 0) {
      V[i, ] <- 0
    } else {
      V[i, ] <- (V[i, ] - row.mean[i]) / (0.5 * row.sd[i])
    }
    V[i, ] <- ifelse(V[i, ] < -6, -6, V[i, ])
    V[i, ] <- ifelse(V[i, ] > 6, 6, V[i, ])
  }

  mycol <- c("#0000FF", "#0000FF", "#4040FF", "#7070FF", "#8888FF", "#A9A9FF", "#D5D5FF", "#EEE5EE", "#FFAADA", "#FF9DB0", "#FF7080", "#FF5A5A", "#FF4040", "#FF0D1D", "#FF0000") # blue-pinkogram colors. The first and last are the colors to indicate the class vector (phenotype). This is the 1998-vintage, pre-gene cluster, original pinkogram color map

  mid.range.V <- mean(range(V)) - 0.1


  heatm <- matrix(0, nrow = n.rows + 1, ncol = n.cols)
  heatm[1:n.rows, ] <- V[seq(n.rows, 1, -1), ]
  heatm[n.rows + 1, ] <- ifelse(col.labels == 0, 7, -7)
  #        par(mar = c(5, 10, 5, 5))
  image(1:n.cols, 1:(n.rows + 1), t(heatm), col = mycol, axes = FALSE, main = main, xlab = xlab, ylab = ylab)

  if (length(row.names) > 1) {
    numC <- nchar(row.names)
    size.row.char <- 35 / (n.rows + 5)
    size.col.char <- 25 / (n.cols + 5)
    for (i in 1:n.rows) {
      row.names[i] <- substr(row.names[i], 1, 10)
    }
    row.names <- c(row.names[seq(n.rows, 1, -1)], "Class")
    axis(2, at = 1:(n.rows + 1), labels = row.names, adj = 0.5, tick = FALSE, las = 1, cex.axis = size.row.char, font.axis = 2, line = -1)
  }

  if (length(col.names) > 1) {
    axis(1, at = 1:n.cols, labels = col.names, tick = FALSE, las = 3, cex.axis = size.col.char, font.axis = 2, line = -1)
  }

  C <- split(col.labels, col.labels)
  class1.size <- length(C[[1]])
  class2.size <- length(C[[2]])
  axis(3, at = c(floor(class1.size / 2), class1.size + floor(class2.size / 2)), labels = col.classes, tick = FALSE, las = 1, cex.axis = 1.25, font.axis = 2, line = -1)

  return()
}

GSEA.Res2Frame <- function(filename = "NULL") {
  ## Reads a gene expression dataset in RES format and converts it into an R data frame
  #
  # The Broad Institute
  # SOFTWARE COPYRIGHT NOTICE AGREEMENT
  # This software and its documentation are copyright 2003 by the
  # Broad Institute/Massachusetts Institute of Technology.
  # All rights are reserved.
  #
  # This software is supplied without any warranty or guaranteed support
  # whatsoever. Neither the Broad Institute nor MIT can be responsible for
  # its use, misuse, or functionality.

  header.cont <- readLines(filename, n = 1)
  temp <- unlist(strsplit(header.cont, "\t"))
  colst <- length(temp)
  header.labels <- temp[seq(3, colst, 2)]
  ds <- read.delim(filename, header = F, row.names = 2, sep = "\t", skip = 3, blank.lines.skip = T, comment.char = "", as.is = T)
  colst <- length(ds[1, ])
  cols <- (colst - 1) / 2
  rows <- length(ds[, 1])
  A <- matrix(nrow = rows - 1, ncol = cols)
  A <- ds[1:rows, seq(2, colst, 2)]
  table1 <- data.frame(A)
  names(table1) <- header.labels
  return(table1)
}

MSIG.Res2Frame <- function(filename = "NULL") {
  header.cont <- readLines(filename, n = 1)
  temp <- unlist(strsplit(header.cont, "\t"))
  colst <- length(temp)
  header.labels <- temp[seq(3, colst, 2)]
  ds <- read.delim(filename, header = F, row.names = 2, sep = "\t", skip = 3, blank.lines.skip = T, comment.char = "", as.is = T)
  colst <- length(ds[1, ])
  cols <- (colst - 1) / 2
  rows <- length(ds[, 1])
  A <- matrix(nrow = rows - 1, ncol = cols)
  A <- ds[1:rows, seq(2, colst, 2)]
  table1 <- data.frame(A)
  names(table1) <- header.labels
  names <- header.labels
  row.names <- row.names(A)
  descs <- ds[, 1]
  return(list(ds = table1, row.names = row.names, descs = descs, names = names))
}

GSEA.Gct2Frame <- function(filename = "NULL") {
  #
  # Reads a gene expression dataset in GCT format and converts it into an R data frame
  #
  # The Broad Institute
  # SOFTWARE COPYRIGHT NOTICE AGREEMENT
  # This software and its documentation are copyright 2003 by the
  # Broad Institute/Massachusetts Institute of Technology.
  # All rights are reserved.
  #
  # This software is supplied without any warranty or guaranteed support
  # whatsoever. Neither the Broad Institute nor MIT can be responsible for
  # its use, misuse, or functionality.
  ds <- read.delim(filename, header = T, sep = "\t", skip = 2, row.names = 1, blank.lines.skip = T, comment.char = "", as.is = T)
  ds <- ds[-1]
  return(ds)
}


MSIG.Gct2Frame <- function(filename = "NULL") {
  #
  # Reads a gene expression dataset in GCT format and converts it into an R data frame
  #
  # The Broad Institute
  # SOFTWARE COPYRIGHT NOTICE AGREEMENT
  # This software and its documentation are copyright 2003 by the
  # Broad Institute/Massachusetts Institute of Technology.
  # All rights are reserved.
  #
  # This software is supplied without any warranty or guaranteed support
  # whatsoever. Neither the Broad Institute nor MIT can be responsible for
  # its use, misuse, or functionality.

  ds <- read.delim(filename, header = T, sep = "\t", skip = 2, row.names = 1, blank.lines.skip = T, comment.char = "", as.is = T, na.strings = "")
  descs <- ds[, 1]
  ds <- ds[-1]
  row.names <- row.names(ds)
  names <- names(ds)
  return(list(ds = ds, row.names = row.names, descs = descs, names = names))
}

GSEA.ReadClsFile <- function(file = "NULL") {
  #
  # Reads a class vector CLS file and defines phenotype and class labels vectors for the samples in a gene expression file (RES or GCT format)
  #
  # The Broad Institute
  # SOFTWARE COPYRIGHT NOTICE AGREEMENT
  # This software and its documentation are copyright 2003 by the
  # Broad Institute/Massachusetts Institute of Technology.
  # All rights are reserved.
  #
  # This software is supplied without any warranty or guaranteed support
  # whatsoever. Neither the Broad Institute nor MIT can be responsible for
  # its use, misuse, or functionality.

  cls.cont <- readLines(file)
  num.lines <- length(cls.cont)
  class.list <- unlist(strsplit(cls.cont[[3]], " "))
  s <- length(class.list)
  t <- table(class.list)
  l <- length(t)
  phen <- vector(length = l, mode = "character")
  phen.label <- vector(length = l, mode = "numeric")
  class.v <- vector(length = s, mode = "numeric")
  for (i in 1:l) {
    phen[i] <- noquote(names(t)[i])
    phen.label[i] <- i - 1
  }
  for (i in 1:s) {
    for (j in 1:l) {
      if (class.list[i] == phen[j]) {
        class.v[i] <- phen.label[j]
      }
    }
  }
  return(list(phen = phen, class.v = class.v))
}

GSEA.Threshold <- function(V, thres, ceil) {
  #
  # Threshold and ceiling pre-processing for gene expression matrix
  #
  # The Broad Institute
  # SOFTWARE COPYRIGHT NOTICE AGREEMENT
  # This software and its documentation are copyright 2003 by the
  # Broad Institute/Massachusetts Institute of Technology.
  # All rights are reserved.
  #
  # This software is supplied without any warranty or guaranteed support
  # whatsoever. Neither the Broad Institute nor MIT can be responsible for
  # its use, misuse, or functionality.

  V[V < thres] <- thres
  V[V > ceil] <- ceil
  return(V)
}

GSEA.VarFilter <- function(V, fold, delta, gene.names = "") {
  #
  # Variation filter pre-processing for gene expression matrix
  #
  # The Broad Institute
  # SOFTWARE COPYRIGHT NOTICE AGREEMENT
  # This software and its documentation are copyright 2003 by the
  # Broad Institute/Massachusetts Institute of Technology.
  # All rights are reserved.
  #
  # This software is supplied without any warranty or guaranteed support
  # whatsoever. Neither the Broad Institute nor MIT can be responsible for
  # its use, misuse, or functionality.

  cols <- length(V[1, ])
  rows <- length(V[, 1])
  row.max <- apply(V, MARGIN = 1, FUN = max)
  row.min <- apply(V, MARGIN = 1, FUN = min)
  flag <- array(dim = rows)
  flag <- (row.max / row.min >= fold) & (row.max - row.min >= delta)
  size <- sum(flag)
  B <- matrix(0, nrow = size, ncol = cols)
  j <- 1
  if (length(gene.names) == 1) {
    for (i in 1:rows) {
      if (flag[i]) {
        B[j, ] <- V[i, ]
        j <- j + 1
      }
    }
    return(B)
  } else {
    new.list <- vector(mode = "character", length = size)
    for (i in 1:rows) {
      if (flag[i]) {
        B[j, ] <- V[i, ]
        new.list[j] <- gene.names[i]
        j <- j + 1
      }
    }
    return(list(V = B, new.list = new.list))
  }
}

MSIG.VarFilter <- function(V, fold, delta, gene.names = "", gene.descs = "") {
  # Variation filter pre-processing for gene expression matrix

  cols <- length(V[1, ])
  rows <- length(V[, 1])
  row.max <- apply(V, MARGIN = 1, FUN = max)
  row.min <- apply(V, MARGIN = 1, FUN = min)
  flag <- array(dim = rows)
  flag <- (row.max / row.min >= fold) & (row.max - row.min >= delta)
  size <- sum(flag)
  B <- matrix(0, nrow = size, ncol = cols)
  j <- 1
  if (length(gene.names) == 1) {
    for (i in 1:rows) {
      if (flag[i]) {
        B[j, ] <- V[i, ]
        j <- j + 1
      }
    }
    return(B)
  } else {
    new.gene.names <- vector(mode = "character", length = size)
    new.gene.descs <- vector(mode = "character", length = size)
    for (i in 1:rows) {
      if (flag[i]) {
        B[j, ] <- V[i, ]
        new.gene.names[j] <- gene.names[i]
        new.gene.descs[j] <- gene.descs[i]
        j <- j + 1
      }
    }
    return(list(V = B, new.gene.names = new.gene.names, new.gene.descs = new.gene.descs, locations = flag))
  }
}

GSEA.NormalizeRows <- function(V) {
  #
  # Stardardize rows of a gene expression matrix
  #
  # The Broad Institute
  # SOFTWARE COPYRIGHT NOTICE AGREEMENT
  # This software and its documentation are copyright 2003 by the
  # Broad Institute/Massachusetts Institute of Technology.
  # All rights are reserved.
  #
  # This software is supplied without any warranty or guaranteed support
  # whatsoever. Neither the Broad Institute nor MIT can be responsible for
  # its use, misuse, or functionality.

  row.mean <- apply(V, MARGIN = 1, FUN = mean)
  row.sd <- apply(V, MARGIN = 1, FUN = sd)

  row.n <- length(V[, 1])
  for (i in 1:row.n) {
    if (row.sd[i] == 0) {
      V[i, ] <- 0
    } else {
      V[i, ] <- (V[i, ] - row.mean[i]) / row.sd[i]
    }
  }
  return(V)
}

GSEA.NormalizeCols <- function(V) {
  #
  # Stardardize columns of a gene expression matrix
  #
  # The Broad Institute
  # SOFTWARE COPYRIGHT NOTICE AGREEMENT
  # This software and its documentation are copyright 2003 by the
  # Broad Institute/Massachusetts Institute of Technology.
  # All rights are reserved.
  #
  # This software is supplied without any warranty or guaranteed support
  # whatsoever. Neither the Broad Institute nor MIT can be responsible for
  # its use, misuse, or functionality.

  col.mean <- apply(V, MARGIN = 2, FUN = mean)
  col.sd <- apply(V, MARGIN = 2, FUN = sd)
  col.n <- length(V[1, ])
  for (i in 1:col.n) {
    if (col.sd[i] == 0) {
      V[i, ] <- 0
    } else {
      V[, i] <- (V[, i] - col.mean[i]) / col.sd[i]
    }
  }
  return(V)
}

GSEA.NormalizeCols.Rank <- function(V) {
  #
  cols <- length(V[1, ])
  rows <- length(V[, 1])
  for (j in 1:cols) { # column rank normalization
    V[, j] <- rank(V[, j], ties.method = "average")
  }

  return(V)
}


MSIG.NormalizeCols.Rank <- function(V) {
  cols <- length(V[1, ])
  rows <- length(V[, 1])
  for (j in 1:cols) { # column rank normalization
    V[, j] <- rank(V[, j], ties.method = "average")
  }

  return(V)
}

MSIG.NormalizeCols.Rescale <- function(V) {
  epsilon <- 0.00001
  cols <- length(V[1, ])
  for (j in 1:cols) { # column rank normalization
    max.v <- max(V[, j])
    min.v <- min(V[, j])
    V[, j] <- (V[, j] - min.v + epsilon) / (max.v - min.v)
  }

  return(V)
}

# end of auxiliary functions

# ----------------------------------------------------------------------------------------

GSEA <- function(
    input.ds,
    input.cls,
    gene.ann = "",
    gs.db,
    gs.ann = "",
    output.directory = "",
    doc.string = "GSEA.analysis",
    non.interactive.run = F,
    reshuffling.type = "sample.labels",
    nperm = 1000,
    weighted.score.type = 1,
    nom.p.val.threshold = -1,
    fwer.p.val.threshold = -1,
    fdr.q.val.threshold = 0.25,
    topgs = 10,
    adjust.FDR.q.val = F,
    gs.size.threshold.min = 25,
    gs.size.threshold.max = 500,
    reverse.sign = F,
    preproc.type = 0,
    random.seed = 123456,
    perm.type = 0,
    fraction = 1.0,
    replace = F) {
  # This is a methodology for the analysis of global molecular profiles called Gene Set Enrichment Analysis (GSEA). It determines
  # whether an a priori defined set of genes shows statistically significant, concordant differences between two biological
  # states (e.g. phenotypes). GSEA operates on all genes from an experiment, rank ordered by the signal to noise ratio and
  # determines whether members of an a priori defined gene set are nonrandomly distributed towards the top or bottom of the
  # list and thus may correspond to an important biological process. To assess significance the program uses an empirical
  # permutation procedure to test deviation from random that preserves correlations between genes.
  #
  # For details see Subramanian et al 2005
  #
  # Inputs:
  #   input.ds: Input gene expression Affymetrix dataset file in RES or GCT format
  #   input.cls:  Input class vector (phenotype) file in CLS format
  #   gene.ann.file: Gene microarray annotation file (Affymetrix Netaffyx *.csv format) (default: none)
  #   gs.file: Gene set database in GMT format
  #   output.directory: Directory where to store output and results (default: .)
  #   doc.string:  Documentation string used as a prefix to name result files (default: "GSEA.analysis")
  #   non.interactive.run: Run in interactive (i.e. R GUI) or batch (R command line) mode (default: F)
  #   reshuffling.type: Type of permutation reshuffling: "sample.labels" or "gene.labels" (default: "sample.labels")
  #   nperm: Number of random permutations (default: 1000)
  #   weighted.score.type: Enrichment correlation-based weighting: 0=no weight (KS), 1=standard weigth, 2 = over-weigth (default: 1)
  #   nom.p.val.threshold: Significance threshold for nominal p-vals for gene sets (default: -1, no thres)
  #   fwer.p.val.threshold: Significance threshold for FWER p-vals for gene sets (default: -1, no thres)
  #   fdr.q.val.threshold: Significance threshold for FDR q-vals for gene sets (default: 0.25)
  #   topgs: Besides those passing test, number of top scoring gene sets used for detailed reports (default: 10)
  #   adjust.FDR.q.val: Adjust the FDR q-vals (default: F)
  #   gs.size.threshold.min: Minimum size (in genes) for database gene sets to be considered (default: 25)
  #   gs.size.threshold.max: Maximum size (in genes) for database gene sets to be considered (default: 500)
  #   reverse.sign: Reverse direction of gene list (pos. enrichment becomes negative, etc.) (default: F)
  #   preproc.type: Preprocessing normalization: 0=none, 1=col(z-score)., 2=col(rank) and row(z-score)., 3=col(rank). (default: 0)
  #   random.seed: Random number generator seed. (default: 123456)
  #   perm.type: Permutation type: 0 = unbalanced, 1 = balanced. For experts only (default: 0)
  #   fraction: Subsampling fraction. Set to 1.0 (no resampling). For experts only (default: 1.0)
  #   replace: Resampling mode (replacement or not replacement). For experts only (default: F)
  #
  #   Output:
  #    The results of the method are stored in the "output.directory" specified by the user as part of the input parameters.
  #      The results files are:
  #    - Two tab-separated global result text files (one for each phenotype). These files are labeled according to the doc
  #      string prefix and the phenotype name from the CLS file: <doc.string>.results.report.<phenotype>.txt
  #    - One set of global plots. They include a.- gene list correlation profile, b.- global observed and null densities, c.- heat map
  #      for the entire sorted dataset, and d.- p-values vs. NES plot. These plots are in a single JPEG file named
  #      <doc.string>.global.plots.<phenotype>.jpg. When the program is run interactively these plots appear on a window in the R GUI.
  #    - A variable number of tab-separated gene result text files according to how many sets pass any of the significance thresholds
  #      ("nom.p.val.threshold," "fwer.p.val.threshold," and "fdr.q.val.threshold") and how many are specified in the "topgs"
  #      parameter. These files are named: <doc.string>.<gene set name>.report.txt.
  #   - A variable number of gene set plots (one for each gene set report file). These plots include a.- Gene set running enrichment
  #      "mountain" plot, b.- gene set null distribution and c.- heat map for genes in the gene set. These plots are stored in a
  #      single JPEG file named <doc.string>.<gene set name>.jpg.
  #  The format (columns) for the global result files is as follows.
  #  GS : Gene set name.
  # SIZE : Size of the set in genes.
  # SOURCE : Set definition or source.
  # ES : Enrichment score.
  # NES : Normalized (multiplicative rescaling) normalized enrichment score.
  # NOM p-val : Nominal p-value (from the null distribution of the gene set).
  # FDR q-val: False discovery rate q-values
  # FWER p-val: Family wise error rate p-values.
  # Tag %: Percent of gene set before running enrichment peak.
  # Gene %: Percent of gene list before running enrichment peak.
  # Signal : enrichment signal strength.
  # FDR (median): FDR q-values from the median of the null distributions.
  # glob.p.val: P-value using a global statistic (number of sets above the set's NES).
  #
  # The rows are sorted by the NES values (from maximum positive or negative NES to minimum)
  #
  # The format (columns) for the gene set result files is as follows.
  #
  # #: Gene number in the (sorted) gene set
  # PROBE_ID : gene name. For example the accession number in the dataset.
  # SYMBOL : gene symbol from the gene annotation file.
  # DESC : gene description (title) from the gene annotation file.
  # LIST LOC : location of the gene in the sorted gene list.
  # S2N : signal to noise ratio (correlation) of the gene in the gene list.
  # RES : value of the running enrichment score at the gene location.
  # CORE_ENRICHMENT: is this gene is the "core enrichment" section of the list? Yes or No variable specifying in the gene location is before (positive ES) or after (negative ES) the running enrichment peak.
  #
  # The rows are sorted by the gene location in the gene list.
  # The function call to GSEA returns a  two element list containing the two global result reports as data frames ($report1, $report2).
  #
  # results1: Global output report for first phenotype
  # result2:  Global putput report for second phenotype
  #
  # The Broad Institute
  # SOFTWARE COPYRIGHT NOTICE AGREEMENT
  # This software and its documentation are copyright 2003 by the
  # Broad Institute/Massachusetts Institute of Technology.
  # All rights are reserved.
  #
  # This software is supplied without any warranty or guaranteed support
  # whatsoever. Neither the Broad Institute nor MIT can be responsible for
  # its use, misuse, or functionality.

  print(" *** Running GSEA Analysis...")

  # Copy input parameters to log file


  if (output.directory != "") {
    filename <- paste(output.directory, doc.string, "_params.txt", sep = "", collapse = "")

    time.string <- as.character(as.POSIXlt(Sys.time(), "GMT"))
    write(paste("Run of GSEA on ", time.string), file = filename)

    if (is.data.frame(input.ds)) {
      #      write(paste("input.ds=", quote(input.ds), sep=" "), file=filename, append=T)
    } else {
      write(paste("input.ds=", input.ds, sep = " "), file = filename, append = T)
    }
    if (is.list(input.cls)) {
      #      write(paste("input.cls=", input.cls, sep=" "), file=filename, append=T)
    } else {
      write(paste("input.cls=", input.cls, sep = " "), file = filename, append = T)
    }
    if (is.data.frame(gene.ann)) {
      #    write(paste("gene.ann =", gene.ann, sep=" "), file=filename, append=T)
    } else {
      write(paste("gene.ann =", gene.ann, sep = " "), file = filename, append = T)
    }
    if (regexpr(pattern = ".gmt", gs.db[1]) == -1) {
      #   write(paste("gs.db=", gs.db, sep=" "), file=filename, append=T)
    } else {
      write(paste("gs.db=", gs.db, sep = " "), file = filename, append = T)
    }
    if (is.data.frame(gs.ann)) {
      #    write(paste("gene.ann =", gene.ann, sep=" "), file=filename, append=T)
    } else {
      write(paste("gs.ann =", gs.ann, sep = " "), file = filename, append = T)
    }
    write(paste("output.directory =", output.directory, sep = " "), file = filename, append = T)
    write(paste("doc.string = ", doc.string, sep = " "), file = filename, append = T)
    write(paste("non.interactive.run =", non.interactive.run, sep = " "), file = filename, append = T)
    write(paste("reshuffling.type =", reshuffling.type, sep = " "), file = filename, append = T)
    write(paste("nperm =", nperm, sep = " "), file = filename, append = T)
    write(paste("weighted.score.type =", weighted.score.type, sep = " "), file = filename, append = T)
    write(paste("nom.p.val.threshold =", nom.p.val.threshold, sep = " "), file = filename, append = T)
    write(paste("fwer.p.val.threshold =", fwer.p.val.threshold, sep = " "), file = filename, append = T)
    write(paste("fdr.q.val.threshold =", fdr.q.val.threshold, sep = " "), file = filename, append = T)
    write(paste("topgs =", topgs, sep = " "), file = filename, append = T)
    write(paste("adjust.FDR.q.val =", adjust.FDR.q.val, sep = " "), file = filename, append = T)
    write(paste("gs.size.threshold.min =", gs.size.threshold.min, sep = " "), file = filename, append = T)
    write(paste("gs.size.threshold.max =", gs.size.threshold.max, sep = " "), file = filename, append = T)
    write(paste("reverse.sign =", reverse.sign, sep = " "), file = filename, append = T)
    write(paste("preproc.type =", preproc.type, sep = " "), file = filename, append = T)
    write(paste("random.seed =", random.seed, sep = " "), file = filename, append = T)
    write(paste("perm.type =", perm.type, sep = " "), file = filename, append = T)
    write(paste("fraction =", fraction, sep = " "), file = filename, append = T)
    write(paste("replace =", replace, sep = " "), file = filename, append = T)
  }

  # Start of GSEA methodology

  if (.Platform$OS.type == "windows") {
    memory.limit(6000000000)
    memory.limit()
    #      print(c("Start memory size=",  memory.size()))
  }

  # Read input data matrix

  set.seed(seed = random.seed, kind = NULL)
  adjust.param <- 0.5

  gc()

  time1 <- proc.time()

  if (is.data.frame(input.ds)) {
    dataset <- input.ds
  } else {
    if (regexpr(pattern = ".gct", input.ds) == -1) {
      dataset <- GSEA.Res2Frame(filename = input.ds)
    } else {
      dataset <- GSEA.Gct2Frame(filename = input.ds)
    }
  }
  gene.labels <- row.names(dataset)
  sample.names <- names(dataset)
  A <- data.matrix(dataset)
  dim(A)
  cols <- length(A[1, ])
  rows <- length(A[, 1])

  # preproc.type control the type of pre-processing: threshold, variation filter, normalization

  if (preproc.type == 1) { # Column normalize (Z-score)
    A <- GSEA.NormalizeCols(A)
  } else if (preproc.type == 2) { # Column (rank) and row (Z-score) normalize
    for (j in 1:cols) { # column rank normalization
      A[, j] <- rank(A[, j])
    }
    A <- GSEA.NormalizeRows(A)
  } else if (preproc.type == 3) { # Column (rank) norm.
    for (j in 1:cols) { # column rank normalization
      A[, j] <- rank(A[, j])
    }
  }

  # Read input class vector

  if (is.list(input.cls)) {
    CLS <- input.cls
  } else {
    CLS <- GSEA.ReadClsFile(file = input.cls)
  }
  class.labels <- CLS$class.v
  class.phen <- CLS$phen

  if (reverse.sign == T) {
    phen1 <- class.phen[2]
    phen2 <- class.phen[1]
  } else {
    phen1 <- class.phen[1]
    phen2 <- class.phen[2]
  }

  # sort samples according to phenotype

  col.index <- order(class.labels, decreasing = F)
  class.labels <- class.labels[col.index]
  sample.names <- sample.names[col.index]
  for (j in 1:rows) {
    A[j, ] <- A[j, col.index]
  }
  names(A) <- sample.names

  # Read input gene set database

  if (regexpr(pattern = ".gmt", gs.db[1]) == -1) {
    temp <- gs.db
  } else {
    temp <- readLines(gs.db)
  }

  max.Ng <- length(temp)
  temp.size.G <- vector(length = max.Ng, mode = "numeric")
  for (i in 1:max.Ng) {
    temp.size.G[i] <- length(unlist(strsplit(temp[[i]], "\t"))) - 2
  }

  max.size.G <- max(temp.size.G)
  gs <- matrix(rep("null", max.Ng * max.size.G), nrow = max.Ng, ncol = max.size.G)
  temp.names <- vector(length = max.Ng, mode = "character")
  temp.desc <- vector(length = max.Ng, mode = "character")
  gs.count <- 1
  for (i in 1:max.Ng) {
    gene.set.size <- length(unlist(strsplit(temp[[i]], "\t"))) - 2
    gs.line <- noquote(unlist(strsplit(temp[[i]], "\t")))
    gene.set.name <- gs.line[1]
    gene.set.desc <- gs.line[1]
    gene.set.tags <- vector(length = gene.set.size, mode = "character")
    for (j in 1:gene.set.size) {
      gene.set.tags[j] <- gs.line[j + 2]
    }
    existing.set <- is.element(gene.set.tags, gene.labels)
    set.size <- length(existing.set[existing.set == T])
    if ((set.size < gs.size.threshold.min) || (set.size > gs.size.threshold.max)) next
    temp.size.G[gs.count] <- set.size
    gs[gs.count, ] <- c(gene.set.tags[existing.set], rep("null", max.size.G - temp.size.G[gs.count]))
    temp.names[gs.count] <- gene.set.name
    temp.desc[gs.count] <- gene.set.desc
    gs.count <- gs.count + 1
  }
  Ng <- gs.count - 1
  gs.names <- vector(length = Ng, mode = "character")
  gs.desc <- vector(length = Ng, mode = "character")
  size.G <- vector(length = Ng, mode = "numeric")
  gs.names <- temp.names[1:Ng]
  gs.desc <- temp.desc[1:Ng]
  size.G <- temp.size.G[1:Ng]

  N <- length(A[, 1])
  Ns <- length(A[1, ])

  print(c("Number of genes:", N))
  print(c("Number of Gene Sets:", Ng))
  print(c("Number of samples:", Ns))
  print(c("Original number of Gene Sets:", max.Ng))
  print(c("Maximum gene set size:", max.size.G))

  # Read gene and gene set annotations if gene annotation file was provided

  all.gene.descs <- vector(length = N, mode = "character")
  all.gene.symbols <- vector(length = N, mode = "character")
  all.gs.descs <- vector(length = Ng, mode = "character")

  if (is.data.frame(gene.ann)) {
    temp <- gene.ann
    a.size <- length(temp[, 1])
    print(c("Number of gene annotation file entries:", a.size))
    accs <- as.character(temp[, 1])
    locs <- match(gene.labels, accs)
    all.gene.descs <- as.character(temp[locs, "Gene.Title"])
    all.gene.symbols <- as.character(temp[locs, "Gene.Symbol"])
    rm(temp)
  } else if (gene.ann == "") {
    for (i in 1:N) {
      all.gene.descs[i] <- " -- "
      all.gene.symbols[i] <- " -- "
    }
  } else {
    temp <- read.delim(gene.ann, header = T, sep = ",", comment.char = "", as.is = T)
    a.size <- length(temp[, 1])
    print(c("Number of gene annotation file entries:", a.size))
    accs <- as.character(temp[, 1])
    locs <- match(gene.labels, accs)
    all.gene.descs <- as.character(temp[locs, "Gene.Title"])
    all.gene.symbols <- as.character(temp[locs, "Gene.Symbol"])
    rm(temp)
  }

  if (is.data.frame(gs.ann)) {
    temp <- gs.ann
    a.size <- length(temp[, 1])
    print(c("Number of gene set annotation file entries:", a.size))
    accs <- as.character(temp[, 1])
    locs <- match(gs.names, accs)
    all.gs.descs <- as.character(temp[locs, "SOURCE"])
    rm(temp)
  } else if (gs.ann == "") {
    for (i in 1:N) {
      all.gene.descs[i] <- " -- "
      all.gene.symbols[i] <- " -- "
    }
  } else {
    temp <- read.delim(gs.ann, header = T, sep = "\t", comment.char = "", as.is = T)
    a.size <- length(temp[, 1])
    print(c("Number of gene set annotation file entries:", a.size))
    accs <- as.character(temp[, 1])
    locs <- match(gs.names, accs)
    all.gs.descs <- as.character(temp[locs, "SOURCE"])
    rm(temp)
  }


  Obs.indicator <- matrix(nrow = Ng, ncol = N)
  Obs.RES <- matrix(nrow = Ng, ncol = N)

  Obs.ES <- vector(length = Ng, mode = "numeric")
  Obs.arg.ES <- vector(length = Ng, mode = "numeric")
  Obs.ES.norm <- vector(length = Ng, mode = "numeric")

  time2 <- proc.time()

  # GSEA methodology

  # Compute observed and random permutation gene rankings

  obs.s2n <- vector(length = N, mode = "numeric")
  signal.strength <- vector(length = Ng, mode = "numeric")
  tag.frac <- vector(length = Ng, mode = "numeric")
  gene.frac <- vector(length = Ng, mode = "numeric")
  coherence.ratio <- vector(length = Ng, mode = "numeric")
  obs.phi.norm <- matrix(nrow = Ng, ncol = nperm)
  correl.matrix <- matrix(nrow = N, ncol = nperm)
  obs.correl.matrix <- matrix(nrow = N, ncol = nperm)
  order.matrix <- matrix(nrow = N, ncol = nperm)
  obs.order.matrix <- matrix(nrow = N, ncol = nperm)

  nperm.per.call <- 200
  n.groups <- nperm %/% nperm.per.call
  n.rem <- nperm %% nperm.per.call
  n.perms <- c(rep(nperm.per.call, n.groups), n.rem)
  n.ends <- cumsum(n.perms)
  n.starts <- n.ends - n.perms + 1

  if (n.rem == 0) {
    n.tot <- n.groups
  } else {
    n.tot <- n.groups + 1
  }

  for (nk in 1:n.tot) {
    call.nperm <- n.perms[nk]

    print(paste("Computing ranked list for actual and permuted phenotypes.......permutations: ", n.starts[nk], "--", n.ends[nk], sep = " "))

    O <- GSEA.GeneRanking(A, class.labels, gene.labels, call.nperm, permutation.type = perm.type, sigma.correction = "GeneCluster", fraction = fraction, replace = replace, reverse.sign = reverse.sign)
    gc()

    order.matrix[, n.starts[nk]:n.ends[nk]] <- O$order.matrix
    obs.order.matrix[, n.starts[nk]:n.ends[nk]] <- O$obs.order.matrix
    correl.matrix[, n.starts[nk]:n.ends[nk]] <- O$s2n.matrix
    obs.correl.matrix[, n.starts[nk]:n.ends[nk]] <- O$obs.s2n.matrix
    rm(O)
  }

  obs.s2n <- apply(obs.correl.matrix, 1, median) # using median to assign enrichment scores
  obs.index <- order(obs.s2n, decreasing = T)
  obs.s2n <- sort(obs.s2n, decreasing = T)

  obs.gene.labels <- gene.labels[obs.index]
  obs.gene.descs <- all.gene.descs[obs.index]
  obs.gene.symbols <- all.gene.symbols[obs.index]

  for (r in 1:nperm) {
    correl.matrix[, r] <- correl.matrix[order.matrix[, r], r]
  }
  for (r in 1:nperm) {
    obs.correl.matrix[, r] <- obs.correl.matrix[obs.order.matrix[, r], r]
  }

  gene.list2 <- obs.index
  for (i in 1:Ng) {
    print(paste("Computing observed enrichment for gene set:", i, gs.names[i], sep = " "))
    gene.set <- gs[i, gs[i, ] != "null"]
    gene.set2 <- vector(length = length(gene.set), mode = "numeric")
    gene.set2 <- match(gene.set, gene.labels)
    GSEA.results <- GSEA.EnrichmentScore(gene.list = gene.list2, gene.set = gene.set2, weighted.score.type = weighted.score.type, correl.vector = obs.s2n)
    Obs.ES[i] <- GSEA.results$ES
    Obs.arg.ES[i] <- GSEA.results$arg.ES
    Obs.RES[i, ] <- GSEA.results$RES
    Obs.indicator[i, ] <- GSEA.results$indicator
    if (Obs.ES[i] >= 0) { # compute signal strength
      tag.frac[i] <- sum(Obs.indicator[i, 1:Obs.arg.ES[i]]) / size.G[i]
      gene.frac[i] <- Obs.arg.ES[i] / N
    } else {
      tag.frac[i] <- sum(Obs.indicator[i, Obs.arg.ES[i]:N]) / size.G[i]
      gene.frac[i] <- (N - Obs.arg.ES[i] + 1) / N
    }
    signal.strength[i] <- tag.frac[i] * (1 - gene.frac[i]) * (N / (N - size.G[i]))
  }

  # Compute enrichment for random permutations

  phi <- matrix(nrow = Ng, ncol = nperm)
  phi.norm <- matrix(nrow = Ng, ncol = nperm)
  obs.phi <- matrix(nrow = Ng, ncol = nperm)

  if (reshuffling.type == "sample.labels") { # reshuffling phenotype labels
    for (i in 1:Ng) {
      print(paste("Computing random permutations' enrichment for gene set:", i, gs.names[i], sep = " "))
      gene.set <- gs[i, gs[i, ] != "null"]
      gene.set2 <- vector(length = length(gene.set), mode = "numeric")
      gene.set2 <- match(gene.set, gene.labels)
      for (r in 1:nperm) {
        gene.list2 <- order.matrix[, r]
        GSEA.results <- GSEA.EnrichmentScore(gene.list = gene.list2, gene.set = gene.set2, weighted.score.type = weighted.score.type, correl.vector = correl.matrix[, r])
        phi[i, r] <- GSEA.results$ES
        obs.gene.list2 <- obs.order.matrix[, r]
        GSEA.results <- GSEA.EnrichmentScore(gene.list = obs.gene.list2, gene.set = gene.set2, weighted.score.type = weighted.score.type, correl.vector = obs.correl.matrix[, r])
        obs.phi[i, r] <- GSEA.results$ES
      }
      gc()
    }
  } else if (reshuffling.type == "gene.labels") { # reshuffling gene labels
    for (i in 1:Ng) {
      gene.set <- gs[i, gs[i, ] != "null"]
      gene.set2 <- vector(length = length(gene.set), mode = "numeric")
      gene.set2 <- match(gene.set, gene.labels)
      for (r in 1:nperm) {
        reshuffled.gene.labels <- sample(1:rows)
        GSEA.results <- GSEA.EnrichmentScore(gene.list = reshuffled.gene.labels, gene.set = gene.set2, weighted.score.type = weighted.score.type, correl.vector = obs.s2n)
        phi[i, r] <- GSEA.results$ES
        obs.gene.list2 <- obs.order.matrix[, r]
        GSEA.results <- GSEA.EnrichmentScore(gene.list = obs.gene.list2, gene.set = gene.set2, weighted.score.type = weighted.score.type, correl.vector = obs.correl.matrix[, r])
        obs.phi[i, r] <- GSEA.results$ES
      }
      gc()
    }
  }

  # Compute 3 types of p-values

  # Find nominal p-values


  print("Computing nominal p-values...")

  p.vals <- matrix(0, nrow = Ng, ncol = 2)

  for (i in 1:Ng) {
    pos.phi <- phi[i, phi[i, ] >= 0]
    neg.phi <- phi[i, phi[i, ] < 0]
    if (Obs.ES[i] >= 0) {
      p.vals[i, 1] <- sum(pos.phi >= Obs.ES[i]) / length(pos.phi)
      p.vals[i, 1] <- signif(p.vals[i, 1], digits = 5)
    } else {
      p.vals[i, 1] <- sum(neg.phi < Obs.ES[i]) / length(neg.phi)
      p.vals[i, 1] <- signif(p.vals[i, 1], digits = 5)
    }
  }

  # Find effective size

  erf <- function(x) {
    2 * pnorm(sqrt(2) * x)
  }

  KS.mean <- function(N) { # KS mean as a function of set size N
    S <- 0
    for (k in -100:100) {
      if (k == 0) next
      S <- S + 4 * (-1)**(k + 1) * (0.25 * exp(-2 * k * k * N) - sqrt(2 * pi) * erf(sqrt(2 * N) * k) / (16 * k * sqrt(N)))
    }
    return(abs(S))
  }

  KS.mean.table <- vector(length = 5000, mode = "numeric")

  for (i in 1:5000) {
    KS.mean.table[i] <- KS.mean(i)
  }

  KS.size <- vector(length = Ng, mode = "numeric")

  # Rescaling normalization for each gene set null

  print("Computing rescaling normalization for each gene set null...")

  for (i in 1:Ng) {
    pos.phi <- phi[i, phi[i, ] >= 0]
    neg.phi <- phi[i, phi[i, ] < 0]
    pos.m <- mean(pos.phi)
    neg.m <- mean(abs(neg.phi))
    if (Obs.ES[i] >= 0) {
      KS.size[i] <- which.min(abs(KS.mean.table - pos.m))
    } else {
      KS.size[i] <- which.min(abs(KS.mean.table - neg.m))
    }
    pos.phi <- pos.phi / pos.m
    neg.phi <- neg.phi / neg.m
    for (j in 1:nperm) {
      if (phi[i, j] >= 0) {
        phi.norm[i, j] <- phi[i, j] / pos.m
      } else {
        phi.norm[i, j] <- phi[i, j] / neg.m
      }
    }
    for (j in 1:nperm) {
      if (obs.phi[i, j] >= 0) {
        obs.phi.norm[i, j] <- obs.phi[i, j] / pos.m
      } else {
        obs.phi.norm[i, j] <- obs.phi[i, j] / neg.m
      }
    }
    if (Obs.ES[i] >= 0) {
      Obs.ES.norm[i] <- Obs.ES[i] / pos.m
    } else {
      Obs.ES.norm[i] <- Obs.ES[i] / neg.m
    }
  }

  # Compute FWER p-vals

  print("Computing FWER p-values...")

  max.ES.vals.p <- NULL
  max.ES.vals.n <- NULL
  for (j in 1:nperm) {
    pos.phi <- phi.norm[phi.norm[, j] >= 0, j]
    neg.phi <- phi.norm[phi.norm[, j] < 0, j]
    if (length(pos.phi) > 0) {
      max.ES.vals.p <- c(max.ES.vals.p, max(pos.phi))
    }
    if (length(neg.phi) > 0) {
      max.ES.vals.n <- c(max.ES.vals.n, min(neg.phi))
    }
  }
  for (i in 1:Ng) {
    if (Obs.ES.norm[i] >= 0) {
      p.vals[i, 2] <- sum(max.ES.vals.p >= Obs.ES.norm[i]) / length(max.ES.vals.p)
    } else {
      p.vals[i, 2] <- sum(max.ES.vals.n < Obs.ES.norm[i]) / length(max.ES.vals.n)
    }
    p.vals[i, 2] <- signif(p.vals[i, 2], digits = 4)
  }

  # Compute FDRs

  print("Computing FDR q-values...")

  NES <- vector(length = Ng, mode = "numeric")
  phi.norm.mean <- vector(length = Ng, mode = "numeric")
  obs.phi.norm.mean <- vector(length = Ng, mode = "numeric")
  phi.norm.median <- vector(length = Ng, mode = "numeric")
  obs.phi.norm.median <- vector(length = Ng, mode = "numeric")
  phi.norm.mean <- vector(length = Ng, mode = "numeric")
  obs.phi.mean <- vector(length = Ng, mode = "numeric")
  FDR.mean <- vector(length = Ng, mode = "numeric")
  FDR.median <- vector(length = Ng, mode = "numeric")
  phi.norm.median.d <- vector(length = Ng, mode = "numeric")
  obs.phi.norm.median.d <- vector(length = Ng, mode = "numeric")

  Obs.ES.index <- order(Obs.ES.norm, decreasing = T)
  Orig.index <- seq(1, Ng)
  Orig.index <- Orig.index[Obs.ES.index]
  Orig.index <- order(Orig.index, decreasing = F)
  Obs.ES.norm.sorted <- Obs.ES.norm[Obs.ES.index]
  gs.names.sorted <- gs.names[Obs.ES.index]

  for (k in 1:Ng) {
    NES[k] <- Obs.ES.norm.sorted[k]
    count.col <- vector(length = nperm, mode = "numeric")
    obs.count.col <- vector(length = nperm, mode = "numeric")
    for (i in 1:nperm) {
      if (NES[k] >= 0) {
        count.col.norm <- sum(phi.norm[, i] >= 0)
        obs.count.col.norm <- sum(obs.phi.norm[, i] >= 0)
        count.col[i] <- sum(phi.norm[, i] >= NES[k]) / count.col.norm
        obs.count.col[i] <- sum(obs.phi.norm[, i] >= NES[k]) / obs.count.col.norm
      } else {
        count.col.norm <- sum(phi.norm[, i] < 0)
        obs.count.col.norm <- sum(obs.phi.norm[, i] < 0)
        count.col[i] <- sum(phi.norm[, i] <= NES[k]) / count.col.norm
        obs.count.col[i] <- sum(obs.phi.norm[, i] <= NES[k]) / obs.count.col.norm
      }
    }
    phi.norm.mean[k] <- mean(count.col)
    obs.phi.norm.mean[k] <- mean(obs.count.col)
    phi.norm.median[k] <- median(count.col)
    obs.phi.norm.median[k] <- median(obs.count.col)
    FDR.mean[k] <- ifelse(phi.norm.mean[k] / obs.phi.norm.mean[k] < 1, phi.norm.mean[k] / obs.phi.norm.mean[k], 1)
    FDR.median[k] <- ifelse(phi.norm.median[k] / obs.phi.norm.median[k] < 1, phi.norm.median[k] / obs.phi.norm.median[k], 1)
  }

  # adjust q-values

  if (adjust.FDR.q.val == T) {
    pos.nes <- length(NES[NES >= 0])
    min.FDR.mean <- FDR.mean[pos.nes]
    min.FDR.median <- FDR.median[pos.nes]
    for (k in seq(pos.nes - 1, 1, -1)) {
      if (FDR.mean[k] < min.FDR.mean) {
        min.FDR.mean <- FDR.mean[k]
      }
      if (min.FDR.mean < FDR.mean[k]) {
        FDR.mean[k] <- min.FDR.mean
      }
    }

    neg.nes <- pos.nes + 1
    min.FDR.mean <- FDR.mean[neg.nes]
    min.FDR.median <- FDR.median[neg.nes]
    for (k in seq(neg.nes + 1, Ng)) {
      if (FDR.mean[k] < min.FDR.mean) {
        min.FDR.mean <- FDR.mean[k]
      }
      if (min.FDR.mean < FDR.mean[k]) {
        FDR.mean[k] <- min.FDR.mean
      }
    }
  }

  obs.phi.norm.mean.sorted <- obs.phi.norm.mean[Orig.index]
  phi.norm.mean.sorted <- phi.norm.mean[Orig.index]
  FDR.mean.sorted <- FDR.mean[Orig.index]
  FDR.median.sorted <- FDR.median[Orig.index]

  #   Compute global statistic

  glob.p.vals <- vector(length = Ng, mode = "numeric")
  NULL.pass <- vector(length = nperm, mode = "numeric")
  OBS.pass <- vector(length = nperm, mode = "numeric")

  for (k in 1:Ng) {
    NES[k] <- Obs.ES.norm.sorted[k]
    if (NES[k] >= 0) {
      for (i in 1:nperm) {
        NULL.pass[i] <- sum(phi.norm[, i] >= NES[k]) / sum(phi.norm[, i] >= 0)
        OBS.pass[i] <- sum(obs.phi.norm[, i] >= NES[k]) / sum(obs.phi.norm[, i] >= 0)
      }
    } else {
      for (i in 1:nperm) {
        NULL.pass[i] <- sum(phi.norm[, i] <= NES[k]) / sum(phi.norm[, i] < 0)
        OBS.pass[i] <- sum(obs.phi.norm[, i] <= NES[k]) / sum(obs.phi.norm[, i] < 0)
      }
    }
    glob.p.vals[k] <- sum(NULL.pass >= mean(OBS.pass)) / nperm
  }
  glob.p.vals.sorted <- glob.p.vals[Orig.index]

  # Produce results report

  print("Producing result tables and plots...")

  Obs.ES <- signif(Obs.ES, digits = 5)
  Obs.ES.norm <- signif(Obs.ES.norm, digits = 5)
  p.vals <- signif(p.vals, digits = 4)
  signal.strength <- signif(signal.strength, digits = 3)
  tag.frac <- signif(tag.frac, digits = 3)
  gene.frac <- signif(gene.frac, digits = 3)
  FDR.mean.sorted <- signif(FDR.mean.sorted, digits = 5)
  FDR.median.sorted <- signif(FDR.median.sorted, digits = 5)
  glob.p.vals.sorted <- signif(glob.p.vals.sorted, digits = 5)

  report <- data.frame(cbind(gs.names, size.G, all.gs.descs, Obs.ES, Obs.ES.norm, p.vals[, 1], FDR.mean.sorted, p.vals[, 2], tag.frac, gene.frac, signal.strength, FDR.median.sorted, glob.p.vals.sorted))
  names(report) <- c("GS", "SIZE", "SOURCE", "ES", "NES", "NOM p-val", "FDR q-val", "FWER p-val", "Tag %", "Gene %", "Signal", "FDR (median)", "glob.p.val")
  #       print(report)
  report2 <- report
  report.index2 <- order(Obs.ES.norm, decreasing = T)
  for (i in 1:Ng) {
    report2[i, ] <- report[report.index2[i], ]
  }
  report3 <- report
  report.index3 <- order(Obs.ES.norm, decreasing = F)
  for (i in 1:Ng) {
    report3[i, ] <- report[report.index3[i], ]
  }
  phen1.rows <- length(Obs.ES.norm[Obs.ES.norm >= 0])
  phen2.rows <- length(Obs.ES.norm[Obs.ES.norm < 0])
  report.phen1 <- report2[1:phen1.rows, ]
  report.phen2 <- report3[1:phen2.rows, ]

  if (output.directory != "") {
    if (phen1.rows > 0) {
      filename <- paste(output.directory, doc.string, ".results.report.", phen1, ".txt", sep = "", collapse = "")
      write.table(report.phen1, file = filename, quote = F, row.names = F, sep = "\t")
    }
    if (phen2.rows > 0) {
      filename <- paste(output.directory, doc.string, ".results.report.", phen2, ".txt", sep = "", collapse = "")
      write.table(report.phen2, file = filename, quote = F, row.names = F, sep = "\t")
    }
  }

  # Global plots

  if (output.directory != "") {
    if (non.interactive.run == F) {
      if (.Platform$OS.type == "windows") {
        glob.filename <- paste(output.directory, doc.string, ".global.plots", sep = "", collapse = "")
        x11(width = 10, height = 10)
      } else if (.Platform$OS.type == "unix") {
        glob.filename <- paste(output.directory, doc.string, ".global.plots.pdf", sep = "", collapse = "")
        pdf(file = glob.filename, height = 10, width = 10)
      }
    } else {
      if (.Platform$OS.type == "unix") {
        glob.filename <- paste(output.directory, doc.string, ".global.plots.pdf", sep = "", collapse = "")
        pdf(file = glob.filename, height = 10, width = 10)
      } else if (.Platform$OS.type == "windows") {
        glob.filename <- paste(output.directory, doc.string, ".global.plots.pdf", sep = "", collapse = "")
        pdf(file = glob.filename, height = 10, width = 10)
      }
    }
  }

  nf <- layout(matrix(c(1, 2, 3, 4), 2, 2, byrow = T), c(1, 1), c(1, 1), TRUE)

  # plot S2N correlation profile

  location <- 1:N
  max.corr <- max(obs.s2n)
  min.corr <- min(obs.s2n)

  x <- plot(location, obs.s2n, ylab = "Signal to Noise Ratio (S2N)", xlab = "Gene List Location", main = "Gene List Correlation (S2N) Profile", type = "l", lwd = 2, cex = 0.9, col = 1)
  for (i in seq(1, N, 20)) {
    lines(c(i, i), c(0, obs.s2n[i]), lwd = 3, cex = 0.9, col = colors()[12]) # shading of correlation plot
  }
  x <- points(location, obs.s2n, type = "l", lwd = 2, cex = 0.9, col = 1)
  lines(c(1, N), c(0, 0), lwd = 2, lty = 1, cex = 0.9, col = 1) # zero correlation horizontal line
  temp <- order(abs(obs.s2n), decreasing = T)
  arg.correl <- temp[N]
  lines(c(arg.correl, arg.correl), c(min.corr, 0.7 * max.corr), lwd = 2, lty = 3, cex = 0.9, col = 1) # zero correlation vertical line

  area.bias <- signif(100 * (sum(obs.s2n[1:arg.correl]) + sum(obs.s2n[arg.correl:N])) / sum(abs(obs.s2n[1:N])), digits = 3)
  area.phen <- ifelse(area.bias >= 0, phen1, phen2)
  delta.string <- paste("Corr. Area Bias to \"", area.phen, "\" =", abs(area.bias), "%", sep = "", collapse = "")
  zero.crossing.string <- paste("Zero Crossing at location ", arg.correl, " (", signif(100 * arg.correl / N, digits = 3), " %)")
  leg.txt <- c(delta.string, zero.crossing.string)
  legend(x = N / 10, y = max.corr, bty = "n", bg = "white", legend = leg.txt, cex = 0.9)

  leg.txt <- paste("\"", phen1, "\" ", sep = "", collapse = "")
  text(x = 1, y = -0.05 * max.corr, adj = c(0, 1), labels = leg.txt, cex = 0.9)

  leg.txt <- paste("\"", phen2, "\" ", sep = "", collapse = "")
  text(x = N, y = 0.05 * max.corr, adj = c(1, 0), labels = leg.txt, cex = 0.9)

  if (Ng > 1) { # make these plots only if there are multiple gene sets.

    # compute plots of actual (weighted) null and observed

    phi.densities.pos <- matrix(0, nrow = 512, ncol = nperm)
    phi.densities.neg <- matrix(0, nrow = 512, ncol = nperm)
    obs.phi.densities.pos <- matrix(0, nrow = 512, ncol = nperm)
    obs.phi.densities.neg <- matrix(0, nrow = 512, ncol = nperm)
    phi.density.mean.pos <- vector(length = 512, mode = "numeric")
    phi.density.mean.neg <- vector(length = 512, mode = "numeric")
    obs.phi.density.mean.pos <- vector(length = 512, mode = "numeric")
    obs.phi.density.mean.neg <- vector(length = 512, mode = "numeric")
    phi.density.median.pos <- vector(length = 512, mode = "numeric")
    phi.density.median.neg <- vector(length = 512, mode = "numeric")
    obs.phi.density.median.pos <- vector(length = 512, mode = "numeric")
    obs.phi.density.median.neg <- vector(length = 512, mode = "numeric")
    x.coor.pos <- vector(length = 512, mode = "numeric")
    x.coor.neg <- vector(length = 512, mode = "numeric")

    for (i in 1:nperm) {
      pos.phi <- phi.norm[phi.norm[, i] >= 0, i]
      if (length(pos.phi) > 2) {
        temp <- density(pos.phi, adjust = adjust.param, n = 512, from = 0, to = 3.5)
      } else {
        temp <- list(x = 3.5 * (seq(1, 512) - 1) / 512, y = rep(0.001, 512))
      }
      phi.densities.pos[, i] <- temp$y
      norm.factor <- sum(phi.densities.pos[, i])
      phi.densities.pos[, i] <- phi.densities.pos[, i] / norm.factor
      if (i == 1) {
        x.coor.pos <- temp$x
      }

      neg.phi <- phi.norm[phi.norm[, i] < 0, i]
      if (length(neg.phi) > 2) {
        temp <- density(neg.phi, adjust = adjust.param, n = 512, from = -3.5, to = 0)
      } else {
        temp <- list(x = 3.5 * (seq(1, 512) - 1) / 512, y = rep(0.001, 512))
      }
      phi.densities.neg[, i] <- temp$y
      norm.factor <- sum(phi.densities.neg[, i])
      phi.densities.neg[, i] <- phi.densities.neg[, i] / norm.factor
      if (i == 1) {
        x.coor.neg <- temp$x
      }
      pos.phi <- obs.phi.norm[obs.phi.norm[, i] >= 0, i]
      if (length(pos.phi) > 2) {
        temp <- density(pos.phi, adjust = adjust.param, n = 512, from = 0, to = 3.5)
      } else {
        temp <- list(x = 3.5 * (seq(1, 512) - 1) / 512, y = rep(0.001, 512))
      }
      obs.phi.densities.pos[, i] <- temp$y
      norm.factor <- sum(obs.phi.densities.pos[, i])
      obs.phi.densities.pos[, i] <- obs.phi.densities.pos[, i] / norm.factor

      neg.phi <- obs.phi.norm[obs.phi.norm[, i] < 0, i]
      if (length(neg.phi) > 2) {
        temp <- density(neg.phi, adjust = adjust.param, n = 512, from = -3.5, to = 0)
      } else {
        temp <- list(x = 3.5 * (seq(1, 512) - 1) / 512, y = rep(0.001, 512))
      }
      obs.phi.densities.neg[, i] <- temp$y
      norm.factor <- sum(obs.phi.densities.neg[, i])
      obs.phi.densities.neg[, i] <- obs.phi.densities.neg[, i] / norm.factor
    }
    phi.density.mean.pos <- apply(phi.densities.pos, 1, mean)
    phi.density.mean.neg <- apply(phi.densities.neg, 1, mean)

    obs.phi.density.mean.pos <- apply(obs.phi.densities.pos, 1, mean)
    obs.phi.density.mean.neg <- apply(obs.phi.densities.neg, 1, mean)

    phi.density.median.pos <- apply(phi.densities.pos, 1, median)
    phi.density.median.neg <- apply(phi.densities.neg, 1, median)

    obs.phi.density.median.pos <- apply(obs.phi.densities.pos, 1, median)
    obs.phi.density.median.neg <- apply(obs.phi.densities.neg, 1, median)

    x <- c(x.coor.neg, x.coor.pos)
    x.plot.range <- range(x)
    y1 <- c(phi.density.mean.neg, phi.density.mean.pos)
    y2 <- c(obs.phi.density.mean.neg, obs.phi.density.mean.pos)
    y.plot.range <- c(-0.3 * max(c(y1, y2)), max(c(y1, y2)))
    plot(x, y1, xlim = x.plot.range, ylim = 1.5 * y.plot.range, type = "l", lwd = 2, col = 2, xlab = "NES", ylab = "P(NES)", main = "Global Observed and Null Densities (Area Normalized)")

    y1.point <- y1[seq(1, length(x), 2)]
    y2.point <- y2[seq(2, length(x), 2)]
    x1.point <- x[seq(1, length(x), 2)]
    x2.point <- x[seq(2, length(x), 2)]

    #     for (i in 1:length(x1.point)) {
    #       lines(c(x1.point[i], x1.point[i]), c(0, y1.point[i]), lwd = 3, cex = 0.9, col = colors()[555]) # shading
    #     }
    #
    #     for (i in 1:length(x2.point)) {
    #       lines(c(x2.point[i], x2.point[i]), c(0, y2.point[i]), lwd = 3, cex = 0.9, col = colors()[29]) # shading
    #     }

    points(x, y1, type = "l", lwd = 2, col = colors()[555])
    points(x, y2, type = "l", lwd = 2, col = colors()[29])

    for (i in 1:Ng) {
      col <- ifelse(Obs.ES.norm[i] > 0, 2, 3)
      lines(c(Obs.ES.norm[i], Obs.ES.norm[i]), c(-0.2 * max(c(y1, y2)), 0), lwd = 1, lty = 1, col = 1)
    }
    leg.txt <- paste("Neg. ES: \"", phen2, " \" ", sep = "", collapse = "")
    text(x = x.plot.range[1], y = -0.25 * max(c(y1, y2)), adj = c(0, 1), labels = leg.txt, cex = 0.9)
    leg.txt <- paste(" Pos. ES: \"", phen1, "\" ", sep = "", collapse = "")
    text(x = x.plot.range[2], y = -0.25 * max(c(y1, y2)), adj = c(1, 1), labels = leg.txt, cex = 0.9)

    leg.txt <- c("Null Density", "Observed Density", "Observed NES values")
    c.vec <- c(colors()[555], colors()[29], 1)
    lty.vec <- c(1, 1, 1)
    lwd.vec <- c(2, 2, 2)
    legend(x = 0, y = 1.5 * y.plot.range[2], bty = "n", bg = "white", legend = leg.txt, lty = lty.vec, lwd = lwd.vec, col = c.vec, cex = 0.9)

    B <- A[obs.index, ]
    GSEA.HeatMapPlot(V = B, col.labels = class.labels, col.classes = class.phen, main = "Heat Map for Genes in Dataset")

    # p-vals plot
    nom.p.vals <- p.vals[Obs.ES.index, 1]
    FWER.p.vals <- p.vals[Obs.ES.index, 2]
    plot.range <- 1.25 * range(NES)
    plot(NES, FDR.mean, ylim = c(0, 1), xlim = plot.range, col = 1, bg = 1, type = "p", pch = 22, cex = 0.75, xlab = "NES", main = "p-values vs. NES", ylab = "p-val/q-val")
    points(NES, nom.p.vals, type = "p", col = 2, bg = 2, pch = 22, cex = 0.75)
    points(NES, FWER.p.vals, type = "p", col = colors()[577], bg = colors()[577], pch = 22, cex = 0.75)
    leg.txt <- c("Nominal p-value", "FWER p-value", "FDR q-value")
    c.vec <- c(2, colors()[577], 1)
    pch.vec <- c(22, 22, 22)
    legend(x = -0.5, y = 0.5, bty = "n", bg = "white", legend = leg.txt, pch = pch.vec, col = c.vec, pt.bg = c.vec, cex = 0.9)
    lines(c(min(NES), max(NES)), c(nom.p.val.threshold, nom.p.val.threshold), lwd = 1, lty = 2, col = 2)
    lines(c(min(NES), max(NES)), c(fwer.p.val.threshold, fwer.p.val.threshold), lwd = 1, lty = 2, col = colors()[577])
    lines(c(min(NES), max(NES)), c(fdr.q.val.threshold, fdr.q.val.threshold), lwd = 1, lty = 2, col = 1)

    if (non.interactive.run == F) {
      if (.Platform$OS.type == "windows") {
        savePlot(filename = glob.filename, type = "jpeg", device = dev.cur())
      } else if (.Platform$OS.type == "unix") {
        dev.off()
      }
    } else {
      dev.off()
    }
  } # if Ng > 1

  #----------------------------------------------------------------------------
  # Produce report for each gene set passing the nominal, FWER or FDR test or the top topgs in each side

  if (topgs > floor(Ng / 2)) {
    topgs <- floor(Ng / 2)
  }

  for (i in 1:Ng) {
    if ((p.vals[i, 1] <= nom.p.val.threshold) ||
      (p.vals[i, 2] <= fwer.p.val.threshold) ||
      (FDR.mean.sorted[i] <= fdr.q.val.threshold) ||
      (is.element(i, c(Obs.ES.index[1:topgs], Obs.ES.index[(Ng - topgs + 1):Ng])))) {
      #  produce report per gene set

      kk <- 1
      gene.number <- vector(length = size.G[i], mode = "character")
      gene.names <- vector(length = size.G[i], mode = "character")
      gene.symbols <- vector(length = size.G[i], mode = "character")
      gene.descs <- vector(length = size.G[i], mode = "character")
      gene.list.loc <- vector(length = size.G[i], mode = "numeric")
      core.enrichment <- vector(length = size.G[i], mode = "character")
      gene.s2n <- vector(length = size.G[i], mode = "numeric")
      gene.RES <- vector(length = size.G[i], mode = "numeric")
      rank.list <- seq(1, N)

      if (Obs.ES[i] >= 0) {
        set.k <- seq(1, N, 1)
      } else {
        set.k <- seq(N, 1, -1)
      }

      for (k in set.k) {
        if (Obs.indicator[i, k] == 1) {
          gene.number[kk] <- kk
          gene.names[kk] <- obs.gene.labels[k]
          gene.symbols[kk] <- substr(obs.gene.symbols[k], 1, 15)
          gene.descs[kk] <- substr(obs.gene.descs[k], 1, 40)
          gene.list.loc[kk] <- k
          gene.s2n[kk] <- signif(obs.s2n[k], digits = 3)
          gene.RES[kk] <- signif(Obs.RES[i, k], digits = 3)
          if (Obs.ES[i] >= 0) {
            core.enrichment[kk] <- ifelse(gene.list.loc[kk] <= Obs.arg.ES[i], "YES", "NO")
          } else {
            core.enrichment[kk] <- ifelse(gene.list.loc[kk] > Obs.arg.ES[i], "YES", "NO")
          }
          kk <- kk + 1
        }
      }

      gene.report <- data.frame(cbind(gene.number, gene.names, gene.symbols, gene.descs, gene.list.loc, gene.s2n, gene.RES, core.enrichment))
      names(gene.report) <- c("#", "PROBE_ID", "SYMBOL", "DESC", "LIST LOC", "S2N", "RES", "CORE_ENRICHMENT")

      #       print(gene.report)

      if (output.directory != "") {
        filename <- paste(output.directory, doc.string, ".", gs.names[i], ".report.txt", sep = "", collapse = "")
        write.table(gene.report, file = filename, quote = F, row.names = F, sep = "\t")


        if (non.interactive.run == F) {
          if (.Platform$OS.type == "windows") {
            gs.filename <- paste(output.directory, doc.string, ".", gs.names[i], sep = "", collapse = "")
            x11(width = 18, height = 6)
          } else if (.Platform$OS.type == "unix") {
            gs.filename <- paste(output.directory, doc.string, ".", gs.names[i], ".pdf", sep = "", collapse = "")
            pdf(file = gs.filename, height = 6, width = 18)
          }
        } else {
          if (.Platform$OS.type == "unix") {
            gs.filename <- paste(output.directory, doc.string, ".", gs.names[i], ".pdf", sep = "", collapse = "")
            pdf(file = gs.filename, height = 6, width = 18)
          } else if (.Platform$OS.type == "windows") {
            gs.filename <- paste(output.directory, doc.string, ".", gs.names[i], ".pdf", sep = "", collapse = "")
            pdf(file = gs.filename, height = 6, width = 18)
          }
        }
      }

      nf <- layout(matrix(c(1, 2, 3), 1, 3, byrow = T), 1, c(1, 1, 1), TRUE)
      ind <- 1:N
      min.RES <- min(Obs.RES[i, ])
      max.RES <- max(Obs.RES[i, ])
      if (max.RES < 0.3) max.RES <- 0.3
      if (min.RES > -0.3) min.RES <- -0.3
      delta <- (max.RES - min.RES) * 0.50
      min.plot <- min.RES - 2 * delta
      max.plot <- max.RES
      max.corr <- max(obs.s2n)
      min.corr <- min(obs.s2n)
      Obs.correl.vector.norm <- (obs.s2n - min.corr) / (max.corr - min.corr) * 1.25 * delta + min.plot
      zero.corr.line <- (-min.corr / (max.corr - min.corr)) * 1.25 * delta + min.plot
      col <- ifelse(Obs.ES[i] > 0, 2, 4)

      # Running enrichment plot

      sub.string <- paste("Number of genes: ", N, " (in list), ", size.G[i], " (in gene set)", sep = "", collapse = "")

      main.string <- paste("Gene Set ", i, ":", gs.names[i])
      plot(ind, Obs.RES[i, ], main = main.string, sub = sub.string, xlab = "Gene List Index", ylab = "Running Enrichment Score (RES)", xlim = c(1, N), ylim = c(min.plot, max.plot), type = "l", lwd = 2, cex = 1, col = col)
      for (j in seq(1, N, 20)) {
        lines(c(j, j), c(zero.corr.line, Obs.correl.vector.norm[j]), lwd = 1, cex = 1, col = colors()[12]) # shading of correlation plot
      }
      lines(c(1, N), c(0, 0), lwd = 1, lty = 2, cex = 1, col = 1) # zero RES line
      lines(c(Obs.arg.ES[i], Obs.arg.ES[i]), c(min.plot, max.plot), lwd = 1, lty = 3, cex = 1, col = col) # max enrichment vertical line
      for (j in 1:N) {
        if (Obs.indicator[i, j] == 1) {
          lines(c(j, j), c(min.plot + 1.25 * delta, min.plot + 1.75 * delta), lwd = 1, lty = 1, cex = 1, col = 1) # enrichment tags
        }
      }
      lines(ind, Obs.correl.vector.norm, type = "l", lwd = 1, cex = 1, col = 1)
      lines(c(1, N), c(zero.corr.line, zero.corr.line), lwd = 1, lty = 1, cex = 1, col = 1) # zero correlation horizontal line
      temp <- order(abs(obs.s2n), decreasing = T)
      arg.correl <- temp[N]
      lines(c(arg.correl, arg.correl), c(min.plot, max.plot), lwd = 1, lty = 3, cex = 1, col = 3) # zero crossing correlation vertical line

      leg.txt <- paste("\"", phen1, "\" ", sep = "", collapse = "")
      text(x = 1, y = min.plot, adj = c(0, 0), labels = leg.txt, cex = 1.0)

      leg.txt <- paste("\"", phen2, "\" ", sep = "", collapse = "")
      text(x = N, y = min.plot, adj = c(1, 0), labels = leg.txt, cex = 1.0)

      adjx <- ifelse(Obs.ES[i] > 0, 0, 1)

      leg.txt <- paste("Peak at ", Obs.arg.ES[i], sep = "", collapse = "")
      text(x = Obs.arg.ES[i], y = min.plot + 1.8 * delta, adj = c(adjx, 0), labels = leg.txt, cex = 1.0)

      leg.txt <- paste("Zero crossing at ", arg.correl, sep = "", collapse = "")
      text(x = arg.correl, y = min.plot + 1.95 * delta, adj = c(adjx, 0), labels = leg.txt, cex = 1.0)

      # nominal p-val histogram

      sub.string <- paste("ES =", signif(Obs.ES[i], digits = 3), " NES =", signif(Obs.ES.norm[i], digits = 3), "Nom. p-val=", signif(p.vals[i, 1], digits = 3), "FWER=", signif(p.vals[i, 2], digits = 3), "FDR=", signif(FDR.mean.sorted[i], digits = 3))
      temp <- density(phi[i, ], adjust = adjust.param)
      x.plot.range <- range(temp$x)
      y.plot.range <- c(-0.125 * max(temp$y), 1.5 * max(temp$y))
      plot(temp$x, temp$y, type = "l", sub = sub.string, xlim = x.plot.range, ylim = y.plot.range, lwd = 2, col = 2, main = "Gene Set Null Distribution", xlab = "ES", ylab = "P(ES)")
      x.loc <- which.min(abs(temp$x - Obs.ES[i]))
      lines(c(Obs.ES[i], Obs.ES[i]), c(0, temp$y[x.loc]), lwd = 2, lty = 1, cex = 1, col = 1)
      lines(x.plot.range, c(0, 0), lwd = 1, lty = 1, cex = 1, col = 1)

      leg.txt <- c("Gene Set Null Density", "Observed Gene Set ES value")
      c.vec <- c(2, 1)
      lty.vec <- c(1, 1)
      lwd.vec <- c(2, 2)
      legend(x = -0.2, y = y.plot.range[2], bty = "n", bg = "white", legend = leg.txt, lty = lty.vec, lwd = lwd.vec, col = c.vec, cex = 1.0)

      leg.txt <- paste("Neg. ES \"", phen2, "\" ", sep = "", collapse = "")
      text(x = x.plot.range[1], y = -0.1 * max(temp$y), adj = c(0, 0), labels = leg.txt, cex = 1.0)
      leg.txt <- paste(" Pos. ES: \"", phen1, "\" ", sep = "", collapse = "")
      text(x = x.plot.range[2], y = -0.1 * max(temp$y), adj = c(1, 0), labels = leg.txt, cex = 1.0)

      # create pinkogram for each gene set

      kk <- 1

      pinko <- matrix(0, nrow = size.G[i], ncol = cols)
      pinko.gene.names <- vector(length = size.G[i], mode = "character")
      for (k in 1:rows) {
        if (Obs.indicator[i, k] == 1) {
          pinko[kk, ] <- A[obs.index[k], ]
          pinko.gene.names[kk] <- obs.gene.symbols[k]
          kk <- kk + 1
        }
      }
      GSEA.HeatMapPlot(V = pinko, row.names = pinko.gene.names, col.labels = class.labels, col.classes = class.phen, col.names = sample.names, main = " Heat Map for Genes in Gene Set", xlab = " ", ylab = " ")

      if (non.interactive.run == F) {
        if (.Platform$OS.type == "windows") {
          savePlot(filename = gs.filename, type = "jpeg", device = dev.cur())
        } else if (.Platform$OS.type == "unix") {
          dev.off()
        }
      } else {
        dev.off()
      }
    } # if p.vals thres
  } # loop over gene sets


  return(list(report1 = report.phen1, report2 = report.phen2))
} # end of definition of GSEA.analysis

write.cls <- function(class.v, phen, filename) {
  f <- file(filename, "w")
  n <- length(phen)
  l <- length(class.v)
  cat(l, n, "1", "\n", file = f, append = TRUE, sep = " ")
  cat("#", phen, "\n", file = f, append = TRUE, sep = " ")
  class.v <- phen[class.v]
  cat(class.v, "\n", file = f, append = TRUE, sep = " ")
  close(f)
}

write.cls.2 <- function(class.v, phen, filename) {
  f <- file(filename, "w")
  n <- length(phen)
  l <- length(class.v)
  cat(l, n, "1", "\n", file = f, append = TRUE, sep = " ")
  cat("#", unlist(phen), "\n", file = f, append = TRUE, sep = " ")
  if (is.vector(class.v)) {
    class.v <- phen[class.v]
    cat(class.v, "\n", file = f, append = TRUE, sep = " ")
  } else {
    class.list <- matrix(0, nrow = length(class.v[, 1]), ncol = length(class.v[1, ]))
    for (i in 1:length(class.v[, 1])) {
      class.list[i, ] <- unlist(phen[[i]])[class.v[i, ]]
      cat(class.list[i, ], "\n", file = f, append = TRUE, sep = " ")
    }
  }
  close(f)
}

write.gct <- function(gct.data.frame, descs = "", filename) {
  f <- file(filename, "w")
  cat("#1.2", "\n", file = f, append = TRUE, sep = "")
  cat(dim(gct.data.frame)[1], "\t", dim(gct.data.frame)[2], "\n", file = f, append = TRUE, sep = "")
  cat("Name", "\t", file = f, append = TRUE, sep = "")
  cat("Description", file = f, append = TRUE, sep = "")

  names <- names(gct.data.frame)
  cat("\t", names[1], file = f, append = TRUE, sep = "")

  if (length(names) > 1) {
    for (j in 2:length(names)) {
      cat("\t", names[j], file = f, append = TRUE, sep = "")
    }
  }
  cat("\n", file = f, append = TRUE, sep = "\t")

  oldWarn <- options(warn = -1)
  m <- matrix(nrow = dim(gct.data.frame)[1], ncol = dim(gct.data.frame)[2] + 2)
  m[, 1] <- row.names(gct.data.frame)
  if (length(descs) > 1) {
    m[, 2] <- descs
  } else {
    m[, 2] <- row.names(gct.data.frame)
  }
  index <- 3
  for (i in 1:dim(gct.data.frame)[2]) {
    m[, index] <- gct.data.frame[, i]
    index <- index + 1
  }
  write.table(m, file = f, append = TRUE, quote = FALSE, sep = "\t", eol = "\n", col.names = FALSE, row.names = FALSE)
  close(f)
  options(warn = 0)
}

matrix.abs.plot <- function(V, axes = F, log = F, norm = T, transpose = T, matrix.order = T, max.v = 1, min.v = 0, main = " ", sub = " ", xlab = " ", ylab = "  ") {
  rows <- length(V[, 1])
  cols <- length(V[1, ])
  if (log == T) {
    V <- log(V)
  }
  B <- matrix(0, nrow = rows, ncol = cols)
  for (i in 1:rows) {
    for (j in 1:cols) {
      if (matrix.order == T) {
        k <- rows - i + 1
      } else {
        k <- i
      }
      if (norm == T) {
        if ((max.v == 1) && (min.v == 0)) {
          max.val <- max(V)
          min.val <- min(V)
        } else {
          max.val <- max.v
          min.val <- min.v
        }
      }
      B[k, j] <- max.val - V[i, j] + min.val
    }
  }
  if (transpose == T) {
    B <- t(B)
  }
  if (norm == T) {
    image(z = B, zlim = c(min.val, max.val), axes = axes, col = rainbow(100, s = 1.0, v = 0.75, start = 0.0, end = 0.75, gamma = 1.5), main = main, sub = sub, xlab = xlab, ylab = ylab)
  } else {
    image(z = B, axes = axes, col = rainbow(100, s = 1, v = 0.6, start = 0.1, end = 0.9, gamma = 1), main = main, sub = sub, xlab = xlab, ylab = ylab)
  }
  return(list(B, max.val, min.val))
}



MSIG.HeatMapPlot <- function(
    V,
    row.names = "NA",
    col.labels = "NA",
    col.classes = "NA",
    col.names = "NA",
    main = " ",
    sub = " ",
    xlab = " ",
    ylab = " ",
    row.norm = TRUE,
    cmap.type = 1) # 1 = pinkogram, 2 = scale of blues
{
  #
  # Plots a heatmap "pinkogram" of a gene expression matrix including phenotype vector and gene, sample and phenotype labels
  #
  # The Broad Institute
  # SOFTWARE COPYRIGHT NOTICE AGREEMENT
  # This software and its documentation are copyright 2003 by the
  # Broad Institute/Massachusetts Institute of Technology.
  # All rights are reserved.
  #
  # This software is supplied without any warranty or guaranteed support
  # whatsoever. Neither the Broad Institute nor MIT can be responsible for
  # its use, misuse, or functionality.

  n.rows <- length(V[, 1])
  n.cols <- length(V[1, ])

  if (row.norm == TRUE) {
    row.mean <- apply(V, MARGIN = 1, FUN = mean)
    row.sd <- apply(V, MARGIN = 1, FUN = sd)
    row.n <- length(V[, 1])
    for (i in 1:n.rows) {
      if (row.sd[i] == 0) {
        V[i, ] <- 0
      } else {
        V[i, ] <- (V[i, ] - row.mean[i]) / (0.333 * row.sd[i])
      }
      V[i, ] <- ifelse(V[i, ] < -6, -6, V[i, ])
      V[i, ] <- ifelse(V[i, ] > 6, 6, V[i, ])
    }
  }

  if (cmap.type == 1) {
    mycol <- c("#BBBBBB", "#0000FF", "#4040FF", "#7070FF", "#8888FF", "#A9A9FF", "#D5D5FF", "#EEE5EE", "#FFAADA", "#FF9DB0", "#FF7080", "#FF5A5A", "#FF4040", "#FF0D1D", "#333333") # blue-pinkogram colors. The first and last are the colors to indicate the class vector (phenotype). This is the 1998-vintage, pre-gene cluster, original pinkogram color map
    ncolors <- length(mycol) - 2
  } else if (cmap.type >= 2) {
    library("RColorBrewer")
    #          mycol <- c("#FFA500", brewer.pal(9, "Purples"), "#2E8657")
    mycol <- c("#FFA500", "#FCFBFD", "#F4F2F8", "#F8F7FB", "#EFEDF5", "#E1E1EF", "#E8E7F2", "#DADAEB", "#C6C7E1", "#D0D1E6", "#BCBDDC", "#A8A6CF", "#B2B2D6", "#9E9AC8", "#8A87BF", "#9491C4", "#807DBA", "#7260AB", "#796FB3", "#6A51A3", "#5C3596", "#63439D", "#54278F", "#460D83", "#4D1A89", "#3F007D", "#2E8657")
    ncolors <- length(mycol) - 2
  }

  heatm <- matrix(0, nrow = n.rows + 1, ncol = n.cols)
  heatm[1:n.rows, ] <- V[seq(n.rows, 1, -1), ]

  #           heatm[n.rows + 1,] <- ifelse(col.labels %% 2 == 0, 7, -7)
  maxv <- max(V)
  minv <- min(V)
  rangev <- maxv - minv
  heatm[n.rows + 1, ] <- ifelse(col.labels %% 2 == 0, maxv + (rangev / (ncolors - 1)), minv - (rangev / (ncolors - 1)))

  par(mar = c(10, 20, 10, 4))
  image(1:n.cols, 1:(n.rows + 1), t(heatm), col = mycol, axes = FALSE, main = main, sub = sub, xlab = xlab, ylab = ylab)

  if (length(row.names) > 1) {
    numC <- nchar(row.names)
    size.row.char <- 20 / (n.rows + 15)
    size.col.char <- 35 / (n.cols + 15)
    for (i in 1:n.rows) {
      row.names[i] <- substr(row.names[i], 1, 35)
    }
    row.names <- c(row.names[seq(n.rows, 1, -1)], "Class")
    axis(2, at = 1:(n.rows + 1), labels = row.names, adj = 0.5, tick = FALSE, las = 1, cex.axis = size.row.char, font.axis = 2, line = -1)
  }

  if (length(col.names) > 1) {
    axis(1, at = 1:n.cols, labels = col.names, tick = FALSE, las = 3, cex.axis = size.col.char, font.axis = 2, line = -1)
  }

  if (length(col.labels) > 1) {
    tab <- table(col.labels)
    zinit <- cumsum(tab) - tab + 1
    zfinal <- cumsum(tab)
    locations <- ceiling((zinit + zfinal) / 2)
    size.class.char <- 13 / (length(tab) + 5)
    #           axis(3, at=locations, labels=col.classes, tick=FALSE, las = 1, cex.axis=size.class.char, font.axis=2, line=-1)
    axis(3, at = locations, labels = col.classes, tick = FALSE, las = 1, cex.axis = size.class.char, font.axis = 2, line = -1)
  }

  return()
}

MSIG.HeatMapPlot.2 <- function(
    V,
    row.names = "NA",
    col.labels = "NA",
    col.names = "NA",
    col.symbols = "NA",
    main = " ",
    sub = " ",
    xlab = " ",
    ylab = " ",
    row.norm = TRUE,
    cmap.type = 1) # 1 = pinkogram, 2 = scale of blues
{
  #
  # Plots a heatmap "pinkogram" of a gene expression matrix including phenotype vector and gene, sample and phenotype labels
  #
  # The Broad Institute
  # SOFTWARE COPYRIGHT NOTICE AGREEMENT
  # This software and its documentation are copyright 2003 by the
  # Broad Institute/Massachusetts Institute of Technology.
  # All rights are reserved.
  #
  # This software is supplied without any warranty or guaranteed support
  # whatsoever. Neither the Broad Institute nor MIT can be responsible for
  # its use, misuse, or functionality.

  n.rows <- length(V[, 1])
  n.cols <- length(V[1, ])

  if (row.norm == TRUE) {
    row.mean <- apply(V, MARGIN = 1, FUN = mean)
    row.sd <- apply(V, MARGIN = 1, FUN = sd)
    row.n <- length(V[, 1])
    for (i in 1:n.rows) {
      if (row.sd[i] == 0) {
        V[i, ] <- 0
      } else {
        V[i, ] <- (V[i, ] - row.mean[i]) / (0.333 * row.sd[i])
      }
      V[i, ] <- ifelse(V[i, ] < -6, -6, V[i, ])
      V[i, ] <- ifelse(V[i, ] > 6, 6, V[i, ])
    }
  }

  if (cmap.type == 1) {
    mycol <- c("#BBBBBB", "#0000FF", "#4040FF", "#7070FF", "#8888FF", "#A9A9FF", "#D5D5FF", "#EEE5EE", "#FFAADA", "#FF9DB0", "#FF7080", "#FF5A5A", "#FF4040", "#FF0D1D", "#333333") # blue-pinkogram colors. The first and last are the colors to indicate the class vector (phenotype). This is the 1998-vintage, pre-gene cluster, original pinkogram color map
    ncolors <- length(mycol) - 2
  } else if (cmap.type >= 2) {
    library("RColorBrewer")
    #          mycol <- c("#FFA500", brewer.pal(9, "Purples"), "#2E8657")
    mycol <- c("#FFA500", "#FCFBFD", "#F4F2F8", "#F8F7FB", "#EFEDF5", "#E1E1EF", "#E8E7F2", "#DADAEB", "#C6C7E1", "#D0D1E6", "#BCBDDC", "#A8A6CF", "#B2B2D6", "#9E9AC8", "#8A87BF", "#9491C4", "#807DBA", "#7260AB", "#796FB3", "#6A51A3", "#5C3596", "#63439D", "#54278F", "#460D83", "#4D1A89", "#3F007D", "#2E8657")
    ncolors <- length(mycol) - 2
  }

  heatm <- matrix(0, nrow = n.rows + 1, ncol = n.cols)
  heatm[1:n.rows, ] <- V[seq(n.rows, 1, -1), ]

  maxv <- max(V)
  minv <- min(V)
  rangev <- maxv - minv


  col.num <- vector(length = length(col.labels), mode = "numeric")

  if (length(col.labels) > 1) {
    current.num <- 0
    col.num[1] <- current.num
    class.labels <- col.labels[1]
    current.class <- col.labels[1]
    locations <- 1
    for (i in 2:length(col.labels)) {
      if (current.class != col.labels[i]) {
        current.class <- col.labels[i]
        class.labels <- c(class.labels, current.class)
        locations <- c(locations, i)
        current.num <- current.num + 1
      }
      col.num[i] <- current.num
    }
    L <- length(class.labels)
    locations[1:(L - 1)] <- (locations[1:(L - 1)] + (locations[2:L] - 1)) / 2.0
    locations[L] <- (locations[L] + length(col.labels) - 1) / 2.0
  } else {
    col.num <- rep(1, n.cols)
  }

  heatm[n.rows + 1, ] <- ifelse(col.num %% 2 == 0, maxv + (rangev / (ncolors - 1)), minv - (rangev / (ncolors - 1)))

  par(mar = c(10, 14, 10, 4))
  image(1:n.cols, 1:(n.rows + 1), t(heatm), col = mycol, axes = FALSE, main = main, sub = sub, xlab = xlab, ylab = ylab)

  if (length(row.names) > 1) {
    numC <- nchar(row.names)
    size.row.char <- 25 / (n.rows + 15)
    size.col.char <- 25 / (n.cols + 15)
    for (i in 1:n.rows) {
      row.names[i] <- substr(row.names[i], 1, 25)
    }
    row.names <- c(row.names[seq(n.rows, 1, -1)], "Class")
    axis(2, at = 1:(n.rows + 1), labels = row.names, adj = 0.5, tick = FALSE, las = 1, cex.axis = size.row.char, font.axis = 2, line = -1)
  }

  if (length(col.names) > 1) {
    if (length(col.symbols) > 1) {
      col.names <- paste(col.names, col.symbols)
      axis(1, at = 1:n.cols, labels = col.names, tick = FALSE, las = 3, cex.axis = size.col.char, font.axis = 2, line = -1)
    } else {
      axis(1, at = 1:n.cols, labels = col.names, tick = FALSE, las = 3, cex.axis = size.col.char, font.axis = 2, line = -1)
    }
  }

  if (length(col.labels) > 1) {
    size.class.char <- 25 / (L + 15)
    axis(3, at = locations, labels = class.labels, tick = FALSE, las = 1, cex.axis = size.class.char, font.axis = 2, line = -1, padj = 1.1)
  }

  return()
}

erf <- function(x) {
  2 * pnorm(sqrt(2) * x)
}


ReadClsFile <- function(file = "NULL") {
  #
  # Reads a class vector CLS file and defines phenotype and class labels vectors for the samples in a gene expression file (RES or GCT format)
  #
  # The Broad Institute
  # SOFTWARE COPYRIGHT NOTICE AGREEMENT
  # This software and its documentation are copyright 2003 by the
  # Broad Institute/Massachusetts Institute of Technology.
  # All rights are reserved.
  #
  # This software is supplied without any warranty or guaranteed support
  # whatsoever. Neither the Broad Institute nor MIT can be responsible for
  # its use, misuse, or functionality.

  cls.cont <- readLines(file)
  num.lines <- length(cls.cont)
  class.list <- unlist(strsplit(cls.cont[[3]], " "))
  s <- length(class.list)
  t <- table(class.list)
  l <- length(t)
  phen <- vector(length = l, mode = "character")
  class.v <- vector(length = s, mode = "numeric")

  current.label <- class.list[1]
  current.number <- 1
  class.v[1] <- current.number
  phen[1] <- current.label
  phen.count <- 1

  if (length(class.list) > 1) {
    for (i in 2:s) {
      if (class.list[i] == current.label) {
        class.v[i] <- current.number
      } else {
        phen.count <- phen.count + 1
        current.number <- current.number + 1
        current.label <- class.list[i]
        phen[phen.count] <- current.label
        class.v[i] <- current.number
      }
    }
  }
  return(list(phen = phen, class.v = class.v))
}

MSIG.ReadClsFile <- function(file = "NULL") {
  #
  # Reads a class vector CLS file and defines phenotype and class labels vectors (numeric and character) for the samples in a gene expression file (RES or GCT format)
  #
  # The Broad Institute
  # SOFTWARE COPYRIGHT NOTICE AGREEMENT
  # This software and its documentation are copyright 2003 by the
  # Broad Institute/Massachusetts Institute of Technology.
  # All rights are reserved.
  #
  # This software is supplied without any warranty or guaranteed support
  # whatsoever. Neither the Broad Institute nor MIT can be responsible for
  # its use, misuse, or functionality.

  cls.cont <- readLines(file)
  num.lines <- length(cls.cont)
  class.list <- unlist(strsplit(cls.cont[[3]], " "))
  s <- length(class.list)
  t <- table(class.list)
  l <- length(t)
  phen <- vector(length = l, mode = "character")
  class.v <- vector(length = s, mode = "numeric")

  current.label <- class.list[1]
  current.number <- 1
  class.v[1] <- current.number
  phen[1] <- current.label
  phen.count <- 1

  if (length(class.list) > 1) {
    for (i in 2:s) {
      if (class.list[i] == current.label) {
        class.v[i] <- current.number
      } else {
        phen.count <- phen.count + 1
        current.number <- current.number + 1
        current.label <- class.list[i]
        phen[phen.count] <- current.label
        class.v[i] <- current.number
      }
    }
  }
  return(list(phen = phen, class.v = class.v, class.list = class.list))
}

SamplePlot2 <- function(S1 = "NULL", S2 = "NULL", title = " ", phen = "NULL", class.v = "NULL", second.phen = "NULL", second.class.v = "NULL") {
  num.samples <- length(S1)
  min.S1 <- min(S1)
  max.S1 <- max(S1)
  min.S2 <- min(S2)
  max.S2 <- max(S2)
  limit.S1 <- max(S1) + 0.35 * (max(S1) - min(S1))
  limit.S2 <- max(S2) + 0.35 * (max(S2) - min(S2))
  mid.S1 <- min.S1 + 0.5 * (max.S1 - min.S1)
  plot(S1, S2, ylim = c(min.S2, limit.S2), type = "n", main = title)
  offset <- 20
  for (j in 1:num.samples) {
    color.code <- class.v[j] + 2
    pch.code <- second.class.v[j] + offset
    points(S1[j], S2[j], pch = pch.code, type = "p", cex = 2.00, bg = color.code, col = 1)
  }
  leg.txt <- phen
  n.phen <- length(phen)
  p.vec <- rep(22, n.phen)
  c.vec <- 2:(n.phen + 1)
  legend(x = min.S1, y = limit.S2, legend = leg.txt, bty = "n", pch = p.vec, bg = "white", pt.bg = c.vec, col = 1, cex = 2.00)

  leg2.txt <- second.phen
  n.phen2 <- length(second.phen)
  p.vec2 <- 1:n.phen2
  p.vec2 <- p.vec2 + offset
  c.vec2 <- rep(1, n.phen2)
  legend(x = mid.S1, y = limit.S2, legend = leg2.txt, bty = "n", pch = p.vec2, bg = "white", pt.bg = c.vec2, col = c.vec2, cex = 2.00)
}

SamplePlot <- function(S1 = "NULL", S2 = "NULL", title = " ", phen = "NULL", class.v = "NULL") {
  num.samples <- length(S1)
  min.S1 <- min(S1)
  min.S2 <- min(S2)
  max.S2 <- max(S2)
  limit.S1 <- max(S1) + 0.25 * (max(S1) - min(S1))
  limit.S2 <- max(S2) + 0.25 * (max(S2) - min(S2))
  plot(S1, S2, ylim = c(min.S2, limit.S2), type = "n", main = title)
  for (j in 1:num.samples) {
    if (min(class.v) == 0) {
      color.code <- class.v[j] + 1
    } else {
      color.code <- class.v[j]
    }
    points(S1[j], S2[j], pch = 22, type = "p", cex = 1.5, bg = color.code, col = color.code)
  }
  leg.txt <- phen
  n.phen <- length(phen)
  p.vec <- rep(22, n.phen)
  c.vec <- 1:n.phen
  legend(x = min.S1, y = limit.S2, legend = leg.txt, bty = "n", pch = p.vec, bg = "white", pt.bg = c.vec, col = c.vec, cex = 1.00)
}


MSIG.Select.Top.Markers <- function(
    input.ds,
    input.cls,
    output.marker.report,
    output.marker.file,
    output.marker.gene.set.file = "",
    output.marker.plot,
    up.and.down.markers = F,
    topgs = 10,
    seed = 1234,
    non.interactive.run = F) {
  print("Running MSIG.Select.Top.Markers...")

  if (output.marker.gene.set.file != "") {
    gs.file <- file(output.marker.gene.set.file, "w")
  }

  # Feature selection using projected dataset

  if (regexpr(pattern = ".gct", input.ds) == -1) {
    dataset <- GSEA.Res2Frame(filename = input.ds)
    gs.names <- row.names(dataset)
    sample.names <- names(dataset)
    m <- data.matrix(dataset)
  } else {
    dataset <- MSIG.Gct2Frame(filename = input.ds)
    gs.names <- dataset$row.names
    gs.descs <- dataset$descs
    sample.names <- dataset$names
    m <- data.matrix(dataset$ds)
  }

  m1 <- m

  dim(m)
  Ns <- length(m[1, ])
  Ng <- length(m[, 1])

  CLS <- ReadClsFile(file = input.cls)
  class.labels <- CLS$class.v
  class.phen <- CLS$phen
  class.list <- CLS$class.list

  #  Perform one vs all selection of topgs for each class

  if (up.and.down.markers == F) {
    topgs <- ifelse(topgs > floor(Ng / length(class.phen)), floor(Ng / length(class.phen)), topgs)
  } else {
    topgs <- ifelse(topgs > floor(Ng / (2 * length(class.phen))), floor(Ng / (2 * length(class.phen))), topgs)
  }

  if (up.and.down.markers == F) {
    sample.molsig.sorted.subset <- matrix(0, nrow = length(class.phen) * topgs, ncol = Ns)
    sample.molsig.sorted.subset.gs <- vector(length = length(class.phen) * topgs, mode = "character")
    sample.molsig.sorted.s2n <- vector(length = length(class.phen) * topgs, mode = "character")
    sample.molsig.sorted.class <- vector(length = length(class.phen) * topgs, mode = "character")
  } else {
    sample.molsig.sorted.subset <- matrix(0, nrow = length(class.phen) * 2 * topgs, ncol = Ns)
    sample.molsig.sorted.subset.gs <- vector(length = length(class.phen) * 2 * topgs, mode = "character")
    sample.molsig.sorted.s2n <- vector(length = length(class.phen) * 2 * topgs, mode = "character")
    sample.molsig.sorted.class <- vector(length = length(class.phen) * 2 * topgs, mode = "character")
  }
  for (k in 1:length(class.phen)) {
    class.k.labels <- ifelse(class.labels == k, 0, 1)

    col.index <- order(class.k.labels, decreasing = F)
    class.k.labels <- class.k.labels[col.index]
    #      sample.names <- sample.names[col.index]
    for (j in 1:Ng) {
      m1[j, ] <- m[j, col.index]
    }
    names(m1) <- sample.names

    #      print(c("k=", k, " labels=", class.k.labels))

    set.seed(seed)

    O <- GSEA.GeneRanking(m1, class.k.labels, gene.labels, 1, permutation.type = 0, sigma.correction = "GeneCluster", fraction = 1, replace = F, reverse.sign = F)
    order.matrix <- O$order.matrix
    obs.order.matrix <- O$obs.order.matrix
    correl.matrix <- O$s2n.matrix
    obs.correl.matrix <- O$obs.s2n.matrix

    rm(O)

    obs.s2n.orig <- apply(obs.correl.matrix, 1, median) # using median to assign enrichment scores
    obs.index <- order(obs.s2n.orig, decreasing = T)
    obs.s2n <- sort(obs.s2n.orig, decreasing = T)
    sample.molsig.sorted <- m[obs.index, ]
    gs.names.sorted <- gs.names[obs.index]

    if (up.and.down.markers == F) {
      start <- (k - 1) * topgs + 1
      end <- k * topgs
      sample.molsig.sorted.subset[start:end, ] <- sample.molsig.sorted[1:topgs, ]
      sample.molsig.sorted.subset.gs[start:end] <- gs.names.sorted[1:topgs]
      sample.molsig.sorted.s2n[start:end] <- signif(obs.s2n[1:topgs], digits = 3)
      sample.molsig.sorted.class[start:end] <- class.phen[k]
    } else {
      start <- (k - 1) * 2 * topgs + 1

      sample.molsig.sorted.subset[start:(start + topgs - 1), ] <- sample.molsig.sorted[1:topgs, ]
      sample.molsig.sorted.subset[(start + topgs):(start + 2 * topgs - 1), ] <- sample.molsig.sorted[seq(Ng, Ng - topgs + 1, -1), ]
      sample.molsig.sorted.subset.gs[start:(start + topgs - 1)] <- gs.names.sorted[1:topgs]
      sample.molsig.sorted.subset.gs[(start + topgs):(start + 2 * topgs - 1)] <- gs.names.sorted[seq(Ng, Ng - topgs + 1, -1)]

      sample.molsig.sorted.s2n[start:(start + topgs - 1)] <- signif(obs.s2n[1:topgs], digits = 3)
      sample.molsig.sorted.s2n[(start + topgs):(start + 2 * topgs - 1)] <- signif(obs.s2n[seq(Ng, Ng - topgs + 1, -1)], digits = 3)

      sample.molsig.sorted.class[start:(start + 2 * topgs - 1)] <- class.phen[k]
    }

    # saving top markers as a gene set

    if (output.marker.gene.set.file != "") {
      if (up.and.down.markers == F) {
        gene.set <- paste(gs.names.sorted[1:topgs], sep = "\t")
        gene.set.name <- paste("Markers_of_", class.phen[k], sep = "")
        gene.set.desc <- paste("Top markers of phenotype: ", class.phen[k], sep = "")
        cat(gene.set.name, gene.set.desc, gene.set, "\n", file = gs.file, append = TRUE, sep = "\t")
      } else {
        gene.set <- paste(gs.names.sorted[1:topgs], sep = "\t")
        gene.set.name <- paste("UP_Markers_of_", class.phen[k], sep = "")
        gene.set.desc <- paste("Top UP markers of phenotype: ", class.phen[k], sep = "")
        cat(gene.set.name, gene.set.desc, gene.set, "\n", file = gs.file, append = TRUE, sep = "\t")

        gene.set <- paste(gs.names.sorted[seq(Ng, Ng - topgs + 1, -1)], sep = "\t")
        gene.set.name <- paste("DOWN_Markers_of_", class.phen[k], sep = "")
        gene.set.desc <- paste("Top DOWN markers of phenotype: ", class.phen[k], sep = "")
        cat(gene.set.name, gene.set.desc, gene.set, "\n", file = gs.file, append = TRUE, sep = "\t")
      }
    }
  } # for

  if (output.marker.gene.set.file != "") {
    close(gs.file)
  }

  if (non.interactive.run == F) {
    if (.Platform$OS.type == "windows") {
      glob.filename <- output.marker.plot
      x11(height = 22, width = 40)
    } else if (.Platform$OS.type == "unix") {
      glob.filename <- paste(output.marker.plot, ".pdf", sep = "", collapse = "")
      pdf(file = glob.filename, height = 22, width = 40)
    }
  } else {
    if (.Platform$OS.type == "unix") {
      glob.filename <- paste(output.marker.plot, ".pdf", sep = "", collapse = "")
      pdf(file = glob.filename, height = 22, width = 40)
    } else if (.Platform$OS.type == "windows") {
      glob.filename <- paste(output.marker.plot, ".pdf", sep = "", collapse = "")
      pdf(file = glob.filename, height = 22, width = 40)
    }
  }

  nf <- layout(matrix(c(1, 2), 1, 2, byrow = T), widths = c(6, 1), heights = 1, respect = FALSE)

  c1 <- c("grey3", "cadetblue", "darkgreen", "chartreuse2", "red1", "darkred", "orange", "blue2", "lightblue", "pink", "coral")



  MSIG.HeatMapPlot.3(V = sample.molsig.sorted.subset, row.names = sample.molsig.sorted.subset.gs, col.labels = class.labels, col.classes = class.phen, phen.cmap = c1[1:length(class.phen)], col.names = sample.names, main = "Top Markers -- Heat Map", xlab = " ", ylab = " ", sub = " ", row.norm = T, cmap.type = 4)

  # legend

  leg.txt <- class.phen
  n.phen <- length(class.phen)
  p.vec <- rep(22, n.phen)
  c.vec <- c1[1:n.phen]
  par(mar = c(0, 0, 0, 0))
  plot(c(0, 0), c(1, 1), xlim = c(0, 1), ylim = c(0, 1), axes = F, type = "n", xlab = "", ylab = "")
  legend(x = 0, y = 1, legend = leg.txt, bty = "n", xjust = 0, yjust = 1, pch = p.vec, pt.bg = c.vec, col = "black", cex = 1.1, pt.cex = 1.1)

  if (non.interactive.run == F) {
    if (.Platform$OS.type == "windows") {
      savePlot(filename = glob.filename, type = "jpeg", device = dev.cur())
    } else if (.Platform$OS.type == "unix") {
      dev.off()
    }
  } else {
    dev.off()
  }

  print("Saving markers dataset...")

  # saving markers report

  n.markers <- length(sample.molsig.sorted.subset.gs)
  report <- data.frame(cbind(sample.molsig.sorted.subset.gs, sample.molsig.sorted.s2n, sample.molsig.sorted.class))
  names(report) <- c("Name", "S2N", "class")
  row.names(report) <- seq(1, n.markers)
  write.table(report, file = output.marker.report, quote = F, sep = "\t")

  # saving markers gct file

  V <- data.frame(sample.molsig.sorted.subset)
  names(V) <- sample.names
  #   row.names(V) <- paste(seq(1, n.markers), sample.molsig.sorted.subset.gs)
  row.names(V) <- sample.molsig.sorted.subset.gs
  write.gct(gct.data.frame = V, descs = seq(1, n.markers), filename = output.marker.file)
}

NMF.div <- function(V, k, maxniter = 2000, seed = 123456, stopconv = 40, stopfreq = 10) {
  N <- length(V[, 1])
  M <- length(V[1, ])
  set.seed(seed)
  W <- matrix(runif(N * k), nrow = N, ncol = k) # Initialize W and H with random numbers
  H <- matrix(runif(k * M), nrow = k, ncol = M)
  VP <- matrix(nrow = N, ncol = M)
  error.v <- vector(mode = "numeric", length = maxniter)
  new.membership <- vector(mode = "numeric", length = M)
  old.membership <- vector(mode = "numeric", length = M)
  no.change.count <- 0
  eps <- .Machine$double.eps
  for (t in 1:maxniter) {
    VP <- W %*% H
    W.t <- t(W)
    H <- H * (W.t %*% (V / VP)) + eps
    norm <- apply(W, MARGIN = 2, FUN = sum)
    for (i in 1:k) {
      H[i, ] <- H[i, ] / norm[i]
    }
    VP <- W %*% H
    H.t <- t(H)
    W <- W * ((V / VP) %*% H.t) + eps
    norm <- apply(H, MARGIN = 1, FUN = sum)
    for (i in 1:k) {
      W[, i] <- W[, i] / norm[i]
    }
    error.v[t] <- sum(V * log((V + eps) / (VP + eps)) - V + VP) / (M * N)
    if (t %% stopfreq == 0) {
      for (j in 1:M) {
        class <- order(H[, j], decreasing = T)
        new.membership[j] <- class[1]
      }
      if (sum(new.membership == old.membership) == M) {
        no.change.count <- no.change.count + 1
      } else {
        no.change.count <- 0
      }
      if (no.change.count == stopconv) break
      old.membership <- new.membership
    }
  }
  return(list(W = W, H = H, t = t, error.v = error.v))
}

NMF <- function(V, k, maxniter = 2000, seed = 123456, stopconv = 40, stopfreq = 10) {
  N <- length(V[, 1])
  M <- length(V[1, ])
  set.seed(seed)
  W <- matrix(runif(N * k), nrow = N, ncol = k) # Initialize W and H with random numbers
  H <- matrix(runif(k * M), nrow = k, ncol = M)
  VP <- matrix(nrow = N, ncol = M)
  error.v <- vector(mode = "numeric", length = maxniter)
  new.membership <- vector(mode = "numeric", length = M)
  old.membership <- vector(mode = "numeric", length = M)
  eps <- .Machine$double.eps
  for (t in 1:maxniter) {
    VP <- W %*% H
    H <- H * (crossprod(W, V) / crossprod(W, VP)) + eps
    VP <- W %*% H
    H.t <- t(H)
    W <- W * (V %*% H.t) / (VP %*% H.t) + eps
    error.v[t] <- sqrt(sum((V - VP)^2)) / (N * M)
    if (t %% stopfreq == 0) {
      for (j in 1:M) {
        class <- order(H[, j], decreasing = T)
        new.membership[j] <- class[1]
      }
      if (sum(new.membership == old.membership) == M) {
        no.change.count <- no.change.count + 1
      } else {
        no.change.count <- 0
      }
      if (no.change.count == stopconv) break
      old.membership <- new.membership
    }
  }
  return(list(W = W, H = H, t = t, error.v = error.v))
}

# Computes a sparse NMF decomposition
# Adapted from an original C++ program by Yuan Gao

SNMF <- function(V, k, maxniter = 2000, seed = 123456, lambda = 1, stopconv = 40, stopfreq = 10) {
  N <- length(V[, 1])
  M <- length(V[1, ])
  set.seed(seed)
  W <- matrix(runif(N * k), nrow = N, ncol = k) # Initialize W and H with random numbers
  H <- matrix(runif(k * M), nrow = k, ncol = M)
  VP <- matrix(nrow = N, ncol = M)
  error.v <- vector(mode = "numeric", length = maxniter)
  new.membership <- vector(mode = "numeric", length = M)
  old.membership <- vector(mode = "numeric", length = M)
  eps <- .Machine$double.eps
  for (t in 1:maxniter) {
    T <- crossprod(W, W) + lambda * diag(k)
    for (j in 1:M) {
      b <- crossprod(W, V[, j])
      H[, j] <- solve(T, b)
    }
    H[H < 0] <- 0
    H.t <- t(H)
    VH <- V %*% H.t
    HH <- H %*% H.t
    WH <- W %*% HH
    W <- W * VH / WH + eps
    VP <- W %*% H
    error.v[t] <- sqrt(sum((V - VP)^2)) / (N * M)
    if (t %% stopfreq == 0) {
      for (j in 1:M) {
        class <- order(H[, j], decreasing = T)
        new.membership[j] <- class[1]
      }
      if (sum(new.membership == old.membership) == M) {
        no.change.count <- no.change.count + 1
      } else {
        no.change.count <- 0
      }
      if (no.change.count == stopconv) break

      old.membership <- new.membership
    }
  }
  return(list(W = W, H = H, t = t, error.v = error.v))
}

NSNMF.div <- function(V, k, theta, maxniter = 2000, seed = 123456, stopconv = 40, stopfreq = 10) {
  # Non-smooth NMF from Carmona-Saez et al

  N <- length(V[, 1])
  M <- length(V[1, ])
  set.seed(seed)
  W <- matrix(runif(N * k), nrow = N, ncol = k) # Initialize W and H with random numbers
  H <- matrix(runif(k * M), nrow = k, ncol = M)
  VP <- matrix(nrow = N, ncol = M)
  error.v <- vector(mode = "numeric", length = maxniter)
  new.membership <- vector(mode = "numeric", length = M)
  old.membership <- vector(mode = "numeric", length = M)
  no.change.count <- 0
  eps <- .Machine$double.eps

  # smoothing matrix
  SS <- theta / k * (rep(1, k) %*% t(rep(1, k)))
  for (i in 1:k) {
    SS[i, i] <- SS[i, i] + (1 - theta)
  }

  for (t in 1:maxniter) {
    W <- W %*% SS
    VP <- W %*% H
    W.t <- t(W)
    H <- H * (W.t %*% (V / VP)) + eps
    norm <- apply(W, MARGIN = 2, FUN = sum)
    for (i in 1:k) {
      H[i, ] <- H[i, ] / norm[i]
    }
    H <- SS %*% H
    VP <- W %*% H
    H.t <- t(H)
    W <- W * ((V / VP) %*% H.t) + eps
    norm <- apply(H, MARGIN = 1, FUN = sum)
    for (i in 1:k) {
      W[, i] <- W[, i] / norm[i]
    }
    error.v[t] <- sum(V * log((V + eps) / (VP + eps)) - V + VP) / (M * N)
    if (t %% stopfreq == 0) {
      for (j in 1:M) {
        class <- order(H[, j], decreasing = T)
        new.membership[j] <- class[1]
      }
      if (sum(new.membership == old.membership) == M) {
        no.change.count <- no.change.count + 1
      } else {
        no.change.count <- 0
      }
      if (no.change.count == stopconv) break
      old.membership <- new.membership
    }
  }
  return(list(W = W, H = H, t = t, error.v = error.v))
}


metagene.plot <- function(H, main = " ", sub = " ", xlab = "samples ", ylab = "amplitude") {
  k <- length(H[, 1])
  S <- length(H[1, ])
  index <- 1:S
  maxval <- max(H)
  minval <- min(H)
  plot(index, H[1, ], xlim = c(1, S), ylim = c(minval, maxval), main = main, sub = sub, ylab = ylab, xlab = xlab, type = "n")
  for (i in 1:k) {
    lines(index, H[i, ], type = "l", col = i, lwd = 2)
  }
}

MSIG.Preprocess.Dataset <- function(
    input.ds,
    output.ds,
    thres = NULL,
    ceil = NULL,
    shift = NULL,
    fold = NULL,
    delta = NULL,
    normalization = NULL,
    cntrl.genes = NULL) {
  print(c("Running MSIG.Preprocess.Dataset... on:", input.ds))
  print(c("output file:", output.ds))
  print(c("normalization =", normalization))

  # Read dataset

  dataset <- MSIG.Gct2Frame(filename = input.ds)
  m <- data.matrix(dataset$ds)
  gs.names <- dataset$row.names
  gs.descs <- dataset$descs
  sample.names <- dataset$names

  # threshold, ceiling and shift

  if (!is.null(thres)) {
    m[m < thres] <- thres
  }
  if (!is.null(ceil)) {
    m[m > ceil] <- ceil
  }
  if (!is.null(shift)) {
    m <- m + shift
  }

  # identify and save control genes

  if (!is.null(cntrl.genes)) {
    gene.names2 <- intersect(cntrl.genes, gs.names)
    locs <- match(gene.names2, gs.names, nomatch = 0)
    msig.cntrl <- m[locs, ]
    msig.cntrl.genes <- gs.names[locs]
    msig.cntrl.descs <- gs.descs[locs]
    m <- m[-locs, ]
    gs.names <- gs.names[-locs]
    gs.descs <- gs.descs[-locs]
  }

  # variation filter

  if ((!is.null(fold)) && (!is.null(delta))) {
    temp <- MSIG.VarFilter(V = m, fold = fold, delta = delta, gene.names = gs.names, gene.descs = gs.descs)
    m <- temp$V
    gs.names <- temp$new.gene.names
    gs.descs <- temp$new.gene.descs
    dim(m)
  }

  # restore control genes

  if (!is.null(cntrl.genes)) {
    m <- rbind(m, msig.cntrl)
    gs.names <- c(gs.names, msig.cntrl.genes)
    gs.descs <- c(gs.descs, msig.cntrl.descs)
  }

  # normalization

  if (!is.null(normalization)) {
    if (normalization == 1) {
      m <- MSIG.NormalizeCols.Rank(m)
    } else if (normalization == 2) {
      m <- MSIG.NormalizeCols.Rank(m) / length(m[, 1])
    } else if (normalization == 3) {
      m <- GSEA.NormalizeCols(m) + 3
      m <- GSEA.Threshold(m, 0.001, 100000)
    } else if (normalization == 4) {
      m <- MSIG.NormalizeCols.Rank(m) / length(m[, 1])
    } else if (normalization == 5) {
      m <- MSIG.NormalizeCols.Rescale(m)
    } else if (normalization == 6) {
      cols <- length(m[1, ])
      for (j in 1:cols) { # column rank normalization from 0 to N - 1
        m[, j] <- rank(m[, j], ties.method = "average") - 1
      }
      m <- 10000 * m / (length(m[, 1]) - 1)
    } else if (normalization == 7) {
      m <- ((100 * MSIG.NormalizeCols.Rank(m)) %/% length(m[, 1]) + 1)
    } else if (normalization == 8) {
      row.mean <- apply(m, MARGIN = 1, FUN = mean)
      for (i in 1:length(m[, 1])) {
        m[i, ] <- m[i, ] / row.mean[i]
      }
    }
  }

  V <- data.frame(m)
  names(V) <- sample.names
  row.names(V) <- gs.names
  write.gct(gct.data.frame = V, descs = gs.descs, filename = output.ds)
}

MSIG.Extract.Factors <- function(
    input.ds,
    input.cls = "",
    output.W.file,
    output.H.file,
    k.proj = 2,
    alg = "NMF.div", # decomposition algorithm: NMF.div, NMF, SNMF, NSNMF.div or PCA
    niter = 1000,
    seed = 1234,
    theta = 0,
    sort.factors = F) {
  # start of methodology

  print(c("Running MSIG.Extract.Factors... on", input.ds))


  # Read dataset

  dataset <- MSIG.Gct2Frame(filename = input.ds)
  m <- data.matrix(dataset$ds)
  gs.names <- dataset$row.names
  gs.descs <- dataset$descs
  sample.names <- dataset$names

  dim(m)

  Ns <- length(m[1, ])
  Ng <- length(m[, 1])

  if (input.cls != "") {
    CLS <- ReadClsFile(file = input.cls)
    class.labels <- CLS$class.v
    class.phen <- CLS$phen
  } else {
    class.labels <- rep(1, Ns)
    class.phen <- "Samples"
  }

  if (alg == "PCA") { # PCA projection
    svd.proj <- svd(m, nv = Ns, nu = Ns)
    H.full <- diag(x = svd.proj$d, nrow = Ns, ncol = Ns) %*% t(svd.proj$v)
    H <- H.full[1:k.proj, ]
    W <- svd.proj$u[, 1:k.proj]
  } else if (alg == "NMF.div") { # NMF divergence
    NMF.out <- NMF.div(V = m, k = k.proj, maxniter = niter, seed = seed, stopconv = 40, stopfreq = 10)
    H <- NMF.out$H
    W <- NMF.out$W
  } else if (alg == "NMF") { # NMF Euclidean
    NMF.out <- NMF(V = m, k = k.proj, maxniter = niter, seed = seed, stopconv = 40, stopfreq = 10)
    H <- NMF.out$H
    W <- NMF.out$W
  } else if (alg == "SNMF") { # Sparse NMF Euclidean
    NMF.out <- SNMF(V = m, k = k.proj, maxniter = niter, seed = seed, lambda = 1, stopconv = 40, stopfreq = 10)
    H <- NMF.out$H
    W <- NMF.out$W
  } else if (alg == "NSNMF.div") { # non-smooth NMF
    NMF.out <- NSNMF.div(V = m, k = k.proj, theta = theta, maxniter = niter, seed = seed, stopconv = 40, stopfreq = 10)
    H <- NMF.out$H
    W <- NMF.out$W
  } else {
    stop(c("unknown algorithm:", alg))
  }

  # sort W columns and H rows to make similar projections get similar ordering

  if (sort.factors == T) {
    dist.matrix <- dist(t(W))
    HC <- hclust(dist.matrix, method = "complete")
    W <- W[, HC$order]
    H <- H[HC$order, ]
  }

  factor.names <- paste("F", seq(1, k.proj), sep = "")
  factor.descs <- paste("NMF Extracted Factor Number ", seq(1, k.proj), sep = "")

  # save extracted factors datasets W and H

  V <- data.frame(W)
  names(V) <- factor.names
  row.names(V) <- gs.names
  write.gct(gct.data.frame = V, descs = gs.descs, filename = output.W.file)

  V <- data.frame(H)
  names(V) <- sample.names
  row.names(V) <- factor.names
  write.gct(gct.data.frame = V, descs = factor.descs, filename = output.H.file)
}


MSIG.Projection.Plots <- function(
    input.ds,
    input.cls = "",
    output.2D.sammon.file,
    output.2D.sammon.plot,
    output.3D.sammon.file,
    output.3D.sammon.plot,
    output.heatmap.plot,
    output.hclust.plot,
    title = "",
    niter = 1000,
    seed = 1234,
    non.interactive.run = F,
    col = c("grey3", "cadetblue", "darkgreen", "chartreuse2", "red1", "darkred", "orange", "blue2", "lightblue", "pink", "coral")) {
  print("Running MSIG.Projection.Plots...")

  library("scatterplot3d")
  #   library(MASS)

  set.seed(seed = seed, kind = NULL)

  # Read dataset

  dataset <- MSIG.Gct2Frame(filename = input.ds)
  m <- data.matrix(dataset$ds)
  gs.names <- dataset$row.names
  gs.descs <- dataset$descs
  sample.names <- dataset$names

  dim(m)

  Ns <- length(m[1, ])
  k.proj <- length(m[, 1])

  if (input.cls != "") {
    CLS <- ReadClsFile(file = input.cls)
    class.labels <- CLS$class.v
    class.phen <- CLS$phen
  } else {
    class.labels <- rep(1, Ns)
    class.phen <- "Samples"
  }

  # Sammon map plots

  # 2D

  sammon.map <- sammon(d = dist(t(m)), k = 2, niter = niter, trace = TRUE, magic = 0.2, tol = 1e-6)

  S1 <- sammon.map$points[, 1]
  S2 <- sammon.map$points[, 2]

  # Normalize them between 0 and 1

  S1 <- (S1 - min(S1)) / (max(S1) - min(S1))
  S2 <- (S2 - min(S2)) / (max(S2) - min(S2))

  V <- data.frame(cbind(S1, S2))
  names(V) <- c("S1", "S2")
  row.names(V) <- sample.names
  write.gct(gct.data.frame = V, filename = output.2D.sammon.file)

  c0 <- col
  c1 <- colors()[match(c0, colors())]
  color <- c1[class.labels]

  if (non.interactive.run == F) {
    if (.Platform$OS.type == "windows") {
      plot.filename <- output.2D.sammon.plot
      x11(height = 8, width = 14)
    } else if (.Platform$OS.type == "unix") {
      plot.filename <- paste(output.2D.sammon.plot, ".pdf", sep = "", collapse = "")
      pdf(file = plot.filename, height = 8, width = 14)
    }
  } else {
    if (.Platform$OS.type == "unix") {
      plot.filename <- paste(output.2D.sammon.plot, ".pdf", sep = "", collapse = "")
      pdf(file = plot.filename, height = 8, width = 14)
    } else if (.Platform$OS.type == "windows") {
      plot.filename <- paste(output.2D.sammon.plot, ".pdf", sep = "", collapse = "")
      pdf(file = plot.filename, height = 8, width = 14)
    }
  }

  nf <- layout(matrix(c(1, 2), 1, 2, byrow = T), widths = c(5, 1), heights = 1, respect = FALSE)

  # 1st subplot

  num.samples <- length(S1)
  plot(S1, S2, xlim = c(0, 1), ylim = c(0, 1), type = "n", main = paste(title, " -- 2D Sammon Map", sep = ""), sub = input.ds)
  for (j in 1:num.samples) {
    if (min(class.labels) == 0) {
      color.code <- c1[class.labels[j] + 1]
    } else {
      color.code <- c1[class.labels[j]]
    }
    points(S1[j], S2[j], pch = 21, type = "p", cex = 2.5, bg = color.code, col = "black")
  }

  # 2nd subplot: legend

  leg.txt <- class.phen
  n.phen <- length(class.phen)
  p.vec <- rep(21, n.phen)
  c.vec <- c1[1:n.phen]

  par(mar = c(0, 0, 0, 0))
  plot(c(0, 0), c(1, 1), xlim = c(0, 1), ylim = c(0, 1), axes = F, type = "n", xlab = "", ylab = "")
  legend(x = 0, y = 1, legend = leg.txt, bty = "n", xjust = 0, yjust = 1, pch = p.vec, pt.bg = c.vec, col = "black", cex = 1.5, pt.cex = 2)

  if (non.interactive.run == F) {
    if (.Platform$OS.type == "windows") {
      savePlot(filename = plot.filename, type = "jpeg", device = dev.cur())
    } else if (.Platform$OS.type == "unix") {
      dev.off()
    }
  } else {
    dev.off()
  }

  # 3D plots

  sammon.map <- sammon(d = dist(t(m)), k = 3, niter = niter, trace = TRUE, magic = 0.2, tol = 1e-6)

  S1 <- sammon.map$points[, 3]
  S2 <- sammon.map$points[, 2]
  S3 <- sammon.map$points[, 1]

  # Normalize them between 0 and 1

  S1 <- (S1 - min(S1)) / (max(S1) - min(S1))
  S2 <- (S2 - min(S2)) / (max(S2) - min(S2))
  S3 <- (S3 - min(S3)) / (max(S3) - min(S3))

  V <- data.frame(cbind(S1, S2, S3))
  names(V) <- c("S1", "S2", "S3")
  row.names(V) <- sample.names
  write.gct(gct.data.frame = V, filename = output.3D.sammon.file)

  color <- c1[class.labels]

  if (non.interactive.run == F) {
    if (.Platform$OS.type == "windows") {
      plot.filename <- output.3D.sammon.plot
      x11(height = 14, width = 22)
    } else if (.Platform$OS.type == "unix") {
      plot.filename <- paste(output.3D.sammon.plot, ".pdf", sep = "", collapse = "")
      pdf(file = plot.filename, height = 14, width = 22)
    }
  } else {
    if (.Platform$OS.type == "unix") {
      plot.filename <- paste(output.3D.sammon.plot, ".pdf", sep = "", collapse = "")
      pdf(file = plot.filename, height = 14, width = 22)
    } else if (.Platform$OS.type == "windows") {
      plot.filename <- paste(output.2D.sammon.plot, ".pdf", sep = "", collapse = "")
      pdf(file = plot.filename, height = 14, width = 22)
    }
  }

  nf <- layout(matrix(c(1, 2, 3), 1, 3, byrow = T), widths = c(3, 3, 1), heights = 1, respect = FALSE)

  # angle 15 good for marker plot
  # angle 100, 135

  # Subplot # 1

  x <- scatterplot3d(S1, S2, S3, xlab = "F1", ylab = "F2", zlab = "F3", type = "p", color, angle = 45, pch = 20, main = paste(title, " -- 3D Sammon Map", sep = ""), sub = " ", cex.symbols = 1)
  for (j in 1:num.samples) {
    if (min(class.labels) == 0) {
      color.code <- c1[class.labels[j] + 1]
    } else {
      color.code <- c1[class.labels[j]]
    }
    cex <- 3 * (max(S2) - S2[j]) / (max(S2) - min(S2)) + 1
    x$points3d(S1[j], S2[j], S3[j], col = "black", pch = 21, bg = color.code, cex = cex)
  }

  # Subplot # 2 (reverse S1 and S2 axes)

  x <- scatterplot3d(S2, S1, S3, xlab = "F2", ylab = "F1", zlab = "F3", type = "p", color, angle = 45, pch = 20, main = paste(title, " -- 3D Sammon Map", sep = ""), sub = input.ds, cex.symbols = 1)
  for (j in 1:num.samples) {
    if (min(class.labels) == 0) {
      color.code <- c1[class.labels[j] + 1]
    } else {
      color.code <- c1[class.labels[j]]
    }
    cex <- 3 * (max(S1) - S1[j]) / (max(S1) - min(S1)) + 1
    x$points3d(S2[j], S1[j], S3[j], col = "black", pch = 21, bg = color.code, cex = cex)
  }

  # Subplot # 3 legeng

  leg.txt <- class.phen
  n.phen <- length(class.phen)
  p.vec <- rep(21, n.phen)
  c.vec <- c1[1:n.phen]
  par(mar = c(0, 0, 0, 0))
  plot(c(0, 0), c(1, 1), xlim = c(0, 1), ylim = c(0, 1), axes = F, type = "n", xlab = "", ylab = "")
  legend(x = 0, y = 1, legend = leg.txt, bty = "n", xjust = 0, yjust = 1, pch = p.vec, pt.bg = c.vec, col = "black", cex = 2, pt.cex = 2.5)

  if (non.interactive.run == F) {
    if (.Platform$OS.type == "windows") {
      savePlot(filename = plot.filename, type = "jpeg", device = dev.cur())
    } else if (.Platform$OS.type == "unix") {
      dev.off()
    }
  } else {
    dev.off()
  }

  # Heat map plot

  if (non.interactive.run == F) {
    if (.Platform$OS.type == "windows") {
      plot.filename <- output.heatmap.plot
      x11(height = 14, width = 15)
    } else if (.Platform$OS.type == "unix") {
      plot.filename <- paste(output.heatmap.plot, ".pdf", sep = "", collapse = "")
      pdf(file = plot.filename, height = 14, width = 15)
    }
  } else {
    if (.Platform$OS.type == "unix") {
      plot.filename <- paste(output.heatmap.plot, ".pdf", sep = "", collapse = "")
      pdf(file = plot.filename, height = 14, width = 15)
    } else if (.Platform$OS.type == "windows") {
      plot.filename <- paste(output.heatmap.plot, ".pdf", sep = "", collapse = "")
      pdf(file = plot.filename, height = 14, width = 15)
    }
  }

  MSIG.HeatMapPlot(V = m, row.names = gs.names, col.labels = class.labels, col.classes = class.phen, col.names = sample.names, main = paste(title, " -- Heat Map", sep = ""), sub = " ", xlab = " ", ylab = " ")

  if (non.interactive.run == F) {
    if (.Platform$OS.type == "windows") {
      savePlot(filename = plot.filename, type = "jpeg", device = dev.cur())
    } else if (.Platform$OS.type == "unix") {
      dev.off()
    }
  } else {
    dev.off()
  }

  dist.matrix <- dist(t(m))
  HC <- hclust(dist.matrix, method = "complete")

  if (non.interactive.run == F) {
    if (.Platform$OS.type == "windows") {
      plot.filename <- output.hclust.plot
      x11(height = 15, width = 15)
    } else if (.Platform$OS.type == "unix") {
      plot.filename <- paste(output.hclust.plot, ".pdf", sep = "", collapse = "")
      pdf(file = plot.filename, height = 15, width = 15)
    }
  } else {
    if (.Platform$OS.type == "unix") {
      plot.filename <- paste(output.hclust.plot, ".pdf", sep = "", collapse = "")
      pdf(file = plot.filename, height = 15, width = 15)
    } else if (.Platform$OS.type == "windows") {
      plot.filename <- paste(output.hclust.plot, ".pdf", sep = "", collapse = "")
      pdf(file = plot.filename, height = 15, width = 15)
    }
  }

  plot(HC, xlab = "samples", cex = 0.75, labels = class.phen[class.labels], col = "blue", main = paste(title, " -- Hierarchical Clustering", sep = ""))

  if (non.interactive.run == F) {
    if (.Platform$OS.type == "windows") {
      savePlot(filename = plot.filename, type = "jpeg", device = dev.cur())
    } else if (.Platform$OS.type == "unix") {
      dev.off()
    }
  } else {
    dev.off()
  }
}

MSIG.Factors.Project.old <- function(
    input.ds,
    factors.ds,
    postprojnorm = TRUE,
    output.file) {
  library(MASS)

  # start of methodology

  print(c("Running MSIG.Factors.Project... on: ", input.ds))

  # Read input dataset

  dataset <- MSIG.Gct2Frame(filename = input.ds)
  m <- data.matrix(dataset$ds)
  gs.names <- dataset$row.names
  gs.descs <- dataset$descs
  sample.names <- dataset$names

  # Read factors dataset

  dataset <- MSIG.Gct2Frame(filename = factors.ds)
  W <- data.matrix(dataset$ds)
  W.row.names <- dataset$row.names
  W.row.descs <- dataset$descs
  W.names <- dataset$names

  # Match features to first dataset and create matching m2 dataset

  overlap <- intersect(gs.names, W.row.names)

  print(c("Size of Input dataset=", length(gs.names), " genes"))
  print(c("Size of W matrix (rows)=", length(W.row.names), " genes"))
  print(c("Size of overlap=", length(overlap), " genes"))

  locations.m <- match(overlap, gs.names, nomatch = 0)
  m2 <- m[locations.m, ]

  locations.W <- match(overlap, W.row.names, nomatch = 0)
  W2 <- W[locations.W, ]

  # Project input dataset using factors input

  H <- ginv(W2) %*% m2

  max.H <- max(H)
  min.H <- min(H)

  H <- (H - min.H) / (max.H - min.H)

  #   H <- ifelse(H < 0, 0, H)

  # Normalize projected dataset to the unit hypersphere

  if (postprojnorm == TRUE) {
    n.col <- length(H[1, ])
    for (i in 1:n.col) {
      S.2 <- sqrt(sum(H[, i] * H[, i]))
      #        S.2 <- sum(H[,i])
      H[, i] <- H[, i] / S.2
    }
  }

  # Save projected dataset

  V <- data.frame(H)
  names(V) <- sample.names
  row.names(V) <- W.names
  write.gct(gct.data.frame = V, filename = output.file)
}

MSIG.Factors.Project <- function(
    input.ds,
    factors.ds,
    postprojnorm = TRUE,
    output.file,
    method = "pseudo-inverse") { # method: pseudo-inverse, nnls-solver

  library(MASS)

  # start of methodology

  print(c("Running MSIG.Factors.Project... on: ", input.ds))

  # Read input dataset

  dataset <- MSIG.Gct2Frame(filename = input.ds)
  m <- data.matrix(dataset$ds)
  gs.names <- dataset$row.names
  gs.descs <- dataset$descs
  sample.names <- dataset$names

  # Read factors dataset

  dataset <- MSIG.Gct2Frame(filename = factors.ds)
  W <- data.matrix(dataset$ds)
  W.row.names <- dataset$row.names
  W.row.descs <- dataset$descs
  W.names <- dataset$names

  # Match features to first dataset and create matching m2 dataset

  overlap <- intersect(gs.names, W.row.names)

  print(c("Size of Input dataset=", length(gs.names), " genes"))
  print(c("Size of W matrix (rows)=", length(W.row.names), " genes"))
  print(c("Size of overlap=", length(overlap), " genes"))

  locations.m <- match(overlap, gs.names, nomatch = 0)
  m2 <- m[locations.m, ]

  locations.W <- match(overlap, W.row.names, nomatch = 0)
  W2 <- W[locations.W, ]

  if (method == "pseudo-inverse") {
    # Project input dataset using factors input

    H <- ginv(W2) %*% m2

    # three potential ways to deal with negative values created in the approximated projection
    #
    # I:
    #   max.H <- max(H)
    #   min.H <- min(H)
    #   H <- (H - min.H)/(max.H - min.H)
    #
    # II:
    #   H <- ifelse(H < 0, 0, H)
    #
    # III:
    #  n.col <- length(H[1,])
    #  for (i in 1:n.col) {
    #        max.H <- max(H[,i])
    #        min.H <- min(H[,i])
    #        H[,i] <- (H[,i] - min.H)/(max.H - min.H)
    #  }

    print(c("projecting using pseudo-inverse..."))
  } else if (method == "nnls-solver") { # using a non-negative least square solver

    H <- matrix(0, nrow = length(W2[1, ]), ncol = length(m2[1, ]))

    for (i in 1:length(m2[1, ])) {
      H[, i] <- nnls.fit(W2, m2[, i], wsqrt = 1, eps = 0, rank.tol = 1e-07)
    }

    print(c("projecting using NNLS solver..."))
  } else {
    stop("unknown method")
  }

  # Normalize projected dataset to the unit hypersphere

  if (postprojnorm == TRUE) {
    n.col <- length(H[1, ])
    for (i in 1:n.col) {
      S.2 <- sqrt(sum(H[, i] * H[, i]))
      #        S.2 <- sum(H[,i])
      H[, i] <- H[, i] / S.2
    }
  }


  # Save projected dataset

  V <- data.frame(H)
  names(V) <- sample.names
  row.names(V) <- W.names
  write.gct(gct.data.frame = V, filename = output.file)
}

MSIG.Match.and.Merge <- function(
    input1.ds,
    input1.cls = "",
    input2.ds,
    input2.cls = "",
    output.ds,
    output.cls = "") {
  # start of methodology

  print(c("Running MSIG.Match.and.Merge... on: ", input1.ds, " ", input2.ds))

  # Read input datasets

  dataset1 <- MSIG.Gct2Frame(filename = input1.ds)
  m1 <- data.matrix(dataset1$ds)
  gs.names1 <- dataset1$row.names
  gs.descs1 <- dataset1$descs
  sample.names1 <- dataset1$names

  dataset2 <- MSIG.Gct2Frame(filename = input2.ds)
  m2 <- data.matrix(dataset2$ds)
  gs.names2 <- dataset2$row.names
  gs.descs2 <- dataset2$descs
  sample.names2 <- dataset2$names

  # Read CLS files

  if ((input1.cls != "") & (input2.cls != "")) {
    CLS1 <- ReadClsFile(file = input1.cls)
    class.labels1 <- CLS1$class.v
    class.phen1 <- CLS1$phen

    CLS2 <- ReadClsFile(file = input2.cls)
    class.labels2 <- CLS2$class.v
    class.phen2 <- CLS2$phen
  }

  # Match features to first dataset and create matching m2 dataset

  gs.names3 <- intersect(gs.names1, gs.names2)

  locations1 <- match(gs.names3, gs.names1, nomatch = 0)
  m1 <- m1[locations1, ]
  gs.descs1 <- gs.descs1[locations1]

  locations2 <- match(gs.names3, gs.names2, nomatch = 0)
  m2 <- m2[locations2, ]
  gs.descs2 <- gs.descs2[locations2]

  # Merge datasets

  m3 <- cbind(m1, m2)
  sample.names3 <- c(sample.names1, sample.names2)

  if ((input1.cls != "") & (input2.cls != "")) {
    class.labels3 <- c(class.labels1, class.labels2 + length(class.phen1))
    class.phen3 <- c(class.phen1, class.phen2)
  }

  # Save datasets

  V <- data.frame(m3)
  row.names(V) <- gs.names3
  write.gct(gct.data.frame = V, descs = gs.descs1, filename = output.ds)

  if ((input1.cls != "") & (input2.cls != "")) {
    write.cls(class.v = class.labels3, phen = class.phen3, filename = output.cls)
  }
}

MSIG.Reorder.Dataset <- function(
    input.ds,
    input.cls,
    phen.reorder,
    output.ds,
    output.cls) {
  # Read input datasets

  dataset <- MSIG.Gct2Frame(filename = input.ds)
  m <- data.matrix(dataset$ds)
  gs.names <- dataset$row.names
  gs.descs <- dataset$descs
  sample.names <- dataset$names

  # Read CLS files

  CLS <- ReadClsFile(file = input.cls)
  class.labels <- CLS$class.v
  class.phen <- CLS$phen

  # Sort samples according to new phenotype order

  new.phen <- class.phen[phen.reorder]

  new.labels <- vector(length(class.labels), mode = "numeric")
  for (i in 1:length(class.labels)) {
    new.labels[i] <- phen.reorder[class.labels[i]]
  }

  col.index <- order(new.labels, decreasing = F)
  new.labels <- new.labels[col.index]
  sample.names <- sample.names[col.index]

  rows <- length(m[, 1])
  for (j in 1:rows) {
    m[j, ] <- m[j, col.index]
  }
  names(m) <- sample.names

  # Save datasets

  V <- data.frame(m)
  names(V) <- sample.names
  row.names(V) <- gs.names
  write.gct(gct.data.frame = V, descs = gs.descs, filename = output.ds)

  write.cls(class.v = new.labels, phen = new.phen, filename = output.cls)
}

MSIG.Resort.Rows <- function(
    input.ds,
    new.order,
    new.row.labs = NULL,
    new.row.descs = NULL,
    output.ds) {
  # Read input datasets

  dataset <- MSIG.Gct2Frame(filename = input.ds)
  m <- data.matrix(dataset$ds)
  gs.names <- dataset$row.names
  gs.descs <- dataset$descs
  sample.names <- dataset$names

  m <- m[new.order, ]

  if (length(new.row.labs) == 1) {
    gs.names <- gs.names[new.order]
    gs.descs <- gs.descs[new.order]
  } else {
    gs.names <- new.row.labs
    gs.names <- new.row.descs
  }

  # Save dataset

  V <- data.frame(m)
  names(V) <- sample.names
  row.names(V) <- gs.names
  write.gct(gct.data.frame = V, descs = gs.descs, filename = output.ds)
}

MSIG.GSEA.Project <- function(
    input.ds,
    input.cls = "",
    gs.db,
    output.gct.file,
    output.plot.file,
    non.interactive.run = F,
    nperm = 25,
    projection.type = 0,
    weighted.score.type = 1,
    topgs = 10,
    preproc.type = 0,
    gs.size.threshold.min = 5,
    gs.size.threshold.max = 100000,
    reverse.sign = F,
    seed = 3338) {
  print(" *** Running GSEA projection...")

  if (.Platform$OS.type == "windows") {
    #      memory.limit(6000000000)
    #      memory.limit()
  }

  # Read input data matrix

  set.seed(seed = seed, kind = NULL)

  if (is.data.frame(input.ds)) {
    dataset <- input.ds
  } else {
    if (regexpr(pattern = ".gct", input.ds) == -1) {
      dataset <- GSEA.Res2Frame(filename = input.ds)
    } else {
      dataset <- GSEA.Gct2Frame(filename = input.ds)
    }
  }
  gene.labels <- row.names(dataset)
  sample.names <- names(dataset)
  A <- data.matrix(dataset)
  dim(A)
  cols <- length(A[1, ])
  rows <- length(A[, 1])

  #  Preproc.type control the type of pre-processing: threshold, variation filter, normalization

  if (preproc.type == 1) { # Column normalize (Z-score)
    A <- GSEA.Threshold(A, 0, 100000)
    A <- GSEA.NormalizeCols(A)
  } else if (preproc.type == 2) { # Column (rank)
    A <- GSEA.Threshold(A, 0, 100000)
    for (j in 1:cols) { # column rank normalization
      A[, j] <- rank(A[, j])
    }
  } else if (preproc.type == 3) { # No normalization
  }

  if (input.cls == "") {
    class.labels <- rep(0, cols)
    class.phen <- "   "
  } else {
    CLS <- ReadClsFile(file = input.cls)
    class.labels <- CLS$class.v
    class.phen <- CLS$phen
  }

  # Read input gene set database

  if (regexpr(pattern = ".gmt", gs.db[1]) == -1) {
    temp <- gs.db
  } else {
    temp <- readLines(gs.db)
  }

  max.Ng <- length(temp)
  temp.size.G <- vector(length = max.Ng, mode = "numeric")
  for (i in 1:max.Ng) {
    temp.size.G[i] <- length(unlist(strsplit(temp[[i]], "\t"))) - 2
  }

  max.size.G <- max(temp.size.G)
  gs <- matrix(rep("null", max.Ng * max.size.G), nrow = max.Ng, ncol = max.size.G)
  temp.names <- vector(length = max.Ng, mode = "character")
  temp.desc <- vector(length = max.Ng, mode = "character")
  gs.count <- 1
  for (i in 1:max.Ng) {
    gene.set.size <- length(unlist(strsplit(temp[[i]], "\t"))) - 2
    gs.line <- noquote(unlist(strsplit(temp[[i]], "\t")))
    gene.set.name <- gs.line[1]
    gene.set.desc <- gs.line[1]
    gene.set.tags <- vector(length = gene.set.size, mode = "character")
    for (j in 1:gene.set.size) {
      gene.set.tags[j] <- gs.line[j + 2]
    }
    existing.set <- is.element(gene.set.tags, gene.labels)
    set.size <- length(existing.set[existing.set == T])
    if ((set.size < gs.size.threshold.min) || (set.size > gs.size.threshold.max)) next
    temp.size.G[gs.count] <- set.size
    gs[gs.count, ] <- c(gene.set.tags[existing.set], rep("null", max.size.G - temp.size.G[gs.count]))
    temp.names[gs.count] <- gene.set.name
    temp.desc[gs.count] <- gene.set.desc
    gs.count <- gs.count + 1
  }

  Ng <- gs.count - 1
  gs.names <- vector(length = Ng, mode = "character")
  gs.desc <- vector(length = Ng, mode = "character")
  size.G <- vector(length = Ng, mode = "numeric")
  gs.names <- temp.names[1:Ng]
  gs.desc <- temp.desc[1:Ng]
  size.G <- temp.size.G[1:Ng]

  N <- length(A[, 1])
  Ns <- length(A[1, ])

  print(c("Number of genes:", N))
  print(c("Number of Gene Sets:", Ng))
  print(c("Number of samples:", Ns))
  print(c("Original number of Gene Sets:", max.Ng))
  print(c("Maximum gene set size:", max.size.G))

  sample.ES <- matrix(0, nrow = Ng, ncol = Ns)
  sample.NES <- matrix(0, nrow = Ng, ncol = Ns)
  sample.p.val <- matrix(0, nrow = Ng, ncol = Ns)

  correl.vector <- vector(length = N, mode = "numeric")

  order.matrix <- matrix(nrow = N, ncol = Ns)

  all.gs.descs <- vector(length = Ng, mode = "character")

  # Compute ES score for each gene set and sample

  phi <- array(0, c(Ns, Ng, nperm))
  for (sample.index in 1:Ns) {
    order.matrix[, sample.index] <- order(A[, sample.index], decreasing = T)
    gene.list2 <- order.matrix[, sample.index]
    for (i in 1:Ng) {
      print(paste("Computing observed enrichment for gene set:", i, "in sample:", sample.index, sep = " "))
      gene.set <- gs[i, gs[i, ] != "null"]
      gene.set2 <- vector(length = length(gene.set), mode = "numeric")
      gene.set2 <- match(gene.set, gene.labels)

      if (weighted.score.type == 0) {
        correl.vector <- rep(1, N)
      } else if (weighted.score.type > 0) {
        correl.vector <- A[gene.list2, sample.index]
      }

      if (projection.type == 0) { # GSEA Enrichment (KS-like)
        GSEA.results <- GSEA.EnrichmentScore2(gene.list = gene.list2, gene.set = gene.set2, weighted.score.type = weighted.score.type, correl.vector = correl.vector)
      } else if (projection.type == 1) { # GSEA area under RES
        GSEA.results <- GSEA.EnrichmentScore3(gene.list = gene.list2, gene.set = gene.set2, weighted.score.type = weighted.score.type, correl.vector = correl.vector)
      } else if (projection.type == 2) { # Wilcox score
        GSEA.results <- Wilcox.Score(gene.list = gene.list2, gene.set = gene.set2)
      }
      sample.ES[i, sample.index] <- GSEA.results$ES

      for (r in 1:nperm) {
        reshuffled.gene.labels <- sample(1:rows)
        if (weighted.score.type == 0) {
          correl.vector <- rep(1, N)
        } else if (weighted.score.type == 1) {
          correl.vector <- A[reshuffled.gene.labels, sample.index]
        }
        if (projection.type == 0) {
          GSEA.results <- GSEA.EnrichmentScore2(gene.list = reshuffled.gene.labels, gene.set = gene.set2, weighted.score.type = weighted.score.type, correl.vector = correl.vector)
        } else if (projection.type == 1) { # GSEA area under RES
          GSEA.results <- GSEA.EnrichmentScore3(gene.list = reshuffled.gene.labels, gene.set = gene.set2, weighted.score.type = weighted.score.type, correl.vector = correl.vector)
        } else if (projection.type == 2) { # Wilcox score
          GSEA.results <- Wilcox.Score(gene.list = reshuffled.gene.labels, gene.set = gene.set2)
        }
        phi[sample.index, i, r] <- GSEA.results$ES
      }

      #         if (projection.type == 0) {
      if (sample.ES[i, sample.index] >= 0) {
        pos.phi <- phi[sample.index, i, phi[sample.index, i, ] >= 0]
        if (length(pos.phi) == 0) {
          pos.phi <- 0.5
        }
        pos.m <- mean(pos.phi)
        sample.NES[i, sample.index] <- sample.ES[i, sample.index] / pos.m
        s <- sum(pos.phi >= sample.ES[i, sample.index]) / length(pos.phi)
        s <- ifelse(s == 0, 1 / nperm, s)
        sample.p.val[i, sample.index] <- 1 / s
      } else {
        neg.phi <- phi[sample.index, i, phi[sample.index, i, ] < 0]
        if (length(neg.phi) == 0) {
          neg.phi <- 0.5
        }
        neg.m <- mean(neg.phi)
        sample.NES[i, sample.index] <- sample.ES[i, sample.index] / abs(neg.m)
        s <- sum(neg.phi <= sample.ES[i, sample.index]) / length(neg.phi)
        s <- ifelse(s == 0, 1 / nperm, s)
        sample.p.val[i, sample.index] <- -1 / s
        #            }
      }
    }
  }


  # Produce plot

  #   MSIG.HeatMapPlot.5(V = t(sample.ES), row.names = sample.names, col.labels = rep(1, length(gs.names)), col.names = gs.names, main = "Single Sample GSEA", xlab=" ", ylab=" ", sub = gct.file, row.norm = T,  cmap.type = 4, rotated.col.labels = F)

  phen.cmap <- c("grey3", "cadetblue", "darkgreen", "chartreuse2", "red1", "darkred", "orange", "blue2", "lightblue", "pink", "coral")
  MSIG.HeatMapPlot.4(V = sample.ES, row.names = gs.names, col.labels = class.labels, col.classes = class.phen, phen.cmap = phen.cmap, col.names = sample.names, main = "Single Sample GSEA Projection", xlab = " ", ylab = " ", sub = gct.file, row.norm = T, cmap.type = 4, char.rescale = 0.75)

  savePlot(filename = output.plot.file, type = "jpeg", device = dev.cur())

  # save enrichment projection

  print("Saving projected dataset...")

  V <- data.frame(sample.NES)
  names(V) <- sample.names
  row.names(V) <- gs.names

  write.gct(gct.data.frame = V, descs = all.gs.descs, filename = output.gct.file)
}

MSIG.Subset.Dataset <- function(
    input.ds,
    input.cls = NULL,
    column.subset = "ALL", # subset of column numbers or names (or phenotypes)
    column.sel.type = "samples", # "samples" or "phenotype"
    row.subset = "ALL", # subset of row numbers or names
    output.ds,
    output.cls = NULL) {
  # start of methodology

  print(c("Running MSIG.Subset.Dataset... on GCT file:", input.ds))
  print(c("Running MSIG.Subset.Dataset... on CLS file:", input.cls))

  # Read input datasets

  dataset <- MSIG.Gct2Frame(filename = input.ds)
  m <- data.matrix(dataset$ds)
  gs.names <- dataset$row.names
  gs.descs <- dataset$descs
  sample.names <- dataset$names

  # Read CLS file

  if (!is.null(input.cls)) {
    CLS <- MSIG.ReadClsFile(file = input.cls)
    class.labels <- CLS$class.v
    class.phen <- CLS$phen
    class.list <- CLS$class.list
  }


  # Select desired column subset

  if (column.sel.type == "samples") {
    if (column.subset[1] == "ALL") {
      m2 <- m
      sample.names2 <- sample.names
      if (!is.null(input.cls)) {
        class.labels2 <- class.labels
      }
    } else {
      m2 <- m[, column.subset]
      if (is.numeric(column.subset[1])) {
        sample.names2 <- sample.names[column.subset]
        if (!is.null(input.cls)) {
          class.labels2 <- class.labels[column.subset]
        }
      } else {
        #            locations <- match(sample.names, column.subset)
        locations <- match(column.subset, sample.names)
        sample.names2 <- sample.names[locations]
        if (!is.null(input.cls)) {
          class.labels2 <- class.labels[locations]
        }
      }
    }
  } else if (column.sel.type == "phenotype") {
    locations <- !is.na(match(class.list, column.subset))
    sample.names2 <- sample.names[locations]
    m2 <- m[, locations]
    if (!is.null(input.cls)) {
      class.labels2 <- class.labels[locations]
    }
  }

  if (row.subset[1] == "ALL") {
    m3 <- m2
    gs.names2 <- gs.names
    gs.descs2 <- gs.descs
  } else {
    m3 <- m2[row.subset, ]
    if (is.numeric(row.subset[1])) {
      gs.names2 <- gs.names[row.subset]
      gs.descs2 <- gs.descs[row.subset]
    } else {
      locations <- match(row.subset, gs.names)
      gs.names2 <- gs.names[locations]
      gs.descs2 <- gs.descs[locations]
    }
  }

  # Save datasets

  V <- data.frame(m3)
  names(V) <- sample.names2
  row.names(V) <- gs.names2
  write.gct(gct.data.frame = V, descs = gs.descs2, filename = output.ds)

  if (!is.null(input.cls)) {
    write.cls(class.v = class.labels2, phen = class.phen, filename = output.cls)
  }
}

MSIG.Relabel.Phenotype <- function(
    input.cls,
    renaming.list, # renaming list:  list(old.name1 = new.name1, old.name2 = new.name2,...)
    output.cls) {
  CLS <- MSIG.ReadClsFile(file = input.cls)
  class.labels <- CLS$class.v
  class.phen <- CLS$phen
  class.list <- CLS$class.list

  # find matches and replace values

  for (k in 1:length(class.phen)) {
    value <- class.phen[k]
    for (i in 1:length(renaming.list)) {
      if (names(renaming.list[i]) == value) {
        class.phen[k] <- renaming.list[i][[1]]
      }
    }
  }

  # Save dataset

  write.cls(class.v = class.labels, phen = class.phen, filename = output.cls)
}


MSIG.Sample.Dataset <- function(
    input.ds,
    input.cls = "",
    column.subset.fraction = "ALL", # percentage [0, 1] of columns or "ALL"
    row.subset.fraction = "ALL", # percentage [0, 1] of rows or "ALL"
    output.ds,
    output.cls) {
  # start of methodology

  print("Running MSIG.Sample.Dataset...")

  # Read input datasets

  dataset <- MSIG.Gct2Frame(filename = input.ds)
  m <- data.matrix(dataset$ds)
  gs.names <- dataset$row.names
  gs.descs <- dataset$descs
  sample.names <- dataset$names

  # Read CLS file

  if (input.cls != "") {
    CLS <- ReadClsFile(file = input.cls)
    class.labels <- CLS$class.v
    class.phen <- CLS$phen
  }

  # Select desired column subset

  if (column.subset.fraction[1] == "ALL") {
    m2 <- m
    sample.names2 <- sample.names
    if (input.cls != "") {
      class.labels2 <- class.labels
    }
  } else {
    ncol <- length(m[1, ])
    column.subset <- sample(x = seq(1, ncol), size = ceiling(column.subset.fraction * ncol), replace = FALSE)
    m2 <- m[, column.subset]
    sample.names2 <- sample.names[column.subset]
    if (input.cls != "") {
      class.labels2 <- class.labels[column.subset]
    }
  }
  if (row.subset.fraction[1] == "ALL") {
    m3 <- m2
    gs.names2 <- gs.names
    gs.descs2 <- gs.descs
  } else {
    nrow <- length(m[, 1])
    row.subset <- sample(x = seq(1, nrow), size = ceiling(row.subset.fraction * nrow), replace = FALSE)
    m3 <- m2[row.subset, ]
    gs.names2 <- gs.names[row.subset]
    gs.descs2 <- gs.descs[row.subset]
  }

  # Save datasets

  V <- data.frame(m3)
  names(V) <- sample.names2
  row.names(V) <- gs.names2
  write.gct(gct.data.frame = V, descs = gs.descs2, filename = output.ds)

  if (input.cls != "") {
    write.cls(class.v = class.labels2, phen = class.phen, filename = output.cls)
  }
}



ConsPlot <- function(V, col.names, main = " ", sub = " ", xlab = " ", ylab = " ") {
  # Plots a heatmap plot of a consensus matrix

  cols <- length(V[1, ])
  B <- matrix(0, nrow = cols, ncol = cols)
  max.val <- max(V)
  min.val <- min(V)
  for (i in 1:cols) {
    for (j in 1:cols) {
      k <- cols - i + 1
      B[k, j] <- max.val - V[i, j] + min.val
    }
  }



  #     col.map <- c(rainbow(100, s = 1.0, v = 0.75, start = 0.0, end = 0.75, gamma = 1.5), "#BBBBBB", "#333333", "#FFFFFF")
  col.map <- rev(c("#0000FF", "#4040FF", "#7070FF", "#8888FF", "#A9A9FF", "#D5D5FF", "#EEE5EE", "#FFAADA", "#FF9DB0", "#FF7080", "#FF5A5A", "#FF4040", "#FF0D1D"))

  #     max.size <- max(nchar(col.names))
  par(mar = c(5, 20, 20, 1))
  image(1:cols, 1:cols, t(B), col = col.map, axes = FALSE, main = main, sub = sub, xlab = xlab, ylab = ylab)

  for (i in 1:cols) {
    col.names[i] <- substr(col.names[i], 1, 25)
  }
  col.names2 <- rev(col.names)

  size.col.char <- ifelse(cols < 20, 1.25, 1.25 * sqrt(20 / cols))

  axis(2, at = 1:cols, labels = col.names2, adj = 0.5, tick = FALSE, las = 1, cex.axis = size.col.char, font.axis = 1, line = -1)
  axis(3, at = 1:cols, labels = col.names, adj = 1, tick = FALSE, las = 3, cex.axis = size.col.char, font.axis = 1, line = -1)

  return()
}


MSIG.Projection.Plots.2 <- function(
    input.ds,
    input.cls = "",
    model.set,
    output.2D.proj.file,
    output.2D.proj.plot,
    output.3D.proj.file,
    output.3D.1.proj.plot,
    output.3D.2.proj.plot,
    output.3D.3.proj.plot,
    output.heatmap.plot,
    output.hclust.plot,
    use.feature.names = FALSE,
    use.biplot = TRUE,
    title = "",
    seed = 1234,
    non.interactive.run = F,
    heatmap.row.norm = T,
    heatmap.cmap.type = 1,
    col = c("grey3", "cadetblue", "darkgreen", "chartreuse2", "red1", "darkred", "orange", "blue2", "lightblue", "pink", "coral")) {
  print(c("Running MSIG.Projection.Plots... on:", input.ds))

  library("scatterplot3d")
  library(MASS)

  set.seed(seed = seed, kind = NULL)

  # Read dataset

  dataset <- MSIG.Gct2Frame(filename = input.ds)
  m <- data.matrix(dataset$ds)
  gs.names <- dataset$row.names
  gs.descs <- dataset$descs
  sample.names <- dataset$names

  dim(m)

  Ns <- length(m[1, ])
  k.proj <- length(m[, 1])

  if (input.cls != "") {
    CLS <- ReadClsFile(file = input.cls)
    class.labels <- CLS$class.v
    class.phen <- CLS$phen
  } else {
    class.labels <- rep(1, Ns)
    class.phen <- "Samples"
  }

  # Separate data into train and test pieces

  m.train <- m[, model.set]
  num.samples.train <- length(model.set)
  sample.names.train <- sample.names[model.set]
  if (input.cls != "") {
    class.labels.train <- class.labels[model.set]
  }
  m.test <- m[, -model.set]
  sample.names.test <- sample.names[-model.set]
  if (input.cls != "") {
    class.labels.test <- class.labels[-model.set]
  }



  #   sammon.map <-  sammon(d = dist(t(m)), k = 2, niter = niter, trace = TRUE, magic = 0.2, tol = 1e-6)
  #  S1 <- sammon.map$points[, 1]
  #   S2 <- sammon.map$points[, 2]

  #   pca <- princomp(t(m.train), cor=T)
  #   S1 <- pca$scores[,1]
  #   S2 <- pca$scores[,2]
  #   S3 <- pca$scores[,3]
  #   X1 <- pca$loadings[,1]
  #  X2 <- pca$loadings[,2]
  #   X3 <- pca$loadings[,3]

  pca <- prcomp(t(m.train), retx = TRUE, center = TRUE, scale. = TRUE)

  S1 <- pca$x[, 1]
  S2 <- pca$x[, 2]
  S3 <- pca$x[, 3]
  X1 <- pca$rotation[, 1]
  X2 <- pca$rotation[, 2]
  X3 <- pca$rotation[, 3]


  # 2D plots

  max.S <- max(sqrt(S1 * S1 + S2 * S2))
  max.X <- max(sqrt(X1 * X1 + X2 * X2))
  X1 <- max.S * X1 / max.X
  X2 <- max.S * X2 / max.X
  max.A <- max(max.S, max.X)


  c0 <- col
  c1 <- colors()[match(c0, colors())]
  color <- c1[class.labels]

  if (non.interactive.run == F) {
    if (.Platform$OS.type == "windows") {
      plot.filename <- output.2D.proj.plot
      x11(height = 14, width = 22)
    } else if (.Platform$OS.type == "unix") {
      plot.filename <- paste(output.2D.proj.plot, ".pdf", sep = "", collapse = "")
      pdf(file = plot.filename, height = 14, width = 22)
    }
  } else {
    if (.Platform$OS.type == "unix") {
      plot.filename <- paste(output.2D.proj.plot, ".pdf", sep = "", collapse = "")
      pdf(file = plot.filename, height = 14, width = 22)
    } else if (.Platform$OS.type == "windows") {
      plot.filename <- paste(output.2D.proj.plot, ".pdf", sep = "", collapse = "")
      pdf(file = plot.filename, height = 14, width = 22)
    }
  }

  nf <- layout(matrix(c(1, 2, 3), 1, 3, byrow = T), widths = c(3, 3, 1), heights = 1, respect = FALSE)

  # 1st subplot

  plot(S1, S2, xlim = c(-max.A, max.A), ylim = c(-max.A, max.A), type = "n", main = paste(title, " -- Model Samples Biplot", sep = ""), sub = input.ds)

  for (j in 1:num.samples.train) {
    if (min(class.labels) == 0) {
      color.code <- c1[class.labels[j] + 1]
    } else {
      color.code <- c1[class.labels[j]]
    }
    points(S1[j], S2[j], pch = 21, type = "p", cex = 2.5, bg = color.code, col = "black")
  }

  if (use.biplot == TRUE) {
    for (j in 1:k.proj) {
      x.coor <- X1[j] * 0.925
      y.coor <- X2[j] * 0.925
      arrows(0, 0, x.coor, y.coor, lwd = 1, length = 0.15, angle = 20, col = "black")

      if (use.feature.names == FALSE) {
        leg.txt <- paste("F", j, sep = "")
      } else {
        leg.txt <- gs.names[j]
      }
      text(X1[j], X2[j], labels = leg.txt, adj = NULL, pos = NULL, offset = 0.5, vfont = NULL, cex = 2, col = "black")
    }

    ang <- vector(length = k.proj, mode = "numeric")
    for (j in 1:k.proj) {
      ang[j] <- ifelse(atan2(X2[j], X1[j]) > 0, atan2(X2[j], X1[j]), 2 * pi + atan2(X2[j], X1[j]))
    }

    ang.index <- order(ang, decreasing = F)
    ang2 <- ang[ang.index]

    for (j in 1:k.proj) {
      if (j == k.proj) {
        angle.in.between <- (ang2[1] - ang2[j] - 2 * pi) / 2 + ang2[j] - 2 * pi
      } else {
        angle.in.between <- (ang2[j + 1] - ang2[j]) / 2 + ang2[j]
      }
      x <- max.S * cos(angle.in.between)
      y <- max.S * sin(angle.in.between)
      arrows(0, 0, x, y, lwd = 4, length = 0, lty = 3, col = "grey50")
    }
  }

  # 2nd subplot Project non-model (test) data

  test.scores <- predict(pca, t(m.test))

  #   S1 <- c(pca$scores[,1], test.scores[,1])
  #   S2 <- c(pca$scores[,2], test.scores[,2])
  #   S3 <- c(pca$scores[,3], test.scores[,3])

  S1 <- c(pca$x[, 1], test.scores[, 1])
  S2 <- c(pca$x[, 2], test.scores[, 2])
  S3 <- c(pca$x[, 3], test.scores[, 3])

  max.S <- max(sqrt(S1 * S1 + S2 * S2))
  max.X <- max(sqrt(X1 * X1 + X2 * X2))
  X1 <- max.S * X1 / max.X
  X2 <- max.S * X2 / max.X
  num.samples <- length(S1)

  plot(S1, S2, xlim = c(-max.A, max.A), ylim = c(-max.A, max.A), type = "n", main = paste(title, " -- All Samples Biplot", sep = ""), sub = input.ds)

  for (j in 1:num.samples) {
    if (min(class.labels) == 0) {
      color.code <- c1[class.labels[j] + 1]
    } else {
      color.code <- c1[class.labels[j]]
    }
    points(S1[j], S2[j], pch = 21, type = "p", cex = 2.5, bg = color.code, col = "black")
  }

  if (use.biplot == TRUE) {
    for (j in 1:k.proj) {
      x.coor <- X1[j] * 0.925
      y.coor <- X2[j] * 0.925
      arrows(0, 0, x.coor, y.coor, lwd = 1, length = 0.15, angle = 20, col = "black")
      if (use.feature.names == FALSE) {
        leg.txt <- paste("F", j, sep = "")
      } else {
        leg.txt <- gs.names[j]
      }
      text(X1[j], X2[j], labels = leg.txt, adj = NULL, pos = NULL, offset = 0.5, vfont = NULL, cex = 2, col = "black")
    }

    ang <- vector(length = k.proj, mode = "numeric")
    for (j in 1:k.proj) {
      ang[j] <- ifelse(atan2(X2[j], X1[j]) > 0, atan2(X2[j], X1[j]), 2 * pi + atan2(X2[j], X1[j]))
    }

    ang.index <- order(ang, decreasing = F)
    ang2 <- ang[ang.index]

    for (j in 1:k.proj) {
      if (j == k.proj) {
        angle.in.between <- (ang2[1] - ang2[j] - 2 * pi) / 2 + ang2[j] - 2 * pi
      } else {
        angle.in.between <- (ang2[j + 1] - ang2[j]) / 2 + ang2[j]
      }
      x <- max.S * cos(angle.in.between)
      y <- max.S * sin(angle.in.between)
      arrows(0, 0, x, y, lwd = 4, length = 0, lty = 3, col = "grey50")
    }
  }


  #   test.scores2 <- predict(pca, W)
  #   S1.W <- test.scores2[,1]
  #   S2.W <- test.scores2[,2]
  #   scales.S1.W <- (S1.W - min(S1.W))/(max(S1.W) - min(S1.W))
  #   S1.W <- max.S * S1.W/max.S.W
  #   S2.W <- max.S * S2.W/max.S.W
  #
  #   for (j in 1:num.features/100) {
  #      points(S1.W[j], S2.W[j], type = "p", cex = 2, pch = 20, col = "black")
  #   }

  # 3nd subplot: legend

  leg.txt <- class.phen
  n.phen <- length(class.phen)
  p.vec <- rep(21, n.phen)
  c.vec <- c1[1:n.phen]

  par(mar = c(0, 0, 0, 0))
  plot(c(0, 0), c(1, 1), xlim = c(0, 1), ylim = c(0, 1), axes = F, type = "n", xlab = "", ylab = "")
  legend(x = 0, y = 1, legend = leg.txt, bty = "n", xjust = 0, yjust = 1, pch = p.vec, pt.bg = c.vec, col = "black", cex = 2, pt.cex = 3)

  if (non.interactive.run == F) {
    if (.Platform$OS.type == "windows") {
      savePlot(filename = plot.filename, type = "jpeg", device = dev.cur())
    } else if (.Platform$OS.type == "unix") {
      dev.off()
    }
  } else {
    dev.off()
  }

  # save data

  #   V <- data.frame(cbind(S1, S2))
  #   names(V) <- c("S1", "S2")
  #   row.names(V) <- sample.names
  #   write.gct(gct.data.frame = V, filename = output.2D.sammon.file)


  # 3D plots

  #   sammon.map <-  sammon(d = dist(t(m)), k = 3, niter = niter, trace = TRUE, magic = 0.2, tol = 1e-6)
  #   S1 <- sammon.map$points[, 3]
  #   S2 <- sammon.map$points[, 2]
  #   S3 <- sammon.map$points[, 1]

  # Normalize them between 0 and 1

  #   S1 <- (S1 - min(S1))/(max(S1) - min(S1))
  #   S2 <- (S2 - min(S2))/(max(S2) - min(S2))
  #   S3 <- (S3 - min(S3))/(max(S3) - min(S3))

  max.S <- max(sqrt(S1 * S1 + S2 * S2 + S3 * S3))
  max.X <- max(sqrt(X1 * X1 + X2 * X2 + X3 * X3))

  X1 <- max.S * X1 / max.X
  X2 <- max.S * X2 / max.X
  X3 <- max.S * X3 / max.X
  max.A <- max(max.S, max.X)

  color <- c1[class.labels]

  if (non.interactive.run == F) {
    if (.Platform$OS.type == "windows") {
      plot.filename <- output.3D.1.proj.plot
      x11(height = 18, width = 22)
    } else if (.Platform$OS.type == "unix") {
      plot.filename <- paste(output.3D.1.proj.plot, ".pdf", sep = "", collapse = "")
      pdf(file = plot.filename, height = 18, width = 22)
    }
  } else {
    if (.Platform$OS.type == "unix") {
      plot.filename <- paste(output.3D.1.proj.plot, ".pdf", sep = "", collapse = "")
      pdf(file = plot.filename, height = 18, width = 22)
    } else if (.Platform$OS.type == "windows") {
      plot.filename <- paste(output.3D.1.proj.plot, ".pdf", sep = "", collapse = "")
      pdf(file = plot.filename, height = 18, width = 22)
    }
  }

  #   nf <- layout(matrix(c(1,2,3,4), 2, 2, byrow=T), widths = c(1, 1), heights = c(1, 1), respect = FALSE)

  # angle 15 good for marker plot
  # angle 100, 135

  # Subplot # 1


  nf <- layout(matrix(c(1, 2), 1, 2, byrow = T), widths = c(3, 1), heights = 1, respect = FALSE)

  x <- scatterplot3d(S1, S2, S3, xlab = "F1", ylab = "F2", zlab = "F3", type = "n", angle = 45, pch = 20, main = paste(title, " -- 3D Biplot", sep = ""), sub = " ", cex.symbols = 1)
  for (j in 1:num.samples) {
    if (min(class.labels) == 0) {
      color.code <- c1[class.labels[j] + 1]
    } else {
      color.code <- c1[class.labels[j]]
    }
    cex <- 2.5 * (max(S2) - S2[j]) / (max(S2) - min(S2)) + 1.5
    x$points3d(S1[j], S2[j], S3[j], col = "black", pch = 21, bg = color.code, cex = cex)
  }

  if (use.biplot == TRUE) {
    origin.3D <- x$xyz.convert(0, 0, 0)
    for (j in 1:k.proj) {
      x.coor <- X1[j] * 0.925
      y.coor <- X2[j] * 0.925
      z.coor <- X3[j] * 0.925
      end.point.3D <- x$xyz.convert(x.coor, y.coor, z.coor)
      arrows(origin.3D$x, origin.3D$y, end.point.3D$x, end.point.3D$y, lwd = 2, length = 0.15, angle = 20, col = "grey")

      if (use.feature.names == FALSE) {
        leg.txt <- paste("F", j, sep = "")
      } else {
        leg.txt <- gs.names[j]
      }

      x.coor <- X1[j]
      y.coor <- X2[j]
      z.coor <- X3[j]
      end.point.3D <- x$xyz.convert(x.coor, y.coor, z.coor)
      text(end.point.3D$x, end.point.3D$y, labels = leg.txt, adj = NULL, pos = NULL, offset = 0.5, vfont = NULL, cex = 1.25, col = "grey")
    }
  }



  leg.txt <- class.phen
  n.phen <- length(class.phen)
  p.vec <- rep(21, n.phen)
  c.vec <- c1[1:n.phen]
  par(mar = c(0, 0, 0, 0))
  plot(c(0, 0), c(1, 1), xlim = c(0, 1), ylim = c(0, 1), axes = F, type = "n", xlab = "", ylab = "")
  legend(x = 0, y = 1, legend = leg.txt, bty = "n", xjust = 0, yjust = 1, pch = p.vec, pt.bg = c.vec, col = "black", cex = 2, pt.cex = 2.5)

  if (non.interactive.run == F) {
    if (.Platform$OS.type == "windows") {
      savePlot(filename = plot.filename, type = "jpeg", device = dev.cur())
    } else if (.Platform$OS.type == "unix") {
      dev.off()
    }
  } else {
    dev.off()
  }

  # Subplot # 2 (reverse S2 and S3 axes)

  S3 <- -S3
  X3 <- -X3


  if (non.interactive.run == F) {
    if (.Platform$OS.type == "windows") {
      plot.filename <- output.3D.2.proj.plot
      x11(height = 18, width = 22)
    } else if (.Platform$OS.type == "unix") {
      plot.filename <- paste(output.3D.2.proj.plot, ".pdf", sep = "", collapse = "")
      pdf(file = plot.filename, height = 18, width = 22)
    }
  } else {
    if (.Platform$OS.type == "unix") {
      plot.filename <- paste(output.3D.2.proj.plot, ".pdf", sep = "", collapse = "")
      pdf(file = plot.filename, height = 18, width = 22)
    } else if (.Platform$OS.type == "windows") {
      plot.filename <- paste(output.3D.2.proj.plot, ".pdf", sep = "", collapse = "")
      pdf(file = plot.filename, height = 18, width = 22)
    }
  }

  nf <- layout(matrix(c(1, 2), 1, 2, byrow = T), widths = c(3, 1), heights = 1, respect = FALSE)

  x <- scatterplot3d(S1, S3, S2, xlab = "F1", ylab = "F3", zlab = "F2", type = "n", angle = 45, pch = 20, main = paste(title, " -- 3D Biplot", sep = ""), sub = " ", cex.symbols = 1)
  for (j in 1:num.samples) {
    if (min(class.labels) == 0) {
      color.code <- c1[class.labels[j] + 1]
    } else {
      color.code <- c1[class.labels[j]]
    }
    cex <- 2.5 * (max(S3) - S3[j]) / (max(S3) - min(S3)) + 1.5
    x$points3d(S1[j], S3[j], S2[j], col = "black", pch = 21, bg = color.code, cex = cex)
  }

  if (use.biplot == TRUE) {
    origin.3D <- x$xyz.convert(0, 0, 0)
    for (j in 1:k.proj) {
      x.coor <- X1[j] * 0.925
      y.coor <- X3[j] * 0.925
      z.coor <- X2[j] * 0.925
      end.point.3D <- x$xyz.convert(x.coor, y.coor, z.coor)
      arrows(origin.3D$x, origin.3D$y, end.point.3D$x, end.point.3D$y, lwd = 2, length = 0.15, angle = 20, col = "grey")

      if (use.feature.names == FALSE) {
        leg.txt <- paste("F", j, sep = "")
      } else {
        leg.txt <- gs.names[j]
      }

      x.coor <- X1[j]
      y.coor <- X3[j]
      z.coor <- X2[j]
      end.point.3D <- x$xyz.convert(x.coor, y.coor, z.coor)
      text(end.point.3D$x, end.point.3D$y, labels = leg.txt, adj = NULL, pos = NULL, offset = 0.5, vfont = NULL, cex = 1.25, col = "grey")
    }
  }



  leg.txt <- class.phen
  n.phen <- length(class.phen)
  p.vec <- rep(21, n.phen)
  c.vec <- c1[1:n.phen]
  par(mar = c(0, 0, 0, 0))
  plot(c(0, 0), c(1, 1), xlim = c(0, 1), ylim = c(0, 1), axes = F, type = "n", xlab = "", ylab = "")
  legend(x = 0, y = 1, legend = leg.txt, bty = "n", xjust = 0, yjust = 1, pch = p.vec, pt.bg = c.vec, col = "black", cex = 2, pt.cex = 2.5)

  if (non.interactive.run == F) {
    if (.Platform$OS.type == "windows") {
      savePlot(filename = plot.filename, type = "jpeg", device = dev.cur())
    } else if (.Platform$OS.type == "unix") {
      dev.off()
    }
  } else {
    dev.off()
  }

  # Subplot # 3 (reverse S2 and S3 axes)

  S1 <- -S1
  X1 <- -X1


  if (non.interactive.run == F) {
    if (.Platform$OS.type == "windows") {
      plot.filename <- output.3D.3.proj.plot
      x11(height = 18, width = 22)
    } else if (.Platform$OS.type == "unix") {
      plot.filename <- paste(output.3D.3.proj.plot, ".pdf", sep = "", collapse = "")
      pdf(file = plot.filename, height = 18, width = 22)
    }
  } else {
    if (.Platform$OS.type == "unix") {
      plot.filename <- paste(output.3D.3.proj.plot, ".pdf", sep = "", collapse = "")
      pdf(file = plot.filename, height = 18, width = 22)
    } else if (.Platform$OS.type == "windows") {
      plot.filename <- paste(output.3D.3.proj.plot, ".pdf", sep = "", collapse = "")
      pdf(file = plot.filename, height = 18, width = 22)
    }
  }


  nf <- layout(matrix(c(1, 2), 1, 2, byrow = T), widths = c(3, 1), heights = 1, respect = FALSE)

  x <- scatterplot3d(S2, S1, S3, xlab = "F2", ylab = "F1", zlab = "F3", type = "n", angle = 45, pch = 20, main = paste(title, " -- 3D Biplot", sep = ""), sub = " ", cex.symbols = 1)
  for (j in 1:num.samples) {
    if (min(class.labels) == 0) {
      color.code <- c1[class.labels[j] + 1]
    } else {
      color.code <- c1[class.labels[j]]
    }
    cex <- 2.5 * (max(S1) - S1[j]) / (max(S1) - min(S1)) + 1.5
    x$points3d(S2[j], S1[j], S3[j], col = "black", pch = 21, bg = color.code, cex = cex)
  }

  if (use.biplot == TRUE) {
    origin.3D <- x$xyz.convert(0, 0, 0)
    for (j in 1:k.proj) {
      x.coor <- X2[j] * 0.925
      y.coor <- X1[j] * 0.925
      z.coor <- X3[j] * 0.925
      end.point.3D <- x$xyz.convert(x.coor, y.coor, z.coor)
      arrows(origin.3D$x, origin.3D$y, end.point.3D$x, end.point.3D$y, lwd = 2, length = 0.15, angle = 20, col = "grey")

      if (use.feature.names == FALSE) {
        leg.txt <- paste("F", j, sep = "")
      } else {
        leg.txt <- gs.names[j]
      }

      x.coor <- X2[j]
      y.coor <- X1[j]
      z.coor <- X3[j]
      end.point.3D <- x$xyz.convert(x.coor, y.coor, z.coor)
      text(end.point.3D$x, end.point.3D$y, labels = leg.txt, adj = NULL, pos = NULL, offset = 0.5, vfont = NULL, cex = 1.25, col = "grey")
    }
  }

  # Subplot # 4 legeng

  ### for paper

  leg.txt <- class.phen
  n.phen <- length(class.phen)
  p.vec <- rep(21, n.phen)
  c.vec <- c1[1:n.phen]
  par(mar = c(0, 0, 0, 0))
  plot(c(0, 0), c(1, 1), xlim = c(0, 1), ylim = c(0, 1), axes = F, type = "n", xlab = "", ylab = "")
  legend(x = 0, y = 1, legend = leg.txt, bty = "n", xjust = 0, yjust = 1, pch = p.vec, pt.bg = c.vec, col = "black", cex = 2, pt.cex = 2.5)

  if (non.interactive.run == F) {
    if (.Platform$OS.type == "windows") {
      savePlot(filename = plot.filename, type = "jpeg", device = dev.cur())
    } else if (.Platform$OS.type == "unix") {
      dev.off()
    }
  } else {
    dev.off()
  }



  #   V <- data.frame(cbind(S1, S2, S3))
  #   names(V) <- c("S1", "S2", "S3")
  #   row.names(V) <- sample.names
  #   write.gct(gct.data.frame = V, filename = output.3D.sammon.file)


  # Heat map plot


  height <- ifelse(k.proj > 50, 15, 0.20 * k.proj + 4.8)

  if (non.interactive.run == F) {
    if (.Platform$OS.type == "windows") {
      plot.filename <- output.heatmap.plot
      x11(height = height, width = 15)
    } else if (.Platform$OS.type == "unix") {
      plot.filename <- paste(output.heatmap.plot, ".pdf", sep = "", collapse = "")
      pdf(file = plot.filename, height = height, width = 15)
    }
  } else {
    if (.Platform$OS.type == "unix") {
      plot.filename <- paste(output.heatmap.plot, ".pdf", sep = "", collapse = "")
      pdf(file = plot.filename, height = height, width = 15)
    } else if (.Platform$OS.type == "windows") {
      plot.filename <- paste(output.heatmap.plot, ".pdf", sep = "", collapse = "")
      pdf(file = plot.filename, height = height, width = 15)
    }
  }

  MSIG.HeatMapPlot(V = m, row.names = gs.names, col.labels = class.labels, col.classes = class.phen, col.names = sample.names, main = paste(title, " -- Heat Map", sep = ""), sub = " ", xlab = " ", ylab = " ", row.norm = heatmap.row.norm, cmap.type = heatmap.cmap.type)

  if (non.interactive.run == F) {
    if (.Platform$OS.type == "windows") {
      savePlot(filename = plot.filename, type = "jpeg", device = dev.cur())
    } else if (.Platform$OS.type == "unix") {
      dev.off()
    }
  } else {
    dev.off()
  }

  dist.matrix <- dist(t(m))
  HC <- hclust(dist.matrix, method = "complete")

  if (non.interactive.run == F) {
    if (.Platform$OS.type == "windows") {
      plot.filename <- output.hclust.plot
      x11(height = 15, width = 15)
    } else if (.Platform$OS.type == "unix") {
      plot.filename <- paste(output.hclust.plot, ".pdf", sep = "", collapse = "")
      pdf(file = plot.filename, height = 15, width = 15)
    }
  } else {
    if (.Platform$OS.type == "unix") {
      plot.filename <- paste(output.hclust.plot, ".pdf", sep = "", collapse = "")
      pdf(file = plot.filename, height = 15, width = 15)
    } else if (.Platform$OS.type == "windows") {
      plot.filename <- paste(output.hclust.plot, ".pdf", sep = "", collapse = "")
      pdf(file = plot.filename, height = 15, width = 15)
    }
  }

  plot(HC, xlab = "samples", cex = 0.75, labels = class.phen[class.labels], col = "blue", main = paste(title, " -- Hierarchical Clustering", sep = ""))

  if (non.interactive.run == F) {
    if (.Platform$OS.type == "windows") {
      savePlot(filename = plot.filename, type = "jpeg", device = dev.cur())
    } else if (.Platform$OS.type == "unix") {
      dev.off()
    }
  } else {
    dev.off()
  }



  # Compute cramer test for each pair of classes

  #   library(cramer)
  #   n.phen <- length(class.phen)
  #   cramer.p.vals <- matrix(0.5, nrow=n.phen, ncol=n.phen)
  #   for (i in 1:n.phen) {
  #      for (j in 1:n.phen) {
  #         if (i <= j)
  #            next
  #         else {
  #            A <- t(m[,class.labels == i])
  #            B <- t(m[,class.labels == j])
  #         }
  #        ct <- cramer.test(x = A, y = B)
  #        cramer.p.vals[i, j] <- ct$p.value
  #      }
  #   }
  #
  # for (i in 1:n.phen) {
  #    for (j in 1:n.phen) {
  #          if (i < j) {
  #             cramer.p.vals[i, j] <- cramer.p.vals[j, i]
  #          }
  #    }
  #  }
  #
  #
  #  cramer.table <- data.frame(cramer.p.vals)
  #  names(cramer.table) <- class.phen
  #  row.names(cramer.table) <- class.phen
  #  print(cramer.table)
  ##  cm <- log(cramer.p.vals + 0.000001)
  ##  cm <- ifelse(cramer.p.vals < 0.05, cramer.p.vals, 1)
  #   cm <- cramer.p.vals
  #  print(cm)
  #  x11()
  #
  #     col.map <- rev(c("#0000FF", "#4040FF", "#7070FF", "#8888FF", "#A9A9FF", "#D5D5FF", "#EEE5EE", "#FFAADA", "#FF9DB0", "#FF7080", "#FF5A5A", "#FF4040", "#FF0D1D", "#000000"))
  #
  ##  cm <- ifelse(cm == min(cm),  cm, cm)
  #
  #
  #     V <- cm
  #     col.names <- class.phen
  #
  #     cols <- length(V[1,])
  #     B <- matrix(0, nrow=cols, ncol=cols)
  #     max.val <- max(V)
  #     min.val <- min(V)
  #     for (i in 1:cols) {
  #         for (j in 1:cols) {
  #             k <- cols - i + 1
  # 	     B[k, j] <-  max.val - V[i, j] + min.val
  #          }
  #     }
  #
  #
  #     par(mar = c(5, 5, 5, 5))
  #     image(1:cols, 1:cols, t(B), col = col.map, axes=FALSE, main="", sub="", xlab= "", ylab="")
  #
  #     for (i in 1:cols) {
  #         col.names[i]  <- substr(col.names[i], 1, 25)
  #     }
  #     col.names2 <- rev(col.names)
  #
  #     size.col.char <- ifelse(cols < 20, 1.25, 1.25*sqrt(20/cols))
  #
  #     axis(2, at=1:cols, labels=col.names2, adj= 0.5, tick=FALSE, las = 1, cex.axis=size.col.char, font.axis=1, line=-1)
  #     axis(3, at=1:cols, labels=col.names, adj= 1, tick=FALSE, las = 3, cex.axis=size.col.char, font.axis=1, line=-1)
  #
  #
  ## Compute inter and intra class distances and ratios
  #
  #   x11()
  #   n.phen <- length(class.phen)
  #   class.dist <- matrix(0, nrow=n.phen, ncol=n.phen)
  #   for (i in 1:n.phen) {
  #      for (j in 1:n.phen) {
  ##         if (i <= j)
  ##            next
  ##         else {
  #            A <- t(m[,class.labels == i])
  #            B <- t(m[,class.labels == j])
  ##         }
  #         A <- t(m[,class.labels == i])
  #         B <- t(m[,class.labels == j])
  #         mean_A <- apply(A, 2, mean)
  #         mean_B <- apply(B, 2, mean)
  #         D_AB <- as.matrix(dist(rbind(mean_A, mean_B, method = "euclidean")))[2,1]
  #         D_A_inner <- as.matrix(dist(A, method = "euclidean"))
  #         D_B_inner <- as.matrix(dist(B, method = "euclidean"))
  #         mean_D_A_inner <-  mean(D_A_inner[D_A_inner != 0])
  #         mean_D_B_inner <-  mean(D_B_inner[D_B_inner != 0])
  #         class.dist[i, j] <- D_AB / sqrt(mean_D_A_inner**2 + mean_D_B_inner**2)
  ##         class.dist[i, j] <- D_AB / (mean_D_A_inner + mean_D_B_inner)
  #      }
  #   }
  #
  #  class.dist <- max(class.dist) - class.dist
  #
  #  class.dist <- ifelse(class.dist == max(class.dist), (1 + 1/length(col.map)) * class.dist, class.dist)
  #
  #  class.dist.table <- data.frame(class.dist)
  #  names(class.dist.table) <- class.phen
  #  row.names(class.dist.table) <- class.phen
  #  print(class.dist.table)
  ##  cm <- log(class.dist.table + 0.000001)
  #
  #     V <- class.dist.table
  #     col.names <- class.phen
  #
  #     cols <- length(V[1,])
  #     B <- matrix(0, nrow=cols, ncol=cols)
  #     max.val <- max(V)
  #     min.val <- min(V)
  #     for (i in 1:cols) {
  #         for (j in 1:cols) {
  #             k <- cols - i + 1
  # 	     B[k, j] <-  max.val - V[i, j] + min.val
  #          }
  #     }
  #
  #     col.map <- rev(c("#FFFFFF", "#0000FF", "#4040FF", "#7070FF", "#8888FF", "#A9A9FF", "#D5D5FF", "#EEE5EE", "#FFAADA", "#FF9DB0", "#FF7080", "#FF5A5A", "#FF4040", "#FF0D1D", "#000000"))
  #     par(mar = c(5, 5, 5, 5))
  #     image(1:cols, 1:cols, t(B), col = col.map, axes=FALSE, main="", sub="", xlab= "", ylab="")
  ##     image(1:cols, 1:cols, t(B), col = topo.colors(12), axes=FALSE, main="", sub="", xlab= "", ylab="")
  #
  #     for (i in 1:cols) {
  #         col.names[i]  <- substr(col.names[i], 1, 25)
  #     }
  #     col.names2 <- rev(col.names)
  #
  #     size.col.char <- ifelse(cols < 20, 1.25, 1.25*sqrt(20/cols))
  #
  #     axis(2, at=1:cols, labels=col.names2, adj= 0.5, tick=FALSE, las = 1, cex.axis=size.col.char, font.axis=1, line=-1)
  #     axis(3, at=1:cols, labels=col.names, adj= 1, tick=FALSE, las = 3, cex.axis=size.col.char, font.axis=1, line=-1)
}



MSIG.Extract.Features <- function(
    input.ds,
    features,
    output.ds) {
  print(c("Running MSIG.Extract.Features... on dataset:", input.ds))

  # read dataset

  dataset <- MSIG.Gct2Frame(filename = input.ds)
  A <- data.matrix(dataset$ds)
  gs.names <- dataset$row.names
  gs.descs <- dataset$descs
  sample.names <- dataset$names

  dim(A)
  cols <- length(A[1, ])
  rows <- length(A[, 1])

  # read gene set of features

  temp <- readLines(features)

  gene.set.size <- length(unlist(strsplit(temp, "\t"))) - 2
  gs.line <- noquote(unlist(strsplit(temp, "\t")))
  gene.set.name <- gs.line[1]
  gene.set.desc <- gs.line[1]
  gene.set.tags <- vector(length = gene.set.size, mode = "character")
  for (j in 1:gene.set.size) {
    gene.set.tags[j] <- gs.line[j + 2]
  }

  # extract those genes in the gene set

  locations <- match(gene.set.tags, gs.names, nomatch = 0)
  A2 <- A[locations, ]
  A2 <- data.frame(A2)
  names(A2) <- sample.names
  row.names(A2) <- gene.set.tags
  write.gct(A2, filename = output.ds)
}

GSEA.write.gct <- function(gct, filename) {
  f <- file(filename, "w")
  cat("#1.2", "\n", file = f, append = TRUE, sep = "")
  cat(dim(gct)[1], "\t", dim(gct)[2], "\n", file = f, append = TRUE, sep = "")
  cat("Name", "\t", file = f, append = TRUE, sep = "")
  cat("Description", file = f, append = TRUE, sep = "")
  names <- names(gct)
  cat("\t", names[1], file = f, append = TRUE, sep = "")
  for (j in 2:length(names)) {
    cat("\t", names[j], file = f, append = TRUE, sep = "")
  }
  cat("\n", file = f, append = TRUE, sep = "\t")
  oldWarn <- options(warn = -1)
  m <- matrix(nrow = dim(gct)[1], ncol = dim(gct)[2] + 2)
  m[, 1] <- row.names(gct)
  m[, 2] <- row.names(gct)
  index <- 3
  for (i in 1:dim(gct)[2]) {
    m[, index] <- gct[, i]
    index <- index + 1
  }
  write.table(m, file = f, append = TRUE, quote = FALSE, sep = "\t", eol = "\n", col.names = FALSE, row.names = FALSE)
  close(f)
  options(warn = 0)
  return(gct)
}


GSEA.ConsPlot <- function(V, col.names, main = " ", sub = " ", xlab = " ", ylab = " ") {
  # Plots a heatmap plot of a consensus matrix

  cols <- length(V[1, ])
  B <- matrix(0, nrow = cols, ncol = cols)
  max.val <- max(V)
  min.val <- min(V)
  for (i in 1:cols) {
    for (j in 1:cols) {
      k <- cols - i + 1
      B[k, j] <- max.val - V[i, j] + min.val
    }
  }



  #     col.map <- c(rainbow(100, s = 1.0, v = 0.75, start = 0.0, end = 0.75, gamma = 1.5), "#BBBBBB", "#333333", "#FFFFFF")
  col.map <- rev(c("#0000FF", "#4040FF", "#7070FF", "#8888FF", "#A9A9FF", "#D5D5FF", "#EEE5EE", "#FFAADA", "#FF9DB0", "#FF7080", "#FF5A5A", "#FF4040", "#FF0D1D"))

  #     max.size <- max(nchar(col.names))
  par(mar = c(5, 15, 15, 2))
  image(1:cols, 1:cols, t(B), col = col.map, axes = FALSE, main = main, sub = sub, xlab = xlab, ylab = ylab)

  for (i in 1:cols) {
    col.names[i] <- substr(col.names[i], 1, 25)
  }
  col.names2 <- rev(col.names)

  size.col.char <- ifelse(cols < 15, 1, sqrt(15 / cols))

  axis(2, at = 1:cols, labels = col.names2, adj = 0.5, tick = FALSE, las = 1, cex.axis = size.col.char, font.axis = 1, line = -1)
  axis(3, at = 1:cols, labels = col.names, adj = 1, tick = FALSE, las = 3, cex.axis = size.col.char, font.axis = 1, line = -1)

  return()
}

GSEA.HeatMapPlot2 <- function(V, row.names = "NA", col.names = "NA", main = " ", sub = " ", xlab = " ", ylab = " ") {
  #
  # Plots a heatmap of a matrix

  n.rows <- length(V[, 1])
  n.cols <- length(V[1, ])

  #        mycol <- c("#0000FF", "#4040FF", "#7070FF", "#8888FF", "#A9A9FF", "#D5D5FF", "#EEE5EE", "#FFAADA", "#FF9DB0", "#FF7080", "#FF5A5A", "#FF4040", "#FF0D1D") # blue-pinkogram colors. This is the 1998-vintage, pre-gene cluster, original pinkogram color map
  mycol <- rev(rainbow(100, s = 1.0, v = 0.75, start = 0.0, end = 0.75, gamma = 1.5))

  heatm <- matrix(0, nrow = n.rows, ncol = n.cols)
  heatm[1:n.rows, ] <- V[seq(n.rows, 1, -1), ]

  par(mar = c(5, 20, 5, 2))
  image(1:n.cols, 1:n.rows, t(heatm), col = mycol, axes = FALSE, main = main, sub = sub, xlab = xlab, ylab = ylab)

  if (length(row.names) > 1) {
    size.row.char <- ifelse(n.rows < 15, 1, sqrt(15 / n.rows))
    size.col.char <- ifelse(n.cols < 25, 1, sqrt(25 / n.cols))
    for (i in 1:n.rows) {
      row.names[i] <- substr(row.names[i], 1, 40)
    }
    row.names <- row.names[seq(n.rows, 1, -1)]
    axis(2, at = 1:n.rows, labels = row.names, adj = 0.5, tick = FALSE, las = 1, cex.axis = size.row.char, font.axis = 1, line = -1)
  }

  if (length(col.names) > 1) {
    axis(1, at = 1:n.cols, labels = col.names, tick = FALSE, las = 3, cex.axis = size.col.char, font.axis = 2, line = -1)
  }

  return()
}


GSEA.Analyze.Sets <- function(
    directory,
    topgs = "",
    non.interactive.run = F) {
  file.list <- list.files(directory)
  files <- file.list[regexpr(pattern = ".report.", file.list) > 1]
  max.sets <- length(files)

  set.table <- matrix(nrow = max.sets, ncol = 5)

  for (i in 1:max.sets) {
    temp1 <- strsplit(files[i], split = ".report.")
    temp2 <- strsplit(temp1[[1]][1], split = ".")
    s <- length(temp2[[1]])
    prefix.name <- paste(temp2[[1]][1:(s - 1)], sep = "", collapse = "")
    set.name <- temp2[[1]][s]
    temp3 <- strsplit(temp1[[1]][2], split = ".")
    phenotype <- temp3[[1]][1]
    seq.number <- temp3[[1]][2]
    dataset <- paste(temp2[[1]][1:(s - 1)], sep = "", collapse = ".")

    set.table[i, 1] <- files[i]

    set.table[i, 3] <- phenotype
    set.table[i, 4] <- as.numeric(seq.number)
    set.table[i, 5] <- dataset

    set.table[i, 2] <- paste(set.name, dataset, sep = "", collapse = "")
  }

  print(c("set name=", prefix.name))
  doc.string <- prefix.name

  set.table <- noquote(set.table)
  phen.order <- order(set.table[, 3], decreasing = T)
  set.table <- set.table[phen.order, ]
  phen1 <- names(table(set.table[, 3]))[1]
  phen2 <- names(table(set.table[, 3]))[2]
  set.table.phen1 <- set.table[set.table[, 3] == phen1, ]
  set.table.phen2 <- set.table[set.table[, 3] == phen2, ]

  seq.order <- order(as.numeric(set.table.phen1[, 4]), decreasing = F)
  set.table.phen1 <- set.table.phen1[seq.order, ]
  seq.order <- order(as.numeric(set.table.phen2[, 4]), decreasing = F)
  set.table.phen2 <- set.table.phen2[seq.order, ]

  #   max.sets.phen1 <- length(set.table.phen1[,1])
  #   max.sets.phen2 <- length(set.table.phen2[,1])

  if (topgs == "") {
    max.sets.phen1 <- length(set.table.phen1[, 1])
    max.sets.phen2 <- length(set.table.phen2[, 1])
  } else {
    max.sets.phen1 <- ifelse(topgs > length(set.table.phen1[, 1]), length(set.table.phen1[, 1]), topgs)
    max.sets.phen2 <- ifelse(topgs > length(set.table.phen2[, 1]), length(set.table.phen2[, 1]), topgs)
  }

  # Analysis for phen1

  leading.lists <- NULL
  for (i in 1:max.sets.phen1) {
    inputfile <- paste(directory, set.table.phen1[i, 1], sep = "", collapse = "")
    gene.set <- read.table(file = inputfile, sep = "\t", header = T, comment.char = "", as.is = T)
    leading.set <- as.vector(gene.set[gene.set[, "CORE_ENRICHMENT"] == "YES", "SYMBOL"])
    leading.lists <- c(leading.lists, list(leading.set))
    if (i == 1) {
      all.leading.genes <- leading.set
    } else {
      all.leading.genes <- union(all.leading.genes, leading.set)
    }
  }
  max.genes <- length(all.leading.genes)
  M <- matrix(0, nrow = max.sets.phen1, ncol = max.genes)
  for (i in 1:max.sets.phen1) {
    M[i, ] <- sign(match(all.leading.genes, as.vector(leading.lists[[i]]), nomatch = 0)) # notice that the sign is 0 (no tag) or 1 (tag)
  }

  Inter <- matrix(0, nrow = max.sets.phen1, ncol = max.sets.phen1)
  for (i in 1:max.sets.phen1) {
    for (j in 1:max.sets.phen1) {
      Inter[i, j] <- length(intersect(leading.lists[[i]], leading.lists[[j]])) / length(union(leading.lists[[i]], leading.lists[[j]]))
    }
  }

  Itable <- data.frame(Inter)
  names(Itable) <- set.table.phen1[1:max.sets.phen1, 2]
  row.names(Itable) <- set.table.phen1[1:max.sets.phen1, 2]

  if (non.interactive.run == F) {
    if (.Platform$OS.type == "windows") {
      filename <- paste(directory, doc.string, ".leading.overlap.", phen1, sep = "", collapse = "")
      x11(height = 15, width = 15)
    } else if (.Platform$OS.type == "unix") {
      filename <- paste(directory, doc.string, ".leading.overlap.", phen1, ".pdf", sep = "", collapse = "")
      pdf(file = filename, height = 15, width = 15)
    }
  } else {
    if (.Platform$OS.type == "unix") {
      filename <- paste(directory, doc.string, ".leading.overlap.", phen1, ".pdf", sep = "", collapse = "")
      pdf(file = filename, height = 15, width = 15)
    } else if (.Platform$OS.type == "windows") {
      filename <- paste(directory, doc.string, ".leading.overlap.", phen1, ".pdf", sep = "", collapse = "")
      pdf(file = filename, height = 15, width = 15)
    }
  }

  GSEA.ConsPlot(Itable, col.names = set.table.phen1[1:max.sets.phen1, 2], main = " ", sub = paste("Leading Subsets Overlap ", doc.string, " - ", phen1, sep = ""), xlab = " ", ylab = " ")

  if (non.interactive.run == F) {
    if (.Platform$OS.type == "windows") {
      savePlot(filename = filename, type = "jpeg", device = dev.cur())
    } else if (.Platform$OS.type == "unix") {
      dev.off()
    }
  } else {
    dev.off()
  }

  # Save leading subsets in a GCT file

  D.phen1 <- data.frame(M)
  names(D.phen1) <- all.leading.genes
  row.names(D.phen1) <- set.table.phen1[1:max.sets.phen1, 2]
  output <- paste(directory, doc.string, ".", phen1, ".leading.genes.gct", sep = "")
  GSEA.write.gct(D.phen1, filename = output)

  # Save leading subsets as a single gene set in a .gmt file

  row.header <- paste(doc.string, ".", phen1, ".all.leading.genes", sep = "")
  output.line <- paste(all.leading.genes, sep = "\t", collapse = "\t")
  output.line <- paste(row.header, row.header, output.line, sep = "\t", collapse = "")
  output <- paste(directory, doc.string, ".", phen1, ".all.leading.genes.gmt", sep = "")
  write(noquote(output.line), file = output, ncolumns = length(output.line))

  if (non.interactive.run == F) {
    if (.Platform$OS.type == "windows") {
      filename <- paste(directory, doc.string, ".leading.assignment.", phen1, sep = "", collapse = "")
      x11(height = 12, width = 17)
    } else if (.Platform$OS.type == "unix") {
      filename <- paste(directory, doc.string, ".leading.assignment.", phen1, ".pdf", sep = "", collapse = "")
      pdf(file = filename, height = 12, width = 17)
    }
  } else {
    if (.Platform$OS.type == "unix") {
      filename <- paste(directory, doc.string, ".leading.assignment.", phen1, ".pdf", sep = "", collapse = "")
      pdf(file = filename, height = 12, width = 17)
    } else if (.Platform$OS.type == "windows") {
      filename <- paste(directory, doc.string, ".leading.assignment.", phen1, ".pdf", sep = "", collapse = "")
      pdf(file = filename, height = 12, width = 17)
    }
  }

  GSEA.HeatMapPlot2(V = data.matrix(D.phen1), row.names = row.names(D.phen1), col.names = names(D.phen1), main = "Leading Subsets Assignment", sub = paste(doc.string, " - ", phen1, sep = ""), xlab = " ", ylab = " ")

  if (non.interactive.run == F) {
    if (.Platform$OS.type == "windows") {
      savePlot(filename = filename, type = "jpeg", device = dev.cur())
    } else if (.Platform$OS.type == "unix") {
      dev.off()
    }
  } else {
    dev.off()
  }

  DT1.phen1 <- data.matrix(t(D.phen1))
  DT2.phen1 <- data.frame(DT1.phen1)
  names(DT2.phen1) <- set.table.phen1[1:max.sets.phen1, 2]
  row.names(DT2.phen1) <- all.leading.genes
  #   GSEA.write.gct(DT2.phen1, filename=outputfile2.phen1)

  # Analysis for phen2

  leading.lists <- NULL
  for (i in 1:max.sets.phen2) {
    inputfile <- paste(directory, set.table.phen2[i, 1], sep = "", collapse = "")
    gene.set <- read.table(file = inputfile, sep = "\t", header = T, comment.char = "", as.is = T)
    leading.set <- as.vector(gene.set[gene.set[, "CORE_ENRICHMENT"] == "YES", "SYMBOL"])
    leading.lists <- c(leading.lists, list(leading.set))
    if (i == 1) {
      all.leading.genes <- leading.set
    } else {
      all.leading.genes <- union(all.leading.genes, leading.set)
    }
  }
  max.genes <- length(all.leading.genes)
  M <- matrix(0, nrow = max.sets.phen2, ncol = max.genes)
  for (i in 1:max.sets.phen2) {
    M[i, ] <- sign(match(all.leading.genes, as.vector(leading.lists[[i]]), nomatch = 0)) # notice that the sign is 0 (no tag) or 1 (tag)
  }

  Inter <- matrix(0, nrow = max.sets.phen2, ncol = max.sets.phen2)
  for (i in 1:max.sets.phen2) {
    for (j in 1:max.sets.phen2) {
      Inter[i, j] <- length(intersect(leading.lists[[i]], leading.lists[[j]])) / length(union(leading.lists[[i]], leading.lists[[j]]))
    }
  }

  Itable <- data.frame(Inter)
  names(Itable) <- set.table.phen2[1:max.sets.phen2, 2]
  row.names(Itable) <- set.table.phen2[1:max.sets.phen2, 2]

  if (non.interactive.run == F) {
    if (.Platform$OS.type == "windows") {
      filename <- paste(directory, doc.string, ".leading.overlap.", phen2, sep = "", collapse = "")
      x11(height = 15, width = 15)
    } else if (.Platform$OS.type == "unix") {
      filename <- paste(directory, doc.string, ".leading.overlap.", phen2, ".pdf", sep = "", collapse = "")
      pdf(file = filename, height = 15, width = 15)
    }
  } else {
    if (.Platform$OS.type == "unix") {
      filename <- paste(directory, doc.string, ".leading.overlap.", phen2, ".pdf", sep = "", collapse = "")
      pdf(file = filename, height = 15, width = 15)
    } else if (.Platform$OS.type == "windows") {
      filename <- paste(directory, doc.string, ".leading.overlap.", phen2, ".pdf", sep = "", collapse = "")
      pdf(file = filename, height = 15, width = 15)
    }
  }

  GSEA.ConsPlot(Itable, col.names = set.table.phen2[1:max.sets.phen2, 2], main = " ", sub = paste("Leading Subsets Overlap ", doc.string, " - ", phen2, sep = ""), xlab = " ", ylab = " ")

  if (non.interactive.run == F) {
    if (.Platform$OS.type == "windows") {
      savePlot(filename = filename, type = "jpeg", device = dev.cur())
    } else if (.Platform$OS.type == "unix") {
      dev.off()
    }
  } else {
    dev.off()
  }

  # Save leading subsets in a GCT file

  D.phen2 <- data.frame(M)
  names(D.phen2) <- all.leading.genes
  row.names(D.phen2) <- set.table.phen2[1:max.sets.phen2, 2]
  output <- paste(directory, doc.string, ".", phen2, ".leading.genes.gct", sep = "")
  GSEA.write.gct(D.phen2, filename = output)

  # Save primary subsets as a single gene set in a .gmt file

  row.header <- paste(doc.string, ".", phen2, ".all.leading.genes", sep = "")
  output.line <- paste(all.leading.genes, sep = "\t", collapse = "\t")
  output.line <- paste(row.header, row.header, output.line, sep = "\t", collapse = "")
  output <- paste(directory, doc.string, ".", phen2, ".all.leading.genes.gmt", sep = "")
  write(noquote(output.line), file = output, ncolumns = length(output.line))

  if (non.interactive.run == F) {
    if (.Platform$OS.type == "windows") {
      filename <- paste(directory, doc.string, ".leading.assignment.", phen2, sep = "", collapse = "")
      x11(height = 12, width = 17)
    } else if (.Platform$OS.type == "unix") {
      filename <- paste(directory, doc.string, ".leading.assignment.", phen2, ".pdf", sep = "", collapse = "")
      pdf(file = filename, height = 12, width = 17)
    }
  } else {
    if (.Platform$OS.type == "unix") {
      filename <- paste(directory, doc.string, ".leading.assignment.", phen2, ".pdf", sep = "", collapse = "")
      pdf(file = filename, height = 12, width = 17)
    } else if (.Platform$OS.type == "windows") {
      filename <- paste(directory, doc.string, ".leading.assignment.", phen2, ".pdf", sep = "", collapse = "")
      pdf(file = filename, height = 12, width = 17)
    }
  }

  GSEA.HeatMapPlot2(V = data.matrix(D.phen2), row.names = row.names(D.phen2), col.names = names(D.phen2), main = "Leading Subsets Assignment", sub = paste(doc.string, " - ", phen2, sep = ""), xlab = " ", ylab = " ")

  if (non.interactive.run == F) {
    if (.Platform$OS.type == "windows") {
      savePlot(filename = filename, type = "jpeg", device = dev.cur())
    } else if (.Platform$OS.type == "unix") {
      dev.off()
    }
  } else {
    dev.off()
  }

  DT1.phen2 <- data.matrix(t(D.phen2))
  DT2.phen2 <- data.frame(DT1.phen2)
  names(DT2.phen2) <- set.table.phen2[1:max.sets.phen2, 2]
  row.names(DT2.phen2) <- all.leading.genes
  #   GSEA.write.gct(DT2.phen2, filename=outputfile2.phen2)

  # Resort columns and rows for phen1

  A <- data.matrix(D.phen1)
  A.row.names <- row.names(D.phen1)
  A.names <- names(D.phen1)

  # Max.genes

  init <- 1
  for (k in 1:max.sets.phen1) {
    end <- which.max(cumsum(A[k, ]))
    if (end - init > 1) {
      B <- A[, init:end]
      B.names <- A.names[init:end]
      dist.matrix <- dist(t(B))
      HC <- hclust(dist.matrix, method = "average")
      B <- B[, HC$order] + 0.2 * (k %% 2)
      A[, init:end] <- B
      A.names[init:end] <- B.names[HC$order]
      init <- end + 1
    }
  }

  #   x11(width=14, height=10)
  #   GSEA.HeatMapPlot2(V = A, row.names = A.row.names, col.names = A.names, sub = "  ", main = paste("Primary Sets Assignment - ", doc.string, " - ", phen1, sep=""), xlab=" ", ylab=" ")

  dist.matrix <- dist(A)
  HC <- hclust(dist.matrix, method = "average")
  A <- A[HC$order, ]
  A.row.names <- A.row.names[HC$order]

  if (non.interactive.run == F) {
    if (.Platform$OS.type == "windows") {
      filename <- paste(directory, doc.string, ".leading.assignment.clustered.", phen1, sep = "", collapse = "")
      x11(height = 12, width = 17)
    } else if (.Platform$OS.type == "unix") {
      filename <- paste(directory, doc.string, ".leading.assignment.clustered.", phen1, ".pdf", sep = "", collapse = "")
      pdf(file = filename, height = 12, width = 17)
    }
  } else {
    if (.Platform$OS.type == "unix") {
      filename <- paste(directory, doc.string, ".leading.assignment.clustered.", phen1, ".pdf", sep = "", collapse = "")
      pdf(file = filename, height = 12, width = 17)
    } else if (.Platform$OS.type == "windows") {
      filename <- paste(directory, doc.string, ".leading.assignment.clustered.", phen1, ".pdf", sep = "", collapse = "")
      pdf(file = filename, height = 12, width = 17)
    }
  }

  GSEA.HeatMapPlot2(V = A, row.names = A.row.names, col.names = A.names, main = "Leading Subsets Assignment (clustered)", sub = paste(doc.string, " - ", phen1, sep = ""), xlab = " ", ylab = " ")


  if (non.interactive.run == F) {
    if (.Platform$OS.type == "windows") {
      savePlot(filename = filename, type = "jpeg", device = dev.cur())
    } else if (.Platform$OS.type == "unix") {
      dev.off()
    }
  } else {
    dev.off()
  }

  # resort columns and rows for phen2

  A <- data.matrix(D.phen2)
  A.row.names <- row.names(D.phen2)
  A.names <- names(D.phen2)

  # Max.genes

  init <- 1
  for (k in 1:max.sets.phen2) {
    end <- which.max(cumsum(A[k, ]))
    if (end - init > 1) {
      B <- A[, init:end]
      B.names <- A.names[init:end]
      dist.matrix <- dist(t(B))
      HC <- hclust(dist.matrix, method = "average")
      B <- B[, HC$order] + 0.2 * (k %% 2)
      A[, init:end] <- B
      A.names[init:end] <- B.names[HC$order]
      init <- end + 1
    }
  }

  #  x11(width=14, height=10)
  #  GESA.HeatMapPlot2(V = A, row.names = A.row.names, col.names = A.names, sub = "  ", main = paste("Primary Sets Assignment - ", doc.string, " - ", phen2, sep=""), xlab=" ", ylab=" ")

  dist.matrix <- dist(A)
  HC <- hclust(dist.matrix, method = "average")
  A <- A[HC$order, ]
  A.row.names <- A.row.names[HC$order]

  if (non.interactive.run == F) {
    if (.Platform$OS.type == "windows") {
      filename <- paste(directory, doc.string, ".leading.assignment.clustered.", phen2, sep = "", collapse = "")
      x11(height = 12, width = 17)
    } else if (.Platform$OS.type == "unix") {
      filename <- paste(directory, doc.string, ".leading.assignment.clustered.", phen2, ".pdf", sep = "", collapse = "")
      pdf(file = filename, height = 12, width = 17)
    }
  } else {
    if (.Platform$OS.type == "unix") {
      filename <- paste(directory, doc.string, ".leading.assignment.clustered.", phen2, ".pdf", sep = "", collapse = "")
      pdf(file = filename, height = 12, width = 17)
    } else if (.Platform$OS.type == "windows") {
      filename <- paste(directory, doc.string, ".leading.assignment.clustered.", phen2, ".pdf", sep = "", collapse = "")
      pdf(file = filename, height = 12, width = 17)
    }
  }

  GSEA.HeatMapPlot2(V = A, row.names = A.row.names, col.names = A.names, main = "Leading Subsets Assignment (clustered)", sub = paste(doc.string, " - ", phen2, sep = ""), xlab = " ", ylab = " ")

  if (non.interactive.run == F) {
    if (.Platform$OS.type == "windows") {
      savePlot(filename = filename, type = "jpeg", device = dev.cur())
    } else if (.Platform$OS.type == "unix") {
      dev.off()
    }
  } else {
    dev.off()
  }
}

MSIG.Replace.Missing.Values <- function(
    input.ds,
    output.ds,
    miss.val.replacement) {
  # Read dataset

  dataset <- MSIG.Gct2Frame(filename = input.ds)
  data <- dataset$ds
  gs.names <- dataset$row.names
  gs.descs <- dataset$descs
  sample.names <- dataset$names

  # replace missing values

  N <- length(data[, 1])
  M <- length(data[1, ])

  for (i in 1:N) {
    for (j in 1:M) {
      if (is.na(data[i, j])) {
        data[i, j] <- miss.val.replacement
      }
    }
  }

  m <- data.matrix(data)

  # save file

  V <- data.frame(m)
  names(V) <- sample.names
  row.names(V) <- gs.names
  write.gct(gct.data.frame = V, descs = gs.descs, filename = output.ds)
}



MSIG.2D.Plot <- function(
    input.ds,
    input.cls = "",
    output.2D.proj.plot,
    output.heatmap.plot,
    title = "",
    non.interactive.run = F,
    heatmap.row.norm = F,
    heatmap.cmap.type = 1,
    col = c("grey3", "cadetblue", "darkgreen", "chartreuse2", "red1", "darkred", "orange", "blue2", "lightblue", "pink", "coral")) {
  print("Running MSIG.2D.Plot...")

  print(heatmap.row.norm)

  # Read dataset

  dataset <- MSIG.Gct2Frame(filename = input.ds)
  m <- data.matrix(dataset$ds)
  gs.names <- dataset$row.names
  gs.descs <- dataset$descs
  sample.names <- dataset$names

  dim(m)

  Ns <- length(m[1, ])
  k.proj <- length(m[, 1])

  if (input.cls != "") {
    CLS <- ReadClsFile(file = input.cls)
    class.labels <- CLS$class.v
    class.phen <- CLS$phen
  } else {
    class.labels <- rep(1, Ns)
    class.phen <- "Samples"
  }

  # 2D plots

  S1 <- as.real(m[1, ])
  S2 <- as.real(m[2, ])

  range.S1 <- range(S1)
  range.S2 <- range(S2)

  c0 <- col
  c1 <- colors()[match(c0, colors())]
  color <- c1[class.labels]

  if (non.interactive.run == F) {
    if (.Platform$OS.type == "windows") {
      plot.filename <- output.2D.proj.plot
      x11(height = 14, width = 22)
    } else if (.Platform$OS.type == "unix") {
      plot.filename <- paste(output.2D.proj.plot, ".pdf", sep = "", collapse = "")
      pdf(file = plot.filename, height = 14, width = 22)
    }
  } else {
    if (.Platform$OS.type == "unix") {
      plot.filename <- paste(output.2D.proj.plot, ".pdf", sep = "", collapse = "")
      pdf(file = plot.filename, height = 14, width = 22)
    } else if (.Platform$OS.type == "windows") {
      plot.filename <- paste(output.2D.proj.plot, ".pdf", sep = "", collapse = "")
      pdf(file = plot.filename, height = 14, width = 22)
    }
  }

  nf <- layout(matrix(c(1, 2), 1, 2, byrow = T), widths = c(3, 2), heights = 1, respect = FALSE)

  # 2 D Plot

  plot(S1, S2, xlim = range.S1, ylim = range.S2, type = "n", main = paste(title, " -- 2D Plot", sep = ""), sub = input.ds)

  for (j in 1:Ns) {
    if (min(class.labels) == 0) {
      color.code <- c1[class.labels[j] + 1]
    } else {
      color.code <- c1[class.labels[j]]
    }
    points(S1[j], S2[j], pch = 22, type = "p", cex = 2, bg = color.code, col = "black")
  }

  # legend

  leg.txt <- class.phen
  n.phen <- length(class.phen)
  p.vec <- rep(21, n.phen)
  c.vec <- c1[1:n.phen]

  par(mar = c(0, 0, 0, 0))
  plot(c(0, 0), c(1, 1), xlim = c(0, 1), ylim = c(0, 1), axes = F, type = "n", xlab = "", ylab = "")
  legend(x = 0, y = 1, legend = leg.txt, bty = "n", xjust = 0, yjust = 1, pch = p.vec, pt.bg = c.vec, col = "black", cex = 2, pt.cex = 3)

  if (non.interactive.run == F) {
    if (.Platform$OS.type == "windows") {
      savePlot(filename = plot.filename, type = "jpeg", device = dev.cur())
    } else if (.Platform$OS.type == "unix") {
      dev.off()
    }
  } else {
    dev.off()
  }


  # Heat map plot


  height <- ifelse(k.proj > 50, 15, 0.20 * k.proj + 4.8)

  if (non.interactive.run == F) {
    if (.Platform$OS.type == "windows") {
      plot.filename <- output.heatmap.plot
      x11(height = height, width = 15)
    } else if (.Platform$OS.type == "unix") {
      plot.filename <- paste(output.heatmap.plot, ".pdf", sep = "", collapse = "")
      pdf(file = plot.filename, height = height, width = 15)
    }
  } else {
    if (.Platform$OS.type == "unix") {
      plot.filename <- paste(output.heatmap.plot, ".pdf", sep = "", collapse = "")
      pdf(file = plot.filename, height = height, width = 15)
    } else if (.Platform$OS.type == "windows") {
      plot.filename <- paste(output.heatmap.plot, ".pdf", sep = "", collapse = "")
      pdf(file = plot.filename, height = height, width = 15)
    }
  }

  MSIG.HeatMapPlot(V = m, row.names = gs.names, col.labels = class.labels, col.classes = class.phen, col.names = sample.names, main = paste(title, " -- Heat Map", sep = ""), sub = " ", xlab = " ", ylab = " ", row.norm = heatmap.row.norm, cmap.type = heatmap.cmap.type)

  if (non.interactive.run == F) {
    if (.Platform$OS.type == "windows") {
      savePlot(filename = plot.filename, type = "jpeg", device = dev.cur())
    } else if (.Platform$OS.type == "unix") {
      dev.off()
    }
  } else {
    dev.off()
  }
}


MSIG.File.to.HeatMap <- function(
    input.ds,
    input.cls = "",
    output.heatmap.plot,
    output.heatmap.sorted.plot,
    title = "",
    non.interactive.run = F,
    heatmap.row.norm = F,
    phen.cmap = c("grey3", "cadetblue", "darkgreen", "chartreuse2", "red1", "darkred", "orange", "blue2", "lightblue", "pink", "coral"),
    heatmap.cmap.type = 1,
    char.rescale = 0.8) {
  print("Running MSIG.File.to.HeatMap...")

  print(heatmap.row.norm)

  # Read dataset

  dataset <- MSIG.Gct2Frame(filename = input.ds)
  m <- data.matrix(dataset$ds)
  gs.names <- dataset$row.names
  gs.descs <- dataset$descs
  sample.names <- dataset$names

  dim(m)

  Ns <- length(m[1, ])
  k.proj <- length(m[, 1])

  if (input.cls != "") {
    CLS <- ReadClsFile(file = input.cls)
    class.labels <- CLS$class.v
    class.phen <- CLS$phen
  } else {
    class.labels <- rep(1, Ns)
    class.phen <- "Samples"
  }

  # Heat map plot

  height <- ifelse(k.proj >= 25, 25, k.proj * 0.8 + 5)
  x11(width = 30, height = height)
  nf <- layout(matrix(c(1, 2), 1, 2, byrow = T), widths = c(4, 1), respect = FALSE)

  MSIG.HeatMapPlot.3(V = m, row.names = gs.names, col.labels = class.labels, col.classes = class.phen, phen.cmap = phen.cmap, col.names = sample.names, main = paste(title, " -- Heat Map", sep = ""), xlab = " ", ylab = " ", sub = " ", row.norm = heatmap.row.norm, cmap.type = heatmap.cmap.type, char.rescale = char.rescale)
  leg.txt <- class.phen
  p.vec <- rep(21, 21)
  c.vec <- phen.cmap
  par(mar = c(0, 0, 0, 0))
  plot(c(0, 0), c(1, 1), xlim = c(0, 1), ylim = c(0, 1), axes = F, type = "n", xlab = "", ylab = "")
  legend(x = 0, y = 0.8, legend = leg.txt, bty = "n", xjust = 0, yjust = 1, pch = p.vec, pt.bg = c.vec, col = "black", cex = 1.0, pt.cex = 1.5)

  savePlot(filename = output.heatmap.plot, type = "jpeg", device = dev.cur())

  height <- ifelse(k.proj >= 25, 25, k.proj * 0.8 + 5)
  x11(width = 30, height = height)
  nf <- layout(matrix(c(1, 2), 1, 2, byrow = T), widths = c(4, 1), respect = FALSE)

  dist.matrix <- dist(t(m))
  HC <- hclust(dist.matrix, method = "complete")
  m <- m[, HC$order]
  sample.names <- sample.names[HC$order]
  class.labels <- class.labels[HC$order]

  dist.matrix <- dist(m)
  HC <- hclust(dist.matrix, method = "complete")
  m <- m[HC$order, ]
  gs.names <- gs.names[HC$order]

  MSIG.HeatMapPlot.3(V = m, row.names = gs.names, col.labels = class.labels, col.classes = class.phen, phen.cmap = phen.cmap, col.names = sample.names, main = paste(title, " -- Heat Map (sorted)", sep = ""), xlab = " ", ylab = " ", sub = " ", row.norm = heatmap.row.norm, cmap.type = heatmap.cmap.type, char.rescale = char.rescale)
  leg.txt <- class.phen
  p.vec <- rep(21, 21)
  c.vec <- phen.cmap
  par(mar = c(0, 0, 0, 0))
  plot(c(0, 0), c(1, 1), xlim = c(0, 1), ylim = c(0, 1), axes = F, type = "n", xlab = "", ylab = "")
  legend(x = 0, y = 0.8, legend = leg.txt, bty = "n", xjust = 0, yjust = 1, pch = p.vec, pt.bg = c.vec, col = "black", cex = 1.0, pt.cex = 1.5)

  savePlot(filename = output.heatmap.sorted.plot, type = "jpeg", device = dev.cur())
}


MSIG.File.to.HeatMap.2 <- function(
    input.ds,
    input.cls = "",
    output.heatmap.plot,
    output.heatmap.sorted.plot,
    output.heatmap.sorted.2.plot,
    title = "",
    non.interactive.run = F,
    heatmap.row.norm = F,
    phen.cmap = c("grey3", "cadetblue", "darkgreen", "chartreuse2", "red1", "darkred", "orange", "blue2", "lightblue", "pink", "coral"),
    heatmap.cmap.type = 1,
    char.rescale = 1.0) {
  print("Running MSIG.File.to.HeatMap...")

  print(heatmap.row.norm)

  # Read dataset

  dataset <- MSIG.Gct2Frame(filename = input.ds)
  m <- data.matrix(dataset$ds)
  gs.names <- dataset$row.names
  gs.descs <- dataset$descs
  sample.names <- dataset$names

  dim(m)

  Ns <- length(m[1, ])
  k.proj <- length(m[, 1])

  if (input.cls != "") {
    CLS <- ReadClsFile(file = input.cls)
    class.labels <- CLS$class.v
    class.phen <- CLS$phen
  } else {
    class.labels <- rep(1, Ns)
    class.phen <- "Samples"
  }

  MSIG.HeatMapPlot.4(V = m, row.names = gs.names, col.labels = class.labels, col.classes = class.phen, phen.cmap = phen.cmap, col.names = sample.names, main = title, xlab = " ", ylab = " ", sub = "heat map ", row.norm = heatmap.row.norm, cmap.type = heatmap.cmap.type, char.rescale = char.rescale)

  savePlot(filename = output.heatmap.plot, type = "jpeg", device = dev.cur())

  # Sorted heat map plot

  dist.matrix <- dist(t(m))
  HC <- hclust(dist.matrix, method = "complete")
  m1 <- m[, HC$order]
  sample.names1 <- sample.names[HC$order]
  class.labels1 <- class.labels[HC$order]

  dist.matrix <- dist(m1)
  HC <- hclust(dist.matrix, method = "complete")
  m1 <- m1[HC$order, ]
  gs.names1 <- gs.names[HC$order]

  MSIG.HeatMapPlot.4(V = m1, row.names = gs.names1, col.labels = class.labels1, col.classes = class.phen, phen.cmap = phen.cmap, col.names = sample.names1, main = title, xlab = " ", ylab = " ", sub = "sorted heat map ", row.norm = heatmap.row.norm, cmap.type = heatmap.cmap.type, char.rescale = char.rescale)

  savePlot(filename = output.heatmap.sorted.plot, type = "jpeg", device = dev.cur())

  # heat map sorted inside each phenotype

  dist.matrix <- dist(m)
  HC <- hclust(dist.matrix, method = "complete")
  m2 <- m[HC$order, ]
  gs.names2 <- gs.names[HC$order]
  sample.names2 <- sample.names
  max.classes <- max(class.labels)

  for (k in 1:max.classes) {
    m3 <- m2[, class.labels == k]
    sn <- sample.names2[class.labels == k]
    dist.matrix <- dist(t(m3))
    HC <- hclust(dist.matrix, method = "complete")
    m3 <- m3[, HC$order]
    sn <- sn[HC$order]
    m2[, class.labels == k] <- m3
    sample.names2[class.labels == k] <- sn
  }

  MSIG.HeatMapPlot.4(V = m2, row.names = gs.names2, col.labels = class.labels, col.classes = class.phen, phen.cmap = phen.cmap, col.names = sample.names2, main = title, xlab = " ", ylab = " ", sub = "sorted heat map (inside phenotype)", row.norm = heatmap.row.norm, cmap.type = heatmap.cmap.type, char.rescale = char.rescale)

  savePlot(filename = output.heatmap.sorted.2.plot, type = "jpeg", device = dev.cur())
}

MSIG.Factors.Project.2 <- function(
    input.ds,
    factors.ds,
    inverse.type = "transpose", # "transpose" or "Penrose-Moore"
    weighting.type = -1, # -1 = dot product (NMF), 0 = rank (KS-like), 1 = weighted with rank and data entry
    postprojnorm = FALSE,
    output.file) {
  library(MASS)

  # start of methodology

  print("Running MSIG.Factors.Project..")

  # Read input dataset

  dataset <- MSIG.Gct2Frame(filename = input.ds)
  m <- data.matrix(dataset$ds)
  gs.names <- dataset$row.names
  gs.descs <- dataset$descs
  sample.names <- dataset$names

  # Read factors dataset

  dataset <- MSIG.Gct2Frame(filename = factors.ds)
  W <- data.matrix(dataset$ds)
  W.row.names <- dataset$row.names
  W.row.descs <- dataset$descs
  W.names <- dataset$names

  # Match features to first dataset and create matching m2 dataset

  overlap <- intersect(gs.names, W.row.names)

  locations <- match(overlap, W.row.names, nomatch = 0)
  W2 <- W[locations, ]

  # redefine m to take into account weighting

  if (inverse.type == "transpose") {
    if (weighting.type == -1) { # dot product
      m2 <- m[locations, ]
    } else if (weighting.type == 0) { # weighted by rank KS-like
      m.rank <- MSIG.NormalizeCols.Rank(m)
      m2 <- m.rank[locations, ]
    } else if (weighting.type == 0) { # weighted by rank and value
      m.rank <- MSIG.NormalizeCols.Rank(m)
      maxm <- max(m)
      minm <- min(m)
      m.weight <- (m - min(m)) / (maxm - minm)
      m2 <- m.rank[locations, ] * m.weight[locations, ]
    }
  }

  # Project input dataset using factors input

  if (inverse.type == "Penrose-Moore") {
    H <- ginv(W2) %*% m2
  } else if (inverse.type == "transpose") {
    H <- t(W2) %*% m2
  }

  # Normalize projected dataset to the unit hypersphere

  if (postprojnorm == TRUE) {
    n.col <- length(H[1, ])
    for (i in 1:n.col) {
      S.2 <- sqrt(sum(H[, i] * H[, i]))
      norm <- 1 / S.2
      H[, i] <- H[, i] / norm
    }
  }

  # Save projected dataset

  V <- data.frame(H)
  names(V) <- sample.names
  row.names(V) <- W.names
  write.gct(gct.data.frame = V, descs = gs.descs, filename = output.file)
}

MSIG.Create.W.from.Sets <- function(
    gs.db,
    output.file,
    non.interactive.run = F,
    gs.size.threshold.min = 5,
    gs.size.threshold.max = 100000) {
  print(" *** Running MSIG.Create.W.from.Sets...")

  if (.Platform$OS.type == "windows") {
    memory.limit(6000000000)
    memory.limit()
  }

  # Read input gene set database

  if (regexpr(pattern = ".gmt", gs.db[1]) == -1) {
    temp <- gs.db
  } else {
    temp <- readLines(gs.db)
  }

  max.Ng <- length(temp)
  temp.size.G <- vector(length = max.Ng, mode = "numeric")
  for (i in 1:max.Ng) {
    temp.size.G[i] <- length(unlist(strsplit(temp[[i]], "\t"))) - 2
  }

  max.size.G <- max(temp.size.G)
  gs <- matrix(rep("null", max.Ng * max.size.G), nrow = max.Ng, ncol = max.size.G)
  temp.names <- vector(length = max.Ng, mode = "character")
  temp.desc <- vector(length = max.Ng, mode = "character")
  gs.count <- 1
  for (i in 1:max.Ng) {
    gene.set.size <- length(unlist(strsplit(temp[[i]], "\t"))) - 2
    gs.line <- noquote(unlist(strsplit(temp[[i]], "\t")))
    gene.set.name <- gs.line[1]
    gene.set.desc <- gs.line[1]
    gene.set.tags <- vector(length = gene.set.size, mode = "character")
    for (j in 1:gene.set.size) {
      gene.set.tags[j] <- gs.line[j + 2]
    }
    set.size <- length(gene.set.tags)
    if ((set.size < gs.size.threshold.min) || (set.size > gs.size.threshold.max)) next
    temp.size.G[gs.count] <- set.size
    gs[gs.count, ] <- c(gene.set.tags, rep("null", max.size.G - temp.size.G[gs.count]))
    temp.names[gs.count] <- gene.set.name
    temp.desc[gs.count] <- gene.set.desc
    gs.count <- gs.count + 1
  }

  Ng <- gs.count - 1
  gs.names <- vector(length = Ng, mode = "character")
  gs.desc <- vector(length = Ng, mode = "character")
  size.G <- vector(length = Ng, mode = "numeric")
  gs.names <- temp.names[1:Ng]
  gs.desc <- temp.desc[1:Ng]
  size.G <- temp.size.G[1:Ng]

  print(c("Number of Gene Sets:", Ng))
  print(c("Original number of Gene Sets:", max.Ng))
  print(c("Maximum gene set size:", max.size.G))

  # Compute W from gene sets

  all.genes <- NULL
  for (i in 1:Ng) {
    gene.set <- gs[i, gs[i, ] != "null"]
    all.genes <- union(all.genes, gene.set)
  }

  W <- matrix(0, nrow = length(all.genes), ncol = Ng)
  for (i in 1:Ng) {
    print(paste("Computing W column for gene set:", i, sep = " "))
    gene.set <- gs[i, gs[i, ] != "null"]
    W[, i] <- as.real(sign(match(all.genes, gene.set, nomatch = 0))) # notice that the sign is 0 (no tag) or 1 (tag)
  }

  V <- data.frame(W)
  names(V) <- gs.names
  row.names(V) <- all.genes

  write.gct(gct.data.frame = V, filename = output.file)
}

MSIG.Evaluate.Projection <- function(
    input.ds,
    input.cls,
    model.set,
    prediction.results.file,
    prediction.matrix.file,
    train.pred.plot,
    test.pred.plot,
    pred.2D.plot,
    col = c("grey3", "cadetblue", "darkgreen", "chartreuse2", "red1", "darkred", "orange", "blue2", "lightblue", "pink", "coral"),
    non.interactive.run = F,
    use.feature.names = F,
    nchar.phen = 3,
    high.conf.thres = 0.7) {
  print(c("Running MSIG.Evaluate.Projection... on:", input.ds))

  library(e1071)
  library(tree)

  # Read dataset

  dataset <- MSIG.Gct2Frame(filename = input.ds)
  m <- data.matrix(dataset$ds)
  gs.names <- dataset$row.names
  gs.descs <- dataset$descs
  sample.names <- dataset$names

  dim(m)

  Ns <- length(m[1, ])
  N <- length(m[, 1])

  CLS <- MSIG.ReadClsFile(file = input.cls)
  class.labels <- CLS$class.v
  class.list <- CLS$class.list
  class.phen <- CLS$phen

  num.classes <- length(class.phen)

  print("Reading dataset completed...")

  # Use first nchar.phen characters of phenotype class to define new phenotypes

  class.list2 <- vector(length = Ns, mode = "numeric")
  for (i in 1:Ns) {
    class.list2[i] <- substr(class.list[i], 1, nchar.phen)
  }
  class.phen2 <- vector(length = num.classes, mode = "character")
  for (i in 1:num.classes) {
    class.phen2[i] <- substr(class.phen[i], 1, 3)
  }
  true.num.classes <- length(table(class.phen2))

  class.labels2 <- match(class.list2, class.phen2)

  # Separate data into train and test pieces

  m.train <- m[, model.set]
  n.train <- length(model.set)
  num.samples.train <- n.train
  sample.names.train <- as.factor(sample.names[model.set])
  class.list.train <- class.list2[model.set]
  class.phen.train <- unique(class.list.train)
  class.labels.train <- class.labels2[model.set]

  if (Ns - length(model.set) > 0) {
    m.test <- as.matrix(m[, -model.set])
    n.test <- length(m.test[1, ])
    sample.names.test <- as.factor(sample.names[-model.set])
    class.list.test <- class.list2[-model.set]
    class.phen.test <- unique(class.list.test)
    class.labels.test <- class.labels2[-model.set]
  }


  # Build SVM and tree models

  print("Building SVM model...")

  svm.model <- svm(x = t(m.train), y = class.list.train, scale = F, tolerance = 0.0001, type = "C-classification", kernel = "linear", cost = 3, probability = T)

  #   tree.data <- data.frame(as.factor(class.list), t(h.df))
  #   tree.model <- tree(formula = as.factor.class.list. ~ ., data = tree.data, split = "deviance")


  #   tree.data <- data.frame(cbind(class.list.train, t(m.train)))
  #   tree.data
  #   print(tree.data)
  #   tree.model <- tree(formula = class.list.train ~ ., data = tree.data)
  #   print(summary(tree.model))
  #   x11(height = 15, width = 15)
  #   plot(tree.model)
  #   text(tree.model)

  print("Computing train set predictions...")

  train.pred <- predict(object = svm.model, newdata = t(m.train), decision.values = T, probability = T)
  dec.vals.train <- attr(train.pred, "decision.values")
  prob.train <- signif(attr(train.pred, "probabilities"), digits = 2)
  confidence.vector <- apply(prob.train, 1, max)
  confidence.call <- ifelse(confidence.vector < high.conf.thres, " L ", " H ")
  error.call <- ifelse(class.list.train == as.character(train.pred), "   ", " * ")
  col.symbols.train <- paste(confidence.call, error.call)

  class.names <- names(data.frame(prob.train))
  bscore <- vector(length = n.train, mode = "numeric")
  for (i in 1:n.train) {
    bscore[i] <- 0
    for (k in 1:length(prob.train[1, ])) {
      if (class.list.train[i] == class.names[k]) {
        bscore[i] <- bscore[i] + (1 - prob.train[i, k])^2
      } else {
        bscore[i] <- bscore[i] + prob.train[i, k]^2
      }
    }
    bscore[i] <- signif(bscore[i], digits = 2)
  }
  Brier.train <- signif(mean(bscore), digits = 2)

  train.results <- data.frame(cbind(as.character(sample.names.train), class.list.train, as.character(train.pred), error.call, confidence.call, prob.train, bscore))
  names(train.results)[1] <- "Train Sample Name"
  names(train.results)[2] <- "Actual"
  names(train.results)[3] <- "Predicted"
  names(train.results)[4] <- "Error (*)"
  names(train.results)[5] <- "Conf (H/L)"

  names(train.results)[6 + length(class.phen.train)] <- "Brier score"
  print(train.results)
  print(c("Brier score (Train) = ", Brier.train))

  write("Training Results \n", file = prediction.results.file, append = F)
  write.table(train.results, file = prediction.results.file, append = T, quote = F, row.names = F, sep = "\t")

  write(c("\n\n Brier score (Train) = ", Brier.train), file = prediction.results.file, append = T)

  conf.table.train <- table(class.list.train, train.pred)
  print(conf.table.train)
  write("\n\n Confusion Matrix (Train) \n", file = prediction.results.file, append = T)
  write.table(conf.table.train, file = prediction.results.file, append = T, quote = F, row.names = F, sep = "\t")
  height <- ifelse(length(class.phen.train) > 50, 15, 0.30 * length(class.phen.train) + 4.8)
  if (non.interactive.run == F) {
    if (.Platform$OS.type == "windows") {
      plot.filename <- train.pred.plot
      x11(height = height, width = 15)
    } else if (.Platform$OS.type == "unix") {
      plot.filename <- paste(train.pred.plot, ".pdf", sep = "", collapse = "")
      pdf(file = plot.filename, height = height, width = 15)
    }
  } else {
    if (.Platform$OS.type == "unix") {
      plot.filename <- paste(train.pred.plot, ".pdf", sep = "", collapse = "")
      pdf(file = plot.filename, height = height, width = 15)
    } else if (.Platform$OS.type == "windows") {
      plot.filename <- paste(train.pred.plot, ".pdf", sep = "", collapse = "")
      pdf(file = plot.filename, height = height, width = 15)
    }
  }

  MSIG.HeatMapPlot.2(V = t(prob.train), row.names = class.names, col.symbols = col.symbols.train, col.labels = class.list.train, col.names = as.character(sample.names.train), main = "Train Samples Predictions", sub = " ", xlab = " ", ylab = " ", row.norm = F, cmap.type = 2)
  #   MSIG.HeatMapPlot.3(V = t(prob.train), row.names = class.names, col.labels = class.list.train, col.names =as.character(sample.names.train), phen.cmap = col[1:length(class.names)], main= "Train Samples Predictions", sub = " ", xlab=" ", ylab=" ", row.norm = F,  cmap.type = 2)

  if (non.interactive.run == F) {
    if (.Platform$OS.type == "windows") {
      savePlot(filename = train.pred.plot, type = "jpeg", device = dev.cur())
    } else if (.Platform$OS.type == "unix") {
      dev.off()
    }
  } else {
    dev.off()
  }

  print("Building SVM model completed. Predicting test data...")

  if (Ns - length(model.set) > 0) {
    test.pred <- predict(object = svm.model, newdata = t(m.test), decision.values = T, probability = T)
    dec.vals.test <- attr(test.pred, "decision.values")
    prob.test <- signif(attr(test.pred, "probabilities"), digits = 2)
    confidence.vector <- apply(prob.test, 1, max)
    confidence.call <- ifelse(confidence.vector < high.conf.thres, " L ", " H ")
    error.call <- ifelse(class.list.test == as.character(test.pred), "   ", " * ")
    col.symbols.test <- paste(confidence.call, error.call)

    class.names <- names(data.frame(prob.test))
    bscore <- vector(length = n.test, mode = "numeric")
    for (i in 1:n.test) {
      bscore[i] <- 0
      for (k in 1:length(prob.test[1, ])) {
        if (class.list.test[i] == class.names[k]) {
          bscore[i] <- bscore[i] + (1 - prob.test[i, k])^2
        } else {
          bscore[i] <- bscore[i] + prob.test[i, k]^2
        }
      }
      bscore[i] <- signif(bscore[i], digits = 2)
    }
    Brier.test <- signif(mean(bscore), digits = 2)


    test.results <- data.frame(cbind(as.character(sample.names.test), class.list.test, as.character(test.pred), error.call, confidence.call, prob.test, bscore))
    names(test.results)[1] <- "Test Sample Name"
    names(test.results)[2] <- "Actual"
    names(test.results)[3] <- "Predicted"
    names(test.results)[4] <- "Error (*)"
    names(test.results)[5] <- "Conf (H/L)"

    #      names(test.results)[6 + length(class.phen.test)] <- "Brier score"
    names(test.results)[6 + length(class.phen.train)] <- "Brier score"
    print(test.results)
    print(c("Brier score (Test) = ", Brier.test))

    write("\n Test Results \n", file = prediction.results.file, append = T)
    write.table(test.results, file = prediction.results.file, append = T, quote = F, row.names = F, sep = "\t")

    write(c("\n\n Brier score (Test) = ", Brier.test), file = prediction.results.file, append = T)

    conf.table.test <- table(class.list.test, test.pred)
    print(conf.table.test)
    write("\n\n Confusion Matrix (Test) \n", file = prediction.results.file, append = T)
    write.table(conf.table.test, file = prediction.results.file, append = T, quote = F, row.names = F, sep = "\t")
    #      height <- ifelse(length(class.phen.test) > 50, 15, 0.20*length(class.phen.test) + 4.8)
    height <- ifelse(length(class.phen.train) > 50, 15, 0.30 * length(class.phen.train) + 4.8)
    if (non.interactive.run == F) {
      if (.Platform$OS.type == "windows") {
        plot.filename <- test.pred.plot
        x11(height = height, width = 15)
      } else if (.Platform$OS.type == "unix") {
        plot.filename <- paste(test.pred.plot, ".pdf", sep = "", collapse = "")
        pdf(file = plot.filename, height = height, width = 15)
      }
    } else {
      if (.Platform$OS.type == "unix") {
        plot.filename <- paste(test.pred.plot, ".pdf", sep = "", collapse = "")
        pdf(file = plot.filename, height = height, width = 15)
      } else if (.Platform$OS.type == "windows") {
        plot.filename <- paste(test.pred.plot, ".pdf", sep = "", collapse = "")
        pdf(file = plot.filename, height = height, width = 15)
      }
    }

    #      MSIG.HeatMapPlot.2(V = t(prob.test), row.names = names(test.results)[seq(6, 6 + length(class.phen.test) - 1)], col.symbols = col.symbols.test, col.labels = class.list.test, col.names =as.character(sample.names.test), main= "Test Samples Predictions", sub = " ", xlab=" ", ylab=" ", row.norm = T,  cmap.type = 2)
    MSIG.HeatMapPlot.2(V = t(prob.test), row.names = names(test.results)[seq(6, 6 + length(class.phen.train) - 1)], col.symbols = col.symbols.test, col.labels = class.list.test, col.names = as.character(sample.names.test), main = "Test Samples Predictions", sub = " ", xlab = " ", ylab = " ", row.norm = F, cmap.type = 2)

    if (non.interactive.run == F) {
      if (.Platform$OS.type == "windows") {
        savePlot(filename = test.pred.plot, type = "jpeg", device = dev.cur())
      } else if (.Platform$OS.type == "unix") {
        dev.off()
      }
    } else {
      dev.off()
    }
  }


  # Save predictions for all classes in a gct and cls files

  print(c("dim train=", dim(prob.train)))
  print(c("dim test=", dim(prob.test)))

  V <- cbind(t(prob.train), t(prob.test))
  W <- data.frame(V)
  names(W) <- c(as.character(sample.names.train), as.character(sample.names.test))
  row.names(W) <- names(train.results)[seq(6, 6 + length(class.phen.train) - 1)]

  print(W)
  write.gct(gct.data.frame = W, descs = row.names(W), filename = prediction.matrix.file)


  print("Done predicting test data...")

  #################################################################################
  # Contour plot results

  pca <- prcomp(t(m.train), retx = TRUE, center = TRUE, scale. = TRUE)

  S1 <- pca$x[, 1]
  S2 <- pca$x[, 2]
  S3 <- pca$x[, 3]
  X1 <- pca$rotation[, 1]
  X2 <- pca$rotation[, 2]
  X3 <- pca$rotation[, 3]

  # 2D plots

  c0 <- col
  c1 <- colors()[match(c0, colors())]
  color <- c1[class.labels]

  height <- 15
  if (non.interactive.run == F) {
    if (.Platform$OS.type == "windows") {
      plot.filename <- pred.2D.plot
      x11(height = height, width = 15)
    } else if (.Platform$OS.type == "unix") {
      plot.filename <- paste(pred.2D.plot, ".pdf", sep = "", collapse = "")
      pdf(file = plot.filename, height = height, width = 15)
    }
  } else {
    if (.Platform$OS.type == "unix") {
      plot.filename <- paste(pred.2D.plot, ".pdf", sep = "", collapse = "")
      pdf(file = plot.filename, height = height, width = 15)
    } else if (.Platform$OS.type == "windows") {
      plot.filename <- paste(pred.2D.plot, ".pdf", sep = "", collapse = "")
      pdf(file = plot.filename, height = height, width = 15)
    }
  }

  nf <- layout(matrix(c(1, 2), 1, 2, byrow = T), widths = c(3, 1), heights = 1, respect = FALSE)

  if (Ns - length(model.set) > 0) {
    test.scores <- predict(pca, t(m.test))

    S1 <- c(pca$x[, 1], test.scores[, 1])
    S2 <- c(pca$x[, 2], test.scores[, 2])
    S3 <- c(pca$x[, 3], test.scores[, 3])
  }

  max.S <- max(sqrt(S1 * S1 + S2 * S2))
  max.X <- max(sqrt(X1 * X1 + X2 * X2))
  X1 <- max.S * X1 / max.X
  X2 <- max.S * X2 / max.X
  max.A <- max(max.S, max.X)
  num.samples <- length(S1)

  plot(S1, S2, xlim = c(-max.A, max.A), ylim = c(-max.A, max.A), type = "n", main = "  ", sub = input.ds)

  for (j in 1:num.samples) {
    if (min(class.labels) == 0) {
      color.code <- c1[class.labels[j] + 1]
    } else {
      color.code <- c1[class.labels[j]]
    }
    points(S1[j], S2[j], pch = 21, type = "p", cex = 2, bg = color.code, col = "black")
  }

  for (j in 1:N) {
    x.coor <- X1[j] * 0.925
    y.coor <- X2[j] * 0.925
    arrows(0, 0, x.coor, y.coor, lwd = 1, length = 0.15, angle = 20, col = "grey50")
    if (use.feature.names == FALSE) {
      leg.txt <- paste("F", j, sep = "")
    } else {
      leg.txt <- gs.names[j]
    }
    text(X1[j], X2[j], labels = leg.txt, adj = NULL, pos = NULL, offset = 0.5, vfont = NULL, cex = 1, col = "grey50")
  }


  # compute grid

  max.points <- 1000
  max.d <- apply(m, 1, max)
  grid <- matrix(0, nrow = N, ncol = max.points)
  for (i in 1:max.points) {
    for (j in 1:N) {
      grid[j, i] <- runif(1, 0, 1.2 * max.d[j])
    }
  }

  grid.pred <- predict(object = svm.model, newdata = t(grid), decision.values = T, probability = T)
  prob.test <- attr(grid.pred, "probabilities")

  grid.df <- data.frame(t(grid))
  names(grid.df) <- gs.names

  test.scores2 <- predict(pca, grid.df)

  #  test.scores2 <- predict(pca, t(grid))

  S1 <- test.scores2[, 1]
  S2 <- test.scores2[, 2]

  for (j in 1:max.points) {
    if (min(class.labels) == 0) {
      color.code <- c1[class.labels[match(grid.pred[j], class.list2)] + 1]
    } else {
      color.code <- c1[class.labels[match(grid.pred[j], class.list2)]]
    }
    points(S1[j], S2[j], pch = 21, type = "p", cex = 1, bg = color.code, col = color.code)
  }

  # legend

  leg.txt <- class.phen
  n.phen <- length(class.phen)
  p.vec <- rep(21, n.phen)
  c.vec <- c1[1:n.phen]

  par(mar = c(0, 0, 0, 0))
  plot(c(0, 0), c(1, 1), xlim = c(0, 1), ylim = c(0, 1), axes = F, type = "n", xlab = "", ylab = "")
  legend(x = 0, y = 1, legend = leg.txt, bty = "n", xjust = 0, yjust = 1, pch = p.vec, pt.bg = c.vec, col = "black", cex = 2, pt.cex = 3)

  if (non.interactive.run == F) {
    if (.Platform$OS.type == "windows") {
      savePlot(filename = pred.2D.plot, type = "jpeg", device = dev.cur())
    } else if (.Platform$OS.type == "unix") {
      dev.off()
    }
  } else {
    dev.off()
  }

  # }


  MSIG.HeatMapPlot.3(V = t(prob.train), row.names = class.names, col.labels = class.labels.train, col.names = as.character(error.call), col.classes = class.names, phen.cmap = col[1:length(class.names)], main = "Train Samples Predictions", sub = " ", xlab = " ", ylab = " ", row.norm = F, cmap.type = 2)

  #      MSIG.HeatMapPlot.5(V = t(prob.train), row.names = class.names, col.labels = class.labels.train, col.names =as.character(error.call), col.classes = class.names, main= "Train Samples Predictions", sub = " ", xlab=" ", ylab=" ", row.norm = F, cmap.type = 2)

  # legend

  leg.txt <- class.phen
  n.phen <- length(class.phen)
  p.vec <- symbs[1:n.phen]
  c.vec <- col[1:n.phen]

  par(mar = c(0, 0, 0, 0))
  plot(c(0, 0), c(1, 1), xlim = c(0, 1), ylim = c(0, 1), axes = F, type = "n", xlab = "", ylab = "")
  legend(x = 0, y = 1, legend = leg.txt, bty = "n", xjust = 0, yjust = 1, pch = p.vec, pt.bg = c.vec, col = "black", cex = symbol.scaling * 1.25, pt.cex = symbol.scaling * 2.5)

  if (non.interactive.run == F) {
    if (.Platform$OS.type == "windows") {
      savePlot(filename = train.pred.plot, type = "jpeg", device = dev.cur())
    } else if (.Platform$OS.type == "unix") {
      dev.off()
    }
  } else {
    dev.off()
  }

  print("Building SVM model completed. Predicting test data...")

  if (Ns - length(model.set) > 0) {
    test.pred <- predict(object = svm.model, newdata = t(m.test), decision.values = T, probability = T)
    dec.vals.test <- attr(test.pred, "decision.values")
    prob.test <- signif(attr(test.pred, "probabilities"), digits = 2)
    confidence.vector <- vector(length = n.test, mode = "numeric")
    bscore <- vector(length = n.test, mode = "numeric")
    max.k <- length(prob.train[1, ])
    random.pred.conf <- ((max.k - 1) / max.k)^2 + (max.k - 1) * (1 / max.k)^2
    for (ii in 1:n.test) {
      probs <- sort(prob.test[ii, ], decreasing = T)
      confidence.vector[ii] <- 1 - ((1 - probs[1])^2 + sum(probs[2:max.k]^2)) / random.pred.conf
      confidence.vector[ii] <- signif(confidence.vector[ii], digits = 3)
      if (class.list.test[ii] == as.character(test.pred[ii])) {
        bscore[ii] <- signif((1 - probs[1])^2, digits = 2)
      } else {
        bscore[ii] <- signif(probs[1]^2, digits = 2)
      }
    }

    confidence.call <- ifelse(confidence.vector >= high.conf.thres, " H ", " L ")
    error.call <- ifelse(class.list.test == as.character(test.pred), "   ", " * ")
    no.call <- ifelse(confidence.vector >= high.conf.thres, 0, 1)
    real.error <- ifelse(((no.call == 0) & (error.call == " * ")), 1, 0)
    correct.call <- ifelse(((no.call == 0) & (error.call == "   ")), 1, 0)
    col.symbols.test <- paste(confidence.call, error.call)
    class.names <- names(data.frame(prob.test))
    Brier.test <- signif(mean(bscore), digits = 2)

    test.results <- data.frame(cbind(as.character(sample.names.test), class.list.test, as.character(test.pred), error.call, confidence.call, confidence.vector, no.call, real.error, correct.call, prob.test, bscore))
    names(test.results)[1] <- "Test Sample Name"
    names(test.results)[2] <- "Actual"
    names(test.results)[3] <- "Predicted"
    names(test.results)[4] <- "Error (*)"
    names(test.results)[5] <- "Conf (H/L)"
    names(test.results)[6] <- "Conf"
    names(test.results)[7] <- "No Call"
    names(test.results)[8] <- "Real Error"
    names(test.results)[9] <- "Correct Call"

    names(test.results)[10 + length(class.phen.train)] <- "Brier score"
    #      print(test.results)
    print(c("Brier score (Test) = ", Brier.test))

    write("\n Test Results \n", file = prediction.results.file, append = T)
    write.table(test.results, file = prediction.results.file, append = T, quote = F, row.names = T, sep = "\t")

    write(c("\n\n Brier score (Test) = ", Brier.test), file = prediction.results.file, append = T)

    no.call.list <- split(no.call, class.list.test)
    real.error.list <- split(real.error, class.list.test)
    correct.call.list <- split(correct.call, class.list.test)
    count.class <- c(sapply(no.call.list, length), length(no.call))
    no.call.class <- c(sapply(no.call.list, sum), sum(no.call))
    real.error.class <- c(sapply(real.error.list, sum), sum(real.error))
    correct.call.class <- c(sapply(correct.call.list, sum), sum(correct.call))

    no.call.class.pct <- no.call.class / count.class
    real.error.class.pct <- real.error.class / count.class
    correct.call.class.pct <- correct.call.class / count.class

    perf.table.test <- data.frame(cbind(count.class, no.call.class, no.call.class.pct, real.error.class, real.error.class.pct, correct.call.class, correct.call.class.pct))
    names(perf.table.test) <- c("Count", "No Call", "No Call (%)", "Real Error", "Real Error (%)", "Correct Call", "Correct Call (%)")
    row.names(perf.table.test) <- c(names(no.call.list), "Total")
    write.table(perf.table.test, file = prediction.results.file, append = T, quote = F, row.names = T, sep = "\t")
    print(perf.table.test)

    conf.table.test <- table(class.list.test, test.pred)
    print(conf.table.test)
    write("\n\n Confusion Matrix (Test) \n", file = prediction.results.file, append = T)
    write.table(conf.table.test, file = prediction.results.file, append = T, quote = F, row.names = T, sep = "\t")
    height <- ifelse(length(class.phen.train) > 50, 20, 0.30 * length(class.phen.train) + 10)

    if (non.interactive.run == F) {
      if (.Platform$OS.type == "windows") {
        plot.filename <- test.pred.plot
        x11(height = height, width = 40)
      } else if (.Platform$OS.type == "unix") {
        plot.filename <- paste(test.pred.plot, ".pdf", sep = "", collapse = "")
        pdf(file = plot.filename, height = height, width = 40)
      }
    } else {
      if (.Platform$OS.type == "unix") {
        plot.filename <- paste(test.pred.plot, ".pdf", sep = "", collapse = "")
        pdf(file = plot.filename, height = height, width = 40)
      } else if (.Platform$OS.type == "windows") {
        plot.filename <- paste(test.pred.plot, ".pdf", sep = "", collapse = "")
        pdf(file = plot.filename, height = height, width = 40)
      }
    }

    nf <- layout(matrix(c(1, 2), 1, 2, byrow = T), widths = c(5, 1), heights = 1, respect = FALSE)

    MSIG.HeatMapPlot.3(V = t(prob.test), row.names = names(test.results)[seq(10, 10 + length(class.phen.train) - 1)], col.labels = class.labels.test, col.names = as.character(error.call), col.classes = class.names, phen.cmap = col[1:length(class.names)], main = "Test Samples Predictions", sub = " ", xlab = " ", ylab = " ", row.norm = F, cmap.type = 2)


    # legend

    leg.txt <- class.phen
    n.phen <- length(class.phen)
    p.vec <- symbs[1:n.phen]
    c.vec <- col[1:n.phen]

    par(mar = c(0, 0, 0, 0))
    plot(c(0, 0), c(1, 1), xlim = c(0, 1), ylim = c(0, 1), axes = F, type = "n", xlab = "", ylab = "")
    legend(x = 0, y = 1, legend = leg.txt, bty = "n", xjust = 0, yjust = 1, pch = p.vec, pt.bg = c.vec, col = "black", cex = symbol.scaling * 1.25, pt.cex = symbol.scaling * 2.5)

    if (non.interactive.run == F) {
      if (.Platform$OS.type == "windows") {
        savePlot(filename = test.pred.plot, type = "jpeg", device = dev.cur())
      } else if (.Platform$OS.type == "unix") {
        dev.off()
      }
    } else {
      dev.off()
    }
  }


  # Save predictions for all classes in a gct and cls files

  #   print(c("dim train=", dim(prob.train)))
  #   print(c("dim test=", dim(prob.test)))

  if (Ns - length(model.set) > 0) {
    V <- cbind(t(prob.train), t(prob.test))
    W <- data.frame(V)
    names(W) <- c(as.character(sample.names.train), as.character(sample.names.test))
    row.names(W) <- names(train.results)[seq(7, 7 + length(class.phen.train) - 1)]
  } else {
    V <- t(prob.train)
    W <- data.frame(V)
    names(W) <- as.character(sample.names.train)
    row.names(W) <- names(train.results)[seq(7, 7 + length(class.phen.train) - 1)]
  }
  #   print(W)
  write.gct(gct.data.frame = W, descs = row.names(W), filename = prediction.matrix.file)

  print("Done predicting test data...")

  #################################################################################
  # Contour plot results

  pca <- prcomp(t(m.train), retx = TRUE, center = T, scale. = T)
  S1 <- pca$x[, 1]
  S2 <- pca$x[, 2]
  S3 <- pca$x[, 3]
  X1 <- pca$rotation[, 1]
  X2 <- pca$rotation[, 2]
  X3 <- pca$rotation[, 3]

  row.mean <- apply(m.train, MARGIN = 1, FUN = mean)
  row.sd <- apply(m.train, MARGIN = 1, FUN = sd)

  # 2D plots

  c0 <- col
  #   c1 <- colors()[match(c0, colors())]
  c1 <- col
  color <- c1[class.labels]

  height <- 25
  if (non.interactive.run == F) {
    if (.Platform$OS.type == "windows") {
      plot.filename <- pred.2D.plot
      x11(height = height, width = 30)
    } else if (.Platform$OS.type == "unix") {
      plot.filename <- paste(pred.2D.plot, ".pdf", sep = "", collapse = "")
      pdf(file = plot.filename, height = height, width = 30)
    }
  } else {
    if (.Platform$OS.type == "unix") {
      plot.filename <- paste(pred.2D.plot, ".pdf", sep = "", collapse = "")
      pdf(file = plot.filename, height = height, width = 30)
    } else if (.Platform$OS.type == "windows") {
      plot.filename <- paste(pred.2D.plot, ".pdf", sep = "", collapse = "")
      pdf(file = plot.filename, height = height, width = 30)
    }
  }

  nf <- layout(matrix(c(1, 2), 1, 2, byrow = T), widths = c(3, 1), heights = 1, respect = FALSE)

  if (Ns - length(model.set) > 0) {
    test.scores <- predict(pca, t(m.test))
    S1 <- c(pca$x[, 1], test.scores[, 1])
    S2 <- c(pca$x[, 2], test.scores[, 2])
    S3 <- c(pca$x[, 3], test.scores[, 3])
  }

  max.S <- max(sqrt(S1 * S1 + S2 * S2))
  max.X <- max(sqrt(X1 * X1 + X2 * X2))
  X1 <- max.S * X1 / max.X
  X2 <- max.S * X2 / max.X
  X3 <- max.S * X3 / max.X
  max.A <- max(max.S, max.X)
  num.samples <- length(S1)

  plot(S1, S2, xlim = c(-max.A, max.A), ylim = c(-max.A, max.A), type = "n", main = "  ", sub = input.ds)

  for (j in 1:num.samples) {
    if (min(class.labels) == 0) {
      symb <- symbs[class.labels[j] + 1]
      color.code <- c1[class.labels[j] + 1]
    } else {
      symb <- symbs[class.labels[j]]
      color.code <- c1[class.labels[j]]
    }
    points(S1[j], S2[j], pch = symb, type = "p", cex = symbol.scaling * 3, bg = color.code, col = "black")
    #      if (j <= length(model.set)) {
    #         points(S1[j], S2[j], pch=22, type="p", cex = 2, bg = color.code, col = "black")
    #       } else {
    #         points(S1[j], S2[j], pch=21, type="p", cex = 1.5, bg = color.code, col = "black")
    #       }
  }

  for (j in 1:N) {
    x.coor <- X1[j] * 0.925
    y.coor <- X2[j] * 0.925
    arrows(0, 0, x.coor, y.coor, lwd = 1, length = 0.15, angle = 20, col = "grey50")
    if (use.feature.names == FALSE) {
      leg.txt <- paste("F", j, sep = "")
    } else {
      leg.txt <- gs.names[j]
    }
    text(X1[j], X2[j], labels = leg.txt, adj = NULL, pos = NULL, offset = 0.5, vfont = NULL, cex = 1, col = "grey50")
  }

  # compute grid

  points.axis <- 200
  x <- vector(length = points.axis, mode = "numeric")
  pca.x <- matrix(0, nrow = 2, ncol = points.axis * points.axis)

  for (i in 1:points.axis) {
    x[i] <- -max.A + i * 2 * max.A / points.axis
  }
  for (i in 1:points.axis) {
    for (j in 1:points.axis) {
      index.point <- i + (j - 1) * points.axis
      pca.x[1, index.point] <- -max.A + i * 2 * max.A / points.axis
      pca.x[2, index.point] <- -max.A + j * 2 * max.A / points.axis
    }
  }

  #   I <- pca$rotation %*% t(pca$x)
  #   I[1,1]*row.sd[1] + row.mean[1]

  grid.H <- pca$rotation[, 1:2] %*% pca.x
  for (i in 1:N) {
    grid.H[i, ] <- grid.H[i, ] * row.sd[i] + row.mean[i]
  }

  grid.pred <- predict(object = svm.model, newdata = t(grid.H), decision.values = T, probability = T)

  #   print(c("grid.pred:", grid.pred))

  prob.test <- attr(grid.pred, "probabilities")

  z <- matrix(0, nrow = points.axis, ncol = points.axis)
  z.class <- array(dim = c(length(class.phen.train), points.axis, points.axis))
  max.k <- length(prob.train[1, ])
  random.pred.conf <- ((max.k - 1) / max.k)^2 + (max.k - 1) * (1 / max.k)^2
  for (i in 1:points.axis) {
    for (j in 1:points.axis) {
      index.point <- i + (j - 1) * points.axis
      probs <- sort(prob.test[index.point, ], decreasing = T)
      z[i, j] <- 1 - ((1 - probs[1])^2 + sum(probs[2:max.k]^2)) / random.pred.conf
      for (k in 1:length(class.phen.train)) {
        if (probs[1] == prob.test[index.point, k]) {
          z.class[k, i, j] <- z[i, j]
        } else {
          z.class[k, i, j] <- 0
        }
      }

      #       probs <- sort(prob.test[index.point,], decreasing=T)
      #       z[i, j] <- probs[1] - probs[2]
      #       for (k in 1:length(class.phen.train)) {
      #          if (probs[1] == prob.test[index.point, k]) {
      #             z.class[k, i, j] <- z[i, j]
      #           } else {
      #             z.class[k, i, j] <- 0
      #           }
      #       }
    }
  }

  #   contour(x, x, z.class[3,,], nlevels = 50, col="black", lwd=2, add=T)

  library("RColorBrewer")
  if (length(levels) > 1) {
    contour(x, x, z, levels = levels, col = brewer.pal(n = 9, name = "Greys")[7], add = T)
  } else {
    contour(x, x, z, nlevels = nlevels, col = brewer.pal(n = 9, name = "Greys")[7], add = T)
  }
  contour(x, x, z, levels = 0.01, col = brewer.pal(n = 9, name = "Greys")[7], lwd = 2, add = T)
  for (k in 1:length(class.phen.train)) {
    contour(x, x, z.class[k, , ], levels = high.conf.thres, col = col[k], lwd = 2, add = T)
    #       contour(x, x, z.class[k,,], nlevels = 20, col=col[k], lwd=2, add=T)
  }

  #   z <- matrix(0, nrow=points.axis, ncol=points.axis)
  #     for (k in 1:length(class.phen.train)) {
  #     for (i in 1:points.axis) {
  #       for (j in 1:points.axis) {
  #        index.point <- i + (j - 1) * points.axis
  #         z[i, j] <- prob.test[index.point, k]
  #       }
  #     }
  #    contour(x, x, z, nlevels = 20, col=c1[k], add=T)
  #   }

  #   print("pca.x:")
  #   print(pca.x)

  #   for (j in seq(1, points.axis*points.axis)) {
  #      if (min(class.labels) == 0) {
  #          color.code <- c1[class.labels[match(grid.pred[j], class.list2)] + 1]
  #      } else {

  #          color.code <- c1[class.labels[match(grid.pred[j], class.list2)]]
  #      }
  #      points(pca.x[1, j], pca.x[2, j], pch=21, type="p", cex = 1, bg = color.code, col = color.code)
  #   }

  #   print(c("class.phen=", class.phen, " class.phen.train=", class.phen.train))
  #   print(c("prob.test:", prob.test))

  # legend

  leg.txt <- class.phen
  n.phen <- length(class.phen)
  #   p.vec <- c(rep(22, length(class.phen.train)), rep(21, n.phen - length(class.phen.train)))
  p.vec <- symbs[1:n.phen]
  c.vec <- c1[1:n.phen]

  par(mar = c(0, 0, 0, 0))
  plot(c(0, 0), c(1, 1), xlim = c(0, 1), ylim = c(0, 1), axes = F, type = "n", xlab = "", ylab = "")
  legend(x = 0, y = 1, legend = leg.txt, bty = "n", xjust = 0, yjust = 1, pch = p.vec, pt.bg = c.vec, col = "black", cex = symbol.scaling * 1.15, pt.cex = symbol.scaling * 3)

  if (non.interactive.run == F) {
    if (.Platform$OS.type == "windows") {
      savePlot(filename = pred.2D.plot, type = "jpeg", device = dev.cur())
    } else if (.Platform$OS.type == "unix") {
      dev.off()
    }
  } else {
    dev.off()
  }
}

MSIG.StJude_to_GCT <- function(
    input.ds,
    output.ds,
    output.cls) {
  # start of methodology

  print(c("Running MSIG.StJude_to_GCT... on:", input.ds))

  # Read input datasets


  filename <-
    ds <- read.delim(filename, header = T, sep = "\t", skip = 1, row.names = 1, blank.lines.skip = T, comment.char = "", as.is = T, na.strings = "")
  descs <- ds[, 1]
  ds <- ds[-1]
  row.names <- row.names(ds)
  names <- names(ds)
  return(list(ds = ds, row.names = row.names, descs = descs, names = names))

  dataset <- MSIG.Gct2Frame(filename = input.ds)
  m <- data.matrix(dataset$ds)
  gs.names <- dataset$row.names
  gs.descs <- dataset$descs
  sample.names <- dataset$names

  # Read CLS file

  if (input.cls != "") {
    CLS <- ReadClsFile(file = input.cls)
    class.labels <- CLS$class.v
    class.phen <- CLS$phen
  }

  # Select desired column subset

  if (column.subset[1] == "ALL") {
    m2 <- m
    sample.names2 <- sample.names
    if (input.cls != "") {
      class.labels2 <- class.labels
    }
  } else {
    m2 <- m[, column.subset]
    if (is.numeric(column.subset[1])) {
      sample.names2 <- sample.names[column.subset]
      if (input.cls != "") {
        class.labels2 <- class.labels[column.subset]
      }
    } else {
      locations <- match(column.subset, sample.names)
      sample.names2 <- sample.names[locations]
      if (input.cls != "") {
        class.labels2 <- class.labels[locations]
      }
    }
  }

  if (row.subset[1] == "ALL") {
    m3 <- m2
    gs.names2 <- gs.names
    gs.descs2 <- gs.descs
  } else {
    m3 <- m2[row.subset, ]
    if (is.numeric(row.subset[1])) {
      gs.names2 <- gs.names[row.subset]
      gs.descs2 <- gs.descs[row.subset]
    } else {
      locations <- match(row.subset, gs.names)
      gs.names2 <- gs.names[locations]
      gs.descs2 <- gs.descs[locations]
    }
  }

  # Save datasets

  V <- data.frame(m3)
  names(V) <- sample.names2
  row.names(V) <- gs.names2
  write.gct(gct.data.frame = V, descs = gs.descs2, filename = output.ds)

  if (input.cls != "") {
    write.cls(class.v = class.labels2, phen = class.phen, filename = output.cls)
  }
}



MSIG.W.to.gmt <- function(
    input.ds,
    threshold = 1.0,
    output.gmt) {
  dataset <- MSIG.Gct2Frame(filename = input.ds)
  m <- data.matrix(dataset$ds)
  N <- length(m[, 1])
  M <- length(m[1, ])
  gs.names <- dataset$row.names
  gs.descs <- dataset$descs
  sample.names <- dataset$names


  for (i in 1:M) {
    gene.set <- gs.names[m[, i] >= threshold]

    row.header <- paste("Factor_", i, sep = "")
    output.line <- paste(gene.set, sep = "\t", collapse = "\t")
    output.line <- paste(row.header, row.header, output.line, sep = "\t", collapse = "")

    print(paste("factor =", i, " length=", length(gene.set), sep = ""))

    if (i == 1) {
      write(noquote(output.line), file = output.gmt, append = F, ncolumns = length(gene.set) + 2)
    } else {
      write(noquote(output.line), file = output.gmt, append = T, ncolumns = length(gene.set) + 2)
    }
  }
}

MSIG.StJude2gct <- function(input.ds, output.ds) {
  # converts a file from St Jude's format to GCT

  header.cont <- readLines(input.ds, n = 1)
  temp <- unlist(strsplit(header.cont, "\t"))
  colst <- length(temp)
  header.labels <- temp[seq(3, colst, 3)]
  ds <- read.delim(input.ds, header = F, row.names = 1, sep = "\t", skip = 2, blank.lines.skip = T, comment.char = "", as.is = T)
  colst <- length(ds[1, ])
  cols <- (colst - 1) / 3
  rows <- length(ds[, 1])
  A <- matrix(nrow = rows - 1, ncol = cols)
  A <- ds[1:rows, seq(2, colst, 3)]
  descs <- ds[, 1]
  table1 <- data.frame(A)
  names(table1) <- header.labels

  write.gct(gct.data.frame = table1, descs = descs, output.ds)
}

MSIG.Projection.Plots.3 <- function(
    input.ds,
    input.cls = "",
    model.set = "ALL",
    output.2D.proj.file,
    output.2D.proj.plot,
    output.3D.proj.file,
    output.3D.1.proj.plot,
    output.3D.2.proj.plot,
    output.3D.3.proj.plot,
    output.heatmap.plot,
    output.heatmap.sorted.plot,
    output.heatmap.sorted.2.plot,
    output.hclust.plot,
    use.feature.names = FALSE,
    use.biplot = TRUE,
    title = "",
    seed = 1234,
    non.interactive.run = F,
    heatmap.row.norm = T,
    heatmap.cmap.type = 1,
    symbol.scaling = 1,
    col = c("greny3", "cadetblue", "darkgreen", "chartreuse2", "red1", "darkred", "orange", "blue2", "lightblue", "pink", "coral"),
    symbs = c(22, 21, 20, 23, 24, 25, 21, 20, 23, 24, 25)) {
  print(c("Running MSIG.Projection.Plots... on:", input.ds))

  library("scatterplot3d")
  library(MASS)

  set.seed(seed = seed, kind = NULL)

  # Read dataset

  dataset <- MSIG.Gct2Frame(filename = input.ds)
  m <- data.matrix(dataset$ds)
  gs.names <- dataset$row.names
  gs.descs <- dataset$descs
  sample.names <- dataset$names

  dim(m)

  Ns <- length(m[1, ])
  k.proj <- length(m[, 1])

  if (input.cls != "") {
    CLS <- ReadClsFile(file = input.cls)
    class.labels <- CLS$class.v
    class.phen <- CLS$phen
  } else {
    class.labels <- rep(1, Ns)
    class.phen <- "Samples"
  }

  # Separate data into train and test pieces

  if (model.set == "ALL") {
    model.set <- seq(1, Ns)
  }

  m.train <- as.matrix(m[, model.set])
  num.samples.train <- length(model.set)
  sample.names.train <- sample.names[model.set]
  if (input.cls != "") {
    class.labels.train <- class.labels[model.set]
  }
  m.test <- as.matrix(m[, -model.set])
  sample.names.test <- sample.names[-model.set]
  if (input.cls != "") {
    class.labels.test <- class.labels[-model.set]
  }

  pca <- prcomp(t(m.train), retx = TRUE, center = TRUE, scale. = TRUE)

  S1 <- pca$x[, 1]
  S2 <- pca$x[, 2]
  S3 <- pca$x[, 3]
  X1 <- pca$rotation[, 1]
  X2 <- pca$rotation[, 2]
  X3 <- pca$rotation[, 3]

  # 2D plots

  max.S <- max(sqrt(S1 * S1 + S2 * S2))
  max.X <- max(sqrt(X1 * X1 + X2 * X2))
  X1 <- max.S * X1 / max.X
  X2 <- max.S * X2 / max.X
  max.A <- max(max.S, max.X)


  c0 <- col
  c1 <- col
  #   c1 <- colors()[match(c0, colors())]
  color <- c1[class.labels]

  if (non.interactive.run == F) {
    if (.Platform$OS.type == "windows") {
      plot.filename <- output.2D.proj.plot
      x11(height = 20, width = 30)
    } else if (.Platform$OS.type == "unix") {
      plot.filename <- paste(output.2D.proj.plot, ".pdf", sep = "", collapse = "")
      pdf(file = plot.filename, height = 20, width = 30)
    }
  } else {
    if (.Platform$OS.type == "unix") {
      plot.filename <- paste(output.2D.proj.plot, ".pdf", sep = "", collapse = "")
      pdf(file = plot.filename, height = 20, width = 30)
    } else if (.Platform$OS.type == "windows") {
      plot.filename <- paste(output.2D.proj.plot, ".pdf", sep = "", collapse = "")
      pdf(file = plot.filename, height = 20, width = 30)
    }
  }

  nf <- layout(matrix(c(1, 2, 3), 1, 3, byrow = T), widths = c(3, 3, 1), heights = 1, respect = FALSE)

  # 1st subplot

  plot(S1, S2, xlim = c(-max.A, max.A), ylim = c(-max.A, max.A), type = "n", main = paste(title, " -- Model Samples Biplot", sep = ""), sub = input.ds)

  for (j in 1:num.samples.train) {
    if (min(class.labels) == 0) {
      color.code <- c1[class.labels[j] + 1]
      symb <- symbs[class.labels[j] + 1]
    } else {
      color.code <- c1[class.labels[j]]
      symb <- symbs[class.labels[j]]
    }
    #         points(S1[j], S2[j], pch=22, type="p", cex = symbol.scaling*3, bg = color.code, col = "black")
    points(S1[j], S2[j], pch = symb, type = "p", cex = symbol.scaling * 3, bg = color.code, col = "black")
  }

  if (use.biplot == TRUE) {
    for (j in 1:k.proj) {
      x.coor <- X1[j] * 0.925
      y.coor <- X2[j] * 0.925
      arrows(0, 0, x.coor, y.coor, lwd = 1, length = 0.15, angle = 20, col = "black")

      if (use.feature.names == FALSE) {
        leg.txt <- paste("F", j, sep = "")
      } else {
        leg.txt <- gs.names[j]
      }
      text(X1[j], X2[j], labels = leg.txt, adj = NULL, pos = NULL, offset = 0.5, vfont = NULL, cex = symbol.scaling * 2, col = "black")
    }

    ang <- vector(length = k.proj, mode = "numeric")
    for (j in 1:k.proj) {
      ang[j] <- ifelse(atan2(X2[j], X1[j]) > 0, atan2(X2[j], X1[j]), 2 * pi + atan2(X2[j], X1[j]))
    }

    ang.index <- order(ang, decreasing = F)
    ang2 <- ang[ang.index]

    for (j in 1:k.proj) {
      if (j == k.proj) {
        angle.in.between <- (ang2[1] - ang2[j] - 2 * pi) / 2 + ang2[j] - 2 * pi
      } else {
        angle.in.between <- (ang2[j + 1] - ang2[j]) / 2 + ang2[j]
      }
      x <- max.S * cos(angle.in.between)
      y <- max.S * sin(angle.in.between)
      arrows(0, 0, x, y, lwd = 4, length = 0, lty = 3, col = "grey50")
    }
  }

  # 2nd subplot Project non-model (test) data

  test.scores <- predict(pca, t(m.test))

  S1 <- c(pca$x[, 1], test.scores[, 1])
  S2 <- c(pca$x[, 2], test.scores[, 2])
  S3 <- c(pca$x[, 3], test.scores[, 3])

  max.S <- max(sqrt(S1 * S1 + S2 * S2))
  max.X <- max(sqrt(X1 * X1 + X2 * X2))
  X1 <- max.S * X1 / max.X
  X2 <- max.S * X2 / max.X
  num.samples <- length(S1)

  plot(S1, S2, xlim = c(-max.A, max.A), ylim = c(-max.A, max.A), type = "n", main = paste(title, " -- All Samples Biplot", sep = ""), sub = input.ds)

  for (j in 1:num.samples) {
    if (min(class.labels) == 0) {
      symb <- symbs[class.labels[j] + 1]
      color.code <- c1[class.labels[j] + 1]
    } else {
      symb <- symbs[class.labels[j]]
      color.code <- c1[class.labels[j]]
    }
    points(S1[j], S2[j], pch = symb, type = "p", cex = symbol.scaling * 3, bg = color.code, col = "black")
    #      if (j <= num.samples.train) {
    #         points(S1[j], S2[j], pch=22, type="p", cex = symbol.scaling*3, bg = color.code, col = "black")
    #       } else {
    #         points(S1[j], S2[j], pch=21, type="p", cex = symbol.scaling*2.5, bg = color.code, col = "black")
    #       }
  }

  if (use.biplot == TRUE) {
    for (j in 1:k.proj) {
      x.coor <- X1[j] * 0.925
      y.coor <- X2[j] * 0.925
      arrows(0, 0, x.coor, y.coor, lwd = 1, length = 0.15, angle = 20, col = "black")
      if (use.feature.names == FALSE) {
        leg.txt <- paste("F", j, sep = "")
      } else {
        leg.txt <- gs.names[j]
      }
      text(X1[j], X2[j], labels = leg.txt, adj = NULL, pos = NULL, offset = 0.5, vfont = NULL, cex = symbol.scaling * 2, col = "black")
    }

    ang <- vector(length = k.proj, mode = "numeric")
    for (j in 1:k.proj) {
      ang[j] <- ifelse(atan2(X2[j], X1[j]) > 0, atan2(X2[j], X1[j]), 2 * pi + atan2(X2[j], X1[j]))
    }

    ang.index <- order(ang, decreasing = F)
    ang2 <- ang[ang.index]

    for (j in 1:k.proj) {
      if (j == k.proj) {
        angle.in.between <- (ang2[1] - ang2[j] - 2 * pi) / 2 + ang2[j] - 2 * pi
      } else {
        angle.in.between <- (ang2[j + 1] - ang2[j]) / 2 + ang2[j]
      }
      x <- max.S * cos(angle.in.between)
      y <- max.S * sin(angle.in.between)
      arrows(0, 0, x, y, lwd = 4, length = 0, lty = 3, col = "grey50")
    }
  }

  # 3nd subplot: legend

  class.phen.train <- unique(class.labels.train)

  leg.txt <- class.phen
  n.phen <- length(class.phen)
  #   p.vec <- c(rep(22, length(class.phen.train)), rep(21, n.phen - length(class.phen.train)))
  p.vec <- symbs[1:n.phen]
  #   p.vec <- rep(21, n.phen)
  c.vec <- c1[1:n.phen]

  par(mar = c(0, 0, 0, 0))
  plot(c(0, 0), c(1, 1), xlim = c(0, 1), ylim = c(0, 1), axes = F, type = "n", xlab = "", ylab = "")
  legend(x = 0, y = 1, legend = leg.txt, bty = "n", xjust = 0, yjust = 1, pch = p.vec, pt.bg = c.vec, col = "black", cex = symbol.scaling * 1.35, pt.cex = symbol.scaling * 3)

  if (non.interactive.run == F) {
    if (.Platform$OS.type == "windows") {
      savePlot(filename = plot.filename, type = "jpeg", device = dev.cur())
    } else if (.Platform$OS.type == "unix") {
      dev.off()
    }
  } else {
    dev.off()
  }

  # save data

  # 3D plots

  # Normalize them between 0 and 1

  max.S <- max(sqrt(S1 * S1 + S2 * S2 + S3 * S3))
  max.X <- max(sqrt(X1 * X1 + X2 * X2 + X3 * X3))

  X1 <- max.S * X1 / max.X
  X2 <- max.S * X2 / max.X
  X3 <- max.S * X3 / max.X
  max.A <- max(max.S, max.X)

  color <- c1[class.labels]

  if (non.interactive.run == F) {
    if (.Platform$OS.type == "windows") {
      plot.filename <- output.3D.1.proj.plot
      x11(height = 20, width = 25)
    } else if (.Platform$OS.type == "unix") {
      plot.filename <- paste(output.3D.1.proj.plot, ".pdf", sep = "", collapse = "")
      pdf(file = plot.filename, height = 20, width = 25)
    }
  } else {
    if (.Platform$OS.type == "unix") {
      plot.filename <- paste(output.3D.1.proj.plot, ".pdf", sep = "", collapse = "")
      pdf(file = plot.filename, height = 20, width = 25)
    } else if (.Platform$OS.type == "windows") {
      plot.filename <- paste(output.3D.1.proj.plot, ".pdf", sep = "", collapse = "")
      pdf(file = plot.filename, height = 20, width = 25)
    }
  }

  # Subplot # 1

  nf <- layout(matrix(c(1, 2), 1, 2, byrow = T), widths = c(3, 1), heights = 1, respect = FALSE)

  x <- scatterplot3d(S1, S2, S3, xlab = "F1", ylab = "F2", zlab = "F3", type = "n", angle = 45, pch = 20, main = paste(title, " -- 3D Biplot", sep = ""), sub = " ", cex.symbols = symbol.scaling * 1)
  for (j in 1:num.samples) {
    if (min(class.labels) == 0) {
      color.code <- c1[class.labels[j] + 1]
      symb <- symbs[class.labels[j] + 1]
    } else {
      symb <- symbs[class.labels[j]]
      color.code <- c1[class.labels[j]]
    }
    cex <- 2.5 * (max(S2) - S2[j]) / (max(S2) - min(S2)) + 1.5
    x$points3d(S1[j], S2[j], S3[j], col = "black", pch = symb, bg = color.code, cex = 1 * symbol.scaling * cex)

    #      if (j <= num.samples.train) {
    #         x$points3d(S1[j], S2[j], S3[j], col="black", pch = 22, bg = color.code, cex=0.8*symbol.scaling*cex)
    #       } else {
    #         x$points3d(S1[j], S2[j], S3[j], col="black", pch = 21, bg = color.code, cex=0.8*symbol.scaling*cex)
    #       }
  }

  if (use.biplot == TRUE) {
    origin.3D <- x$xyz.convert(0, 0, 0)
    for (j in 1:k.proj) {
      x.coor <- X1[j] * 0.925
      y.coor <- X2[j] * 0.925
      z.coor <- X3[j] * 0.925
      end.point.3D <- x$xyz.convert(x.coor, y.coor, z.coor)
      arrows(origin.3D$x, origin.3D$y, end.point.3D$x, end.point.3D$y, lwd = 2, length = 0.15, angle = 20, col = "grey")

      if (use.feature.names == FALSE) {
        leg.txt <- paste("F", j, sep = "")
      } else {
        leg.txt <- gs.names[j]
      }

      x.coor <- X1[j]
      y.coor <- X2[j]
      z.coor <- X3[j]
      end.point.3D <- x$xyz.convert(x.coor, y.coor, z.coor)
      text(end.point.3D$x, end.point.3D$y, labels = leg.txt, adj = NULL, pos = NULL, offset = 0.5, vfont = NULL, cex = symbol.scaling * 2, col = "grey")
    }
  }

  leg.txt <- class.phen
  n.phen <- length(class.phen)
  #   p.vec <- rep(21, n.phen)
  #   p.vec <- c(rep(22, length(class.phen.train)), rep(21, n.phen - length(class.phen.train)))
  p.vec <- symbs[1:n.phen]
  c.vec <- c1[1:n.phen]
  par(mar = c(0, 0, 0, 0))
  plot(c(0, 0), c(1, 1), xlim = c(0, 1), ylim = c(0, 1), axes = F, type = "n", xlab = "", ylab = "")
  legend(x = 0, y = 1, legend = leg.txt, bty = "n", xjust = 0, yjust = 1, pch = p.vec, pt.bg = c.vec, col = "black", cex = symbol.scaling * 1.20, pt.cex = symbol.scaling * 1.5)

  if (non.interactive.run == F) {
    if (.Platform$OS.type == "windows") {
      savePlot(filename = plot.filename, type = "jpeg", device = dev.cur())
    } else if (.Platform$OS.type == "unix") {
      dev.off()
    }
  } else {
    dev.off()
  }

  # Subplot # 2 (reverse S2 and S3 axes)

  S3 <- -S3
  X3 <- -X3

  if (non.interactive.run == F) {
    if (.Platform$OS.type == "windows") {
      plot.filename <- output.3D.2.proj.plot
      x11(height = 20, width = 25)
    } else if (.Platform$OS.type == "unix") {
      plot.filename <- paste(output.3D.2.proj.plot, ".pdf", sep = "", collapse = "")
      pdf(file = plot.filename, height = 20, width = 25)
    }
  } else {
    if (.Platform$OS.type == "unix") {
      plot.filename <- paste(output.3D.2.proj.plot, ".pdf", sep = "", collapse = "")
      pdf(file = plot.filename, height = 20, width = 25)
    } else if (.Platform$OS.type == "windows") {
      plot.filename <- paste(output.3D.2.proj.plot, ".pdf", sep = "", collapse = "")
      pdf(file = plot.filename, height = 20, width = 25)
    }
  }

  nf <- layout(matrix(c(1, 2), 1, 2, byrow = T), widths = c(3, 1), heights = 1, respect = FALSE)

  x <- scatterplot3d(S1, S3, S2, xlab = "F1", ylab = "F3", zlab = "F2", type = "n", angle = 45, pch = 20, main = paste(title, " -- 3D Biplot", sep = ""), sub = " ", cex.symbols = symbol.scaling * 1)
  for (j in 1:num.samples) {
    if (min(class.labels) == 0) {
      symb <- symbs[class.labels[j] + 1]
      color.code <- c1[class.labels[j] + 1]
    } else {
      symb <- symbs[class.labels[j]]
      color.code <- c1[class.labels[j]]
    }
    cex <- 2.5 * (max(S3) - S3[j]) / (max(S3) - min(S3)) + 1.5
    x$points3d(S1[j], S3[j], S2[j], col = "black", pch = symb, bg = color.code, cex = 1 * symbol.scaling * cex)
    #      if (j <= num.samples.train) {
    #          x$points3d(S1[j], S3[j], S2[j], col="black", pch = 22, bg = color.code, cex=0.8*symbol.scaling*cex)
    #      } else {
    #          x$points3d(S1[j], S3[j], S2[j], col="black", pch = 21, bg = color.code, cex=0.8*symbol.scaling*cex)
    #      }
  }

  if (use.biplot == TRUE) {
    origin.3D <- x$xyz.convert(0, 0, 0)
    for (j in 1:k.proj) {
      x.coor <- X1[j] * 0.925
      y.coor <- X3[j] * 0.925
      z.coor <- X2[j] * 0.925
      end.point.3D <- x$xyz.convert(x.coor, y.coor, z.coor)
      arrows(origin.3D$x, origin.3D$y, end.point.3D$x, end.point.3D$y, lwd = 2, length = 0.15, angle = 20, col = "grey")

      if (use.feature.names == FALSE) {
        leg.txt <- paste("F", j, sep = "")
      } else {
        leg.txt <- gs.names[j]
      }

      x.coor <- X1[j]
      y.coor <- X3[j]
      z.coor <- X2[j]
      end.point.3D <- x$xyz.convert(x.coor, y.coor, z.coor)
      text(end.point.3D$x, end.point.3D$y, labels = leg.txt, adj = NULL, pos = NULL, offset = 0.5, vfont = NULL, cex = symbol.scaling * 2, col = "grey")
    }
  }

  leg.txt <- class.phen
  n.phen <- length(class.phen)
  #   p.vec <- rep(21, n.phen)
  #   p.vec <- c(rep(22, length(class.phen.train)), rep(21, n.phen - length(class.phen.train)))
  p.vec <- symbs[1:n.phen]
  c.vec <- c1[1:n.phen]
  par(mar = c(0, 0, 0, 0))
  plot(c(0, 0), c(1, 1), xlim = c(0, 1), ylim = c(0, 1), axes = F, type = "n", xlab = "", ylab = "")
  legend(x = 0, y = 1, legend = leg.txt, bty = "n", xjust = 0, yjust = 1, pch = p.vec, pt.bg = c.vec, col = "black", cex = symbol.scaling * 1.20, pt.cex = symbol.scaling * 1.5)

  if (non.interactive.run == F) {
    if (.Platform$OS.type == "windows") {
      savePlot(filename = plot.filename, type = "jpeg", device = dev.cur())
    } else if (.Platform$OS.type == "unix") {
      dev.off()
    }
  } else {
    dev.off()
  }

  # Subplot # 3 (reverse S2 and S3 axes)

  S1 <- -S1
  X1 <- -X1


  if (non.interactive.run == F) {
    if (.Platform$OS.type == "windows") {
      plot.filename <- output.3D.3.proj.plot
      x11(height = 20, width = 25)
    } else if (.Platform$OS.type == "unix") {
      plot.filename <- paste(output.3D.3.proj.plot, ".pdf", sep = "", collapse = "")
      pdf(file = plot.filename, height = 20, width = 25)
    }
  } else {
    if (.Platform$OS.type == "unix") {
      plot.filename <- paste(output.3D.3.proj.plot, ".pdf", sep = "", collapse = "")
      pdf(file = plot.filename, height = 20, width = 25)
    } else if (.Platform$OS.type == "windows") {
      plot.filename <- paste(output.3D.3.proj.plot, ".pdf", sep = "", collapse = "")
      pdf(file = plot.filename, height = 20, width = 25)
    }
  }


  nf <- layout(matrix(c(1, 2), 1, 2, byrow = T), widths = c(3, 1), heights = 1, respect = FALSE)

  x <- scatterplot3d(S2, S1, S3, xlab = "F2", ylab = "F1", zlab = "F3", type = "n", angle = 45, pch = 20, main = paste(title, " -- 3D Biplot", sep = ""), sub = " ", cex.symbols = symbol.scaling * 1)
  for (j in 1:num.samples) {
    if (min(class.labels) == 0) {
      symb <- symbs[class.labels[j] + 1]
      color.code <- c1[class.labels[j] + 1]
    } else {
      symb <- symbs[class.labels[j]]
      color.code <- c1[class.labels[j]]
    }
    cex <- 2.5 * (max(S1) - S1[j]) / (max(S1) - min(S1)) + 1.5
    x$points3d(S2[j], S1[j], S3[j], col = "black", pch = symb, bg = color.code, cex = 1 * symbol.scaling * cex)
    #      if (j <= num.samples.train) {
    #         x$points3d(S2[j], S1[j], S3[j], col="black", pch = 22, bg = color.code, cex=0.8*symbol.scaling*cex)
    #       } else {
    #         x$points3d(S2[j], S1[j], S3[j], col="black", pch = 21, bg = color.code, cex=0.8*symbol.scaling*cex)
    #       }
  }

  if (use.biplot == TRUE) {
    origin.3D <- x$xyz.convert(0, 0, 0)
    for (j in 1:k.proj) {
      x.coor <- X2[j] * 0.925
      y.coor <- X1[j] * 0.925
      z.coor <- X3[j] * 0.925
      end.point.3D <- x$xyz.convert(x.coor, y.coor, z.coor)
      arrows(origin.3D$x, origin.3D$y, end.point.3D$x, end.point.3D$y, lwd = 2, length = 0.15, angle = 20, col = "grey")

      if (use.feature.names == FALSE) {
        leg.txt <- paste("F", j, sep = "")
      } else {
        leg.txt <- gs.names[j]
      }

      x.coor <- X2[j]
      y.coor <- X1[j]
      z.coor <- X3[j]
      end.point.3D <- x$xyz.convert(x.coor, y.coor, z.coor)
      text(end.point.3D$x, end.point.3D$y, labels = leg.txt, adj = NULL, pos = NULL, offset = 0.5, vfont = NULL, cex = symbol.scaling * 2, col = "grey")
    }
  }

  # Subplot # 4 legeng

  leg.txt <- class.phen
  n.phen <- length(class.phen)
  #   p.vec <- rep(21, n.phen)
  #   p.vec <- c(rep(22, length(class.phen.train)), rep(21, n.phen - length(class.phen.train)))
  p.vec <- symbs[1:n.phen]
  c.vec <- c1[1:n.phen]
  par(mar = c(0, 0, 0, 0))
  plot(c(0, 0), c(1, 1), xlim = c(0, 1), ylim = c(0, 1), axes = F, type = "n", xlab = "", ylab = "")
  legend(x = 0, y = 1, legend = leg.txt, bty = "n", xjust = 0, yjust = 1, pch = p.vec, pt.bg = c.vec, col = "black", cex = symbol.scaling * 1.20, pt.cex = 1.5)

  if (non.interactive.run == F) {
    if (.Platform$OS.type == "windows") {
      savePlot(filename = plot.filename, type = "jpeg", device = dev.cur())
    } else if (.Platform$OS.type == "unix") {
      dev.off()
    }
  } else {
    dev.off()
  }

  # Heat map plot

  height <- ifelse(k.proj > 50, 20, 0.50 * k.proj + 7)
  #      height <- ifelse(k.proj >= 25, 25,k.proj*0.8 + 5)


  if (non.interactive.run == F) {
    if (.Platform$OS.type == "windows") {
      plot.filename <- output.heatmap.plot
      #           x11(height = height, width = 25)
    } else if (.Platform$OS.type == "unix") {
      plot.filename <- paste(output.heatmap.plot, ".pdf", sep = "", collapse = "")
      pdf(file = plot.filename, height = height, width = 25)
    }
  } else {
    if (.Platform$OS.type == "unix") {
      plot.filename <- paste(output.heatmap.plot, ".pdf", sep = "", collapse = "")
      pdf(file = plot.filename, height = height, width = 25)
    } else if (.Platform$OS.type == "windows") {
      plot.filename <- paste(output.heatmap.plot, ".pdf", sep = "", collapse = "")
      pdf(file = plot.filename, height = height, width = 25)
    }
  }
  #   nf <- layout(matrix(c(1,2), 1, 2, byrow=T), widths = c(6, 1), heights = 1, respect = FALSE)

  #   MSIG.HeatMapPlot.3(V = m, row.names = gs.names, col.labels = class.labels, col.classes = class.phen, phen.cmap = c1[1:n.phen], col.names = sample.names, main= paste(title, " -- Heat Map ", sep=""), sub = " ", xlab=" ", ylab=" ", row.norm = heatmap.row.norm,  cmap.type = heatmap.cmap.type, char.rescale = 0.8)

  #   MSIG.HeatMapPlot.4(V = m, row.names = gs.names, col.labels = class.labels, col.classes = class.phen, phen.cmap = c1[1:n.phen], col.names = sample.names, main= paste(title, " -- Heat Map ", sep=""), sub = " ", xlab=" ", ylab=" ", row.norm = heatmap.row.norm,  cmap.type = heatmap.cmap.type, char.rescale = 1)

  MSIG.HeatMapPlot.5(V = m, row.names = gs.names, col.labels = class.labels, col.classes = class.phen, col.names = sample.names, main = paste(title, " -- Heat Map ", sep = ""), sub = " ", xlab = " ", ylab = " ", row.norm = heatmap.row.norm, cmap.type = heatmap.cmap.type, char.rescale = 1)

  # legend

  #   leg.txt <- class.phen
  #   n.phen <- length(class.phen)
  #   p.vec <- c(rep(22, length(class.phen.train)), rep(21, n.phen - length(class.phen.train)))
  #   c.vec <- c1[1:n.phen]
  #   par(mar = c(0, 0, 0, 0))
  #   plot(c(0,0), c(1, 1), xlim = c(0, 1), ylim = c(0, 1), axes=F, type="n", xlab = "", ylab="")
  #   legend(x=0, y=1, legend=leg.txt, bty="n", xjust=0, yjust= 1, pch = p.vec, pt.bg = c.vec, col = "black", cex = symbol.scaling*1.2, pt.cex=symbol.scaling*2)

  if (non.interactive.run == F) {
    if (.Platform$OS.type == "windows") {
      savePlot(filename = plot.filename, type = "jpeg", device = dev.cur())
    } else if (.Platform$OS.type == "unix") {
      dev.off()
    }
  } else {
    dev.off()
  }

  # H. Tree plot

  dist.matrix <- dist(t(m))
  HC <- hclust(dist.matrix, method = "complete")

  if (non.interactive.run == F) {
    if (.Platform$OS.type == "windows") {
      plot.filename <- output.hclust.plot
      x11(height = 20, width = 25)
    } else if (.Platform$OS.type == "unix") {
      plot.filename <- paste(output.hclust.plot, ".pdf", sep = "", collapse = "")
      pdf(file = plot.filename, height = 20, width = 25)
    }
  } else {
    if (.Platform$OS.type == "unix") {
      plot.filename <- paste(output.hclust.plot, ".pdf", sep = "", collapse = "")
      pdf(file = plot.filename, height = 20, width = 25)
    } else if (.Platform$OS.type == "windows") {
      plot.filename <- paste(output.hclust.plot, ".pdf", sep = "", collapse = "")
      pdf(file = plot.filename, height = 20, width = 25)
    }
  }

  #   plot(HC, xlab="samples", cex = symbol.scaling*0.75, labels = class.phen[class.labels], col = "blue", main = paste(title, " -- Hierarchical Clustering", sep=""))

  HC$labels <- class.phen[class.labels]
  dhc <- as.dendrogram(HC, hang = 0.01, edge.root = T, dLeaf = 2)
  local({
    colLab <<- function(n) {
      if (is.leaf(n)) {
        a <- attributes(n)
        i <<- i + 1
        attr(n, "nodePar") <-
          c(a$nodePar, list(lab.col = mycols[i], pch = c(0, 0), col = c(0, 0), bg = c(0, 0), cex = c(0.8, 0.8), lab.font = i %% 1))
      }
      n
    }
    mycols <- col[class.labels[HC$order]]
    i <- 0
  })
  dL <- dendrapply(dhc, colLab)
  plot(dL, cex = 1, edge.root = T) ## --> colored labels!


  if (non.interactive.run == F) {
    if (.Platform$OS.type == "windows") {
      savePlot(filename = plot.filename, type = "jpeg", device = dev.cur())
    } else if (.Platform$OS.type == "unix") {
      dev.off()
    }
  } else {
    dev.off()
  }

  # Sorted heat map plot

  dist.matrix <- dist(t(m))
  HC <- hclust(dist.matrix, method = "complete")
  m2 <- m[, HC$order]
  sample.names2 <- sample.names[HC$order]
  class.labels2 <- class.labels[HC$order]

  dist.matrix <- dist(m)
  HC <- hclust(dist.matrix, method = "complete")
  m2 <- m2[HC$order, ]
  gs.names2 <- gs.names[HC$order]

  height <- ifelse(k.proj > 50, 20, 0.50 * k.proj + 7)

  if (non.interactive.run == F) {
    if (.Platform$OS.type == "windows") {
      plot.filename <- output.heatmap.sorted.plot
      #           x11(height = height, width = 25)
    } else if (.Platform$OS.type == "unix") {
      plot.filename <- paste(output.heatmap.sorted.plot, ".pdf", sep = "", collapse = "")
      pdf(file = plot.filename, height = height, width = 25)
    }
  } else {
    if (.Platform$OS.type == "unix") {
      plot.filename <- paste(output.heatmap.sorted.plot, ".pdf", sep = "", collapse = "")
      pdf(file = plot.filename, height = height, width = 25)
    } else if (.Platform$OS.type == "windows") {
      plot.filename <- paste(output.heatmap.sorted.plot, ".pdf", sep = "", collapse = "")
      pdf(file = plot.filename, height = height, width = 25)
    }
  }

  #   nf <- layout(matrix(c(1,2), 1, 2, byrow=T), widths = c(6, 1), heights = 1, respect = FALSE)

  #   MSIG.HeatMapPlot.4(V = m2, row.names = gs.names2, col.labels = class.labels2, col.classes = class.phen, phen.cmap = c1[1:n.phen], col.names = sample.names2, main= paste(title, " -- Heat Map (sorted)", sep=""), sub = " ", xlab=" ", ylab=" ", row.norm = heatmap.row.norm,  cmap.type = heatmap.cmap.type, char.rescale = 1)

  MSIG.HeatMapPlot.5(V = m2, row.names = gs.names2, col.labels = class.labels2, col.classes = class.phen, col.names = sample.names2, main = paste(title, " -- Heat Map (sorted)", sep = ""), sub = " ", xlab = " ", ylab = " ", row.norm = heatmap.row.norm, cmap.type = heatmap.cmap.type, char.rescale = 1)

  # legeng

  #   leg.txt <- class.phen
  #   n.phen <- length(class.phen)
  #   p.vec <- c(rep(22, length(class.phen.train)), rep(21, n.phen - length(class.phen.train)))
  #  c.vec <- c1[1:n.phen]
  #   par(mar = c(0, 0, 0, 0))
  #   plot(c(0,0), c(1, 1), xlim = c(0, 1), ylim = c(0, 1), axes=F, type="n", xlab = "", ylab="")
  #   legend(x=0, y=1, legend=leg.txt, bty="n", xjust=0, yjust= 1, pch = p.vec, pt.bg = c.vec, col = "black", cex = symbol.scaling*1.2, pt.cex=2)

  if (non.interactive.run == F) {
    if (.Platform$OS.type == "windows") {
      savePlot(filename = plot.filename, type = "jpeg", device = dev.cur())
    } else if (.Platform$OS.type == "unix") {
      dev.off()
    }
  } else {
    dev.off()
  }

  # heat map sorted inside each phenotype

  #   HC <- hclust(dist.matrix, method="complete")
  #   m2 <- m[HC$order, ]
  #  gs.names2 <- gs.names[HC$order]
  m2 <- m
  gs.names2 <- gs.names
  sample.names2 <- sample.names

  max.classes <- max(class.labels)

  for (k in 1:max.classes) {
    if (sum(class.labels == k) > 1) {
      m3 <- m2[, class.labels == k]
      sn <- sample.names2[class.labels == k]
      dist.matrix <- dist(t(m3))
      HC <- hclust(dist.matrix, method = "complete")
      m3 <- m3[, HC$order]
      sn <- sn[HC$order]
      m2[, class.labels == k] <- m3
      sample.names2[class.labels == k] <- sn
    }
  }

  height <- ifelse(k.proj > 50, 20, 0.50 * k.proj + 7)

  if (non.interactive.run == F) {
    if (.Platform$OS.type == "windows") {
      plot.filename <- output.heatmap.sorted.2.plot
      #           x11(height = height, width = 25)
    } else if (.Platform$OS.type == "unix") {
      plot.filename <- paste(output.heatmap.sorted.2.plot, ".pdf", sep = "", collapse = "")
      pdf(file = plot.filename, height = height, width = 25)
    }
  } else {
    if (.Platform$OS.type == "unix") {
      plot.filename <- paste(output.heatmap.sorted.2.plot, ".pdf", sep = "", collapse = "")
      pdf(file = plot.filename, height = height, width = 25)
    } else if (.Platform$OS.type == "windows") {
      plot.filename <- paste(output.heatmap.sorted.2.plot, ".pdf", sep = "", collapse = "")
      pdf(file = plot.filename, height = height, width = 25)
    }
  }

  #  nf <- layout(matrix(c(1,2), 1, 2, byrow=T), widths = c(6, 1), heights = 1, respect = FALSE)

  #   MSIG.HeatMapPlot.4(V = m2, row.names = gs.names2, col.labels = class.labels, col.classes = class.phen, phen.cmap = c1[1:n.phen], col.names = sample.names2, main= paste(title, " -- Heat Map (sorted inside phenotype class)", sep=""), sub = " ", xlab=" ", ylab=" ", row.norm = heatmap.row.norm,  cmap.type = heatmap.cmap.type, char.rescale = 1)

  MSIG.HeatMapPlot.5(V = m2, row.names = gs.names2, col.labels = class.labels, col.classes = class.phen, col.names = sample.names2, main = paste(title, " -- Heat Map (sorted inside phenotype class)", sep = ""), sub = " ", xlab = " ", ylab = " ", row.norm = heatmap.row.norm, cmap.type = heatmap.cmap.type, char.rescale = 1)

  # legeng

  #   leg.txt <- class.phen
  #   n.phen <- length(class.phen)
  #   p.vec <- c(rep(22, length(class.phen.train)), rep(21, n.phen - length(class.phen.train)))
  #   c.vec <- c1[1:n.phen]
  #   par(mar = c(0, 0, 0, 0))
  #  plot(c(0,0), c(1, 1), xlim = c(0, 1), ylim = c(0, 1), axes=F, type="n", xlab = "", ylab="")
  #   legend(x=0, y=1, legend=leg.txt, bty="n", xjust=0, yjust= 1, pch = p.vec, pt.bg = c.vec, col = "black", cex = symbol.scaling*1.2, pt.cex=2)

  if (non.interactive.run == F) {
    if (.Platform$OS.type == "windows") {
      savePlot(filename = plot.filename, type = "jpeg", device = dev.cur())
    } else if (.Platform$OS.type == "unix") {
      dev.off()
    }
  } else {
    dev.off()
  }
}

MSIG.HeatMapPlot.3 <- function(
    V,
    row.names = "NA",
    col.labels = "NA",
    col.classes = "NA",
    phen.cmap = "NA",
    col.names = "NA",
    main = " ",
    sub = " ",
    xlab = " ",
    ylab = " ",
    row.norm = TRUE,
    char.rescale = 1.0,
    cmap.type = 1, # 1 = vintage pinkogram, 2 = scale of blues, 3 = high-resolution pinkogram for probabilities [0, 1], 4 = high-resolution pinkogram for general values, 5 = color map for normalized enrichment scores, 6 = GenePattern color map
    max.v = "NA") {
  #
  # Plots a heatmap "pinkogram" of a gene expression matrix including phenotype vector and gene, sample and phenotype labels
  #
  # The Broad Institute
  # SOFTWARE COPYRIGHT NOTICE AGREEMENT
  # This software and its documentation are copyright 2003 by the
  # Broad Institute/Massachusetts Institute of Technology.
  # All rights are reserved.
  #
  # This software is supplied without any warranty or guaranteed support
  # whatsoever. Neither the Broad Institute nor MIT can be responsible for
  # its use, misuse, or functionality.

  n.rows <- length(V[, 1])
  n.cols <- length(V[1, ])

  if ((cmap.type == 3) | (cmap.type == 5)) {
    row.norm <- F
  }

  if (row.norm == TRUE) {
    row.mean <- apply(V, MARGIN = 1, FUN = mean)
    row.sd <- apply(V, MARGIN = 1, FUN = sd)
    row.n <- length(V[, 1])
    for (i in 1:n.rows) {
      if (row.sd[i] == 0) {
        V[i, ] <- 0
      } else {
        V[i, ] <- (V[i, ] - row.mean[i]) / (0.333 * row.sd[i])
      }
      V[i, ] <- ifelse(V[i, ] < -6, -6, V[i, ])
      V[i, ] <- ifelse(V[i, ] > 6, 6, V[i, ])
    }
  }

  if (cmap.type == 1) {
    mycol <- c(
      "#0000FF", "#4040FF", "#7070FF", "#8888FF", "#A9A9FF", "#D5D5FF", "#EEE5EE", "#FFAADA", "#FF9DB0", "#FF7080",
      "#FF5A5A", "#FF4040", "#FF0D1D"
    ) # blue-pinkogram colors. This is the 1998-vintage, pre-gene cluster, original pinkogram color map
  } else if (cmap.type == 2) {
    mycol <- c(
      "#FCFBFD", "#F4F2F8", "#F8F7FB", "#EFEDF5", "#E1E1EF", "#E8E7F2", "#DADAEB", "#C6C7E1", "#D0D1E6", "#BCBDDC", "#A8A6CF",
      "#B2B2D6", "#9E9AC8", "#8A87BF", "#9491C4", "#807DBA", "#7260AB", "#796FB3", "#6A51A3", "#5C3596", "#63439D", "#54278F", "#460D83", "#4D1A89", "#3F007D"
    )
  } else if ((cmap.type == 3) | (cmap.type == 4) | (cmap.type == 5)) {
    mycol <- vector(length = 512, mode = "numeric")

    for (k in 1:256) {
      mycol[k] <- rgb(255, k - 1, k - 1, maxColorValue = 255)
    }
    for (k in 257:512) {
      mycol[k] <- rgb(511 - (k - 1), 511 - (k - 1), 255, maxColorValue = 255)
    }
    mycol <- rev(mycol)
  } else if (cmap.type == 6) {
    mycol <- c("#4500AD", "#2700D1", "#6B58EF", "#8888FF", "#C7C1FF", "#D5D5FF", "#FFC0E5", "#FF8989", "#FF7080", "#FF5A5A", "#EF4040", "#D60C00") # blue-pinkogram colors. This is the GenePattern ComparativeMarker Selection pinkogram color map
  }

  ncolors <- length(mycol)
  mycol <- c(mycol, phen.cmap[1:length(col.classes)])

  if (cmap.type == 5) {
    if (max.v == "NA") {
      max.v <- max(max(V), -min(V))
    }
    V <- ceiling(ncolors * (V - (-max.v)) / (1.001 * (max.v - (-max.v))))
  } else {
    V <- ceiling(ncolors * (V - min(V)) / (1.001 * (max(V) - min(V))))
  }

  heatm <- matrix(0, nrow = n.rows + 1, ncol = n.cols)
  heatm[1:n.rows, ] <- V[seq(n.rows, 1, -1), ]
  heatm[n.rows + 1, ] <- ncolors + col.labels

  if (cmap.type == 2) {
    par(mar = c(3, 7, 3, 1))
  } else {
    par(mar = c(4, 15, 4, 1))
  }

  print(c("range=", range(V)))
  if (cmap.type == 5) {
    image(1:n.cols, 1:(n.rows + 1), t(heatm), zlim = c(0, ncolors + max(col.labels)), col = mycol, axes = FALSE, main = main, sub = sub, xlab = xlab, ylab = ylab)
  } else {
    image(1:n.cols, 1:(n.rows + 1), t(heatm), col = mycol, axes = FALSE, main = main, sub = sub, xlab = xlab, ylab = ylab)
  }

  if (row.names[1] != "NA") {
    numC <- nchar(row.names)
    size.row.char <- char.rescale * 30 / (n.rows + 15)
    size.col.char <- char.rescale * 30 / (n.cols + 15)
    for (i in 1:n.rows) {
      row.names[i] <- substr(row.names[i], 1, 30)
    }
    row.names <- c(row.names[seq(n.rows, 1, -1)], "Class")
    axis(2, at = 1:(n.rows + 1), labels = row.names, adj = 0.5, tick = FALSE, las = 1, cex.axis = size.row.char, font.axis = 2, line = -1)
  }

  if (col.names[1] != "NA") {
    axis(1, at = 1:n.cols, labels = col.names, tick = FALSE, las = 3, cex.axis = size.col.char, font.axis = 2, line = -1)
  }

  return()
}


MSIG.HeatMapPlot.4 <- function(
    V,
    row.names = "NA",
    row.names2 = "NA",
    col.labels = "NA",
    col.labels2 = "NA",
    col.classes = "NA",
    phen.cmap = "NA",
    col.names = "NA",
    main = " ",
    sub = " ",
    xlab = " ",
    ylab = " ",
    row.norm = TRUE,
    char.rescale = 1.0,
    cmap.type = 1, # 1 = vintage pinkogram, 2 = scale of blues, 3 = high-resolution pinkogram for probabilities [0, 1], 4 = high-resolution pinkogram for general values, 5 = color map for normalized enrichment scores, 6 = scale of red purples, 7 = scale of Oranges, 8 = scale of Greens, 9 = scale of Blues
    max.v = "NA") {
  #
  # Plots a heatmap "pinkogram" of a gene expression matrix including phenotype vector and gene, sample and phenotype labels
  #
  # The Broad Institute
  # SOFTWARE COPYRIGHT NOTICE AGREEMENT
  # This software and its documentation are copyright 2003 by the
  # Broad Institute/Massachusetts Institute of Technology.
  # All rights are reserved.
  #
  # This software is supplied without any warranty or guaranteed support
  # whatsoever. Neither the Broad Institute nor MIT can be responsible for
  # its use, misuse, or functionality.

  n.rows <- length(V[, 1])
  n.cols <- length(V[1, ])
  V1 <- matrix(0, nrow = n.rows, ncol = n.cols)

  if ((cmap.type == 3) | (cmap.type == 5)) {
    row.norm <- F
  }

  if (row.norm == TRUE) {
    row.mean <- apply(V, MARGIN = 1, FUN = mean)
    row.sd <- apply(V, MARGIN = 1, FUN = sd)
    row.n <- length(V[, 1])
    for (i in 1:n.rows) {
      if (row.sd[i] == 0) {
        V1[i, ] <- 0
      } else {
        V1[i, ] <- (V[i, ] - row.mean[i]) / (0.333 * row.sd[i])
      }
      V1[i, ] <- ifelse(V1[i, ] < -6, -6, V1[i, ])
      V1[i, ] <- ifelse(V1[i, ] > 6, 6, V1[i, ])
    }
  } else {
    V1 <- V
  }

  if (cmap.type == 1) {
    mycol <- c(
      "#0000FF", "#4040FF", "#7070FF", "#8888FF", "#A9A9FF", "#D5D5FF", "#EEE5EE", "#FFAADA",
      "#FF9DB0", "#FF7080",
      "#FF5A5A", "#FF4040", "#FF0D1D"
    ) # blue-pinkogram colors. This is the 1998-vintage, pre-gene cluster, original pinkogram color map
  } else if (cmap.type == 2) {
    mycol <- c(
      "#FCFBFD", "#F4F2F8", "#F8F7FB", "#EFEDF5", "#E1E1EF", "#E8E7F2", "#DADAEB", "#C6C7E1", "#D0D1E6",
      "#BCBDDC", "#A8A6CF",
      "#B2B2D6", "#9E9AC8", "#8A87BF", "#9491C4", "#807DBA", "#7260AB", "#796FB3", "#6A51A3", "#5C3596",
      "#63439D", "#54278F", "#460D83", "#4D1A89", "#3F007D"
    )
  } else if (cmap.type == 6) {
    mycol <- c(
      "#FFF7F3", "#FDE0DD", "#FCC5C0", "#FA9FB5", "#F768A1", "#DD3497", "#AE017E",
      "#7A0177", "#49006A"
    )
  } else if (cmap.type == 7) {
    mycol <- c(
      "#FFF5EB", "#FEE6CE", "#FDD0A2", "#FDAE6B", "#FD8D3C", "#F16913", "#D94801",
      "#A63603", "#7F2704"
    )
  } else if (cmap.type == 8) {
    mycol <- c(
      "#F7FCF5", "#E5F5E0", "#C7E9C0", "#A1D99B", "#74C476", "#41AB5D", "#238B45",
      "#006D2C", "#00441B"
    )
  } else if (cmap.type == 9) {
    mycol <- c(
      "#F7FBFF", "#DEEBF7", "#C6DBEF", "#9ECAE1", "#6BAED6", "#4292C6", "#2171B5",
      "#08519C", "#08306B"
    )
  } else if ((cmap.type == 3) | (cmap.type == 4) | (cmap.type == 5)) {
    mycol <- vector(length = 512, mode = "numeric")

    for (k in 1:256) {
      mycol[k] <- rgb(255, k - 1, k - 1, maxColorValue = 255)
    }
    for (k in 257:512) {
      mycol[k] <- rgb(511 - (k - 1), 511 - (k - 1), 255, maxColorValue = 255)
    }
    mycol <- rev(mycol)
  }

  ncolors <- length(mycol)
  mycol <- c(mycol, phen.cmap[1:length(col.classes)])

  if (cmap.type == 5) {
    if (max.v == "NA") {
      max.v <- max(max(V1), -min(V1))
    }
    V2 <- ceiling(ncolors * (V1 - (-max.v)) / (1.001 * (max.v - (-max.v))))
  } else {
    V2 <- ceiling(ncolors * (V1 - min(V1)) / (1.001 * (max(V1) - min(V1))))
  }

  heatm <- matrix(0, nrow = n.rows + 1, ncol = n.cols)
  heatm[1:n.rows, ] <- V2[seq(n.rows, 1, -1), ]
  heatm[n.rows + 1, ] <- ncolors + col.labels

  #        height <- ifelse(n.rows >= 25, 25, n.rows*0.8 + 5)
  height <- ifelse(n.rows >= 9, 9, n.rows * 0.44 + 5)
  x11(width = 14, height = height)
  nf <- layout(matrix(c(1, 2, 3, 0), 2, 2, byrow = T), widths = c(10, 2), heights = c(5, 1), respect = FALSE)

  par(mar = c(3, 8, 3, 13))

  #       if (cmap.type == 5) {
  image(1:n.cols, 1:(n.rows + 1), t(heatm), zlim = c(0, ncolors + max(col.labels)), col = mycol, axes = FALSE, main = main, sub = sub, xlab = xlab, ylab = ylab)
  #         } else {
  #           image(1:n.cols, 1:(n.rows + 1), t(heatm), col=mycol, axes=FALSE, main=main, sub = sub, xlab= xlab, ylab=ylab)
  #         }

  # add lines to separate phenotypes or subgroups

  if (col.labels2[1] != "NA") {
    groups <- split(col.labels2, col.labels2)
    len.vec <- lapply(groups, length)
    plot.div <- c(0.51, cumsum(len.vec) + 0.5)
    for (i in plot.div) {
      lines(c(i, i), c(0, n.cols), lwd = 2, cex = 0.9, col = "black")
    }
    lines(c(0.51, n.cols + 0.49), c(0.51, 0.51), lwd = 2, cex = 0.9, col = "black")
    lines(c(0.51, n.cols + 0.49), c(n.rows + 1.48, n.rows + 1.48), lwd = 2, cex = 0.9, col = "black")
  }
  if (row.names[1] != "NA") {
    numC <- nchar(row.names)
    size.row.char <- char.rescale * 25 / (n.rows + 15)
    size.col.char <- char.rescale * 35 / (n.cols + 15)
    for (i in 1:n.rows) {
      row.names[i] <- substr(row.names[i], 1, 30)
      row.names[i] <- paste(row.names[i], " ", sep = "")
    }
    row.names <- c(row.names[seq(n.rows, 1, -1)], "Class")
    axis(2, at = 1:(n.rows + 1), labels = row.names, adj = 0.5, tick = FALSE, las = 1, cex.axis = size.row.char, font.axis = 2, line = -1)
  }

  if (row.names2[1] != "NA") {
    numC <- nchar(row.names)
    size.row.char <- char.rescale * 18 / (n.rows + 15)
    for (i in 1:n.rows) {
      row.names2[i] <- substr(row.names2[i], 1, 50)
      row.names2[i] <- paste(" ", row.names2[i], sep = "")
    }
    row.names2 <- c(row.names2[seq(n.rows, 1, -1)], "     ")
    axis(4, at = 1:(n.rows + 1), labels = row.names2, adj = 0.5, tick = FALSE, las = 1, cex.axis = size.row.char, font.axis = 2, line = -1)
  }

  if (col.names[1] != "NA") {
    axis(1, at = 1:n.cols, labels = col.names, tick = FALSE, las = 3, cex.axis = size.col.char, font.axis = 2, line = -1)
  }


  # Phenotype Legend

  #       par(mar = c(15, 5, 15, 5))
  leg.txt <- col.classes
  p.vec <- rep(22, 22)
  c.vec <- phen.cmap
  par(mar = c(0, 0, 0, 0))
  plot(c(0, 0), c(1, 1), xlim = c(0, 1), ylim = c(0, 1), axes = F, type = "n", xlab = "", ylab = "")
  legend(x = 0, y = 0.8, legend = leg.txt, bty = "n", xjust = 0, yjust = 1, pch = p.vec, pt.bg = c.vec, col = "black", cex = char.rescale * 1.2, pt.cex = char.rescale * 2)

  # Color map legend

  print(c("range V=", range(V)))
  print(c("range V1=", range(V1)))
  print(c("range V2=", range(V2)))

  #       par(mar = c(8, 2, 8, 2))
  par(mar = c(2, 8, 2, 8))
  num.v <- 20
  #        if (cmap.type == 5) {
  range.v <- range(V2)
  incr <- (range.v[1] - range.v[2]) / (num.v - 1)
  heatm.v <- matrix(rev(seq(range.v[2], range.v[1], incr)), nrow = num.v, ncol = 1)
  #          image(1:1, 1:num.v, t(heatm.v), zlim = c(0, ncolors + max(col.labels)), col=mycol, axes=FALSE,
  #                main="Color \n Legend ", sub = " ", xlab= xlab, ylab=ylab)
  image(1:num.v, 1:1, heatm.v,
    zlim = c(0, ncolors + max(col.labels)), col = mycol, axes = FALSE,
    main = " ", sub = " ", xlab = ylab, ylab = xlab
  )
  range.v <- range(V1)
  incr <- (range.v[1] - range.v[2]) / (num.v - 1)
  heatm.v2 <- matrix(signif(rev(seq(range.v[2], range.v[1], incr)), digits = 2), nrow = num.v, ncol = 1)
  print(c("heatm.v2=", heatm.v2))
  axis(3, at = 1:num.v, labels = heatm.v2, adj = 0.5, tick = FALSE, las = 1, cex.axis = 0.5 * char.rescale, font.axis = 1)

  #         } else {
  #          range.v <- range(V2)
  #          incr <-  (range.v[1] - range.v[2])/(num.v - 1)
  #          heatm.v <- matrix(rev(seq(range.v[2], range.v[1], incr)), nrow=num.v, ncol=1)
  #          image(1:1, 1:num.v, t(heatm.v), col=mycol, axes=FALSE, main="Color \n Legend ", sub = " ", xlab= xlab, ylab=ylab)
  #          range.v <- range(V1)
  #          incr <-  (range.v[1] - range.v[2])/(num.v - 1)
  #          heatm.v <- matrix(rev(signif(seq(range.v[2], range.v[1], incr), digits=3)), nrow=num.v, ncol=1)
  #          axis(2, at=1:num.v, labels=heatm.v, adj= 0.5, tick=FALSE, las = 1, cex.axis=size.row.char, font.axis=2, line=-1)
  #         }

  return()
}

MSIG.Signature.Plot <- function(
    V,
    row.names = "NA",
    col.labels = "NA",
    col.classes = "NA",
    phen.cmap = "NA",
    col.names = "NA",
    main = " ",
    sub = " ",
    xlab = " ",
    ylab = " ",
    row.norm = TRUE,
    char.rescale = 1.0,
    max.v = "NA",
    seed = 1729) {
  n.rows <- length(V[, 1])
  n.cols <- length(V[1, ])
  V1 <- matrix(0, nrow = n.rows, ncol = n.cols)

  if (row.norm == TRUE) {
    row.mean <- apply(V, MARGIN = 1, FUN = mean)
    row.sd <- apply(V, MARGIN = 1, FUN = sd)
    row.n <- length(V[, 1])
    for (i in 1:n.rows) {
      if (row.sd[i] == 0) {
        V1[i, ] <- 0
      } else {
        V1[i, ] <- (V[i, ] - row.mean[i]) / row.sd[i]
      }
      V1[i, ] <- ifelse(V1[i, ] < -4, -4, V1[i, ])
      V1[i, ] <- ifelse(V1[i, ] > 4, 4, V1[i, ])
    }
  } else {
    V1 <- V
  }
  mycol <- vector(length = 512, mode = "numeric")

  for (k in 1:256) {
    mycol[k] <- rgb(255, k - 1, k - 1, maxColorValue = 255)
  }
  for (k in 257:512) {
    mycol[k] <- rgb(511 - (k - 1), 511 - (k - 1), 255, maxColorValue = 255)
  }
  mycol <- rev(mycol)
  ncolors <- length(mycol)
  mycol <- c(mycol, phen.cmap[1:length(col.classes)])

  V2 <- ceiling(ncolors * (V1 - min(V1)) / (1.001 * (max(V1) - min(V1))))

  heatm <- matrix(0, nrow = n.rows + 1, ncol = n.cols)
  heatm[1:n.rows, ] <- V2[seq(n.rows, 1, -1), ]

  heatm[n.rows + 1, ] <- ncolors + col.labels
  height <- ifelse(n.rows >= 25, 25, n.rows * 0.8 + 5)

  x11(width = 19, height = 11)
  nf <- layout(matrix(c(1, 2, 3), 1, 3, byrow = T), widths = c(4, 4, 1), heights = 1, respect = FALSE)

  par(mar = c(8, 9, 4, 4))
  image(1:n.cols, 1:(n.rows + 1), t(heatm),
    zlim = c(0, ncolors + max(col.labels)), col = mycol, axes = FALSE, main = main,
    sub = sub, xlab = xlab, ylab = ylab
  )

  if (row.names[1] != "NA") {
    numC <- nchar(row.names)
    size.row.char <- char.rescale * 35 / (n.rows + 15)
    size.col.char <- char.rescale * 30 / (n.cols + 15)
    for (i in 1:n.rows) {
      row.names[i] <- substr(row.names[i], 1, 30)
    }
    row.names <- c(row.names[seq(n.rows, 1, -1)], "Class")
    axis(2,
      at = 1:(n.rows + 1), labels = row.names, adj = 0.5, tick = FALSE, las = 1, cex.axis = size.row.char,
      font.axis = 2, line = -1
    )
  }
  if (col.names[1] != "NA") {
    axis(1, at = 1:n.cols, labels = col.names, tick = FALSE, las = 3, cex.axis = size.col.char, font.axis = 2, line = -1)
  }

  par(mar = c(8, 2, 4, 2))

  library(RColorBrewer)
  col.list <- c(
    brewer.pal(n = 9, name = "Oranges"), brewer.pal(n = 9, name = "Reds"), brewer.pal(n = 9, name = "Blues"),
    brewer.pal(n = 9, name = "Greens"), brewer.pal(n = 9, name = "Greys"), brewer.pal(n = 9, name = "Purples"),
    brewer.pal(n = 9, name = "Pastel1"), brewer.pal(n = 9, name = "Set1"), brewer.pal(n = 8, name = "Dark2"),
    brewer.pal(n = 8, name = "Accent"), brewer.pal(n = 10, name = "Spectral"), brewer.pal(n = 9, name = "Set1"),
    brewer.pal(n = 8, name = "Set2"), brewer.pal(n = 8, name = "Set3"), brewer.pal(n = 8, name = "BuGn"),
    brewer.pal(n = 8, name = "GnBu"), brewer.pal(n = 8, name = "OrRd"), brewer.pal(n = 8, name = "PuBu"),
    brewer.pal(n = 8, name = "PuBuGn"), brewer.pal(n = 8, name = "RdPu"), brewer.pal(n = 8, name = "YlGn"),
    brewer.pal(n = 8, name = "YlOrBr"), brewer.pal(n = 8, name = "YlOrRd")
  )
  set.seed(seed)
  sig.col <- sample(col.list, size = 1000, replace = T)

  h.range <- range(apply(V2, MARGIN = 2, FUN = sum))
  h.size <- h.range[2] / 30

  V3 <- rbind(rep(h.size, n.cols), apply(V2, MARGIN = 2, FUN = rev))
  barplot(V3,
    main = main, font.axis = 1.25, cex.lab = 1.5, cex.axis = 1.25, cex.names = 1.25, width = 1, space = 0,
    border = 1, col = sig.col
  )

  barplot(rep(h.size, n.cols),
    cex.lab = 1.5, cex.axis = 1.25, cex.names = 1.25, width = 1, space = 0,
    border = 0, col = phen.cmap[col.labels], add = T
  )

  leg.txt <- col.classes
  p.vec <- rep(21, 21)
  c.vec <- phen.cmap
  par(mar = c(0, 0, 0, 0))
  plot(c(0, 0), c(1, 1), xlim = c(0, 1), ylim = c(0, 1), axes = F, type = "n", xlab = "", ylab = "")
  legend(
    x = 0, y = 0.8, legend = leg.txt, bty = "n", xjust = 0, yjust = 1, pch = p.vec, pt.bg = c.vec, col = "black",
    cex = char.rescale * 1.3, pt.cex = char.rescale * 2.5
  )

  return()
}


MSIG.Score.Plot <- function(
    z,
    main = "",
    phen.cmap,
    char.rescale = 1,
    col.classes,
    col.labels,
    create.legend = T,
    create.window = T,
    xlab = " ",
    ylab = " ") {
  size <- length(z)

  if (create.window == T) {
    x11(width = 19, height = 11)
  }
  if ((create.window == T && create.legend == T)) {
    nf <- layout(matrix(c(1, 2), 1, 2, byrow = T), widths = c(6, 1), heights = 1, respect = FALSE)
  }
  barplot(z,
    xlab = xlab, ylab = ylab, main = main, font.axis = 1, cex.lab = 1, cex.axis = 1, cex.names = 1,
    width = 1, space = 0, col = phen.cmap[col.labels]
  )
  if (create.legend == T) {
    leg.txt <- col.classes
    p.vec <- rep(21, 21)
    c.vec <- phen.cmap
    par(mar = c(0, 0, 0, 0))
    plot(c(0, 0), c(1, 1), xlim = c(0, 1), ylim = c(0, 1), axes = F, type = "n", xlab = "", ylab = "")
    legend(
      x = 0, y = 0.8, legend = leg.txt, bty = "n", xjust = 0, yjust = 1, pch = p.vec, pt.bg = c.vec, col = "black",
      cex = char.rescale, pt.cex = char.rescale * 2
    )
  }
}

MSIG.SplitTrainTest <- function(V, class.v = 0, fraction = 0.50) {
  cols <- length(V[1, ])
  rows <- length(V[, 1])
  train.size <- ceiling(cols * fraction)
  test.size <- cols - train.size
  index <- sample(1:cols)
  train <- matrix(0, nrow = rows, ncol = train.size)
  test <- matrix(0, nrow = rows, ncol = test.size)
  for (i in 1:train.size) {
    train[, i] <- V[, index[i]]
  }
  for (i in 1:test.size) {
    test[, i] <- V[, index[train.size + i]]
  }
  if (length(class.v) > 0) { # there is a class vector
    cls.train <- vector(mode = typeof(class.v), length = train.size)
    cls.test <- vector(mode = typeof(class.v), length = test.size)
    for (i in 1:train.size) {
      cls.train[i] <- class.v[index[i]]
    }
    for (i in 1:test.size) {
      cls.test[i] <- class.v[index[train.size + i]]
    }
    return(list(train = train, test = test, cls.train = cls.train, cls.test = cls.test))
  } else {
    return(list(train = train, test = test))
  }
}


MPM <- function(
    model.dataset.table,
    test.datasets.table,
    norm,
    identifier,
    k.proj,
    alg = "NMF.div",
    niter = 1500,
    seed = 123,
    nchar.phen = 2,
    postprojnorm = TRUE,
    use.biplot = TRUE,
    non.interactive.run = FALSE,
    heatmap.row.norm = FALSE,
    heatmap.cmap.type = 3,
    use.feature.names = FALSE,
    high.conf.thres = 0.5,
    output.dir,
    col = c("green", "blue", "pink", "red", "orange", "red4", "steelblue2", "violet"),
    symbs = c(22, 21, 20, 23, 24, 25, 21, 20, 23, 24, 25),
    symbol.scaling = 1,
    levels = NULL,
    nlevels = 20,
    kernel = "radial",
    cost = 1,
    gamma = 0.05,
    theta = 0,
    model.set.refinement = T) {
  # -------------------------------------------------- Projection Methodology

  O <- MSIG.Subset.Dataset(
    input.ds = model.dataset.table$gct.file,
    input.cls = model.dataset.table$cls.file,
    column.subset = model.dataset.table$column.subset,
    column.sel.type = model.dataset.table$column.sel.type,
    row.subset = "ALL",
    output.ds = paste(output.dir, "temp1.gct", sep = ""),
    output.cls = paste(output.dir, identifier, ".model_set.1.cls", sep = "")
  )

  O <- MSIG.Preprocess.Dataset(
    input.ds = paste(output.dir, "temp1.gct", sep = ""),
    output.ds = paste(output.dir, identifier, ".model_set.1.gct", sep = ""),
    thres = model.dataset.table$thres,
    ceil = model.dataset.table$ceil,
    fold = model.dataset.table$fold,
    delta = model.dataset.table$delta,
    normalization = norm
  )

  O <- MSIG.Extract.Factors(
    input.ds = paste(output.dir, identifier, ".model_set.1.gct", sep = ""),
    input.cls = paste(output.dir, identifier, ".model_set.1.cls", sep = ""),
    output.W.file = paste(output.dir, identifier, ".model.W.gct", sep = ""),
    output.H.file = paste(output.dir, identifier, ".model.H.gct", sep = ""),
    k.proj = k.proj,
    alg = alg,
    niter = niter,
    seed = seed,
    theta = theta
  )

  O <- MSIG.Factors.Project(
    input.ds = paste(output.dir, identifier, ".model_set.1.gct", sep = ""),
    factors.ds = paste(output.dir, identifier, ".model.W.gct", sep = ""),
    postprojnorm = postprojnorm,
    output.file = paste(output.dir, identifier, ".model_dataset.H.gct", sep = "")
  )

  CLS <- MSIG.ReadClsFile(file = paste(output.dir, identifier, ".model_set.1.cls", sep = ""))
  model.size <- length(CLS$class.v)

  O <- MSIG.Evaluate.Projection.2(
    input.ds = paste(output.dir, identifier, ".model_dataset.H.gct", sep = ""),
    input.cls = paste(output.dir, identifier, ".model_set.1.cls", sep = ""),
    model.set = seq(1, model.size),
    prediction.results.file = paste(output.dir, identifier, ".prelim.pred.txt", sep = ""),
    prediction.matrix.file = paste(output.dir, identifier, ".prelim.pred.gct", sep = ""),
    train.pred.plot = paste(output.dir, identifier, ".prelim.train.pred", sep = ""),
    test.pred.plot = paste(output.dir, identifier, ".prelim.test.pred", sep = ""),
    pred.2D.plot = paste(output.dir, identifier, ".prelim.2D.pred", sep = ""),
    col = col,
    use.feature.names = use.feature.names,
    nchar.phen = nchar.phen,
    high.conf.thres = high.conf.thres,
    symbs = symbs,
    symbol.scaling = symbol.scaling,
    levels = levels,
    nlevels = nlevels,
    kernel = kernel,
    cost = cost,
    gamma = gamma
  )


  input.txt <- paste(output.dir, identifier, ".prelim.pred.txt", sep = "")
  pred.table <- read.table(file = input.txt, skip = 2, nrow = model.size, sep = "\t", header = T, comment.char = "", as.is = T)
  conf.list <- pred.table["Conf..H.L."]
  actual.list <- pred.table["Actual"]
  predicted <- pred.table["Predicted"]
  sample.select <- ((conf.list == " H ") & (actual.list == predicted))

  if (model.set.refinement == T) {
    high.conf.set <- seq(1, model.size)[sample.select]
  } else {
    high.conf.set <- seq(1, model.size)
  }

  #  high.conf.set <- seq(1, model.size)[((conf.list == " H ") & (error.list != " * "))]
  #     high.conf.set <- seq(1, model.size)[conf.list == " H "]
  print(c("pred table from file=", pred.table))
  print(c("model sizet=", model.size))
  print(c("high.conf.set=", high.conf.set))
  print(c("sample.select", sample.select))
  print(paste("Original: ", model.size, " (samples); New: ", length(high.conf.set), " (samples); diff: ", model.size - length(high.conf.set), sep = " "))

  O <- MSIG.Subset.Dataset(
    input.ds = paste(output.dir, identifier, ".model_set.1.gct", sep = ""),
    input.cls = paste(output.dir, identifier, ".model_set.1.cls", sep = ""),
    column.subset = high.conf.set,
    row.subset = "ALL",
    output.ds = paste(output.dir, identifier, ".model_set.2.gct", sep = ""),
    output.cls = paste(output.dir, identifier, ".model_set.2.cls", sep = "")
  )

  O <- MSIG.Subset.Dataset(
    input.ds = model.dataset.table$gct.file,
    input.cls = model.dataset.table$cls.file,
    column.subset = high.conf.set,
    row.subset = "ALL",
    output.ds = paste(output.dir, identifier, ".model_set.0.gct", sep = ""),
    output.cls = paste(output.dir, identifier, ".model_set.0.cls", sep = "")
  )

  O <- MSIG.Extract.Factors(
    input.ds = paste(output.dir, identifier, ".model_set.2.gct", sep = ""),
    input.cls = paste(output.dir, identifier, ".model_set.2.cls", sep = ""),
    output.W.file = paste(output.dir, identifier, ".model.W.gct", sep = ""),
    output.H.file = paste(output.dir, identifier, ".model.H.gct", sep = ""),
    k.proj = k.proj,
    alg = alg,
    niter = niter,
    seed = seed,
    theta = theta
  )

  O <- MSIG.Factors.Project(
    input.ds = paste(output.dir, identifier, ".model_set.2.gct", sep = ""),
    factors.ds = paste(output.dir, identifier, ".model.W.gct", sep = ""),
    postprojnorm = postprojnorm,
    output.file = paste(output.dir, identifier, ".model_dataset.H.gct", sep = "")
  )

  O <- MSIG.Subset.Dataset(
    input.ds = paste(output.dir, identifier, ".model_dataset.H.gct", sep = ""),
    input.cls = paste(output.dir, identifier, ".model_set.2.cls", sep = ""),
    column.subset = "ALL",
    column.sel.type = "samples",
    row.subset = "ALL",
    output.ds = paste(output.dir, identifier, ".all.H.gct", sep = ""),
    output.cls = paste(output.dir, identifier, ".all.H.cls", sep = "")
  )

  O <- MSIG.Subset.Dataset(
    input.ds = paste(output.dir, identifier, ".model_set.2.gct", sep = ""),
    input.cls = paste(output.dir, identifier, ".model_set.2.cls", sep = ""),
    column.subset = "ALL",
    column.sel.type = "samples",
    row.subset = "ALL",
    output.ds = paste(output.dir, identifier, ".all.gct", sep = ""),
    output.cls = paste(output.dir, identifier, ".all.cls", sep = "")
  )

  # Pre-process and project all the test datasets

  if (!is.null(test.datasets.table)) {
    max.files <- length(test.datasets.table)
    for (ds in 1:max.files) {
      print(c("Processing test file: ", test.datasets.table[[ds]]$gct.file))

      O <- MSIG.Subset.Dataset(
        input.ds = test.datasets.table[[ds]]$gct.file,
        input.cls = test.datasets.table[[ds]]$cls.file,
        column.subset = test.datasets.table[[ds]]$column.subset,
        column.sel.type = test.datasets.table[[ds]]$column.sel.type,
        row.subset = "ALL",
        output.ds = paste(output.dir, "temp1.gct", sep = ""),
        output.cls = paste(output.dir, "temp1.cls", sep = "")
      )

      O <- MSIG.Preprocess.Dataset(
        input.ds = paste(output.dir, "temp1.gct", sep = ""),
        output.ds = paste(output.dir, "temp2.gct", sep = ""),
        thres = test.datasets.table[[ds]]$thres,
        ceil = test.datasets.table[[ds]]$ceil,
        normalization = norm
      )

      O <- MSIG.Factors.Project(
        input.ds = paste(output.dir, "temp2.gct", sep = ""),
        factors.ds = paste(output.dir, identifier, ".model.W.gct", sep = ""),
        postprojnorm = postprojnorm,
        output.file = paste(output.dir, "temp3.gct", sep = "")
      )

      O <- MSIG.Match.and.Merge(
        input1.ds = paste(output.dir, identifier, ".all.H.gct", sep = ""),
        input1.cls = paste(output.dir, identifier, ".all.H.cls", sep = ""),
        input2.ds = paste(output.dir, "temp3.gct", sep = ""),
        input2.cls = paste(output.dir, "temp1.cls", sep = ""),
        output.ds = paste(output.dir, identifier, ".all.H.gct", sep = ""),
        output.cls = paste(output.dir, identifier, ".all.H.cls", sep = "")
      )

      O <- MSIG.Match.and.Merge(
        input1.ds = paste(output.dir, identifier, ".all.gct", sep = ""),
        input1.cls = paste(output.dir, identifier, ".all.cls", sep = ""),
        input2.ds = paste(output.dir, "temp2.gct", sep = ""),
        input2.cls = paste(output.dir, "temp1.cls", sep = ""),
        output.ds = paste(output.dir, identifier, ".all.gct", sep = ""),
        output.cls = paste(output.dir, identifier, ".all.cls", sep = "")
      )

      print(c("Done processing test file: ", test.datasets.table[[ds]]$gct.file))
    }
  }
  # Plots for NMP

  CLS <- MSIG.ReadClsFile(file = paste(output.dir, identifier, ".model_set.2.cls", sep = ""))
  model.size <- length(CLS$class.v)

  O <- MSIG.Projection.Plots.3(
    input.ds = paste(output.dir, identifier, ".all.H.gct", sep = ""),
    input.cls = paste(output.dir, identifier, ".all.H.cls", sep = ""),
    model.set = seq(1, model.size),
    output.2D.proj.file = paste(output.dir, identifier, ".2D.proj.gct", sep = ""),
    output.2D.proj.plot = paste(output.dir, identifier, ".2D.proj", sep = ""),
    output.3D.proj.file = paste(output.dir, identifier, ".3D.proj.gct", sep = ""),
    output.3D.1.proj.plot = paste(output.dir, identifier, ".3D.1.proj", sep = ""),
    output.3D.2.proj.plot = paste(output.dir, identifier, ".3D.2.proj", sep = ""),
    output.3D.3.proj.plot = paste(output.dir, identifier, ".3D.3.proj", sep = ""),
    output.heatmap.plot = paste(output.dir, identifier, ".heatmap", sep = ""),
    output.heatmap.sorted.plot = paste(output.dir, identifier, ".heatmap.sorted", sep = ""),
    output.heatmap.sorted.2.plot = paste(output.dir, identifier, ".heatmap.sorted.2", sep = ""),
    output.hclust.plot = paste(output.dir, identifier, ".hclust", sep = ""),
    use.biplot = use.biplot,
    title = identifier,
    seed = seed,
    non.interactive.run = non.interactive.run,
    heatmap.row.norm = heatmap.row.norm,
    heatmap.cmap.type = heatmap.cmap.type,
    symbol.scaling = symbol.scaling,
    col = col,
    symbs = symbs
  )

  # Evaluate projection

  O <- MSIG.Evaluate.Projection.2(
    input.ds = paste(output.dir, identifier, ".all.H.gct", sep = ""),
    input.cls = paste(output.dir, identifier, ".all.H.cls", sep = ""),
    model.set = seq(1, model.size),
    prediction.results.file = paste(output.dir, identifier, ".pred.txt", sep = ""),
    prediction.matrix.file = paste(output.dir, identifier, ".pred.gct", sep = ""),
    train.pred.plot = paste(output.dir, identifier, ".train.pred", sep = ""),
    test.pred.plot = paste(output.dir, identifier, ".test.pred", sep = ""),
    pred.2D.plot = paste(output.dir, identifier, ".2D.pred", sep = ""),
    col = col,
    use.feature.names = use.feature.names,
    nchar.phen = nchar.phen,
    high.conf.thres = high.conf.thres,
    symbol.scaling = symbol.scaling,
    symbs = symbs,
    levels = levels,
    nlevels = nlevels,
    kernel = kernel,
    cost = cost,
    gamma = gamma
  )

  # Compute hierarchical clustering

  input.ds <- paste(output.dir, identifier, ".all.gct", sep = "")
  input.cls <- paste(output.dir, identifier, ".all.cls", sep = "")
  dataset <- MSIG.Gct2Frame(filename = input.ds)
  m.ds <- data.matrix(dataset$ds)
  N.ds <- length(m.ds[, 1])
  M.ds <- length(m.ds[1, ])
  ds.names <- dataset$row.names
  ds.descs <- dataset$descs
  ds.sample.names <- dataset$names

  # Read class labels

  CLS <- MSIG.ReadClsFile(file = input.cls)
  class.labels <- CLS$class.v
  class.phen <- CLS$phen
  class.list <- CLS$class.list

  class.labels <- match(class.list, class.phen)
  class.phen <- unique(class.list)

  # Compute hierarchical tree clustering

  dist.matrix <- dist(t(m.ds))
  HC <- hclust(dist.matrix, method = "complete")

  x11(height = 24, width = 30)
  plot(HC, xlab = "samples", cex = 0.7, labels = class.list, col = "blue", main = " Hierarchical Clustering (original data)")

  plot.filename <- paste(output.dir, identifier, ".htree", sep = "")
  savePlot(filename = plot.filename, type = "jpeg", device = dev.cur())

  # Read projected dataset

  input.ds <- paste(output.dir, identifier, ".all.H.gct", sep = "")
  input.cls <- paste(output.dir, identifier, ".all.H.cls", sep = "")
  dataset <- MSIG.Gct2Frame(filename = input.ds)
  m.ds <- data.matrix(dataset$ds)
  N.ds <- length(m.ds[, 1])
  M.ds <- length(m.ds[1, ])
  ds.names <- dataset$row.names
  ds.descs <- dataset$descs
  ds.sample.names <- dataset$names

  # Compute hierarchical tree clustering

  dist.matrix <- dist(t(m.ds))
  HC <- hclust(dist.matrix, method = "complete")

  x11(height = 24, width = 30)
  plot(HC, xlab = "samples", cex = 0.7, labels = class.list, col = "blue", main = " Hierarchical Clustering (projected data)")

  plot.filename <- paste(output.dir, identifier, ".H.htree", sep = "")
  savePlot(filename = plot.filename, type = "jpeg", device = dev.cur())

  # Compute class membership

  membership <- vector(length = M.ds, mode = "numeric")

  for (j in 1:M.ds) { # Find membership
    membership[j] <- order(m.ds[, j], decreasing = T)
  }

  mem.order <- order(membership, decreasing = F)
  membership.sorted <- membership[mem.order]
  ds.sample.names <- paste(class.list, ds.sample.names, sep = "_")
  ds.sample.names.sorted <- ds.sample.names[mem.order]
  class.list.sorted <- class.list[mem.order]

  mem.table <- data.frame(cbind(class.list, ds.sample.names, membership, rep(" ", M.ds), class.list.sorted, ds.sample.names.sorted, membership.sorted))
  row.names(mem.table) <- seq(1, M.ds)
  names(mem.table) <- c("Phen", "Sample Names", "Membership", " ", "Phen Sorted", "Sample Names Sorted", "Membership Sorted")

  mem.filename <- paste(output.dir, identifier, ".H.mem.txt", sep = "")

  write.table(file = mem.filename, mem.table, quote = F, sep = "\t")

  table(class.list.sorted, membership.sorted)
}

MSIG.Match.and.Select <- function(
    input1.ds,
    input2.ds,
    output.ds) {
  # Match the genes of the first dataset on the second and select those rows from the second

  # start of methodology

  print(c("Running MSIG.Match.and.Select... on: ", input1.ds, " ", input2.ds))

  # Read input datasets

  dataset1 <- MSIG.Gct2Frame(filename = input1.ds)
  m1 <- data.matrix(dataset1$ds)
  gs.names1 <- dataset1$row.names
  gs.descs1 <- dataset1$descs
  sample.names1 <- dataset1$names

  dataset2 <- MSIG.Gct2Frame(filename = input2.ds)
  m2 <- data.matrix(dataset2$ds)
  gs.names2 <- dataset2$row.names
  gs.descs2 <- dataset2$descs
  sample.names2 <- dataset2$names

  # Match features to first dataset and create matching m2 dataset

  gs.names3 <- intersect(gs.names1, gs.names2)

  locations2 <- match(gs.names3, gs.names2, nomatch = 0)
  gs.names2 <- gs.names2[locations2]
  gs.descs2 <- gs.descs2[locations2]
  m2 <- m2[locations2, ]

  # Save dataset

  V <- data.frame(m2)
  names(V) <- sample.names2
  row.names(V) <- gs.names2
  write.gct(gct.data.frame = V, descs = gs.descs2, filename = output.ds)
}

MPM.2 <- function(
    model.dataset.table,
    test.datasets.table,
    identifier,
    k.proj,
    alg = "NMF.div",
    niter = 1500,
    seed = 123,
    nchar.phen = 2,
    postprojnorm = TRUE,
    use.biplot = TRUE,
    non.interactive.run = FALSE,
    heatmap.row.norm = FALSE,
    heatmap.cmap.type = 3,
    use.feature.names = FALSE,
    high.conf.thres = 0.5,
    output.dir,
    col = c("green", "blue", "pink", "red", "orange", "red4", "steelblue2", "violet"),
    symbs = c(22, 21, 20, 23, 24, 25, 21, 20, 23, 24, 25),
    symbol.scaling = 1,
    levels = NULL,
    nlevels = 20,
    kernel = "radial",
    cost = 1,
    gamma = 0.05,
    theta = 0,
    model.set.refinement = T,
    produce.contours = T) {
  print(c(model.dataset.table))
  print(c(test.datasets.table))


  # -------------------------------------------------- Projection Methodology

  set.seed(seed = seed, kind = NULL)

  O <- MSIG.Subset.Dataset(
    input.ds = model.dataset.table$gct.file,
    input.cls = model.dataset.table$cls.file,
    column.subset = model.dataset.table$column.subset,
    column.sel.type = model.dataset.table$column.sel.type,
    row.subset = "ALL",
    output.ds = paste(output.dir, "temp1.gct", sep = ""),
    output.cls = paste(output.dir, identifier, ".model_set.1.cls", sep = "")
  )

  O <- MSIG.Preprocess.Dataset(
    input.ds = paste(output.dir, "temp1.gct", sep = ""),
    output.ds = paste(output.dir, identifier, ".model_set.1.gct", sep = ""),
    thres = model.dataset.table$thres,
    ceil = model.dataset.table$ceil,
    fold = model.dataset.table$fold,
    delta = model.dataset.table$delta,
    normalization = model.dataset.table$norm
  )

  O <- MSIG.Extract.Factors(
    input.ds = paste(output.dir, identifier, ".model_set.1.gct", sep = ""),
    input.cls = paste(output.dir, identifier, ".model_set.1.cls", sep = ""),
    output.W.file = paste(output.dir, identifier, ".model.W.gct", sep = ""),
    output.H.file = paste(output.dir, identifier, ".model.H.gct", sep = ""),
    k.proj = k.proj,
    alg = alg,
    niter = niter,
    seed = seed,
    theta = theta,
    sort.factors = T
  )

  O <- MSIG.Factors.Project(
    input.ds = paste(output.dir, identifier, ".model_set.1.gct", sep = ""),
    factors.ds = paste(output.dir, identifier, ".model.W.gct", sep = ""),
    postprojnorm = postprojnorm,
    output.file = paste(output.dir, identifier, ".model_dataset.H.gct", sep = "")
  )

  CLS <- MSIG.ReadClsFile(file = paste(output.dir, identifier, ".model_set.1.cls", sep = ""))
  model.size <- length(CLS$class.v)

  O <- MSIG.Evaluate.Projection.2(
    input.ds = paste(output.dir, identifier, ".model_dataset.H.gct", sep = ""),
    input.cls = paste(output.dir, identifier, ".model_set.1.cls", sep = ""),
    model.set = seq(1, model.size),
    prediction.results.file = paste(output.dir, identifier, ".prelim.pred.txt", sep = ""),
    prediction.matrix.file = paste(output.dir, identifier, ".prelim.pred.gct", sep = ""),
    train.pred.plot = paste(output.dir, identifier, ".prelim.train.pred", sep = ""),
    test.pred.plot = paste(output.dir, identifier, ".prelim.test.pred", sep = ""),
    pred.2D.plot = paste(output.dir, identifier, ".prelim.2D.pred", sep = ""),
    col = col,
    use.feature.names = use.feature.names,
    nchar.phen = nchar.phen,
    high.conf.thres = high.conf.thres,
    symbs = symbs,
    symbol.scaling = symbol.scaling,
    levels = levels,
    nlevels = nlevels,
    kernel = kernel,
    cost = cost,
    gamma = gamma,
    produce.contours = F
  )

  input.txt <- paste(output.dir, identifier, ".prelim.pred.txt", sep = "")
  pred.table <- read.table(file = input.txt, skip = 2, nrow = model.size, sep = "\t", header = T, comment.char = "", as.is = T)
  conf.list <- pred.table["Conf..H.L."]
  actual.list <- pred.table["Actual"]
  predicted <- pred.table["Predicted"]
  sample.select <- ((conf.list == " H ") & (actual.list == predicted))

  if (model.set.refinement == T) {
    high.conf.set <- seq(1, model.size)[sample.select]
  } else {
    high.conf.set <- seq(1, model.size)
  }

  print(c("pred table from file=", pred.table))
  print(c("model sizet=", model.size))
  print(c("high.conf.set=", high.conf.set))
  print(c("sample.select", sample.select))
  print(paste("Original: ", model.size, " (samples); New: ", length(high.conf.set), " (samples); diff: ", model.size - length(high.conf.set), sep = " "))

  O <- MSIG.Subset.Dataset(
    input.ds = paste(output.dir, identifier, ".model_set.1.gct", sep = ""),
    input.cls = paste(output.dir, identifier, ".model_set.1.cls", sep = ""),
    column.subset = high.conf.set,
    row.subset = "ALL",
    output.ds = paste(output.dir, identifier, ".model_set.2.gct", sep = ""),
    output.cls = paste(output.dir, identifier, ".model_set.2.cls", sep = "")
  )

  O <- MSIG.Subset.Dataset(
    input.ds = model.dataset.table$gct.file,
    input.cls = model.dataset.table$cls.file,
    column.subset = high.conf.set,
    row.subset = "ALL",
    output.ds = paste(output.dir, identifier, ".model_set.0.gct", sep = ""),
    output.cls = paste(output.dir, identifier, ".model_set.0.cls", sep = "")
  )

  O <- MSIG.Extract.Factors(
    input.ds = paste(output.dir, identifier, ".model_set.2.gct", sep = ""),
    input.cls = paste(output.dir, identifier, ".model_set.2.cls", sep = ""),
    output.W.file = paste(output.dir, identifier, ".model.W.gct", sep = ""),
    output.H.file = paste(output.dir, identifier, ".model.H.gct", sep = ""),
    k.proj = k.proj,
    alg = alg,
    niter = niter,
    seed = seed,
    theta = theta,
    sort.factors = T
  )

  O <- MSIG.Factors.Project(
    input.ds = paste(output.dir, identifier, ".model_set.2.gct", sep = ""),
    factors.ds = paste(output.dir, identifier, ".model.W.gct", sep = ""),
    postprojnorm = postprojnorm,
    output.file = paste(output.dir, identifier, ".model_dataset.H.gct", sep = "")
  )

  O <- MSIG.Subset.Dataset(
    input.ds = paste(output.dir, identifier, ".model_dataset.H.gct", sep = ""),
    input.cls = paste(output.dir, identifier, ".model_set.2.cls", sep = ""),
    column.subset = "ALL",
    column.sel.type = "samples",
    row.subset = "ALL",
    output.ds = paste(output.dir, identifier, ".all.H.gct", sep = ""),
    output.cls = paste(output.dir, identifier, ".all.H.cls", sep = "")
  )

  O <- MSIG.Subset.Dataset(
    input.ds = paste(output.dir, identifier, ".model_set.2.gct", sep = ""),
    input.cls = paste(output.dir, identifier, ".model_set.2.cls", sep = ""),
    column.subset = "ALL",
    column.sel.type = "samples",
    row.subset = "ALL",
    output.ds = paste(output.dir, identifier, ".all.gct", sep = ""),
    output.cls = paste(output.dir, identifier, ".all.cls", sep = "")
  )

  # Pre-process and project all the test datasets

  if (!is.null(test.datasets.table)) {
    max.files <- length(test.datasets.table)
    for (ds in 1:max.files) {
      print(c("Processing test file: ", test.datasets.table[[ds]]$gct.file))

      O <- MSIG.Subset.Dataset(
        input.ds = test.datasets.table[[ds]]$gct.file,
        input.cls = test.datasets.table[[ds]]$cls.file,
        column.subset = test.datasets.table[[ds]]$column.subset,
        column.sel.type = test.datasets.table[[ds]]$column.sel.type,
        row.subset = "ALL",
        output.ds = paste(output.dir, "temp1.gct", sep = ""),
        output.cls = paste(output.dir, "temp1.cls", sep = "")
      )

      O <- MSIG.Preprocess.Dataset(
        input.ds = paste(output.dir, "temp1.gct", sep = ""),
        output.ds = paste(output.dir, "temp2.gct", sep = ""),
        thres = test.datasets.table[[ds]]$thres,
        ceil = test.datasets.table[[ds]]$ceil,
        normalization = "NULL"
      )

      #### New methodology: match first before normalization

      O <- MSIG.Match.and.Select(
        input1.ds = paste(output.dir, identifier, ".model_set.2.gct", sep = ""),
        input2.ds = paste(output.dir, "temp2.gct", sep = ""),
        output.ds = paste(output.dir, "temp3.gct", sep = "")
      )

      O <- MSIG.Preprocess.Dataset(
        input.ds = paste(output.dir, "temp3.gct", sep = ""),
        output.ds = paste(output.dir, "temp4.gct", sep = ""),
        normalization = test.datasets.table[[ds]]$norm
      )

      ####

      O <- MSIG.Factors.Project(
        input.ds = paste(output.dir, "temp4.gct", sep = ""),
        factors.ds = paste(output.dir, identifier, ".model.W.gct", sep = ""),
        postprojnorm = postprojnorm,
        output.file = paste(output.dir, "temp5.gct", sep = "")
      )

      O <- MSIG.Match.and.Merge(
        input1.ds = paste(output.dir, identifier, ".all.H.gct", sep = ""),
        input1.cls = paste(output.dir, identifier, ".all.H.cls", sep = ""),
        input2.ds = paste(output.dir, "temp5.gct", sep = ""),
        input2.cls = paste(output.dir, "temp1.cls", sep = ""),
        output.ds = paste(output.dir, identifier, ".all.H.gct", sep = ""),
        output.cls = paste(output.dir, identifier, ".all.H.cls", sep = "")
      )

      O <- MSIG.Match.and.Merge(
        input1.ds = paste(output.dir, identifier, ".all.gct", sep = ""),
        input1.cls = paste(output.dir, identifier, ".all.cls", sep = ""),
        input2.ds = paste(output.dir, "temp4.gct", sep = ""),
        input2.cls = paste(output.dir, "temp1.cls", sep = ""),
        output.ds = paste(output.dir, identifier, ".all.gct", sep = ""),
        output.cls = paste(output.dir, identifier, ".all.cls", sep = "")
      )

      print(c("Done processing test file: ", test.datasets.table[[ds]]$gct.file))
    }
  }

  # Plots for NMP

  CLS <- MSIG.ReadClsFile(file = paste(output.dir, identifier, ".model_set.2.cls", sep = ""))
  model.size <- length(CLS$class.v)

  O <- MSIG.Projection.Plots.3(
    input.ds = paste(output.dir, identifier, ".all.H.gct", sep = ""),
    input.cls = paste(output.dir, identifier, ".all.H.cls", sep = ""),
    model.set = seq(1, model.size),
    output.2D.proj.file = paste(output.dir, identifier, ".2D.proj.gct", sep = ""),
    output.2D.proj.plot = paste(output.dir, identifier, ".2D.proj", sep = ""),
    output.3D.proj.file = paste(output.dir, identifier, ".3D.proj.gct", sep = ""),
    output.3D.1.proj.plot = paste(output.dir, identifier, ".3D.1.proj", sep = ""),
    output.3D.2.proj.plot = paste(output.dir, identifier, ".3D.2.proj", sep = ""),
    output.3D.3.proj.plot = paste(output.dir, identifier, ".3D.3.proj", sep = ""),
    output.heatmap.plot = paste(output.dir, identifier, ".heatmap", sep = ""),
    output.heatmap.sorted.plot = paste(output.dir, identifier, ".heatmap.sorted", sep = ""),
    output.heatmap.sorted.2.plot = paste(output.dir, identifier, ".heatmap.sorted.2", sep = ""),
    output.hclust.plot = paste(output.dir, identifier, ".hclust", sep = ""),
    use.biplot = use.biplot,
    title = identifier,
    seed = seed,
    non.interactive.run = non.interactive.run,
    heatmap.row.norm = heatmap.row.norm,
    heatmap.cmap.type = heatmap.cmap.type,
    symbol.scaling = symbol.scaling,
    col = col,
    symbs = symbs
  )

  # Evaluate projection

  O <- MSIG.Evaluate.Projection.2(
    input.ds = paste(output.dir, identifier, ".all.H.gct", sep = ""),
    input.cls = paste(output.dir, identifier, ".all.H.cls", sep = ""),
    model.set = seq(1, model.size),
    prediction.results.file = paste(output.dir, identifier, ".pred.txt", sep = ""),
    prediction.matrix.file = paste(output.dir, identifier, ".pred.gct", sep = ""),
    train.pred.plot = paste(output.dir, identifier, ".train.pred", sep = ""),
    test.pred.plot = paste(output.dir, identifier, ".test.pred", sep = ""),
    pred.2D.plot = paste(output.dir, identifier, ".2D.pred", sep = ""),
    col = col,
    use.feature.names = use.feature.names,
    nchar.phen = nchar.phen,
    high.conf.thres = high.conf.thres,
    symbol.scaling = symbol.scaling,
    symbs = symbs,
    levels = levels,
    nlevels = nlevels,
    kernel = kernel,
    cost = cost,
    gamma = gamma,
    produce.contours = produce.contours
  )

  # Compute hierarchical clustering

  input.ds <- paste(output.dir, identifier, ".all.gct", sep = "")
  input.cls <- paste(output.dir, identifier, ".all.cls", sep = "")
  dataset <- MSIG.Gct2Frame(filename = input.ds)
  m.ds <- data.matrix(dataset$ds)
  N.ds <- length(m.ds[, 1])
  M.ds <- length(m.ds[1, ])
  ds.names <- dataset$row.names
  ds.descs <- dataset$descs
  ds.sample.names <- dataset$names

  # Read class labels

  CLS <- MSIG.ReadClsFile(file = input.cls)
  class.labels <- CLS$class.v
  class.phen <- CLS$phen
  class.list <- CLS$class.list

  class.labels <- match(class.list, class.phen)
  class.phen <- unique(class.list)

  # Compute hierarchical tree clustering

  dist.matrix <- dist(t(m.ds))
  HC <- hclust(dist.matrix, method = "complete")

  x11(height = 24, width = 30)

  #  plot(HC, xlab="samples", cex = 0.7, labels = class.list, col = "blue", main = " Hierarchical Clustering (original data)")

  #       HC$labels <- class.phen[class.labels]
  HC$labels <- class.list
  dhc <- as.dendrogram(HC, hang = 0.01, edge.root = T, dLeaf = 2, edgePar = list(col = c("blue", "green"), lty = c(1, 1), lwd = c(2, 2), t.col = c(1, 1)))
  local({
    colLab <<- function(n) {
      if (is.leaf(n)) {
        a <- attributes(n)
        i <<- i + 1
        attr(n, "nodePar") <-
          c(a$nodePar, list(lab.col = mycols[i], pch = c(0, 0), col = c(0, 0), bg = c(0, 0), cex = c(0.8, 0.8), lab.font = i %% 1))
      }
    }
    mycols <- col[class.labels[HC$order]]
    i <- 0
  })
  dL <- dendrapply(dhc, colLab)
  #     plot(dL, cex=1, edge.root = T, main = " Hierarchical Clustering (original data)", xlab = "samples") ## --> colored labels!

  plot.filename <- paste(output.dir, identifier, ".htree", sep = "")
  #   savePlot(filename = plot.filename, type ="jpeg", device = dev.cur())

  # Read projected dataset

  input.ds <- paste(output.dir, identifier, ".all.H.gct", sep = "")
  input.cls <- paste(output.dir, identifier, ".all.H.cls", sep = "")
  dataset <- MSIG.Gct2Frame(filename = input.ds)
  m.ds <- data.matrix(dataset$ds)
  N.ds <- length(m.ds[, 1])
  M.ds <- length(m.ds[1, ])
  ds.names <- dataset$row.names
  ds.descs <- dataset$descs
  ds.sample.names <- dataset$names

  # Compute hierarchical tree clustering

  dist.matrix <- dist(t(m.ds))
  HC <- hclust(dist.matrix, method = "complete")

  x11(height = 24, width = 30)
  #   plot(HC, xlab="samples", cex = 0.7, labels = class.list, col = "blue", main = " Hierarchical Clustering (projected data)")

  HC$labels <- class.list
  dhc <- as.dendrogram(HC, hang = 0.01, edge.root = T, dLeaf = 2, edgePar = list(col = c("blue", "green"), lty = c(1, 1), lwd = c(2, 2), t.col = c(1, 1)))
  local({
    colLab <<- function(n) {
      if (is.leaf(n)) {
        a <- attributes(n)
        i <<- i + 1
        attr(n, "nodePar") <-
          c(a$nodePar, list(lab.col = mycols[i], pch = c(0, 0), col = c(0, 0), bg = c(0, 0), cex = c(0.8, 0.8), lab.font = i %% 1))
      }
    }
    mycols <- col[class.labels[HC$order]]
    i <- 0
  })
  dL <- dendrapply(dhc, colLab)
  #     plot(dL, cex=1, edge.root = T, main = " Hierarchical Clustering (projected data)", xlab = "samples") ## --> colored labels!


  #   plot.filename <- paste(output.dir, identifier, ".H.htree", sep="")
  savePlot(filename = plot.filename, type = "jpeg", device = dev.cur())

  # Compute class membership

  membership <- vector(length = M.ds, mode = "numeric")

  for (j in 1:M.ds) { # Find membership
    membership[j] <- order(m.ds[, j], decreasing = T)
  }

  mem.order <- order(membership, decreasing = F)
  membership.sorted <- membership[mem.order]
  ds.sample.names <- paste(class.list, ds.sample.names, sep = "_")
  ds.sample.names.sorted <- ds.sample.names[mem.order]
  class.list.sorted <- class.list[mem.order]

  mem.table <- data.frame(cbind(class.list, ds.sample.names, membership, rep(" ", M.ds), class.list.sorted, ds.sample.names.sorted, membership.sorted))
  row.names(mem.table) <- seq(1, M.ds)
  names(mem.table) <- c("Phen", "Sample Names", "Membership", " ", "Phen Sorted", "Sample Names Sorted", "Membership Sorted")

  mem.filename <- paste(output.dir, identifier, ".H.mem.txt", sep = "")

  write.table(file = mem.filename, mem.table, quote = F, sep = "\t")

  table(class.list.sorted, membership.sorted)
}


Analyze.factors.2 <- function(
    input.ds,
    input.cls,
    output.dir,
    identifier,
    model.size,
    seed = 123,
    nchar.phen = 2,
    use.biplot = TRUE,
    non.interactive.run = FALSE,
    heatmap.row.norm = FALSE,
    heatmap.cmap.type = 3,
    use.feature.names = FALSE,
    high.conf.thres = 0.5,
    col = c("green", "blue", "pink", "red", "orange", "red4", "steelblue2", "violet"),
    symbs = c(22, 21, 20, 23, 24, 25, 21, 20, 23, 24, 25),
    symbol.scaling = 1,
    levels = NULL,
    nlevels = 20,
    kernel = "radial",
    cost = 1,
    gamma = 0.05,
    theta = 0,
    model.set.refinement = T,
    produce.contours = T) {
  print(c("Running MSIG.Analyze.factors.2... on: ", input.ds, " ", input.cls))

  # Plots for NMP

  # CLS <- MSIG.ReadClsFile(file =  paste(output.dir, identifier, ".model_set.2.cls", sep=""))
  # model.size <- length(CLS$class.v)

  # O <- MSIG.Projection.Plots.3(
  #     input.ds =                    input.ds,
  #     input.cls =                   input.cls,
  #        model.set =                   seq(1, model.size),
  #        output.2D.proj.file =         paste(output.dir, identifier, ".2D.proj.gct", sep=""),
  #        output.2D.proj.plot =         paste(output.dir, identifier, ".2D.proj", sep=""),
  #        output.3D.proj.file =         paste(output.dir, identifier, ".3D.proj.gct", sep=""),
  #        output.3D.1.proj.plot =       paste(output.dir, identifier, ".3D.1.proj", sep=""),
  #        output.3D.2.proj.plot =       paste(output.dir, identifier, ".3D.2.proj", sep=""),
  #       output.3D.3.proj.plot =       paste(output.dir, identifier, ".3D.3.proj", sep=""),
  #        output.heatmap.plot =         paste(output.dir, identifier, ".heatmap", sep=""),
  #        output.heatmap.sorted.plot =  paste(output.dir, identifier, ".heatmap.sorted", sep=""),
  #        output.heatmap.sorted.2.plot =  paste(output.dir, identifier, ".heatmap.sorted.2", sep=""),
  #        output.hclust.plot =          paste(output.dir, identifier, ".hclust", sep=""),
  #        use.biplot =                  use.biplot,
  #       title =                       identifier,
  #        seed =                        seed,
  #       non.interactive.run =         non.interactive.run,
  #        heatmap.row.norm =            heatmap.row.norm,
  #        heatmap.cmap.type =           heatmap.cmap.type,
  #        symbol.scaling =              symbol.scaling,
  #        col =                         col,
  #        symbs =                       symbs)

  # Evaluate projection

  O <- MSIG.Evaluate.Projection.2(
    input.ds = input.ds,
    input.cls = input.cls,
    model.set = seq(1, model.size),
    prediction.results.file = paste(output.dir, identifier, ".pred.txt", sep = ""),
    prediction.matrix.file = paste(output.dir, identifier, ".pred.gct", sep = ""),
    train.pred.plot = paste(output.dir, identifier, ".train.pred", sep = ""),
    test.pred.plot = paste(output.dir, identifier, ".test.pred", sep = ""),
    pred.2D.plot = paste(output.dir, identifier, ".2D.pred", sep = ""),
    col = col,
    use.feature.names = use.feature.names,
    nchar.phen = nchar.phen,
    high.conf.thres = high.conf.thres,
    symbol.scaling = symbol.scaling,
    symbs = symbs,
    levels = levels,
    nlevels = nlevels,
    kernel = kernel,
    cost = cost,
    gamma = gamma,
    produce.contours = produce.contours
  )

  # Compute hierarchical clustering

  dataset <- MSIG.Gct2Frame(filename = input.ds)
  m.ds <- data.matrix(dataset$ds)
  N.ds <- length(m.ds[, 1])
  M.ds <- length(m.ds[1, ])
  ds.names <- dataset$row.names
  ds.descs <- dataset$descs
  ds.sample.names <- dataset$names

  # Read class labels

  CLS <- MSIG.ReadClsFile(file = input.cls)
  class.labels <- CLS$class.v
  class.phen <- CLS$phen
  class.list <- CLS$class.list

  class.labels <- match(class.list, class.phen)
  class.phen <- unique(class.list)

  # Compute hierarchical tree clustering

  dist.matrix <- dist(t(m.ds))
  HC <- hclust(dist.matrix, method = "complete")

  x11(height = 24, width = 30)

  #  plot(HC, xlab="samples", cex = 0.7, labels = class.list, col = "blue", main = " Hierarchical Clustering (original data)")

  #       HC$labels <- class.phen[class.labels]
  HC$labels <- class.list
  dhc <- as.dendrogram(HC, hang = 0.01, edge.root = T, dLeaf = 2, edgePar = list(col = c("blue", "green"), lty = c(1, 1), lwd = c(2, 2), t.col = c(1, 1)))
  local({
    colLab <<- function(n) {
      if (is.leaf(n)) {
        a <- attributes(n)
        i <<- i + 1
        attr(n, "nodePar") <-
          c(a$nodePar, list(lab.col = mycols[i], pch = c(0, 0), col = c(0, 0), bg = c(0, 0), cex = c(0.8, 0.8), lab.font = i %% 1))
      }
      n
    }
    mycols <- col[class.labels[HC$order]]
    i <- 0
  })
  dL <- dendrapply(dhc, colLab)
  plot(dL, cex = 1, edge.root = T, main = " Hierarchical Clustering (original data)", xlab = "samples") ## --> colored labels!

  plot.filename <- paste(output.dir, identifier, ".htree", sep = "")
  savePlot(filename = plot.filename, type = "jpeg", device = dev.cur())

  # Read projected dataset

  dataset <- MSIG.Gct2Frame(filename = input.ds)
  m.ds <- data.matrix(dataset$ds)
  N.ds <- length(m.ds[, 1])
  M.ds <- length(m.ds[1, ])
  ds.names <- dataset$row.names
  ds.descs <- dataset$descs
  ds.sample.names <- dataset$names

  # Compute hierarchical tree clustering

  dist.matrix <- dist(t(m.ds))
  HC <- hclust(dist.matrix, method = "complete")

  x11(height = 24, width = 30)
  #   plot(HC, xlab="samples", cex = 0.7, labels = class.list, col = "blue", main = " Hierarchical Clustering (projected data)")

  HC$labels <- class.list
  dhc <- as.dendrogram(HC, hang = 0.01, edge.root = T, dLeaf = 2, edgePar = list(col = c("blue", "green"), lty = c(1, 1), lwd = c(2, 2), t.col = c(1, 1)))
  local({
    colLab <<- function(n) {
      if (is.leaf(n)) {
        a <- attributes(n)
        i <<- i + 1
        attr(n, "nodePar") <-
          c(a$nodePar, list(lab.col = mycols[i], pch = c(0, 0), col = c(0, 0), bg = c(0, 0), cex = c(0.8, 0.8), lab.font = i %% 1))
      }
      n
    }
    mycols <- col[class.labels[HC$order]]
    i <- 0
  })
  dL <- dendrapply(dhc, colLab)
  plot(dL, cex = 1, edge.root = T, main = " Hierarchical Clustering (projected data)", xlab = "samples") ## --> colored labels!


  plot.filename <- paste(output.dir, identifier, ".H.htree", sep = "")
  savePlot(filename = plot.filename, type = "jpeg", device = dev.cur())

  # Compute class membership

  membership <- vector(length = M.ds, mode = "numeric")

  for (j in 1:M.ds) { # Find membership
    membership[j] <- order(m.ds[, j], decreasing = T)
  }

  mem.order <- order(membership, decreasing = F)
  membership.sorted <- membership[mem.order]
  ds.sample.names <- paste(class.list, ds.sample.names, sep = "_")
  ds.sample.names.sorted <- ds.sample.names[mem.order]
  class.list.sorted <- class.list[mem.order]

  mem.table <- data.frame(cbind(class.list, ds.sample.names, membership, rep(" ", M.ds), class.list.sorted, ds.sample.names.sorted, membership.sorted))
  row.names(mem.table) <- seq(1, M.ds)
  names(mem.table) <- c("Phen", "Sample Names", "Membership", " ", "Phen Sorted", "Sample Names Sorted", "Membership Sorted")

  mem.filename <- paste(output.dir, identifier, ".H.mem.txt", sep = "")

  write.table(file = mem.filename, mem.table, quote = F, sep = "\t")

  table(class.list.sorted, membership.sorted)
}

MSIG.Match.and.Merge.2 <- function(
    input1.ds,
    input1.cls = "",
    input2.ds,
    input2.cls = "",
    output.ds,
    output.cls,
    mode = "intersection") { # gene set mode = intersection (default) or union

  # start of methodology

  print(c("Running MSIG.Match.and.Merge... on: ", input1.ds, " ", input2.ds))

  # Read input datasets

  dataset1 <- MSIG.Gct2Frame(filename = input1.ds)
  m1 <- data.matrix(dataset1$ds)
  gs.names1 <- dataset1$row.names
  gs.descs1 <- dataset1$descs
  sample.names1 <- dataset1$names
  N1 <- length(m1[1, ])
  Ng1 <- length(m1[, 1])

  print(c("dataset 1:", N1, " samples"))
  print(c("dataset 1:", Ng1, " genes"))

  dataset2 <- MSIG.Gct2Frame(filename = input2.ds)
  m2 <- data.matrix(dataset2$ds)
  gs.names2 <- dataset2$row.names
  gs.descs2 <- dataset2$descs
  sample.names2 <- dataset2$names
  N2 <- length(m2[1, ])
  Ng2 <- length(m2[, 1])

  print(c("dataset 2:", N2, " samples"))
  print(c("dataset 2:", Ng2, " genes"))

  # Read CLS files

  if ((input1.cls != "") & (input2.cls != "")) {
    CLS1 <- ReadClsFile(file = input1.cls)
    class.labels1 <- CLS1$class.v
    class.phen1 <- CLS1$phen

    CLS2 <- ReadClsFile(file = input2.cls)
    class.labels2 <- CLS2$class.v
    class.phen2 <- CLS2$phen
  }

  # Match features to first dataset and create matching m2 dataset

  sample.names3 <- c(sample.names1, sample.names2)
  N3 <- N1 + N2

  if (mode == "intersection") {
    gs.names3 <- intersect(gs.names1, gs.names2)

    locations1 <- match(gs.names3, gs.names1, nomatch = 0)
    m1 <- m1[locations1, ]
    gs.descs1 <- gs.descs1[locations1]

    locations2 <- match(gs.names3, gs.names2, nomatch = 0)
    m2 <- m2[locations2, ]
    gs.descs2 <- gs.descs2[locations2]
    m3 <- cbind(m1, m2)
  } else if (mode == "union") {
    gs.names3 <- union(gs.names1, gs.names2)
    M3 <- length(gs.names3)

    m3 <- matrix(0, nrow = M3, ncol = N3)

    locations1 <- match(gs.names3, gs.names1, nomatch = 0)
    locations2 <- match(gs.names3, gs.names2, nomatch = 0)

    for (i in 1:M3) {
      if (locations1[i] != 0) {
        m3[i, 1:N1] <- m1[locations1[i], ]
      }
      if (locations2[i] != 0) {
        m3[i, (N1 + 1):N3] <- m2[locations2[i], ]
      }
    }
  } else {
    stop(c("unknown mode", mode))
  }


  if ((input1.cls != "") & (input2.cls != "")) {
    class.labels3 <- c(class.labels1, class.labels2 + length(class.phen1))
    class.phen3 <- c(class.phen1, class.phen2)
  }

  # Save datasets

  N3 <- length(m3[1, ])
  Ng3 <- length(m3[, 1])

  print(c("dataset 3:", N3, " samples"))
  print(c("dataset 3:", Ng3, " genes"))

  V <- data.frame(m3)
  names(V) <- sample.names3
  row.names(V) <- gs.names3
  gs.descs1 <- gs.names3
  write.gct(gct.data.frame = V, descs = gs.descs1, filename = output.ds)

  if ((input1.cls != "") & (input2.cls != "")) {
    write.cls(class.v = class.labels3, phen = class.phen3, filename = output.cls)
  }
}

Analyze.factors <- function(
    input.ds,
    input.cls,
    output.dir,
    identifier,
    model.size,
    seed = 123,
    nchar.phen = 2,
    use.biplot = TRUE,
    non.interactive.run = FALSE,
    heatmap.row.norm = FALSE,
    heatmap.cmap.type = 3,
    use.feature.names = FALSE,
    high.conf.thres = 0.5,
    col = c("green", "blue", "pink", "red", "orange", "red4", "steelblue2", "violet"),
    symbs = c(22, 21, 20, 23, 24, 25, 21, 20, 23, 24, 25),
    symbol.scaling = 1,
    levels = NULL,
    nlevels = 20,
    kernel = "radial",
    cost = 1,
    gamma = 0.05,
    theta = 0,
    model.set.refinement = T) {
  # This is very much like the plots and analysis second part of MPM.2

  # CLS <- MSIG.ReadClsFile(file =  input.cls)
  # model.size <- length(CLS$class.v)

  O <- MSIG.Projection.Plots.3(
    input.ds = input.ds,
    input.cls = input.cls,
    model.set = seq(1, model.size),
    output.2D.proj.file = paste(output.dir, identifier, ".2D.proj.gct", sep = ""),
    output.2D.proj.plot = paste(output.dir, identifier, ".2D.proj", sep = ""),
    output.3D.proj.file = paste(output.dir, identifier, ".3D.proj.gct", sep = ""),
    output.3D.1.proj.plot = paste(output.dir, identifier, ".3D.1.proj", sep = ""),
    output.3D.2.proj.plot = paste(output.dir, identifier, ".3D.2.proj", sep = ""),
    output.3D.3.proj.plot = paste(output.dir, identifier, ".3D.3.proj", sep = ""),
    output.heatmap.plot = paste(output.dir, identifier, ".heatmap", sep = ""),
    output.heatmap.sorted.plot = paste(output.dir, identifier, ".heatmap.sorted", sep = ""),
    output.heatmap.sorted.2.plot = paste(output.dir, identifier, ".heatmap.sorted.2", sep = ""),
    output.hclust.plot = paste(output.dir, identifier, ".hclust", sep = ""),
    use.biplot = use.biplot,
    title = identifier,
    seed = seed,
    non.interactive.run = non.interactive.run,
    heatmap.row.norm = heatmap.row.norm,
    heatmap.cmap.type = heatmap.cmap.type,
    symbol.scaling = symbol.scaling,
    col = col,
    symbs = symbs
  )

  # Evaluate projection

  O <- MSIG.Evaluate.Projection.2(
    input.ds = input.ds,
    input.cls = input.cls,
    model.set = seq(1, model.size),
    prediction.results.file = paste(output.dir, identifier, ".pred.txt", sep = ""),
    prediction.matrix.file = paste(output.dir, identifier, ".pred.gct", sep = ""),
    train.pred.plot = paste(output.dir, identifier, ".train.pred", sep = ""),
    test.pred.plot = paste(output.dir, identifier, ".test.pred", sep = ""),
    pred.2D.plot = paste(output.dir, identifier, ".2D.pred", sep = ""),
    col = col,
    use.feature.names = use.feature.names,
    nchar.phen = nchar.phen,
    high.conf.thres = high.conf.thres,
    symbol.scaling = symbol.scaling,
    symbs = symbs,
    levels = levels,
    nlevels = nlevels,
    kernel = kernel,
    cost = cost,
    gamma = gamma
  )

  # Compute hierarchical clustering

  dataset <- MSIG.Gct2Frame(filename = input.ds)
  m.ds <- data.matrix(dataset$ds)
  N.ds <- length(m.ds[, 1])
  M.ds <- length(m.ds[1, ])
  ds.names <- dataset$row.names
  ds.descs <- dataset$descs
  ds.sample.names <- dataset$names

  # Read class labels

  CLS <- MSIG.ReadClsFile(file = input.cls)
  class.labels <- CLS$class.v
  class.phen <- CLS$phen
  class.list <- CLS$class.list

  class.labels <- match(class.list, class.phen)
  class.phen <- unique(class.list)

  # Compute hierarchical tree clustering

  dist.matrix <- dist(t(m.ds))
  HC <- hclust(dist.matrix, method = "complete")

  x11(height = 24, width = 30)
  plot(HC, xlab = "samples", cex = 0.7, labels = class.list, col = "black", main = " Hierarchical Clustering (original data)")

  plot.filename <- paste(output.dir, identifier, ".htree", sep = "")
  savePlot(filename = plot.filename, type = "jpeg", device = dev.cur())

  # Read projected dataset

  dataset <- MSIG.Gct2Frame(filename = input.ds)
  m.ds <- data.matrix(dataset$ds)
  N.ds <- length(m.ds[, 1])
  M.ds <- length(m.ds[1, ])
  ds.names <- dataset$row.names
  ds.descs <- dataset$descs
  ds.sample.names <- dataset$names

  # Compute hierarchical tree clustering

  dist.matrix <- dist(t(m.ds))
  HC <- hclust(dist.matrix, method = "complete")

  x11(height = 24, width = 30)
  plot(HC, xlab = "samples", cex = 0.7, labels = class.list, col = "blue", main = " Hierarchical Clustering (projected data)")

  plot.filename <- paste(output.dir, identifier, ".H.htree", sep = "")
  savePlot(filename = plot.filename, type = "jpeg", device = dev.cur())

  # Compute class membership

  membership <- vector(length = M.ds, mode = "numeric")

  for (j in 1:M.ds) { # Find membership
    membership[j] <- order(m.ds[, j], decreasing = T)
  }

  mem.order <- order(membership, decreasing = F)
  membership.sorted <- membership[mem.order]
  ds.sample.names <- paste(class.list, ds.sample.names, sep = "_")
  ds.sample.names.sorted <- ds.sample.names[mem.order]
  class.list.sorted <- class.list[mem.order]

  mem.table <- data.frame(cbind(class.list, ds.sample.names, membership, rep(" ", M.ds), class.list.sorted, ds.sample.names.sorted, membership.sorted))
  row.names(mem.table) <- seq(1, M.ds)
  names(mem.table) <- c("Phen", "Sample Names", "Membership", " ", "Phen Sorted", "Sample Names Sorted", "Membership Sorted")

  mem.filename <- paste(output.dir, identifier, ".H.mem.txt", sep = "")

  write.table(file = mem.filename, mem.table, quote = F, sep = "\t")

  table(class.list.sorted, membership.sorted)
}


MSIG.Reconstruct.Dataset <- function(
    input.H.ds,
    input.W.ds,
    output.file) {
  library(MASS)

  # start of methodology

  print(c("Running MSIG.Reconstruct.Dataset... on: ", input.H.ds, input.W.ds))

  # Read input datasets

  dataset <- MSIG.Gct2Frame(filename = input.W.ds)
  W <- data.matrix(dataset$ds)
  W.row.names <- dataset$row.names
  W.row.descs <- dataset$descs
  W.names <- dataset$names

  dataset <- MSIG.Gct2Frame(filename = input.H.ds)
  H <- data.matrix(dataset$ds)
  H.row.names <- dataset$row.names
  H.row.descs <- dataset$descs
  H.names <- dataset$names

  # Project input dataset using factors input

  A <- W %*% H

  # Save reconstructed dataset

  V <- data.frame(A)
  names(V) <- H.names
  row.names(V) <- W.row.names
  write.gct(gct.data.frame = V, filename = output.file)
}

MSIG.HeatMapPlot.5 <- function(
    V,
    row.names = "NA",
    col.labels = "NA",
    col.classes = "NA",
    col.names = "NA",
    main = " ",
    sub = " ",
    xlab = " ",
    ylab = " ",
    row.norm = TRUE,
    char.rescale = 1.0,
    cmap.type = 1, # 1 = vintage pinkogram, 2 = scale of grays, 3 = high-resolution pinkogram for probabilities [0, 1], 4 = high-resolution pinkogram for general values, 5 = color map for normalized enrichment scores, 6 = color map for metagene factors
    max.v = "NA",
    rotated.col.labels = F,
    create.legend = T,
    create.window = T) {
  #
  # Plots a heatmap "pinkogram" of a gene expression matrix including phenotype vector and gene, sample and phenotype labels
  #
  # The Broad Institute
  # SOFTWARE COPYRIGHT NOTICE AGREEMENT
  # This software and its documentation are copyright 2003 by the
  # Broad Institute/Massachusetts Institute of Technology.
  # All rights are reserved.
  #
  # This software is supplied without any warranty or guaranteed support
  # whatsoever. Neither the Broad Institute nor MIT can be responsible for
  # its use, misuse, or functionality.

  n.rows <- length(V[, 1])
  n.cols <- length(V[1, ])
  n.phen <- length(col.classes)
  V1 <- matrix(0, nrow = n.rows, ncol = n.cols)

  if ((cmap.type == 3) | (cmap.type == 5)) {
    row.norm <- F
  }

  if (row.norm == TRUE) {
    row.mean <- apply(V, MARGIN = 1, FUN = mean)
    row.sd <- apply(V, MARGIN = 1, FUN = sd)
    row.n <- length(V[, 1])
    for (i in 1:n.rows) {
      if (row.sd[i] == 0) {
        V1[i, ] <- 0
      } else {
        V1[i, ] <- (V[i, ] - row.mean[i]) / (0.333 * row.sd[i])
      }
      V1[i, ] <- ifelse(V1[i, ] < -6, -6, V1[i, ])
      V1[i, ] <- ifelse(V1[i, ] > 6, 6, V1[i, ])
    }
  } else {
    V1 <- V
  }

  if (cmap.type == 1) {
    mycol <- c(
      "#0000FF", "#4040FF", "#7070FF", "#8888FF", "#A9A9FF", "#D5D5FF", "#EEE5EE", "#FFAADA", "#FF9DB0", "#FF7080",
      "#FF5A5A", "#FF4040", "#FF0D1D"
    ) # blue-pinkogram colors. This is the 1998-vintage, pre-gene cluster, original pinkogram color map
  } else if (cmap.type == 2) {
    mycol <- vector(length = 256, mode = "numeric")
    for (k in 1:256) {
      red <- (k - 1) * 0.80 + 50
      green <- (k - 1) * 0.68 + 80
      blue <- (k - 1) * 0.60 + 100
      mycol[k] <- rgb(red, green, blue, maxColorValue = 255)
    }
    mycol <- rev(mycol)

    #          mycol <- c("#FCFBFD","#F4F2F8","#F8F7FB","#EFEDF5","#E1E1EF","#E8E7F2","#DADAEB","#C6C7E1","#D0D1E6","#BCBDDC","#A8A6CF",
    #                        "#B2B2D6","#9E9AC8","#8A87BF","#9491C4","#807DBA","#7260AB","#796FB3","#6A51A3","#5C3596","#63439D","#54278F","#460D83","#4D1A89","#3F007D")
  } else if ((cmap.type == 3) | (cmap.type == 4) | (cmap.type == 5)) {
    mycol <- vector(length = 512, mode = "numeric")

    for (k in 1:256) {
      mycol[k] <- rgb(255, k - 1, k - 1, maxColorValue = 255)
    }
    for (k in 257:512) {
      mycol[k] <- rgb(511 - (k - 1), 511 - (k - 1), 255, maxColorValue = 255)
    }
    mycol <- rev(mycol)
  } else if (cmap.type == 6) {
    mycol <- vector(length = 256, mode = "numeric")
    max.pos.V1 <- max(V1)
    if (min(V1) < 0) {
      min.neg.V1 <- min(V1)
    } else {
      min.neg.V1 <- 0
    }
    neg.k <- ceiling(256 * (abs(min.neg.V1) / (max.pos.V1 - min.neg.V1)))
    #            print(c("neg.k=", neg.k, " max.pos.V1=", max.pos.V1, "min.neg.V1=", min.neg.V1))
    for (k in 1:neg.k) {
      max.red <- 255 - (255 - 50) * abs(min.neg.V1)
      min.red <- 255
      red <- max.red + (min.red - max.red) * (k - 1) / (neg.k - 1)
      max.green <- (255 - (255 - 50) * abs(min.neg.V1))
      min.green <- 255
      green <- max.green + (min.green - max.green) * (k - 1) / (neg.k - 1)
      max.blue <- (255 - (255 - 200) * abs(min.neg.V1))
      min.blue <- 255
      blue <- max.blue + (min.blue - max.blue) * (k - 1) / (neg.k - 1)
      mycol[k] <- rgb(red, green, blue, maxColorValue = 255)
    }
    for (k in (neg.k + 1):256) {
      max.red <- 205
      min.red <- 255
      red <- min.red - (min.red - max.red) * (k - (neg.k + 1)) / (256 - (neg.k + 1))
      max.green <- 50
      min.green <- 255
      green <- min.green - (min.green - max.green) * (k - (neg.k + 1)) / (256 - (neg.k + 1))
      max.blue <- 50
      min.blue <- 255
      blue <- min.blue - (min.blue - max.blue) * (k - (neg.k + 1)) / (256 - (neg.k + 1))
      mycol[k] <- rgb(red, green, blue, maxColorValue = 255)
    }
  }

  ncolors <- length(mycol)

  if (cmap.type == 5) {
    if (max.v == "NA") {
      max.v <- max(max(V1), -min(V1))
    }
    V2 <- ceiling(ncolors * (V1 - (-max.v)) / (1.001 * (max.v - (-max.v))))
  } else if (cmap.type == 3) {
    V2 <- ceiling(ncolors * (V1 / 1.001))
  } else {
    V2 <- ceiling(ncolors * (V1 - min(V1)) / (1.001 * (max(V1) - min(V1))))
  }

  heatm <- matrix(0, nrow = n.rows, ncol = n.cols)
  heatm[1:n.rows, ] <- V2[, seq(n.cols, 1, -1)]
  ###
  col.labels <- col.labels[seq(n.cols, 1, -1)]
  if (col.names[1] != "NA") {
    col.names <- col.names[seq(n.cols, 1, -1)]
  }
  height <- ifelse(n.rows >= 25, 25, n.rows * 0.8 + 2)
  #        x11(width=24, height=14, rescale="fit")


  if (create.window == T) {
    x11(width = 31, height = 19)
  }
  if ((create.window == T && create.legend == T)) {
    nf <- layout(matrix(c(1, 2), 1, 2, byrow = T), widths = c(7, 1), respect = FALSE)
  }
  n.rows2 <- ifelse(n.rows < 4, 4, n.rows)
  a <- -12 / 16
  b <- 27
  margin <- a * n.rows2 + b
  margin <- ifelse(margin < 2, 2, ifelse(margin > 24, 24, margin))
  #       margin <- ifelse(col.names[1] != "NA", margin + 10, margin)
  #       par(mar = c(2, margin, 4, 10))

  if (rotated.col.labels == F) {
    par(mar = c(4, margin, 4, 10))
  } else {
    par(mar = c(4, margin, 10, 10))
  }
  #        par(mar = c(2, 24, 4, 10))
  #       par(mar = c(2, 24, 4, 14))

  image(1:n.rows, 1:n.cols, heatm, zlim = c(0, ncolors), col = mycol, axes = FALSE, main = main, sub = sub, xlab = xlab, ylab = ylab)

  if (row.names[1] != "NA") {
    numC <- nchar(row.names)
    size.row.char <- char.rescale * 15 / (n.rows + 15)
    size.col.char <- char.rescale * 30 / (n.cols + 15)
    size.lab.char <- char.rescale * 30 / (n.phen + 15)
    for (i in 1:n.rows) {
      row.names[i] <- substr(row.names[i], 1, 20)
    }
    if (rotated.col.labels == F) {
      axis(3, at = 1:n.rows, labels = row.names, adj = 1, tick = FALSE, las = 1, cex.axis = size.row.char, font.axis = 1, line = -0.75)
    } else {
      axis(3, at = 1:n.rows, labels = row.names, adj = 1, tick = FALSE, las = 3, cex.axis = size.row.char, font.axis = 1, line = -0.75)
    }
  }

  #       print(c(" col.names=", col.names))

  if (col.names[1] != "NA") {
    axis(2, at = 1:n.cols, labels = col.names, tick = FALSE, las = 1, cex.axis = size.col.char, font.axis = 2.5, line = -1)
  }

  ####
  for (i in 1:(n.rows)) {
    lines(x = c(i + 0.5, i + 0.5), y = c(0.5, n.cols + 0.5), type = "l", lwd = 1, col = "black")
  }

  boundaries <- cumsum(sapply(split(rep(1, n.cols), col.labels), sum))
  boundaries <- n.cols - boundaries

  loc.phen <- vector(length = n.phen, mode = "numeric")
  for (i in 1:(n.phen)) {
    lines(x = c(1 - 0.5, n.rows + 0.5), y = c(boundaries[i] + 0.5, boundaries[i] + 0.5), type = "l", lwd = 1, col = "black")
    if (i > 1) {
      loc.phen[i] <- mean(c(boundaries[i - 1], boundaries[i]))
    } else {
      loc.phen[i] <- mean(c(n.cols, boundaries[i]))
    }
    #          print(c(i, " line=", boundaries[i], " label=", loc.phen[i]))
  }
  axis(4, at = loc.phen, labels = col.classes, tick = FALSE, las = 1, cex.axis = size.lab.char, font.axis = 1, line = -0.75)


  lines(x = c(0.50, n.rows + 0.50), y = c(0.50, 0.50), type = "l", lwd = 1, col = "black")
  lines(x = c(0.50, n.rows + 0.50), y = c(n.cols + 0.50, n.cols + 0.50), type = "l", lwd = 1, col = "black")
  lines(x = c(0.50, 0.50), y = c(0.50, n.cols + 0.50), type = "l", lwd = 1, col = "black")
  lines(x = c(n.rows + 0.50, n.rows + 0.50), y = c(0.50, n.cols + 0.50), type = "l", lwd = 1, col = "black")

  ####

  # Color map legend

  if (create.legend == T) {
    #       print(c("range V=", range(V)))
    #       print(c("range V1=", range(V1)))
    #       print(c("range V2=", range(V2)))

    par(mar = c(10, 2, 10, 2))
    num.v <- 20
    if (cmap.type == 3) {
      range.v <- c(0, ncolors)
    } else {
      range.v <- range(V2)
    }
    incr <- (range.v[1] - range.v[2]) / (num.v - 1)
    heatm.v <- matrix(rev(seq(range.v[2], range.v[1], incr)), nrow = num.v, ncol = 1)
    image(1:1, 1:num.v, t(heatm.v), zlim = c(0, ncolors + max(col.labels)), col = mycol, axes = FALSE, sub = "Color \n Legend ", main = " ", xlab = xlab, ylab = ylab)
    if (cmap.type == 3) {
      range.v <- c(0, 1)
    } else {
      range.v <- range(V1)
    }
    incr <- (range.v[1] - range.v[2]) / (num.v - 1)
    heatm.v2 <- matrix(signif(rev(seq(range.v[2], range.v[1], incr)), digits = 3), nrow = num.v, ncol = 1)
    #          print(c("heatm.v2=", heatm.v2))
    axis(2, at = 1:num.v, labels = heatm.v2, adj = 0.5, tick = FALSE, las = 1, cex.axis = char.rescale * 0.6, font.axis = 1.25, line = -0.8)

    lines(x = c(0.5, 1.5), y = c(0.5, 0.5), type = "l", lwd = 1, col = "black")
    lines(x = c(0.5, 1.5), y = c(num.v + 0.5, num.v + 0.5), type = "l", lwd = 1, col = "black")

    lines(x = c(0.6, 0.6), y = c(0.5, num.v + 0.5), type = "l", lwd = 1, col = "black")
    lines(x = c(1.4, 1.4), y = c(0.5, num.v + 0.5), type = "l", lwd = 1, col = "black")
  }

  return()
}

MSIG.Evaluate.Projection.2 <- function(
    input.ds,
    input.cls,
    model.set,
    prediction.results.file,
    prediction.matrix.file,
    train.pred.plot,
    test.pred.plot,
    pred.2D.plot,
    col = c("grey3", "cadetblue", "darkgreen", "chartreuse2", "red1", "darkred", "orange", "blue2", "lightblue", "pink", "coral"),
    symbs = c(22, 21, 20, 23, 24, 25, 21, 20, 23, 24, 25),
    non.interactive.run = F,
    use.feature.names = F,
    nchar.phen = 3,
    high.conf.thres = 0.75,
    symbol.scaling = symbol.scaling,
    levels = NULL,
    nlevels = 10,
    kernel = "radial",
    cost = 5,
    gamma = 0.05,
    produce.contours = T) {
  print(c("Running MSIG.Evaluate.Projection2... on:", input.ds))

  library(e1071)
  library(tree)

  # Read dataset

  dataset <- MSIG.Gct2Frame(filename = input.ds)
  m <- data.matrix(dataset$ds)
  max.m <- max(m)
  m <- m / max.m
  gs.names <- dataset$row.names
  gs.descs <- dataset$descs
  sample.names <- dataset$names

  dim(m)

  Ns <- length(m[1, ])
  N <- length(m[, 1])

  CLS <- MSIG.ReadClsFile(file = input.cls)
  class.labels <- CLS$class.v
  class.list <- CLS$class.list
  class.phen <- CLS$phen

  num.classes <- length(class.phen)

  print("Reading dataset completed...")

  # Use first nchar.phen characters of phenotype class to define new phenotypes

  class.list2 <- vector(length = Ns, mode = "character")
  for (i in 1:Ns) {
    class.list2[i] <- substr(class.list[i], 1, nchar.phen)
  }
  class.phen2 <- vector(length = num.classes, mode = "character")
  for (i in 1:num.classes) {
    class.phen2[i] <- substr(class.phen[i], 1, nchar.phen)
  }
  true.num.classes <- length(table(class.phen2))

  class.labels2 <- match(class.list2, class.phen2)

  # Separate data into train and test pieces

  m.train <- m[, model.set]
  n.train <- length(model.set)
  num.samples.train <- n.train
  sample.names.train <- as.factor(sample.names[model.set])
  class.list.train <- class.list2[model.set]
  class.phen.train <- unique(class.list.train)
  class.labels.train <- class.labels2[model.set]
  orig.class.labels.train <- class.labels[model.set]

  if (Ns - length(model.set) > 0) {
    m.test <- as.matrix(m[, -model.set])
    n.test <- length(m.test[1, ])
    sample.names.test <- as.factor(sample.names[-model.set])
    class.list.test <- class.list2[-model.set]
    class.phen.test <- unique(class.list.test)
    class.labels.test <- class.labels2[-model.set]
  }

  # Build SVM and tree models

  print("Building SVM model...")

  #  svm.model <- svm(x = t(m.train), y = class.list.train, scale = F, tolerance = 0.0001, type = "C-classification", kernel = "linear", cost = 3, probability = T)

  #   svm.model <- svm(x = t(m.train), y = class.list.train, scale = F, tolerance = 0.0001, type = "C-classification", kernel = "radial", cost = 3, gamma = 25.0, probability = T)

  #      svm.model <- svm(x = t(m.train), y = class.list.train, scale = F, tolerance = 0.00001, type = "C-classification", kernel = "radial", cost = 5, gamma = 2.5, probability = T)

  #      svm.model <- svm(x = t(m.train), y = class.list.train, scale = F, tolerance = 0.00001, type = "C-classification", kernel = "radial", cost = 5, gamma = 1, probability = T)

  #      svm.model <- svm(x = t(m.train), y = class.list.train, scale = F, tolerance = 0.00001, type = "C-classification", kernel = kernel, cost = cost, gamma = gamma, probability = T)

  one.over <- function(x) {
    return(100 / length(x))
  }
  class.number.list <- split(rep(1, length(class.list.train)), class.list.train)
  class.weights <- sapply(class.number.list, one.over)
  print(c("class.weights=", class.weights))

  #         svm.model <- svm(x = t(m.train), y = class.list.train, scale = F, tolerance = 0.00001, type = "C-classification", kernel = kernel, cost = cost, gamma = gamma, probability = T, class.weights = class.weights)

  svm.model <- svm(x = t(m.train), y = class.list.train, scale = F, tolerance = 0.00001, type = "C-classification", kernel = kernel, cost = cost, gamma = gamma, probability = T)

  #               svm.model <- svm(x = t(m.train), y = class.list.train, scale = F, tolerance = 0.00001, type = "C-classification", kernel = kernel, cost = cost, gamma = gamma, probability = T, cross=n.train)


  #   tree.data <- data.frame(as.factor(class.list), t(h.df))
  #   tree.model <- tree(formula = as.factor.class.list. ~ ., data = tree.data, split = "deviance")
  #   tree.data <- data.frame(cbind(class.list.train, t(m.train)))
  #   tree.data
  #   print(tree.data)
  #   tree.model <- tree(formula = class.list.train ~ ., data = tree.data)
  #   print(summary(tree.model))
  #   x11(height = 15, width = 15)
  #   plot(tree.model)
  #   text(tree.model)

  print("Computing train set predictions...")

  train.pred <- predict(object = svm.model, newdata = t(m.train), decision.values = T, probability = T)

  dec.vals.train <- attr(train.pred, "decision.values")
  prob.train <- signif(attr(train.pred, "probabilities"), digits = 2)
  confidence.vector <- vector(length = n.train, mode = "numeric")
  bscore <- vector(length = n.train, mode = "numeric")
  max.k <- length(prob.train[1, ])
  random.pred.conf <- ((max.k - 1) / max.k)^2 + (max.k - 1) * (1 / max.k)^2
  for (ii in 1:n.train) {
    probs <- sort(prob.train[ii, ], decreasing = T)
    confidence.vector[ii] <- 1 - ((1 - probs[1])^2 + sum(probs[2:max.k]^2)) / random.pred.conf
    confidence.vector[ii] <- signif(confidence.vector[ii], digits = 3)
    if (class.list.train[ii] == as.character(train.pred[ii])) {
      bscore[ii] <- signif((1 - probs[1])^2, digits = 2)
    } else {
      bscore[ii] <- signif(probs[1]^2, digits = 2)
    }
  }
  confidence.call <- ifelse(confidence.vector >= high.conf.thres, " H ", " L ")
  error.call <- ifelse(class.list.train == as.character(train.pred), "   ", " * ")
  no.call <- ifelse(confidence.vector >= high.conf.thres, 0, 1)
  real.error <- ifelse(((no.call == 0) & (error.call == " * ")), 1, 0)
  correct.call <- ifelse(((no.call == 0) & (error.call == "   ")), 1, 0)

  col.symbols.train <- paste(confidence.call, error.call)
  class.names <- names(data.frame(prob.train))
  Brier.train <- signif(mean(bscore), digits = 2)

  train.results <- data.frame(cbind(as.character(sample.names.train), class.list.train, as.character(train.pred), error.call, confidence.call, confidence.vector, no.call, real.error, correct.call, prob.train, bscore))
  names(train.results)[1] <- "Train Sample Name"
  names(train.results)[2] <- "Actual"
  names(train.results)[3] <- "Predicted"
  names(train.results)[4] <- "Error (*)"
  names(train.results)[5] <- "Conf (H/L)"
  names(train.results)[6] <- "Conf"
  names(train.results)[7] <- "No Call"
  names(train.results)[8] <- "Real Error"
  names(train.results)[9] <- "Correct Call"

  names(train.results)[10 + length(class.phen.train)] <- "Brier score"
  #   print(train.results)
  print(c("Brier score (Train) = ", Brier.train))

  write("Training Results \n", file = prediction.results.file, append = F)
  write.table(train.results, file = prediction.results.file, append = T, quote = F, row.names = F, sep = "\t")

  write(c("\n\n Brier score (Train) = ", Brier.train), file = prediction.results.file, append = T)

  no.call.list <- split(no.call, class.list.train)
  real.error.list <- split(real.error, class.list.train)
  correct.call.list <- split(correct.call, class.list.train)
  count.class <- c(sapply(no.call.list, length), length(no.call))
  no.call.class <- c(sapply(no.call.list, sum), sum(no.call))
  real.error.class <- c(sapply(real.error.list, sum), sum(real.error))
  correct.call.class <- c(sapply(correct.call.list, sum), sum(correct.call))
  train.pred.high.conf <- ifelse(no.call == 0, as.character(train.pred), "-- no call")
  #   print(c("train.pred.high.conf =", train.pred.high.conf))

  no.call.class.pct <- no.call.class / count.class
  real.error.class.pct <- real.error.class / count.class
  correct.call.class.pct <- correct.call.class / count.class
  perf.table.train <- data.frame(cbind(c(names(no.call.list), "Total"), count.class, no.call.class, no.call.class.pct, real.error.class, real.error.class.pct, correct.call.class, correct.call.class.pct))
  names(perf.table.train) <- c("Class", "Count", "No Call", "No Call (%)", "Real Error", "Real Error (%)", "Correct Call", "Correct Call (%)")
  write.table(perf.table.train, file = prediction.results.file, append = T, quote = F, row.names = F, sep = "\t")
  print(perf.table.train)

  conf.table.train <- table(class.list.train, train.pred.high.conf)
  conf.table.train <- data.frame(cbind(row.names(conf.table.train), conf.table.train))
  print(conf.table.train)
  write("\n\n Confusion Matrix (Train) \n", file = prediction.results.file, append = T)
  write.table(conf.table.train, file = prediction.results.file, append = T, quote = F, row.names = F, sep = "\t")
  height <- ifelse(length(class.phen.train) > 50, 20, 0.30 * length(class.phen.train) + 10)

  if (non.interactive.run == F) {
    if (.Platform$OS.type == "windows") {
      plot.filename <- train.pred.plot
      x11(height = height, width = 40)
    } else if (.Platform$OS.type == "unix") {
      plot.filename <- paste(train.pred.plot, ".pdf", sep = "", collapse = "")
      pdf(file = plot.filename, height = height, width = 40)
    }
  } else {
    if (.Platform$OS.type == "unix") {
      plot.filename <- paste(train.pred.plot, ".pdf", sep = "", collapse = "")
      pdf(file = plot.filename, height = height, width = 40)
    } else if (.Platform$OS.type == "windows") {
      plot.filename <- paste(train.pred.plot, ".pdf", sep = "", collapse = "")
      pdf(file = plot.filename, height = height, width = 40)
    }
  }

  nf <- layout(matrix(c(1, 2), 1, 2, byrow = T), widths = c(3, 1), heights = 1, respect = FALSE)

  #   par(mar = c(0, 0, 0, 0))
  #   MSIG.HeatMapPlot.3(V = t(prob.train), row.names = class.names, col.labels = class.labels.train, col.names =as.character(error.call), col.classes = class.names, phen.cmap = col[1:length(class.names)], main= "Train Samples Predictions", sub = " ", xlab=" ", ylab=" ", row.norm = F, cmap.type = 2)
  MSIG.HeatMapPlot.3(V = t(prob.train), row.names = class.names, col.labels = orig.class.labels.train, col.names = as.character(error.call), col.classes = class.names, phen.cmap = col[1:length(class.names)], main = "Train Samples Predictions", sub = " ", xlab = " ", ylab = " ", row.norm = F, cmap.type = 2)

  # legend

  leg.txt <- class.phen
  n.phen <- length(class.phen)
  p.vec <- symbs[1:n.phen]
  c.vec <- col[1:n.phen]

  par(mar = c(0, 0, 0, 0))
  plot(c(0, 0), c(1, 1), xlim = c(0, 1), ylim = c(0, 1), axes = F, type = "n", xlab = "", ylab = "")
  legend(x = 0, y = 1, legend = leg.txt, bty = "n", xjust = 0, yjust = 1, pch = p.vec, pt.bg = c.vec, col = "black", cex = symbol.scaling * 1.25, pt.cex = symbol.scaling * 2.5)

  if (non.interactive.run == F) {
    if (.Platform$OS.type == "windows") {
      savePlot(filename = train.pred.plot, type = "jpeg", device = dev.cur())
    } else if (.Platform$OS.type == "unix") {
      dev.off()
    }
  } else {
    dev.off()
  }

  print("Building SVM model completed. Predicting test data...")

  if (Ns - length(model.set) > 0) {
    test.pred <- predict(object = svm.model, newdata = t(m.test), decision.values = T, probability = T)
    dec.vals.test <- attr(test.pred, "decision.values")
    prob.test <- signif(attr(test.pred, "probabilities"), digits = 2)
    confidence.vector <- vector(length = n.test, mode = "numeric")
    bscore <- vector(length = n.test, mode = "numeric")
    max.k <- length(prob.train[1, ])
    random.pred.conf <- ((max.k - 1) / max.k)^2 + (max.k - 1) * (1 / max.k)^2
    for (ii in 1:n.test) {
      probs <- sort(prob.test[ii, ], decreasing = T)
      confidence.vector[ii] <- 1 - ((1 - probs[1])^2 + sum(probs[2:max.k]^2)) / random.pred.conf
      confidence.vector[ii] <- signif(confidence.vector[ii], digits = 3)
      if (class.list.test[ii] == as.character(test.pred[ii])) {
        bscore[ii] <- signif((1 - probs[1])^2, digits = 2)
      } else {
        bscore[ii] <- signif(probs[1]^2, digits = 2)
      }
    }

    confidence.call <- ifelse(confidence.vector >= high.conf.thres, " H ", " L ")
    error.call <- ifelse(class.list.test == as.character(test.pred), "   ", " * ")
    no.call <- ifelse(confidence.vector >= high.conf.thres, 0, 1)
    real.error <- ifelse(((no.call == 0) & (error.call == " * ")), 1, 0)
    correct.call <- ifelse(((no.call == 0) & (error.call == "   ")), 1, 0)
    col.symbols.test <- paste(confidence.call, error.call)
    class.names <- names(data.frame(prob.test))
    Brier.test <- signif(mean(bscore), digits = 2)

    test.results <- data.frame(cbind(as.character(sample.names.test), class.list.test, as.character(test.pred), error.call, confidence.call, confidence.vector, no.call, real.error, correct.call, prob.test, bscore))
    names(test.results)[1] <- "Test Sample Name"
    names(test.results)[2] <- "Actual"
    names(test.results)[3] <- "Predicted"
    names(test.results)[4] <- "Error (*)"
    names(test.results)[5] <- "Conf (H/L)"
    names(test.results)[6] <- "Conf"
    names(test.results)[7] <- "No Call"
    names(test.results)[8] <- "Real Error"
    names(test.results)[9] <- "Correct Call"

    names(test.results)[10 + length(class.phen.train)] <- "Brier score"
    #      print(test.results)
    print(c("Brier score (Test) = ", Brier.test))

    write("\n Test Results \n", file = prediction.results.file, append = T)
    write.table(test.results, file = prediction.results.file, append = T, quote = F, row.names = F, sep = "\t")

    write(c("\n\n Brier score (Test) = ", Brier.test), file = prediction.results.file, append = T)

    no.call.list <- split(no.call, class.list.test)
    real.error.list <- split(real.error, class.list.test)
    correct.call.list <- split(correct.call, class.list.test)
    count.class <- c(sapply(no.call.list, length), length(no.call))
    no.call.class <- c(sapply(no.call.list, sum), sum(no.call))
    real.error.class <- c(sapply(real.error.list, sum), sum(real.error))
    correct.call.class <- c(sapply(correct.call.list, sum), sum(correct.call))
    test.pred.high.conf <- ifelse(no.call == 0, as.character(test.pred), "-- no call")
    #   print(c("test.pred.high.conf =", test.pred.high.conf))

    no.call.class.pct <- no.call.class / count.class
    real.error.class.pct <- real.error.class / count.class
    correct.call.class.pct <- correct.call.class / count.class

    perf.table.test <- data.frame(cbind(c(names(no.call.list), "Total"), count.class, no.call.class, no.call.class.pct, real.error.class, real.error.class.pct, correct.call.class, correct.call.class.pct))
    names(perf.table.test) <- c("Class", "Count", "No Call", "No Call (%)", "Real Error", "Real Error (%)", "Correct Call", "Correct Call (%)")
    write.table(perf.table.test, file = prediction.results.file, append = T, quote = F, row.names = F, sep = "\t")
    print(perf.table.test)

    conf.table.test <- table(class.list.test, test.pred.high.conf)
    conf.table.test <- data.frame(cbind(row.names(conf.table.test), conf.table.test))
    print(conf.table.test)
    write("\n\n Confusion Matrix (Test) \n", file = prediction.results.file, append = T)
    write.table(conf.table.test, file = prediction.results.file, append = T, quote = F, row.names = F, sep = "\t")
    height <- ifelse(length(class.phen.train) > 50, 20, 0.30 * length(class.phen.train) + 10)

    if (non.interactive.run == F) {
      if (.Platform$OS.type == "windows") {
        plot.filename <- test.pred.plot
        x11(height = height, width = 40)
      } else if (.Platform$OS.type == "unix") {
        plot.filename <- paste(test.pred.plot, ".pdf", sep = "", collapse = "")
        pdf(file = plot.filename, height = height, width = 40)
      }
    } else {
      if (.Platform$OS.type == "unix") {
        plot.filename <- paste(test.pred.plot, ".pdf", sep = "", collapse = "")
        pdf(file = plot.filename, height = height, width = 40)
      } else if (.Platform$OS.type == "windows") {
        plot.filename <- paste(test.pred.plot, ".pdf", sep = "", collapse = "")
        pdf(file = plot.filename, height = height, width = 40)
      }
    }

    nf <- layout(matrix(c(1, 2), 1, 2, byrow = T), widths = c(5, 1), heights = 1, respect = FALSE)

    MSIG.HeatMapPlot.3(V = t(prob.test), row.names = names(test.results)[seq(10, 10 + length(class.phen.train) - 1)], col.labels = class.labels.test, col.names = as.character(error.call), col.classes = class.names, phen.cmap = col[1:length(class.names)], main = "Test Samples Predictions", sub = " ", xlab = " ", ylab = " ", row.norm = F, cmap.type = 2)


    # legend

    leg.txt <- class.phen
    n.phen <- length(class.phen)
    p.vec <- symbs[1:n.phen]
    c.vec <- col[1:n.phen]

    par(mar = c(0, 0, 0, 0))
    plot(c(0, 0), c(1, 1), xlim = c(0, 1), ylim = c(0, 1), axes = F, type = "n", xlab = "", ylab = "")
    legend(x = 0, y = 1, legend = leg.txt, bty = "n", xjust = 0, yjust = 1, pch = p.vec, pt.bg = c.vec, col = "black", cex = symbol.scaling * 1.25, pt.cex = symbol.scaling * 2.5)

    if (non.interactive.run == F) {
      if (.Platform$OS.type == "windows") {
        savePlot(filename = test.pred.plot, type = "jpeg", device = dev.cur())
      } else if (.Platform$OS.type == "unix") {
        dev.off()
      }
    } else {
      dev.off()
    }
  }


  # Save predictions for all classes in a gct and cls files

  #   print(c("dim train=", dim(prob.train)))
  #   print(c("dim test=", dim(prob.test)))

  if (Ns - length(model.set) > 0) {
    V <- cbind(t(prob.train), t(prob.test))
    W <- data.frame(V)
    names(W) <- c(as.character(sample.names.train), as.character(sample.names.test))
    row.names(W) <- names(train.results)[seq(7, 7 + length(class.phen.train) - 1)]
  } else {
    V <- t(prob.train)
    W <- data.frame(V)
    names(W) <- as.character(sample.names.train)
    row.names(W) <- names(train.results)[seq(7, 7 + length(class.phen.train) - 1)]
  }
  #   print(W)
  write.gct(gct.data.frame = W, descs = row.names(W), filename = prediction.matrix.file)

  print("Done predicting test data...")


  if (produce.contours == T) {
    #################################################################################
    # Contour plot results

    pca <- prcomp(t(m.train), retx = TRUE, center = T, scale. = T)
    S1 <- pca$x[, 1]
    S2 <- pca$x[, 2]
    S3 <- pca$x[, 3]
    X1 <- pca$rotation[, 1]
    X2 <- pca$rotation[, 2]
    X3 <- pca$rotation[, 3]

    row.mean <- apply(m.train, MARGIN = 1, FUN = mean)
    row.sd <- apply(m.train, MARGIN = 1, FUN = sd)

    # 2D plots

    c0 <- col
    #   c1 <- colors()[match(c0, colors())]
    c1 <- col
    color <- c1[class.labels]

    height <- 25
    if (non.interactive.run == F) {
      if (.Platform$OS.type == "windows") {
        plot.filename <- pred.2D.plot
        x11(height = height, width = 30)
      } else if (.Platform$OS.type == "unix") {
        plot.filename <- paste(pred.2D.plot, ".pdf", sep = "", collapse = "")
        pdf(file = plot.filename, height = height, width = 30)
      }
    } else {
      if (.Platform$OS.type == "unix") {
        plot.filename <- paste(pred.2D.plot, ".pdf", sep = "", collapse = "")
        pdf(file = plot.filename, height = height, width = 30)
      } else if (.Platform$OS.type == "windows") {
        plot.filename <- paste(pred.2D.plot, ".pdf", sep = "", collapse = "")
        pdf(file = plot.filename, height = height, width = 30)
      }
    }

    nf <- layout(matrix(c(1, 2), 1, 2, byrow = T), widths = c(3, 1), heights = 1, respect = FALSE)

    if (Ns - length(model.set) > 0) {
      test.scores <- predict(pca, t(m.test))
      S1 <- c(pca$x[, 1], test.scores[, 1])
      S2 <- c(pca$x[, 2], test.scores[, 2])
      S3 <- c(pca$x[, 3], test.scores[, 3])
    }

    max.S <- max(sqrt(S1 * S1 + S2 * S2))
    max.X <- max(sqrt(X1 * X1 + X2 * X2))
    X1 <- max.S * X1 / max.X
    X2 <- max.S * X2 / max.X
    X3 <- max.S * X3 / max.X
    max.A <- max(max.S, max.X)
    num.samples <- length(S1)

    plot(S1, S2, xlim = c(-max.A, max.A), ylim = c(-max.A, max.A), type = "n", main = "  ", sub = input.ds)

    for (j in 1:num.samples) {
      if (min(class.labels) == 0) {
        symb <- symbs[class.labels[j] + 1]
        color.code <- c1[class.labels[j] + 1]
      } else {
        symb <- symbs[class.labels[j]]
        color.code <- c1[class.labels[j]]
      }
      points(S1[j], S2[j], pch = symb, type = "p", cex = symbol.scaling * 3, bg = color.code, col = "black")
      #      if (j <= length(model.set)) {
      #         points(S1[j], S2[j], pch=22, type="p", cex = 2, bg = color.code, col = "black")
      #       } else {
      #         points(S1[j], S2[j], pch=21, type="p", cex = 1.5, bg = color.code, col = "black")
      #       }
    }

    for (j in 1:N) {
      x.coor <- X1[j] * 0.925
      y.coor <- X2[j] * 0.925
      arrows(0, 0, x.coor, y.coor, lwd = 1, length = 0.15, angle = 20, col = "grey50")
      if (use.feature.names == FALSE) {
        leg.txt <- paste("F", j, sep = "")
      } else {
        leg.txt <- gs.names[j]
      }
      text(X1[j], X2[j], labels = leg.txt, adj = NULL, pos = NULL, offset = 0.5, vfont = NULL, cex = 1, col = "grey50")
    }

    # compute grid

    points.axis <- 200
    x <- vector(length = points.axis, mode = "numeric")
    pca.x <- matrix(0, nrow = 2, ncol = points.axis * points.axis)

    for (i in 1:points.axis) {
      x[i] <- -max.A + i * 2 * max.A / points.axis
    }
    for (i in 1:points.axis) {
      for (j in 1:points.axis) {
        index.point <- i + (j - 1) * points.axis
        pca.x[1, index.point] <- -max.A + i * 2 * max.A / points.axis
        pca.x[2, index.point] <- -max.A + j * 2 * max.A / points.axis
      }
    }

    #   I <- pca$rotation %*% t(pca$x)
    #   I[1,1]*row.sd[1] + row.mean[1]

    grid.H <- pca$rotation[, 1:2] %*% pca.x
    for (i in 1:N) {
      grid.H[i, ] <- grid.H[i, ] * row.sd[i] + row.mean[i]
    }

    grid.pred <- predict(object = svm.model, newdata = t(grid.H), decision.values = T, probability = T)

    #   print(c("grid.pred:", grid.pred))

    prob.test <- attr(grid.pred, "probabilities")

    z <- matrix(0, nrow = points.axis, ncol = points.axis)
    z.class <- array(dim = c(length(class.phen.train), points.axis, points.axis))
    max.k <- length(prob.train[1, ])
    random.pred.conf <- ((max.k - 1) / max.k)^2 + (max.k - 1) * (1 / max.k)^2
    for (i in 1:points.axis) {
      for (j in 1:points.axis) {
        index.point <- i + (j - 1) * points.axis
        probs <- sort(prob.test[index.point, ], decreasing = T)
        z[i, j] <- 1 - ((1 - probs[1])^2 + sum(probs[2:max.k]^2)) / random.pred.conf
        for (k in 1:length(class.phen.train)) {
          if (probs[1] == prob.test[index.point, k]) {
            z.class[k, i, j] <- z[i, j]
          } else {
            z.class[k, i, j] <- 0
          }
        }

        #       probs <- sort(prob.test[index.point,], decreasing=T)
        #       z[i, j] <- probs[1] - probs[2]
        #       for (k in 1:length(class.phen.train)) {
        #          if (probs[1] == prob.test[index.point, k]) {
        #             z.class[k, i, j] <- z[i, j]
        #           } else {
        #             z.class[k, i, j] <- 0
        #           }
        #       }
      }
    }

    #   contour(x, x, z.class[3,,], nlevels = 50, col="black", lwd=2, add=T)

    library("RColorBrewer")
    if (length(levels) > 1) {
      contour(x, x, z, levels = levels, col = brewer.pal(n = 9, name = "Greys")[7], add = T)
    } else {
      contour(x, x, z, nlevels = nlevels, col = brewer.pal(n = 9, name = "Greys")[7], add = T)
    }
    contour(x, x, z, levels = 0.01, col = brewer.pal(n = 9, name = "Greys")[7], lwd = 2, add = T)
    for (k in 1:length(class.phen.train)) {
      contour(x, x, z.class[k, , ], levels = high.conf.thres, col = col[k], lwd = 2, add = T)
      #       contour(x, x, z.class[k,,], nlevels = 20, col=col[k], lwd=2, add=T)
    }

    #   z <- matrix(0, nrow=points.axis, ncol=points.axis)
    #     for (k in 1:length(class.phen.train)) {
    #     for (i in 1:points.axis) {
    #       for (j in 1:points.axis) {
    #        index.point <- i + (j - 1) * points.axis
    #         z[i, j] <- prob.test[index.point, k]
    #       }
    #     }
    #    contour(x, x, z, nlevels = 20, col=c1[k], add=T)
    #   }

    #   print("pca.x:")
    #   print(pca.x)

    #   for (j in seq(1, points.axis*points.axis)) {
    #      if (min(class.labels) == 0) {
    #          color.code <- c1[class.labels[match(grid.pred[j], class.list2)] + 1]
    #      } else {

    #          color.code <- c1[class.labels[match(grid.pred[j], class.list2)]]
    #      }
    #      points(pca.x[1, j], pca.x[2, j], pch=21, type="p", cex = 1, bg = color.code, col = color.code)
    #   }

    #   print(c("class.phen=", class.phen, " class.phen.train=", class.phen.train))
    #   print(c("prob.test:", prob.test))

    # legend

    leg.txt <- class.phen
    n.phen <- length(class.phen)
    #   p.vec <- c(rep(22, length(class.phen.train)), rep(21, n.phen - length(class.phen.train)))
    p.vec <- symbs[1:n.phen]
    c.vec <- c1[1:n.phen]

    par(mar = c(0, 0, 0, 0))
    plot(c(0, 0), c(1, 1), xlim = c(0, 1), ylim = c(0, 1), axes = F, type = "n", xlab = "", ylab = "")
    legend(x = 0, y = 1, legend = leg.txt, bty = "n", xjust = 0, yjust = 1, pch = p.vec, pt.bg = c.vec, col = "black", cex = symbol.scaling * 1.15, pt.cex = symbol.scaling * 3)

    if (non.interactive.run == F) {
      if (.Platform$OS.type == "windows") {
        savePlot(filename = pred.2D.plot, type = "jpeg", device = dev.cur())
      } else if (.Platform$OS.type == "unix") {
        dev.off()
      }
    } else {
      dev.off()
    }
  }
}


nnls.fit <- function(x, y, wsqrt = 1, eps = 0, rank.tol = 1e-07) {
  ## Purpose: Nonnegative Least Squares (similar to the S-Plus function
  ## with the same name) with the help of the R-library quadprog
  ## ------------------------------------------------------------------------
  ## Attention:
  ## - weights are square roots of usual weights
  ## - the constraint is coefficient>=eps
  ## ------------------------------------------------------------------------
  ## Author: Marcel Wolbers, July 99
  ##
  ## ========================================================================
  require("quadprog")
  m <- NCOL(x)
  if (length(eps) == 1) eps <- rep(eps, m)
  x <- x * wsqrt
  y <- y * wsqrt
  #  sometimes a rescaling of x and y helps (if solve.QP.compact fails otherwise)
  xscale <- apply(abs(x), 2, mean)
  yscale <- mean(abs(y))
  x <- t(t(x) / xscale)
  y <- y / yscale
  Rinv <- backsolve(qr.R(qr(x)), diag(m))
  cf <- solve.QP.compact(
    Dmat = Rinv, dvec = t(x) %*% y, Amat = rbind(rep(1, m)),
    Aind = rbind(rep(1, m), 1:m), bvec = eps * xscale / yscale,
    factorized = TRUE
  )$sol
  cf <- cf * yscale / xscale # scale back
  cf
}

MSIG.replace.row.col <- function(input.ds, output.ds, mode = "row", number, values) {
  dataset <- MSIG.Gct2Frame(filename = input.ds)
  m <- data.matrix(dataset$ds)


  if (mode == "row") {
    m[number, ] <- ifelse(length(values) == 1, rep(values, length(m[number, ])), values)
  } else if (mode == "col") {
    m[, number] <- ifelse(length(values) == 1, rep(values, length(m[, number])), values)
  } else {
    stop(c("unknown mode:", mode))
  }

  V <- data.frame(m)
  names(V) <- dataset$names
  row.names(V) <- dataset$row.names

  write.gct(gct.data.frame = V, descs = dataset$descs, filename = output.ds)
}


MSIG.apply.model <- function(
    gct.file,
    cls.file,
    phen.annot.file = NULL,
    output.dir,
    database.dir,
    identifiers,
    column.subset = "ALL",
    column.sel.type = "samples",
    thres = "NULL",
    ceil = "NULL",
    shift = "NULL",
    fold = 1,
    delta = 0,
    norm = 6,
    no.call.range.max = NULL,
    no.call.range.min = NULL) {
  ## Read test set

  print(c("Processing test file: ", gct.file))

  O <- MSIG.Subset.Dataset(
    input.ds = gct.file,
    input.cls = cls.file,
    column.subset = column.subset,
    column.sel.type = column.sel.type,
    row.subset = "ALL",
    output.ds = paste(output.dir, "temp2.gct", sep = ""),
    output.cls = paste(output.dir, "temp2.cls", sep = "")
  )

  O <- MSIG.Preprocess.Dataset(
    input.ds = paste(output.dir, "temp2.gct", sep = ""),
    output.ds = paste(output.dir, "temp3.gct", sep = ""),
    thres = thres,
    ceil = ceil,
    normalization = "NULL"
  )

  dataset <- MSIG.Gct2Frame(filename = paste(output.dir, "temp3.gct", sep = ""))
  m.test <- data.matrix(dataset$ds)
  gs.names.test <- dataset$row.names
  gs.descs.test <- dataset$descs
  sample.names.test <- dataset$names

  Ns.test <- length(m.test[1, ])
  Ng.test <- length(m.test[, 1])

  CLS <- MSIG.ReadClsFile(file = paste(output.dir, "temp2.cls", sep = ""))
  class.labels.test <- CLS$class.v
  class.phen.test <- CLS$phen
  class.list.test <- CLS$class.list

  # loop over signatures

  for (sig in identifiers) {
    # Read parameters and signatures (gmt, gct and cls files)

    filename <- paste(database.dir, sig, ".msig.params", sep = "")
    temp <- readLines(filename)
    seed <- as.numeric(noquote(unlist(strsplit(temp[[1]], "\t")))[2])
    topgs <- as.numeric(noquote(unlist(strsplit(temp[[2]], "\t")))[2])
    link.function <- unlist(strsplit(temp[[3]], "\t"))[2]
    model.type <- unlist(strsplit(temp[[4]], "\t"))[2]
    burnin.iter <- as.numeric(noquote(unlist(strsplit(temp[[5]], "\t")))[2])
    mcmc.iter <- as.numeric(noquote(unlist(strsplit(temp[[6]], "\t")))[2])
    col.target <- unlist(strsplit(temp[[7]], "\t"))[2]
    col.control <- unlist(strsplit(temp[[8]], "\t"))[2]
    no.call.r.max <- as.numeric(noquote(unlist(strsplit(temp[[9]], "\t")))[2])
    no.call.r.min <- as.numeric(noquote(unlist(strsplit(temp[[10]], "\t")))[2])
    beta0.train <- as.numeric(noquote(unlist(strsplit(temp[[11]], "\t")))[2])
    beta1.train <- as.numeric(noquote(unlist(strsplit(temp[[12]], "\t")))[2])
    target.class <- unlist(strsplit(temp[[13]], "\t"))[2]

    c1 <- c(col.target, col.control) # color for target class is always first

    if (is.null(no.call.range.max)) {
      no.call.range.max <- no.call.r.max
    }
    if (is.null(no.call.range.min)) {
      no.call.range.min <- no.call.r.min
    }

    filename <- paste(database.dir, sig, ".msig.gct", sep = "")
    dataset <- MSIG.Gct2Frame(filename = filename)
    sample.molsig.sorted.subset <- dataset$ds
    Ns <- length(sample.molsig.sorted.subset[1, ])
    msize.all <- length(sample.molsig.sorted.subset[, 1])
    sample.molsig.sorted.subset.gs <- dataset$row.names
    sample.names <- dataset$names

    filename <- paste(database.dir, sig, ".msig.gct", sep = "")
    dataset <- MSIG.Gct2Frame(filename = filename)
    sample.molsig.sorted.subset <- dataset$ds
    Ns <- length(sample.molsig.sorted.subset[1, ])
    msize.all <- length(sample.molsig.sorted.subset[, 1])
    sample.molsig.sorted.subset.gs <- dataset$row.names
    sample.names <- dataset$names

    filename <- paste(database.dir, sig, ".msig.cls", sep = "")
    CLS <- MSIG.ReadClsFile(file = filename)
    class.labels <- CLS$class.v
    class.phen <- CLS$phen
    class.list <- CLS$class.list

    # rename phenotypes according to target class and set everything else as control

    for (i in 1:length(class.list)) {
      if (class.list[i] == target.class) {
        class.labels[i] <- 1
      } else {
        class.list[i] <- "CNTL"
        class.labels[i] <- 0
      }
    }

    print(c("Target class:", target.class))
    print(c("Class labels:", class.labels))

    col.index <- order(class.labels, decreasing = T)
    for (j in 1:msize.all) {
      sample.molsig.sorted.subset[j, ] <- sample.molsig.sorted.subset[j, col.index]
    }
    sample.names <- sample.names[col.index]
    class.labels <- class.labels[col.index]
    class.list <- class.list[col.index]
    class.phen <- c(target.class, "CNTL")
    control.class <- "CNTL"

    # match test set and model

    gs.names2 <- intersect(sample.molsig.sorted.subset.gs, gs.names.test)
    locations <- match(gs.names2, gs.names.test, nomatch = 0)
    m.test2 <- m.test[locations, ]
    locations2 <- match(gs.names2, sample.molsig.sorted.subset.gs)
    m.train <- sample.molsig.sorted.subset[locations2, ]

    print(c("Matched signature and test set: overlap=", length(gs.names2), " Total original signature size= ", length(sample.molsig.sorted.subset.gs)))

    # define signature model

    msize <- length(locations)
    sig.matrix <- array(0, dim = c(msize, Ns))
    sig.matrix.test <- array(0, dim = c(msize, Ns.test))
    for (k in 1:Ns) {
      sig.matrix[, k] <- rank(m.train[, k], ties.method = "average")
    }
    for (k in 1:Ns.test) {
      sig.matrix.test[, k] <- rank(m.test2[, k], ties.method = "average")
    }

    sig.matrix.all <- cbind(sig.matrix, sig.matrix.test)
    sample.names.all <- c(sample.names, sample.names.test)

    MSIG.HeatMapPlot.5(V = t(sig.matrix.all), row.names = sample.names.all, col.labels = rep(1, msize), col.classes = "C", col.names = gs.names2, main = paste(sig, gct.file, sep = " / "), xlab = " ", ylab = " ", row.norm = F, cmap.type = 2)

    t.class.point <- apply(sig.matrix[, class.list == target.class], MARGIN = 1, FUN = mean)
    c.class.point <- apply(sig.matrix[, class.list == control.class], MARGIN = 1, FUN = mean)

    d.t.class <- vector(length = Ns, mode = "numeric")
    d.c.class <- vector(length = Ns, mode = "numeric")
    d.c.t.class <- sum(abs(t.class.point - c.class.point))
    x <- vector(length = Ns, mode = "numeric")
    y <- vector(length = Ns, mode = "numeric")

    d.t.class.test <- vector(length = Ns.test, mode = "numeric")
    d.c.class.test <- vector(length = Ns.test, mode = "numeric")
    x.test <- vector(length = Ns.test, mode = "numeric")
    y.test <- vector(length = Ns.test, mode = "numeric")

    for (i in 1:Ns) {
      d.t.class[i] <- sum(abs(t.class.point - sig.matrix[, i])) / d.c.t.class
      d.c.class[i] <- sum(abs(c.class.point - sig.matrix[, i])) / d.c.t.class
      x[i] <- (d.t.class[i]^2 - d.c.class[i]^2 - 1) / (-2)
      y[i] <- sqrt(d.c.class[i]^2 - x[i]^2)
    }

    # Create regression model using x as the input variable

    print(c("Creating regression signature model using overlap..."))

    target.var <- ifelse(class.list == target.class, 1, 0)

    if (model.type == "Bayesian") {
      if (link.function == "logit") {
        reg.model <- MCMClogit(target.var ~ x, burnin = burnin.iter, mcmc = mcmc.iter, bayes.resid = T) # Logit
      } else if (link.function == "probit") {
        reg.model <- MCMCprobit(target.var ~ x, burnin = burnin.iter, mcmc = mcmc.iter, bayes.resid = T) # Probit
      } else {
        stop("Unknown link function")
      }
    } else if (model.type == "Classic") {
      if (link.function == "logit") {
        reg.model <- glm(target.var ~ x, family = binomial("logit")) # Logit
      } else if (link.function == "probit") {
        reg.model <- glm(target.var ~ x, family = binomial("probit")) # Probit
      } else {
        stop("Unknown link function")
      }
    } else {
      stop("Unknown model type")
    }

    if (model.type == "Bayesian") {
      beta0 <- reg.model[, 1]
      beta1 <- reg.model[, 2]
      print(c("beta0=", beta0, " beta1=", beta1))
      prob.i <- matrix(0, nrow = Ns, ncol = 3)
    } else if (model.type == "Classic") {
      beta0 <- reg.model[[1]][1]
      beta1 <- reg.model[[1]][2]
      print(c("beta0=", beta0, " beta1=", beta1))
      prob.i <- matrix(0, nrow = Ns, ncol = 3)
    } else {
      stop("Unknown model type")
    }

    print(c("beta0 train=", beta0.train, " beta0=", beta0))
    print(c("beta1 train=", beta1.train, " beta1=", beta1))

    xmin <- min(x)
    xmax <- max(x)
    range.x <- xmax - xmin
    prob.m <- matrix(0, nrow = 1000, ncol = 3)
    x.m <- vector(length = 1000, mode = "numeric")
    for (k in 1:1000) {
      x.m[k] <- xmin + k * (range.x / 1000)
      if (link.function == "logit") {
        p.vec <- (exp(beta0 + beta1 * x.m[k]) / (1 + exp(beta0 + beta1 * x.m[k]))) # Logit
      } else if (link.function == "probit") {
        p.vec <- (erf(beta0 + beta1 * x.m[k]) + 1) / 2 # Probit
      } else {
        nstop("Unknown link function")
      }
      prob.m[k, 1] <- quantile(p.vec, probs = 0.5)
      prob.m[k, 2] <- quantile(p.vec, probs = 0.05)
      prob.m[k, 3] <- quantile(p.vec, probs = 0.95)
    }
    istar <- which.min(abs(0.5 - prob.m[, 1]))
    istar <- xmin + istar * (range.x / 1000)

    for (i in 1:Ns.test) {
      d.t.class.test[i] <- sum(abs(t.class.point - sig.matrix.test[, i])) / d.c.t.class
      d.c.class.test[i] <- sum(abs(c.class.point - sig.matrix.test[, i])) / d.c.t.class
      x.test[i] <- (d.t.class.test[i]^2 - d.c.class.test[i]^2 - 1) / (-2)
      y.test[i] <- sqrt(d.c.class.test[i]^2 - x.test[i]^2)
    }

    # plot results for this test dataset

    x.range <- range(c(x, x.test, 0, 1))
    y.range <- range(c(y, y.test, 0))

    x11(height = 24, width = 30)
    plot(x, y, xlim = x.range, ylim = y.range, type = "n", main = sig, sub = gct.file)
    points(0, 0, cex = 2, pch = 21, col = 1, bg = 3)
    points(1, 0, cex = 2, pch = 21, col = 1, bg = 2)
    points(x[class.list == control.class], y[class.list == control.class], cex = 1, pch = 21, col = 1, bg = 3)
    points(x[class.list == target.class], y[class.list == target.class], cex = 1, pch = 21, col = 1, bg = 2)
    k <- 1
    for (i in class.list.test) {
      points(x.test[class.list.test == i], y.test[class.list.test == i], cex = 1, pch = 22, col = 1, bg = k %% 5)
      k <- k + 1
    }

    prob.i.test <- matrix(0, nrow = Ns.test, ncol = 3)

    for (i in 1:Ns.test) {
      if (link.function == "logit") {
        p.vec.test <- (exp(beta0 + beta1 * x.test[i]) / (1 + exp(beta0 + beta1 * x.test[i]))) # Logit
      } else if (link.function == "probit") {
        p.vec.test <- (erf(beta0 + beta1 * x.test[i]) + 1) / 2 # Probit
      } else {
        stop("Unknown link function")
      }
      prob.i.test[i, 1] <- quantile(p.vec.test, probs = 0.5)
      prob.i.test[i, 2] <- quantile(p.vec.test, probs = 0.05)
      prob.i.test[i, 3] <- quantile(p.vec.test, probs = 0.95)
    }

    x.index <- order(x.test, decreasing = F)
    x.order.test <- x.test[x.index]
    prob.i.order.test <- prob.i.test[x.index, ]
    class.list.test.order <- class.list.test[x.index]

    x11(height = 7, width = 9.5)
    nf <- layout(matrix(c(1, 2), 1, 2, byrow = T), widths = c(3.75, 1), heights = 1, respect = FALSE)

    plot(x.order.test, prob.i.order.test[, 1], sub = gct.file, pch = 20, ylim = c(-0.05, 1.07), main = sig, xlim = c(-0.1, 1.1), col = 0, cex.axis = 1.35, cex = 3, cex.lab = 1.35, xlab = "Activation Index", ylab = "Probability")
    points(x.m, prob.m[, 1], type = "l", lwd = 2, col = 1, lty = 1, cex = 1)
    points(x.m, prob.m[, 2], type = "l", col = 4, lty = 1, cex = 1)
    points(x.m, prob.m[, 3], type = "l", col = 4, lty = 1, cex = 1)

    arrows(x.order.test, prob.i.order.test[, 2], x.order.test, prob.i.order.test[, 3], col = 4, angle = 90, code = 3, length = 0.0)
    range.x <- range(x.order.test)
    points(range.x, c(0.5, 0.5), type = "l", lty = 3, col = 1, lwd = 2)
    points(c(istar, istar), c(-0.07, 1.07), type = "l", lty = 3, col = 1, lwd = 2)
    k <- 1
    for (i in class.list.test) {
      points(x.order.test[class.list.test.order == i], prob.i.order.test[class.list.test.order == i, 1], pch = 21, bg = k %% 5, col = 1, cex = 2)
      k <- k + 1
    }
    leg.txt <- unique(class.list.test.order)
    p.vec <- rep(21, length(unique(class.list.test.order)))
    c.vec <- rep(seq(1, 5), length(unique(class.list.test.order)))
    par(mar = c(0, 0, 0, 0))
    plot(c(0, 0), c(1, 1), xlim = c(0, 1), ylim = c(0, 1), axes = F, type = "n", xlab = "", ylab = "")
    legend(x = 0, y = 0.8, legend = leg.txt, bty = "n", xjust = 0, yjust = 1, pch = p.vec, pt.bg = c.vec, col = "black", cex = 1.2, pt.cex = 2)

    activation.indicator <- ifelse(prob.i.test[, 1] >= 0.5, 1, 0)
    activation.indicator <- ifelse((prob.i.test[, 1] >= no.call.range.max) | (prob.i.test[, 1] <= no.call.range.min), activation.indicator, 0.5)

    # read phenotype annotation and append them to test set results

    if (!is.null(phen.annot.file)) {
      filename <- phen.annot.file
      dataset <- MSIG.Gct2Frame(filename = filename)
      phen.annot <- data.matrix(dataset$ds)
      phen.annot.gs <- dataset$row.names

      for (i in 1:length(phen.annot[, 1])) {
        phen.annot[i, ] <- (phen.annot[i, ] - min(phen.annot[i, ])) / (max(phen.annot[i, ]) - min(phen.annot[i, ]))
      }

      z <- rbind(prob.i.test[, 1], activation.indicator, phen.annot)
      p.lab <- c(paste("P(", sig, ")", sep = ""), paste("A(", sig, ")", sep = ""), phen.annot.gs)
    } else {
      z <- rbind(prob.i.test[, 1], activation.indicator)
      p.lab <- c(paste("P(", sig, ")", sep = ""), paste("A(", sig, ")", sep = ""))
    }

    MSIG.HeatMapPlot.5(V = z, row.names = p.lab, col.labels = class.labels.test, col.classes = class.phen.test, col.names = sample.names.test, main = paste(sig, " Activation on Test", sep = ""), xlab = " ", ylab = " ", sub = gct.file, row.norm = F, cmap.type = 3, rotated.col.labels = T)

    if (sig == identifiers[[1]]) {
      z.all <- prob.i.test[, 1]
      z.act.all <- activation.indicator
      if (!is.null(phen.annot.file)) {
        phen.annot.all <- phen.annot
        phen.annot.gs.all <- phen.annot.gs
      }
      p.lab.all <- paste("P(", sig, ")", sep = "")
      p.act.lab.all <- paste("A(", sig, ")", sep = "")
    } else {
      z.all <- rbind(z.all, prob.i.test[, 1])
      z.act.all <- rbind(z.act.all, activation.indicator)
      p.lab.all <- c(p.lab.all, paste("P(", sig, ")", sep = ""))
      p.act.lab.all <- c(p.act.lab.all, paste("A(", sig, ")", sep = ""))
    }
  }

  if (!is.null(phen.annot.file)) {
    z.all <- rbind(z.all, phen.annot)
    z.act.all <- rbind(z.act.all, phen.annot)
    p.lab.all <- c(p.lab.all, phen.annot.gs)
    p.act.lab.all <- c(p.act.lab.all, phen.annot.gs)
  }

  print(c("dim z.all=", dim(z.all)))

  MSIG.HeatMapPlot.5(V = z.all, row.names = p.lab.all, col.labels = class.labels.test, col.classes = class.phen.test, col.names = sample.names.test, main = " ", xlab = " ", ylab = " ", sub = gct.file, row.norm = F, cmap.type = 2, rotated.col.labels = T)

  MSIG.HeatMapPlot.5(V = z.act.all, row.names = p.act.lab.all, col.labels = class.labels.test, col.classes = class.phen.test, col.names = sample.names.test, main = " ", xlab = " ", ylab = " ", sub = gct.file, row.norm = F, cmap.type = 2, rotated.col.labels = T)
}


MSIG.Create.model <- function(
    gct.file,
    cls.file,
    output.dir,
    database.dir,
    identifier,
    target.class,
    column.subset = "ALL",
    column.sel.type = "samples",
    thres = "NULL",
    ceil = "NULL",
    shift = "NULL",
    fold = 1,
    delta = 0,
    norm = 6,
    seed = 1234,
    topgs = 25,
    link.function = "logit", #  "probit" or "logit"
    model.type = "Classic", #  "Bayesian" or "Classic"
    burnin.iter = 5000,
    mcmc.iter = 25000,
    col.target = "darkgreen",
    col.control = "yellow",
    no.call.range.min = 0.3,
    no.call.range.max = 0.7) {
  c1 <- c(col.target, col.control) # color for target class is always first

  # start of methodology
  set.seed(seed)

  # Re-order/subset dataset

  O <- MSIG.Subset.Dataset(
    input.ds = gct.file,
    input.cls = cls.file,
    column.subset = column.subset,
    column.sel.type = column.sel.type,
    row.subset = "ALL",
    output.ds = paste(output.dir, "temp1.gct", sep = ""),
    output.cls = paste(output.dir, "temp1.cls", sep = "")
  )

  O <- MSIG.Preprocess.Dataset(
    input.ds = paste(output.dir, "temp1.gct", sep = ""),
    output.ds = paste(output.dir, "temp2.gct", sep = ""),
    thres = thres,
    ceil = ceil,
    fold = fold,
    delta = delta,
    normalization = norm
  )

  # Read dataset

  dataset <- MSIG.Gct2Frame(filename = paste(output.dir, "temp2.gct", sep = ""))
  m <- data.matrix(dataset$ds)
  gs.names <- dataset$row.names
  gs.descs <- dataset$descs
  sample.names <- dataset$names

  dim(m)

  Ns <- length(m[1, ])
  Ng <- length(m[, 1])

  CLS <- MSIG.ReadClsFile(file = paste(output.dir, "temp1.cls", sep = ""))
  class.labels <- CLS$class.v
  class.phen <- CLS$phen
  class.list <- CLS$class.list
  if (is.na(match(target.class, class.phen))) stop(c("target class is not phenotype in:", model.cls))

  # rename phenotypes according to target class and set everything else as control

  print("Renaming phenotypes...")

  for (i in 1:length(class.list)) {
    if (class.list[i] == target.class) {
      class.labels[i] <- 1
    } else {
      class.list[i] <- "CNTL"
      class.labels[i] <- 0
    }
  }

  col.index <- order(class.labels, decreasing = T)
  for (j in 1:Ng) {
    m[j, ] <- m[j, col.index]
  }
  sample.names <- sample.names[col.index]
  class.labels <- class.labels[col.index]
  class.list <- class.list[col.index]
  class.phen <- c(target.class, "CNTL")
  control.class <- "CNTL"

  ##### Marker selection

  print("Executing marker selection...")

  topgs <- ifelse(topgs > floor(Ng / length(class.phen)), floor(Ng / length(class.phen)), topgs)
  sample.molsig.sorted.subset <- matrix(0, nrow = length(class.phen) * topgs, ncol = Ns)
  sample.molsig.sorted.subset.gs <- vector(length = length(class.phen) * topgs, mode = "character")
  sample.molsig.sorted.s2n <- vector(length = length(class.phen) * topgs, mode = "character")
  sample.molsig.sorted.class <- vector(length = length(class.phen) * topgs, mode = "character")
  sample.molsig.sorted.index <- vector(length = length(class.phen) * topgs, mode = "character")

  num.k <- 1
  for (k in class.phen) {
    print(c("Executing marker selection for class:", k))

    class.k.labels <- ifelse(class.list == k, 0, 1)
    col.index <- order(class.k.labels, decreasing = F)
    m1 <- m
    for (j in 1:Ng) {
      m1[j, ] <- m[j, col.index]
    }
    names(m1) <- sample.names[col.index]
    class.k.labels <- class.k.labels[col.index]

    O <- GSEA.GeneRanking(m1, class.k.labels, gs.names, 1, permutation.type = 0, sigma.correction = "GeneCluster", fraction = 1, replace = F, reverse.sign = F)
    order.matrix <- O$order.matrix
    obs.order.matrix <- O$obs.order.matrix
    correl.matrix <- O$s2n.matrix
    obs.correl.matrix <- O$obs.s2n.matrix

    rm(O)

    obs.s2n.orig <- apply(obs.correl.matrix, 1, median) # using median to assign enrichment scores
    obs.index <- order(obs.s2n.orig, decreasing = T)
    obs.s2n <- sort(obs.s2n.orig, decreasing = T)
    sample.molsig.sorted <- m[obs.index, ]
    gs.names.sorted <- gs.names[obs.index]

    start <- (num.k - 1) * topgs + 1
    end <- num.k * topgs
    sample.molsig.sorted.subset[start:end, ] <- sample.molsig.sorted[1:topgs, ]
    sample.molsig.sorted.subset.gs[start:end] <- gs.names.sorted[1:topgs]
    sample.molsig.sorted.s2n[start:end] <- signif(obs.s2n[1:topgs], digits = 3)
    sample.molsig.sorted.class[start:end] <- class.phen[num.k]
    sample.molsig.sorted.index[start:end] <- obs.index[1:topgs]

    num.k <- num.k + 1
  }

  #    x11(height = 24, width = 30)
  #   nf <- layout(matrix(c(1,2), 1, 2, byrow=T), widths = c(6, 1), heights = 1, respect = FALSE)

  c1 <- c("grey3", "cadetblue", "darkgreen", "chartreuse2", "red1", "darkred", "orange", "blue2", "lightblue", "pink", "coral")

  print("Making heat maps...")

  MSIG.HeatMapPlot.5(V = t(sample.molsig.sorted.subset), row.names = sample.names, col.labels = c(rep(1, topgs), rep(0, topgs)), col.classes = c("C1", "C2"), col.names = sample.molsig.sorted.subset.gs, main = "Original Expression (Norm. Signature genes)", xlab = " ", ylab = " ", sub = gct.file, row.norm = T, cmap.type = 4, rotated.col.labels = T)

  #   MSIG.HeatMapPlot.5(V = t(log(sample.molsig.sorted.subset)), row.names = sample.names, col.labels = c(rep(1, topgs), rep(0, topgs)), col.classes = c("C1", "C2"), col.names = sample.molsig.sorted.subset.gs, main = "Original Expression (Signature genes)", xlab=" ", ylab=" ", sub = gct.file, row.norm = F,  cmap.type = 2, rotated.col.labels = T)

  msize <- 2 * topgs
  sig.matrix <- array(0, dim = c(msize, Ns))
  for (k in 1:Ns) {
    sig.matrix[, k] <- rank(sample.molsig.sorted.subset[, k], ties.method = "average")
  }

  MSIG.HeatMapPlot.5(V = t(sig.matrix), row.names = sample.names, col.labels = c(rep(1, topgs), rep(0, topgs)), col.classes = c("C1", "C2"), col.names = sample.molsig.sorted.subset.gs, main = "Signature -- Training", xlab = " ", ylab = " ", sub = gct.file, row.norm = F, cmap.type = 2, rotated.col.labels = T)

  print("Defining signature...")

  t.class.point <- apply(sig.matrix[, class.list == target.class], MARGIN = 1, FUN = mean)
  c.class.point <- apply(sig.matrix[, class.list == control.class], MARGIN = 1, FUN = mean)

  d.t.class <- vector(length = Ns, mode = "numeric")
  d.c.class <- vector(length = Ns, mode = "numeric")
  d.c.t.class <- sum(abs(t.class.point - c.class.point))
  x <- vector(length = Ns, mode = "numeric")
  y <- vector(length = Ns, mode = "numeric")

  for (i in 1:Ns) {
    d.t.class[i] <- sum(abs(t.class.point - sig.matrix[, i])) / d.c.t.class
    d.c.class[i] <- sum(abs(c.class.point - sig.matrix[, i])) / d.c.t.class
    x[i] <- (d.t.class[i]^2 - d.c.class[i]^2 - 1) / (-2)
    y[i] <- sqrt(d.c.class[i]^2 - x[i]^2)
  }

  x.range <- range(c(x, 0, 1))
  y.range <- range(c(y, 0))

  x11(height = 24, width = 30)
  plot(x, y, xlim = x.range, ylim = y.range, type = "n")
  points(0, 0, cex = 2, pch = 21, col = 1, bg = 3)
  points(1, 0, cex = 2, pch = 21, col = 1, bg = 2)
  points(x[class.list == control.class], y[class.list == control.class], cex = 1, pch = 21, col = 1, bg = 3)
  points(x[class.list == target.class], y[class.list == target.class], cex = 1, pch = 21, col = 1, bg = 2)

  # Create regression model using x as the input variable

  target.var <- ifelse(class.list == target.class, 1, 0)

  if (model.type == "Bayesian") {
    if (link.function == "logit") {
      reg.model <- MCMClogit(target.var ~ x, burnin = burnin.iter, mcmc = mcmc.iter, bayes.resid = T) # Logit
    } else if (link.function == "probit") {
      reg.model <- MCMCprobit(target.var ~ x, burnin = burnin.iter, mcmc = mcmc.iter, bayes.resid = T) # Probit
    } else {
      stop("Unknown link function")
    }
  } else if (model.type == "Classic") {
    if (link.function == "logit") {
      reg.model <- glm(target.var ~ x, family = binomial("logit")) # Logit
    } else if (link.function == "probit") {
      reg.model <- glm(target.var ~ x, family = binomial("probit")) # Probit
    } else {
      stop("Unknown link function")
    }
  } else {
    stop("Unknown model type")
  }

  if (model.type == "Bayesian") {
    beta0 <- reg.model[, 1]
    beta1 <- reg.model[, 2]
    print(c("beta0=", beta0, " beta1=", beta1))
    prob.i <- matrix(0, nrow = Ns, ncol = 3)
  } else if (model.type == "Classic") {
    beta0 <- reg.model[[1]][1]
    beta1 <- reg.model[[1]][2]
    print(c("beta0=", beta0, " beta1=", beta1))
    prob.i <- matrix(0, nrow = Ns, ncol = 3)
  } else {
    stop("Unknown model type")
  }

  for (i in 1:Ns) {
    if (link.function == "logit") {
      p.vec <- (exp(beta0 + beta1 * x[i]) / (1 + exp(beta0 + beta1 * x[i]))) # Logit
    } else if (link.function == "probit") {
      p.vec <- (erf(beta0 + beta1 * x[i]) + 1) / 2 # Probit
    } else {
      stop("Unknown link function")
    }
    prob.i[i, 1] <- quantile(p.vec, probs = 0.5)
    prob.i[i, 2] <- quantile(p.vec, probs = 0.05)
    prob.i[i, 3] <- quantile(p.vec, probs = 0.95)
  }
  xmin <- min(x)
  xmax <- max(x)
  range.x <- xmax - xmin
  prob.m <- matrix(0, nrow = 1000, ncol = 3)
  x.m <- vector(length = 1000, mode = "numeric")
  for (k in 1:1000) {
    x.m[k] <- xmin + k * (range.x / 1000)
    if (link.function == "logit") {
      p.vec <- (exp(beta0 + beta1 * x.m[k]) / (1 + exp(beta0 + beta1 * x.m[k]))) # Logit
    } else if (link.function == "probit") {
      p.vec <- (erf(beta0 + beta1 * x.m[k]) + 1) / 2 # Probit
    } else {
      stop("Unknown link function")
    }
    prob.m[k, 1] <- quantile(p.vec, probs = 0.5)
    prob.m[k, 2] <- quantile(p.vec, probs = 0.05)
    prob.m[k, 3] <- quantile(p.vec, probs = 0.95)
  }
  istar <- which.min(abs(0.5 - prob.m[, 1]))
  istar <- xmin + istar * (range.x / 1000)

  # save prob.m in file

  # V <- data.frame(rbind(x.m, prob.m[,1], prob.m[,2], prob.m[,3]))
  # names(V) <- seq(1, length(x.m))
  # row.names(V) <- c("x.m", "prob.m_0.5", "prob.m_0.05", "prob.m_0.95")
  # descs <- c("x.m", "prob.m_0.5", "prob.m_0.05", "prob.m_0.95")
  # write.gct(gct.data.frame = V, descs = descs, filename = paste(data.dir, file.iden, ".train.set.prob.model.gct", sep=""))

  x.index <- order(x, decreasing = F)
  x.order <- x[x.index]
  prob.i.order <- prob.i[x.index, ]
  target.var.order <- target.var[x.index]
  target.var.order <- ifelse(target.var.order == 1, c1[1], c1[2])
  target.var <- ifelse(target.var == 1, c1[1], c1[2])

  # save probabilities in file

  # V <- data.frame(rbind(class.labels, z.vector, rep(istar, length(z.vector)), prob.i[,1], prob.i[,2], prob.i[,3]))
  # names(V) <- sample.names
  # row.names(V) <- c("class", "zscore", "istar", "prob_0.5", "prob_0.05", "prob_0.95")
  # descs <- c("class", "zscore", "istar", "prob_0.5", "prob_0.05", "prob_0.95")
  # write.gct(gct.data.frame = V, descs = descs, filename = paste(data.dir, file.iden, ".train.set.probs.gct", sep=""))

  # Plot posterior probabilities

  x11(height = 7, width = 9.5)
  nf <- layout(matrix(c(1, 2), 1, 2, byrow = T), widths = c(3.75, 1), heights = 1, respect = FALSE)

  plot(x.order, prob.i.order[, 1], sub = gct.file, pch = 20, ylim = c(-0.2, 1.07), xlim = c(-0.1, 1.1), col = 0, cex = 2, xlab = "Activation Index", ylab = "Probability")
  points(x.m, prob.m[, 1], type = "l", lwd = 2, col = 1, lty = 1, cex = 1)
  points(x.m, prob.m[, 2], type = "l", col = 4, lty = 1, cex = 1)
  points(x.m, prob.m[, 3], type = "l", col = 4, lty = 1, cex = 1)
  arrows(x.order, prob.i.order[, 2], x.order, prob.i.order[, 3], col = 4, angle = 90, code = 3, length = 0.0)

  range.x <- range(x.order)
  points(range.x, c(0.5, 0.5), type = "l", lty = 3, col = 1, lwd = 2)
  points(range.x, c(-.15, -0.15), type = "l", lty = 1, col = 1, lwd = 2)
  points(c(istar, istar), c(-0.07, 1.07), type = "l", lty = 3, col = 1, lwd = 2)
  points(x.order, prob.i.order[, 1], pch = 21, bg = target.var.order, col = 1, cex = 2)
  points(x.order, rep(-0.15, length(x.order)), pch = 21, bg = target.var.order, col = 1, cex = 2)

  leg.txt <- class.phen
  p.vec <- rep(21, 21)
  c.vec <- c1
  par(mar = c(0, 0, 0, 0))
  plot(c(0, 0), c(1, 1), xlim = c(0, 1), ylim = c(0, 1), axes = F, type = "n", xlab = "", ylab = "")
  legend(x = 0, y = 0.8, legend = leg.txt, bty = "n", xjust = 0, yjust = 1, pch = p.vec, pt.bg = c.vec, col = "black", cex = 1.2, pt.cex = 2)

  activation.indicator <- ifelse(prob.i[, 1] >= 0.5, 1, 0)
  activation.indicator <- ifelse((prob.i[, 1] >= no.call.range.max) | (prob.i[, 1] <= no.call.range.min), activation.indicator, 0.5)
  z <- rbind(prob.i[, 1], activation.indicator)

  print("z:")
  print(z)
  print("class.labels:")
  print(class.labels)
  print("class.phen:")
  print(class.phen)
  print("sample.names:")
  print(sample.names)

  MSIG.HeatMapPlot.5(V = z, row.names = c("P(A)", "A"), col.labels = class.labels, col.classes = class.phen, col.names = sample.names, main = "Activation Assessment (training set)", xlab = " ", ylab = " ", sub = gct.file, row.norm = F, cmap.type = 2)

  # Save parameters and signature as a gene set file plus a gct with the numerical values and a cls with the class vector
  filename <- paste(database.dir, identifier, ".msig.params", sep = "")
  write(paste("seed", seed, sep = "\t"), file = filename, append = F, ncolumns = 3)
  write(paste("topgs", topgs, sep = "\t"), file = filename, append = T, ncolumns = 3)
  write(paste("link.function", link.function, sep = "\t"), file = filename, append = T, ncolumns = 3)
  write(paste("model.type", model.type, sep = "\t"), file = filename, append = T, ncolumns = 3)
  write(paste("burnin.iter", burnin.iter, sep = "\t"), file = filename, append = T, ncolumns = 3)
  write(paste("mcmc.iter", mcmc.iter, sep = "\t"), file = filename, append = T, ncolumns = 3)
  write(paste("col.target", col.target, sep = "\t"), file = filename, append = T, ncolumns = 3)
  write(paste("col.control", col.control, sep = "\t"), file = filename, append = T, ncolumns = 3)
  write(paste("no.call.range.max", no.call.range.max, sep = "\t"), file = filename, append = T, ncolumns = 3)
  write(paste("no.call.range.min", no.call.range.min, sep = "\t"), file = filename, append = T, ncolumns = 3)
  write(paste("beta0", beta0, sep = "\t"), file = filename, append = T, ncolumns = 3)
  write(paste("beta1", beta1, sep = "\t"), file = filename, append = T, ncolumns = 3)
  write(paste("target class", target.class, sep = "\t"), file = filename, append = T, ncolumns = 3)

  lset <- length(sample.molsig.sorted.subset.gs)

  filename <- paste(database.dir, identifier, ".msig.up.gmt", sep = "")
  gene.set <- sample.molsig.sorted.subset.gs[seq(1, lset / 2)]
  output.line <- paste(gene.set, sep = "\t", collapse = "\t")
  gene.set.name <- paste(identifier, ".up", sep = "")
  output.line <- paste(gene.set.name, gene.set.name, output.line, sep = "\t", collapse = "")
  write(noquote(output.line), file = filename, append = F, ncolumns = length(gene.set) + 2)

  filename <- paste(database.dir, identifier, ".msig.dn.gmt", sep = "")
  gene.set <- sample.molsig.sorted.subset.gs[seq(lset / 2 + 1, lset)]
  output.line <- paste(gene.set, sep = "\t", collapse = "\t")
  gene.set.name <- paste(identifier, ".dn", sep = "")
  output.line <- paste(gene.set.name, gene.set.name, output.line, sep = "\t", collapse = "")
  write(noquote(output.line), file = filename, append = F, ncolumns = length(gene.set) + 2)

  filename <- paste(database.dir, identifier, ".msig.gct", sep = "")
  z <- data.frame(sample.molsig.sorted.subset)
  names(z) <- sample.names
  row.names(z) <- sample.molsig.sorted.subset.gs
  write.gct(gct.data.frame = z, descs = sample.molsig.sorted.subset.gs, filename = filename)


  class.labels <- ifelse(class.labels == 1, 1, 2)

  print("class.labels:")
  print(class.labels)

  filename <- paste(database.dir, identifier, ".msig.cls", sep = "")
  write.cls(class.v = class.labels, phen = class.phen, filename = filename)
}

fix.sd <- function(s, m, s.percent = 0.2) { # function to threshold stdev as done in GeneCluster (so as to be
  # able to compare results)
  min.s <- s.percent * abs(m)
  if (min.s < s) {
    min.s <- s
  }
  if (min.s == 0) {
    min.s <- 0.2
  }
  return(min.s)
}

S2N_GC <- function(A, C) { # computes signal to noise ratio for one gene
  x <- split(A, C)
  mean.val <- sapply(x, mean)
  std.val <- sapply(x, sd)
  s2n <- (mean.val[1] - mean.val[2]) / (fix.sd(std.val[1], mean.val[1]) + fix.sd(std.val[2], mean.val[2]))
  return(s2n)
}

S2N <- function(A, C) { # computes robust signal to noise ratio for one gene
  x <- split(A, C)
  m1 <- mean(x[[1]])
  m2 <- mean(x[[2]])
  s1 <- ifelse(length(x[[1]]) > 1, sd(x[[1]]), 0)
  s2 <- ifelse(length(x[[2]]) > 1, sd(x[[2]]), 0)
  s1 <- ifelse(s1 < 0.1 * abs(m1), 0.1 * abs(m1), s1)
  s2 <- ifelse(s2 < 0.1 * abs(m2), 0.1 * abs(m2), s2)
  s2n <- (m1 - m2) / (s1 + s2 + 0.1)
  return(s2n)
}

RS2N <- function(A, C) { # computes robust signal to noise ratio for one gene
  x <- split(A, C)
  m1 <- median(x[[1]])
  m2 <- median(x[[2]])
  s1 <- mad(x[[1]])
  s2 <- mad(x[[2]])
  s1 <- ifelse(s1 < 0.1 * abs(m1), 0.1 * abs(m1), s1)
  s2 <- ifelse(s2 < 0.1 * abs(m2), 0.1 * abs(m2), s2)
  rs2n <- (m1 - m2) / (s1 + s2 + 0.1)
  return(rs2n)
}

MEAN.DIFF <- function(A, C) { # computes the difference in means for one gene
  x <- split(A, C)
  m1 <- mean(x[[1]])
  m2 <- mean(x[[2]])
  mean.diff <- m1 - m2
  return(mean.diff)
}

MEDIAN.DIFF <- function(A, C) { # computes the difference in medians for one gene
  x <- split(A, C)
  m1 <- median(x[[1]])
  m2 <- median(x[[2]])
  median.diff <- m1 - m2
  return(median.diff)
}

Gene.ranking <- function(dataset, class.labels, method = "S2N") {
  # This function ranks the genes according to the standard (S2N), robust (RS2N) signal to noise ratio
  # or difference of means (MEAN.DIFF) or ROC
  rows <- length(dataset[, 1])
  cols <- length(dataset[1, ])
  values.vector <- vector(length = rows, mode = "numeric")
  if (method == "S2N") {
    for (i in 1:rows) {
      #        values.vector[i] <- S2N_GC(dataset[i,], class.labels)
      values.vector[i] <- S2N(dataset[i, ], class.labels)
    }
  } else if (method == "RS2N") {
    for (i in 1:rows) {
      values.vector[i] <- RS2N(dataset[i, ], class.labels)
    }
  } else if (method == "MEAN.DIFF") {
    for (i in 1:rows) {
      values.vector[i] <- MEAN.DIFF(dataset[i, ], class.labels)
    }
  } else if (method == "MEDIAN.DIFF") {
    for (i in 1:rows) {
      values.vector[i] <- MEDIAN.DIFF(dataset[i, ], class.labels)
    }
  } else if (method == "ROC") {
    status <- ifelse(class.labels == class.labels[1], 1, 0)
    for (i in 1:rows) {
      m.score <- as.numeric(dataset[i, ])
      m.score.norm <- (m.score - min(m.score)) / (max(m.score) - min(m.score))
      values.vector[i] <- signif(roc.area(status, m.score.norm)$A, digits = 3)
    }
  } else {
    stop(c("Gene.ranking -- unknown method:", method))
  }
  return(values.vector)
}

fold.changes <- function(dataset, class.labels, method = "MEAN.DIFF", thres = 1) {
  # This function ranks the genes according to the standard (S2N), robust (RS2N) signal to noise ratio
  # or difference of means (MEAN.DIFF) or medians (MEDIAN.DIFF)
  rows <- length(dataset[, 1])
  cols <- length(dataset[1, ])
  fold.changes <- vector(length = rows, mode = "numeric")

  for (i in 1:rows) {
    y <- dataset[i, ]
    y[y < thres] <- thres
    x <- split(y, class.labels)
    if (method == "MEAN.DIFF") {
      m1 <- mean(x[[1]])
      m2 <- mean(x[[2]])
    } else if (method == "MEDIAN.DIFF") {
      m1 <- median(x[[1]])
      m2 <- median(x[[2]])
    } else {
      stop("ERROR in fold.changes(): unknown method")
    }
    if (abs(m2) > 1e-6) {
      fold.changes[i] <- m1 / m2
    } else {
      fold.changes[i] <- 1
    }
  }
  return(fold.changes)
}

consensusNMF <- function(input.ds, k.init, k.final, num.clusterings, maxniter, error.function, rseed = 123456789, directory = "", stopconv = 40, stopfreq = 10, non.interactive.run = F, doc.string = "", ...) {
  #
  #  GenePattern Methodology for:
  #
  #  Metagenes and Molecular Pattern Discovery using Matrix Factorization
  #  Jean-Philippe Brunet, Pablo Tamayo, Todd. R. Golub, and Jill P. Mesirov
  #
  #  Author:  Pablo Tamayo (tamayo@genome.wi.mit.edu)
  #
  #  Based on the original matlab version written by Jean-Philippe Brunet (brunet@broad.mit.edu) and
  #  with additional contributions from: Ted Liefeld (liefeld@broad.mit.edu)
  #  Date:  November 27, 2003
  #
  #  Last change March 3, 2005: modifications to make the output more readable.
  #
  #  Execute from an R console window with this command:
  #  source("<this file>", echo = TRUE)
  #  E.g. someoutput <- mynmf2(input.ds="c:\\nmf\\all_aml.res",k.init=2,k.final=5,num.clusterings=20,maxniter=500)
  #
  #  For details on the method see:
  #
  #  Proc. Natl. Acad. Sci. USA 2004 101: 4164-4169
  #  http://www.broad.mit.edu/cgi-bin/cancer/publications/pub_paper.cgi?mode=view&paper_id=89
  #
  #  Input parameters
  #
  #   input.ds
  #                       input gene expression dataset in GCT or RES format
  #   k.init
  #                       initial value of k
  #   k.final
  #                       final value of k
  #   num.clusterings
  #                       number of NMF clusterings to build consensus matrix
  #   maxniter
  #                       maximum number of NMF iterations
  #   error.function
  #                       NMF error function: "divergence" of "euclidean"
  #   rseed
  #                       random number generator seed
  #   directory
  #                       file directory where to store the result files
  #   stopconv
  #                       how many no change checks are needed to stop NMF iterations (convergence)
  #   stopfreq
  #                       frequency (NMF iterations) of "no change" checks
  #   non.interactive.run
  #                       flag controlling if the plots are produced interatively (Rgui and saved) or only saved in files
  #   doc.string
  #                       prefix to be added to the output files
  #
  #  Output files are (prefix with "doc.string")
  #
  #   params.txt
  #                       run parameters and time of execution
  #   membership.gct
  # 			membership results for samples at all values of K
  #   cophenetic.txt
  # 			cophenetic values for each K
  #   cophenetic.plot.jpeg
  # 			plot of cophenetic for each value of K
  #   consensus.k#.gct (for each value of K)
  # 			consensus matrix for k=#
  #   consensus.plot.k#.jpeg (for each value of K)
  # 			plot of consensus matrix for k=#
  #   graphs.k#.jpeg (for each value of K)

  # save input parameters

  filename <- paste(directory, doc.string, ".params.txt", sep = "", collapse = "")

  time.string <- as.character(as.POSIXlt(Sys.time(), "GMT"))
  write(paste("Run of NMF on ", time.string), file = filename)

  write(paste("input.ds =", input.ds, sep = " "), file = filename, append = T)
  write(paste("k.init = ", k.init, sep = " "), file = filename, append = T)
  write(paste("k.final =", k.final, sep = " "), file = filename, append = T)
  write(paste("num.clusterings =", num.clusterings, sep = " "), file = filename, append = T)
  write(paste("maxniter =", maxniter, sep = " "), file = filename, append = T)
  write(paste("error.function =", error.function, sep = " "), file = filename, append = T)
  write(paste("rseed =", rseed, sep = " "), file = filename, append = T)
  write(paste("directory =", directory, sep = " "), file = filename, append = T)
  write(paste("stopconv =", stopconv, sep = " "), file = filename, append = T)
  write(paste("stopfreq =", stopfreq, sep = " "), file = filename, append = T)
  write(paste("non.interctive.run =", non.interactive.run, sep = " "), file = filename, append = T)
  write(paste("doc.string =", doc.string, sep = " "), file = filename, append = T)


  k.init <- as.integer(k.init)
  k.final <- as.integer(k.final)
  num.clusterings <- as.integer(num.clusterings)
  n.iter <- as.integer(maxniter)
  if (!is.na(rseed)) {
    seed <- as.integer(rseed)
  }


  # library(mva)
  # library(MASS)
  # library(GenePattern)

  D <- CNMF.read.dataset(input.ds)
  A <- data.matrix(D)

  # Threshold negative values to small quantity

  eps <- .Machine$double.eps
  A[A < 0] <- eps



  cols <- length(A[1, ])
  rows <- length(A[, 1])

  col.names <- names(D)

  num.k <- k.final - k.init + 1

  rho <- vector(mode = "numeric", length = num.k)
  k.vector <- vector(mode = "numeric", length = num.k)

  k.index <- 1

  connect.matrix.ordered <- array(0, c(num.k, cols, cols))

  for (k in k.init:k.final) {
    if (non.interactive.run == F) {
      if (.Platform$OS.type == "windows") {
        filename <- paste(directory, doc.string, ".", "graphs.k", k, sep = "", collapse = "")
        #             x11(width = 9, height = 11)
        x11(width = 40, height = 22)
      } else if (.Platform$OS.type == "unix") {
        filename <- paste(directory, doc.string, ".", "graphs.k", k, ".pdf", sep = "", collapse = "")
        pdf(file = filename, width = 9, height = 11)
      }
    } else {
      if (.Platform$OS.type == "unix") {
        filename <- paste(directory, doc.string, ".", "graphs.k", k, ".pdf", sep = "", collapse = "")
        pdf(file = filename, width = 9, height = 11)
      } else if (.Platform$OS.type == "windows") {
        filename <- paste(directory, doc.string, ".", "graphs.k", k, ".pdf", sep = "", collapse = "")
        pdf(file = filename, width = 9, height = 11)
      }
    }

    #   nf <- layout(matrix(c(1,2,3,4,5,6,7,8), 4, 2, byrow=T), c(1, 1, 1, 1), c(1, 1), TRUE)
    nf <- layout(matrix(c(1, 2, 3, 4, 5, 6, 7, 8), 2, 4, byrow = T), c(1, 1), c(1, 1, 1, 1), TRUE)
    assign <- matrix(0, nrow = num.clusterings, ncol = cols)

    for (i in 1:num.clusterings) {
      print(paste("Computing clustering number=", i, " for k=", k, sep = ""))

      if (error.function == "divergence") {
        NMF.out <- NMF.div(V = A, k = k, maxniter = n.iter, seed = seed + i, stopconv = stopconv, stopfreq = stopfreq)
      } else if (error.function == "euclidean") {
        NMF.out <- NMF(V = A, k = k, maxniter = n.iter, seed = seed + i, stopconv = stopconv, stopfreq = stopfreq)
      } else {
        stop(paste("Un-supported error function=", error.function, sep = ""))
      }
      print(paste(NMF.out$t, " NMF iterations performed", sep = ""))

      for (j in 1:cols) { # Find membership
        class <- order(NMF.out$H[, j], decreasing = T)
        assign[i, j] <- class[1]
      }

      if (i == 1) { # Plot example for first clustering iteration
        H.saved <- NMF.out$H
        sub.string <- paste(doc.string, " k=", k, sep = "")
        plot(1:NMF.out$t, NMF.out$error.v[1:NMF.out$t], pch = 20, cex = 1.5, col = 1, xlab = "time", ylab = "NMF error", sub = sub.string, main = paste("Convergence plot k=", k, " example", sep = ""))


        if (rows < 1000) {
          W <- NMF.out$W
        } else {
          W <- NMF.out$W[sample(x = 1:rows, size = 1000), ]
        }
        sub.string <- paste(doc.string, " k=", k, sep = "")
        CNMF.matrix.abs.plot(W, sub = sub.string, log = F, main = "Example W matrix (orig. order)", ylab = "genes", xlab = "metasamples")
        CNMF.matrix.abs.plot(H.saved, sub = sub.string, log = F, main = "Example H matrix (orig. order)", ylab = "metagenes", xlab = "samples")
        CNMF.metagene.plot(H = H.saved, main = "Metagenes Example (orig. order)", sub = sub.string, xlab = "samples", ylab = "metagenes")
      }

      rm(NMF.out)
    } ## end  for (i in 1:num.clusterings)


    # compute consensus matrix
    connect.matrix <- matrix(0, nrow = cols, ncol = cols)

    for (i in 1:num.clusterings) {
      for (j in 1:cols) {
        for (p in 1:cols) {
          if (j != p) {
            if (assign[i, j] == assign[i, p]) {
              connect.matrix[j, p] <- connect.matrix[j, p] + 1
            }
          } else {
            connect.matrix[j, p] <- connect.matrix[j, p] + 1
          }
        }
      }
    }

    connect.matrix <- connect.matrix / num.clusterings

    dist.matrix <- 1 - connect.matrix
    dist.matrix <- as.dist(dist.matrix)
    HC <- hclust(dist.matrix, method = "average")

    dist.coph <- cophenetic(HC)
    k.vector[k.index] <- k
    rho[k.index] <- cor(dist.matrix, dist.coph)
    rho[k.index] <- signif(rho[k.index], digits = 4)

    #     connect.matrix.ordered <- matrix(0, nrow=cols, ncol = cols)

    for (i in 1:cols) {
      for (j in 1:cols) {
        connect.matrix.ordered[k.index, i, j] <- connect.matrix[HC$order[i], HC$order[j]]
      }
    }

    # compute consensus clustering membership

    membership <- cutree(HC, k = k)

    max.k <- max(membership)
    items.names.ordered <- col.names[HC$order]
    membership.ordered <- membership[HC$order]
    results <- data.frame(cbind(membership.ordered, items.names.ordered))

    if (k > k.init) {
      all.membership <- cbind(all.membership, membership)
    } else {
      all.membership <- cbind(membership)
    }

    sub.string <- paste(doc.string, " k=", k, sep = "")
    CNMF.matrix.abs.plot(connect.matrix.ordered[k.index, , ], sub = sub.string, log = F, main = "Ordered Consensus Matrix", ylab = "samples", xlab = "samples")
    plot(HC, xlab = "samples", cex = 0.75, labels = col.names, sub = sub.string, col = "blue", main = paste("Ordered Linkage Tree. Coph=", rho[k.index]))

    resultsGct <- data.frame(membership.ordered)
    row.names(resultsGct) <- items.names.ordered
    filename <- paste(directory, doc.string, ".", "consensus.k.", k, ".gct", sep = "", collapse = "")
    CNMF.write.gct(resultsGct, filename)

    H.sorted <- H.saved[, HC$order]
    sub.string <- paste(doc.string, " k=", k, sep = "")
    CNMF.matrix.abs.plot(H.sorted, sub = sub.string, log = F, main = "Example H matrix (ordered)", ylab = "metagenes", xlab = "samples")
    CNMF.metagene.plot(H = H.sorted, sub = sub.string, main = "Metagenes Example (ordered)", xlab = "samples", ylab = "metagenes")

    if (non.interactive.run == F) {
      if (.Platform$OS.type == "windows") {
        savePlot(filename = filename, type = "jpeg", device = dev.cur())
      } else if (.Platform$OS.type == "unix") {
        dev.off()
      }
    } else {
      dev.off()
    }

    if (non.interactive.run == F) {
      if (.Platform$OS.type == "windows") {
        filename <- paste(directory, doc.string, ".", "consensus.plot.k", k, sep = "", collapse = "")
        #             x11(width = 8.5, height = 11)
        x11(width = 20, height = 20)
      } else if (.Platform$OS.type == "unix") {
        filename <- paste(directory, doc.string, ".", "consensus.plot.k", k, ".pdf", sep = "", collapse = "")
        pdf(file = filename, width = 8.5, height = 11)
      }
    } else {
      if (.Platform$OS.type == "unix") {
        filename <- paste(directory, doc.string, ".", "consensus.plot.k", k, ".pdf", sep = "", collapse = "")
        pdf(file = filename, width = 8.5, height = 11)
      } else if (.Platform$OS.type == "windows") {
        filename <- paste(directory, doc.string, ".", "consensus.plot.k", k, ".pdf", sep = "", collapse = "")
        pdf(file = filename, width = 8.5, height = 11)
      }
    }

    nf <- layout(matrix(c(1), 1, 1, byrow = T), c(1, 1), c(1, 1), TRUE)

    conlabel <- paste("Consensus k =", k, sep = " ", collapse = "")

    sub.string <- paste("Consensus matrix k=", k, "; dataset= ", input.ds, sep = "")
    CNMF.ConsPlot(connect.matrix.ordered[k.index, , ], col.labels = membership.ordered, col.names = items.names.ordered, main = " ", sub = sub.string, xlab = " ", ylab = " ")

    if (non.interactive.run == F) {
      if (.Platform$OS.type == "windows") {
        savePlot(filename = filename, type = "jpeg", device = dev.cur())
      } else if (.Platform$OS.type == "unix") {
        dev.off()
      }
    } else {
      dev.off()
    }

    k.index <- k.index + 1
  } # end of loop over k


  # Save consensus matrices in one file

  if (non.interactive.run == F) {
    if (.Platform$OS.type == "windows") {
      filename <- paste(directory, doc.string, ".", "consensus.all.k.plot", sep = "")
      #             x11(width = 8.5, height = 11)
      x11(width = 30, height = 22)
    } else if (.Platform$OS.type == "unix") {
      filename <- paste(directory, doc.string, ".", "consensus.all.k.plot.pdf", sep = "")
      pdf(file = filename, width = 8.5, height = 11)
    }
  } else {
    if (.Platform$OS.type == "unix") {
      filename <- paste(directory, doc.string, ".", "consensus.all.k.plot.pdf", sep = "")
      pdf(file = filename, width = 8.5, height = 11)
    } else if (.Platform$OS.type == "windows") {
      filename <- paste(directory, doc.string, ".", "consensus.all.k.plot.pdf", sep = "")
      pdf(file = filename, width = 8.5, height = 11)
    }
  }

  nf <- layout(matrix(c(1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16), 4, 4, byrow = T), c(1, 1, 1, 1), c(1, 1, 1, 1), TRUE)

  for (k in 1:num.k) {
    CNMF.matrix.abs.plot(connect.matrix.ordered[k, , ],
      log = F, main = paste("k=", k.vector[k]),
      sub = paste("Cophenetic coef.=", rho[k]), ylab = "samples", xlab = "samples"
    )
  }

  y.range <- c(1 - 2 * (1 - min(rho)), 1)
  plot(k.vector, rho,
    main = "Cophenetic Coefficient", xlim = c(k.init, k.final), ylim = y.range,
    xlab = "k", ylab = "Cophenetic correlation", type = "n"
  )
  lines(k.vector, rho, type = "l", col = "black")
  points(k.vector, rho, pch = 22, type = "p", cex = 1.25, bg = "black", col = "black")

  if (non.interactive.run == F) {
    if (.Platform$OS.type == "windows") {
      savePlot(filename = filename, type = "jpeg", device = dev.cur())
    } else if (.Platform$OS.type == "unix") {
      dev.off()
    }
  } else {
    dev.off()
  }

  if (non.interactive.run == F) {
    if (.Platform$OS.type == "windows") {
      filename <- paste(directory, doc.string, ".", "cophenetic.plot", sep = "")
      #             x11(width = 8.5, height = 11)
      x11(width = 20, height = 20)
    } else if (.Platform$OS.type == "unix") {
      filename <- paste(directory, doc.string, ".", "cophenetic.plot.pdf", sep = "")
      pdf(file = filename, width = 8.5, height = 11)
    }
  } else {
    if (.Platform$OS.type == "unix") {
      filename <- paste(directory, doc.string, ".", "cophenetic.plot.pdf", sep = "")
      pdf(file = filename, width = 8.5, height = 11)
    } else if (.Platform$OS.type == "windows") {
      filename <- paste(directory, doc.string, ".", "cophenetic.plot.pdf", sep = "")
      pdf(file = filename, width = 8.5, height = 11)
    }
  }


  # Write the membership matrix

  resultsmembership <- data.frame(all.membership)
  row.names(resultsmembership) <- col.names

  print("Membership:")

  print(resultsmembership)

  filename <- paste(directory, doc.string, ".", "membership", ".txt", sep = "", collapse = "")

  write.table(resultsmembership, file = filename, append = FALSE, quote = FALSE, sep = "\t", eol = "\n", col.names = FALSE, row.names = FALSE)

  # CNMF.write.gct(resultsmembership , filename)

  y.range <- c(1 - 2 * (1 - min(rho)), 1)
  plot(k.vector, rho, main = "Cophenetic Coefficient", xlim = c(k.init, k.final), ylim = y.range, xlab = "k", ylab = "Cophenetic correlation", type = "n")
  lines(k.vector, rho, type = "l", col = "black")
  points(k.vector, rho, pch = 22, type = "p", cex = 1.25, bg = "black", col = "black")

  if (non.interactive.run == F) {
    if (.Platform$OS.type == "windows") {
      savePlot(filename = filename, type = "jpeg", device = dev.cur())
    } else if (.Platform$OS.type == "unix") {
      dev.off()
    }
  } else {
    dev.off()
  }

  xx <- cbind(k.vector, rho)
  write(xx, file = paste(directory, doc.string, ".", "cophenetic.txt", sep = ""))
}


CNMF.read.dataset <- function(file) {
  result <- regexpr(paste(".gct", "$", sep = ""), tolower(file))
  if (result[[1]] != -1) {
    return(CNMF.read.gct(file))
  }
  result <- regexpr(paste(".res", "$", sep = ""), tolower(file))
  if (result[[1]] != -1) {
    return(CNMF.read.res(file))
  }
  stop("Input is not a res or gct file.")
}

CNMF.matrix.abs.plot <- function(V, axes = F, log = F, norm = T, transpose = T, matrix.order = T, max.v = 1, min.v = 0, main = " ", sub = " ", xlab = " ", ylab = "  ") {
  rows <- length(V[, 1])
  cols <- length(V[1, ])
  if (log == T) {
    V <- log(V)
  }
  B <- matrix(0, nrow = rows, ncol = cols)
  for (i in 1:rows) {
    for (j in 1:cols) {
      if (matrix.order == T) {
        k <- rows - i + 1
      } else {
        k <- i
      }
      if (norm == T) {
        if ((max.v == 1) && (min.v == 0)) {
          max.val <- max(V)
          min.val <- min(V)
        } else {
          max.val <- max.v
          min.val <- min.v
        }
      }
      B[k, j] <- max.val - V[i, j] + min.val
    }
  }
  if (transpose == T) {
    B <- t(B)
  }
  if (norm == T) {
    image(z = B, zlim = c(min.val, max.val), axes = axes, col = rainbow(100, s = 1.0, v = 0.75, start = 0.0, end = 0.75, gamma = 1.5), main = main, sub = sub, xlab = xlab, ylab = ylab)
  } else {
    image(z = B, axes = axes, col = rainbow(100, s = 1, v = 0.6, start = 0.1, end = 0.9, gamma = 1), main = main, sub = sub, xlab = xlab, ylab = ylab)
  }
  return(list(B, max.val, min.val))
}

CNMF.metagene.plot <- function(H, main = " ", sub = " ", xlab = "samples ", ylab = "amplitude") {
  k <- length(H[, 1])
  S <- length(H[1, ])
  index <- 1:S
  maxval <- max(H)
  minval <- min(H)
  plot(index, H[1, ], xlim = c(1, S), ylim = c(minval, maxval), main = main, sub = sub, ylab = ylab, xlab = xlab, type = "n")
  for (i in 1:k) {
    lines(index, H[i, ], type = "l", col = i, lwd = 2)
  }
}


CNMF.ConsPlot <- function(V, col.labels, col.names, main = " ", sub = " ", xlab = " ", ylab = " ") {
  # Plots a heatmap plot of a consensus matrix

  cols <- length(V[1, ])
  B <- matrix(0, nrow = cols, ncol = cols)
  max.val <- max(V)
  min.val <- min(V)
  for (i in 1:cols) {
    for (j in 1:cols) {
      k <- cols - i + 1
      B[k, j] <- max.val - V[i, j] + min.val
    }
  }

  col.names2 <- rev(col.names)
  col.labels2 <- rev(col.labels)
  D <- matrix(0, nrow = (cols + 1), ncol = (cols + 1))

  col.tag <- vector(length = cols, mode = "numeric")
  current.tag <- 0
  col.tag[1] <- current.tag
  for (i in 2:cols) {
    if (col.labels[i] != col.labels[i - 1]) {
      current.tag <- 1 - current.tag
    }
    col.tag[i] <- current.tag
  }
  col.tag2 <- rev(col.tag)
  D[(cols + 1), 2:(cols + 1)] <- ifelse(col.tag %% 2 == 0, 1.02, 1.01)
  D[1:cols, 1] <- ifelse(col.tag2 %% 2 == 0, 1.02, 1.01)
  D[(cols + 1), 1] <- 1.03
  D[1:cols, 2:(cols + 1)] <- B[1:cols, 1:cols]

  col.map <- c(rainbow(100, s = 1.0, v = 0.75, start = 0.0, end = 0.75, gamma = 1.5), "#BBBBBB", "#333333", "#FFFFFF")
  image(1:(cols + 1), 1:(cols + 1), t(D), col = col.map, axes = FALSE, main = main, sub = sub, xlab = xlab, ylab = ylab)
  for (i in 1:cols) {
    col.names[i] <- paste("      ", substr(col.names[i], 1, 12), sep = "")
    col.names2[i] <- paste(substr(col.names2[i], 1, 12), "     ", sep = "")
  }

  axis(2, at = 1:cols, labels = col.names2, adj = 0.5, tick = FALSE, las = 1, cex.axis = 0.50, font.axis = 1, line = -1)
  axis(2, at = 1:cols, labels = col.labels2, adj = 0.5, tick = FALSE, las = 1, cex.axis = 0.65, font.axis = 1, line = -1)

  axis(3, at = 2:(cols + 1), labels = col.names, adj = 1, tick = FALSE, las = 3, cex.axis = 0.50, font.axis = 1, line = -1)
  axis(3, at = 2:(cols + 1), labels = as.character(col.labels), adj = 1, tick = FALSE, las = 1, cex.axis = 0.65, font.axis = 1, line = -1)

  return()
}

CNMF.read.res <- function(filename = "NULL") {
  #
  # Reads a gene expression dataset in RES format and converts it into an R data frame
  #
  header.cont <- readLines(filename, n = 1)
  temp <- unlist(strsplit(header.cont, "\t"))
  colst <- length(temp)
  header.labels <- temp[seq(3, colst, 2)]
  ds <- read.delim(filename, header = F, row.names = 2, sep = "\t", skip = 3, blank.lines.skip = T, comment.char = "", as.is = T)
  colst <- length(ds[1, ])
  cols <- (colst - 1) / 2
  rows <- length(ds[, 1])
  A <- matrix(nrow = rows - 1, ncol = cols)
  A <- ds[1:rows, seq(2, colst, 2)]
  table1 <- data.frame(A)
  names(table1) <- header.labels
  return(table1)
}

CNMF.read.gct <- function(filename = "NULL") {
  #
  # Reads a gene expression dataset in GCT format and converts it into an R data frame
  #
  ds <- read.delim(filename, header = T, sep = "\t", skip = 2, row.names = 1, blank.lines.skip = T, comment.char = "", as.is = T)
  ds <- ds[-1]
  return(ds)
}

CNMF.write.gct <- function(gct, filename) {
  f <- file(filename, "w")
  cat("#1.2", "\n", file = f, append = TRUE, sep = "")
  cat(dim(gct)[1], "     ", dim(gct)[2], "\n", file = f, append = TRUE, sep = "")
  cat("Name", "       ", file = f, append = TRUE, sep = "")
  cat("Description", file = f, append = TRUE, sep = "")
  names <- names(gct)
  cat("       ", names[1], file = f, append = TRUE, sep = "")
  for (j in 2:length(names)) {
    cat("   ", names[j], file = f, append = TRUE, sep = "")
  }
  cat("\n", file = f, append = TRUE, sep = "")
  oldWarn <- options(warn = -1)
  m <- matrix(nrow = dim(gct)[1], ncol = dim(gct)[2] + 2)
  m[, 1] <- row.names(gct)
  m[, 2] <- row.names(gct)
  index <- 3
  for (i in 1:dim(gct)[2]) {
    m[, index] <- gct[, i]
    index <- index + 1
  }
  write.table(m, file = f, append = TRUE, quote = FALSE, sep = "      ", eol = "\n", col.names = FALSE, row.names = FALSE)
  close(f)
  options(warn = 0)
  return(gct)
}

GSEA.EnrichmentScore4 <- function(gene.list, gene.set, statistic = "Kolmogorov-Smirnov", alpha = 1, correl.vector = NULL) {
  #
  # Computes the weighted GSEA score of gene.set in gene.list.
  #
  # This version supports multiple statistics:
  #  "Kolmogorov-Smirnov" "Cramer-von-Mises" "Anderson-Darling" "Zhang_A" "Zhang_C" "Zhang_K"
  #
  # The weighted score type is the exponent of the correlation
  # weight: 0 (unweighted = Kolmogorov-Smirnov), 1 (weighted), and 2 (over-weighted). When the score type is 1 or 2 it is
  # necessary to input the correlation vector with the values in the same order as in the gene list.
  #
  # Inputs:
  #   gene.list: The ordered gene list (e.g. integers indicating the original position in the input dataset)
  #   gene.set: A gene set (e.g. integers indicating the location of those genes in the input dataset)
  #   weighted.score.type: Type of score: weight: 0 (unweighted = Kolmogorov-Smirnov), 1 (weighted), and 2 (over-weighted)
  #  correl.vector: A vector with the coorelations (e.g. signal to noise scores) corresponding to the genes in the gene list
  #
  # Outputs:
  #   ES: Enrichment score (real number between -1 and +1)
  #   arg.ES: Location in gene.list where the peak running enrichment occurs (peak of the "mountain")
  #   RES: Numerical vector containing the running enrichment score for all locations in the gene list
  #   tag.indicator: Binary vector indicating the location of the gene sets (1's) in the gene list
  #
  # The Broad Institute
  # SOFTWARE COPYRIGHT NOTICE AGREEMENT
  # This software and its documentation are copyright 2003 by the
  # Broad Institute/Massachusetts Institute of Technology.
  # All rights are reserved.
  #
  # This software is supplied without any warranty or guaranteed support
  # whatsoever. Neither the Broad Institute nor MIT can be responsible for
  # its use, misuse, or functionality.

  tag.indicator <- sign(match(gene.list, gene.set, nomatch = 0)) # notice that the sign is 0 (no tag) or 1 (tag)
  no.tag.indicator <- 1 - tag.indicator
  N <- length(gene.list)
  Nh <- length(gene.set)
  Nm <- N - Nh
  orig.correl.vector <- correl.vector
  if (alpha == 0) correl.vector <- rep(1, N) # unweighted case
  correl.vector <- abs(correl.vector)^alpha
  sum.correl <- sum(correl.vector[tag.indicator == 1])
  P0 <- no.tag.indicator / Nm
  F0 <- cumsum(P0)
  Pn <- tag.indicator * correl.vector / sum.correl
  Fn <- cumsum(Pn)
  if (statistic == "Kolmogorov-Smirnov") {
    RES <- Fn - F0
    max.ES <- max(RES)
    min.ES <- min(RES)
    if (max.ES > -min.ES) {
      ES <- signif(max.ES, digits = 5)
      arg.ES <- which.max(RES)
    } else {
      ES <- signif(min.ES, digits = 5)
      arg.ES <- which.min(RES)
    }
    return(list(ES = ES, arg.ES = arg.ES, RES = RES, indicator = tag.indicator))
  } else if (statistic == "Cramer-von-Mises") {
    # Based on Choulakian et al Canadian J. of Statistics 22, 125 (1994)
    # Modified to separate positive and negative enrichment parts
    RES <- Fn - F0
    X <- RES^2 * P0
    X_p <- X[RES >= 0]
    X_n <- X[RES < 0]
    ES_p <- sqrt(sum(X_p) / N)
    ES_n <- sqrt(sum(X_n) / N)
    if (ES_p > ES_n) {
      ES <- signif(ES_p, digits = 5)
      arg.ES <- which.min(abs(X - max(X_p)))
    } else {
      ES <- -signif(ES_n, digits = 5)
      arg.ES <- which.min(abs(X - max(X_n)))
    }
    return(list(ES = ES, RES = RES, arg.ES = arg.ES, indicator = tag.indicator))
  } else if (statistic == "Anderson-Darling") {
    # Based on Choulakian et al Canadian J. of Statistics 22, 125 (1994)
    # Modified to separate positive and negative enrichment parts
    RES <- Fn - F0
    F0_factor <- ifelse(F0 < 1 / Nm | F0 > (Nm - 1) / Nm, rep(1, N), F0 * (1 - F0))
    X <- RES^2 * P0 / F0_factor
    X_p <- X[RES >= 0]
    X_n <- X[RES < 0]
    ES_p <- sqrt(sum(X_p) / N)
    ES_n <- sqrt(sum(X_n) / N)
    if (ES_p > ES_n) {
      ES <- signif(ES_p, digits = 5)
      arg.ES <- which.min(abs(X - max(X_p)))
    } else {
      ES <- -signif(ES_n, digits = 5)
      arg.ES <- which.min(abs(X - max(X_n)))
    }
    return(list(ES = ES, RES = RES, arg.ES = arg.ES, indicator = tag.indicator))
  } else if (statistic == "Zhang_A") {
    # Based on Zhang, J. R. Statist. Soc. B 64, Part2, 281 (2002)
    # Modified to separate positive and negative enrichment parts
    RES <- Fn - F0
    Fact1 <- ifelse(F0 < 1 / Nm | Fn < 1 / sum.correl, 0, Fn * log(Fn / F0))
    Fact2 <- ifelse(F0 > (Nm - 1) / Nm | Fn > (sum.correl - 1) / sum.correl, 0, (1 - Fn) * log((1 - Fn) / (1 - F0)))
    Fn_factor <- ifelse(Fn < 1 / sum.correl | Fn > (sum.correl - 1) / sum.correl, rep(1, N), Fn * (1 - Fn))
    G <- (Fact1 + Fact2) * Pn / Fn_factor
    G_p <- G[RES >= 0]
    G_n <- G[RES < 0]
    ES_p <- sum(G_p) / N
    ES_n <- sum(G_n) / N
    if (ES_p > ES_n) {
      ES <- signif(ES_p, digits = 5)
      arg.ES <- which.min(abs(G - max(G_p)))
    } else {
      ES <- -signif(ES_n, digits = 5)
      arg.ES <- which.min(abs(G - max(G_n)))
    }
    return(list(ES = ES, RES = RES, arg.ES = arg.ES, indicator = tag.indicator))
  } else if (statistic == "Zhang_C") {
    # Based on Zhang, J. R. Statist. Soc. B 64, Part2, 281 (2002)
    # Modified to separate positive and negative enrichment parts
    RES <- Fn - F0
    Fact1 <- ifelse(F0 < 1 / Nm | Fn < 1 / sum.correl, 0, Fn * log(Fn / F0))
    Fact2 <- ifelse(F0 > (Nm - 1) / Nm | Fn > (sum.correl - 1) / sum.correl, 0, (1 - Fn) * log((1 - Fn) / (1 - F0)))
    F0_factor <- ifelse(F0 < 1 / Nm | F0 > (Nm - 1) / Nm, rep(1, N), F0 * (1 - F0))
    G <- (Fact1 + Fact2) * P0 / F0_factor
    G_p <- G[RES >= 0]
    G_n <- G[RES < 0]
    ES_p <- sum(G_p) / N
    ES_n <- sum(G_n) / N
    if (ES_p > ES_n) {
      ES <- signif(ES_p, digits = 5)
      arg.ES <- which.min(abs(G - max(G_p)))
    } else {
      ES <- -signif(ES_n, digits = 5)
      arg.ES <- which.min(abs(G - max(G_n)))
    }
    return(list(ES = ES, RES = RES, arg.ES = arg.ES, indicator = tag.indicator))
  } else if (statistic == "Zhang_K") {
    # Based on Zhang, J. R. Statist. Soc. B 64, Part2, 281 (2002)
    # Modified to separate positive and negative enrichment parts
    RES <- Fn - F0
    Fact1 <- ifelse(F0 < 1 / Nm | Fn < 1 / sum.correl, 0, Fn * log(Fn / F0))
    Fact2 <- ifelse(F0 > (Nm - 1) / Nm | Fn > (sum.correl - 1) / sum.correl, 0, (1 - Fn) * log((1 - Fn) / (1 - F0)))
    G <- Fact1 + Fact2
    G_p <- G[RES >= 0]
    G_n <- G[RES < 0]
    ES_p <- max(G_p)
    ES_n <- max(G_n)
    if (ES_p > ES_n) {
      ES <- signif(ES_p, digits = 5)
      arg.ES <- which.min(abs(G - ES_p))
    } else {
      ES <- -signif(ES_n, digits = 5)
      arg.ES <- which.min(abs(G - ES_n))
    }
    return(list(ES = ES, RES = RES, arg.ES = arg.ES, indicator = tag.indicator))
  }
}

MSIG.Signature.Plot.2 <- function(
    V,
    row.names = "NA",
    row.names2 = "NA",
    col.labels = "NA",
    col.classes = "NA",
    phen.cmap = "NA", # this color map contain the colors for the col.classes plus one extra color for the missing values
    col.names = "NA",
    main = " ",
    sub = " ",
    xlab = " ",
    ylab = " ",
    char.rescale = 1.0,
    max.v = "NA",
    seed = 1729) { # This function plots small sets of genes (signatures). It supports missing values (NA)

  n.rows <- length(V[, 1])
  n.cols <- length(V[1, ])
  V[V > 4] <- 4
  V[V < -4] <- -4
  mycol <- vector(length = 512, mode = "numeric")

  for (k in 1:256) {
    mycol[k] <- rgb(255, k - 1, k - 1, maxColorValue = 255)
  }
  for (k in 257:512) {
    mycol[k] <- rgb(511 - (k - 1), 511 - (k - 1), 255, maxColorValue = 255)
  }
  mycol <- rev(mycol)
  ncolors <- length(mycol)
  if (is.na(max(V))) {
    mycol <- c(mycol, phen.cmap[1:(length(col.classes) + 1)])
  } else {
    mycol <- c(mycol, phen.cmap[1:length(col.classes)])
  }

  V2 <- ceiling(ncolors * (V - -4) / (8 + 0.0001))

  #   print(V)
  #   print(V2)

  heatm <- matrix(0, nrow = n.rows + 1, ncol = n.cols)
  heatm[1:n.rows, ] <- V2[seq(n.rows, 1, -1), ]

  for (i in 1:length(heatm[, 1])) { # replace missing values with last color in color map
    for (j in 1:length(heatm[1, ])) {
      heatm[i, j] <- ifelse(is.na(heatm[i, j]), ncolors + max(col.labels) + 1, heatm[i, j])
    }
  }

  heatm[n.rows + 1, ] <- ncolors + col.labels
  height <- ifelse(n.rows >= 25, 25, n.rows * 0.8 + 5)

  x11(width = 31, height = 19)
  nf <- layout(matrix(c(1, 2), 1, 2, byrow = T), widths = c(7, 1), respect = FALSE)

  par(mar = c(5, 7, 4, 7))

  if (is.na(max(V))) {
    image(1:n.cols, 1:(n.rows + 1), t(heatm),
      zlim = c(0, ncolors + max(col.labels) + 1), col = mycol,
      axes = FALSE, main = main, sub = sub, xlab = xlab, ylab = ylab
    )
  } else {
    image(1:n.cols, 1:(n.rows + 1), t(heatm),
      zlim = c(0, ncolors + max(col.labels)), col = mycol,
      axes = FALSE, main = main, sub = sub, xlab = xlab, ylab = ylab
    )
  }

  #    image(1:n.cols, 1:(n.rows + 1), t(heatm), col=mycol, axes=FALSE, main=main,sub = sub, xlab= xlab, ylab=ylab)

  if (row.names[1] != "NA") {
    numC <- nchar(row.names)
    size.row.char <- char.rescale * 30 / (n.rows + 15)
    size.col.char <- char.rescale * 20 / (n.cols + 15)
    for (i in 1:n.rows) {
      row.names[i] <- substr(row.names[i], 1, 30)
    }
    row.names <- c(row.names[seq(n.rows, 1, -1)], "Class")
    axis(2,
      at = 1:(n.rows + 1), labels = row.names, adj = 0.5, tick = FALSE, las = 1, cex.axis = size.row.char,
      font.axis = 2, line = -1
    )
  }
  if (row.names2[1] != "NA") {
    numC <- nchar(row.names2)
    size.row.char <- char.rescale * 35 / (n.rows + 15)
    size.col.char <- char.rescale * 20 / (n.cols + 15)
    for (i in 1:n.rows) {
      row.names2[i] <- substr(row.names2[i], 1, 30)
    }
    row.names2 <- c(row.names2[seq(n.rows, 1, -1)], "Class")
    axis(4,
      at = 1:(n.rows + 1), labels = row.names2, adj = 0.5, tick = FALSE, las = 1, cex.axis = size.row.char,
      font.axis = 2, line = -1
    )
  }
  if (col.names[1] != "NA") {
    axis(1, at = 1:n.cols, labels = col.names, tick = FALSE, las = 3, cex.axis = size.col.char, font.axis = 2, line = -1)
  }

  if (col.classes[1] != "NA") {
    C <- split(col.labels, col.labels)
    class1.size <- length(C[[1]])
    class2.size <- length(C[[2]])
    axis(3,
      at = c(class1.size / 2, class1.size + class2.size / 2), labels = col.classes,
      tick = FALSE, las = 1, cex.axis = 1.25, font.axis = 2, line = -1
    )
  }


  par(mar = c(10, 2, 10, 2))
  num.v <- 20
  range.v <- range(V2, na.rm = T)
  incr <- (range.v[1] - range.v[2]) / (num.v - 1)
  heatm.v <- matrix(rev(seq(range.v[2], range.v[1], incr)), nrow = num.v, ncol = 1)
  if (is.na(max(V))) {
    image(1:1, 1:num.v, t(heatm.v),
      zlim = c(0, ncolors + max(col.labels) + 1), col = mycol, axes = FALSE,
      sub = "Color \n Legend ", main = " ", xlab = xlab, ylab = ylab
    )
  } else {
    image(1:1, 1:num.v, t(heatm.v),
      zlim = c(0, ncolors + max(col.labels)), col = mycol, axes = FALSE,
      sub = "Color \n Legend ", main = " ", xlab = xlab, ylab = ylab
    )
  }
  range.v <- range(V, na.rm = T)
  incr <- (range.v[1] - range.v[2]) / (num.v - 1)
  heatm.v2 <- matrix(signif(rev(seq(range.v[2], range.v[1], incr)), digits = 3), nrow = num.v, ncol = 1)
  #          print(c("heatm.v2=", heatm.v2))
  axis(2,
    at = 1:num.v, labels = heatm.v2, adj = 0.5, tick = FALSE, las = 1, cex.axis = char.rescale * 0.6,
    font.axis = 1.25, line = -0.8
  )
}

MSIG.ReadClsFile2 <- function(file = "NULL") {
  #
  # Reads a class vector CLS file and defines phenotype and class labels vectors (numeric and character) for the samples in a gene expression file (RES or GCT format)
  #
  # The Broad Institute
  # SOFTWARE COPYRIGHT NOTICE AGREEMENT
  # This software and its documentation are copyright 2003 by the
  # Broad Institute/Massachusetts Institute of Technology.
  # All rights are reserved.
  #
  # This software is supplied without any warranty or guaranteed support
  # whatsoever. Neither the Broad Institute nor MIT can be responsible for
  # its use, misuse, or functionality.

  cls.cont <- readLines(file)

  class.list <- unlist(strsplit(cls.cont[[3]], " "))
  phen <- unique(class.list)
  class.v <- match(class.list, phen)
  return(list(phen = phen, class.v = class.v, class.list = class.list))
}


MSIG.Sort.Dataset <- function(
    input.ds,
    input.cls,
    phen.order = NULL, # (optional) order of phenotypes
    output.ds,
    output.cls) {
  # Sorts a CLS and GCT file according to phenotype labels and optionally reorders phenotypes

  # start of methodology

  print(c("Running MSIG.Sort.Dataset... on GCT file:", input.ds))
  print(c("Running MSIG.Sort.Dataset... on CLS file:", input.cls))

  # Read input datasets

  dataset <- MSIG.Gct2Frame(filename = input.ds)
  m <- data.matrix(dataset$ds)
  gs.names <- dataset$row.names
  gs.descs <- dataset$descs
  sample.names <- dataset$names

  # Read CLS file

  CLS <- MSIG.ReadClsFile2(file = input.cls)
  class.labels <- CLS$class.v
  class.phen <- CLS$phen
  class.list <- CLS$class.list

  # Sort according to phenotype

  class.labels.order <- order(class.labels, decreasing = F)
  class.labels <- class.labels[class.labels.order]
  m <- m[, class.labels.order]
  sample.names <- sample.names[class.labels.order]

  # Reorder phenotypes

  if (!is.null(phen.order)) {
    class.phen <- class.phen[phen.order]
    new.labels <- vector(length(class.labels), mode = "numeric")
    for (i in 1:length(class.labels)) {
      new.labels[i] <- phen.order[class.labels[i]]
    }
    col.index <- order(new.labels, decreasing = F)
    class.labels <- new.labels[col.index]
    sample.names <- sample.names[col.index]
    rows <- length(m[, 1])
    for (j in 1:rows) {
      m[j, ] <- m[j, col.index]
    }
  }

  # Save datasets

  V <- data.frame(m)
  names(V) <- sample.names
  row.names(V) <- gs.names
  write.gct(gct.data.frame = V, descs = gs.descs, filename = output.ds)
  write.cls(class.v = class.labels, phen = class.phen, filename = output.cls)
}

GSEA.EnrichmentScore5 <- function(
    #
    # Computes the weighted GSEA score of gene.set in gene.list.
    #
    # This version supports multiple statistics:
    #  "Kolmogorov-Smirnov" "Cramer-von-Mises" "Anderson-Darling" "Zhang_A" "Zhang_C" "Zhang_K" and "area.under.RES"
    #
    # The weighted score type is the exponent of the correlation
    # weight: 0 (unweighted = Kolmogorov-Smirnov), 1 (weighted), and 2 (over-weighted). When the score type is 1 or 2 it is
    # necessary to input the correlation vector with the values in the same order as in the gene list.
    #
    # Inputs:
    #   gene.list: The ordered gene list (e.g. integers indicating the original position in the input dataset)
    #   gene.set: A gene set (e.g. integers indicating the location of those genes in the input dataset)
    #   weighted.score.type: "Kolmogorov-Smirnov" "Cramer-von-Mises" "Anderson-Darling" "Zhang_A" "Zhang_C" "Zhang_K", "area.under.RES", or "Wilcoxon"
    #  correl.vector: A vector with the coorelations (e.g. signal to noise scores) corresponding to the genes in the gene list
    #
    # Outputs:
    #   ES: Enrichment score (e.g. real number between -1 and +1 for KS)
    #   arg.ES: Location in gene.list where the peak running enrichment occurs (peak of the "mountain")
    #   RES: Numerical vector containing the running enrichment score for all locations in the gene list
    #   tag.indicator: Binary vector indicating the location of the gene sets (1's) in the gene list
    #
    # The Broad Institute
    # SOFTWARE COPYRIGHT NOTICE AGREEMENT
    # This software and its documentation are copyright 2008 by the
    # Broad Institute/Massachusetts Institute of Technology.
    # All rights are reserved.
    #
    # This software is supplied without any warranty or guaranteed support
    # whatsoever. Neither the Broad Institute nor MIT can be responsible for
    # its use, misuse, or functionality.
    gene.list, # The ordered gene list (e.g. integers indicating the original position in the input dataset)
    gene.set, # A gene set (e.g. integers indicating the location of those genes in the input dataset)
    statistic = "Kolmogorov-Smirnov", # "Kolmogorov-Smirnov", "Cramer-von-Mises", "Anderson-Darling", "Zhang_A", "Zhang_C",
    #  "Zhang_K", "area.under.RES", or "Wilcoxon"
    alpha = 1, # The weight exponent
    correl.vector = NULL # A correlation vector of genes with phenotype or other appropriate weight
    ) {
  tag.indicator <- sign(match(gene.list, gene.set, nomatch = 0)) # notice that the sign is 0 (no tag) or 1 (tag)
  no.tag.indicator <- 1 - tag.indicator
  N <- length(gene.list)
  Nh <- length(gene.set)
  Nm <- N - Nh
  orig.correl.vector <- correl.vector
  if (alpha == 0) correl.vector <- rep(1, N) # unweighted case
  correl.vector <- abs(correl.vector)^alpha
  sum.correl <- sum(correl.vector[tag.indicator == 1])
  P0 <- no.tag.indicator / Nm
  F0 <- cumsum(P0)
  Pn <- tag.indicator * correl.vector / sum.correl
  Fn <- cumsum(Pn)
  if (statistic == "Kolmogorov-Smirnov") {
    RES <- Fn - F0
    max.ES <- max(RES)
    min.ES <- min(RES)
    if (max.ES > -min.ES) {
      ES <- signif(max.ES, digits = 5)
      arg.ES <- which.max(RES)
    } else {
      ES <- signif(min.ES, digits = 5)
      arg.ES <- which.min(RES)
    }
    return(list(ES = ES, arg.ES = arg.ES, RES = RES, indicator = tag.indicator))
  } else if (statistic == "Cramer-von-Mises") {
    # Based on Choulakian et al Canadian J. of Statistics 22, 125 (1994)
    # Modified to separate positive and negative enrichment parts
    RES <- Fn - F0
    X <- RES^2 * P0
    X_p <- X[RES >= 0]
    X_n <- X[RES < 0]
    ES_p <- sqrt(sum(X_p) / N)
    ES_n <- sqrt(sum(X_n) / N)
    if (ES_p > ES_n) {
      ES <- signif(ES_p, digits = 5)
      arg.ES <- which.min(abs(X - max(X_p)))
    } else {
      ES <- -signif(ES_n, digits = 5)
      arg.ES <- which.min(abs(X - max(X_n)))
    }
    return(list(ES = ES, RES = RES, arg.ES = arg.ES, indicator = tag.indicator))
  } else if (statistic == "Anderson-Darling") {
    # Based on Choulakian et al Canadian J. of Statistics 22, 125 (1994)
    # Modified to separate positive and negative enrichment parts
    RES <- Fn - F0
    F0_factor <- ifelse(F0 < 1 / Nm | F0 > (Nm - 1) / Nm, rep(1, N), F0 * (1 - F0))
    X <- RES^2 * P0 / F0_factor
    X_p <- X[RES >= 0]
    X_n <- X[RES < 0]
    ES_p <- sqrt(sum(X_p) / N)
    ES_n <- sqrt(sum(X_n) / N)
    if (ES_p > ES_n) {
      ES <- signif(ES_p, digits = 5)
      arg.ES <- which.min(abs(X - max(X_p)))
    } else {
      ES <- -signif(ES_n, digits = 5)
      arg.ES <- which.min(abs(X - max(X_n)))
    }
    return(list(ES = ES, RES = RES, arg.ES = arg.ES, indicator = tag.indicator))
  } else if (statistic == "Zhang_A") {
    # Based on Zhang, J. R. Statist. Soc. B 64, Part2, 281 (2002)
    # Modified to separate positive and negative enrichment parts
    RES <- Fn - F0
    Fact1 <- ifelse(F0 < 1 / Nm | Fn < 1 / sum.correl, 0, Fn * log(Fn / F0))
    Fact2 <- ifelse(F0 > (Nm - 1) / Nm | Fn > (sum.correl - 1) / sum.correl, 0, (1 - Fn) * log((1 - Fn) / (1 - F0)))
    Fn_factor <- ifelse(Fn < 1 / sum.correl | Fn > (sum.correl - 1) / sum.correl, rep(1, N), Fn * (1 - Fn))
    G <- (Fact1 + Fact2) * Pn / Fn_factor
    G_p <- G[RES >= 0]
    G_n <- G[RES < 0]
    ES_p <- sum(G_p) / N
    ES_n <- sum(G_n) / N
    if (ES_p > ES_n) {
      ES <- signif(ES_p, digits = 5)
      arg.ES <- which.min(abs(G - max(G_p)))
    } else {
      ES <- -signif(ES_n, digits = 5)
      arg.ES <- which.min(abs(G - max(G_n)))
    }
    return(list(ES = ES, RES = RES, arg.ES = arg.ES, indicator = tag.indicator))
  } else if (statistic == "Zhang_C") {
    # Based on Zhang, J. R. Statist. Soc. B 64, Part2, 281 (2002)
    # Modified to separate positive and negative enrichment parts
    RES <- Fn - F0
    Fact1 <- ifelse(F0 < 1 / Nm | Fn < 1 / sum.correl, 0, Fn * log(Fn / F0))
    Fact2 <- ifelse(F0 > (Nm - 1) / Nm | Fn > (sum.correl - 1) / sum.correl, 0, (1 - Fn) * log((1 - Fn) / (1 - F0)))
    F0_factor <- ifelse(F0 < 1 / Nm | F0 > (Nm - 1) / Nm, rep(1, N), F0 * (1 - F0))
    G <- (Fact1 + Fact2) * P0 / F0_factor
    G_p <- G[RES >= 0]
    G_n <- G[RES < 0]
    ES_p <- sum(G_p) / N
    ES_n <- sum(G_n) / N
    if (ES_p > ES_n) {
      ES <- signif(ES_p, digits = 5)
      arg.ES <- which.min(abs(G - max(G_p)))
    } else {
      ES <- -signif(ES_n, digits = 5)
      arg.ES <- which.min(abs(G - max(G_n)))
    }
    return(list(ES = ES, RES = RES, arg.ES = arg.ES, indicator = tag.indicator))
  } else if (statistic == "Zhang_K") {
    # Based on Zhang, J. R. Statist. Soc. B 64, Part2, 281 (2002)
    # Modified to separate positive and negative enrichment parts
    RES <- Fn - F0
    Fact1 <- ifelse(F0 < 1 / Nm | Fn < 1 / sum.correl, 0, Fn * log(Fn / F0))
    Fact2 <- ifelse(F0 > (Nm - 1) / Nm | Fn > (sum.correl - 1) / sum.correl, 0, (1 - Fn) * log((1 - Fn) / (1 - F0)))
    G <- Fact1 + Fact2
    G_p <- G[RES >= 0]
    G_n <- G[RES < 0]
    ES_p <- max(G_p)
    ES_n <- max(G_n)
    if (ES_p > ES_n) {
      ES <- signif(ES_p, digits = 5)
      arg.ES <- which.min(abs(G - ES_p))
    } else {
      ES <- -signif(ES_n, digits = 5)
      arg.ES <- which.min(abs(G - ES_n))
    }
    return(list(ES = ES, RES = RES, arg.ES = arg.ES, indicator = tag.indicator))
  } else if (statistic == "area.under.RES") {
    # Area under the RES score
    RES <- Fn - F0
    max.ES <- max(RES)
    min.ES <- min(RES)
    if (max.ES > -min.ES) {
      arg.ES <- which.max(RES)
    } else {
      arg.ES <- which.min(RES)
    }
    ES <- sum(RES)
    return(list(ES = ES, arg.ES = arg.ES, RES = RES, indicator = tag.indicator))
  } else if (statistic == "Wilcoxon") {
    # Wilcoxon test score
    library(exactRankTests)
    seq.index <- seq(1, N)
    gene.set.ranks <- seq.index[tag.indicator == 1]
    gene.set.comp.ranks <- seq.index[tag.indicator == 0]
    W <- wilcox.exact(
      x = gene.set.ranks, y = gene.set.comp.ranks, alternative = "two.sided", mu = 0,
      paired = FALSE, exact = F, conf.int = T, conf.level = 0.95
    )
    ES <- log(1 / W$p.value)
    return(list(ES = ES, arg.ES = NULL, RES = NULL, indicator = tag.indicator))
  }
}


GSEA.EnrichmentScore6 <- function(
    #
    # Computes the weighted GSEA score of gene.set in gene.list.
    #
    # This version supports multiple statistics:
    #  "Kolmogorov-Smirnov" "Cramer-von-Mises" "Anderson-Darling" "Zhang_A" "Zhang_C" "Zhang_K" and "area.under.RES"
    #
    # The weighted score type is the exponent of the correlation
    # weight: 0 (unweighted = Kolmogorov-Smirnov), 1 (weighted), and 2 (over-weighted). When the score type is 1 or 2 it is
    # necessary to input the correlation vector with the values in the same order as in the gene list.
    #
    # Inputs:
    #   gene.list: The ordered gene list (e.g. integers indicating the original position in the input dataset)
    #   gene.set: A gene set (e.g. integers indicating the location of those genes in the input dataset)
    #   weighted.score.type: "Kolmogorov-Smirnov" "Cramer-von-Mises" "Anderson-Darling" "Zhang_A" "Zhang_C" "Zhang_K", "area.under.RES", or "Wilcoxon"
    #  correl.vector: A vector with the coorelations (e.g. signal to noise scores) corresponding to the genes in the gene list
    #
    # Outputs:
    #   ES: Enrichment score (e.g. real number between -1 and +1 for KS)
    #   arg.ES: Location in gene.list where the peak running enrichment occurs (peak of the "mountain")
    #   RES: Numerical vector containing the running enrichment score for all locations in the gene list
    #   tag.indicator: Binary vector indicating the location of the gene sets (1's) in the gene list
    #
    # The Broad Institute
    # SOFTWARE COPYRIGHT NOTICE AGREEMENT
    # This software and its documentation are copyright 2008 by the
    # Broad Institute/Massachusetts Institute of Technology.
    # All rights are reserved.
    #
    # This software is supplied without any warranty or guaranteed support
    # whatsoever. Neither the Broad Institute nor MIT can be responsible for
    # its use, misuse, or functionality.
    gene.list, # The ordered gene list (e.g. integers indicating the original position in the input dataset)
    gene.set, # A gene set (e.g. integers indicating the location of those genes in the input dataset)
    statistic = "Kolmogorov-Smirnov", # "Kolmogorov-Smirnov", "Cramer-von-Mises", "Anderson-Darling", "Zhang_A", "Zhang_C",
    #  "Zhang_K", "area.under.RES", or "Wilcoxon"
    alpha = 1, # The weight exponent
    correl.vector = NULL # A correlation vector of genes with phenotype or other appropriate weight
    ) {
  tag.indicator <- sign(match(gene.list, gene.set, nomatch = 0)) # notice that the sign is 0 (no tag) or 1 (tag)
  no.tag.indicator <- 1 - tag.indicator
  N <- length(gene.list)
  Nh <- length(gene.set)
  Nm <- N - Nh
  orig.correl.vector <- correl.vector
  if (alpha == 0) correl.vector <- rep(1, N) # unweighted case
  correl.vector <- abs(correl.vector)^alpha
  sum.correl <- sum(correl.vector[tag.indicator == 1])
  P0 <- no.tag.indicator / Nm
  F0 <- cumsum(P0)
  Pn <- tag.indicator * correl.vector / sum.correl
  Fn <- cumsum(Pn)
  if (statistic == "Kolmogorov-Smirnov") {
    RES <- Fn - F0
    max.ES <- max(RES)
    min.ES <- min(RES)
    if (max.ES > -min.ES) {
      ES <- signif(max.ES, digits = 5)
      arg.ES <- which.max(RES)
    } else {
      ES <- signif(min.ES, digits = 5)
      arg.ES <- which.min(RES)
    }
    return(list(ES = ES, arg.ES = arg.ES, RES = RES, indicator = tag.indicator))
  } else if (statistic == "Cramer-von-Mises.old") {
    # Based on Choulakian et al Canadian J. of Statistics 22, 125 (1994)
    # Modified to separate positive and negative enrichment parts
    RES <- Fn - F0
    X <- RES^2 * P0
    X_p <- X[RES >= 0]
    X_n <- X[RES < 0]
    ES_p <- sqrt(sum(X_p) / N)
    ES_n <- sqrt(sum(X_n) / N)
    if (ES_p > ES_n) {
      ES <- signif(ES_p, digits = 5)
      arg.ES <- which.min(abs(X - max(X_p)))
    } else {
      ES <- -signif(ES_n, digits = 5)
      arg.ES <- which.min(abs(X - max(X_n)))
    }
    return(list(ES = ES, RES = RES, arg.ES = arg.ES, indicator = tag.indicator))
  } else if (statistic == "Cramer-von-Mises") {
    # Based on T. W. Anderson Annals of Mathematical Statistics 1962.
    # Modified to separate positive and negative enrichment parts
    RES <- Fn - F0
    X <- RES^2
    X_p <- X[RES >= 0]
    X_n <- X[RES < 0]
    ES_p <- ((Nh * Nm) / N) * sum(X_p)
    ES_n <- ((Nh * Nm) / N) * sum(X_n)
    if (ES_p > ES_n) {
      ES <- signif(ES_p, digits = 5)
      arg.ES <- which.min(abs(X - max(X_p)))
    } else {
      ES <- -signif(ES_n, digits = 5)
      arg.ES <- which.min(abs(X - max(X_n)))
    }
    return(list(ES = ES, RES = RES, arg.ES = arg.ES, indicator = tag.indicator))
  } else if (statistic == "Anderson-Darling") {
    # Based on Choulakian et al Canadian J. of Statistics 22, 125 (1994)
    # Modified to separate positive and negative enrichment parts
    RES <- Fn - F0
    F0_factor <- ifelse(F0 < 1 / Nm | F0 > (Nm - 1) / Nm, rep(1, N), F0 * (1 - F0))
    X <- RES^2 * P0 / F0_factor
    X_p <- X[RES >= 0]
    X_n <- X[RES < 0]
    ES_p <- sqrt(sum(X_p) / N)
    ES_n <- sqrt(sum(X_n) / N)
    if (ES_p > ES_n) {
      ES <- signif(ES_p, digits = 5)
      arg.ES <- which.min(abs(X - max(X_p)))
    } else {
      ES <- -signif(ES_n, digits = 5)
      arg.ES <- which.min(abs(X - max(X_n)))
    }
    return(list(ES = ES, RES = RES, arg.ES = arg.ES, indicator = tag.indicator))
  } else if (statistic == "Zhang_A") {
    # Based on Zhang, J. R. Statist. Soc. B 64, Part2, 281 (2002)
    # Modified to separate positive and negative enrichment parts
    RES <- Fn - F0
    Fact1 <- ifelse(F0 < 1 / Nm | Fn < 1 / sum.correl, 0, Fn * log(Fn / F0))
    Fact2 <- ifelse(F0 > (Nm - 1) / Nm | Fn > (sum.correl - 1) / sum.correl, 0, (1 - Fn) * log((1 - Fn) / (1 - F0)))
    Fn_factor <- ifelse(Fn < 1 / sum.correl | Fn > (sum.correl - 1) / sum.correl, rep(1, N), Fn * (1 - Fn))
    G <- (Fact1 + Fact2) * Pn / Fn_factor
    G_p <- G[RES >= 0]
    G_n <- G[RES < 0]
    ES_p <- sum(G_p) / N
    ES_n <- sum(G_n) / N
    if (ES_p > ES_n) {
      ES <- signif(ES_p, digits = 5)
      arg.ES <- which.min(abs(G - max(G_p)))
    } else {
      ES <- -signif(ES_n, digits = 5)
      arg.ES <- which.min(abs(G - max(G_n)))
    }
    return(list(ES = ES, RES = RES, arg.ES = arg.ES, indicator = tag.indicator))
  } else if (statistic == "Zhang_C") {
    # Based on Zhang, J. R. Statist. Soc. B 64, Part2, 281 (2002)
    # Modified to separate positive and negative enrichment parts
    RES <- Fn - F0
    #      Fact1 <- ifelse(F0 < 1/Nm | Fn < 1/sum.correl, 0, Fn * log(Fn/F0))
    Fact1 <- ifelse(F0 < 1 / Nm | Fn == 0, rep(0, N), Fn * log(Fn / F0))

    #      Fact2 <- ifelse(F0 > (Nm - 1)/Nm | Fn > (sum.correl - 1)/sum.correl,
    #                            rep(0, N), (1 - Fn) * log( (1 - Fn) / (1 - F0) ))

    Fact2 <- ifelse(F0 > (Nm - 1) / Nm | Fn > 0.999, rep(0, N), (1 - Fn) * log((1 - Fn) / (1 - F0)))

    #      F0_factor <- ifelse(F0 < 1/Nm | F0 > (Nm - 1)/Nm, rep(1, N), F0 * (1 - F0))

    F0_factor <- ifelse(F0 < 1 / Nm | F0 == 1, rep(1, N), F0 * (1 - F0))


    G <- (Fact1 + Fact2) * P0 / F0_factor
    G_p <- G[RES >= 0]
    G_n <- G[RES < 0]
    ES_p <- sum(G_p) / N
    ES_n <- sum(G_n) / N
    if (ES_p > ES_n) {
      ES <- signif(ES_p, digits = 5)
      arg.ES <- which.min(abs(G - max(G_p)))
    } else {
      ES <- -signif(ES_n, digits = 5)
      arg.ES <- which.min(abs(G - max(G_n)))
    }
    return(list(ES = ES, RES = RES, arg.ES = arg.ES, indicator = tag.indicator))
  } else if (statistic == "Zhang_C.test") {
    # Based on Zhang, J. R. Statist. Soc. B 64, Part2, 281 (2002)
    # Modified to separate positive and negative enrichment parts
    RES <- Fn - F0
    Fact1 <- Fn * log(Fn / F0)
    Fact2 <- (1 - Fn) * log((1 - Fn) / (1 - F0))
    F0_factor <- F0 * (1 - F0)
    G <- (Fact1 + Fact2) * P0 / F0_factor
    ES <- signif(sum(G) / N, digits = 5)
    arg.ES <- which.min(abs(G - max(G)))
    return(list(ES = ES, RES = RES, arg.ES = arg.ES, indicator = tag.indicator))
  } else if (statistic == "Zhang_K") {
    # Based on Zhang, J. R. Statist. Soc. B 64, Part2, 281 (2002)
    # Modified to separate positive and negative enrichment parts
    RES <- Fn - F0
    Fact1 <- ifelse(F0 < 1 / Nm | Fn < 1 / sum.correl, 0, Fn * log(Fn / F0))
    Fact2 <- ifelse(F0 > (Nm - 1) / Nm | Fn > (sum.correl - 1) / sum.correl, 0, (1 - Fn) * log((1 - Fn) / (1 - F0)))
    G <- Fact1 + Fact2
    G_p <- G[RES >= 0]
    G_n <- G[RES < 0]
    ES_p <- max(G_p)
    ES_n <- max(G_n)
    if (ES_p > ES_n) {
      ES <- signif(ES_p, digits = 5)
      arg.ES <- which.min(abs(G - ES_p))
    } else {
      ES <- -signif(ES_n, digits = 5)
      arg.ES <- which.min(abs(G - ES_n))
    }
    return(list(ES = ES, RES = RES, arg.ES = arg.ES, indicator = tag.indicator))
  } else if (statistic == "area.under.RES") {
    # Area under the RES score
    RES <- Fn - F0
    max.ES <- max(RES)
    min.ES <- min(RES)
    if (max.ES > -min.ES) {
      arg.ES <- which.max(RES)
    } else {
      arg.ES <- which.min(RES)
    }
    ES <- sum(RES)
    return(list(ES = ES, arg.ES = arg.ES, RES = RES, indicator = tag.indicator))
  } else if (statistic == "Wilcoxon") {
    # Wilcoxon test score
    library(exactRankTests)
    seq.index <- seq(1, N)
    gene.set.ranks <- seq.index[tag.indicator == 1]
    gene.set.comp.ranks <- seq.index[tag.indicator == 0]
    W <- wilcox.exact(
      x = gene.set.ranks, y = gene.set.comp.ranks, alternative = "two.sided", mu = 0,
      paired = FALSE, exact = F, conf.int = T, conf.level = 0.95
    )
    ES <- log(1 / W$p.value)
    return(list(ES = ES, arg.ES = NULL, RES = NULL, indicator = tag.indicator))
  }
}

MSIG.HeatMapPlot.6 <- function(
    V,
    row.names = "NA",
    row.names2 = "NA",
    col.labels = "NA",
    col.labels2 = "NA",
    col.classes = "NA",
    phen.cmap = "NA",
    col.names = "NA",
    phen.names = "NA",
    main = " ",
    sub = " ",
    xlab = " ",
    ylab = " ",
    row.norm = TRUE,
    char.rescale = 0.85,
    cmap.type = 1, # 1 = vintage pinkogram, 2 = scale of blues, 3 = high-resolution pinkogram for scores or probabilities [0, 1], 4 = high-resolution pinkogram for general values, 5 = color map for normalized enrichment scores, 6 = scale of red purples, 7 = scale of Oranges, 8 = scale of Greens, 9 = scale of Blues
    max.v = "NA",
    legend = T) {
  #
  # Plots a heatmap "pinkogram" of a gene expression matrix including phenotype vector and gene, sample and phenotype labels
  #
  # The Broad Institute
  # SOFTWARE COPYRIGHT NOTICE AGREEMENT
  # This software and its documentation are copyright 2003 by the
  # Broad Institute/Massachusetts Institute of Technology.
  # All rights are reserved.
  #
  # This software is supplied without any warranty or guaranteed support
  # whatsoever. Neither the Broad Institute nor MIT can be responsible for
  # its use, misuse, or functionality.

  n.rows <- length(V[, 1])
  n.cols <- length(V[1, ])
  V1 <- matrix(0, nrow = n.rows, ncol = n.cols)

  #       if ((cmap.type == 5) | (cmap.type == 3)) {
  if (cmap.type == 5) {
    row.norm <- F
  }

  if (row.norm == TRUE) {
    row.mean <- apply(V, MARGIN = 1, FUN = mean)
    row.sd <- apply(V, MARGIN = 1, FUN = sd)
    row.n <- length(V[, 1])
    for (i in 1:n.rows) {
      if (row.sd[i] == 0) {
        V1[i, ] <- 0
      } else {
        V1[i, ] <- (V[i, ] - row.mean[i]) / (0.333 * row.sd[i])
      }
      V1[i, ] <- ifelse(V1[i, ] < -4, -4, V1[i, ])
      V1[i, ] <- ifelse(V1[i, ] > 4, 4, V1[i, ])
    }
  } else {
    V1 <- V
  }

  if (cmap.type == 1) {
    mycol <- c(
      "#0000FF", "#4040FF", "#7070FF", "#8888FF", "#A9A9FF", "#D5D5FF", "#EEE5EE", "#FFAADA",
      "#FF9DB0", "#FF7080",
      "#FF5A5A", "#FF4040", "#FF0D1D"
    ) # blue-pinkogram colors. This is the 1998-vintage,
    # pre-gene cluster, original pinkogram color map
  } else if (cmap.type == 2) {
    violet.palette <- colorRampPalette(c("#400030", "white"), space = "rgb")
    mycol <- rev(violet.palette(20))

    #          mycol <- c("#FCFBFD","#F4F2F8","#F8F7FB","#EFEDF5","#E1E1EF","#E8E7F2","#DADAEB","#C6C7E1","#D0D1E6",
    #                        "#BCBDDC","#A8A6CF",
    #                        "#B2B2D6","#9E9AC8","#8A87BF","#9491C4","#807DBA","#7260AB","#796FB3","#6A51A3","#5C3596",
    #                        "#63439D","#54278F","#460D83","#4D1A89","#3F007D")
  } else if (cmap.type == 6) {
    mycol <- c(
      "#FFF7F3", "#FDE0DD", "#FCC5C0", "#FA9FB5", "#F768A1", "#DD3497", "#AE017E",
      "#7A0177", "#49006A"
    )
  } else if (cmap.type == 7) {
    mycol <- c(
      "#FFF5EB", "#FEE6CE", "#FDD0A2", "#FDAE6B", "#FD8D3C", "#F16913", "#D94801",
      "#A63603", "#7F2704"
    )
  } else if (cmap.type == 8) {
    mycol <- c(
      "#F7FCF5", "#E5F5E0", "#C7E9C0", "#A1D99B", "#74C476", "#41AB5D", "#238B45",
      "#006D2C", "#00441B"
    )
  } else if (cmap.type == 9) {
    mycol <- c(
      "#F7FBFF", "#DEEBF7", "#C6DBEF", "#9ECAE1", "#6BAED6", "#4292C6", "#2171B5",
      "#08519C", "#08306B"
    )
  } else if ((cmap.type == 3) | (cmap.type == 4) | (cmap.type == 5)) {
    mycol <- vector(length = 512, mode = "numeric")

    for (k in 1:256) {
      mycol[k] <- rgb(255, k - 1, k - 1, maxColorValue = 255)
    }
    for (k in 257:512) {
      mycol[k] <- rgb(511 - (k - 1), 511 - (k - 1), 255, maxColorValue = 255)
    }
    mycol <- rev(mycol)
  }

  ncolors <- length(mycol)

  if (cmap.type == 5) {
    if (max.v == "NA") {
      max.v <- max(max(V1), -min(V1))
    }
    V2 <- ceiling(ncolors * (V1 - (-max.v)) / (1.001 * (max.v - (-max.v))))
  } else {
    V2 <- ceiling(ncolors * (V1 - min(V1)) / (1.001 * (max(V1) - min(V1))))
  }

  if (col.labels[1] == "NA") {
    heatm <- matrix(0, nrow = n.rows + 1, ncol = n.cols)
    heatm[1:n.rows, ] <- V2[seq(n.rows, 1, -1), ]
    tot.cols <- ncolors
    if (legend == T) {
      nf <- layout(matrix(c(1, 2, 3, 0), 2, 2, byrow = T), widths = c(6, 1), heights = c(10, 1), respect = FALSE)
    } else {
      nf <- layout(matrix(c(1, 2), 2, 1, byrow = T), heights = c(8, 1), respect = FALSE)
    }
    par(mar = c(3, 13, 3, 13))
    mycol <- c(mycol, phen.cmap[1:length(col.classes)])
    image(1:n.cols, 1:n.rows, t(heatm),
      zlim = c(0, tot.cols), col = mycol, axes = FALSE,
      main = main, sub = sub, xlab = xlab, ylab = ylab
    )
    n.rows.phen <- 0
  } else {
    tot.cols <- ncolors
    if (is.vector(col.labels)) {
      heatm <- matrix(0, nrow = n.rows + 1, ncol = n.cols)
      heatm[1:n.rows, ] <- V2[seq(n.rows, 1, -1), ]
      n.rows.phen <- 1
      heatm[n.rows + 1, ] <- tot.cols + col.labels
      cols.row <- length(unique(col.labels))
      tot.cols <- tot.cols + cols.row
      phen.cmap <- phen.cmap[1:cols.row]
    } else {
      n.rows.phen <- length(col.labels[, 1])
      cols.row <- vector(length = n.rows.phen, mode = "numeric")
      heatm <- matrix(0, nrow = n.rows + n.rows.phen, ncol = n.cols)
      heatm[1:n.rows, ] <- V2[seq(n.rows, 1, -1), ]
      for (k in seq(n.rows + n.rows.phen, n.rows + 1, -1)) {
        heatm[k, ] <- tot.cols + col.labels[n.rows + n.rows.phen - k + 1, ]
        cols.row[n.rows + n.rows.phen - k + 1] <- length(unique(col.labels[n.rows + n.rows.phen - k + 1, ]))
        tot.cols <- tot.cols + cols.row[n.rows + n.rows.phen - k + 1]
        #                 print(c("col:", k, ":", tot.cols + col.labels[n.rows + n.rows.phen - k + 1,], "tot.cols:", tot.cols))
      }
      phen.cmap <- phen.cmap[1:sum(unlist(lapply(col.classes, length)))]
    }
    if (legend == T) {
      nf <- layout(matrix(c(1, 2, 3, 0), 2, 2, byrow = T), widths = c(10, 2), heights = c(6, 1), respect = FALSE)
    } else {
      nf <- layout(matrix(c(1, 2), 2, 1, byrow = T), heights = c(5, 1), respect = FALSE)
    }
    par(mar = c(3, 13, 3, 10))
    mycol <- c(mycol, phen.cmap)
    image(1:n.cols, 1:(n.rows + n.rows.phen), t(heatm),
      zlim = c(0, tot.cols), col = mycol, axes = FALSE,
      main = main, sub = sub, xlab = xlab, ylab = ylab
    )
  }

  # Add lines to separate phenotypes or subgroups

  if (col.labels2[1] != "NA") {
    groups <- split(col.labels2, col.labels2)
    len.vec <- lapply(groups, length)
    plot.div <- c(0.51, cumsum(len.vec) + 0.5)
    for (i in plot.div) {
      lines(c(i, i), c(0, n.rows + n.rows.phen + 0.48), lwd = 2, cex = 0.9, col = "black")
    }
    lines(c(0.51, n.cols + 0.49), c(0.51, 0.51), lwd = 2, cex = 0.9, col = "black")
    lines(c(0.51, n.cols + 0.49), c(n.rows + n.rows.phen + 0.48, n.rows + n.rows.phen + 0.48),
      lwd = 2,
      cex = 0.9, col = "black"
    )
    lines(c(0.51, n.cols + 0.49), c(n.rows + 0.50, n.rows + 0.50),
      lwd = 2,
      cex = 0.9, col = "black"
    )
  }
  if (row.names[1] != "NA") {
    numC <- nchar(row.names)
    size.row.char <- char.rescale * 25 / (n.rows + 20)
    for (i in 1:n.rows) {
      row.names[i] <- substr(row.names[i], 1, 40)
      row.names[i] <- paste(row.names[i], " ", sep = "")
    }
    if (phen.names[1] == "NA") {
      head.names <- paste("Class", seq(n.rows.phen, 1, -1))
    } else {
      head.names <- as.character(rev(phen.names))
    }
    row.names <- c(row.names[seq(n.rows, 1, -1)], head.names)
    #            print(paste("n.rows:", n.rows))
    #            print(paste("Phen names:", phen.names))
    #            print(paste("Head names:", head.names))
    #            print(paste("Row names:", row.names))
    axis(2,
      at = 1:(n.rows + n.rows.phen), labels = row.names, adj = 0.5, tick = FALSE, las = 1, cex.axis = size.row.char,
      font.axis = 2, line = -1
    )
  }

  if (row.names2[1] != "NA") {
    numC <- nchar(row.names2)
    size.row.char <- char.rescale * 25 / (n.rows + 20)
    for (i in 1:n.rows) {
      row.names2[i] <- substr(row.names2[i], 1, 40)
      row.names2[i] <- paste(" ", row.names2[i], sep = "")
    }
    row.names2 <- c(row.names2[seq(n.rows, 1, -1)], "     ")
    axis(4,
      at = 1:(n.rows + 1), labels = row.names2, adj = 0.5, tick = FALSE, las = 1, cex.axis = size.row.char,
      font.axis = 2, line = -1
    )
  }

  if (col.names[1] != "NA") {
    size.col.char <- char.rescale * 20 / (n.cols + 25)
    axis(1, at = 1:n.cols, labels = col.names, tick = FALSE, las = 3, cex.axis = size.col.char, font.axis = 2, line = -1)
  }

  # Phenotype Legend

  #      print("--------------------------------------------------------------------------------------------")
  #      print(c("col.classes:", col.classes))
  #      print(c("phen.cmap:", phen.cmap))
  #      print(c("cols.row:", cols.row))

  if (legend == T) {
    leg.txt <- NULL
    p.vec <- NULL
    c.vec <- NULL
    c2.vec <- NULL
    ind <- 1

    for (i in 1:n.rows.phen) {
      if (is.vector(col.labels)) {
        phen.v <- as.character(col.classes)
      } else {
        phen.v <- as.character(col.classes[[i]])
      }
      leg.txt <- c(leg.txt, as.character(rev(head.names)[i]), phen.v, "  ")
      p.vec <- c(p.vec, rep(22, cols.row[i] + 2))
      c.vec <- c(c.vec, "#FFFFFF", phen.cmap[ind:(ind + cols.row[i] - 1)], "#FFFFFF")
      c2.vec <- c(c2.vec, "#FFFFFF", rep("black", cols.row[i]), "#FFFFFF")
      ind <- ind + cols.row[i]
    }
    par(mar = c(1, 0, 1, 0))
    plot(c(0, 0), c(1, 1), xlim = c(0, 1), ylim = c(0, 1), axes = F, type = "n", xlab = "", ylab = "")
    legend(
      x = 0, y = 1, legend = leg.txt, bty = "n", xjust = 0, yjust = 1, pch = p.vec, pt.bg = c.vec, col = c2.vec,
      cex = char.rescale, pt.cex = char.rescale * 1.7
    )
  }

  # Color map legend

  #       print(c("range V=", range(V)))
  #       print(c("range V1=", range(V1)))
  #       print(c("range V2=", range(V2)))

  par(mar = c(2, 12, 2, 12))
  num.v <- 20
  range.v <- range(V2)
  incr <- (range.v[1] - range.v[2]) / (num.v - 1)
  heatm.v <- matrix(rev(seq(range.v[2], range.v[1], incr)), nrow = num.v, ncol = 1)
  image(1:num.v, 1:1, heatm.v,
    zlim = c(0, tot.cols), col = mycol, axes = FALSE,
    main = " ", sub = " ", xlab = ylab, ylab = xlab
  )
  range.v <- range(V1)
  incr <- (range.v[1] - range.v[2]) / (num.v - 1)
  heatm.v2 <- matrix(signif(rev(seq(range.v[2], range.v[1], incr)), digits = 2), nrow = num.v, ncol = 1)
  #          print(c("heatm.v2=", heatm.v2))
  axis(3, at = 1:num.v, labels = heatm.v2, adj = 0.5, tick = FALSE, las = 1, cex.axis = 0.5 * char.rescale, font.axis = 1)

  return()
}

MSIG.ReadPhenFile <- function(file = "NULL") {
  #
  # Reads a matrix of class vectors from a CLS file and defines phenotype and class labels vectors
  #  (numeric and character) for the samples in a gene expression file (RES or GCT format)
  #
  # The Broad Institute
  # SOFTWARE COPYRIGHT NOTICE AGREEMENT
  # This software and its documentation are copyright 2003 by the
  # Broad Institute/Massachusetts Institute of Technology.
  # All rights are reserved.
  #
  # This software is supplied without any warranty or guaranteed support
  # whatsoever. Neither the Broad Institute nor MIT can be responsible for
  # its use, misuse, or functionality.

  cls.cont <- readLines(file)
  num.lines <- length(cls.cont)
  temp <- unlist(strsplit(cls.cont[[1]], " "))
  if (length(temp) == 3) {
    phen.names <- NULL
    col.phen <- NULL
  } else {
    l.phen.names <- match("phen.names:", temp)
    l.col.phen <- match("col.phen:", temp)
    phen.names <- temp[(l.phen.names + 1):(l.col.phen - 1)]
    col.phen <- temp[(l.col.phen + 1):length(temp)]
  }
  temp <- unlist(strsplit(cls.cont[[2]], " "))
  phen.list <- temp[2:length(temp)]

  for (k in 1:(num.lines - 2)) {
    temp <- unlist(strsplit(cls.cont[[k + 2]], " "))
    if (k == 1) {
      len <- length(temp)
      class.list <- matrix(0, nrow = num.lines - 2, ncol = len)
      class.v <- matrix(0, nrow = num.lines - 2, ncol = len)
      phen <- NULL
    }
    class.list[k, ] <- temp
    classes <- unique(temp)
    class.v[k, ] <- match(temp, classes)
    phen[[k]] <- classes
  }
  if (num.lines == 3) {
    class.list <- as.vector(class.list)
    class.v <- as.vector(class.v)
    phen <- unlist(phen)
  }
  return(list(
    phen.list = phen.list, phen = phen, phen.names = phen.names, col.phen = col.phen,
    class.v = class.v, class.list = class.list
  ))
}

MSIG.Subset.Dataset.2 <- function(
    input.ds,
    input.cls = NULL,
    column.subset = "ALL", # subset of column numbers or names (or phenotypes)
    column.sel.type = "samples", # "samples" or "phenotype"
    row.subset = "ALL", # subset of row numbers or names
    output.ds,
    output.cls = NULL) {
  # start of methodology

  print(c("Running MSIG.Subset.Dataset... on GCT file:", input.ds))
  print(c("Running MSIG.Subset.Dataset... on CLS file:", input.cls))

  # Read input datasets

  dataset <- MSIG.Gct2Frame(filename = input.ds)
  m <- data.matrix(dataset$ds)
  gs.names <- dataset$row.names
  gs.descs <- dataset$descs
  sample.names <- dataset$names

  # Read CLS file

  if (!is.null(input.cls)) {
    CLS <- MSIG.ReadPhenFile(file = input.cls)
    class.labels <- CLS$class.v
    class.phen <- CLS$phen
    class.list <- CLS$class.list
  }

  # Select desired column subset

  if (column.sel.type == "samples") {
    if (column.subset[1] == "ALL") {
      m2 <- m
      sample.names2 <- sample.names
      if (!is.null(input.cls)) {
        class.labels2 <- class.labels
      }
    } else {
      if (is.numeric(column.subset[1])) {
        m2 <- m[, column.subset]
        sample.names2 <- sample.names[column.subset]
        if (!is.null(input.cls)) {
          if (is.vector(class.labels)) {
            class.labels2 <- class.labels[column.subset]
          } else {
            class.labels2 <- class.labels[, column.subset]
          }
        }
      } else {
        locations <- !is.na(match(sample.names, column.subset))
        sample.names2 <- sample.names[locations]
        m2 <- m[, locations]
        if (!is.null(input.cls)) {
          if (is.vector(class.labels)) {
            class.labels2 <- class.labels[locations]
          } else {
            class.labels2 <- class.labels[, locations]
          }
        }
      }
    }
  } else if (column.sel.type == "phenotype") {
    locations <- !is.na(match(class.list, column.subset))
    sample.names2 <- sample.names[locations]
    m2 <- m[, locations]
    if (!is.null(input.cls)) {
      if (is.vector(class.labels)) {
        class.labels2 <- class.labels[locations]
      } else {
        class.labels2 <- class.labels[, locations]
      }
    }
  }

  if (row.subset[1] == "ALL") {
    m3 <- m2
    gs.names2 <- gs.names
    gs.descs2 <- gs.descs
  } else {
    locations <- !is.na(match(gs.names, row.subset))
    m3 <- m2[locations, ]
    gs.names2 <- gs.names[locations]
    gs.descs2 <- gs.descs[locations]
  }

  # Save datasets

  V <- data.frame(m3)
  names(V) <- sample.names2
  row.names(V) <- gs.names2
  write.gct(gct.data.frame = V, descs = gs.descs2, filename = output.ds)

  if (!is.null(input.cls)) {
    write.cls.2(class.v = class.labels2, phen = class.phen, filename = output.cls)
  }
}


MSIG.Create.CLS.from.Table <- function(
    file.gct,
    table.txt,
    output.gct = NULL,
    output.cls,
    sort.by = NULL, # column in table used for sorting (e.g. 1, 2 etc.  NULL = no sorting
    then.sort.by = NULL, # column in table used to break ties in first sort.by NULL = no sorting
    only.matches = T,
    rename = T) # create cls (and gct) with only matching entries
{
  # Read input dataset

  dataset1 <- MSIG.Gct2Frame(filename = file.gct)
  m <- data.matrix(dataset1$ds)
  gene.names <- dataset1$row.names
  gene.decs <- dataset1$descs
  sample.names.gct <- dataset1$names
  Ns <- length(sample.names.gct)

  # Read Table

  tab <- read.delim(table.txt, header = T, row.names = 1, sep = "\t", skip = 0, blank.lines.skip = T, comment.char = "", as.is = T)
  sample.names.tab <- row.names(tab)
  phen.names <- names(tab)
  print(c("phen names:", phen.names))

  # Match sample names in GCT file to sample names in table

  locs <- match(sample.names.gct, sample.names.tab, nomatch = NULL)
  print(c("locs:", locs))
  tab2 <- matrix(0, nrow = length(locs), ncol = length(tab[1, ]))
  #  tab2 <- data.matrix(tab[locs,])
  tab2 <- tab[locs, ]
  row.names.tab2 <- row.names(tab)[locs]
  cls.table <- as.matrix(t(tab2))

  print(c("tab2", tab2))

  print(c("sample.names.gct:", sample.names.gct))
  print(c("sample.names.tab:", sample.names.tab))
  print(c("Total samples GCT file:", length(sample.names.gct)))
  print(c("Total samples table file:", length(sample.names.tab)))


  if (only.matches) {
    temp <- substr(row.names.tab2, 1, 2)
    #    locs <-  temp != "NA"
    locs <- !is.na(temp)
    print(c("locs:", locs))
    if (is.vector(cls.table)) {
      cls.table <- cls.table[locs]
      print(c("Number of matches:", length(cls.table)))
    } else {
      cls.table <- cls.table[, locs]
      print(c("Number of matches:", length(cls.table[1, ])))
    }
    if (!is.null(output.gct)) {
      m <- m[, locs]
      sample.names.gct <- sample.names.gct[locs]
    }
  }

  print(c("matching.names: (after)", sample.names.gct))


  if (is.vector(cls.table)) {
    name <- phen.names
    for (j in 1:length(cls.table)) {
      if (rename == T) {
        if (is.na(cls.table[j])) {
          cls.table[j] <- "UNK"
        } else if (cls.table[j] == 1) {
          cls.table[j] <- name
        } else {
          cls.table[j] <- "WT"
        }
      }
    }
    class.phen <- unique(cls.table)

    n <- length(class.phen)
    l <- length(cls.table)
    cat(l, n, "1", "\n", file = output.cls, append = FALSE, sep = " ")
    cat("#", class.phen, "\n", file = output.cls, append = TRUE, sep = " ")
    cat(cls.table, "\n", file = output.cls, append = TRUE, sep = " ")
  } else {
    class.phen <- NULL
    for (i in 1:length(cls.table[, 1])) {
      name <- row.names(cls.table)[i]
      for (j in 1:length(cls.table[1, ])) {
        if (rename == T) {
          if (is.na(cls.table[i, j])) {
            cls.table[i, j] <- "UNK"
          } else if (cls.table[i, j] == 1) {
            cls.table[i, j] <- name
          } else {
            cls.table[i, j] <- "WT"
          }
        }
      }
      class.phen[[i]] <- unique(cls.table[i, ])
    }
    n <- length(class.phen)
    l <- length(cls.table[1, ])
    cat(l, n, "1", "\n", file = output.cls, append = FALSE, sep = " ")
    cat("#", unlist(class.phen), "\n", file = output.cls, append = TRUE, sep = " ")
    for (i in 1:length(cls.table[, 1])) {
      cat(cls.table[i, ], "\n", file = output.cls, append = TRUE, sep = " ")
    }
  }

  # Sort if required

  if (!is.null(sort.by)) {
    print(cls.table)
    if (!is.null(then.sort.by)) {
      new.order <- order(cls.table[sort.by, ], cls.table[then.sort.by, ], decreasing = F)
    } else {
      new.order <- order(cls.table[sort.by, ], decreasing = F)
    }
    if (is.vector(cls.table[, new.order])) {
      cls.table <- cls.table[new.order]
    } else {
      cls.table <- cls.table[, new.order]
    }
    print(cls.table)
    print(dim(cls.table))
    sample.names.gct <- sample.names.gct[new.order]
    if (!is.null(output.gct)) {
      m <- m[, new.order]
    }
  }

  if (!is.null(output.gct)) {
    V <- data.frame(m)
    names(V) <- sample.names.gct
    row.names(V) <- gene.names
    write.gct(gct.data.frame = V, descs = gene.decs, filename = output.gct)
  }
}

MSIG.Match.Tables <- function(
    input1.txt, # input table1
    input2.txt, # input table2
    output.txt, # output merged table
    match.rows = T, # match rows or columns
    mode = "intersection") { # mode = "intersection", "union" or "match.to.first"

  # Read input tables

  tab1 <- read.delim(input1.txt,
    header = T, row.names = 1, sep = "\t", skip = 0,
    blank.lines.skip = T, comment.char = "", as.is = T
  )
  tab2 <- read.delim(input2.txt,
    header = T, row.names = 1, sep = "\t", skip = 0,
    blank.lines.skip = T, comment.char = "", as.is = T
  )

  if (match.rows == T) {
    entries.tab1 <- row.names(tab1)
    N1 <- length(entries.tab1)
    M1 <- length(names(tab1))
    entries.tab2 <- row.names(tab2)
    N2 <- length(entries.tab2)
    M2 <- length(names(tab2))

    print(c("table 1:", N1, " entries"))
    print(c("table 1:", M1, " columns"))

    print(c("table 2:", N2, " entries"))
    print(c("table 2:", M2, " columns"))

    # Match features to first dataset and create matching m2 dataset

    if (mode == "intersection") {
      entries.tab3 <- intersect(entries.tab1, entries.tab2)
      N3 <- length(entries.tab3)
      M3 <- M1 + M2

      locations1 <- match(entries.tab3, entries.tab1, nomatch = 0)
      tab1 <- tab1[locations1, ]
      locations2 <- match(entries.tab3, entries.tab2, nomatch = 0)
      tab2 <- tab2[locations2, ]
      tab3 <- cbind(tab1, tab2)
    } else if (mode == "union") {
      entries.tab3 <- union(entries.tab1, entries.tab2)
      N3 <- length(entries.tab3)
      M3 <- M1 + M2

      m3 <- matrix("UNK", nrow = N3, ncol = M3)
      tab3 <- matrix("NA", nrow = N3, ncol = M3)

      locations1 <- match(entries.tab3, entries.tab1, nomatch = 0)
      locations2 <- match(entries.tab3, entries.tab2, nomatch = 0)

      for (i in 1:N3) {
        if (locations1[i] != 0) {
          tab3[i, 1:M1] <- as.matrix(tab1[locations1[i], ])
        }
        if (locations2[i] != 0) {
          tab3[i, (M1 + 1):M3] <- as.matrix(tab2[locations2[i], ])
        }
      }
      colnames(tab3) <- c(colnames(tab1), colnames(tab2))
    } else if (mode == "match.to.first") {
      N3 <- length(entries.tab1)
      M3 <- M1 + M2

      m3 <- matrix("UNK", nrow = N3, ncol = M3)
      tab3 <- matrix("NA", nrow = N3, ncol = M3)
      locations2 <- match(entries.tab1, entries.tab2, nomatch = 0)

      for (i in 1:N3) {
        tab3[i, 1:M1] <- as.matrix(tab1[i, ])
        if (locations2[i] != 0) {
          tab3[i, (M1 + 1):M3] <- as.matrix(tab2[locations2[i], ])
        }
      }
      colnames(tab3) <- c(colnames(tab1), colnames(tab2))
      row.names(tab3) <- entries.tab1
    } else {
      stop(c("unknown mode", mode))
    }

    # Save table

    print(c("table 3:", N3, " entries"))
    print(c("table 3:", M3, " columns"))

    col.names <- paste(colnames(tab3), collapse = "\t")
    col.names <- paste("SAMPLE", col.names, sep = "\t")
    write(noquote(col.names), file = output.txt, append = F, ncolumns = length(col.names))
    write.table(tab3, file = output.txt, quote = F, col.names = F, row.names = T, append = T, sep = "\t")
  } else { # match columns

    # Match features to first dataset and create matching m2 dataset

    entries.tab1 <- colnames(tab1)
    M1 <- length(entries.tab1)
    N1 <- length(row.names(tab1))
    entries.tab2 <- colnames(tab2)
    M2 <- length(entries.tab2)
    N2 <- length(row.names(tab2))

    print(c("table 1:", N1, " rows"))
    print(c("table 1:", M1, " columns"))

    print(c("table 2:", N2, " rows"))
    print(c("table 2:", M2, " columns"))

    if (mode == "intersection") {
      entries.tab3 <- intersect(entries.tab1, entries.tab2)
      M3 <- length(entries.tab3)
      N3 <- N1 + N2

      locations1 <- match(entries.tab3, entries.tab1, nomatch = 0)
      tab1 <- tab1[, locations1]
      locations2 <- match(entries.tab3, entries.tab2, nomatch = 0)
      tab2 <- tab2[, locations2]
      tab3 <- rbind(tab1, tab2)
    } else if (mode == "union") {
      entries.tab3 <- union(entries.tab1, entries.tab2)
      M3 <- length(entries.tab3)
      N3 <- N1 + N2

      m3 <- matrix("UNK", nrow = N3, ncol = M3)
      tab3 <- matrix("NA", nrow = N3, ncol = M3)

      locations1 <- match(entries.tab3, entries.tab1, nomatch = 0)
      locations2 <- match(entries.tab3, entries.tab2, nomatch = 0)

      for (i in 1:M3) {
        if (locations1[i] != 0) {
          tab3[1:N1, i] <- as.matrix(tab1[, locations1[i]])
        }
        if (locations2[i] != 0) {
          tab3[(N1 + 1):N3, i] <- as.matrix(tab2[, locations2[i]])
        }
      }
      row.names(tab3) <- c(row.names(tab1), row.names(tab2))
    } else if (mode == "match.to.first") {
      M3 <- length(entries.tab1)
      N3 <- N1 + N2

      m3 <- matrix("UNK", nrow = N3, ncol = M3)
      tab3 <- matrix("NA", nrow = N3, ncol = M3)
      locations2 <- match(entries.tab1, entries.tab2, nomatch = 0)

      for (i in 1:M3) {
        tab3[1:N1, i] <- as.matrix(tab1[, i])
        if (locations2[i] != 0) {
          tab3[(N1 + 1):N3, i] <- as.matrix(tab2[, locations2[i]])
        }
      }
      row.names(tab3) <- c(row.names(tab1), row.names(tab2))
      colnames(tab3) <- entries.tab1
    } else {
      stop(c("unknown mode", mode))
    }

    # Save table

    print(c("table 3:", N3, " entries"))
    print(c("table 3:", M3, " columns"))

    col.names <- paste(colnames(tab3), collapse = "\t")
    col.names <- paste("SAMPLE", col.names, sep = "\t")
    write(noquote(col.names), file = output.txt, append = F, ncolumns = length(col.names))
    write.table(tab3, file = output.txt, quote = F, col.names = F, row.names = T, append = T, sep = "\t")
  }
}

MSIG.Select.UP.DN.Markers <- function(
    input.ds,
    input.cls,
    output.marker.file,
    output.marker.plot,
    output.marker.gct,
    num.of.markers = 10,
    disc.metric = "S2N.DIFF") # discrimination metric: "S2N", "RS2N", "MEAN.DIFF", "MEDIAN.DIFF"
{
  print("Running MSIG.Select.UP.DN.Markers...")

  # Feature selection using projected dataset

  if (regexpr(pattern = ".gct", input.ds) == -1) {
    dataset <- GSEA.Res2Frame(filename = input.ds)
    gs.names <- row.names(dataset)
    sample.names <- names(dataset)
    m <- data.matrix(dataset)
  } else {
    dataset <- MSIG.Gct2Frame(filename = input.ds)
    gs.names <- dataset$row.names
    gs.descs <- dataset$descs
    sample.names <- dataset$names
    m <- data.matrix(dataset$ds)
  }

  m1 <- m

  dim(m)
  Ns <- length(m[1, ])
  Ng <- length(m[, 1])

  CLS <- ReadClsFile(file = input.cls)
  class.labels <- CLS$class.v
  class.phen <- CLS$phen
  class.list <- CLS$class.list

  #  Perform one vs all selection of num.of.markers for each class

  num.of.markers <- ifelse(num.of.markers > floor(Ng / (2 * length(class.phen))), floor(Ng / (2 * length(class.phen))), num.of.markers)
  sample.molsig.sorted.subset <- matrix(0, nrow = length(class.phen) * 2 * num.of.markers, ncol = Ns)
  sample.molsig.sorted.subset.gs <- vector(length = length(class.phen) * 2 * num.of.markers, mode = "character")
  sample.molsig.sorted.s2n <- vector(length = length(class.phen) * 2 * num.of.markers, mode = "character")
  sample.molsig.sorted.class <- vector(length = length(class.phen) * 2 * num.of.markers, mode = "character")
  class.k.labels <- class.labels
  col.index <- order(class.k.labels, decreasing = F)
  class.k.labels <- class.k.labels[col.index]
  for (j in 1:Ng) {
    m1[j, ] <- m[j, col.index]
  }
  names(m1) <- sample.names

  print("Executing marker selection...")

  if (disc.metric == "S2N") { # rank genes using signal to noise ratio (S2N)
    obs.s2n <- Gene.ranking(m1, class.k.labels, method = "S2N")
    fold.changes <- fold.changes(m1, class.k.labels, method = "MEAN.DIFF", thres = 1)
  } else if (disc.metric == "RS2N") { # rank genes using robust signal to noise ratio (RS2N)
    obs.s2n <- Gene.ranking(m1, class.k.labels, method = "RS2N")
    fold.changes <- fold.changes(m1, class.k.labels, method = "MEDIAN.DIFF", thres = 1)
  } else if (disc.metric == "MEAN.DIFF") {
    obs.s2n <- Gene.ranking(m1, class.k.labels, method = "MEAN.DIFF")
    fold.changes <- fold.changes(m1, class.k.labels, method = "MEAN.DIFF", thres = 1)
  } else if (disc.metric == "MEDIAN.DIFF") {
    obs.s2n <- Gene.ranking(m1, class.k.labels, method = "MEDIAN.DIFF")
    fold.changes <- fold.changes(m1, class.k.labels, method = "MEDIAN.DIFF", thres = 1)
  }

  obs.index <- order(obs.s2n, decreasing = T)
  obs.s2n <- sort(obs.s2n, decreasing = T)
  sample.molsig.sorted <- m[obs.index, ]
  gs.names.sorted <- gs.names[obs.index]
  msig.up <- m1[obs.index[1:num.of.markers], ]
  msig.up.genes <- gs.names[obs.index[1:num.of.markers]]
  msig.up.genes <- paste(msig.up.genes, signif(fold.changes[obs.index[1:num.of.markers]], digits = 3), sep = "_")
  msig.dn <- m1[obs.index[seq(Ng, Ng - num.of.markers + 1, -1)], ]
  msig.dn.genes <- gs.names[obs.index[seq(Ng, Ng - num.of.markers + 1, -1)]]
  msig.dn.genes <- paste(msig.dn.genes, signif(1 / fold.changes[obs.index[seq(Ng, Ng - num.of.markers + 1, -1)]], digits = 3), sep = "_")

  descs.up <- paste("UP", seq(1, num.of.markers), sep = "_")
  descs.dn <- paste("DN", seq(1, num.of.markers), sep = "_")
  descs <- c(descs.up, descs.dn)

  s2n.up <- signif(obs.s2n[1:num.of.markers], digits = 4)
  s2n.dn <- signif(obs.s2n[seq(Ng, Ng - num.of.markers + 1, -1)], digits = 4)
  s2n <- c(s2n.up, s2n.dn)

  markers <- data.frame(cbind(descs, c(msig.up.genes, msig.dn.genes), s2n))
  names(markers) <- c("UP/DN Rank", "Gene", disc.metric)

  descs.up <- paste(descs.up, s2n.up, sep = "_")
  descs.dn <- paste(descs.dn, s2n.dn, sep = "_")

  print(markers)
  write.table(markers, file = output.marker.file, quote = F, row.names = F, sep = "\t")

  c1 <- c("red", "blue")

  # saving markers gct file

  V <- data.frame(rbind(msig.up, msig.dn))
  names(V) <- sample.names
  row.names(V) <- c(msig.up.genes, msig.dn.genes)
  write.gct(gct.data.frame = V, descs = descs, filename = output.marker.gct)

  msig.up.genes <- paste(msig.up.genes, descs.up, sep = "     ")
  glob.filename <- paste(output.marker.plot, ".UP", sep = "")
  x11(height = 9, width = 12)
  nf <- layout(matrix(c(1, 2), 1, 2, byrow = T), widths = c(6, 1), heights = 1, respect = FALSE)
  MSIG.HeatMapPlot.3(
    V = msig.up, row.names = msig.up.genes, col.labels = class.labels,
    col.classes = class.phen, phen.cmap = c1[1:length(class.phen)], col.names = sample.names,
    main = "Top UP Markers", xlab = " ", ylab = " ", sub = " ", row.norm = T, cmap.type = 4, char.rescale = 1.25
  )

  # legend

  leg.txt <- class.phen
  n.phen <- length(class.phen)
  p.vec <- rep(22, n.phen)
  c.vec <- c1[1:n.phen]
  par(mar = c(0, 0, 0, 0))
  plot(c(0, 0), c(1, 1), xlim = c(0, 1), ylim = c(0, 1), axes = F, type = "n", xlab = "", ylab = "")
  legend(x = 0, y = 1, legend = leg.txt, bty = "n", xjust = 0, yjust = 1, pch = p.vec, pt.bg = c.vec, col = "black", cex = 1.5, pt.cex = 1.5)
  savePlot(filename = glob.filename, type = "jpeg", device = dev.cur())


  msig.dn.genes <- paste(msig.dn.genes, descs.dn, sep = "     ")
  glob.filename <- paste(output.marker.plot, ".DN", sep = "")
  x11(height = 9, width = 12)
  nf <- layout(matrix(c(1, 2), 1, 2, byrow = T), widths = c(6, 1), heights = 1, respect = FALSE)
  MSIG.HeatMapPlot.3(
    V = msig.dn, row.names = msig.dn.genes, col.labels = class.labels,
    col.classes = class.phen, phen.cmap = c1[1:length(class.phen)], col.names = sample.names,
    main = "Top DN Markers", xlab = " ", ylab = " ", sub = " ", row.norm = T, cmap.type = 4, char.rescale = 1.25
  )

  # legend

  leg.txt <- class.phen
  n.phen <- length(class.phen)
  p.vec <- rep(22, n.phen)
  c.vec <- c1[1:n.phen]
  par(mar = c(0, 0, 0, 0))
  plot(c(0, 0), c(1, 1), xlim = c(0, 1), ylim = c(0, 1), axes = F, type = "n", xlab = "", ylab = "")
  legend(x = 0, y = 1, legend = leg.txt, bty = "n", xjust = 0, yjust = 1, pch = p.vec, pt.bg = c.vec, col = "black", cex = 1.5, pt.cex = 1.5)
  savePlot(filename = glob.filename, type = "jpeg", device = dev.cur())
}

MSIG.Make.Biplot <- function(
    m.train,
    labels.train,
    phen.train,
    m.test = NULL,
    labels.test = NULL,
    phen.test = NULL,
    row.names = NULL,
    cols = NULL) {
  if (is.null(cols)) cols <- sample(colors(), 20)

  pca <- prcomp(t(m.train), retx = TRUE, center = TRUE, scale. = TRUE)
  S1 <- pca$x[, 1]
  S2 <- pca$x[, 2]
  X1 <- pca$rotation[, 1]
  X2 <- pca$rotation[, 2]
  max.S <- max(sqrt(S1 * S1 + S2 * S2))
  max.X <- max(sqrt(X1 * X1 + X2 * X2))
  X1 <- max.S * X1 / max.X
  X2 <- max.S * X2 / max.X
  max.A <- max(max.S, max.X)
  color <- cols[class.labels]
  num.samples <- length(S1)
  num.rows <- length(m.train[, 1])
  class.labels <- labels.train
  if (min(class.labels) == 0) class.labels <- class.labels + 1
  class.phen <- phen.train

  if (!is.null(m.test)) {
    test.scores <- predict(pca, t(m.test))
    S1 <- c(pca$x[, 1], test.scores[, 1])
    S2 <- c(pca$x[, 2], test.scores[, 2])
    max.S <- max(sqrt(S1 * S1 + S2 * S2))
    max.X <- max(sqrt(X1 * X1 + X2 * X2))
    X1 <- max.S * X1 / max.X
    X2 <- max.S * X2 / max.X
    num.samples <- length(S1)
    if (min(labels.test) == 0) labels.test <- labels.test + 1
    class.labels <- c(labels.train, labels.test + max(labels.train))
    class.phen <- c(class.phen, phen.test)
  }

  x11(height = 9, width = 14)
  nf <- layout(matrix(c(1, 2), 1, 2, byrow = T), widths = c(5, 1), heights = 1, respect = FALSE)

  plot(S1, S2, xlim = c(-max.A, max.A), ylim = c(-max.A, max.A), type = "n", main = "Biplot", sub = input.ds)
  for (j in 1:num.samples) {
    if (!is.null(m.test)) {
      if (j <= length(labels.train)) {
        pch <- 22
      } else {
        pch <- 21
      }
    } else {
      pch <- 22
    }
    points(S1[j], S2[j], pch = pch, type = "p", cex = 2.0, bg = cols[class.labels[j]], col = "black")
  }
  for (j in 1:num.rows) {
    x.coor <- X1[j] * 0.925
    y.coor <- X2[j] * 0.925
    arrows(0, 0, x.coor, y.coor, lwd = 1, length = 0.15, angle = 20, col = "black")
    leg.txt <- ifelse(!is.null(row.names), row.names[j], paste("F", j, sep = ""))
    text(X1[j], X2[j], labels = leg.txt, adj = NULL, pos = NULL, offset = 0.5, vfont = NULL, cex = 1, col = "grey")
  }
  ang <- vector(length = k.proj, mode = "numeric")
  for (j in 1:k.proj) {
    ang[j] <- ifelse(atan2(X2[j], X1[j]) > 0, atan2(X2[j], X1[j]), 2 * pi + atan2(X2[j], X1[j]))
  }
  ang.index <- order(ang, decreasing = F)
  ang2 <- ang[ang.index]
  for (j in 1:num.rows) {
    if (j == num.rows) {
      angle.in.between <- (ang2[1] - ang2[j] - 2 * pi) / 2 + ang2[j] - 2 * pi
    } else {
      angle.in.between <- (ang2[j + 1] - ang2[j]) / 2 + ang2[j]
    }
    x <- max.S * cos(angle.in.between)
    y <- max.S * sin(angle.in.between)
    arrows(0, 0, x, y, lwd = 4, length = 0, lty = 3, col = "grey50")
  }

  if (!is.null(m.test)) {
    leg.txt <- c("Train:", phen.train, "Test:", phen.test)
    p.vec <- c(rep(22, length(phen.train) + 1), rep(21, length(phen.test) + 1))
    c.vec <- c("white", cols[1:length(phen.train)], "white", cols[seq(length(phen.train) + 1, length(class.phen))])
    b.vec <- c("white", rep("black", length(phen.train)), "white", rep("black", length(phen.test)))
  } else {
    leg.txt <- class.phen
    n.phen <- length(class.phen)
    p.vec <- rep(21, n.phen)
    c.vec <- cols[1:n.phen]
    b.vec <- rep("black", n.phen)
  }
  par(mar = c(0, 0, 0, 0))
  plot(c(0, 0), c(1, 1), xlim = c(0, 1), ylim = c(0, 1), axes = F, type = "n", xlab = "", ylab = "")
  legend(x = 0, y = 1, legend = leg.txt, bty = "n", xjust = 0, yjust = 1, pch = p.vec, pt.bg = c.vec, col = b.vec, cex = 1.5, pt.cex = 2)
}

Read.GeneSets.db <- function(
    gs.db,
    thres.min = 2,
    thres.max = 2000,
    gene.names = NULL) {
  temp <- readLines(gs.db)
  max.Ng <- length(temp)
  temp.size.G <- vector(length = max.Ng, mode = "numeric")
  for (i in 1:max.Ng) {
    temp.size.G[i] <- length(unlist(strsplit(temp[[i]], "\t"))) - 2
  }
  max.size.G <- max(temp.size.G)
  gs <- matrix(rep("null", max.Ng * max.size.G), nrow = max.Ng, ncol = max.size.G)
  temp.names <- vector(length = max.Ng, mode = "character")
  temp.desc <- vector(length = max.Ng, mode = "character")
  gs.count <- 1
  for (i in 1:max.Ng) {
    gene.set.size <- length(unlist(strsplit(temp[[i]], "\t"))) - 2
    gs.line <- noquote(unlist(strsplit(temp[[i]], "\t")))
    gene.set.name <- gs.line[1]
    gene.set.desc <- gs.line[2]
    gene.set.tags <- vector(length = gene.set.size, mode = "character")
    for (j in 1:gene.set.size) {
      gene.set.tags[j] <- gs.line[j + 2]
    }
    if (is.null(gene.names)) {
      existing.set <- rep(TRUE, length(gene.set.tags))
    } else {
      existing.set <- is.element(gene.set.tags, gene.names)
    }
    set.size <- length(existing.set[existing.set == T])
    if ((set.size < thres.min) || (set.size > thres.max)) next
    temp.size.G[gs.count] <- set.size
    gs[gs.count, ] <- c(gene.set.tags[existing.set], rep("null", max.size.G - temp.size.G[gs.count]))
    temp.names[gs.count] <- gene.set.name
    temp.desc[gs.count] <- gene.set.desc
    gs.count <- gs.count + 1
  }
  Ng <- gs.count - 1
  gs.names <- vector(length = Ng, mode = "character")
  gs.desc <- vector(length = Ng, mode = "character")
  size.G <- vector(length = Ng, mode = "numeric")

  gs.names <- temp.names[1:Ng]
  gs.desc <- temp.desc[1:Ng]
  size.G <- temp.size.G[1:Ng]

  return(list(N.gs = Ng, gs = gs, gs.names = gs.names, gs.desc = gs.desc, size.G = size.G, max.N.gs = max.Ng))
}

OPAM.Projection <- function(
    data.array,
    gene.names,
    n.cols,
    n.rows,
    weight = 0,
    statistic = "Kolmogorov-Smirnov", # "Kolmogorov-Smirnov", # "Kolmogorov-Smirnov", "Cramer-von-Mises",
    # "Anderson-Darling", "Zhang_A", "Zhang_C", "Zhang_K",
    # "area.under.RES", or "Wilcoxon"
    gene.set,
    nperm = 200) {
  ES.vector <- vector(length = n.cols)
  NES.vector <- vector(length = n.cols)
  p.val.vector <- vector(length = n.cols)
  correl.vector <- vector(length = n.rows, mode = "numeric")

  # Compute ES score for signatures in each sample

  #   print("Computing GSEA.....")
  phi <- array(0, c(n.cols, nperm))
  for (sample.index in 1:n.cols) {
    gene.list <- order(data.array[, sample.index], decreasing = T)

    #      print(paste("Computing observed enrichment for UP signature in sample:", sample.index, sep=" "))
    gene.set2 <- match(gene.set, gene.names)

    if (weight == 0) {
      correl.vector <- rep(1, n.rows)
    } else if (weight > 0) {
      correl.vector <- data.array[gene.list, sample.index]
    }
    GSEA.results <- GSEA.EnrichmentScore5(
      gene.list = gene.list, gene.set = gene.set2,
      statistic = statistic, alpha = weight, correl.vector = correl.vector
    )
    ES.vector[sample.index] <- GSEA.results$ES

    if (nperm == 0) {
      NES.vector[sample.index] <- ES.vector[sample.index]
      p.val.vector[sample.index] <- 1
    } else {
      for (r in 1:nperm) {
        reshuffled.gene.labels <- sample(1:n.rows)
        if (weight == 0) {
          correl.vector <- rep(1, n.rows)
        } else if (weight > 0) {
          correl.vector <- data.array[reshuffled.gene.labels, sample.index]
        }
        GSEA.results <- GSEA.EnrichmentScore5(
          gene.list = reshuffled.gene.labels, gene.set = gene.set2,
          statistic = statistic, alpha = weight, correl.vector = correl.vector
        )
        phi[sample.index, r] <- GSEA.results$ES
      }
      if (ES.vector[sample.index] >= 0) {
        pos.phi <- phi[sample.index, phi[sample.index, ] >= 0]
        if (length(pos.phi) == 0) pos.phi <- 0.5
        pos.m <- mean(pos.phi)
        NES.vector[sample.index] <- ES.vector[sample.index] / pos.m
        s <- sum(pos.phi >= ES.vector[sample.index]) / length(pos.phi)
        p.val.vector[sample.index] <- ifelse(s == 0, 1 / nperm, s)
      } else {
        neg.phi <- phi[sample.index, phi[sample.index, ] < 0]
        if (length(neg.phi) == 0) neg.phi <- 0.5
        neg.m <- mean(neg.phi)
        NES.vector[sample.index] <- ES.vector[sample.index] / abs(neg.m)
        s <- sum(neg.phi <= ES.vector[sample.index]) / length(neg.phi)
        p.val.vector[sample.index] <- ifelse(s == 0, 1 / nperm, s)
      }
    }
  }
  return(list(ES.vector = ES.vector, NES.vector = NES.vector, p.val.vector = p.val.vector))
} # end of OPAM.Projection

OPAM.write.param.line <- function(
    param,
    param.name,
    file,
    append = T) {
  #  param <- eval(parse(text=param.name))
  if (typeof(param) == "character") {
    output.line <- paste("'", param, "'", sep = "", collapse = ",")
  } else {
    output.line <- paste(param, collapse = ",")
  }
  output.line
  if (length(param) > 1) {
    output.line <- paste(param.name, paste("c(", noquote(output.line), ")", sep = ""), sep = "\t")
  } else {
    output.line <- paste(param.name, noquote(output.line), sep = "\t")
  }
  output.line
  write(output.line, file = file, append = append, ncolumns = length(param) + 1)
}

OPAM.apply.model <- function(
    input.ds,
    input.cls,
    models.dir,
    models,
    raw.score.outfile,
    norm.score.outfile,
    model.score.outfile,
    prob.outfile,
    gmt.file = NULL) { #----------------------------------------------------------------------------------------

  # Load libraries
  erf <- function(x) 2 * pnorm(x * sqrt(2)) - 1
  library(gtools)
  library(verification)
  library(ROCR)
  library(MASS)
  library(RColorBrewer)
  library(MCMCpack)

  # Read test dataset
  dataset <- MSIG.Gct2Frame(filename = input.ds) # Read gene expression dataset (GCT format)
  m <- data.matrix(dataset$ds)
  gene.names <- dataset$row.names
  gene.descs <- dataset$descs
  sample.names <- dataset$names
  Ns <- length(m[1, ])
  Ng <- length(m[, 1])
  temp <- strsplit(input.ds, split = "/") # Extract test file name
  s <- length(temp[[1]])
  test.file.name <- temp[[1]][s]
  temp <- strsplit(test.file.name, split = ".gct")
  test.file.prefix <- temp[[1]][1]

  # Test set color map
  c.test <- c(brewer.pal(n = 7, name = "Set2"), brewer.pal(n = 7, name = "Dark2"), brewer.pal(n = 7, name = "Set1"))

  if (!is.na(input.cls)) { # Read phenotype if CLS file is provided
    CLS <- MSIG.ReadClsFile(file = input.cls) # Read phenotype file (CLS format)
    class.labels <- CLS$class.v
    class.phen <- CLS$phen
    class.list <- CLS$class.list
  }

  # Loop over models

  n.models <- length(models)

  score.matrix <- matrix(0, nrow = n.models, ncol = Ns)
  norm.score.matrix <- matrix(0, nrow = n.models, ncol = Ns)
  model.score.matrix <- matrix(0, nrow = n.models, ncol = Ns)
  probability.matrix <- matrix(0, nrow = n.models, ncol = Ns)
  models.descs <- NULL


  for (model.i in 1:n.models) {
    # Read parameters from model file
    m.file <- paste(models.dir, models[model.i], ".mod", sep = "")
    file.content <- readLines(m.file, n = -1)
    len <- length(file.content)
    for (i in 1:len) {
      temp <- unlist(strsplit(file.content[[i]], "\t"))
      len.param <- length(temp)
      if (len.param == 2) {
        param.vals <- temp[2]
      } else {
        param.vals <- paste(noquote(temp[2:len.param]), collapse = ",")
        param.vals <- paste("c(", param.vals, ")", sep = "")
      }
      assignment.string <- paste(noquote(temp[1]), " <- ", param.vals, sep = "")
      print(c("executing assigment:", assignment.string))
      eval(parse(text = assignment.string))
    }


    if (!is.null(gmt.file)) { # save genes in a gmt file
      m.name <- paste(model.name, "_UP", sep = "")
      genes.string <- paste(msig.up.genes, sep = "\t", collapse = "\t")
      output.line <- paste(m.name, model.description, genes.string, sep = "\t", collapse = "")
      if (model.i == 1) {
        write(noquote(output.line), file = gmt.file, append = F, ncolumns = length(msig.up.genes3) + 2)
      } else {
        write(noquote(output.line), file = gmt.file, append = T, ncolumns = length(msig.up.genes3) + 2)
      }
      m.name <- paste(model.name, "_DN", sep = "")
      genes.string <- paste(msig.dn.genes, sep = "\t", collapse = "\t")
      output.line <- paste(m.name, model.description, genes.string, sep = "\t", collapse = "")
      write(noquote(output.line), file = gmt.file, append = T, ncolumns = length(msig.up.genes3) + 2)
    }

    # Set parameters
    set.seed(random.seed)

    # Sample normalization
    if (sample.norm.type == "rank") {
      for (j in 1:Ns) { # column rank normalization
        m[, j] <- rank(m[, j], ties.method = "average")
      }
      m <- 10000 * m / Ng
    } else if (sample.norm.type == "log.rank") {
      for (j in 1:Ns) { # column rank normalization
        m[, j] <- rank(m[, j], ties.method = "average")
      }
      m <- log(10000 * m / Ng + exp(1))
    } else if (sample.norm.type == "log") {
      m[m < 1] <- 1
      m <- log(m + exp(1))
    }

    # Control signature normalization
    if (!is.na(msig.cntrl.genes)) {
      gene.names.int <- intersect(msig.cntrl.genes, gene.names)
      locs <- match(gene.names.int, gene.names, nomatch = 0)
      msig.cntrl <- m[locs, ]
      msig.cntrl.genes <- gene.names[locs]
      msig.cntrl.descs <- gene.descs[locs]
      msig.cntrl.size <- length(locs)
      if (msig.cntrl.size < 1) {
        msig.cntrl.center <- rep(1, Ns)
      } else if (msig.cntrl.size == 1) {
        msig.cntrl.center <- msig.cntrl
      } else if (msig.cntrl.size > 1) msig.cntrl.center <- apply(msig.cntrl, MARGIN = 2, FUN = mean)
      for (i in 1:Ng) {
        m[i, ] <- m[i, ] / msig.cntrl.center
      }
    }

    # Obtain UP3 & DN3 signatures
    gene.names.int <- intersect(msig.up.genes3, gene.names)
    locs <- match(gene.names.int, gene.names, nomatch = 0)
    msig.up.test <- m[locs, ]
    msig.up.genes.test <- gene.names[locs]
    msig.up.descs.test <- gene.descs[locs]
    msig.up.size.test <- length(locs)

    gene.names.int <- intersect(msig.dn.genes3, gene.names)
    locs <- match(gene.names.int, gene.names, nomatch = 0)
    msig.dn.test <- m[locs, ]
    msig.dn.genes.test <- gene.names[locs]
    msig.dn.descs.test <- gene.descs[locs]
    msig.dn.size.test <- length(locs)

    # Plot signatures2
    x11(height = 9, width = 12)
    MSIG.HeatMapPlot.3(
      V = msig.up.test, row.names = msig.up.genes.test, col.labels = class.labels, col.classes = class.phen,
      phen.cmap = c.test, col.names = sample.names, main = paste(model.name, " UP signature test"),
      xlab = " ", ylab = " ", sub = " ", row.norm = T, cmap.type = 4, char.rescale = 1
    )

    x11(height = 9, width = 12)
    MSIG.HeatMapPlot.3(
      V = msig.dn.test, row.names = msig.dn.genes.test, col.labels = class.labels, col.classes = class.phen,
      phen.cmap = c.test, col.names = sample.names, main = paste(model.name, " DN signature test"),
      xlab = " ", ylab = " ", sub = " ", row.norm = T, cmap.type = 4, char.rescale = 1
    )

    if (!is.na(msig.cntrl.genes)) {
      x11(height = 9, width = 12)
      MSIG.HeatMapPlot.3(
        V = msig.cntrl.test, row.names = msig.cntrl.genes.test, col.labels = class.labels,
        col.classes = class.phen, phen.cmap = c.test, col.names = sample.names,
        main = paste(model.name, " CNTRL signature"), xlab = " ", ylab = " ", sub = " ", row.norm = T,
        cmap.type = 4, char.rescale = 1
      )
    }

    # Project test dataset
    OPAM <- OPAM.Projection(m, gene.names, Ns, Ng, weight, statistic, msig.up.genes.test, nperm = nperm)
    score.up <- OPAM$ES.vector
    OPAM <- OPAM.Projection(m, gene.names, Ns, Ng, weight, statistic, msig.dn.genes.test, nperm = nperm)
    score.dn <- OPAM$ES.vector
    score <- score.up - score.dn

    x11(width = 14, height = 9)
    nf <- layout(matrix(c(1, 2, 3, 0, 4, 0), 3, 2, byrow = T), widths = c(5, 1), heights = c(1, 1, 1), respect = FALSE)
    par(mar = c(2, 4, 2, 4))
    barplot(score.up,
      main = paste(model.name, " OPAM Score UP (test)", sep = ""), font.axis = 1.25, cex.lab = 1.5,
      cex.axis = 1.25, cex.names = 1.25, width = 1, space = 0, col = c.test[class.labels]
    )
    leg.txt <- class.phen
    p.vec <- rep(22, length(leg.txt))
    par(mar = c(0, 0, 0, 0))
    plot(c(0, 0), c(1, 1), xlim = c(0, 1), ylim = c(0, 1), axes = F, type = "n", xlab = "", ylab = "")
    legend(
      x = 0, y = 0.8, legend = leg.txt, bty = "n", xjust = 0, yjust = 1, pch = p.vec, pt.bg = c.test, col = "black",
      cex = 1.25, pt.cex = 2.5
    )
    par(mar = c(2, 4, 2, 4))
    barplot(score.dn,
      main = paste(model.name, " OPAM Score DOWN (test)", sep = ""), font.axis = 1.25, cex.lab = 1.5,
      cex.axis = 1.25, cex.names = 1.25, width = 1, space = 0, col = c.test[class.labels]
    )
    par(mar = c(2, 4, 2, 4))
    barplot(score,
      main = paste(model.name, " OPAM Total Score (test)", sep = ""), font.axis = 1.25, cex.lab = 1.5,
      cex.axis = 1.25, cex.names = 1.25, width = 1, space = 0, col = c.test[class.labels]
    )

    # Apply MCMC logit or probit model to test dataset
    model.formula <- "beta.0 + beta.1 * score[i]"
    model.formula
    prob.i <- matrix(0, nrow = Ns, ncol = 3)
    model.score <- vector(length = Ns, mode = "numeric")
    for (i in 1:Ns) {
      model.score[i] <- eval(parse(text = model.formula))
      if (link.function == "logit") {
        p.vec <- paste("inv.logit(x=", model.formula, ")", sep = "")
      } else if (link.function == "probit") {
        p.vec <- paste("(erf(", model.formula, ") + 1)/2", sep = "")
      } else {
        stop("Unknown link function")
      }
      val <- eval(parse(text = p.vec))
      prob.i[i, 1] <- quantile(val, probs = 0.5)
      prob.i[i, 2] <- quantile(val, probs = 0.05)
      prob.i[i, 3] <- quantile(val, probs = 0.95)
    }
    probability <- prob.i[, 1]
    xmin <- min(model.score)
    xmax <- max(model.score)
    range.x <- xmax - xmin
    n.points <- 1000
    prob.m <- matrix(0, nrow = n.points, ncol = 3)
    x.m <- vector(length = n.points, mode = "numeric")
    for (k in 1:n.points) {
      x.m[k] <- xmin + k * (range.x / n.points)
      if (link.function == "logit") {
        p.vec <- paste("inv.logit(x=", x.m[k], ")", sep = "")
      } else if (link.function == "probit") {
        p.vec <- paste("(erf(", x.m[k], ") + 1)/2", sep = "")
      } else {
        stop("Unknown link function")
      }
      val <- eval(parse(text = p.vec))
      prob.m[k, 1] <- quantile(val, probs = 0.5)
      prob.m[k, 2] <- quantile(val, probs = 0.05)
      prob.m[k, 3] <- quantile(val, probs = 0.95)
    }
    istar <- which.min(abs(0.5 - prob.m[, 1]))
    istar <- xmin + istar * (range.x / 1000)
    x.index <- order(model.score, decreasing = F)
    x.order <- model.score[x.index]
    prob.i.order <- prob.i[x.index, ]
    target.var.order <- c.test[class.labels[x.index]]
    class.labels.order <- class.labels[x.index]

    # Plot bar graph of z-scores
    boundary <- istar
    pred.class <- ifelse(prob.i.order[, 1] >= 0.5, 2, 1)
    z.range <- range(x.order)
    x11(height = 7, width = 9.5)
    nf <- layout(matrix(c(1, 2), 1, 2, byrow = T), widths = c(3.75, 1), heights = 1, respect = FALSE)
    plot(x.order, prob.i.order[, 1], sub = model.name, pch = 20, col = 0, cex = 2, xlab = "Activation Index", ylab = "Probability")
    points(x.m, prob.m[, 1], type = "l", lwd = 2, col = 1, lty = 1, cex = 1)
    points(x.m, prob.m[, 2], type = "l", col = 4, lty = 1, cex = 1)
    points(x.m, prob.m[, 3], type = "l", col = 4, lty = 1, cex = 1)
    arrows(x.order, prob.i.order[, 2], x.order, prob.i.order[, 3], col = 4, angle = 90, code = 3, length = 0.0)
    range.x <- range(x.order)
    points(range.x, c(0.5, 0.5), type = "l", lty = 3, col = 1, lwd = 2)
    points(range.x, c(-.15, -0.15), type = "l", lty = 1, col = 1, lwd = 2)
    points(c(istar, istar), c(-0.07, 1.07), type = "l", lty = 3, col = 1, lwd = 2)
    points(x.order, prob.i.order[, 1], pch = 21, bg = target.var.order, col = 1, cex = 2)
    points(x.order, rep(-0.15, length(x.order)), pch = 21, bg = target.var.order, col = 1, cex = 2)
    leg.txt <- class.phen
    p.vec <- rep(22, length(leg.txt))
    c.vec <- c.test
    par(mar = c(0, 0, 0, 0))
    plot(c(0, 0), c(1, 1), xlim = c(0, 1), ylim = c(0, 1), axes = F, type = "n", xlab = "", ylab = "")
    legend(
      x = 0, y = 0.8, legend = leg.txt, bty = "n", xjust = 0, yjust = 1, pch = p.vec, pt.bg = c.vec,
      col = "black", cex = 1.25, pt.cex = 2.5
    )

    x11(width = 14, height = 9)
    nf <- layout(matrix(c(1, 2, 3, 4, 5, 6), 3, 2, byrow = T), widths = c(5, 1), heights = c(1, 1, 1), respect = FALSE)
    par(mar = c(4, 7, 4, 7))
    MSIG.Score.Plot(
      z = score, main = paste(model.name, " Model Score (test)"), phen.cmap = c.test,
      char.rescale = 1, col.classes = class.phen, col.labels = class.labels,
      xlab = "Samples", ylab = "Score", create.window = F, create.legend = T
    )

    par(mar = c(4, 7, 4, 7))
    norm.score <- (score - min(score)) / (max(score) - min(score))
    MSIG.Score.Plot(
      z = norm.score, main = paste(model.name, " Normalized Model Score (test)"), phen.cmap = c.test,
      char.rescale = 1, col.classes = class.phen, col.labels = class.labels,
      xlab = "Samples", ylab = "Normalized Score", create.window = F, create.legend = T
    )

    par(mar = c(4, 7, 4, 7))
    MSIG.Score.Plot(
      z = prob.i[, 1], main = paste(model.name, " Probabiliy (test)"), phen.cmap = c.test, char.rescale = 1,
      col.classes = class.phen, col.labels = class.labels,
      xlab = "Samples", ylab = "Probability", create.window = F, create.legend = T
    )

    x11(width = 14, height = 9)
    MSIG.HeatMapPlot.6(
      V = rbind(score, norm.score, model.score, probability), row.names = c(
        "raw.score", "norm.score",
        "model.score", "probability"
      ),
      row.names2 = c(model.name, model.name, model.name, model.name),
      col.labels = class.labels, col.labels2 = class.labels,
      col.classes = class.phen, phen.cmap = c.test, phen.names = model.name,
      col.names = sample.names, main = model.name,
      xlab = "  ", ylab = "  ", sub = "   ", row.norm = T, cmap.type = 3,
      char.rescale = 1, legend = T
    )

    score.matrix[model.i, ] <- score
    norm.score.matrix[model.i, ] <- norm.score
    model.score.matrix[model.i, ] <- model.score
    probability.matrix[model.i, ] <- probability
    models.descs <- c(models.descs, model.description)

    graphics.off()
    gc()
  } # end of loop over models

  # Plot projections for all models

  x11(width = 14, height = 9)
  MSIG.HeatMapPlot.6(
    V = score.matrix, row.names = models, row.names2 = models.descs, col.labels = class.labels,
    col.labels2 = class.labels, col.classes = class.phen, phen.cmap = c.test,
    col.names = sample.names, main = paste(test.file.name, " - Raw Scores"), xlab = " ", ylab = " ",
    sub = "Raw Scores", row.norm = T, cmap.type = 3, char.rescale = 1, legend = T
  )

  x11(width = 14, height = 9)
  MSIG.HeatMapPlot.6(
    V = norm.score.matrix, row.names = models, row.names2 = models.descs, col.labels = class.labels,
    col.labels2 = class.labels, col.classes = class.phen, phen.cmap = c.test,
    col.names = sample.names, main = paste(test.file.name, " - Norm Scores"), xlab = " ", ylab = " ",
    sub = "Norm Scores", row.norm = T, cmap.type = 3, char.rescale = 1, legend = T
  )

  x11(width = 14, height = 9)
  MSIG.HeatMapPlot.6(
    V = model.score.matrix, row.names = models, row.names2 = models.descs, col.labels = class.labels,
    col.labels2 = class.labels, col.classes = class.phen, phen.cmap = c.test,
    col.names = sample.names, main = paste(test.file.name, " - Model Scores"), xlab = " ", ylab = " ",
    sub = "Model Scores", row.norm = T, cmap.type = 3, char.rescale = 1, legend = T
  )

  x11(width = 14, height = 9)
  MSIG.HeatMapPlot.6(
    V = probability.matrix, row.names = models, row.names2 = models.descs, col.labels = class.labels,
    col.labels2 = class.labels, col.classes = class.phen, phen.cmap = c.test,
    col.names = sample.names, main = paste(test.file.name, " - Probabilities"), xlab = " ", ylab = " ",
    sub = "Probabilities", row.norm = T, cmap.type = 3, char.rescale = 1, legend = T
  )

  # Save projections in files

  V.GCT <- data.frame(score.matrix)
  names(V.GCT) <- sample.names
  row.names(V.GCT) <- models
  write.gct(gct.data.frame = V.GCT, descs = models.descs, filename = raw.score.outfile)

  V.GCT <- data.frame(norm.score.matrix)
  names(V.GCT) <- sample.names
  row.names(V.GCT) <- models
  write.gct(gct.data.frame = V.GCT, descs = models.descs, filename = norm.score.outfile)

  V.GCT <- data.frame(model.score.matrix)
  names(V.GCT) <- sample.names
  row.names(V.GCT) <- models
  write.gct(gct.data.frame = V.GCT, descs = models.descs, filename = model.score.outfile)

  V.GCT <- data.frame(probability.matrix)
  names(V.GCT) <- sample.names
  row.names(V.GCT) <- models
  write.gct(gct.data.frame = V.GCT, descs = models.descs, filename = prob.outfile)
} # end of function

OPAM.create.model <- function(
    input.ds,
    input.cls,
    input2.ds,
    input2.cls,
    models.dir,
    target.class,
    target.class2,
    model.name,
    model.description,
    # input.signature.file <- NA    # file with gene sets
    # input.gene.set.up <- NA       # name of UP signature
    # input.gene.set.dn <- NA       # name of DN signature
    sample.norm.type = "rank", # "rank", "log.rank", "log"
    marker.disc = "MEAN.DIFF",
    top.markers.up = 20,
    top.markers.dn = 20,
    top.markers.up2 = 20,
    top.markers.dn2 = 20,
    statistic = "area.under.RES",
    weight = 0.25,
    msig.cntrl.genes = NA,
    random.seed = 12345,
    nperm = 0,
    link.function = "logit",
    burnin.iter = 5000, # number of burnin iteration in MCMC (default: 5000)
    mcmc.iter = 25000) # number of MCMC iterations (default: 25000)
{ #----------------------------------------------------------------------------------------

  # Load libraries

  erf <- function(x) 2 * pnorm(x * sqrt(2)) - 1
  library(gtools)
  library(verification)
  library(ROCR)
  library(MASS)
  library(RColorBrewer)
  library(MCMCpack)

  #  Set parameters
  c1 <- c("black", "lightgrey")
  set.seed(random.seed)
  models.file <- paste(models.dir, "/", model.name, ".mod", sep = "")

  # Read dataset
  dataset <- MSIG.Gct2Frame(filename = input.ds) # Read gene expression dataset (GCT format)
  m <- data.matrix(dataset$ds)
  gene.names <- dataset$row.names
  gene.descs <- dataset$descs
  sample.names <- dataset$names
  Ns <- length(m[1, ])
  Ng <- length(m[, 1])

  # Read phenotype and redefine classes: target.class vs. cntrl
  CLS <- MSIG.ReadClsFile(file = input.cls) # Read phenotype file (CLS format)
  class.labels <- CLS$class.v
  class.phen <- CLS$phen
  class.list <- CLS$class.list
  if (is.na(match(target.class, class.phen))) stop(c("target class is not phenotype in:", input.cls))
  for (i in 1:length(class.list)) class.labels[i] <- ifelse(class.list[i] == target.class, 1, 2)
  col.index <- order(class.labels, decreasing = F)
  for (j in 1:Ng) m[j, ] <- m[j, col.index] # reorder samples with target.class first
  sample.names <- sample.names[col.index]
  class.labels <- class.labels[col.index]
  class.list <- class.list[col.index]
  class.phen <- c(target.class, "CNTL")
  control.class <- "CNTL"

  # Sample normalization
  if (sample.norm.type == "rank") {
    for (j in 1:Ns) { # column rank normalization
      m[, j] <- rank(m[, j], ties.method = "average")
    }
    m <- 10000 * m / Ng
  } else if (sample.norm.type == "log.rank") {
    for (j in 1:Ns) { # column rank normalization
      m[, j] <- rank(m[, j], ties.method = "average")
    }
    m <- log(10000 * m / Ng + exp(1))
  } else if (sample.norm.type == "log") {
    m[m < 1] <- 1
    m <- log(m + exp(1))
  }

  # Control signature normalization
  if (!is.na(msig.cntrl.genes)) {
    gene.names.int <- intersect(msig.cntrl.genes, gene.names)
    locs <- match(gene.names.int, gene.names, nomatch = 0)
    msig.cntrl <- m[locs, ]
    msig.cntrl.genes <- gene.names[locs]
    msig.cntrl.descs <- gene.descs[locs]
    msig.cntrl.size <- length(locs)
    if (msig.cntrl.size < 1) {
      msig.cntrl.center <- rep(1, Ns)
    } else if (msig.cntrl.size == 1) {
      msig.cntrl.center <- msig.cntrl
    } else if (msig.cntrl.size > 1) msig.cntrl.center <- apply(msig.cntrl, MARGIN = 2, FUN = mean)
    for (i in 1:Ng) {
      m[i, ] <- m[i, ] / msig.cntrl.center
    }
  }

  # Obtain UP & DN signatures
  temp <- Gene.ranking(m, class.labels, method = marker.disc)
  gene.index <- order(temp, decreasing = T)
  gene.scores <- temp[gene.index]
  msig.up <- m[gene.index[1:top.markers.up], ]
  msig.up.size <- top.markers.up
  msig.up.genes <- gene.names[gene.index[1:top.markers.up]]
  msig.up.descs <- gene.descs[gene.index[1:top.markers.up]]
  msig.dn <- m[gene.index[seq(Ng, Ng - top.markers.dn + 1, -1)], ]
  msig.dn.size <- top.markers.dn
  msig.dn.genes <- gene.names[gene.index[seq(Ng, Ng - top.markers.dn + 1, -1)]]
  msig.dn.descs <- gene.descs[gene.index[seq(Ng, Ng - top.markers.dn + 1, -1)]]
  print("Signatures (up/dn) created from gene marker selection")
  print(c("msig.up.size:", msig.up.size))
  print(c("msig.up.genes:", msig.up.genes))
  # print(c("msig.up", msig.up))
  print(c("..."))
  print(c("msig.dn.size:", msig.dn.size))
  print(c("msig.dn.genes:", msig.dn.genes))
  # print(c("msig.dn", msig.dn))
  print(c("..."))
  if (!is.na(msig.cntrl.genes)) {
    print(c("msig.cntrl.size:", msig.cntrl.size))
    print(c("msig.cntrl.genes:", msig.cntrl.genes))
    #    print(c("msig.cntrl[1,]", msig.cntrl[1,]))
  }

  # Plot signatures
  x11(height = 9, width = 12)
  MSIG.HeatMapPlot.3(
    V = msig.up, row.names = msig.up.genes, col.labels = class.labels, col.classes = class.phen,
    phen.cmap = c1, col.names = sample.names, main = paste(model.name, " UP signature"),
    xlab = " ", ylab = " ", sub = " ", row.norm = T, cmap.type = 4, char.rescale = 1
  )

  x11(height = 9, width = 12)
  MSIG.HeatMapPlot.3(
    V = msig.dn, row.names = msig.dn.genes, col.labels = class.labels, col.classes = class.phen,
    phen.cmap = c1, col.names = sample.names, main = paste(model.name, " DN signature"),
    xlab = " ", ylab = " ", sub = " ", row.norm = T, cmap.type = 4, char.rescale = 1
  )

  if (!is.na(msig.cntrl.genes)) {
    x11(height = 9, width = 12)
    MSIG.HeatMapPlot.3(
      V = msig.cntrl, row.names = msig.cntrl.genes, col.labels = class.labels, col.classes = class.phen,
      phen.cmap = c1, col.names = sample.names, main = paste(model.name, " CNTRL signature"),
      xlab = " ", ylab = " ", sub = " ", row.norm = T, cmap.type = 4, char.rescale = 1
    )
  }

  ######
  # Read refinement dataset and redefine classes: target.class2 vs. cntrl2

  dataset <- MSIG.Gct2Frame(filename = input2.ds) # Read second gene expression dataset (GCT format)
  m2 <- data.matrix(dataset$ds)
  gene.names2 <- dataset$row.names
  gene.descs2 <- dataset$descs
  sample.names2 <- dataset$names
  Ns2 <- length(m2[1, ])
  Ng2 <- length(m2[, 1])

  # Extract file name2
  temp <- strsplit(input2.ds, split = "/")
  s <- length(temp[[1]])
  temp <- temp[[1]][4]
  temp <- strsplit(temp, split = ".gct")
  results.prefix2 <- temp[[1]][1]
  results.prefix2

  CLS <- MSIG.ReadClsFile(file = input2.cls) # Read phenotype file2 (CLS format)
  class.labels2 <- CLS$class.v
  class.phen2 <- CLS$phen
  class.list2 <- CLS$class.list
  if (is.na(match(target.class2, class.phen2))) stop(c("target class2 is not phenotype in:", input2.cls))
  for (i in 1:length(class.list2)) class.labels2[i] <- ifelse(class.list2[i] == target.class2, 1, 2)
  col.index <- order(class.labels2, decreasing = F)
  for (j in 1:Ng2) m2[j, ] <- m2[j, col.index] # reorder samples with target.class first
  sample.names2 <- sample.names2[col.index]
  class.labels2 <- class.labels2[col.index]
  class.list2 <- class.list2[col.index]
  class.phen2 <- c(target.class2, "CNTL")
  control.class2 <- "CNTL"

  # Sample normalization2
  if (sample.norm.type == "rank") {
    for (j in 1:Ns2) { # column rank normalization
      m2[, j] <- rank(m2[, j], ties.method = "average")
    }
    m2 <- 10000 * m2 / Ng2
  } else if (sample.norm.type == "log.rank") {
    for (j in 1:Ns2) { # column rank normalization
      m2[, j] <- rank(m2[, j], ties.method = "average")
    }
    m2 <- log(10000 * m2 / Ng2 + exp(1))
  } else if (sample.norm.type == "log") {
    m2[m2 < 1] <- 1
    m2 <- log(m2 + exp(1))
  }

  # Control signature normalization2

  if (!is.na(msig.cntrl.genes)) {
    gene.names.int2 <- intersect(msig.cntrl.genes, gene.names2)
    locs <- match(gene.names.int2, gene.names2, nomatch = 0)
    msig.cntrl2 <- m2[locs, ]
    msig.cntrl.genes2 <- gene.names2[locs]
    msig.cntrl.descs2 <- gene.descs2[locs]
    msig.cntrl.size2 <- length(locs)
    if (msig.cntrl.size2 < 1) {
      msig.cntrl.center2 <- rep(1, Ns2)
    } else if (msig.cntrl.size2 == 1) {
      msig.cntrl.center2 <- msig.cntrl2
    } else if (msig.cntrl.size2 > 1) msig.cntrl.center2 <- apply(msig.cntrl2, MARGIN = 2, FUN = mean)
    for (i in 1:Ng2) {
      m2[i, ] <- m2[i, ] / msig.cntrl.center2
    }
  }

  gene.names.int2 <- intersect(msig.up.genes, gene.names2)
  locs <- match(gene.names.int2, gene.names2, nomatch = 0)
  msig.up2 <- m2[locs, ]
  msig.up.genes2 <- gene.names2[locs]
  msig.up.descs2 <- gene.descs2[locs]
  msig.up.size2 <- length(locs)

  gene.names.int2 <- intersect(msig.dn.genes, gene.names2)
  locs <- match(gene.names.int2, gene.names2, nomatch = 0)
  msig.dn2 <- m2[locs, ]
  msig.dn.genes2 <- gene.names2[locs]
  msig.dn.descs2 <- gene.descs2[locs]
  msig.dn.size2 <- length(locs)

  # Plot signatures2
  x11(height = 9, width = 12)
  MSIG.HeatMapPlot.3(
    V = msig.up2, row.names = msig.up.genes2, col.labels = class.labels2, col.classes = class.phen2,
    phen.cmap = c1, col.names = sample.names2, main = paste(model.name, " UP signature2"),
    xlab = " ", ylab = " ", sub = " ", row.norm = T, cmap.type = 4, char.rescale = 1
  )

  x11(height = 9, width = 12)
  MSIG.HeatMapPlot.3(
    V = msig.dn2, row.names = msig.dn.genes2, col.labels = class.labels2, col.classes = class.phen2,
    phen.cmap = c1, col.names = sample.names2, main = paste(model.name, " DN signature2"),
    xlab = " ", ylab = " ", sub = " ", row.norm = T, cmap.type = 4, char.rescale = 1
  )

  if (!is.na(msig.cntrl.genes)) {
    x11(height = 9, width = 12)
    MSIG.HeatMapPlot.3(
      V = msig.cntrl2, row.names = msig.cntrl.genes2, col.labels = class.labels2, col.classes = class.phen2,
      phen.cmap = c1, col.names = sample.names2, main = paste(model.name, " CNTRL signature2"),
      xlab = " ", ylab = " ", sub = " ", row.norm = T, cmap.type = 4, char.rescale = 1
    )
  }

  # Refine signature by second marker selection on dataset2 (using only marker genes)

  m3 <- rbind(msig.up2, msig.dn2)
  Ng3 <- length(m3[, 1])
  gene.names3 <- c(msig.up.genes2, msig.dn.genes2)
  gene.descs3 <- c(msig.up.descs2, msig.dn.descs2)

  # Obtain UP3 & DN3 signatures
  temp <- Gene.ranking(m3, class.labels2, method = marker.disc)
  gene.index <- order(temp, decreasing = T)
  gene.scores3 <- temp[gene.index]
  msig.up3 <- m3[gene.index[1:top.markers.up2], ]
  msig.up.size3 <- top.markers.up2
  msig.up.genes3 <- gene.names3[gene.index[1:top.markers.up2]]
  msig.up.descs3 <- gene.descs3[gene.index[1:top.markers.up2]]
  msig.dn3 <- m3[gene.index[seq(Ng3, Ng3 - top.markers.dn2 + 1, -1)], ]
  msig.dn.size3 <- top.markers.dn2
  msig.dn.genes3 <- gene.names3[gene.index[seq(Ng3, Ng3 - top.markers.dn2 + 1, -1)]]
  msig.dn.descs3 <- gene.descs3[gene.index[seq(Ng3, Ng3 - top.markers.dn2 + 1, -1)]]
  print("Signatures3 (up/dn) created from gene marker selection")
  print(c("msig.up.size3:", msig.up.size3))
  print(c("msig.up.genes3:", msig.up.genes3))
  # print(c("msig.up3", msig.up3))
  print(c("..."))
  print(c("msig.dn.size3:", msig.dn.size3))
  print(c("msig.dn.genes3:", msig.dn.genes3))
  # print(c("msig.dn3", msig.dn3))
  print(c("..."))

  # Plot signatures3
  x11(height = 9, width = 12)
  MSIG.HeatMapPlot.3(
    V = msig.up3, row.names = msig.up.genes3, col.labels = class.labels2, col.classes = class.phen2,
    phen.cmap = c1, col.names = sample.names2, main = paste(model.name, " UP signature3"),
    xlab = " ", ylab = " ", sub = " ", row.norm = T, cmap.type = 4, char.rescale = 1
  )

  x11(height = 9, width = 12)
  MSIG.HeatMapPlot.3(
    V = msig.dn3, row.names = msig.dn.genes3, col.labels = class.labels2, col.classes = class.phen2,
    phen.cmap = c1, col.names = sample.names2, main = paste(model.name, " DN signature3"),
    xlab = " ", ylab = " ", sub = " ", row.norm = T, cmap.type = 4, char.rescale = 1
  )


  # Project refinement dataset
  OPAM <- OPAM.Projection(m2, gene.names2, Ns2, Ng2, weight, statistic, msig.up.genes3, nperm = nperm)
  score.up <- OPAM$ES.vector
  # score.up <- sign(OPAM$ES.vector)/OPAM$p.val.vector
  OPAM <- OPAM.Projection(m2, gene.names2, Ns2, Ng2, weight, statistic, msig.dn.genes3, nperm = nperm)
  score.dn <- OPAM$ES.vector
  # score.dn <- sign(OPAM$ES.vector)/OPAM$p.val.vector
  score <- score.up - score.dn

  x11(width = 14, height = 9)
  nf <- layout(matrix(c(1, 2, 3, 0, 4, 0), 3, 2, byrow = T), widths = c(5, 1), heights = c(1, 1, 1), respect = FALSE)
  par(mar = c(2, 4, 2, 4))
  barplot(score.up,
    main = "OPAM Score UP (refinement)", font.axis = 1.25, cex.lab = 1.5, cex.axis = 1.25, cex.names =
      1.25, width = 1, space = 0, col = c1[class.labels]
  )
  leg.txt <- class.phen
  p.vec <- rep(22, length(leg.txt))
  par(mar = c(0, 0, 0, 0))
  plot(c(0, 0), c(1, 1), xlim = c(0, 1), ylim = c(0, 1), axes = F, type = "n", xlab = "", ylab = "")
  legend(
    x = 0, y = 0.8, legend = leg.txt, bty = "n", xjust = 0, yjust = 1, pch = p.vec, pt.bg = c1, col = "black",
    cex = 1.25, pt.cex = 2.5
  )
  par(mar = c(2, 4, 2, 4))
  barplot(score.dn,
    main = "OPAM Score DOWN (refinement)", font.axis = 1.25, cex.lab = 1.5, cex.axis = 1.25,
    cex.names = 1.25, width = 1, space = 0, col = c1[class.labels]
  )
  par(mar = c(2, 4, 2, 4))
  barplot(score,
    main = "OPAM Total Score (refinement)", font.axis = 1.25, cex.lab = 1.5, cex.axis = 1.25,
    cex.names = 1.25, width = 1, space = 0, col = c1[class.labels]
  )

  # Fit MCMC logit or probit model to model score using the refinement dataset
  target.var <- ifelse(class.list2 == target.class2, 1, 0)
  Bayesian.function <- ifelse(link.function == "logit", "MCMClogit(", "MCMCprobit(")
  model.formula <- paste(Bayesian.function,
    "target.var ~ score,  burnin = burnin.iter, mcmc = mcmc.iter, bayes.resid=T)",
    sep = ""
  )
  model.formula
  reg.model <- eval(parse(text = model.formula))
  beta.0 <- reg.model[, 1]
  beta.1 <- reg.model[, 2]
  model.formula <- "beta.0 + beta.1 * score[i]"
  model.formula
  prob.i <- matrix(0, nrow = Ns2, ncol = 3)
  model.score <- vector(length = Ns2, mode = "numeric")
  for (i in 1:Ns2) {
    model.score[i] <- eval(parse(text = model.formula))
    if (link.function == "logit") {
      p.vec <- paste("inv.logit(x=", model.formula, ")", sep = "")
    } else if (link.function == "probit") {
      p.vec <- paste("(erf(", model.formula, ") + 1)/2", sep = "")
    } else {
      stop("Unknown link function")
    }
    val <- eval(parse(text = p.vec))
    prob.i[i, 1] <- quantile(val, probs = 0.5)
    prob.i[i, 2] <- quantile(val, probs = 0.05)
    prob.i[i, 3] <- quantile(val, probs = 0.95)
  }
  xmin <- min(model.score)
  xmax <- max(model.score)
  range.x <- xmax - xmin
  n.points <- 1000
  prob.m <- matrix(0, nrow = n.points, ncol = 3)
  x.m <- vector(length = n.points, mode = "numeric")
  for (k in 1:n.points) {
    x.m[k] <- xmin + k * (range.x / n.points)
    if (link.function == "logit") {
      p.vec <- paste("inv.logit(x=", x.m[k], ")", sep = "")
    } else if (link.function == "probit") {
      p.vec <- paste("(erf(", x.m[k], ") + 1)/2", sep = "")
    } else {
      stop("Unknown link function")
    }
    val <- eval(parse(text = p.vec))
    prob.m[k, 1] <- quantile(val, probs = 0.5)
    prob.m[k, 2] <- quantile(val, probs = 0.05)
    prob.m[k, 3] <- quantile(val, probs = 0.95)
  }
  istar <- which.min(abs(0.5 - prob.m[, 1]))
  istar <- xmin + istar * (range.x / 1000)
  x.index <- order(model.score, decreasing = F)
  x.order <- model.score[x.index]
  prob.i.order <- prob.i[x.index, ]
  target.var.order <- ifelse(target.var[x.index] == 1, c1[1], c1[2])
  class.labels.order <- class.labels[x.index]

  # Plot bar graph of z-scores
  boundary <- istar
  pred.class <- ifelse(prob.i.order[, 1] >= 0.5, 2, 1)
  z.range <- range(x.order)
  x11(height = 7, width = 9.5)
  nf <- layout(matrix(c(1, 2), 1, 2, byrow = T), widths = c(3.75, 1), heights = 1, respect = FALSE)
  plot(x.order, prob.i.order[, 1], sub = model.name, pch = 20, col = 0, cex = 2, xlab = "Activation Index", ylab = "Probability")
  points(x.m, prob.m[, 1], type = "l", lwd = 2, col = 1, lty = 1, cex = 1)
  points(x.m, prob.m[, 2], type = "l", col = 4, lty = 1, cex = 1)
  points(x.m, prob.m[, 3], type = "l", col = 4, lty = 1, cex = 1)
  arrows(x.order, prob.i.order[, 2], x.order, prob.i.order[, 3], col = 4, angle = 90, code = 3, length = 0.0)
  range.x <- range(x.order)
  points(range.x, c(0.5, 0.5), type = "l", lty = 3, col = 1, lwd = 2)
  points(range.x, c(-.15, -0.15), type = "l", lty = 1, col = 1, lwd = 2)
  points(c(istar, istar), c(-0.07, 1.07), type = "l", lty = 3, col = 1, lwd = 2)
  points(x.order, prob.i.order[, 1], pch = 21, bg = target.var.order, col = 1, cex = 2)
  points(x.order, rep(-0.15, length(x.order)), pch = 21, bg = target.var.order, col = 1, cex = 2)
  leg.txt <- class.phen
  p.vec <- rep(22, length(leg.txt))
  c.vec <- c1
  par(mar = c(0, 0, 0, 0))
  plot(c(0, 0), c(1, 1), xlim = c(0, 1), ylim = c(0, 1), axes = F, type = "n", xlab = "", ylab = "")
  legend(
    x = 0, y = 0.8, legend = leg.txt, bty = "n", xjust = 0, yjust = 1, pch = p.vec, pt.bg = c.vec,
    col = "black", cex = 1.25, pt.cex = 2.5
  )

  x11(width = 14, height = 9)
  nf <- layout(matrix(c(1, 2, 3, 4, 5, 6), 3, 2, byrow = T), widths = c(5, 1), heights = c(1, 1, 1), respect = FALSE)
  par(mar = c(4, 7, 4, 7))
  MSIG.Score.Plot(
    z = score, main = paste(model.name, " Model Score (train)"), phen.cmap = c1,
    char.rescale = 1, col.classes = class.phen, col.labels = class.labels,
    xlab = "Samples", ylab = "Score", create.window = F, create.legend = T
  )

  par(mar = c(4, 7, 4, 7))
  norm.score <- (score - min(score)) / (max(score) - min(score))
  MSIG.Score.Plot(
    z = norm.score, main = paste(model.name, " Normalized Model Score (train)"), phen.cmap = c1,
    char.rescale = 1, col.classes = class.phen, col.labels = class.labels,
    xlab = "Samples", ylab = "Normalized Score", create.window = F, create.legend = T
  )

  par(mar = c(4, 7, 4, 7))
  MSIG.Score.Plot(
    z = prob.i[, 1], main = paste(model.name, " Probabiliy (train)"), phen.cmap = c1, char.rescale = 1,
    col.classes = class.phen, col.labels = class.labels,
    xlab = "Samples", ylab = "Probability", create.window = F, create.legend = T
  )


  MSIG.HeatMapPlot.6(
    V = rbind(score, model.score, t(prob.i[, 1])), row.names = c("raw.score", "model.score", "probability"),
    row.names2 = c(model.name, model.name, model.name),
    col.labels = class.labels2, col.labels2 = class.labels2,
    col.classes = class.phen2, phen.cmap = c1, phen.names = model.name,
    col.names = sample.names2, main = model.name,
    xlab = "  ", ylab = "  ", sub = "   ", row.norm = T, cmap.type = 3,
    char.rescale = 1, legend = T
  )

  # Save parameters in model file for annotation and to be used by apply function

  model.creation.date <- date()
  OPAM.write.param.line(param = model.creation.date, param.name = "model.creation.date", file = models.file, append = F)
  OPAM.write.param.line(param = input.ds, param.name = "input.ds", file = models.file, append = T)
  OPAM.write.param.line(param = input.ds, param.name = "input.ds", file = models.file, append = T)
  OPAM.write.param.line(param = input.cls, param.name = "input.cls", file = models.file, append = T)
  OPAM.write.param.line(param = input2.ds, param.name = "input2.ds", file = models.file, append = T)
  OPAM.write.param.line(param = input2.cls, param.name = "input2.cls", file = models.file, append = T)
  OPAM.write.param.line(param = target.class, param.name = "target.class", file = models.file, append = T)
  OPAM.write.param.line(param = target.class2, param.name = "target.class2", file = models.file, append = T)
  OPAM.write.param.line(param = model.name, param.name = "model.name", file = models.file, append = T)
  OPAM.write.param.line(param = model.description, param.name = "model.description", file = models.file, append = T)
  OPAM.write.param.line(param = sample.norm.type, param.name = "sample.norm.type", file = models.file, append = T)
  OPAM.write.param.line(param = marker.disc, param.name = "marker.disc", file = models.file, append = T)
  OPAM.write.param.line(param = top.markers.up, param.name = "top.markers.up", file = models.file, append = T)
  OPAM.write.param.line(param = top.markers.dn, param.name = "top.markers.dn", file = models.file, append = T)
  OPAM.write.param.line(param = top.markers.up2, param.name = "top.markers.up2", file = models.file, append = T)
  OPAM.write.param.line(param = top.markers.dn2, param.name = "top.markers.dn2", file = models.file, append = T)
  OPAM.write.param.line(param = statistic, param.name = "statistic", file = models.file, append = T)
  OPAM.write.param.line(param = weight, param.name = "weight", file = models.file, append = T)
  OPAM.write.param.line(param = random.seed, param.name = "random.seed", file = models.file, append = T)
  OPAM.write.param.line(param = nperm, param.name = "nperm", file = models.file, append = T)
  OPAM.write.param.line(param = link.function, param.name = "link.function", file = models.file, append = T)
  OPAM.write.param.line(param = c1, param.name = "c1", file = models.file, append = T)
  OPAM.write.param.line(param = msig.cntrl.genes, param.name = "msig.cntrl.genes", file = models.file, append = T)
  OPAM.write.param.line(param = msig.up.genes, param.name = "msig.up.genes", file = models.file, append = T)
  OPAM.write.param.line(param = msig.dn.genes, param.name = "msig.dn.genes", file = models.file, append = T)
  OPAM.write.param.line(param = msig.up.genes2, param.name = "msig.up.genes2", file = models.file, append = T)
  OPAM.write.param.line(param = msig.dn.genes2, param.name = "msig.dn.genes2", file = models.file, append = T)
  OPAM.write.param.line(param = msig.up.genes3, param.name = "msig.up.genes3", file = models.file, append = T)
  OPAM.write.param.line(param = msig.dn.genes3, param.name = "msig.dn.genes3", file = models.file, append = T)
  OPAM.write.param.line(param = beta.0, param.name = "beta.0", file = models.file, append = T)
  OPAM.write.param.line(param = beta.1, param.name = "beta.1", file = models.file, append = T)
} # end of function


OPAM.analyze.projection <- function(
    input.ds,
    input.cls,
    results.dir,
    top.class,
    top.phen,
    normalize.score = T,
    normalization.type = zero.one,
    feature.sel.thres = 0.05,
    markers.num = 5,
    user.colors = NA,
    k.proj = NA,
    markers.metric = "S2N",
    markers.file = NULL) {
  library(gtools)

  library(verification)
  library(ROCR)
  library(MASS)
  library(RColorBrewer)
  library(heatmap.plus)

  dataset <- MSIG.Gct2Frame(filename = input.ds) # Read gene expression dataset (GCT format)
  m <- data.matrix(dataset$ds)
  model.names <- dataset$row.names
  model.descs <- dataset$descs
  sample.names <- dataset$names
  Ns <- length(m[1, ])
  n.models <- length(m[, 1])
  temp <- strsplit(input.ds, split = "/") # Extract test file name
  s <- length(temp[[1]])
  test.file.name <- temp[[1]][s]
  temp <- strsplit(test.file.name, split = ".gct")

  test.file.prefix <- temp[[1]][1]
  char.res <- 0.013 * n.models + 0.65

  # normalize scores
  if (normalize.score == T) {
    if (normalization.type == "zero.one") {
      for (i in 1:n.models) {
        m[i, ] <- (m[i, ] - min(m[i, ])) / (max(m[i, ]) - min(m[i, ]))
      }
    } else if (normalization.type == "z.score") {
      for (i in 1:n.models) {
        m[i, ] <- (m[i, ] - mean(m[i, ])) / sd(m[i, ])
      }
    } else if (normalization.type == "r.z.score") {
      for (i in 1:n.models) {
        m[i, ] <- (m[i, ] - median(m[i, ])) / mad(m[i, ])
      }
    }
  }
  CLS <- MSIG.ReadPhenFile(file = input.cls) # Read phenotype file (CLS format)
  cls.labels <- CLS$class.v
  cls.phen <- CLS$phen
  cls.list <- CLS$class.list

  if (is.vector(cls.labels)) {
    cls.labels2 <- cls.labels
    n.phen <- 1
  } else {
    cls.labels2 <- as.vector(cls.labels[top.phen, ])
    n.phen <- length(cls.labels[, 1])
  }
  if (!is.na(user.colors)) {
    c.test <- user.colors
  } else {
    if (!is.null(CLS$col.phen)) {
      c.test <- CLS$col.phen
    } else {
      c.test <- c(brewer.pal(n = 7, name = "Set2"), brewer.pal(n = 7, name = "Dark2"), brewer.pal(n = 7, name = "Set1"))
    }
  }
  if (!is.null(CLS$phen.names)) {
    phen.names <- CLS$phen.names
  } else {
    phen.names <- "NA"
  }

  print("CLS:")
  print(CLS)

  n.classes <- vector(length = n.phen, mode = "numeric")
  if (n.phen == 1) {
    max.classes <- length(cls.phen)
    n.classes[1] <- max.classes
  } else {
    max.classes <- max(unlist(lapply(cls.phen, FUN = length)))
    for (i in 1:n.phen) {
      n.classes[i] <- length(cls.phen[[i]])
    }
  }

  print("computing ROC...")

  roc.array <- array(data = 0, dim = c(n.models, n.phen, max.classes), dimnames = NULL)
  p.val.array <- array(data = 0, dim = c(n.models, n.phen, max.classes), dimnames = NULL)

  for (i in 1:n.models) {
    for (j in 1:n.phen) {
      for (k in 1:n.classes[[j]]) {
        if (is.vector(cls.labels)) {
          bin.class <- ifelse(cls.labels == k, 1, 0)
        } else {
          bin.class <- ifelse(cls.labels[j, ] == k, 1, 0)
        }
        m.score <- m[i, ]
        m.score.norm <- (m.score - min(m.score)) / (max(m.score) - min(m.score))
        perf.auc <- roc.area(bin.class, m.score.norm)
        roc.array[i, j, k] <- signif(perf.auc$A, digits = 3)
        p.val.array[i, j, k] <- signif(perf.auc$p.value, digits = 3)
      }
    }
  }

  print("computing ROC for top class")

  top.roc.vector <- vector(length = n.models, mode = "numeric")
  top.p.val.vector <- vector(length = n.models, mode = "numeric")
  for (i in 1:n.models) {
    if (is.vector(cls.labels)) {
      bin.class <- ifelse(cls.list == top.class, 1, 0)
    } else {
      bin.class <- ifelse(cls.list[top.phen, ] == top.class, 1, 0)
    }
    if (markers.metric == "ROC") {
      m.score <- m[i, ]
      m.score.norm <- (m.score - min(m.score)) / (max(m.score) - min(m.score))
      perf.auc <- roc.area(bin.class, m.score.norm)
      top.roc.vector[i] <- signif(perf.auc$A, digits = 3)
      top.p.val.vector[i] <- signif(perf.auc$p.value, digits = 3)
    } else if (markers.metric == "S2N") {
      top.roc.vector[i] <- signif(S2N(m[i, ], bin.class), digits = 3)
      temp <- split(m[i, ], bin.class)
      x <- temp[[1]]
      y <- temp[[2]]
      t.test.out <- t.test(x = x, y = y)
      top.p.val.vector[i] <- signif(t.test.out$p.value, digits = 3)
    }
  }

  roc.table <- data.frame(cbind(model.names, model.descs, top.roc.vector, top.p.val.vector))
  roc.order <- order(top.roc.vector, decreasing = T)
  roc.table <- roc.table[roc.order, ]
  names(roc.table) <- c("Feature:", "Description:", markers.metric, "p-value:")

  #   print(roc.table)

  top.pos.features.m <- m[top.p.val.vector <= feature.sel.thres, ]
  top.pos.features <- model.names[top.p.val.vector <= feature.sel.thres]
  top.pos.features.descs <- model.descs[top.p.val.vector <= feature.sel.thres]
  top.pos.features.roc <- top.roc.vector[top.p.val.vector <= feature.sel.thres]
  top.pos.features.p.val <- top.p.val.vector[top.p.val.vector <= feature.sel.thres]

  top.neg.features.m <- m[top.p.val.vector >= 1 - feature.sel.thres, ]
  top.neg.features <- model.names[top.p.val.vector >= 1 - feature.sel.thres]
  top.neg.features.descs <- model.descs[top.p.val.vector >= 1 - feature.sel.thres]
  top.neg.features.roc <- top.roc.vector[top.p.val.vector >= 1 - feature.sel.thres]
  top.neg.features.p.val <- top.p.val.vector[top.p.val.vector >= 1 - feature.sel.thres]

  top.features <- data.frame(rbind(cbind(
    top.pos.features, top.pos.features.descs, top.pos.features.roc,
    top.pos.features.p.val
  ), cbind(
    top.neg.features, top.neg.features.descs,
    top.neg.features.roc, top.neg.features.p.val
  )))
  names(top.features) <- c("Feature:", "Description:", markers.metric, "p-value:")

  top.m <- rbind(top.pos.features.m, top.neg.features.m)
  top.features.order <- order(top.features[, 3], decreasing = T)
  top.features <- top.features[top.features.order, ]
  top.features

  # Heatmap

  print("plot in original order")

  height <- ifelse(n.models + n.phen >= 9, 9, (n.models + n.phen) * 0.44 + 5)

  #   print(c("char.res:", char.res))
  x11(width = 14, height = height)
  MSIG.HeatMapPlot.6(
    V = m, row.names = model.names, row.names2 = model.descs,
    col.labels = cls.labels, col.labels2 = cls.labels2,
    col.classes = cls.phen, phen.cmap = c.test, phen.names = phen.names,
    col.names = sample.names, main = paste(test.file.prefix, "-- Original Order"),
    xlab = "  ", ylab = "  ", sub = "Original order", row.norm = T, cmap.type = 3,
    char.rescale = char.res, legend = T
  )

  filename <- paste(results.dir, test.file.prefix, ".HEATMAP", sep = "")
  #   savePlot(filename = filename, type ="jpeg", device = dev.cur())

  pdf(file = paste(filename, ".pdf", sep = ""), height = 8.5, width = 11)

  MSIG.HeatMapPlot.6(
    V = m, row.names = model.names, row.names2 = model.descs,
    col.labels = cls.labels, col.labels2 = cls.labels2,
    col.classes = cls.phen, phen.cmap = c.test, phen.names = phen.names,
    col.names = sample.names, main = paste(test.file.prefix, "-- Original Order"),
    xlab = "  ", ylab = "  ", sub = "Original order", row.norm = T, cmap.type = 3,
    char.rescale = char.res, legend = T
  )
  dev.off()

  cls.labels.renam <- cls.labels
  if (!is.vector(cls.labels)) {
    max.lab <- max(cls.labels[, 1])
    for (k in 2:length(cls.labels[, 1])) {
      cls.labels.renam[k, ] <- cls.labels[k, ] + max.lab
      max.label <- max(cls.labels[k, ])
    }
  }

  #   print(dim(cls.labels.renam))
  #   print(dim(m))
  #   print("cls.labels.renam:")
  #   print(t(cls.labels.renam))

  # Model order AUC ROC

  row.order <- order(top.roc.vector, decreasing = TRUE)
  V1 <- m[row.order, ]
  model.names2 <- model.names[row.order]
  model.descs2 <- paste(model.descs[row.order], top.roc.vector[row.order], top.p.val.vector[row.order])
  score.index <- order(V1[1, ], decreasing = TRUE)
  cls.labels2.sorted <- cls.labels2[score.index]
  if (is.vector(cls.labels)) {
    cls.labels.sorted <- cls.labels[score.index]
  } else {
    cls.labels.sorted <- cls.labels[, score.index]
  }
  sample.names.sorted <- sample.names[score.index]

  print("plot in top model order. Best match to target class")

  height <- ifelse(n.models + n.phen >= 9, 9, (n.models + n.phen) * 0.44 + 5)
  #   char.res <-  0.0125 * n.models + 0.60
  x11(width = 14, height = height)
  MSIG.HeatMapPlot.6(
    V = V1[, score.index], row.names = model.names2,
    row.names2 = model.descs2, col.labels = cls.labels.sorted,
    col.classes = cls.phen, phen.cmap = c.test, phen.names = phen.names,
    col.names = sample.names.sorted, main = paste(test.file.prefix, "-- Model Order", markers.metric),
    xlab = "  ", ylab = "  ", sub = "Top Model order", row.norm = T, cmap.type = 3,
    char.rescale = char.res, legend = T
  )

  filename <- paste(results.dir, test.file.prefix, ".HEATMAP.MODEL", sep = "")
  #   savePlot(filename = filename, type ="jpeg", device = dev.cur())

  pdf(file = paste(filename, ".pdf", sep = ""), height = 8.5, width = 11)

  MSIG.HeatMapPlot.6(
    V = V1[, score.index], row.names = model.names2,
    row.names2 = model.descs2, col.labels = cls.labels.sorted,
    col.classes = cls.phen, phen.cmap = c.test, phen.names = phen.names,
    col.names = sample.names.sorted, main = paste(test.file.prefix, "-- Top Model Order"),
    xlab = "  ", ylab = "  ", sub = "Model order", row.norm = T, cmap.type = 3,
    char.rescale = char.res, legend = T
  )

  dev.off()

  # Sort top phenotype and then Sort inside each sub-phenotype

  dist.matrix <- dist(m)
  HC <- hclust(dist.matrix, method = "complete")
  V2 <- m[HC$order, ]
  model.names2 <- model.names[HC$order]
  model.descs2 <- model.descs[HC$order]
  top.roc.vector2 <- top.roc.vector[HC$order]
  top.p.val.vector2 <- top.p.val.vector[HC$order]

  sample.names2 <- sample.names
  num.phen <- length(unique(cls.labels2))
  for (k in 1:num.phen) {
    V3 <- V2[, cls.labels2 == k]
    sn <- sample.names2[cls.labels2 == k]
    if (is.vector(cls.labels)) {
      cl <- cls.labels[cls.labels2 == k]
    } else {
      cl <- cls.labels[, cls.labels2 == k]
    }
    dist.matrix <- dist(t(V3))
    HC <- hclust(dist.matrix, method = "complete")
    V3 <- V3[, HC$order]
    sn <- sn[HC$order]
    if (is.vector(cls.labels)) {
      cl <- cl[HC$order]
      cls.labels[cls.labels2 == k] <- cl
    } else {
      cl <- cl[, HC$order]
      cls.labels[, cls.labels2 == k] <- cl
    }
    V2[, cls.labels2 == k] <- V3
    sample.names2[cls.labels2 == k] <- sn
  }

  print("plot sorted inside each phenotype")
  model.descs2 <- paste(model.descs2, top.roc.vector2, top.p.val.vector2)
  height <- ifelse(n.models + n.phen >= 9, 9, (n.models + n.phen) * 0.44 + 5)
  x11(width = 14, height = height)
  #     char.res <-  0.015 * n.models + 0.725
  MSIG.HeatMapPlot.6(
    V = V2, row.names = model.names2, row.names2 = model.descs2,
    col.labels = cls.labels, col.labels2 = cls.labels2,
    col.classes = cls.phen, phen.cmap = c.test, phen.names = phen.names,
    col.names = sample.names2, main = paste(test.file.prefix, "-- Sorted Inside Class"),
    xlab = "  ", ylab = "  ", sub = "Sorted Inside Class", row.norm = T, cmap.type = 3,
    char.rescale = char.res, legend = T
  )

  filename <- paste(results.dir, test.file.prefix, ".HEATMAP.SORT.PHEN", sep = "")
  #     savePlot(filename = filename, type ="jpeg", device = dev.cur())

  pdf(file = paste(filename, ".pdf", sep = ""), height = 8.5, width = 11)

  MSIG.HeatMapPlot.6(
    V = V2, row.names = model.names2, row.names2 = model.descs2,
    col.labels = cls.labels, col.labels2 = cls.labels2,
    col.classes = cls.phen, phen.cmap = c.test, phen.names = phen.names,
    col.names = sample.names2, main = paste(test.file.prefix, "-- Sorted Inside Class"),
    xlab = "  ", ylab = "  ", sub = "Sorted Inside Class", row.norm = T, cmap.type = 3,
    char.rescale = char.res, legend = T
  )
  dev.off()

  # Markers for each class

  print("plot in marker order")

  if (is.vector(cls.labels)) {
    classes <- unique(cls.list)
  } else {
    classes <- unique(cls.list[top.phen, ])
  }
  if (length(classes) > 2) {
    only.up <- T
  } else {
    only.up <- F
  }
  if (markers.num == 0) {
    if (is.vector(cls.labels)) {
      bin.class <- ifelse(cls.list == top.class, 0, 1)
    } else {
      bin.class <- ifelse(cls.list[top.phen, ] == top.class, 0, 1)
    }
    if (markers.metric == "S2N") {
      metric <- Gene.ranking(m, bin.class, method = "RS2N")
    } else if (markers.metric == "ROC") {
      bin.class <- ifelse(bin.class == 1, 0, 1)
      metric <- vector(length = n.models, mode = "numeric")
      for (i in 1:n.models) {
        m.score <- m[i, ]
        m.score.norm <- (m.score - min(m.score)) / (max(m.score) - min(m.score))
        perf.auc <- roc.area(bin.class, m.score.norm)
        metric[i] <- perf.auc$A
      }
    }
    metric <- signif(metric, digits = 3)
    metric.order <- order(metric, decreasing = T)
    markers <- model.names[metric.order]
    markers.descs <- model.descs[metric.order]
    metric.list <- metric[metric.order]
    markers.num <- length(m[, 1]) / 2
    k.class <- rep(top.class, markers.num)
  } else {
    if (length(classes) == 2) classes <- classes[1]
    markers <- NULL
    markers.descs <- NULL
    metric.list <- NULL
    k.class <- NULL
    for (k in classes) {
      if (is.vector(cls.labels)) {
        bin.class <- ifelse(cls.list == k, 0, 1)
      } else {
        bin.class <- ifelse(cls.list[top.phen, ] == k, 0, 1)
      }
      if (markers.metric == "S2N") {
        metric <- Gene.ranking(m, bin.class, method = "RS2N")
      } else if (markers.metric == "ROC") {
        bin.class <- ifelse(bin.class == 1, 0, 1)
        metric <- vector(length = n.models, mode = "numeric")
        for (i in 1:n.models) {
          m.score <- m[i, ]
          m.score.norm <- (m.score - min(m.score)) / (max(m.score) - min(m.score))
          perf.auc <- roc.area(bin.class, m.score.norm)
          metric[i] <- perf.auc$A
        }
      }
      metric <- signif(metric, digits = 2)
      metric.order <- order(metric, decreasing = T)
      if (only.up == T) {
        markers <- c(markers, model.names[metric.order][1:markers.num])
        markers.descs <- c(markers.descs, model.descs[metric.order][1:markers.num])
        metric.list <- c(metric.list, metric[metric.order][1:markers.num])
        k.class <- c(k.class, rep(k, markers.num))
      } else {
        markers <- c(
          markers, model.names[metric.order][1:markers.num],
          model.names[metric.order][(length(model.names) - markers.num + 1):length(model.names)]
        )
        markers.descs <- c(
          markers.descs, model.descs[metric.order][1:markers.num],
          model.descs[metric.order][(length(model.names) - markers.num + 1):length(model.names)]
        )
        metric.list <- c(
          metric.list, metric[metric.order][1:markers.num],
          metric[metric.order][(length(model.names) - markers.num + 1):length(model.names)]
        )
        k.class <- c(k.class, rep(k, markers.num), rep(paste("not", k), markers.num))
      }
    }
  }
  print("markers")
  #   print(markers)
  #   V1 <- m[markers,]
  V3 <- V2[markers, ]
  model.descs2 <- paste(markers.descs, metric.list, k.class)

  height <- ifelse(length(markers) + n.phen >= 9, 9, (length(markers) + n.phen) * 0.44 + 5)
  char.res <- 0.0125 * length(markers) + 0.70

  x11(width = 15, height = height)
  MSIG.HeatMapPlot.6(
    V = V3, row.names = markers,
    row.names2 = model.descs2, col.labels = cls.labels, col.labels2 = cls.labels2,
    col.classes = cls.phen, phen.cmap = c.test, phen.names = phen.names,
    col.names = sample.names2, main = paste(test.file.prefix, "-- Markers Order"),
    xlab = "  ", ylab = "  ", sub = paste("Markers Order (metric: ", metric, ")"), row.norm = T,
    cmap.type = 3, char.rescale = char.res, legend = T
  )

  filename <- paste(results.dir, test.file.prefix, ".HEATMAP.MARKERS.", markers.metric, sep = "")
  #   savePlot(filename = filename, type ="jpeg", device = dev.cur())

  pdf(file = paste(filename, ".pdf", sep = ""), height = 14, width = 11)

  MSIG.HeatMapPlot.6(
    V = V3, row.names = markers,
    row.names2 = model.descs2, col.labels = cls.labels, col.labels2 = cls.labels2,
    col.classes = cls.phen, phen.cmap = c.test, phen.names = phen.names,
    col.names = sample.names2, main = paste(test.file.prefix, "-- Markers Order"),
    xlab = "  ", ylab = "  ", sub = paste("Markers Order (metric: ", metric, ")"), row.norm = T,
    cmap.type = 3, char.rescale = char.res, legend = T
  )
  dev.off()

  names(V3)[seq(1, Ns)] <- sample.names2
  row.names(V3) <- paste(markers, seq(1, length(markers)), sep = "_")
  if (!is.null(markers.file)) {
    write.gct(gct.data.frame = V3, descs = model.descs2, filename = markers.file)
  }

  # NMF projection

  print("NMF projection")

  if (is.na(k.proj)) {
    k.proj <- ifelse(length(model.names2) < 10, 3, ceiling(length(model.names2) / 8))
  }

  NMF.out <- NMF.div(V = V2, k = k.proj, maxniter = 4000, seed = 1234, stopconv = 40, stopfreq = 10)
  H <- NMF.out$H
  W <- NMF.out$W

  k.proj.names <- paste("NMF_", seq(1, k.proj), sep = "")
  #
  height <- ifelse(length(k.proj.names) + n.phen >= 9, 9, (length(markers) + n.phen) * 0.44 + 5)
  char.res <- 0.0125 * length(markers) + 0.70

  print("H matrix")

  x11(width = 15, height = height)

  MSIG.HeatMapPlot.6(
    V = H, row.names = k.proj.names, row.names2 = k.proj.names,
    col.labels = cls.labels, col.labels2 = cls.labels2,
    col.classes = cls.phen, phen.cmap = c.test, phen.names = phen.names,
    col.names = sample.names2, main = paste(test.file.prefix, "-- H matrix"),
    xlab = "  ", ylab = "  ", sub = "H Matrix", row.norm = T, cmap.type = 3,
    char.rescale = char.res, legend = T
  )

  filename <- paste(results.dir, test.file.prefix, ".HEATMAP.SORT.PHEN.H.MATRIX", sep = "")
  #     savePlot(filename = filename, type ="jpeg", device = dev.cur())
  pdf(file = paste(filename, ".pdf", sep = ""), height = 8.5, width = 11)

  MSIG.HeatMapPlot.6(
    V = H, row.names = k.proj.names, row.names2 = k.proj.names,
    col.labels = cls.labels, col.labels2 = cls.labels2,
    col.classes = cls.phen, phen.cmap = c.test, phen.names = phen.names,
    col.names = sample.names2, main = paste(test.file.prefix, "-- H matrix"),
    xlab = "  ", ylab = "  ", sub = "H Matrix", row.norm = T, cmap.type = 3,
    char.rescale = char.res, legend = T
  )
  dev.off()

  print("W matrix")

  height <- ifelse(length(k.proj.names) + n.phen >= 9, 9, (length(markers) + n.phen) * 0.44 + 5)
  char.res <- 0.0125 * length(markers) + 0.70

  x11(width = 15, height = height)
  MSIG.HeatMapPlot.6(
    V = W, row.names = model.names2, row.names2 = model.descs2,
    col.labels = seq(1, k.proj),
    col.classes = k.proj.names, phen.cmap = c.test,
    col.names = k.proj.names, main = paste(test.file.prefix, "-- W matrix"),
    xlab = "  ", ylab = "  ", sub = "W Matrix", row.norm = T, cmap.type = 3,
    char.rescale = char.res, legend = T
  )

  filename <- paste(results.dir, test.file.prefix, ".HEATMAP.SORT.PHEN.W.MATRIX", sep = "")
  #     savePlot(filename = filename, type ="jpeg", device = dev.cur())
  pdf(file = paste(filename, ".pdf", sep = ""), height = 8.5, width = 11)

  MSIG.HeatMapPlot.6(
    V = W, row.names = model.names2, row.names2 = model.descs2,
    col.labels = seq(1, k.proj),
    col.classes = k.proj.names, phen.cmap = c.test,
    col.names = k.proj.names, main = paste(test.file.prefix, "-- W matrix"),
    xlab = "  ", ylab = "  ", sub = "W Matrix", row.norm = T, cmap.type = 3,
    char.rescale = char.res, legend = T
  )

  dev.off()

  # sort entire set

  print("sorting entire set")

  dist.matrix <- dist(m)
  HC <- hclust(dist.matrix, method = "complete")
  V2 <- m[HC$order, ]
  model.names2 <- model.names[HC$order]
  model.descs2 <- model.descs[HC$order]
  sample.names2 <- sample.names
  dist.matrix <- dist(t(V2))
  HC <- hclust(dist.matrix, method = "complete")
  V2 <- V2[, HC$order]
  sample.names.sorted <- sample.names2[HC$order]
  if (is.vector(cls.labels)) {
    cls.labels.sorted <- cls.labels[HC$order]
  } else {
    cls.labels.sorted <- cls.labels[, HC$order]
  }

  print("plot bisorted")

  height <- ifelse(n.models + n.phen >= 9, 9, (n.models + n.phen) * 0.44 + 5)
  x11(width = 14, height = height)
  #     char.res <-  0.015 * n.models + 0.725
  MSIG.HeatMapPlot.6(
    V = V2, row.names = model.names2, row.names2 = model.descs2,
    col.labels = cls.labels.sorted,
    col.classes = cls.phen, phen.cmap = c.test, phen.names = phen.names,
    col.names = sample.names.sorted, main = paste(test.file.prefix, "-- Bisorted"),
    xlab = "  ", ylab = "  ", sub = "Bisorted", row.norm = F, cmap.type = 3,
    char.rescale = char.res, legend = T
  )

  filename <- paste(results.dir, test.file.prefix, ".HEATMAP.BISORT", sep = "")
  #   savePlot(filename = filename, type ="jpeg", device = dev.cur())

  pdf(file = paste(filename, ".pdf", sep = ""), height = 8.5, width = 11)

  MSIG.HeatMapPlot.6(
    V = V2, row.names = model.names2, row.names2 = model.descs2,
    col.labels = cls.labels.sorted,
    col.classes = cls.phen, phen.cmap = c.test, phen.names = phen.names,
    col.names = sample.names.sorted, main = paste(test.file.prefix, "-- Bisorted"),
    xlab = "  ", ylab = "  ", sub = "Bisorted", row.norm = F, cmap.type = 3,
    char.rescale = char.res, legend = T
  )

  dev.off()
}


OPAM.analyze.projection.2 <- function(
    input.ds,
    input.cls,
    results.dir,
    top.class,
    top.phen,
    normalize.score = T,
    normalization.type = zero.one,
    feature.sel.thres = 0.05,
    markers.num = 5,
    user.colors = NA,
    k.proj = NA,
    markers.metric = "ROC", # "ROC" or "T.TEST"
    markers.file = NULL,
    markers.file.cls = NULL) {
  library(gtools)
  library(verification)
  library(ROCR)
  library(MASS)
  library(RColorBrewer)
  library(heatmap.plus)

  dataset <- MSIG.Gct2Frame(filename = input.ds) # Read gene expression dataset (GCT format)
  m <- data.matrix(dataset$ds)
  model.names <- dataset$row.names
  model.descs <- dataset$descs
  dim(m)
  sample.names <- dataset$names
  Ns <- length(m[1, ])
  n.models <- length(m[, 1])
  temp <- strsplit(input.ds, split = "/") # Extract test file name
  s <- length(temp[[1]])
  test.file.name <- temp[[1]][s]
  temp <- strsplit(test.file.name, split = ".gct")

  test.file.prefix <- temp[[1]][1]
  char.res <- 0.013 * n.models + 0.65

  for (i in length(m[, 1])) {
    if (sd(m[i, ]) == 0) {
      val <- m[i, 1]
      m[i, ] <- m[i, ] + runif(n = Ns, min = val - 0.001, max = val + 0.001) # add small noise to flat profiles
    }
  }

  # normalize scores
  if (normalize.score == T) {
    if (normalization.type == "zero.one") {
      for (i in 1:n.models) {
        m[i, ] <- (m[i, ] - min(m[i, ])) / (max(m[i, ]) - min(m[i, ]))
      }
    } else if (normalization.type == "z.score") {
      for (i in 1:n.models) {
        m[i, ] <- (m[i, ] - mean(m[i, ])) / sd(m[i, ])
      }
    } else if (normalization.type == "r.z.score") {
      for (i in 1:n.models) {
        m[i, ] <- (m[i, ] - median(m[i, ])) / mad(m[i, ])
      }
    }
  }
  CLS <- MSIG.ReadPhenFile(file = input.cls) # Read phenotype file (CLS format)
  cls.labels <- CLS$class.v
  cls.phen <- CLS$phen
  cls.list <- CLS$class.list
  print(CLS)

  if (is.vector(cls.labels)) {
    cls.labels2 <- cls.labels
    n.phen <- 1
  } else {
    cls.labels2 <- as.vector(cls.labels[top.phen, ])
    n.phen <- length(cls.labels[, 1])
  }
  if (!is.na(user.colors)) {
    c.test <- user.colors
  } else {
    if (!is.null(CLS$col.phen)) {
      c.test <- CLS$col.phen
    } else {
      c.test <- c(brewer.pal(n = 7, name = "Set1"), brewer.pal(n = 7, name = "Dark2"), brewer.pal(n = 7, name = "Set2"))
    }
  }
  if (!is.null(CLS$phen.names)) {
    phen.names <- CLS$phen.names
  } else {
    phen.names <- "NA"
  }

  print("CLS:")
  print(CLS)

  n.classes <- vector(length = n.phen, mode = "numeric")
  if (n.phen == 1) {
    max.classes <- length(cls.phen)
    n.classes[1] <- max.classes
  } else {
    max.classes <- max(unlist(lapply(cls.phen, FUN = length)))
    for (i in 1:n.phen) {
      n.classes[i] <- length(cls.phen[[i]])
    }
  }

  print("computing ROC...")

  roc.array <- array(data = 0, dim = c(n.models, n.phen, max.classes), dimnames = NULL)
  p.val.array <- array(data = 0, dim = c(n.models, n.phen, max.classes), dimnames = NULL)

  for (i in 1:n.models) {
    for (j in 1:n.phen) {
      for (k in 1:n.classes[[j]]) {
        if (is.vector(cls.labels)) {
          bin.class <- ifelse(cls.labels == k, 1, 0)
        } else {
          bin.class <- ifelse(cls.labels[j, ] == k, 1, 0)
        }
        m.score <- m[i, ]
        m.score.norm <- (m.score - min(m.score)) / (max(m.score) - min(m.score))
        #           perf.auc <- roc.area(bin.class, m.score.norm)
        #           roc.array[i, j, k] <- signif(perf.auc$A, digits=3)
        #           p.val.array[i, j, k] <- signif(perf.auc$p.value, digits=3)
      }
    }
  }

  #   print("computing ROC or T.TEST stats for top class")

  top.roc.vector <- vector(length = n.models, mode = "numeric")
  top.p.val.vector <- vector(length = n.models, mode = "numeric")
  for (i in 1:n.models) {
    if (is.vector(cls.labels)) {
      bin.class <- ifelse(cls.list == top.class, 1, 0)
    } else {
      bin.class <- ifelse(cls.list[top.phen, ] == top.class, 1, 0)
    }
    if (markers.metric == "ROC") {
      m.score <- m[i, ]
      m.score.norm <- (m.score - min(m.score)) / (max(m.score) - min(m.score))
      perf.auc <- roc.area(bin.class, m.score.norm)
      top.roc.vector[i] <- signif(perf.auc$A, digits = 3)
      top.p.val.vector[i] <- signif(perf.auc$p.value, digits = 3)
    } else if (markers.metric == "T.TEST") {
      temp <- split(m[i, ], bin.class)
      x <- temp[[1]]
      y <- temp[[2]]
      top.roc.vector[i] <- signif(t.test(x = x, y = y)$statistic, digits = 3)
      top.p.val.vector[i] <- signif(t.test(x = x, y = y)$p.value, digits = 3)
    }
  }

  roc.table <- data.frame(cbind(model.names, model.descs, top.roc.vector, top.p.val.vector))
  roc.order <- order(top.roc.vector, decreasing = T)
  roc.table <- roc.table[roc.order, ]
  names(roc.table) <- c("Feature:", "Description:", markers.metric, "p-value:")

  print(roc.table)

  top.pos.features.m <- m[top.p.val.vector <= feature.sel.thres, ]
  top.pos.features <- model.names[top.p.val.vector <= feature.sel.thres]
  top.pos.features.descs <- model.descs[top.p.val.vector <= feature.sel.thres]
  top.pos.features.roc <- top.roc.vector[top.p.val.vector <= feature.sel.thres]
  top.pos.features.p.val <- top.p.val.vector[top.p.val.vector <= feature.sel.thres]

  top.neg.features.m <- m[top.p.val.vector >= 1 - feature.sel.thres, ]
  top.neg.features <- model.names[top.p.val.vector >= 1 - feature.sel.thres]
  top.neg.features.descs <- model.descs[top.p.val.vector >= 1 - feature.sel.thres]
  top.neg.features.roc <- top.roc.vector[top.p.val.vector >= 1 - feature.sel.thres]
  top.neg.features.p.val <- top.p.val.vector[top.p.val.vector >= 1 - feature.sel.thres]

  top.features <- data.frame(rbind(cbind(
    top.pos.features, top.pos.features.descs, top.pos.features.roc,
    top.pos.features.p.val
  ), cbind(
    top.neg.features, top.neg.features.descs,
    top.neg.features.roc, top.neg.features.p.val
  )))
  names(top.features) <- c("Feature:", "Description:", markers.metric, "p-value:")

  top.m <- rbind(top.pos.features.m, top.neg.features.m)
  top.features.order <- order(top.features[, 3], decreasing = T)
  top.features <- top.features[top.features.order, ]
  top.features

  # Heatmap

  print("plot in original order")

  height <- ifelse(n.models + n.phen >= 9, 9, (n.models + n.phen) * 0.44 + 5)

  #   print(c("char.res:", char.res))
  x11(width = 14, height = height)
  MSIG.HeatMapPlot.6(
    V = m, row.names = model.names, row.names2 = model.descs,
    col.labels = cls.labels,
    col.classes = cls.phen, phen.cmap = c.test, phen.names = phen.names,
    col.names = sample.names, main = paste(test.file.prefix, " - ", top.class, " - Original Order"),
    xlab = "  ", ylab = "  ", sub = "Original order", row.norm = T, cmap.type = 3,
    char.rescale = char.res, legend = T
  )

  filename <- paste(results.dir, test.file.prefix, ".", top.class, ".HEATMAP", sep = "")
  #   savePlot(filename = filename, type ="jpeg", device = dev.cur())

  pdf(file = paste(filename, ".pdf", sep = ""), height = 8.5, width = 11)

  MSIG.HeatMapPlot.6(
    V = m, row.names = model.names, row.names2 = model.descs,
    col.labels = cls.labels,
    col.classes = cls.phen, phen.cmap = c.test, phen.names = phen.names,
    col.names = sample.names, main = paste(test.file.prefix, " - ", top.class, "- Original Order."),
    xlab = "  ", ylab = "  ", sub = "Original order", row.norm = T, cmap.type = 3,
    char.rescale = char.res, legend = T
  )
  dev.off()

  cls.labels.renam <- cls.labels
  if (!is.vector(cls.labels)) {
    max.lab <- max(cls.labels[, 1])
    for (k in 2:length(cls.labels[, 1])) {
      cls.labels.renam[k, ] <- cls.labels[k, ] + max.lab
      max.label <- max(cls.labels[k, ])
    }
  }

  #   print(dim(cls.labels.renam))
  #   print(dim(m))
  #   print("cls.labels.renam:")
  #   print(t(cls.labels.renam))

  # Sort according to top phenotype

  print("cls.list:")
  print(cls.list)

  phen.index <- order(cls.labels2, decreasing = FALSE)
  if (is.vector(cls.labels)) {
    cls.labels <- cls.labels[phen.index]
    cls.list <- cls.list[phen.index]
  } else {
    cls.labels <- cls.labels[, phen.index]
    cls.list <- cls.list[, phen.index]
  }
  cls.labels2 <- cls.labels2[phen.index]
  sample.names <- sample.names[phen.index]
  m <- m[, phen.index]


  print("phen.index:")
  print(phen.index)

  print("cls.list:")
  print(cls.list)

  print("cls.labels:")
  print(cls.labels)

  print("cls.labels2:")
  print(cls.labels2)

  # Sort inside each sub-phenotype

  dist.matrix <- dist(m)
  HC <- hclust(dist.matrix, method = "complete")
  V2 <- m[HC$order, ]
  model.names2 <- model.names[HC$order]
  model.descs2 <- model.descs[HC$order]
  top.roc.vector2 <- top.roc.vector[HC$order]
  top.p.val.vector2 <- top.p.val.vector[HC$order]

  sample.names2 <- sample.names
  num.phen <- length(unique(cls.labels2))
  for (k in 1:num.phen) {
    V3 <- V2[, cls.labels2 == k]
    sn <- sample.names2[cls.labels2 == k]
    if (is.vector(cls.labels)) {
      cl <- cls.labels[cls.labels2 == k]
    } else {
      cl <- cls.labels[, cls.labels2 == k]
    }
    dist.matrix <- dist(t(V3))
    HC <- hclust(dist.matrix, method = "complete")
    V3 <- V3[, HC$order]
    sn <- sn[HC$order]
    if (is.vector(cls.labels)) {
      cl <- cl[HC$order]
      cls.labels[cls.labels2 == k] <- cl
    } else {
      cl <- cl[, HC$order]
      cls.labels[, cls.labels2 == k] <- cl
    }
    V2[, cls.labels2 == k] <- V3
    sample.names2[cls.labels2 == k] <- sn
  }

  print("plot sorted inside each phenotype")
  model.descs2 <- paste(model.descs2, top.roc.vector2, top.p.val.vector2)
  height <- ifelse(n.models + n.phen >= 9, 9, (n.models + n.phen) * 0.44 + 5)
  x11(width = 14, height = height)
  #     char.res <-  0.015 * n.models + 0.725
  MSIG.HeatMapPlot.6(
    V = V2, row.names = model.names2, row.names2 = model.descs2,
    col.labels = cls.labels,
    col.classes = cls.phen, phen.cmap = c.test, phen.names = phen.names,
    col.names = sample.names2, main = paste(test.file.prefix, " - ", top.class, "- Sorted Inside Class"),
    xlab = "  ", ylab = "  ", sub = "Sorted Inside Class", row.norm = T, cmap.type = 3,
    char.rescale = char.res, legend = T
  )

  filename <- paste(results.dir, test.file.prefix, ".", top.class, ".HEATMAP.SORT.PHEN", sep = "")
  #     savePlot(filename = filename, type ="jpeg", device = dev.cur())

  pdf(file = paste(filename, ".pdf", sep = ""), height = 8.5, width = 11)

  MSIG.HeatMapPlot.6(
    V = V2, row.names = model.names2, row.names2 = model.descs2,
    col.labels = cls.labels,
    col.classes = cls.phen, phen.cmap = c.test, phen.names = phen.names,
    col.names = sample.names2, main = paste(test.file.prefix, " - ", top.class, "- Sorted Inside Class"),
    xlab = "  ", ylab = "  ", sub = "Sorted Inside Class", row.norm = T, cmap.type = 3,
    char.rescale = char.res, legend = T
  )
  dev.off()

  # Markers for each class

  print("plot in marker order")

  if (is.vector(cls.labels)) {
    classes <- unique(cls.list)
  } else {
    classes <- unique(cls.list[top.phen, ])
  }
  if (length(classes) > 2) {
    only.up <- T
  } else {
    only.up <- F
  }
  if (markers.num == 0) { # find and sort all markers for the top class
    if (is.vector(cls.labels)) {
      bin.class <- ifelse(cls.list == top.class, 0, 1)
    } else {
      bin.class <- ifelse(cls.list[top.phen, ] == top.class, 0, 1)
    }
    if (markers.metric == "T.TEST") {
      for (i in 1:n.models) {
        temp <- split(m[i, ], bin.class)
        x <- temp[[1]]
        y <- temp[[2]]
        metric[i] <- signif(t.test(x = x, y = y)$statistic, digits = 3)
        p.val[i] <- signif(t.test(x = x, y = y)$p.value, digits = 3)
      }
    } else if (markers.metric == "ROC") {
      bin.class <- ifelse(bin.class == 1, 0, 1)
      metric <- vector(length = n.models, mode = "numeric")
      p.val <- vector(length = n.models, mode = "numeric")
      for (i in 1:n.models) {
        m.score <- m[i, ]
        m.score.norm <- (m.score - min(m.score)) / (max(m.score) - min(m.score))
        perf.auc <- roc.area(bin.class, m.score.norm)
        metric[i] <- perf.auc$A
        p.val[i] <- signif(perf.auc$p.value, digits = 3)
      }
    }
    metric <- signif(metric, digits = 3)
    metric.order <- order(metric, decreasing = T)
    markers <- model.names[metric.order]
    markers.descs <- model.descs[metric.order]
    metric.list <- metric[metric.order]
    markers.num <- length(m[, 1]) / 2
    k.class <- rep(top.class, markers.num)
  } else {
    if (length(classes) == 2) classes <- classes[1]
    markers <- NULL
    markers.descs <- NULL
    metric.list <- NULL
    p.val.list <- NULL
    k.class <- NULL
    for (k in classes) {
      if (is.vector(cls.labels)) {
        bin.class <- ifelse(cls.list == k, 0, 1)
      } else {
        bin.class <- ifelse(cls.list[top.phen, ] == k, 0, 1)
      }
      if (markers.metric == "T.TEST") {
        metric <- vector(length = n.models, mode = "numeric")
        p.val <- vector(length = n.models, mode = "numeric")
        for (i in 1:n.models) {
          temp <- split(m[i, ], bin.class)
          x <- temp[[1]]
          y <- temp[[2]]
          metric[i] <- signif(t.test(x = x, y = y)$statistic, digits = 3)
          p.val[i] <- signif(t.test(x = x, y = y)$p.value, digits = 3)
        }
      } else if (markers.metric == "ROC") {
        bin.class <- ifelse(bin.class == 1, 0, 1)
        metric <- vector(length = n.models, mode = "numeric")
        p.val <- vector(length = n.models, mode = "numeric")
        for (i in 1:n.models) {
          m.score <- m[i, ]
          m.score.norm <- (m.score - min(m.score)) / (max(m.score) - min(m.score))
          perf.auc <- roc.area(bin.class, m.score.norm)
          metric[i] <- signif(perf.auc$A, digits = 3)
          p.val[i] <- signif(perf.auc$p.value, digits = 3)
        }
      }
      metric.order <- order(metric, decreasing = T)
      if (only.up == T) {
        markers <- c(markers, model.names[metric.order][1:markers.num])
        markers.descs <- c(markers.descs, model.descs[metric.order][1:markers.num])
        metric.list <- c(metric.list, metric[metric.order][1:markers.num])
        p.val.list <- c(p.val.list, p.val[metric.order][1:markers.num])
        k.class <- c(k.class, rep(k, markers.num))
      } else {
        markers <- c(
          markers, model.names[metric.order][1:markers.num],
          model.names[metric.order][(length(model.names) - markers.num + 1):length(model.names)]
        )
        markers.descs <- c(
          markers.descs, model.descs[metric.order][1:markers.num],
          model.descs[metric.order][(length(model.names) - markers.num + 1):length(model.names)]
        )
        metric.list <- c(
          metric.list, metric[metric.order][1:markers.num],
          metric[metric.order][(length(model.names) - markers.num + 1):length(model.names)]
        )
        p.val.list <- c(
          p.val.list, p.val[metric.order][1:markers.num],
          p.val[metric.order][(length(model.names) - markers.num + 1):length(model.names)]
        )
        k.class <- c(k.class, rep(k, markers.num), rep(paste("not", k), markers.num))
      }
    }
  }
  print("markers")
  #   print(markers)
  #   V1 <- m[markers,]
  V3 <- V2[markers, ]
  model.descs2 <- paste(markers.descs, metric.list, p.val.list, k.class)

  height <- ifelse(length(markers) + n.phen >= 9, 9, (length(markers) + n.phen) * 0.44 + 5)
  char.res <- 0.0125 * length(markers) + 0.70

  x11(width = 15, height = height)
  MSIG.HeatMapPlot.6(
    V = V3, row.names = markers,
    row.names2 = model.descs2, col.labels = cls.labels,
    col.classes = cls.phen, phen.cmap = c.test, phen.names = phen.names,
    col.names = sample.names2, main = paste(test.file.prefix, " - ", top.class, "- Markers Order - ", markers.metric),
    xlab = "  ", ylab = "  ", row.norm = T,
    cmap.type = 3, char.rescale = char.res, legend = T
  )

  filename <- paste(results.dir, test.file.prefix, ".", top.class, ".HEATMAP.MARKERS.", markers.metric, sep = "")
  #   savePlot(filename = filename, type ="jpeg", device = dev.cur())

  pdf(file = paste(filename, ".pdf", sep = ""), height = 8.5, width = 11)

  MSIG.HeatMapPlot.6(
    V = V3, row.names = markers,
    row.names2 = model.descs2, col.labels = cls.labels,
    col.classes = cls.phen, phen.cmap = c.test, phen.names = phen.names,
    col.names = sample.names2, main = paste(test.file.prefix, " - ", top.class, "- Markers Order - ", markers.metric),
    xlab = "  ", ylab = "  ", row.norm = T,
    cmap.type = 3, char.rescale = char.res, legend = T
  )
  dev.off()

  V3 <- data.frame(V3)
  colnames(V3) <- sample.names2
  row.names(V3) <- paste(markers, seq(1, length(markers)), sep = "_")
  if (!is.null(markers.file)) {
    write.gct(gct.data.frame = V3, descs = model.descs2, filename = markers.file)
  }
  if (!is.null(markers.file.cls)) {
    write.cls.2(class.v = cls.labels, phen = cls.phen, filename = markers.file.cls)
  }

  # NMF projection

  print("NMF projection")

  if (is.na(k.proj)) {
    k.proj <- ifelse(length(model.names2) < 10, 3, ceiling(length(model.names2) / 8))
  }

  NMF.out <- NMF.div(V = V2, k = k.proj, maxniter = 4000, seed = 1234, stopconv = 40, stopfreq = 10)
  H <- NMF.out$H
  W <- NMF.out$W

  k.proj.names <- paste("NMF_", seq(1, k.proj), sep = "")
  #
  height <- ifelse(length(k.proj.names) + n.phen >= 9, 9, (length(markers) + n.phen) * 0.44 + 5)
  char.res <- 0.0125 * length(markers) + 0.70

  print("H matrix")

  x11(width = 15, height = height)

  MSIG.HeatMapPlot.6(
    V = H, row.names = k.proj.names, row.names2 = k.proj.names,
    col.labels = cls.labels,
    col.classes = cls.phen, phen.cmap = c.test, phen.names = phen.names,
    col.names = sample.names2, main = paste(test.file.prefix, " - ", top.class, "- H matrix"),
    xlab = "  ", ylab = "  ", sub = "H Matrix", row.norm = T, cmap.type = 3,
    char.rescale = char.res, legend = T
  )

  filename <- paste(results.dir, test.file.prefix, ".", top.class, ".HEATMAP.SORT.PHEN.H.MATRIX", sep = "")
  #     savePlot(filename = filename, type ="jpeg", device = dev.cur())

  pdf(file = paste(filename, ".pdf", sep = ""), height = 8.5, width = 11)

  MSIG.HeatMapPlot.6(
    V = H, row.names = k.proj.names, row.names2 = k.proj.names,
    col.labels = cls.labels,
    col.classes = cls.phen, phen.cmap = c.test, phen.names = phen.names,
    col.names = sample.names2, main = paste(test.file.prefix, " - ", top.class, "- H matrix"),
    xlab = "  ", ylab = "  ", sub = "H Matrix", row.norm = T, cmap.type = 3,
    char.rescale = char.res, legend = T
  )
  dev.off()

  print("W matrix")

  height <- ifelse(length(k.proj.names) + n.phen >= 9, 9, (length(markers) + n.phen) * 0.44 + 5)
  char.res <- 0.0125 * length(markers) + 0.70

  x11(width = 15, height = height)
  MSIG.HeatMapPlot.6(
    V = W, row.names = model.names2, row.names2 = model.descs2,
    col.labels = seq(1, k.proj),
    col.classes = k.proj.names, phen.cmap = seq(1, k.proj),
    col.names = k.proj.names, main = paste(test.file.prefix, " - ", top.class, "- W matrix"),
    xlab = "  ", ylab = "  ", sub = "W Matrix", row.norm = T, cmap.type = 3,
    char.rescale = char.res, legend = T
  )

  filename <- paste(results.dir, test.file.prefix, ".", top.class, ".HEATMAP.SORT.PHEN.W.MATRIX", sep = "")
  #     savePlot(filename = filename, type ="jpeg", device = dev.cur())

  pdf(file = paste(filename, ".pdf", sep = ""), height = 8.5, width = 11)

  MSIG.HeatMapPlot.6(
    V = W, row.names = model.names2, row.names2 = model.descs2,
    col.labels = seq(1, k.proj),
    col.classes = k.proj.names, phen.cmap = c.test,
    col.names = k.proj.names, main = paste(test.file.prefix, " - ", top.class, "- W matrix"),
    xlab = "  ", ylab = "  ", sub = "W Matrix", row.norm = T, cmap.type = 3,
    char.rescale = char.res, legend = T
  )

  dev.off()

  # sort entire set

  print("sorting entire set")

  dist.matrix <- dist(m)
  HC <- hclust(dist.matrix, method = "complete")
  V2 <- m[HC$order, ]
  model.names2 <- model.names[HC$order]
  model.descs2 <- model.descs[HC$order]
  sample.names2 <- sample.names
  dist.matrix <- dist(t(V2))
  HC <- hclust(dist.matrix, method = "complete")
  V2 <- V2[, HC$order]
  sample.names.sorted <- sample.names2[HC$order]
  if (is.vector(cls.labels)) {
    cls.labels.sorted <- cls.labels[HC$order]
  } else {
    cls.labels.sorted <- cls.labels[, HC$order]
  }

  print("plot bisorted")

  height <- ifelse(n.models + n.phen >= 9, 9, (n.models + n.phen) * 0.44 + 5)
  x11(width = 14, height = height)
  #     char.res <-  0.015 * n.models + 0.725
  MSIG.HeatMapPlot.6(
    V = V2, row.names = model.names2, row.names2 = model.descs2,
    col.labels = cls.labels.sorted,
    col.classes = cls.phen, phen.cmap = c.test, phen.names = phen.names,
    col.names = sample.names.sorted, main = paste(test.file.prefix, " - ", top.class, "- Bisorted"),
    xlab = "  ", ylab = "  ", sub = "Bisorted", row.norm = T, cmap.type = 3,
    char.rescale = char.res, legend = T
  )

  filename <- paste(results.dir, test.file.prefix, ".", top.class, ".HEATMAP.BISORT", sep = "")
  #   savePlot(filename = filename, type ="jpeg", device = dev.cur())

  pdf(file = paste(filename, ".pdf", sep = ""), height = 8.5, width = 11)

  MSIG.HeatMapPlot.6(
    V = V2, row.names = model.names2, row.names2 = model.descs2,
    col.labels = cls.labels.sorted,
    col.classes = cls.phen, phen.cmap = c.test, phen.names = phen.names,
    col.names = sample.names.sorted, main = paste(test.file.prefix, " - ", top.class, "- Bisorted"),
    xlab = "  ", ylab = "  ", sub = "Bisorted", row.norm = T, cmap.type = 3,
    char.rescale = char.res, legend = T
  )

  dev.off()
}

MSIG.Define.Dataset.from.Table <- function(
    input.gct,
    table.txt,
    output.gct,
    output.cls,
    prefix_entries = F) {
  # Read input dataset

  library(RColorBrewer)

  dataset1 <- MSIG.Gct2Frame(filename = input.gct)
  m <- data.matrix(dataset1$ds)
  gene.names <- dataset1$row.names
  gene.decs <- dataset1$descs
  sample.names.gct <- dataset1$names
  Ns <- length(sample.names.gct)

  # Read Table


  tab <- read.delim(table.txt, header = T, row.names = 1, sep = "\t", skip = 0, blank.lines.skip = T, comment.char = "", as.is = T)
  sample.names.tab <- row.names(tab)
  phen.names <- names(tab)
  overlap <- intersect(sample.names.tab, sample.names.gct)
  print("sample names GCT")
  print(sample.names.gct)
  print("sample names TAB")
  print(sample.names.tab)

  locs.gct <- match(overlap, sample.names.gct)
  print(match(sample.names.tab, sample.names.gct))
  print(match(sample.names.gct, sample.names.tab))
  locs.tab <- match(overlap, sample.names.tab)
  print(locs.tab)
  print(c("GCT matching set (", length(locs.gct), " samples):", sample.names.gct[locs.gct]))
  print(c("TAB matching set (", length(overlap), " samples):", sample.names.tab[locs.tab]))
  print(c("overlap set (", length(overlap), " samples):", overlap))

  m2 <- m[, locs.gct]
  sample.names.gct <- sample.names.gct[locs.gct]
  sample.names.tab <- sample.names.tab[locs.tab]
  cls.table <- t(tab[locs.tab, ])

  if (prefix_entries == TRUE) {
    for (i in 1:length(cls.table[, 1])) {
      #        cls.table[i,] <- paste(row.names(cls.table)[i], cls.table[i,], sep=".")
      cls.table[i, ] <- paste(colnames(tab)[i], tab[, i], sep = ".")
    }
  }

  if (!is.null(output.gct)) {
    V <- data.frame(m2)
    names(V) <- sample.names.gct
    row.names(V) <- gene.names
    write.gct(gct.data.frame = V, descs = gene.decs, filename = output.gct)
  }

  class.phen <- unique(cls.table)
  n <- length(class.phen)
  l <- length(cls.table[1, ])

  col.list <- c(
    brewer.pal(n = 7, name = "Set2"), brewer.pal(n = 7, name = "Dark2"), brewer.pal(n = 7, name = "Set1"),
    brewer.pal(n = 8, name = "Accent"), brewer.pal(n = 10, name = "Spectral"), brewer.pal(n = 8, name = "Set3"),
    brewer.pal(n = 8, name = "BuGn"),
    brewer.pal(n = 7, name = "Set2"), brewer.pal(n = 7, name = "Dark2"), brewer.pal(n = 7, name = "Set1"),
    brewer.pal(n = 8, name = "Accent"), brewer.pal(n = 10, name = "Spectral"), brewer.pal(n = 8, name = "Set3"),
    brewer.pal(n = 8, name = "BuGn"),
    brewer.pal(n = 7, name = "Set2"), brewer.pal(n = 7, name = "Dark2"), brewer.pal(n = 7, name = "Set1"),
    brewer.pal(n = 8, name = "Accent"), brewer.pal(n = 10, name = "Spectral"), brewer.pal(n = 8, name = "Set3"),
    brewer.pal(n = 8, name = "BuGn")
  )
  num <- 0
  class.order.list <- NULL
  for (i in 1:length(cls.table[, 1])) {
    num <- num + length(unique(cls.table[i, ]))
    class.order.list <- c(class.order.list, unique(cls.table[i, ]))
  }

  phen.names.string <- paste("phen.names:", paste(phen.names, collapse = " "), sep = " ")
  sig.col <- col.list[1:num]
  col.phen.string <- paste("col.phen:", paste(sig.col, collapse = " "), sep = " ")
  cat(paste(l, num, length(cls.table[, 1]), phen.names.string, col.phen.string, sep = " "), "\n",
    file = output.cls, append = FALSE, sep = ""
  )
  cat("# ", paste(class.order.list, collapse = " "), "\n", file = output.cls, append = TRUE, sep = "")
  for (i in 1:length(cls.table[, 1])) {
    cat(paste(cls.table[i, ], collapse = " "), "\n", file = output.cls, append = TRUE, sep = "")
  }
}

OPAM.apply.model.2 <- function(
    input.ds,
    input.cls = NA,
    models.dir,
    models = "ALL",
    raw.score.outfile,
    norm.score.outfile,
    model.score.outfile,
    prob.outfile,
    gmt.file = NULL,
    graphics.off = F) { #----------------------------------------------------------------------------------------

  # Load libraries
  erf <- function(x) 2 * pnorm(x * sqrt(2)) - 1
  library(gtools)
  library(verification)
  library(ROCR)
  library(MASS)
  library(RColorBrewer)
  library(MCMCpack)

  # Read test dataset
  dataset <- MSIG.Gct2Frame(filename = input.ds) # Read gene expression dataset (GCT format)
  m <- data.matrix(dataset$ds)
  gene.names <- dataset$row.names
  gene.descs <- dataset$descs
  sample.names <- dataset$names
  Ns <- length(m[1, ])
  Ng <- length(m[, 1])
  temp <- strsplit(input.ds, split = "/") # Extract test file name
  s <- length(temp[[1]])
  test.file.name <- temp[[1]][s]
  temp <- strsplit(test.file.name, split = ".gct")
  test.file.prefix <- temp[[1]][1]

  # Test set color map
  c.test <- c(brewer.pal(n = 7, name = "Set2"), brewer.pal(n = 7, name = "Dark2"), brewer.pal(n = 7, name = "Set1"))

  if (!is.na(input.cls)) { # Read phenotype if CLS file is provided
    CLS <- MSIG.ReadClsFile(file = input.cls) # Read phenotype file (CLS format)
    class.labels <- CLS$class.v
    class.phen <- CLS$phen
    class.list <- CLS$class.list
  } else {
    class.labels <- rep(1, Ns)
    class.phen <- "UNDEFINED_PHEN"
    class.list <- rep("U", Ns)
  }

  # Loop over models

  if (models[[1]] == "ALL") { # use all models in directory
    file.list <- list.files(models.dir)
    models <- file.list[regexpr(pattern = ".mod", file.list) > 1]
    for (k.model in 1:length(models)) {
      temp <- strsplit(models[k.model], ".mod")
      models[k.model] <- temp[[1]]
    }
    models <- unique(models)
  }

  n.models <- length(models)
  score.matrix <- matrix(0, nrow = n.models, ncol = Ns)
  norm.score.matrix <- matrix(0, nrow = n.models, ncol = Ns)
  model.score.matrix <- matrix(0, nrow = n.models, ncol = Ns)
  probability.matrix <- matrix(0, nrow = n.models, ncol = Ns)
  models.descs <- NULL

  for (model.i in 1:n.models) {
    print(paste(model.i, "File:", models[model.i]))
    # Read parameters from model file
    m.file <- paste(models.dir, models[model.i], ".mod", sep = "")
    con <- file(m.file, "r")
    file.content <- readLines(con, n = -1)
    close(con)
    gc()
    len <- length(file.content)
    for (i in 1:len) {
      temp <- unlist(strsplit(file.content[[i]], "\t"))
      len.param <- length(temp)
      if (len.param == 2) {
        param.vals <- temp[2]
      } else {
        param.vals <- paste(noquote(temp[2:len.param]), collapse = ",")
        param.vals <- paste("c(", param.vals, ")", sep = "")
      }
      assignment.string <- paste(noquote(temp[1]), " <- ", param.vals, sep = "")
      #         print(c("executing assigment:", assignment.string))
      eval(parse(text = assignment.string))
    }

    print(paste("Model:", model.i, " model name:", model.name))


    # Set parameters
    if (!exists("random.seed")) random.seed <- 12345
    set.seed(random.seed)

    # Sample normalization

    if (!exists("sample.norm.type")) sample.norm.type <- "rank"

    if (sample.norm.type == "rank") {
      for (j in 1:Ns) { # column rank normalization
        m[, j] <- rank(m[, j], ties.method = "average")
      }
      m <- 10000 * m / Ng
    } else if (sample.norm.type == "log.rank") {
      for (j in 1:Ns) { # column rank normalization
        m[, j] <- rank(m[, j], ties.method = "average")
      }
      m <- log(10000 * m / Ng + exp(1))
    } else if (sample.norm.type == "log") {
      m[m < 1] <- 1
      m <- log(m + exp(1))
    }

    # Control signature normalization
    if (!exists("msig.cntrl.genes")) msig.cntrl.genes <- NA
    if (!is.na(msig.cntrl.genes)) {
      gene.names.int <- intersect(msig.cntrl.genes, gene.names)
      locs <- match(gene.names.int, gene.names, nomatch = 0)
      msig.cntrl <- m[locs, ]
      msig.cntrl.genes <- gene.names[locs]
      msig.cntrl.descs <- gene.descs[locs]
      msig.cntrl.size <- length(locs)
      if (msig.cntrl.size < 1) {
        msig.cntrl.center <- rep(1, Ns)
      } else if (msig.cntrl.size == 1) {
        msig.cntrl.center <- msig.cntrl
      } else if (msig.cntrl.size > 1) msig.cntrl.center <- apply(msig.cntrl, MARGIN = 2, FUN = mean)
      for (i in 1:Ng) {
        m[i, ] <- m[i, ] / msig.cntrl.center
      }
    }

    # Obtain UP & DN signatures

    if (exists("msig.up.genes3")) msig.up.genes <- msig.up.genes3
    gene.names.int <- intersect(msig.up.genes, gene.names)
    if (length(gene.names.int) < 2) {
      score.matrix[model.i, ] <- norm.score.matrix[model.i, ] <- model.score.matrix[model.i, ] <- probability.matrix[model.i, ] <- rep(0, Ns)
      models.descs <- c(models.descs, model.description)
      rm(
        model.creation.date, input.ds, input.cls, input2.ds, input2.cls, target.class,
        target.class2, model.name, model.description, sample.norm.type, marker.disc,
        top.markers.up, top.markers.dn, top.markers.up2, top.markers.dn2, statistic, weight,
        random.seed, nperm, link.function, c1, msig.cntrl.genes, msig.up.genes, msig.dn.genes,
        msig.up.genes2, msig.dn.genes2, msig.up.genes3, msig.dn.genes3, beta.0, beta.1, score,
        score.up, score.dn
      )
      next
    }

    if (!is.null(gmt.file)) { # save genes in a gmt file
      m.name <- paste(model.name, "_UP", sep = "")
      genes.string <- paste(msig.up.genes, sep = "\t", collapse = "\t")
      output.line <- paste(m.name, model.description, genes.string, sep = "\t", collapse = "")
      if (model.i == 1) {
        write(noquote(output.line), file = gmt.file, append = F, ncolumns = length(msig.up.genes) + 2)
      } else {
        write(noquote(output.line), file = gmt.file, append = T, ncolumns = length(msig.up.genes) + 2)
      }
      if (exists("msig.dn.genes")) {
        m.name <- paste(model.name, "_DN", sep = "")
        genes.string <- paste(msig.dn.genes, sep = "\t", collapse = "\t")
        output.line <- paste(m.name, model.description, genes.string, sep = "\t", collapse = "")
        write(noquote(output.line), file = gmt.file, append = T, ncolumns = length(msig.up.genes) + 2)
      }
    }


    locs <- match(gene.names.int, gene.names, nomatch = 0)
    msig.up.test <- m[locs, ]
    msig.up.genes.test <- gene.names[locs]
    msig.up.descs.test <- gene.descs[locs]
    msig.up.size.test <- length(locs)
    if (graphics.off == F) {
      x11(height = 9, width = 12)
      MSIG.HeatMapPlot.3(
        V = msig.up.test, row.names = msig.up.genes.test, col.labels = class.labels,
        col.classes = class.phen, phen.cmap = c.test, col.names = sample.names,
        main = paste(model.name, " UP signature test"), xlab = " ", ylab = " ", sub = " ",
        row.norm = T, cmap.type = 4, char.rescale = 1
      )
    }
    # Project test dataset
    OPAM <- OPAM.Projection(m, gene.names, Ns, Ng, weight, statistic, msig.up.genes.test, nperm = nperm)
    score.up <- OPAM$ES.vector

    if (exists("msig.dn.genes3")) msig.dn.genes <- msig.dn.genes3
    if (exists("msig.dn.genes")) {
      gene.names.int <- intersect(msig.dn.genes, gene.names)
      if (length(gene.names.int) < 2) {
        score.matrix[model.i, ] <- norm.score.matrix[model.i, ] <- model.score.matrix[model.i, ] <- probability.matrix[model.i, ] <- rep(0, Ns)
        models.descs <- c(models.descs, model.description)
        rm(
          model.creation.date, input.ds, input.cls, input2.ds, input2.cls, target.class,
          target.class2, model.name, model.description, sample.norm.type, marker.disc,
          top.markers.up, top.markers.dn, top.markers.up2, top.markers.dn2, statistic, weight,
          random.seed, nperm, link.function, c1, msig.cntrl.genes, msig.up.genes, msig.dn.genes,
          msig.up.genes2, msig.dn.genes2, msig.up.genes3, msig.dn.genes3, beta.0, beta.1, score,
          score.up, score.dn
        )
        next
      }
      locs <- match(gene.names.int, gene.names, nomatch = 0)
      msig.dn.test <- m[locs, ]
      msig.dn.genes.test <- gene.names[locs]
      msig.dn.descs.test <- gene.descs[locs]
      msig.dn.size.test <- length(locs)
      if (graphics.off == F) {
        x11(height = 9, width = 12)
        MSIG.HeatMapPlot.3(
          V = msig.dn.test, row.names = msig.dn.genes.test, col.labels = class.labels,
          col.classes = class.phen, phen.cmap = c.test, col.names = sample.names,
          main = paste(model.name, " DN signature test"), xlab = " ", ylab = " ", sub = " ",
          row.norm = T, cmap.type = 4, char.rescale = 1
        )
      }
      OPAM <- OPAM.Projection(m, gene.names, Ns, Ng, weight, statistic, msig.dn.genes.test, nperm = nperm)
      score.dn <- OPAM$ES.vector
    }

    if (!is.na(msig.cntrl.genes)) {
      if (graphics.off == F) {
        x11(height = 9, width = 12)
        MSIG.HeatMapPlot.3(
          V = msig.cntrl.test, row.names = msig.cntrl.genes.test, col.labels = class.labels,
          col.classes = class.phen, phen.cmap = c.test, col.names = sample.names,
          main = paste(model.name, " CNTRL signature"), xlab = " ", ylab = " ", sub = " ", row.norm = T,
          cmap.type = 4, char.rescale = 1
        )
      }
    }

    if (exists("msig.dn.genes")) {
      score <- score.up - score.dn
    } else {
      score <- score.up
    }

    if (graphics.off == F) {
      x11(width = 14, height = 9)
      if (exists("msig.dn.genes")) {
        nf <- layout(matrix(c(1, 2, 3, 0, 4, 0), 3, 2, byrow = T), widths = c(5, 1), heights = c(1, 1, 1), respect = FALSE)
      } else {
        nf <- layout(matrix(c(1, 2), 1, 2, byrow = T), widths = c(5, 1), heights = 1, respect = FALSE)
      }

      par(mar = c(2, 4, 2, 4))
      barplot(score.up,
        main = paste(model.name, " OPAM Score UP (test)", sep = ""), font.axis = 1.25, cex.lab = 1.5,
        cex.axis = 1.25, cex.names = 1.25, width = 1, space = 0, col = c.test[class.labels]
      )
      leg.txt <- class.phen
      p.vec <- rep(22, length(leg.txt))
      par(mar = c(0, 0, 0, 0))
      plot(c(0, 0), c(1, 1), xlim = c(0, 1), ylim = c(0, 1), axes = F, type = "n", xlab = "", ylab = "")
      legend(
        x = 0, y = 0.8, legend = leg.txt, bty = "n", xjust = 0, yjust = 1, pch = p.vec, pt.bg = c.test, col = "black",
        cex = 1.25, pt.cex = 2.5
      )

      if (exists("msig.dn.genes")) {
        par(mar = c(2, 4, 2, 4))
        barplot(score.dn,
          main = paste(model.name, " OPAM Score DOWN (test)", sep = ""), font.axis = 1.25, cex.lab = 1.5,
          cex.axis = 1.25, cex.names = 1.25, width = 1, space = 0, col = c.test[class.labels]
        )

        par(mar = c(2, 4, 2, 4))
        barplot(score,
          main = paste(model.name, " OPAM Total Score (test)", sep = ""), font.axis = 1.25, cex.lab = 1.5,
          cex.axis = 1.25, cex.names = 1.25, width = 1, space = 0, col = c.test[class.labels]
        )
      }
    }

    if (!exists("beta.0")) beta.0 <- 0
    if (!exists("beta.1")) beta.1 <- 1
    if (!exists("link.function")) link.function <- "logit"

    # Apply MCMC logit or probit model to test dataset

    model.formula <- "beta.0 + beta.1 * score[i]"
    model.formula
    prob.i <- matrix(0, nrow = Ns, ncol = 3)
    model.score <- vector(length = Ns, mode = "numeric")
    for (i in 1:Ns) {
      model.score[i] <- eval(parse(text = model.formula))
      if (link.function == "logit") {
        p.vec <- paste("inv.logit(x=", model.formula, ")", sep = "")
      } else if (link.function == "probit") {
        p.vec <- paste("(erf(", model.formula, ") + 1)/2", sep = "")
      } else {
        stop("Unknown link function")
      }
      val <- eval(parse(text = p.vec))
      prob.i[i, 1] <- quantile(val, probs = 0.5)
      prob.i[i, 2] <- quantile(val, probs = 0.05)
      prob.i[i, 3] <- quantile(val, probs = 0.95)
    }
    probability <- prob.i[, 1]
    xmin <- min(model.score)
    xmax <- max(model.score)
    range.x <- xmax - xmin
    n.points <- 1000
    prob.m <- matrix(0, nrow = n.points, ncol = 3)
    x.m <- vector(length = n.points, mode = "numeric")
    for (k in 1:n.points) {
      x.m[k] <- xmin + k * (range.x / n.points)
      if (link.function == "logit") {
        p.vec <- paste("inv.logit(x=", x.m[k], ")", sep = "")
      } else if (link.function == "probit") {
        p.vec <- paste("(erf(", x.m[k], ") + 1)/2", sep = "")
      } else {
        stop("Unknown link function")
      }
      val <- eval(parse(text = p.vec))
      prob.m[k, 1] <- quantile(val, probs = 0.5)
      prob.m[k, 2] <- quantile(val, probs = 0.05)
      prob.m[k, 3] <- quantile(val, probs = 0.95)
    }
    istar <- which.min(abs(0.5 - prob.m[, 1]))
    istar <- xmin + istar * (range.x / 1000)
    x.index <- order(model.score, decreasing = F)
    x.order <- model.score[x.index]
    prob.i.order <- prob.i[x.index, ]
    target.var.order <- c.test[class.labels[x.index]]
    class.labels.order <- class.labels[x.index]


    # Plot bar graph of z-scores
    boundary <- istar
    pred.class <- ifelse(prob.i.order[, 1] >= 0.5, 2, 1)
    z.range <- range(x.order)
    norm.score <- (score - min(score)) / (max(score) - min(score))
    if (graphics.off == F) {
      x11(height = 7, width = 9.5)
      nf <- layout(matrix(c(1, 2), 1, 2, byrow = T), widths = c(3.75, 1), heights = 1, respect = FALSE)
      plot(x.order, prob.i.order[, 1], sub = model.name, pch = 20, col = 0, cex = 2, xlab = "Activation Index", ylab = "Probability")
      points(x.m, prob.m[, 1], type = "l", lwd = 2, col = 1, lty = 1, cex = 1)
      points(x.m, prob.m[, 2], type = "l", col = 4, lty = 1, cex = 1)
      points(x.m, prob.m[, 3], type = "l", col = 4, lty = 1, cex = 1)
      arrows(x.order, prob.i.order[, 2], x.order, prob.i.order[, 3], col = 4, angle = 90, code = 3, length = 0.0)
      range.x <- range(x.order)
      points(range.x, c(0.5, 0.5), type = "l", lty = 3, col = 1, lwd = 2)
      points(range.x, c(-.15, -0.15), type = "l", lty = 1, col = 1, lwd = 2)
      points(c(istar, istar), c(-0.07, 1.07), type = "l", lty = 3, col = 1, lwd = 2)
      points(x.order, prob.i.order[, 1], pch = 21, bg = target.var.order, col = 1, cex = 2)
      points(x.order, rep(-0.15, length(x.order)), pch = 21, bg = target.var.order, col = 1, cex = 2)
      leg.txt <- class.phen
      p.vec <- rep(22, length(leg.txt))
      c.vec <- c.test
      par(mar = c(0, 0, 0, 0))
      plot(c(0, 0), c(1, 1), xlim = c(0, 1), ylim = c(0, 1), axes = F, type = "n", xlab = "", ylab = "")
      legend(
        x = 0, y = 0.8, legend = leg.txt, bty = "n", xjust = 0, yjust = 1, pch = p.vec, pt.bg = c.vec,
        col = "black", cex = 1.25, pt.cex = 2.5
      )

      x11(width = 14, height = 9)
      nf <- layout(matrix(c(1, 2, 3, 4, 5, 6), 3, 2, byrow = T), widths = c(5, 1), heights = c(1, 1, 1), respect = FALSE)
      par(mar = c(4, 7, 4, 7))
      MSIG.Score.Plot(
        z = score, main = paste(model.name, " Model Score (test)"), phen.cmap = c.test,
        char.rescale = 1, col.classes = class.phen, col.labels = class.labels,
        xlab = "Samples", ylab = "Score", create.window = F, create.legend = T
      )

      par(mar = c(4, 7, 4, 7))

      MSIG.Score.Plot(
        z = norm.score, main = paste(model.name, " Normalized Model Score (test)"), phen.cmap = c.test,
        char.rescale = 1, col.classes = class.phen, col.labels = class.labels,
        xlab = "Samples", ylab = "Normalized Score", create.window = F, create.legend = T
      )

      par(mar = c(4, 7, 4, 7))
      MSIG.Score.Plot(
        z = prob.i[, 1], main = paste(model.name, " Probabiliy (test)"), phen.cmap = c.test, char.rescale = 1,
        col.classes = class.phen, col.labels = class.labels,
        xlab = "Samples", ylab = "Probability", create.window = F, create.legend = T
      )

      x11(width = 14, height = 9)
      MSIG.HeatMapPlot.6(
        V = rbind(score, norm.score, model.score, probability), row.names = c(
          "raw.score", "norm.score",
          "model.score", "probability"
        ),
        row.names2 = c(model.name, model.name, model.name, model.name),
        col.labels = class.labels, col.labels2 = class.labels,
        col.classes = class.phen, phen.cmap = c.test, phen.names = model.name,
        col.names = sample.names, main = model.name,
        xlab = "  ", ylab = "  ", sub = "   ", row.norm = T, cmap.type = 3,
        char.rescale = 1, legend = T
      )
    }
    score.matrix[model.i, ] <- score
    norm.score.matrix[model.i, ] <- norm.score
    model.score.matrix[model.i, ] <- model.score
    probability.matrix[model.i, ] <- probability
    models.descs <- c(models.descs, model.description)
    rm(
      model.creation.date, input.ds, input.cls, input2.ds, input2.cls, target.class,
      target.class2, model.name, model.description, sample.norm.type, marker.disc,
      top.markers.up, top.markers.dn, top.markers.up2, top.markers.dn2, statistic, weight,
      random.seed, nperm, link.function, c1, msig.cntrl.genes, msig.up.genes, msig.dn.genes,
      msig.up.genes2, msig.dn.genes2, msig.up.genes3, msig.dn.genes3, beta.0, beta.1, score,
      score.up, score.dn
    )

    if (graphics.off == F) {
      if (model.i %% 5 == 0) {
        graphics.off()
        gc()
      }
    }
  } # end of loop over models

  # Plot projections for all models

  if (FALSE) {
    '
   x11(width = 14, height = 9)
   MSIG.HeatMapPlot.6(V = score.matrix, row.names = models, row.names2 = models.descs, col.labels = class.labels,
                      col.labels2 = class.labels, col.classes = class.phen, phen.cmap = c.test,
                      col.names = sample.names, main = paste(test.file.name, " - Raw Scores"), xlab=" ", ylab=" ",
                      sub = "Raw Scores", row.norm = T, cmap.type = 3, char.rescale = 1, legend=T)

   x11(width = 14, height = 9)
   MSIG.HeatMapPlot.6(V = norm.score.matrix, row.names = models, row.names2 = models.descs, col.labels = class.labels,
                      col.labels2 = class.labels, col.classes = class.phen, phen.cmap = c.test,
                      col.names = sample.names, main = paste(test.file.name, " - Norm Scores"), xlab=" ", ylab=" ",
                      sub = "Norm Scores", row.norm = T, cmap.type = 3, char.rescale = 1, legend=T)

   x11(width = 14, height = 9)
   MSIG.HeatMapPlot.6(V = model.score.matrix, row.names = models, row.names2 = models.descs, col.labels = class.labels,
                      col.labels2 = class.labels, col.classes = class.phen, phen.cmap = c.test,
                      col.names = sample.names, main = paste(test.file.name, " - Model Scores"), xlab=" ", ylab=" ",
                      sub = "Model Scores", row.norm = T, cmap.type = 3, char.rescale = 1, legend=T)

   x11(width = 14, height = 9)
   MSIG.HeatMapPlot.6(V = probability.matrix, row.names = models, row.names2 = models.descs, col.labels = class.labels,
                      col.labels2 = class.labels, col.classes = class.phen, phen.cmap = c.test,
                      col.names = sample.names, main = paste(test.file.name, " - Probabilities"), xlab=" ", ylab=" ",
                      sub = "Probabilities", row.norm = T, cmap.type = 3, char.rescale = 1, legend=T)
   '
  }

  # Save projections in files

  V.GCT <- data.frame(score.matrix)
  names(V.GCT) <- sample.names
  row.names(V.GCT) <- models
  write.gct(gct.data.frame = V.GCT, descs = models.descs, filename = raw.score.outfile)

  V.GCT <- data.frame(norm.score.matrix)
  names(V.GCT) <- sample.names
  row.names(V.GCT) <- models
  write.gct(gct.data.frame = V.GCT, descs = models.descs, filename = norm.score.outfile)

  V.GCT <- data.frame(model.score.matrix)
  names(V.GCT) <- sample.names
  row.names(V.GCT) <- models
  write.gct(gct.data.frame = V.GCT, descs = models.descs, filename = model.score.outfile)

  V.GCT <- data.frame(probability.matrix)
  names(V.GCT) <- sample.names
  row.names(V.GCT) <- models
  write.gct(gct.data.frame = V.GCT, descs = models.descs, filename = prob.outfile)
} # end of function

OPAM.create.MSigDB.models <- function(
    gs.db,
    models.dir,
    gs.size.threshold.min = 10,
    gs.size.threshold.max = 500,
    source = "MSigDB") {
  temp <- readLines(gs.db)
  max.Ng <- length(temp)
  temp.size.G <- vector(length = max.Ng, mode = "numeric")
  for (i in 1:max.Ng) {
    temp.size.G[i] <- length(unlist(strsplit(temp[[i]], "\t"))) - 2
  }
  max.size.G <- max(temp.size.G)
  gs <- matrix(rep("null", max.Ng * max.size.G), nrow = max.Ng, ncol = max.size.G)
  temp.names <- vector(length = max.Ng, mode = "character")
  temp.desc <- vector(length = max.Ng, mode = "character")
  gs.count <- 1
  for (i in 1:max.Ng) {
    gene.set.size <- length(unlist(strsplit(temp[[i]], "\t"))) - 2
    gs.line <- noquote(unlist(strsplit(temp[[i]], "\t")))
    gene.set.name <- gs.line[1]
    gene.set.desc <- gs.line[2]
    gene.set.tags <- vector(length = gene.set.size, mode = "character")
    for (j in 1:gene.set.size) {
      gene.set.tags[j] <- gs.line[j + 2]
    }
    set.size <- length(gene.set.tags)
    if ((set.size < gs.size.threshold.min) || (set.size > gs.size.threshold.max)) next
    temp.size.G[gs.count] <- set.size
    gs[gs.count, ] <- c(gene.set.tags, rep("null", max.size.G - temp.size.G[gs.count]))
    temp.names[gs.count] <- gene.set.name
    temp.desc[gs.count] <- gene.set.desc
    gs.count <- gs.count + 1
  }
  Ng <- gs.count - 1
  gs.names <- vector(length = Ng, mode = "character")
  gs.desc <- vector(length = Ng, mode = "character")
  size.G <- vector(length = Ng, mode = "numeric")
  gs.names <- temp.names[1:Ng]
  gs.desc <- temp.desc[1:Ng]
  size.G <- temp.size.G[1:Ng]

  print(c("Number of Gene Sets:", Ng))
  print(c("Original number of Gene Sets:", max.Ng))
  print(c("Maximum gene set size:", max.size.G))

  for (i in 1:Ng) {
    print(paste("Creating model for gene set:", i, gs.names[i], sep = " "))
    gene.set <- gs[i, gs[i, ] != "null"]
    gene.set.string <- paste("c('", paste(gene.set, collapse = "','"), "')\n", sep = "")
    gene.set.name <- gs.names[i]
    gene.set.desc <- gs.desc[i]
    m.file <- paste(models.dir, gene.set.name, ".mod", sep = "")
    cat("model.creation.date", paste("'", date(), "'\n", sep = ""), file = m.file, append = FALSE, sep = "\t")
    cat("model.name", paste("'", gene.set.name, "'\n", sep = ""), file = m.file, append = TRUE, sep = "\t")
    cat("model.description", paste("'", source, "'\n", sep = ""), file = m.file, append = TRUE, sep = "\t")
    cat("sample.norm.type", "'rank'\n", file = m.file, append = TRUE, sep = "\t")
    cat("statistic", "'area.under.RES'\n", file = m.file, append = TRUE, sep = "\t")
    cat("weight", paste(0.25, "\n", sep = ""), file = m.file, append = TRUE, sep = "\t")
    cat("random.seed", paste(12345, "\n", sep = ""), file = m.file, append = TRUE, sep = "\t")
    cat("nperm", paste(0, "\n", sep = ""), file = m.file, append = TRUE, sep = "\t")
    cat("link.function", "'logit'\n", file = m.file, append = TRUE, sep = "\t")
    cat("c1", "c('black','lightgrey')\n", file = m.file, append = TRUE, sep = "\t")
    cat("msig.up.genes", gene.set.string, file = m.file, append = TRUE, sep = "\t")
  }
}

OPAM.compute.models.overlap <- function(
    model,
    model.set = NA,
    models.dir,
    thres = 0.05,
    size.lim = 500,
    results.file,
    produce.overlap.models = F) { #----------------------------------------------------------------------------------------
  file.list <- list.files(models.dir)
  files <- file.list[regexpr(pattern = ".mod", file.list) > 1]

  if (!is.na(model.set[1])) {
    model.set <- c(model, model.set)
    model.set <- paste(model.set, ".mod", sep = "")
    files <- intersect(files, model.set)
  }
  print(c("files=", files))
  max.models <- length(files)
  model.sizes <- NULL
  model.genes <- NULL
  model.names <- NULL
  num.models <- 0

  for (model.i in 1:max.models) {
    # Read parameters from model file

    m.file <- paste(models.dir, files[model.i], sep = "")
    #      print(paste(model.i, "Reading file:", m.file))
    file.content <- readLines(m.file, n = -1)
    len <- length(file.content)
    for (i in 1:len) {
      temp <- unlist(strsplit(file.content[[i]], "\t"))
      len.param <- length(temp)
      if (len.param == 2) {
        param.vals <- temp[2]
      } else {
        param.vals <- paste(noquote(temp[2:len.param]), collapse = ",")
        param.vals <- paste("c(", param.vals, ")", sep = "")
      }
      assignment.string <- paste(noquote(temp[1]), " <- ", param.vals, sep = "")
      #         print(c("executing assigment:", assignment.string))
      eval(parse(text = assignment.string))
    }
    if (exists("msig.up.genes3")) msig.up.genes <- msig.up.genes3
    m.size <- length(msig.up.genes)
    if (m.size <= size.lim) {
      model.names <- c(model.names, model.name)
      model.sizes <- c(model.sizes, m.size)
      model.genes <- rbind(model.genes, c(msig.up.genes, rep(NA, size.lim - m.size)))
      num.models <- num.models + 1
    }
    rm(model.name, m.size, msig.up.genes, msig.up.genes3)
  } # loop over models

  print(c("model:", model, " model.names:", model.names))
  loc <- match(model, model.names)
  if (is.na(loc)) stop(paste("model:", model, " is not in models directory: ", models.dir))
  m.size <- model.sizes[loc]
  m.genes <- model.genes[loc, 1:m.size]
  model.overlap.size <- vector(length = num.models, mode = "numeric")
  model.overlap.genes <- matrix(NA, nrow = num.models, ncol = size.lim)
  model.overlap.signif <- vector(length = num.models, mode = "numeric")

  for (i in 1:num.models) {
    overlap <- intersect(m.genes, model.genes[i, 1:model.sizes[i]])
    if (length(overlap) == 0) {
      model.overlap.size[i] <- 0
    } else {
      model.overlap.size[i] <- length(overlap)
      model.overlap.genes[i, 1:length(overlap)] <- overlap
      model.overlap.signif[i] <- (length(overlap) / m.size) * (length(overlap) / model.sizes[i])
    }
  }

  save.flag <- ifelse(model.overlap.signif >= thres, 1, 0)
  total.save <- sum(save.flag)
  results.tab <- cbind(
    model.names, model.sizes, rep(m.size, num.models), model.overlap.size,
    model.overlap.signif, save.flag
  )
  ind <- order(model.overlap.signif, decreasing = T)
  results.tab <- results.tab[ind, ]
  colnames(results.tab) <- c("Name", "Size", "Model Size", "Overlap", "Score", "Save Model")
  results.file <- paste(models.dir, "OVER_", model, ".txt", sep = "")
  write.table(results.tab, file = results.file, quote = F, row.names = F, sep = "\t")
  print(noquote(results.tab[1:(min(50, num.models)), ]))

  if (produce.overlap.models == T) {
    for (k in 1:total.save) {
      if (model == model.names[ind[k]]) next
      print(paste("Creating model for model:", model.names[ind[k]], sep = " "))
      gene.set <- model.overlap.genes[ind[k], 1:model.overlap.size[ind[k]]]
      gene.set.string <- paste("c('", paste(gene.set, collapse = "','"), "')\n", sep = "")
      overlap.model.name <- paste("OVER_", model, "_", model.names[ind[k]], sep = "")
      print(paste("Saving model: ", overlap.model.name))
      print(gene.set.string)
      print("------------")
      m.file <- paste(models.dir, overlap.model.name, ".mod", sep = "")
      cat("model.creation.date", paste("'", date(), "'\n", sep = ""), file = m.file, append = FALSE, sep = "\t")
      cat("model.name", paste("'", overlap.model.name, "'\n", sep = ""), file = m.file, append = TRUE, sep = "\t")
      cat("model.description", paste("'Overlap Model'\n", sep = ""), file = m.file, append = TRUE, sep = "\t")
      cat("sample.norm.type", "'rank'\n", file = m.file, append = TRUE, sep = "\t")
      cat("statistic", "'area.under.RES'\n", file = m.file, append = TRUE, sep = "\t")
      cat("weight", paste(0.25, "\n", sep = ""), file = m.file, append = TRUE, sep = "\t")
      cat("random.seed", paste(12345, "\n", sep = ""), file = m.file, append = TRUE, sep = "\t")
      cat("nperm", paste(0, "\n", sep = ""), file = m.file, append = TRUE, sep = "\t")
      cat("link.function", "'logit'\n", file = m.file, append = TRUE, sep = "\t")
      cat("c1", "c('black','lightgrey')\n", file = m.file, append = TRUE, sep = "\t")
      cat("msig.up.genes", gene.set.string, file = m.file, append = TRUE, sep = "\t")
    }
  }
} # end of function

make.contin.table <- function(pred, actual) {
  # computes a complete all-class contingency table (even when there are missing classes in predicted entries)

  tab <- table(pred, actual)
  pred.classes <- unique(pred)
  actual.classes <- unique(actual)
  missing.classes <- setdiff(actual.classes, pred.classes)
  tab.rows <- row.names(tab)
  tab <- rbind(tab, matrix(0, nrow = length(missing.classes), ncol = length(tab[1, ])))
  row.names(tab) <- c(tab.rows, missing.classes)
  row.index <- order(row.names(tab), decreasing = F)
  col.index <- order(colnames(tab), decreasing = F)
  tab2 <- tab[row.index, col.index]
  return(tab2)
}


MISG.Res2Gct <- function(
    res.file,
    gct.file) {
  dataset <- MSIG.Res2Frame(filename = res.file) # read RES file
  A <- dataset$ds
  row.names(A) <- dataset$row.names
  colnames(A) <- dataset$names
  descs <- dataset$descs
  write.gct(gct.data.frame = A, descs = descs, filename = gct.file)
}


OPAM.match.projection.to.phenotypes <- function(
    input.ds,
    input.cls,
    results.dir,
    normalize.score = T,
    normalization.type = "zero.one", # "zero.one", "z.score" or "r.z.score"
    markers.num = 5,
    user.colors = NA,
    markers.metric = "ROC", # "ROC" or "T.TEST"
    markers.file = NULL,
    sort.phenotypes = T,
    sort.decreasing = T, # T = decreasing, F = increasing
    sort.expression = T,
    sort.decreasing.genes = T,
    legend = T,
    char.res = 1,
    only.up = F,
    cmap.type = 3,
    row.norm = T) {
  library(gtools)
  library(verification)
  library(ROCR)
  library(MASS)
  library(RColorBrewer)
  library(heatmap.plus)

  dataset <- MSIG.Gct2Frame(filename = input.ds) # Read gene expression dataset (GCT format)
  m <- data.matrix(dataset$ds)
  model.names <- dataset$row.names
  model.descs <- dataset$descs
  Ns <- length(m[1, ])

  for (i in 1:length(m[, 1])) {
    if (sd(m[i, ]) == 0) {
      val <- m[i, 1]
      m[i, ] <- m[i, ] + runif(n = Ns, min = val - 0.001, max = val + 0.001) # add small noise to flat profiles
    }
  }
  dim(m)
  sample.names <- dataset$names

  n.models <- length(m[, 1])
  temp <- strsplit(input.ds, split = "/") # Extract test file name
  s <- length(temp[[1]])
  test.file.name <- temp[[1]][s]
  temp <- strsplit(test.file.name, split = ".gct")
  test.file.prefix <- temp[[1]][1]
  #   char.res <-  0.013 * n.models + 0.65

  # normalize scores

  if (normalize.score == T) {
    if (normalization.type == "zero.one") {
      for (i in 1:n.models) {
        m[i, ] <- (m[i, ] - min(m[i, ])) / (max(m[i, ]) - min(m[i, ]))
      }
    } else if (normalization.type == "z.score") {
      for (i in 1:n.models) {
        m[i, ] <- (m[i, ] - mean(m[i, ])) / sd(m[i, ])
      }
    } else if (normalization.type == "r.z.score") {
      for (i in 1:n.models) {
        m[i, ] <- (m[i, ] - median(m[i, ])) / mad(m[i, ])
      }
    }
  }

  CLS <- MSIG.ReadPhenFile(file = input.cls) # Read phenotype file (CLS format)
  cls.labels <- CLS$class.v
  cls.phen <- CLS$phen
  cls.list <- CLS$class.list

  if (is.vector(cls.labels)) {
    n.phen <- 1
  } else {
    n.phen <- length(cls.labels[, 1])
  }
  if (!is.na(user.colors)) {
    c.test <- user.colors
  } else {
    if (!is.null(CLS$col.phen)) {
      c.test <- CLS$col.phen
    } else {
      c.test <- c(
        brewer.pal(n = 7, name = "Dark2"), brewer.pal(n = 7, name = "Set2"), brewer.pal(n = 7, name = "Set1"),
        brewer.pal(n = 8, name = "Accent"), brewer.pal(n = 9, name = "Spectral"), brewer.pal(n = 8, name = "Set3"),
        brewer.pal(n = 8, name = "BuGn"),
        brewer.pal(n = 7, name = "Set2"), brewer.pal(n = 7, name = "Dark2"), brewer.pal(n = 7, name = "Set1"),
        brewer.pal(n = 8, name = "Accent"), brewer.pal(n = 10, name = "Spectral"), brewer.pal(n = 8, name = "Set3"),
        brewer.pal(n = 8, name = "BuGn"),
        brewer.pal(n = 7, name = "Set2"), brewer.pal(n = 7, name = "Dark2"), brewer.pal(n = 7, name = "Set1"),
        brewer.pal(n = 8, name = "Accent"), brewer.pal(n = 10, name = "Spectral"), brewer.pal(n = 8, name = "Set3"),
        brewer.pal(n = 8, name = "BuGn")
      )
    }
  }


  if (!is.null(CLS$phen.names)) {
    phen.names <- CLS$phen.names
    #      if (is.vector(cls.list)) {
    #         cls.phen <- paste(phen.names, cls.phen, collapse="_")
    #      } else {
    #         for (i in 1:length(cls.phen)) {
    #            for (j in 1:length(cls.phen[[i]])) {
    #               cls.phen[[i]][j] <- paste(phen.names[i], cls.phen[[i]][j], collapse="_")
    #            }
    #         }
    #      }
  } else {
    phen.names <- "NA"
  }

  cls.phen.index <- unlist(cls.phen)
  cls.phen.colors <- c.test[1:length(cls.phen.index)]

  n.classes <- vector(length = n.phen, mode = "numeric")
  if (n.phen == 1) {
    max.classes <- length(cls.phen)
    n.classes[1] <- max.classes
  } else {
    max.classes <- max(unlist(lapply(cls.phen, FUN = length)))
    for (i in 1:n.phen) {
      n.classes[i] <- length(cls.phen[[i]])
    }
  }

  x <- rbind(sample.names, cls.list, cls.labels)
  print("before loop")
  print(x)
  print(cls.phen)
  print(phen.names)

  filename <- paste(results.dir, test.file.prefix, ".PHEN.MARKERS.", markers.metric, sep = "")
  #               pdf(file=glob.filename, height = 10, width = 10)
  pdf(file = paste(filename, ".pdf", sep = ""), height = 8.5, width = 11)


  # Loop over phenotypes

  for (k.phen in 1:n.phen) {
    if (is.vector(cls.labels)) {
      k.phen.labels <- cls.labels
      k.phen.list <- cls.list
    } else {
      k.phen.labels <- as.vector(cls.labels[k.phen, ])
      k.phen.list <- as.vector(cls.list[k.phen, ])
    }

    # Sort according to current phenotype

    if (sort.expression == T) {
      phen.index <- order(k.phen.labels, decreasing = sort.decreasing)
    } else {
      phen.index <- seq(1, length(k.phen.labels))
    }
    if (is.vector(cls.labels)) {
      cls.labels2 <- cls.labels[phen.index]
      cls.list2 <- cls.list[phen.index]
    } else {
      cls.labels2 <- cls.labels[, phen.index]
      cls.list2 <- cls.list[, phen.index]
    }
    k.phen.labels <- k.phen.labels[phen.index]
    k.phen.list <- k.phen.list[phen.index]
    sample.names2 <- sample.names[phen.index]
    m2 <- m[, phen.index]

    x <- rbind(sample.names2, cls.list2, cls.labels2)
    print(paste("inside loop phen=", k.phen))
    print(x)
    print(cls.phen)
    print(phen.names)

    # Markers for each class

    if (is.vector(cls.labels2)) {
      classes <- unique(cls.list2)
    } else {
      classes <- unique(cls.list2[k.phen, ])
    }
    if (length(classes) > 2) {
      k.only.up <- T
    } else {
      k.only.up <- only.up
    }

    if (length(classes) == 2) classes <- classes[1]
    markers <- NULL
    markers.descs <- NULL
    metric.list <- NULL
    p.val.list <- NULL
    k.class <- NULL
    for (k in classes) {
      if (is.vector(cls.labels2)) {
        bin.class <- ifelse(cls.list2 == k, 0, 1)
      } else {
        bin.class <- ifelse(cls.list2[k.phen, ] == k, 0, 1)
      }
      if (markers.metric == "T.TEST") {
        metric <- vector(length = n.models, mode = "numeric")
        p.val <- vector(length = n.models, mode = "numeric")
        for (i in 1:n.models) {
          temp <- split(m2[i, ], bin.class)
          x <- temp[[1]]
          y <- temp[[2]]
          metric[i] <- signif(t.test(x = x, y = y)$statistic, digits = 3)
          p.val[i] <- signif(t.test(x = x, y = y)$p.value, digits = 3)
        }
      } else if (markers.metric == "ROC") {
        bin.class <- ifelse(bin.class == 1, 0, 1)
        metric <- vector(length = n.models, mode = "numeric")
        p.val <- vector(length = n.models, mode = "numeric")
        for (i in 1:n.models) {
          m.score <- m2[i, ]
          m.score.norm <- (m.score - min(m.score)) / (max(m.score) - min(m.score))
          perf.auc <- roc.area(bin.class, m.score.norm)
          metric[i] <- signif(perf.auc$A, digits = 3)
          p.val[i] <- signif(perf.auc$p.value, digits = 3)
        }
      } else if (markers.metric == "MEAN.DIFF") {
        bin.class <- ifelse(bin.class == 1, 0, 1)
        metric <- vector(length = n.models, mode = "numeric")
        p.val <- vector(length = n.models, mode = "numeric")
        for (i in 1:n.models) {
          temp <- split(m2[i, ], bin.class)
          x <- temp[[1]]
          y <- temp[[2]]
          metric[i] <- signif(mean(x) - mean(y), digits = 3)
          p.val[i] <- signif(t.test(x = x, y = y)$p.value, digits = 3)
        }
      }

      if (is.na(sort.decreasing.genes)) {
        metric.order <- seq(1, length(metric))
      } else {
        metric.order <- order(metric, decreasing = sort.decreasing.genes)
      }
      if (only.up == TRUE) {
        if (length(classes) == 2) {
          k.markers.num <- ifelse(markers.num > n.models, n.models, markers.num)
        } else {
          k.markers.num <- ifelse(length(classes) * markers.num > n.models,
            floor(n.models / length(classes)), markers.num
          )
        }
        markers <- c(markers, model.names[metric.order][1:k.markers.num])
        markers.descs <- c(markers.descs, model.descs[metric.order][1:k.markers.num])
        metric.list <- c(metric.list, metric[metric.order][1:k.markers.num])
        p.val.list <- c(p.val.list, p.val[metric.order][1:k.markers.num])
        k.class <- c(k.class, rep(k, k.markers.num))
      } else {
        k.markers.num <- ifelse(length(classes) * markers.num > n.models, floor(n.models / length(classes)),
          markers.num
        )
        markers <- c(
          markers, model.names[metric.order][1:k.markers.num],
          model.names[metric.order][(length(model.names) - k.markers.num + 1):length(model.names)]
        )
        markers.descs <- c(
          markers.descs, model.descs[metric.order][1:k.markers.num],
          model.descs[metric.order][(length(model.names) - k.markers.num + 1):length(model.names)]
        )
        metric.list <- c(
          metric.list, metric[metric.order][1:k.markers.num],
          metric[metric.order][(length(model.names) - k.markers.num + 1):length(model.names)]
        )
        p.val.list <- c(
          p.val.list, p.val[metric.order][1:k.markers.num],
          p.val[metric.order][(length(model.names) - k.markers.num + 1):length(model.names)]
        )
        k.class <- c(k.class, rep(k, k.markers.num), rep(paste("not", k), k.markers.num))
      }
    }

    V3 <- m2[markers, ]
    model.descs2 <- paste(metric.list, p.val.list, k.class, markers.descs)
    height <- ifelse(length(markers) + n.phen >= 9, 10, (length(markers) + n.phen) * 0.44 + 5)
    #      char.res <-  0.0085 * length(markers) + 0.65


    # Sort markers inside each phenotype class

    if (sort.expression == T) {
      for (j in unique(k.phen.labels)) {
        V4 <- V3[, k.phen.labels == j]
        sn <- sample.names2[k.phen.labels == j]
        if (is.vector(cls.labels)) {
          clab <- cls.labels2[k.phen.labels == j]
          clis <- cls.list2[k.phen.labels == j]
        } else {
          clab <- cls.labels2[, k.phen.labels == j]
          clis <- cls.list2[, k.phen.labels == j]
        }
        dist.matrix <- dist(t(V4))
        HC <- hclust(dist.matrix, method = "complete")
        V4 <- V4[, HC$order]
        sn <- sn[HC$order]
        if (is.vector(cls.labels2)) {
          clab <- clab[HC$order]
          clis <- clis[HC$order]
        } else {
          clab <- clab[, HC$order]
          clis <- clis[, HC$order]
        }
        V3[, k.phen.labels == j] <- V4
        sample.names2[k.phen.labels == j] <- sn
        if (is.vector(cls.labels2)) {
          cls.labels2[k.phen.labels == j] <- clab
          cls.list2[k.phen.labels == j] <- clis
        } else {
          cls.labels2[, k.phen.labels == j] <- clab
          cls.list2[, k.phen.labels == j] <- clis
        }
      }
    }
    x <- rbind(sample.names2, cls.list2, cls.labels2)
    print(paste("inside loop after in-class sort phen=", k.phen))
    print(x)
    print(cls.phen)
    print(phen.names)

    # Recompute cls.phen and cls.labels2 as order may have changed

    cls.phen2 <- NULL
    if (is.vector(cls.labels2)) {
      classes <- unique(cls.list2)
      cls.phen2 <- classes
      cls.labels2 <- match(cls.list2, cls.phen2)
    } else {
      for (kk in 1:length(cls.list2[, 1])) {
        classes <- unique(cls.list2[kk, ])
        cls.phen2[[kk]] <- classes
        cls.labels2[kk, ] <- match(cls.list2[kk, ], cls.phen2[[kk]])
      }
    }

    x <- rbind(sample.names2, cls.list2, cls.labels2)
    print(paste("inside loop after cls.phen renorm phen=", k.phen))
    print(cls.phen2)
    print(phen.names)


    library(gmodels)
    if (!is.vector(cls.labels2)) {
      if (sort.phenotypes == T) {
        phen.score <- vector(length = n.phen, mode = "numeric")
        for (k.lab in 1:n.phen) {
          tab <- table(as.vector(cls.list2[k.lab, ]), k.phen.list)
          print(tab)
          #              phen.score[k.lab] <- 1 - chisq.test(tab)$p.value
          #              phen.score[k.lab] <- 1 - fisher.test(tab)$p.value
          CT <- CrossTable(tab, chisq = T)
          phen.score[k.lab] <- CT$chisq$p.value
          print(phen.score[k.lab])
        }
        phen.order <- order(phen.score, decreasing = T)
        print(phen.order)
        cls.labels2 <- cls.labels2[phen.order, ]
        cls.phen2 <- cls.phen2[phen.order]
        phen.names2 <- phen.names[phen.order]
        main.string <- paste(test.file.prefix, " - ", phen.names2[n.phen], markers.metric, " order")
      } else {
        phen.names2 <- phen.names
        main.string <- paste(test.file.prefix, " - ", phen.names2[k.phen], markers.metric, " order")
      }
    } else {
      phen.names2 <- phen.names[1]
      main.string <- paste(test.file.prefix, " - ", phen.names2, markers.metric, " order")
    }

    #     x11(width=15, height=height)


    x <- rbind(sample.names2, cls.list2, cls.labels2)
    print(paste("inside loop after phen sort before figure phen=", k.phen))
    print(x)
    print(cls.phen2)
    print(phen.names2)

    phen.list <- unlist(cls.phen2)
    colors.list <- cls.phen.colors[match(phen.list, cls.phen.index)]

    print(rbind(phen.list, colors.list))

    markers <- paste(markers, seq(1, length(markers)), sep = "_")
    MSIG.HeatMapPlot.7(
      V = V3, row.names = markers,
      row.names2 = model.descs2, col.labels = cls.labels2,
      col.classes = cls.phen2, phen.cmap = colors.list, phen.names = phen.names2,
      col.names = sample.names2, main = main.string, xlab = "  ", ylab = "  ",
      row.norm = row.norm,
      cmap.type = cmap.type, char.rescale = char.res, legend = legend
    )

    V3 <- data.frame(V3)
    colnames(V3) <- sample.names2
    row.names(V3) <- markers

    if (!is.null(markers.file)) {
      write.gct(gct.data.frame = V3, descs = model.descs2, filename = markers.file)
    }
  } # end loop over phenotypes

  dev.off()
}

MSIG.HeatMapPlot.7 <- function(
    V,
    row.names = "NA",
    row.names2 = "NA",
    col.labels = "NA",
    col.labels2 = "NA",
    col.classes = "NA",
    phen.cmap = "NA",
    col.names = "NA",
    phen.names = "NA",
    main = " ",
    sub = " ",
    xlab = " ",
    ylab = " ",
    row.norm = TRUE,
    char.rescale = 0.85,
    cmap.type = 1, # 1 = vintage pinkogram, 2 = scale of blues, 3 = high-resolution pinkogram for scores or probabilities [0, 1], 4 = high-resolution pinkogram for general values, 5 = color map for normalized enrichment scores, 6 = scale of red purples, 7 = scale of Oranges, 8 = scale of Greens, 9 = scale of Blues
    max.v = "NA",
    legend = T) {
  #
  # Plots a heatmap "pinkogram" of a gene expression matrix including phenotype vector and gene, sample and phenotype labels
  #
  # The Broad Institute
  # SOFTWARE COPYRIGHT NOTICE AGREEMENT
  # This software and its documentation are copyright 2003 by the
  # Broad Institute/Massachusetts Institute of Technology.
  # All rights are reserved.
  #
  # This software is supplied without any warranty or guaranteed support
  # whatsoever. Neither the Broad Institute nor MIT can be responsible for
  # its use, misuse, or functionality.

  n.rows <- length(V[, 1])
  n.cols <- length(V[1, ])
  V1 <- matrix(0, nrow = n.rows, ncol = n.cols)

  #       if ((cmap.type == 5) | (cmap.type == 3)) {
  if (cmap.type == 5) {
    row.norm <- F
  }

  if (row.norm == TRUE) {
    row.mean <- apply(V, MARGIN = 1, FUN = mean)
    row.sd <- apply(V, MARGIN = 1, FUN = sd)
    row.n <- length(V[, 1])
    for (i in 1:n.rows) {
      if (row.sd[i] == 0) {
        V1[i, ] <- 0
      } else {
        V1[i, ] <- (V[i, ] - row.mean[i]) / (0.333 * row.sd[i])
      }
      V1[i, ] <- ifelse(V1[i, ] < -4, -4, V1[i, ])
      V1[i, ] <- ifelse(V1[i, ] > 4, 4, V1[i, ])
    }
  } else {
    V1 <- V
  }

  if (cmap.type == 1) {
    mycol <- c(
      "#0000FF", "#4040FF", "#7070FF", "#8888FF", "#A9A9FF", "#D5D5FF", "#EEE5EE", "#FFAADA",
      "#FF9DB0", "#FF7080",
      "#FF5A5A", "#FF4040", "#FF0D1D"
    ) # blue-pinkogram colors. This is the 1998-vintage,
    # pre-gene cluster, original pinkogram color map
  } else if (cmap.type == 2) {
    violet.palette <- colorRampPalette(c("#400030", "white"), space = "rgb")
    mycol <- rev(violet.palette(20))

    #          mycol <- c("#FCFBFD","#F4F2F8","#F8F7FB","#EFEDF5","#E1E1EF","#E8E7F2","#DADAEB","#C6C7E1","#D0D1E6",
    #                        "#BCBDDC","#A8A6CF",
    #                        "#B2B2D6","#9E9AC8","#8A87BF","#9491C4","#807DBA","#7260AB","#796FB3","#6A51A3","#5C3596",
    #                        "#63439D","#54278F","#460D83","#4D1A89","#3F007D")
  } else if (cmap.type == 6) {
    mycol <- c(
      "#FFF7F3", "#FDE0DD", "#FCC5C0", "#FA9FB5", "#F768A1", "#DD3497", "#AE017E",
      "#7A0177", "#49006A"
    )
  } else if (cmap.type == 7) {
    mycol <- c(
      "#FFF5EB", "#FEE6CE", "#FDD0A2", "#FDAE6B", "#FD8D3C", "#F16913", "#D94801",
      "#A63603", "#7F2704"
    )
  } else if (cmap.type == 8) {
    mycol <- c(
      "#F7FCF5", "#E5F5E0", "#C7E9C0", "#A1D99B", "#74C476", "#41AB5D", "#238B45",
      "#006D2C", "#00441B"
    )
  } else if (cmap.type == 9) {
    mycol <- c(
      "#F7FBFF", "#DEEBF7", "#C6DBEF", "#9ECAE1", "#6BAED6", "#4292C6", "#2171B5",
      "#08519C", "#08306B"
    )
  } else if ((cmap.type == 3) | (cmap.type == 4) | (cmap.type == 5)) {
    mycol <- vector(length = 512, mode = "numeric")

    for (k in 1:256) {
      mycol[k] <- rgb(255, k - 1, k - 1, maxColorValue = 255)
    }
    for (k in 257:512) {
      mycol[k] <- rgb(511 - (k - 1), 511 - (k - 1), 255, maxColorValue = 255)
    }
    mycol <- rev(mycol)
  }

  ncolors <- length(mycol)

  if (cmap.type == 5) {
    if (max.v == "NA") {
      max.v <- max(max(V1), -min(V1))
    }
    V2 <- ceiling(ncolors * (V1 - (-max.v)) / (1.001 * (max.v - (-max.v))))
  } else {
    V2 <- ceiling(ncolors * (V1 - min(V1)) / (1.001 * (max(V1) - min(V1))))
  }

  if (col.labels[1] == "NA") {
    heatm <- matrix(0, nrow = n.rows + 1, ncol = n.cols)
    heatm[1:n.rows, ] <- V2[seq(n.rows, 1, -1), ]
    tot.cols <- ncolors
    if (legend == T) {
      nf <- layout(matrix(c(1, 2, 3, 0), 2, 2, byrow = T), widths = c(5, 1), heights = c(10, 1), respect = FALSE)
    } else {
      nf <- layout(matrix(c(1, 2), 2, 1, byrow = T), heights = c(8, 1), respect = FALSE)
    }
    par(mar = c(3, 16, 3, 16))
    mycol <- c(mycol, phen.cmap[1:length(col.classes)])
    image(1:n.cols, 1:n.rows, t(heatm),
      zlim = c(0, tot.cols), col = mycol, axes = FALSE,
      main = main, sub = sub, xlab = xlab, ylab = ylab
    )
    n.rows.phen <- 0
  } else {
    tot.cols <- ncolors
    if (is.vector(col.labels)) {
      heatm <- matrix(0, nrow = n.rows + 1, ncol = n.cols)
      heatm[1:n.rows, ] <- V2[seq(n.rows, 1, -1), ]
      n.rows.phen <- 1
      heatm[n.rows + 1, ] <- tot.cols + col.labels
      cols.row <- length(unique(col.labels))
      tot.cols <- tot.cols + cols.row
      phen.cmap <- phen.cmap[1:cols.row]
    } else {
      n.rows.phen <- length(col.labels[, 1])
      cols.row <- vector(length = n.rows.phen, mode = "numeric")
      heatm <- matrix(0, nrow = n.rows + n.rows.phen, ncol = n.cols)
      heatm[1:n.rows, ] <- V2[seq(n.rows, 1, -1), ]
      for (k in seq(n.rows + n.rows.phen, n.rows + 1, -1)) {
        heatm[k, ] <- tot.cols + col.labels[n.rows + n.rows.phen - k + 1, ]
        cols.row[n.rows + n.rows.phen - k + 1] <- length(unique(col.labels[n.rows + n.rows.phen - k + 1, ]))
        tot.cols <- tot.cols + cols.row[n.rows + n.rows.phen - k + 1]
        #                 print(c("col:", k, ":", tot.cols + col.labels[n.rows + n.rows.phen - k + 1,], "tot.cols:", tot.cols))
      }
      phen.cmap <- phen.cmap[1:sum(unlist(lapply(col.classes, length)))]
    }
    if (legend == T) {
      #              nf <- layout(matrix(c(1, 2, 3, 0), 2, 2, byrow=T), widths = c(10, 2), heights = c(6, 1), respect = FALSE)
      nf <- layout(matrix(c(1, 2, 3), 3, 1, byrow = T), heights = c(8, 4, 1), respect = FALSE)
    } else {
      nf <- layout(matrix(c(1, 2), 2, 1, byrow = T), heights = c(5, 1), respect = FALSE)
    }
    par(mar = c(3, 16, 3, 16))
    mycol <- c(mycol, phen.cmap)
    image(1:n.cols, 1:(n.rows + n.rows.phen), t(heatm),
      zlim = c(0, tot.cols), col = mycol, axes = FALSE,
      main = main, sub = sub, xlab = xlab, ylab = ylab
    )
  }

  # Add lines to separate phenotypes or subgroups

  if (col.labels2[1] != "NA") {
    groups <- split(col.labels2, col.labels2)
    len.vec <- lapply(groups, length)
    plot.div <- c(0.51, cumsum(len.vec) + 0.5)
    for (i in plot.div) {
      lines(c(i, i), c(0, n.rows + n.rows.phen + 0.48), lwd = 2, cex = 0.9, col = "black")
    }
    lines(c(0.51, n.cols + 0.49), c(0.51, 0.51), lwd = 2, cex = 0.9, col = "black")
    lines(c(0.51, n.cols + 0.49), c(n.rows + n.rows.phen + 0.48, n.rows + n.rows.phen + 0.48),
      lwd = 2,
      cex = 0.9, col = "black"
    )
    lines(c(0.51, n.cols + 0.49), c(n.rows + 0.50, n.rows + 0.50),
      lwd = 2,
      cex = 0.9, col = "black"
    )
  }
  if (row.names[1] != "NA") {
    numC <- nchar(row.names)
    size.row.char <- char.rescale * 25 / (n.rows + 20)
    for (i in 1:n.rows) {
      row.names[i] <- substr(row.names[i], 1, 40)
      row.names[i] <- paste(row.names[i], " ", sep = "")
    }
    if (phen.names[1] == "NA") {
      head.names <- paste("Class", seq(n.rows.phen, 1, -1))
    } else {
      head.names <- as.character(rev(phen.names))
    }
    row.names <- c(row.names[seq(n.rows, 1, -1)], head.names)
    #            print(paste("n.rows:", n.rows))
    #            print(paste("Phen names:", phen.names))
    #            print(paste("Head names:", head.names))
    #            print(paste("Row names:", row.names))
    axis(2,
      at = 1:(n.rows + n.rows.phen), labels = row.names, adj = 0.5, tick = FALSE, las = 1, cex.axis = size.row.char,
      font.axis = 2, line = -1
    )
  }

  if (row.names2[1] != "NA") {
    numC <- nchar(row.names2)
    size.row.char <- char.rescale * 25 / (n.rows + 20)
    for (i in 1:n.rows) {
      row.names2[i] <- substr(row.names2[i], 1, 40)
      row.names2[i] <- paste(" ", row.names2[i], sep = "")
    }
    row.names2 <- c(row.names2[seq(n.rows, 1, -1)], "   ")
    axis(4,
      at = 1:(n.rows + 1), labels = row.names2, adj = 0.5, tick = FALSE, las = 1,
      cex.axis = size.row.char, font.axis = 2, line = -1
    )
  }

  if (col.names[1] != "NA") {
    size.col.char <- char.rescale * 20 / (n.cols + 25)
    axis(1, at = 1:n.cols, labels = col.names, tick = FALSE, las = 3, cex.axis = size.col.char, font.axis = 2, line = -1)
  }

  # Phenotype Legend

  #      print("--------------------------------------------------------------------------------------------")
  if (legend == T) {
    leg.txt <- NULL
    p.vec <- NULL
    c.vec <- NULL
    c2.vec <- NULL
    ind <- 1
    par(mar = c(0, 0, 0, 0))
    plot(c(0, 0), c(1, 1), xlim = c(0, 1), ylim = c(-.050, 1.05), axes = F, type = "n", xlab = "", ylab = "")
    for (i in 1:n.rows.phen) {
      if (is.vector(col.labels)) {
        phen.v <- as.character(col.classes)
      } else {
        phen.v <- as.character(col.classes[[i]])
      }
      p.name <- paste(as.character(rev(head.names)[i]), ":   ", sep = "")
      leg.txt <- c(p.name, phen.v)
      p.vec <- rep(22, cols.row[i] + 1)
      c.vec <- c("#FFFFFF", phen.cmap[ind:(ind + cols.row[i] - 1)])
      c2.vec <- c("#FFFFFF", rep("black", cols.row[i]))
      ind <- ind + cols.row[i]
      offset <- 0.07
      legend(
        x = 0, y = 1 - offset * i,
        horiz = T, x.intersp = 0.5, legend = leg.txt, bty = "n", xjust = 0, yjust = 1, pch = p.vec,
        pt.bg = c.vec, col = c2.vec, cex = 1.20, pt.cex = 1.75
      )
    }
  }

  # Color map legend

  #       print(c("range V=", range(V)))
  #       print(c("range V1=", range(V1)))
  #       print(c("range V2=", range(V2)))

  par(mar = c(2, 12, 2, 12))
  num.v <- 20
  range.v <- range(V2)
  incr <- (range.v[1] - range.v[2]) / (num.v - 1)
  heatm.v <- matrix(rev(seq(range.v[2], range.v[1], incr)), nrow = num.v, ncol = 1)
  image(1:num.v, 1:1, heatm.v,
    zlim = c(0, tot.cols), col = mycol, axes = FALSE,
    main = " ", sub = " ", xlab = ylab, ylab = xlab
  )
  range.v <- range(V1)
  incr <- (range.v[1] - range.v[2]) / (num.v - 1)
  heatm.v2 <- matrix(signif(rev(seq(range.v[2], range.v[1], incr)), digits = 2), nrow = num.v, ncol = 1)
  #          print(c("heatm.v2=", heatm.v2))
  axis(3, at = 1:num.v, labels = heatm.v2, adj = 0.5, tick = FALSE, las = 1, cex.axis = 0.5 * char.rescale, font.axis = 1)

  return()
}

OPAM.create.table.from.GISTIC <- function(
    gistic.ds,
    mapping.ds = NULL,
    output.file) {
  gistic <- read.delim(gistic.ds, header = T, sep = "\t", blank.lines.skip = T, comment.char = "", as.is = T)
  temp <- strsplit(gistic.ds, split = "/") # Extract gistic file name
  s <- length(temp[[1]])
  gistic.file.name <- temp[[1]][s]
  temp <- strsplit(gistic.file.name, split = ".txt")
  gistic.file.prefix <- temp[[1]][1]

  if (!is.null(mapping.ds)) { # rename GISTIC columns according to mapping file
    map <- read.delim(mapping.ds, header = T, sep = "\t", blank.lines.skip = T, comment.char = "", as.is = T)
    common.gistic.map <- intersect(colnames(gistic), map[, 1])
    loc.gistic <- match(common.gistic.map, colnames(gistic))
    loc.map <- match(common.gistic.map, map[, 1])
    colnames(gistic)[loc.gistic] <- map[loc.map, 2]
  }

  gistic2 <- NULL
  for (i in 1:length(row.names(gistic))) {
    if (regexpr(pattern = "log2", gistic[i, 1]) != -1) next # exclude log2
    #       if (regexpr(pattern="Broad", gistic[i, 1]) != -1) next  # exclude Broad events
    str <- paste(substr(gistic[i, 1], 1, 3), gistic[i, 2], sep = "_")
    gistic2 <- rbind(cbind(str, gistic[i, ]), gistic2)
  }
  gistic2[gistic2 == 2] <- 1
  gistic3 <- cbind(gistic2[, 1], gistic2[, 11:length(gistic[1, ])])
  colnames(gistic3) <- c("Sample", colnames(gistic3)[2:length(gistic3[1, ])])
  gistic3 <- t(gistic3)
  write.table(gistic3, file = output.file, quote = F, col.names = F, row.names = T, append = F, sep = "\t")
}

OPAM.sort.projection.by.score <- function(
    input.ds,
    input.cls,
    results.dir,
    normalize.score = T,
    normalization.type = "zero.one",
    model,
    user.colors = NA,
    decreasing.order = T,
    output.dataset = NA) {
  library(gtools)
  library(verification)
  library(ROCR)
  library(MASS)
  library(RColorBrewer)
  library(heatmap.plus)

  dataset <- MSIG.Gct2Frame(filename = input.ds) # Read gene expression dataset (GCT format)
  m <- data.matrix(dataset$ds)
  model.names <- dataset$row.names
  model.descs <- dataset$descs
  Ns <- length(m[1, ])
  dim(m)
  sample.names <- dataset$names

  n.models <- length(m[, 1])
  temp <- strsplit(input.ds, split = "/") # Extract test file name
  s <- length(temp[[1]])
  test.file.name <- temp[[1]][s]
  temp <- strsplit(test.file.name, split = ".gct")
  test.file.prefix <- temp[[1]][1]
  char.res <- 0.013 * n.models + 0.65

  # normalize scores

  if (normalize.score == T) {
    if (normalization.type == "zero.one") {
      for (i in 1:n.models) {
        m[i, ] <- (m[i, ] - min(m[i, ])) / (max(m[i, ]) - min(m[i, ]))
      }
    } else if (normalization.type == "z.score") {
      for (i in 1:n.models) {
        m[i, ] <- (m[i, ] - mean(m[i, ])) / sd(m[i, ])
      }
    } else if (normalization.type == "r.z.score") {
      for (i in 1:n.models) {
        m[i, ] <- (m[i, ] - median(m[i, ])) / mad(m[i, ])
      }
    }
  }

  CLS <- MSIG.ReadPhenFile(file = input.cls) # Read phenotype file (CLS format)
  cls.labels <- CLS$class.v
  cls.phen <- CLS$phen
  cls.list <- CLS$class.list

  if (is.vector(cls.labels)) {
    n.phen <- 1
  } else {
    n.phen <- length(cls.labels[, 1])
  }
  if (!is.na(user.colors)) {
    c.test <- user.colors
  } else {
    if (!is.null(CLS$col.phen)) {
      c.test <- CLS$col.phen
    } else {
      c.test <- c(
        brewer.pal(n = 7, name = "Set2"), brewer.pal(n = 7, name = "Dark2"), brewer.pal(n = 7, name = "Set1"),
        brewer.pal(n = 8, name = "Accent"), brewer.pal(n = 9, name = "Spectral"), brewer.pal(n = 8, name = "Set3"),
        brewer.pal(n = 8, name = "BuGn"),
        brewer.pal(n = 7, name = "Set2"), brewer.pal(n = 7, name = "Dark2"), brewer.pal(n = 7, name = "Set1"),
        brewer.pal(n = 8, name = "Accent"), brewer.pal(n = 10, name = "Spectral"), brewer.pal(n = 8, name = "Set3"),
        brewer.pal(n = 8, name = "BuGn"),
        brewer.pal(n = 7, name = "Set2"), brewer.pal(n = 7, name = "Dark2"), brewer.pal(n = 7, name = "Set1"),
        brewer.pal(n = 8, name = "Accent"), brewer.pal(n = 10, name = "Spectral"), brewer.pal(n = 8, name = "Set3"),
        brewer.pal(n = 8, name = "BuGn")
      )
    }
  }


  if (!is.null(CLS$phen.names)) {
    phen.names <- CLS$phen.names
  } else {
    phen.names <- "NA"
  }

  cls.phen.index <- unlist(cls.phen)
  cls.phen.colors <- c.test[1:length(cls.phen.index)]

  n.classes <- vector(length = n.phen, mode = "numeric")
  if (n.phen == 1) {
    max.classes <- length(cls.phen)
    n.classes[1] <- max.classes
  } else {
    max.classes <- max(unlist(lapply(cls.phen, FUN = length)))
    for (i in 1:n.phen) {
      n.classes[i] <- length(cls.phen[[i]])
    }
  }

  filename <- paste(results.dir, test.file.prefix, ".SORT.PROJ", sep = "")
  #    pdf(file=paste(filename, ".pdf", sep=""), height = 11, width = 9)
  pdf(file = paste(filename, ".pdf", sep = ""), height = 8.5, width = 11)
  #   x11(width=12, height=8)

  loc <- match(model, model.names)
  print(c("loc:", loc))
  s.order <- order(m[loc, ], decreasing = decreasing.order)
  m2 <- m[, s.order]

  sample.names2 <- sample.names[s.order]

  if (is.vector(cls.labels)) {
    cls.labels2 <- cls.labels[s.order]
    cls.list2 <- cls.list[s.order]
  } else {
    cls.labels2 <- cls.labels[, s.order]
    cls.list2 <- cls.list[, s.order]
  }
  # Recompute cls.phen and cls.labels2 as order may have changed

  cls.phen2 <- NULL
  if (is.vector(cls.labels)) {
    classes <- unique(cls.list2)
    cls.phen2 <- classes
    cls.labels2 <- match(cls.list2, cls.phen2)
  } else {
    for (kk in 1:length(cls.list2[, 1])) {
      classes <- unique(cls.list2[kk, ])
      cls.phen2[[kk]] <- classes
      cls.labels2[kk, ] <- match(cls.list2[kk, ], cls.phen2[[kk]])
    }
  }


  correl <- cor(t(m2))[, loc]
  m.order <- order(correl, decreasing = T)
  correl2 <- correl[m.order]
  m2 <- m2[m.order, ]
  model.names2 <- model.names[m.order]
  model.descs2 <- paste(model.descs[m.order], signif(correl2, digits = 3))
  phen.list <- unlist(cls.phen2)
  colors.list <- cls.phen.colors[match(phen.list, cls.phen.index)]

  MSIG.HeatMapPlot.7(
    V = m2, row.names = model.names2,
    row.names2 = model.descs2, col.labels = cls.labels2,
    col.classes = cls.phen2, phen.cmap = colors.list, phen.names = phen.names,
    col.names = sample.names2, main = " ", xlab = "  ", ylab = "  ", row.norm = T,
    cmap.type = 3, char.rescale = 1, legend = T
  )

  dev.off()

  if (!is.na(output.dataset)) {
    V.GCT <- m2
    colnames(V.GCT) <- sample.names2
    row.names(V.GCT) <- model.names2
    write.gct(gct.data.frame = V.GCT, descs = model.descs2, filename = output.dataset)
  }
}

OPAM.create.model.2 <- function(
    input.ds,
    input.cls,
    input2.ds,
    input2.cls,
    models.dir,
    target.class,
    target.class2,
    model.name,
    model.description,
    # input.signature.file <- NA    # file with gene sets
    # input.gene.set.up <- NA       # name of UP signature
    # input.gene.set.dn <- NA       # name of DN signature
    sample.norm.type = "rank", # "rank", "log.rank", "log"
    marker.disc = "MEAN.DIFF",
    top.markers.up = 20,
    top.markers.dn = 20,
    top.markers.up2 = 20,
    top.markers.dn2 = 20,
    statistic = "area.under.RES",
    weight = 0.25,
    msig.cntrl.genes = NA,
    random.seed = 12345,
    nperm = 0,
    link.function = "logit",
    burnin.iter = 5000, # number of burnin iteration in MCMC (default: 5000)
    mcmc.iter = 25000) # number of MCMC iterations (default: 25000)
{ #----------------------------------------------------------------------------------------

  # Load libraries

  erf <- function(x) 2 * pnorm(x * sqrt(2)) - 1
  library(gtools)
  library(verification)
  library(ROCR)
  library(MASS)
  library(RColorBrewer)
  library(MCMCpack)

  #  Set parameters
  c1 <- c("black", "lightgrey")
  set.seed(random.seed)
  models.file <- paste(models.dir, "/", model.name, ".mod", sep = "")

  # Read dataset
  dataset <- MSIG.Gct2Frame(filename = input.ds) # Read gene expression dataset (GCT format)
  m <- data.matrix(dataset$ds)
  gene.names <- dataset$row.names
  gene.descs <- dataset$descs
  sample.names <- dataset$names
  Ns <- length(m[1, ])
  Ng <- length(m[, 1])

  # Read phenotype and redefine classes: target.class vs. cntrl
  CLS <- MSIG.ReadClsFile(file = input.cls) # Read phenotype file (CLS format)
  class.labels <- CLS$class.v
  class.phen <- CLS$phen
  class.list <- CLS$class.list
  if (is.na(match(target.class, class.phen))) stop(c("target class is not phenotype in:", input.cls))
  for (i in 1:length(class.list)) class.labels[i] <- ifelse(class.list[i] == target.class, 1, 2)
  col.index <- order(class.labels, decreasing = F)
  for (j in 1:Ng) m[j, ] <- m[j, col.index] # reorder samples with target.class first
  sample.names <- sample.names[col.index]
  class.labels <- class.labels[col.index]
  class.list <- class.list[col.index]
  class.phen <- c(target.class, "CNTL")
  control.class <- "CNTL"

  # Sample normalization
  if (sample.norm.type == "rank") {
    for (j in 1:Ns) { # column rank normalization
      m[, j] <- rank(m[, j], ties.method = "average")
    }
    m <- 10000 * m / Ng
  } else if (sample.norm.type == "log.rank") {
    for (j in 1:Ns) { # column rank normalization
      m[, j] <- rank(m[, j], ties.method = "average")
    }
    m <- log(10000 * m / Ng + exp(1))
  } else if (sample.norm.type == "log") {
    m[m < 1] <- 1
    m <- log(m + exp(1))
  }

  # Control signature normalization
  if (!is.na(msig.cntrl.genes)) {
    gene.names.int <- intersect(msig.cntrl.genes, gene.names)
    locs <- match(gene.names.int, gene.names, nomatch = 0)
    msig.cntrl <- m[locs, ]
    msig.cntrl.genes <- gene.names[locs]
    msig.cntrl.descs <- gene.descs[locs]
    msig.cntrl.size <- length(locs)
    if (msig.cntrl.size < 1) {
      msig.cntrl.center <- rep(1, Ns)
    } else if (msig.cntrl.size == 1) {
      msig.cntrl.center <- msig.cntrl
    } else if (msig.cntrl.size > 1) msig.cntrl.center <- apply(msig.cntrl, MARGIN = 2, FUN = mean)
    for (i in 1:Ng) {
      m[i, ] <- m[i, ] / msig.cntrl.center
    }
  }

  # Obtain UP & DN signatures
  temp <- Gene.ranking(m, class.labels, method = marker.disc)
  gene.index <- order(temp, decreasing = T)
  gene.scores <- temp[gene.index]
  msig.up <- m[gene.index[1:top.markers.up], ]
  msig.up.size <- top.markers.up
  msig.up.genes <- gene.names[gene.index[1:top.markers.up]]
  msig.up.descs <- gene.descs[gene.index[1:top.markers.up]]
  msig.dn <- m[gene.index[seq(Ng, Ng - top.markers.dn + 1, -1)], ]
  msig.dn.size <- top.markers.dn
  msig.dn.genes <- gene.names[gene.index[seq(Ng, Ng - top.markers.dn + 1, -1)]]
  msig.dn.descs <- gene.descs[gene.index[seq(Ng, Ng - top.markers.dn + 1, -1)]]
  print("Signatures (up/dn) created from gene marker selection")
  print(c("msig.up.size:", msig.up.size))
  print(c("msig.up.genes:", msig.up.genes))
  # print(c("msig.up", msig.up))
  print(c("..."))
  print(c("msig.dn.size:", msig.dn.size))
  print(c("msig.dn.genes:", msig.dn.genes))
  # print(c("msig.dn", msig.dn))
  print(c("..."))
  if (!is.na(msig.cntrl.genes)) {
    print(c("msig.cntrl.size:", msig.cntrl.size))
    print(c("msig.cntrl.genes:", msig.cntrl.genes))
    #    print(c("msig.cntrl[1,]", msig.cntrl[1,]))
  }

  # Plot signatures
  x11(height = 9, width = 12)
  MSIG.HeatMapPlot.3(
    V = msig.up, row.names = msig.up.genes, col.labels = class.labels, col.classes = class.phen,
    phen.cmap = c1, col.names = sample.names, main = paste(model.name, " UP signature"),
    xlab = " ", ylab = " ", sub = " ", row.norm = T, cmap.type = 4, char.rescale = 1
  )

  x11(height = 9, width = 12)
  MSIG.HeatMapPlot.3(
    V = msig.dn, row.names = msig.dn.genes, col.labels = class.labels, col.classes = class.phen,
    phen.cmap = c1, col.names = sample.names, main = paste(model.name, " DN signature"),
    xlab = " ", ylab = " ", sub = " ", row.norm = T, cmap.type = 4, char.rescale = 1
  )

  if (!is.na(msig.cntrl.genes)) {
    x11(height = 9, width = 12)
    MSIG.HeatMapPlot.3(
      V = msig.cntrl, row.names = msig.cntrl.genes, col.labels = class.labels, col.classes = class.phen,
      phen.cmap = c1, col.names = sample.names, main = paste(model.name, " CNTRL signature"),
      xlab = " ", ylab = " ", sub = " ", row.norm = T, cmap.type = 4, char.rescale = 1
    )
  }

  ######
  # Read refinement dataset and redefine classes: target.class2 vs. cntrl2

  dataset <- MSIG.Gct2Frame(filename = input2.ds) # Read second gene expression dataset (GCT format)
  m2 <- data.matrix(dataset$ds)
  gene.names2 <- dataset$row.names
  gene.descs2 <- dataset$descs
  sample.names2 <- dataset$names
  Ns2 <- length(m2[1, ])
  Ng2 <- length(m2[, 1])

  # Extract file name2
  temp <- strsplit(input2.ds, split = "/")
  s <- length(temp[[1]])
  temp <- temp[[1]][4]
  temp <- strsplit(temp, split = ".gct")
  results.prefix2 <- temp[[1]][1]
  results.prefix2

  CLS <- MSIG.ReadClsFile(file = input2.cls) # Read phenotype file2 (CLS format)
  class.labels2 <- CLS$class.v
  class.phen2 <- CLS$phen
  class.list2 <- CLS$class.list
  if (is.na(match(target.class2, class.phen2))) stop(c("target class2 is not phenotype in:", input2.cls))
  for (i in 1:length(class.list2)) class.labels2[i] <- ifelse(class.list2[i] == target.class2, 1, 2)
  col.index <- order(class.labels2, decreasing = F)
  for (j in 1:Ng2) m2[j, ] <- m2[j, col.index] # reorder samples with target.class first
  sample.names2 <- sample.names2[col.index]
  class.labels2 <- class.labels2[col.index]
  class.list2 <- class.list2[col.index]
  class.phen2 <- c(target.class2, "CNTL")
  control.class2 <- "CNTL"

  # Sample normalization2
  if (sample.norm.type == "rank") {
    for (j in 1:Ns2) { # column rank normalization
      m2[, j] <- rank(m2[, j], ties.method = "average")
    }
    m2 <- 10000 * m2 / Ng2
  } else if (sample.norm.type == "log.rank") {
    for (j in 1:Ns2) { # column rank normalization
      m2[, j] <- rank(m2[, j], ties.method = "average")
    }
    m2 <- log(10000 * m2 / Ng2 + exp(1))
  } else if (sample.norm.type == "log") {
    m2[m2 < 1] <- 1
    m2 <- log(m2 + exp(1))
  }

  # Control signature normalization2

  if (!is.na(msig.cntrl.genes)) {
    gene.names.int2 <- intersect(msig.cntrl.genes, gene.names2)
    locs <- match(gene.names.int2, gene.names2, nomatch = 0)
    msig.cntrl2 <- m2[locs, ]
    msig.cntrl.genes2 <- gene.names2[locs]
    msig.cntrl.descs2 <- gene.descs2[locs]
    msig.cntrl.size2 <- length(locs)
    if (msig.cntrl.size2 < 1) {
      msig.cntrl.center2 <- rep(1, Ns2)
    } else if (msig.cntrl.size2 == 1) {
      msig.cntrl.center2 <- msig.cntrl2
    } else if (msig.cntrl.size2 > 1) msig.cntrl.center2 <- apply(msig.cntrl2, MARGIN = 2, FUN = mean)
    for (i in 1:Ng2) {
      m2[i, ] <- m2[i, ] / msig.cntrl.center2
    }
  }

  gene.names.int2 <- intersect(msig.up.genes, gene.names2)
  locs <- match(gene.names.int2, gene.names2, nomatch = 0)
  msig.up2 <- m2[locs, ]
  msig.up.genes2 <- gene.names2[locs]
  msig.up.descs2 <- gene.descs2[locs]
  msig.up.size2 <- length(locs)

  gene.names.int2 <- intersect(msig.dn.genes, gene.names2)
  locs <- match(gene.names.int2, gene.names2, nomatch = 0)
  msig.dn2 <- m2[locs, ]
  msig.dn.genes2 <- gene.names2[locs]
  msig.dn.descs2 <- gene.descs2[locs]
  msig.dn.size2 <- length(locs)

  # Plot signatures2
  x11(height = 9, width = 12)
  MSIG.HeatMapPlot.3(
    V = msig.up2, row.names = msig.up.genes2, col.labels = class.labels2, col.classes = class.phen2,
    phen.cmap = c1, col.names = sample.names2, main = paste(model.name, " UP signature2"),
    xlab = " ", ylab = " ", sub = " ", row.norm = T, cmap.type = 4, char.rescale = 1
  )

  x11(height = 9, width = 12)
  MSIG.HeatMapPlot.3(
    V = msig.dn2, row.names = msig.dn.genes2, col.labels = class.labels2, col.classes = class.phen2,
    phen.cmap = c1, col.names = sample.names2, main = paste(model.name, " DN signature2"),
    xlab = " ", ylab = " ", sub = " ", row.norm = T, cmap.type = 4, char.rescale = 1
  )

  if (!is.na(msig.cntrl.genes)) {
    x11(height = 9, width = 12)
    MSIG.HeatMapPlot.3(
      V = msig.cntrl2, row.names = msig.cntrl.genes2, col.labels = class.labels2, col.classes = class.phen2,
      phen.cmap = c1, col.names = sample.names2, main = paste(model.name, " CNTRL signature2"),
      xlab = " ", ylab = " ", sub = " ", row.norm = T, cmap.type = 4, char.rescale = 1
    )
  }

  # Refine signature by second marker selection on dataset2 (using only marker genes)

  m3 <- rbind(msig.up2, msig.dn2)
  Ng3 <- length(m3[, 1])
  gene.names3 <- c(msig.up.genes2, msig.dn.genes2)
  gene.descs3 <- c(msig.up.descs2, msig.dn.descs2)

  # Obtain UP3 & DN3 signatures
  temp <- Gene.ranking(m3, class.labels2, method = marker.disc)
  gene.index <- order(temp, decreasing = T)
  gene.scores3 <- temp[gene.index]
  msig.up3 <- m3[gene.index[1:top.markers.up2], ]
  msig.up.size3 <- top.markers.up2
  msig.up.genes3 <- gene.names3[gene.index[1:top.markers.up2]]
  msig.up.descs3 <- gene.descs3[gene.index[1:top.markers.up2]]
  msig.dn3 <- m3[gene.index[seq(Ng3, Ng3 - top.markers.dn2 + 1, -1)], ]
  msig.dn.size3 <- top.markers.dn2
  msig.dn.genes3 <- gene.names3[gene.index[seq(Ng3, Ng3 - top.markers.dn2 + 1, -1)]]
  msig.dn.descs3 <- gene.descs3[gene.index[seq(Ng3, Ng3 - top.markers.dn2 + 1, -1)]]
  print("Signatures3 (up/dn) created from gene marker selection")
  print(c("msig.up.size3:", msig.up.size3))
  print(c("msig.up.genes3:", msig.up.genes3))
  # print(c("msig.up3", msig.up3))
  print(c("..."))
  print(c("msig.dn.size3:", msig.dn.size3))
  print(c("msig.dn.genes3:", msig.dn.genes3))
  # print(c("msig.dn3", msig.dn3))
  print(c("..."))

  # Plot signatures3
  x11(height = 9, width = 12)
  MSIG.HeatMapPlot.3(
    V = msig.up3, row.names = msig.up.genes3, col.labels = class.labels2, col.classes = class.phen2,
    phen.cmap = c1, col.names = sample.names2, main = paste(model.name, " UP signature3"),
    xlab = " ", ylab = " ", sub = " ", row.norm = T, cmap.type = 4, char.rescale = 1
  )

  x11(height = 9, width = 12)
  MSIG.HeatMapPlot.3(
    V = msig.dn3, row.names = msig.dn.genes3, col.labels = class.labels2, col.classes = class.phen2,
    phen.cmap = c1, col.names = sample.names2, main = paste(model.name, " DN signature3"),
    xlab = " ", ylab = " ", sub = " ", row.norm = T, cmap.type = 4, char.rescale = 1
  )


  # Project refinement dataset
  OPAM <- OPAM.Projection(m2, gene.names2, Ns2, Ng2, weight, statistic, msig.up.genes3, nperm = nperm)
  score.up <- OPAM$ES.vector
  # score.up <- sign(OPAM$ES.vector)/OPAM$p.val.vector
  OPAM <- OPAM.Projection(m2, gene.names2, Ns2, Ng2, weight, statistic, msig.dn.genes3, nperm = nperm)
  score.dn <- OPAM$ES.vector
  # score.dn <- sign(OPAM$ES.vector)/OPAM$p.val.vector
  score <- score.up - score.dn

  x11(width = 14, height = 9)
  nf <- layout(matrix(c(1, 2, 3, 0, 4, 0), 3, 2, byrow = T), widths = c(5, 1), heights = c(1, 1, 1), respect = FALSE)
  par(mar = c(2, 4, 2, 4))
  barplot(score.up,
    main = "OPAM Score UP (refinement)", font.axis = 1.25, cex.lab = 1.5, cex.axis = 1.25, cex.names =
      1.25, width = 1, space = 0, col = c1[class.labels]
  )
  leg.txt <- class.phen
  p.vec <- rep(22, length(leg.txt))
  par(mar = c(0, 0, 0, 0))
  plot(c(0, 0), c(1, 1), xlim = c(0, 1), ylim = c(0, 1), axes = F, type = "n", xlab = "", ylab = "")
  legend(
    x = 0, y = 0.8, legend = leg.txt, bty = "n", xjust = 0, yjust = 1, pch = p.vec, pt.bg = c1, col = "black",
    cex = 1.25, pt.cex = 2.5
  )
  par(mar = c(2, 4, 2, 4))
  barplot(score.dn,
    main = "OPAM Score DOWN (refinement)", font.axis = 1.25, cex.lab = 1.5, cex.axis = 1.25,
    cex.names = 1.25, width = 1, space = 0, col = c1[class.labels]
  )
  par(mar = c(2, 4, 2, 4))
  barplot(score,
    main = "OPAM Total Score (refinement)", font.axis = 1.25, cex.lab = 1.5, cex.axis = 1.25,
    cex.names = 1.25, width = 1, space = 0, col = c1[class.labels]
  )

  # Fit MCMC logit or probit model to model score using the refinement dataset
  target.var <- ifelse(class.list2 == target.class2, 1, 0)
  Bayesian.function <- ifelse(link.function == "logit", "MCMClogit(", "MCMCprobit(")
  model.formula <- paste(Bayesian.function,
    "target.var ~ score,  burnin = burnin.iter, mcmc = mcmc.iter, bayes.resid=T)",
    sep = ""
  )
  model.formula
  reg.model <- eval(parse(text = model.formula))
  beta.0 <- reg.model[, 1]
  beta.1 <- reg.model[, 2]
  model.formula <- "beta.0 + beta.1 * score[i]"
  model.formula
  prob.i <- matrix(0, nrow = Ns2, ncol = 3)
  model.score <- vector(length = Ns2, mode = "numeric")
  for (i in 1:Ns2) {
    model.score[i] <- eval(parse(text = model.formula))
    if (link.function == "logit") {
      p.vec <- paste("inv.logit(x=", model.formula, ")", sep = "")
    } else if (link.function == "probit") {
      p.vec <- paste("(erf(", model.formula, ") + 1)/2", sep = "")
    } else {
      stop("Unknown link function")
    }
    val <- eval(parse(text = p.vec))
    prob.i[i, 1] <- quantile(val, probs = 0.5)
    prob.i[i, 2] <- quantile(val, probs = 0.05)
    prob.i[i, 3] <- quantile(val, probs = 0.95)
  }
  xmin <- min(model.score)
  xmax <- max(model.score)
  range.x <- xmax - xmin
  n.points <- 1000
  prob.m <- matrix(0, nrow = n.points, ncol = 3)
  x.m <- vector(length = n.points, mode = "numeric")
  for (k in 1:n.points) {
    x.m[k] <- xmin + k * (range.x / n.points)
    if (link.function == "logit") {
      p.vec <- paste("inv.logit(x=", x.m[k], ")", sep = "")
    } else if (link.function == "probit") {
      p.vec <- paste("(erf(", x.m[k], ") + 1)/2", sep = "")
    } else {
      stop("Unknown link function")
    }
    val <- eval(parse(text = p.vec))
    prob.m[k, 1] <- quantile(val, probs = 0.5)
    prob.m[k, 2] <- quantile(val, probs = 0.05)
    prob.m[k, 3] <- quantile(val, probs = 0.95)
  }
  istar <- which.min(abs(0.5 - prob.m[, 1]))
  istar <- xmin + istar * (range.x / 1000)
  x.index <- order(model.score, decreasing = F)
  x.order <- model.score[x.index]
  prob.i.order <- prob.i[x.index, ]
  target.var.order <- ifelse(target.var[x.index] == 1, c1[1], c1[2])
  class.labels.order <- class.labels[x.index]

  # Plot bar graph of z-scores
  boundary <- istar
  pred.class <- ifelse(prob.i.order[, 1] >= 0.5, 2, 1)
  z.range <- range(x.order)
  x11(height = 7, width = 9.5)
  nf <- layout(matrix(c(1, 2), 1, 2, byrow = T), widths = c(3.75, 1), heights = 1, respect = FALSE)
  plot(x.order, prob.i.order[, 1], sub = model.name, pch = 20, col = 0, cex = 2, xlab = "Activation Index", ylab = "Probability")
  points(x.m, prob.m[, 1], type = "l", lwd = 2, col = 1, lty = 1, cex = 1)
  points(x.m, prob.m[, 2], type = "l", col = 4, lty = 1, cex = 1)
  points(x.m, prob.m[, 3], type = "l", col = 4, lty = 1, cex = 1)
  arrows(x.order, prob.i.order[, 2], x.order, prob.i.order[, 3], col = 4, angle = 90, code = 3, length = 0.0)
  range.x <- range(x.order)
  points(range.x, c(0.5, 0.5), type = "l", lty = 3, col = 1, lwd = 2)
  points(range.x, c(-.15, -0.15), type = "l", lty = 1, col = 1, lwd = 2)
  points(c(istar, istar), c(-0.07, 1.07), type = "l", lty = 3, col = 1, lwd = 2)
  points(x.order, prob.i.order[, 1], pch = 21, bg = target.var.order, col = 1, cex = 2)
  points(x.order, rep(-0.15, length(x.order)), pch = 21, bg = target.var.order, col = 1, cex = 2)
  leg.txt <- class.phen
  p.vec <- rep(22, length(leg.txt))
  c.vec <- c1
  par(mar = c(0, 0, 0, 0))
  plot(c(0, 0), c(1, 1), xlim = c(0, 1), ylim = c(0, 1), axes = F, type = "n", xlab = "", ylab = "")
  legend(
    x = 0, y = 0.8, legend = leg.txt, bty = "n", xjust = 0, yjust = 1, pch = p.vec, pt.bg = c.vec,
    col = "black", cex = 1.25, pt.cex = 2.5
  )

  x11(width = 14, height = 9)
  nf <- layout(matrix(c(1, 2, 3, 4, 5, 6), 3, 2, byrow = T), widths = c(5, 1), heights = c(1, 1, 1), respect = FALSE)
  par(mar = c(4, 7, 4, 7))
  MSIG.Score.Plot(
    z = score, main = paste(model.name, " Model Score (train)"), phen.cmap = c1,
    char.rescale = 1, col.classes = class.phen, col.labels = class.labels,
    xlab = "Samples", ylab = "Score", create.window = F, create.legend = T
  )

  par(mar = c(4, 7, 4, 7))
  norm.score <- (score - min(score)) / (max(score) - min(score))
  MSIG.Score.Plot(
    z = norm.score, main = paste(model.name, " Normalized Model Score (train)"), phen.cmap = c1,
    char.rescale = 1, col.classes = class.phen, col.labels = class.labels,
    xlab = "Samples", ylab = "Normalized Score", create.window = F, create.legend = T
  )

  par(mar = c(4, 7, 4, 7))
  MSIG.Score.Plot(
    z = prob.i[, 1], main = paste(model.name, " Probabiliy (train)"), phen.cmap = c1, char.rescale = 1,
    col.classes = class.phen, col.labels = class.labels,
    xlab = "Samples", ylab = "Probability", create.window = F, create.legend = T
  )


  MSIG.HeatMapPlot.6(
    V = rbind(score, model.score, t(prob.i[, 1])), row.names = c("raw.score", "model.score", "probability"),
    row.names2 = c(model.name, model.name, model.name),
    col.labels = class.labels2, col.labels2 = class.labels2,
    col.classes = class.phen2, phen.cmap = c1, phen.names = model.name,
    col.names = sample.names2, main = model.name,
    xlab = "  ", ylab = "  ", sub = "   ", row.norm = T, cmap.type = 3,
    char.rescale = 1, legend = T
  )

  # Save parameters in model file for annotation and to be used by apply function

  model.creation.date <- date()
  OPAM.write.param.line(param = model.creation.date, param.name = "model.creation.date", file = models.file, append = F)
  OPAM.write.param.line(param = input.ds, param.name = "input.ds", file = models.file, append = T)
  OPAM.write.param.line(param = input.ds, param.name = "input.ds", file = models.file, append = T)
  OPAM.write.param.line(param = input.cls, param.name = "input.cls", file = models.file, append = T)
  OPAM.write.param.line(param = input2.ds, param.name = "input2.ds", file = models.file, append = T)
  OPAM.write.param.line(param = input2.cls, param.name = "input2.cls", file = models.file, append = T)
  OPAM.write.param.line(param = target.class, param.name = "target.class", file = models.file, append = T)
  OPAM.write.param.line(param = target.class2, param.name = "target.class2", file = models.file, append = T)
  OPAM.write.param.line(param = model.name, param.name = "model.name", file = models.file, append = T)
  OPAM.write.param.line(param = model.description, param.name = "model.description", file = models.file, append = T)
  OPAM.write.param.line(param = sample.norm.type, param.name = "sample.norm.type", file = models.file, append = T)
  OPAM.write.param.line(param = marker.disc, param.name = "marker.disc", file = models.file, append = T)
  OPAM.write.param.line(param = top.markers.up, param.name = "top.markers.up", file = models.file, append = T)
  OPAM.write.param.line(param = top.markers.dn, param.name = "top.markers.dn", file = models.file, append = T)
  OPAM.write.param.line(param = top.markers.up2, param.name = "top.markers.up2", file = models.file, append = T)
  OPAM.write.param.line(param = top.markers.dn2, param.name = "top.markers.dn2", file = models.file, append = T)
  OPAM.write.param.line(param = statistic, param.name = "statistic", file = models.file, append = T)
  OPAM.write.param.line(param = weight, param.name = "weight", file = models.file, append = T)
  OPAM.write.param.line(param = random.seed, param.name = "random.seed", file = models.file, append = T)
  OPAM.write.param.line(param = nperm, param.name = "nperm", file = models.file, append = T)
  OPAM.write.param.line(param = link.function, param.name = "link.function", file = models.file, append = T)
  OPAM.write.param.line(param = c1, param.name = "c1", file = models.file, append = T)
  OPAM.write.param.line(param = msig.cntrl.genes, param.name = "msig.cntrl.genes", file = models.file, append = T)
  OPAM.write.param.line(param = msig.up.genes, param.name = "msig.up.genes", file = models.file, append = T)
  OPAM.write.param.line(param = msig.dn.genes, param.name = "msig.dn.genes", file = models.file, append = T)
  OPAM.write.param.line(param = msig.up.genes2, param.name = "msig.up.genes2", file = models.file, append = T)
  OPAM.write.param.line(param = msig.dn.genes2, param.name = "msig.dn.genes2", file = models.file, append = T)
  OPAM.write.param.line(param = msig.up.genes3, param.name = "msig.up.genes3", file = models.file, append = T)
  OPAM.write.param.line(param = msig.dn.genes3, param.name = "msig.dn.genes3", file = models.file, append = T)
  OPAM.write.param.line(param = beta.0, param.name = "beta.0", file = models.file, append = T)
  OPAM.write.param.line(param = beta.1, param.name = "beta.1", file = models.file, append = T)
} # end of function



##

OPAM.find.split <- function(msig, status) {
  index <- order(msig, decreasing = F)
  msig.s <- msig[index]
  status.s <- status[index]
  error.1 <- cumsum(status.s == 1)
  error.0 <- rev(cumsum(rev(status.s) == 0))
  min.i <- which.min(abs(error.1 - error.0))
  thres <- msig.s[min.i]
  # x11()
  # plot(msig.s, error.1, col=2, ylim = c(0, length(msig)))
  # points(msig.s, error.0, col=3)
  # points(c(thres, thres), c(0, length(msig)), col=1, type="l")

  return(thres)
}

OPAM.find.split.2 <- function(msig, status, method = "lowest.error") {
  if (method == "equal.errors") {
    index <- order(msig, decreasing = F)
    msig.s <- msig[index]
    status.s <- status[index]
    error.1 <- cumsum(status.s == 1)
    error.0 <- rev(cumsum(rev(status.s) == 0))
    min.i <- which.min(abs(error.1 - error.0))
    thres <- msig.s[min.i]
  } else if (method == "lowest.error") {
    index <- order(msig, decreasing = F)
    msig.s <- msig[index]
    status.s <- status[index]
    error.1 <- cumsum(status.s == 1)
    error.0 <- rev(cumsum(rev(status.s) == 0))
    min.i <- which.min(error.1 + error.0)
    thres <- msig.s[min.i]
  } else if (method == "logistic") {
    thres <- 0.5
  }

  return(thres)
}


OPAM.project.dataset <- function(
    input.ds,
    output.ds,
    gene.set.databases,
    gene.set.selection  = "ALL", # "ALL" or list with names of gene sets
    sample.norm.type    = "rank", # "rank", "log" or "log.rank"
    output.score.type   = "ES") # "ES" or "NES"
{ #----------------------------------------------------------------------------------------

  # Load libraries
  library(gtools)
  library(verification)
  library(RColorBrewer)

  # Read input dataset

  dataset <- MSIG.Gct2Frame(filename = input.ds) # Read gene expression dataset (GCT format)
  m <- data.matrix(dataset$ds)
  gene.names <- dataset$row.names
  gene.descs <- dataset$descs
  sample.names <- dataset$names
  Ns <- length(m[1, ])
  Ng <- length(m[, 1])
  temp <- strsplit(input.ds, split = "/") # Extract input file name
  s <- length(temp[[1]])
  input.file.name <- temp[[1]][s]
  temp <- strsplit(input.file.name, split = ".gct")
  input.file.prefix <- temp[[1]][1]

  # Sample normalization

  if (sample.norm.type == "rank") {
    for (j in 1:Ns) { # column rank normalization
      m[, j] <- rank(m[, j], ties.method = "average")
    }
    m <- 10000 * m / Ng
  } else if (sample.norm.type == "log.rank") {
    for (j in 1:Ns) { # column rank normalization
      m[, j] <- rank(m[, j], ties.method = "average")
    }
    m <- log(10000 * m / Ng + exp(1))
  } else if (sample.norm.type == "log") {
    m[m < 1] <- 1
    m <- log(m + exp(1))
  }

  # Read gene set databases

  max.G <- 0
  max.N <- 0
  for (gsdb in gene.set.databases) {
    GSDB <- Read.GeneSets.db(gsdb, thres.min = 2, thres.max = 2000, gene.names = NULL)
    max.G <- max(max.G, max(GSDB$size.G))
    max.N <- max.N + GSDB$N.gs
  }
  N.gs <- 0
  gs <- matrix("null", nrow = max.N, ncol = max.G)
  gs.names <- vector(length = max.N, mode = "character")
  gs.desc <- vector(length = max.N, mode = "character")
  size.G <- vector(length = max.N, mode = "numeric")
  start <- 1
  for (gsdb in gene.set.databases) {
    GSDB <- Read.GeneSets.db(gsdb, thres.min = 2, thres.max = 2000, gene.names = NULL)
    N.gs <- GSDB$N.gs
    gs.names[start:(start + N.gs - 1)] <- GSDB$gs.names
    gs.desc[start:(start + N.gs - 1)] <- GSDB$gs.desc
    size.G[start:(start + N.gs - 1)] <- GSDB$size.G
    gs[start:(start + N.gs - 1), 1:max(GSDB$size.G)] <- GSDB$gs[1:N.gs, 1:max(GSDB$size.G)]
    start <- start + N.gs
  }
  N.gs <- max.N

  # Select desired gene sets

  if (gene.set.selection[1] != "ALL") {
    locs <- match(gene.set.selection, gs.names)
    N.gs <- sum(!is.na(locs))
    gs <- gs[locs, ]
    gs.names <- gs.names[locs]
    gs.desc <- gs.desc[locs]
    size.G <- size.G[locs]
  }

  # Loop over gene sets

  score.matrix <- matrix(0, nrow = N.gs, ncol = Ns)
  for (gs.i in 1:N.gs) {
    gene.set <- gs[gs.i, 1:size.G[gs.i]]
    gene.overlap <- intersect(gene.set, gene.names)
    print(paste(gs.i, "gene set:", gs.names[gs.i], " overlap=", length(gene.overlap)))
    if (length(gene.overlap) == 0) {
      score.matrix[gs.i, ] <- runif(Ns, min = 1E-06, max = 1.1E-06)
      next
    } else {
      gene.set.locs <- match(gene.overlap, gene.set)
      gene.names.locs <- match(gene.overlap, gene.names)
      msig <- m[gene.names.locs, ]
      msig.names <- gene.names[gene.names.locs]
      if (output.score.type == "ES") {
        OPAM <- OPAM.Projection(
          data.array = m, gene.names = gene.names, n.cols = Ns,
          n.rows = Ng, weight = 0.25, statistic = "area.under.RES",
          gene.set = gene.overlap, nperm = 1
        )
        score.matrix[gs.i, ] <- OPAM$ES.vector
      } else if (output.score.type == "NES") {
        OPAM <- OPAM.Projection(
          data.array = m, gene.names = gene.names, n.cols = Ns,
          n.rows = Ng, weight = 0.25, statistic = "area.under.RES",
          gene.set = gene.overlap, nperm = 100
        )
        score.matrix[gs.i, ] <- OPAM$NES.vector
      }
    }
  }

  V.GCT <- data.frame(score.matrix)
  names(V.GCT) <- sample.names
  row.names(V.GCT) <- gs.names
  write.gct(gct.data.frame = V.GCT, descs = gs.desc, filename = output.ds)
} # end of OPAM.project.dataset



OPAM.sort.projection.by.score.2 <- function(
    input.ds,
    input.cls,
    results.dir,
    normalize.score = T,
    normalization.type = "zero.one",
    model,
    target.phen = NA,
    target.class = NA,
    user.colors = NA,
    decreasing.order = T,
    output.dataset = NA,
    char.rescale = 1,
    cmap.type = 3,
    row.norm = T) {
  library(gtools)
  library(verification)
  library(ROCR)
  library(MASS)
  library(RColorBrewer)
  library(heatmap.plus)

  dataset <- MSIG.Gct2Frame(filename = input.ds) # Read gene expression dataset (GCT format)
  m <- data.matrix(dataset$ds)
  model.names <- dataset$row.names
  model.descs <- dataset$descs
  Ns <- length(m[1, ])
  dim(m)
  sample.names <- dataset$names

  n.models <- length(m[, 1])
  temp <- strsplit(input.ds, split = "/") # Extract test file name
  s <- length(temp[[1]])
  test.file.name <- temp[[1]][s]
  temp <- strsplit(test.file.name, split = ".gct")
  test.file.prefix <- temp[[1]][1]
  char.res <- 0.013 * n.models + 0.65

  # normalize scores

  if (normalize.score == T) {
    if (normalization.type == "zero.one") {
      for (i in 1:n.models) {
        m[i, ] <- (m[i, ] - min(m[i, ])) / (max(m[i, ]) - min(m[i, ]))
      }
    } else if (normalization.type == "z.score") {
      for (i in 1:n.models) {
        m[i, ] <- (m[i, ] - mean(m[i, ])) / sd(m[i, ])
      }
    } else if (normalization.type == "r.z.score") {
      for (i in 1:n.models) {
        m[i, ] <- (m[i, ] - median(m[i, ])) / mad(m[i, ])
      }
    }
  }

  CLS <- MSIG.ReadPhenFile(file = input.cls) # Read phenotype file (CLS format)
  cls.labels <- CLS$class.v
  cls.phen <- CLS$phen
  cls.list <- CLS$class.list

  if (is.vector(cls.labels)) {
    n.phen <- 1
  } else {
    n.phen <- length(cls.labels[, 1])
  }
  if (!is.na(user.colors)) {
    c.test <- user.colors
  } else {
    if (!is.null(CLS$col.phen)) {
      c.test <- CLS$col.phen
    } else {
      c.test <- c(
        brewer.pal(n = 7, name = "Set2"), brewer.pal(n = 7, name = "Dark2"), brewer.pal(n = 7, name = "Set1"),
        brewer.pal(n = 8, name = "Accent"), brewer.pal(n = 9, name = "Spectral"), brewer.pal(n = 8, name = "Set3"),
        brewer.pal(n = 8, name = "BuGn"),
        brewer.pal(n = 7, name = "Set2"), brewer.pal(n = 7, name = "Dark2"), brewer.pal(n = 7, name = "Set1"),
        brewer.pal(n = 8, name = "Accent"), brewer.pal(n = 10, name = "Spectral"), brewer.pal(n = 8, name = "Set3"),
        brewer.pal(n = 8, name = "BuGn"),
        brewer.pal(n = 7, name = "Set2"), brewer.pal(n = 7, name = "Dark2"), brewer.pal(n = 7, name = "Set1"),
        brewer.pal(n = 8, name = "Accent"), brewer.pal(n = 10, name = "Spectral"), brewer.pal(n = 8, name = "Set3"),
        brewer.pal(n = 8, name = "BuGn")
      )
    }
  }


  if (!is.null(CLS$phen.names)) {
    phen.names <- CLS$phen.names
  } else {
    phen.names <- "NA"
  }

  cls.phen.index <- unlist(cls.phen)
  cls.phen.colors <- c.test[1:length(cls.phen.index)]

  n.classes <- vector(length = n.phen, mode = "numeric")
  if (n.phen == 1) {
    max.classes <- length(cls.phen)
    n.classes[1] <- max.classes
  } else {
    max.classes <- max(unlist(lapply(cls.phen, FUN = length)))
    for (i in 1:n.phen) {
      n.classes[i] <- length(cls.phen[[i]])
    }
  }

  filename <- paste(results.dir, test.file.prefix, ".SORT.PROJ", sep = "")
  #    pdf(file=paste(filename, ".pdf", sep=""), height = 11, width = 9)
  pdf(file = paste(filename, ".pdf", sep = ""), height = 8.5, width = 11)
  #   x11(width=12, height=8)

  loc <- match(model, model.names)
  print(c("loc:", loc))
  s.order <- order(m[loc, ], decreasing = decreasing.order)
  m2 <- m[, s.order]

  sample.names2 <- sample.names[s.order]

  if (is.vector(cls.labels)) {
    cls.labels2 <- cls.labels[s.order]
    cls.list2 <- cls.list[s.order]
  } else {
    cls.labels2 <- cls.labels[, s.order]
    cls.list2 <- cls.list[, s.order]
  }
  # Recompute cls.phen and cls.labels2 as order may have changed

  cls.phen2 <- NULL
  if (is.vector(cls.labels)) {
    classes <- unique(cls.list2)
    cls.phen2 <- classes
    cls.labels2 <- match(cls.list2, cls.phen2)
  } else {
    for (kk in 1:length(cls.list2[, 1])) {
      classes <- unique(cls.list2[kk, ])
      cls.phen2[[kk]] <- classes
      cls.labels2[kk, ] <- match(cls.list2[kk, ], cls.phen2[[kk]])
    }
  }


  correl <- cor(t(m2))[, loc]
  m.order <- order(correl, decreasing = T)
  correl2 <- correl[m.order]
  m2 <- m2[m.order, ]
  model.names2 <- model.names[m.order]
  model.descs2 <- paste(model.descs[m.order], signif(correl2, digits = 3))
  phen.list <- unlist(cls.phen2)
  colors.list <- cls.phen.colors[match(phen.list, cls.phen.index)]

  if (!is.na(target.phen)) {
    bin.class <- ifelse(cls.list2[target.phen, ] == target.class, 1, 0)
  } else {
    bin.class <- ifelse(cls.list2[1, ] == cls.list2[1, 1], 1, 0)
  }
  for (i in 1:n.models) {
    m.score <- m2[i, ]
    m.score.norm <- (m.score - min(m.score)) / (max(m.score) - min(m.score))
    perf.auc <- roc.area(bin.class, m.score.norm)
    print(paste("ROC=", signif(perf.auc$A, digits = 3), " p-val=", signif(perf.auc$p.value, digits = 3)))
    model.descs2[i] <- paste(signif(perf.auc$A, digits = 3), " (", signif(perf.auc$p.value, digits = 3), ")")
  }

  MSIG.HeatMapPlot.7(
    V = m2, row.names = model.names2,
    row.names2 = model.descs2, col.labels = cls.labels2,
    col.classes = cls.phen2, phen.cmap = colors.list, phen.names = phen.names,
    col.names = sample.names2, main = " ", xlab = "  ", ylab = "  ", row.norm = row.norm,
    cmap.type = cmap.type, char.rescale = char.rescale, legend = T
  )

  dev.off()

  if (!is.na(output.dataset)) {
    V.GCT <- m2
    colnames(V.GCT) <- sample.names2
    row.names(V.GCT) <- model.names2
    write.gct(gct.data.frame = V.GCT, descs = model.descs2, filename = output.dataset)
  }
}

plot.cont.model <- function(
    feature,
    feature.mean,
    feature.sd,
    # feature.min,
    # feature.max,
    num.points.extrapolation = 500,
    x,
    prob.i,
    beta0,
    beta1,
    col.vec,
    color.map = c("red", "green"),
    target.feature) {
  # plot
  range.x <- max(x) - min(x)
  prob.m <- matrix(0, nrow = num.points.extrapolation, ncol = 3)
  x.m <- vector(length = num.points.extrapolation, mode = "numeric")

  for (j in 1:num.points.extrapolation) {
    x.m[j] <- min(x) + j * (range.x / 500)
    p.vec <- exp(beta0 + beta1 * x.m[j]) / (1 + exp(beta0 + beta1 * x.m[j])) # Logit
    prob.m[j, ] <- quantile(p.vec, probs = c(0.05, 0.50, 0.95))
  }
  istar <- min(x) + which.min(abs(0.5 - prob.m[, 2])) * (range.x / num.points.extrapolation)

  x.index <- order(x, decreasing = F)
  x.order <- x[x.index]
  prob.i.order <- prob.i[x.index, ]
  col.vec.order <- col.vec[x.index]
  col.vec.order <- ifelse(col.vec.order == 1, color.map[1], color.map[2])
  #          x11(height = 7, width = 9.5)
  nf <- layout(matrix(c(1, 2), 1, 2, byrow = T), widths = c(3.5, 1), heights = 1, respect = FALSE)
  plot(x.order, prob.i.order[, 2],
    main = feature,
    #     sub=paste("Boundary: ", signif(istar, 3), "(", signif(istar*(feature.max - feature.min) + feature.min, 3), ")"),
    #      sub=paste("Boundary: ", signif(istar, 3), "(", signif(istar*feature.sd + feature.mean, 3), ")"),
    sub = paste("Boundary: ", signif(istar, 3)),
    pch = 20, ylim = c(-0.3, 1.07), col = 0, cex = 2, xlab = "Activation score", ylab = "Probability"
  )
  for (h in 1:length(x.m)) points(c(x.m[h], x.m[h]), c(prob.m[h, 1], prob.m[h, 3]), type = "l", col = "gray90", lty = 1, lwd = 2)
  points(x.m, prob.m[, 2], type = "l", lwd = 2, col = 1, lty = 1, cex = 1)
  points(x.m, prob.m[, 1], type = "l", col = 1, lty = 1, cex = 1)
  points(x.m, prob.m[, 3], type = "l", col = 1, lty = 1, cex = 1)
  range.x <- range(x.order)
  points(range.x, c(0.5, 0.5), type = "l", lty = 3, col = "gray", lwd = 2)
  points(range.x, c(0, 0), type = "l", lty = 1, col = 1, lwd = 1)
  points(c(istar, istar), c(-0.3, 1.07), type = "l", lty = 3, col = "gray", lwd = 2)
  points(x.order, prob.i.order[, 2], pch = 21, bg = col.vec.order, col = 1, cex = 1.5)
  red.points <- x.order[col.vec.order == color.map[1]]
  green.points <- x.order[col.vec.order == color.map[2]]
  points(range.x, c(-.1, -.1), type = "l", lty = 1, col = 1, lwd = 1)
  points(range.x, c(-.2, -.2), type = "l", lty = 1, col = 1, lwd = 1)
  points(red.points, rep(-0.1, length(red.points)), pch = 21, bg = color.map[1], col = 1, cex = 1.5)
  points(green.points, rep(-0.2, length(green.points)), pch = 21, bg = color.map[2], col = 1, cex = 1.5)
  leg.txt <- c(target.class, other.class)
  p.vec <- rep(21, 21)
  c.vec <- color.map
  old.par <- par(no.readonly = TRUE)
  par(mar = c(1, 1, 1, 1))
  plot(c(0, 0), c(1, 1), xlim = c(0, 1), ylim = c(0, 1), axes = F, type = "n", xlab = "", ylab = "")
  text(x = 0.5, y = 0.8, labels = target.feature, cex = 1.3)
  legend(x = 0, y = 0.8, legend = leg.txt, bty = "n", xjust = 0, yjust = 1, pch = p.vec, pt.bg = c.vec, col = "black", cex = 1.2, pt.cex = 1.2)
  par(old.par)
  return(istar)
}

OPAM.apply.model.3 <- function(
    input.ds,
    input.cls = NA,
    models.dir,
    models = "ALL",
    plots.outfile,
    cmap = NA,
    raw.score.outfile,
    norm.score.outfile,
    model.score.outfile,
    prob.outfile,
    gmt.file = NULL,
    graphics.off = F) { #----------------------------------------------------------------------------------------

  # Load libraries
  erf <- function(x) 2 * pnorm(x * sqrt(2)) - 1
  library(gtools)
  library(verification)
  library(ROCR)
  library(MASS)
  library(RColorBrewer)
  library(MCMCpack)

  pdf(file = plots.outfile, height = 8.5, width = 11)

  # Read test dataset
  dataset <- MSIG.Gct2Frame(filename = input.ds) # Read gene expression dataset (GCT format)
  m <- data.matrix(dataset$ds)
  gene.names <- dataset$row.names
  gene.descs <- dataset$descs
  sample.names <- dataset$names
  Ns <- length(m[1, ])
  Ng <- length(m[, 1])
  temp <- strsplit(input.ds, split = "/") # Extract test file name
  s <- length(temp[[1]])
  test.file.name <- temp[[1]][s]
  temp <- strsplit(test.file.name, split = ".gct")
  test.file.prefix <- temp[[1]][1]

  # Test set color map
  if (is.na(cmap)) cmap <- c(brewer.pal(n = 7, name = "Set2"), brewer.pal(n = 7, name = "Dark2"), brewer.pal(n = 7, name = "Set1"))

  if (!is.na(input.cls)) { # Read phenotype if CLS file is provided
    CLS <- MSIG.ReadClsFile(file = input.cls) # Read phenotype file (CLS format)
    class.labels <- CLS$class.v
    class.phen <- CLS$phen
    class.list <- CLS$class.list
  } else {
    class.labels <- rep(1, Ns)
    class.phen <- "UNDEFINED_PHEN"
    class.list <- rep("U", Ns)
  }

  # Loop over models

  if (models[[1]] == "ALL") { # use all models in directory
    file.list <- list.files(models.dir)
    models <- file.list[regexpr(pattern = ".mod", file.list) > 1]
    for (k.model in 1:length(models)) {
      temp <- strsplit(models[k.model], ".mod")
      models[k.model] <- temp[[1]]
    }
    models <- unique(models)
  }

  n.models <- length(models)
  score.matrix <- matrix(0, nrow = n.models, ncol = Ns)
  norm.score.matrix <- matrix(0, nrow = n.models, ncol = Ns)
  model.score.matrix <- matrix(0, nrow = n.models, ncol = Ns)
  probability.matrix <- matrix(0, nrow = n.models, ncol = Ns)
  models.descs <- NULL

  for (model.i in 1:n.models) {
    print(paste(model.i, "File:", models[model.i]))
    # Read parameters from model file
    m.file <- paste(models.dir, models[model.i], ".mod", sep = "")
    con <- file(m.file, "r")
    file.content <- readLines(con, n = -1)
    close(con)
    gc()
    len <- length(file.content)
    for (i in 1:len) {
      temp <- unlist(strsplit(file.content[[i]], "\t"))
      len.param <- length(temp)
      if (len.param == 2) {
        param.vals <- temp[2]
      } else {
        param.vals <- paste(noquote(temp[2:len.param]), collapse = ",")
        param.vals <- paste("c(", param.vals, ")", sep = "")
      }
      assignment.string <- paste(noquote(temp[1]), " <- ", param.vals, sep = "")
      #         print(c("executing assigment:", assignment.string))
      eval(parse(text = assignment.string))
    }

    print(paste("Model:", model.i, " model name:", model.name))


    # Set parameters
    if (!exists("random.seed")) random.seed <- 12345
    set.seed(random.seed)

    # Sample normalization

    if (!exists("sample.norm.type")) sample.norm.type <- "rank"

    if (sample.norm.type == "rank") {
      for (j in 1:Ns) { # column rank normalization
        m[, j] <- rank(m[, j], ties.method = "average")
      }
      m <- 10000 * m / Ng
    } else if (sample.norm.type == "log.rank") {
      for (j in 1:Ns) { # column rank normalization
        m[, j] <- rank(m[, j], ties.method = "average")
      }
      m <- log(10000 * m / Ng + exp(1))
    } else if (sample.norm.type == "log") {
      m[m < 1] <- 1
      m <- log(m + exp(1))
    }

    # Control signature normalization
    if (!exists("msig.cntrl.genes")) msig.cntrl.genes <- NA
    if (!is.na(msig.cntrl.genes)) {
      gene.names.int <- intersect(msig.cntrl.genes, gene.names)
      locs <- match(gene.names.int, gene.names, nomatch = 0)
      msig.cntrl <- m[locs, ]
      msig.cntrl.genes <- gene.names[locs]
      msig.cntrl.descs <- gene.descs[locs]
      msig.cntrl.size <- length(locs)
      if (msig.cntrl.size < 1) {
        msig.cntrl.center <- rep(1, Ns)
      } else if (msig.cntrl.size == 1) {
        msig.cntrl.center <- msig.cntrl
      } else if (msig.cntrl.size > 1) msig.cntrl.center <- apply(msig.cntrl, MARGIN = 2, FUN = mean)
      for (i in 1:Ng) {
        m[i, ] <- m[i, ] / msig.cntrl.center
      }
    }

    # Obtain UP & DN signatures

    if (exists("msig.up.genes3")) msig.up.genes <- msig.up.genes3
    gene.names.int <- intersect(msig.up.genes, gene.names)
    if (length(gene.names.int) < 2) {
      score.matrix[model.i, ] <- norm.score.matrix[model.i, ] <- model.score.matrix[model.i, ] <- probability.matrix[model.i, ] <- rep(0, Ns)
      models.descs <- c(models.descs, model.description)
      rm(
        model.creation.date, input.ds, input.cls, input2.ds, input2.cls, target.class,
        target.class2, model.name, model.description, sample.norm.type, marker.disc,
        top.markers.up, top.markers.dn, top.markers.up2, top.markers.dn2, statistic, weight,
        random.seed, nperm, link.function, c1, msig.cntrl.genes, msig.up.genes, msig.dn.genes,
        msig.up.genes2, msig.dn.genes2, msig.up.genes3, msig.dn.genes3, beta.0, beta.1, score,
        score.up, score.dn
      )
      next
    }

    if (!is.null(gmt.file)) { # save genes in a gmt file
      m.name <- paste(model.name, "_UP", sep = "")
      genes.string <- paste(msig.up.genes, sep = "\t", collapse = "\t")
      output.line <- paste(m.name, model.description, genes.string, sep = "\t", collapse = "")
      if (model.i == 1) {
        write(noquote(output.line), file = gmt.file, append = F, ncolumns = length(msig.up.genes) + 2)
      } else {
        write(noquote(output.line), file = gmt.file, append = T, ncolumns = length(msig.up.genes) + 2)
      }
      if (exists("msig.dn.genes")) {
        m.name <- paste(model.name, "_DN", sep = "")
        genes.string <- paste(msig.dn.genes, sep = "\t", collapse = "\t")
        output.line <- paste(m.name, model.description, genes.string, sep = "\t", collapse = "")
        write(noquote(output.line), file = gmt.file, append = T, ncolumns = length(msig.up.genes) + 2)
      }
    }


    locs <- match(gene.names.int, gene.names, nomatch = 0)
    msig.up.test <- m[locs, ]
    msig.up.genes.test <- gene.names[locs]
    msig.up.descs.test <- gene.descs[locs]
    msig.up.size.test <- length(locs)
    if (graphics.off == F) {
      #         x11(height = 9, width = 12)
      MSIG.HeatMapPlot.3(
        V = msig.up.test, row.names = msig.up.genes.test, col.labels = class.labels,
        col.classes = class.phen, phen.cmap = cmap, col.names = sample.names,
        main = paste(model.name, " UP signature test"), xlab = " ", ylab = " ", sub = " ",
        row.norm = T, cmap.type = 4, char.rescale = 1
      )
    }
    # Project test dataset
    OPAM <- OPAM.Projection(m, gene.names, Ns, Ng, weight, statistic, msig.up.genes.test, nperm = nperm)
    score.up <- OPAM$ES.vector

    if (exists("msig.dn.genes3")) msig.dn.genes <- msig.dn.genes3
    if (exists("msig.dn.genes")) {
      gene.names.int <- intersect(msig.dn.genes, gene.names)
      if (length(gene.names.int) < 2) {
        score.matrix[model.i, ] <- norm.score.matrix[model.i, ] <- model.score.matrix[model.i, ] <- probability.matrix[model.i, ] <- rep(0, Ns)
        models.descs <- c(models.descs, model.description)
        rm(
          model.creation.date, input.ds, input.cls, input2.ds, input2.cls, target.class,
          target.class2, model.name, model.description, sample.norm.type, marker.disc,
          top.markers.up, top.markers.dn, top.markers.up2, top.markers.dn2, statistic, weight,
          random.seed, nperm, link.function, c1, msig.cntrl.genes, msig.up.genes, msig.dn.genes,
          msig.up.genes2, msig.dn.genes2, msig.up.genes3, msig.dn.genes3, beta.0, beta.1, score,
          score.up, score.dn
        )
        next
      }
      locs <- match(gene.names.int, gene.names, nomatch = 0)
      msig.dn.test <- m[locs, ]
      msig.dn.genes.test <- gene.names[locs]
      msig.dn.descs.test <- gene.descs[locs]
      msig.dn.size.test <- length(locs)
      if (graphics.off == F) {
        #            x11(height = 9, width = 12)
        MSIG.HeatMapPlot.3(
          V = msig.dn.test, row.names = msig.dn.genes.test, col.labels = class.labels,
          col.classes = class.phen, phen.cmap = cmap, col.names = sample.names,
          main = paste(model.name, " DN signature test"), xlab = " ", ylab = " ", sub = " ",
          row.norm = T, cmap.type = 4, char.rescale = 1
        )
      }
      OPAM <- OPAM.Projection(m, gene.names, Ns, Ng, weight, statistic, msig.dn.genes.test, nperm = nperm)
      score.dn <- OPAM$ES.vector
    }

    if (!is.na(msig.cntrl.genes)) {
      if (graphics.off == F) {
        #            x11(height = 9, width = 12)
        MSIG.HeatMapPlot.3(
          V = msig.cntrl.test, row.names = msig.cntrl.genes.test, col.labels = class.labels,
          col.classes = class.phen, phen.cmap = cmap, col.names = sample.names,
          main = paste(model.name, " CNTRL signature"), xlab = " ", ylab = " ", sub = " ", row.norm = T,
          cmap.type = 4, char.rescale = 1
        )
      }
    }

    if (exists("msig.dn.genes")) {
      score <- score.up - score.dn
    } else {
      score <- score.up
    }

    if (graphics.off == F) {
      #         x11(width=14, height=9)
      if (exists("msig.dn.genes")) {
        nf <- layout(matrix(c(1, 2, 3, 0, 4, 0), 3, 2, byrow = T), widths = c(5, 1), heights = c(1, 1, 1), respect = FALSE)
      } else {
        nf <- layout(matrix(c(1, 2), 1, 2, byrow = T), widths = c(5, 1), heights = 1, respect = FALSE)
      }

      par(mar = c(2, 4, 2, 4))
      barplot(score.up,
        main = paste(model.name, " OPAM Score UP (test)", sep = ""), font.axis = 1.25, cex.lab = 1.5,
        cex.axis = 1.25, cex.names = 1.25, width = 1, space = 0, col = cmap[class.labels]
      )
      leg.txt <- class.phen
      p.vec <- rep(22, length(leg.txt))
      par(mar = c(0, 0, 0, 0))
      plot(c(0, 0), c(1, 1), xlim = c(0, 1), ylim = c(0, 1), axes = F, type = "n", xlab = "", ylab = "")
      legend(
        x = 0, y = 0.8, legend = leg.txt, bty = "n", xjust = 0, yjust = 1, pch = p.vec, pt.bg = cmap, col = "black",
        cex = 1.25, pt.cex = 2.5
      )

      if (exists("msig.dn.genes")) {
        par(mar = c(2, 4, 2, 4))
        barplot(score.dn,
          main = paste(model.name, " OPAM Score DOWN (test)", sep = ""), font.axis = 1.25, cex.lab = 1.5,
          cex.axis = 1.25, cex.names = 1.25, width = 1, space = 0, col = cmap[class.labels]
        )

        par(mar = c(2, 4, 2, 4))
        barplot(score,
          main = paste(model.name, " OPAM Total Score (test)", sep = ""), font.axis = 1.25, cex.lab = 1.5,
          cex.axis = 1.25, cex.names = 1.25, width = 1, space = 0, col = cmap[class.labels]
        )
      }
    }

    if (!exists("beta.0")) beta.0 <- 0
    if (!exists("beta.1")) beta.1 <- 1
    if (!exists("link.function")) link.function <- "logit"

    # Apply MCMC logit or probit model to test dataset

    model.formula <- "beta.0 + beta.1 * score[i]"
    model.formula
    prob.i <- matrix(0, nrow = Ns, ncol = 3)
    model.score <- vector(length = Ns, mode = "numeric")
    for (i in 1:Ns) {
      model.score[i] <- eval(parse(text = model.formula))
      if (link.function == "logit") {
        p.vec <- paste("inv.logit(x=", model.formula, ")", sep = "")
      } else if (link.function == "probit") {
        p.vec <- paste("(erf(", model.formula, ") + 1)/2", sep = "")
      } else {
        stop("Unknown link function")
      }
      val <- eval(parse(text = p.vec))
      prob.i[i, 1] <- quantile(val, probs = 0.5)
      prob.i[i, 2] <- quantile(val, probs = 0.05)
      prob.i[i, 3] <- quantile(val, probs = 0.95)
    }
    probability <- prob.i[, 1]
    xmin <- min(model.score)
    xmax <- max(model.score)
    range.x <- xmax - xmin
    n.points <- 1000
    prob.m <- matrix(0, nrow = n.points, ncol = 3)
    x.m <- vector(length = n.points, mode = "numeric")
    for (k in 1:n.points) {
      x.m[k] <- xmin + k * (range.x / n.points)
      if (link.function == "logit") {
        p.vec <- paste("inv.logit(x=", x.m[k], ")", sep = "")
      } else if (link.function == "probit") {
        p.vec <- paste("(erf(", x.m[k], ") + 1)/2", sep = "")
      } else {
        stop("Unknown link function")
      }
      val <- eval(parse(text = p.vec))
      prob.m[k, 1] <- quantile(val, probs = 0.5)
      prob.m[k, 2] <- quantile(val, probs = 0.05)
      prob.m[k, 3] <- quantile(val, probs = 0.95)
    }
    istar <- which.min(abs(0.5 - prob.m[, 1]))
    istar <- xmin + istar * (range.x / 1000)
    x.index <- order(model.score, decreasing = F)
    x.order <- model.score[x.index]
    prob.i.order <- prob.i[x.index, ]
    target.var.order <- cmap[class.labels[x.index]]
    class.labels.order <- class.labels[x.index]

    # Plot bar graph of z-scores
    boundary <- istar
    pred.class <- ifelse(prob.i.order[, 1] >= 0.5, 2, 1)
    z.range <- range(x.order)
    norm.score <- (score - min(score)) / (max(score) - min(score))
    if (graphics.off == F) {
      #         x11(height = 7, width = 9.5)
      nf <- layout(matrix(c(1, 2), 1, 2, byrow = T), widths = c(3.75, 1), heights = 1, respect = FALSE)
      plot(x.order, prob.i.order[, 1], sub = model.name, pch = 20, col = 0, cex = 2, xlab = "Activation Index", ylab = "Probability")
      points(x.m, prob.m[, 1], type = "l", lwd = 2, col = 1, lty = 1, cex = 1)
      points(x.m, prob.m[, 2], type = "l", col = 4, lty = 1, cex = 1)
      points(x.m, prob.m[, 3], type = "l", col = 4, lty = 1, cex = 1)
      arrows(x.order, prob.i.order[, 2], x.order, prob.i.order[, 3], col = 4, angle = 90, code = 3, length = 0.0)
      range.x <- range(x.order)
      points(range.x, c(0.5, 0.5), type = "l", lty = 3, col = 1, lwd = 2)
      points(range.x, c(-.15, -0.15), type = "l", lty = 1, col = 1, lwd = 2)
      points(c(istar, istar), c(-0.07, 1.07), type = "l", lty = 3, col = 1, lwd = 2)
      points(x.order, prob.i.order[, 1], pch = 21, bg = target.var.order, col = 1, cex = 2)
      points(x.order, rep(-0.15, length(x.order)), pch = 21, bg = target.var.order, col = 1, cex = 2)
      leg.txt <- class.phen
      p.vec <- rep(22, length(leg.txt))
      c.vec <- cmap
      par(mar = c(0, 0, 0, 0))
      plot(c(0, 0), c(1, 1), xlim = c(0, 1), ylim = c(0, 1), axes = F, type = "n", xlab = "", ylab = "")
      legend(
        x = 0, y = 0.8, legend = leg.txt, bty = "n", xjust = 0, yjust = 1, pch = p.vec, pt.bg = c.vec,
        col = "black", cex = 1.25, pt.cex = 2.5
      )

      #         x11(width = 14, height = 9)
      nf <- layout(matrix(c(1, 2, 3, 4, 5, 6), 3, 2, byrow = T), widths = c(5, 1), heights = c(1, 1, 1), respect = FALSE)
      par(mar = c(4, 7, 4, 7))
      MSIG.Score.Plot(
        z = score, main = paste(model.name, " Model Score (test)"), phen.cmap = cmap,
        char.rescale = 1, col.classes = class.phen, col.labels = class.labels,
        xlab = "Samples", ylab = "Score", create.window = F, create.legend = T
      )

      par(mar = c(4, 7, 4, 7))

      MSIG.Score.Plot(
        z = norm.score, main = paste(model.name, " Normalized Model Score (test)"), phen.cmap = cmap,
        char.rescale = 1, col.classes = class.phen, col.labels = class.labels,
        xlab = "Samples", ylab = "Normalized Score", create.window = F, create.legend = T
      )

      par(mar = c(4, 7, 4, 7))
      MSIG.Score.Plot(
        z = prob.i[, 1], main = paste(model.name, " Probabiliy (test)"), phen.cmap = cmap, char.rescale = 1,
        col.classes = class.phen, col.labels = class.labels,
        xlab = "Samples", ylab = "Probability", create.window = F, create.legend = T
      )

      #         x11(width = 14, height = 9)
      MSIG.HeatMapPlot.6(
        V = rbind(score, norm.score, model.score, probability), row.names = c(
          "raw.score", "norm.score",
          "model.score", "probability"
        ),
        row.names2 = c(model.name, model.name, model.name, model.name),
        col.labels = class.labels, col.labels2 = class.labels,
        col.classes = class.phen, phen.cmap = cmap, phen.names = model.name,
        col.names = sample.names, main = model.name,
        xlab = "  ", ylab = "  ", sub = "   ", row.norm = T, cmap.type = 3,
        char.rescale = 1, legend = T
      )
    }
    score.matrix[model.i, ] <- score
    norm.score.matrix[model.i, ] <- norm.score
    model.score.matrix[model.i, ] <- model.score
    probability.matrix[model.i, ] <- probability
    models.descs <- c(models.descs, model.description)
    rm(
      model.creation.date, input.ds, input.cls, input2.ds, input2.cls, target.class,
      target.class2, model.name, model.description, sample.norm.type, marker.disc,
      top.markers.up, top.markers.dn, top.markers.up2, top.markers.dn2, statistic, weight,
      random.seed, nperm, link.function, c1, msig.cntrl.genes, msig.up.genes, msig.dn.genes,
      msig.up.genes2, msig.dn.genes2, msig.up.genes3, msig.dn.genes3, beta.0, beta.1, score,
      score.up, score.dn
    )

    if (graphics.off == F) {
      if (model.i %% 5 == 0) {
        graphics.off()
        gc()
      }
    }
  } # end of loop over models

  # Plot projections for all models

  #   x11(width = 14, height = 9)
  MSIG.HeatMapPlot.6(
    V = score.matrix, row.names = models, row.names2 = models.descs, col.labels = class.labels,
    col.labels2 = class.labels, col.classes = class.phen, phen.cmap = cmap,
    col.names = sample.names, main = paste(test.file.name, " - Raw Scores"), xlab = " ", ylab = " ",
    sub = "Raw Scores", row.norm = T, cmap.type = 3, char.rescale = 1, legend = T
  )

  #   x11(width = 14, height = 9)
  MSIG.HeatMapPlot.6(
    V = norm.score.matrix, row.names = models, row.names2 = models.descs, col.labels = class.labels,
    col.labels2 = class.labels, col.classes = class.phen, phen.cmap = cmap,
    col.names = sample.names, main = paste(test.file.name, " - Norm Scores"), xlab = " ", ylab = " ",
    sub = "Norm Scores", row.norm = T, cmap.type = 3, char.rescale = 1, legend = T
  )

  #   x11(width = 14, height = 9)
  MSIG.HeatMapPlot.6(
    V = model.score.matrix, row.names = models, row.names2 = models.descs, col.labels = class.labels,
    col.labels2 = class.labels, col.classes = class.phen, phen.cmap = cmap,
    col.names = sample.names, main = paste(test.file.name, " - Model Scores"), xlab = " ", ylab = " ",
    sub = "Model Scores", row.norm = T, cmap.type = 3, char.rescale = 1, legend = T
  )

  #   x11(width = 14, height = 9)
  MSIG.HeatMapPlot.6(
    V = probability.matrix, row.names = models, row.names2 = models.descs, col.labels = class.labels,
    col.labels2 = class.labels, col.classes = class.phen, phen.cmap = cmap,
    col.names = sample.names, main = paste(test.file.name, " - Probabilities"), xlab = " ", ylab = " ",
    sub = "Probabilities", row.norm = T, cmap.type = 3, char.rescale = 1, legend = T
  )

  # Save projections in files

  V.GCT <- data.frame(score.matrix)
  names(V.GCT) <- sample.names
  row.names(V.GCT) <- models
  write.gct(gct.data.frame = V.GCT, descs = models.descs, filename = raw.score.outfile)

  V.GCT <- data.frame(norm.score.matrix)
  names(V.GCT) <- sample.names
  row.names(V.GCT) <- models
  write.gct(gct.data.frame = V.GCT, descs = models.descs, filename = norm.score.outfile)

  V.GCT <- data.frame(model.score.matrix)
  names(V.GCT) <- sample.names
  row.names(V.GCT) <- models
  write.gct(gct.data.frame = V.GCT, descs = models.descs, filename = model.score.outfile)

  V.GCT <- data.frame(probability.matrix)
  names(V.GCT) <- sample.names
  row.names(V.GCT) <- models
  write.gct(gct.data.frame = V.GCT, descs = models.descs, filename = prob.outfile)

  dev.off()
} # end of function


runSsGSEAwithPermutation <- function(profile_data_file, MOD_file, number_perms,
                                     select.mods) {
  # NEW ADD 1/3
  profile_data_path <- dirname(file.path(profile_data_file))
  profile_data_name <- basename(profile_data_file)
  MOD_file <- MOD_file

  # please keep the first and second column as the same with NAME and Description
  data <- read.table(paste(profile_data_path, profile_data_name, sep = "/"), header = T, row.names = 1, skip = 2) # read data from a gct file and ignore first 2 lines and 1 column
  size <- dim(data[-1]) # size of data
  num_rows <- size[1] # number of columns
  num_cols <- size[2] # number of rows
  num_perm <- number_perms # number of permutations
  colnames <- colnames(data) # g et the column names

  print("Data loading finished!")

  line1 <- "#1.2"
  line2 <- paste(num_rows, num_perm, sep = "\t")
  line3 <- paste("X", 1:num_perm, sep = "", collapse = "\t")
  line3 <- paste("NAME", "Description", line3, sep = "\t")
  write.table(line1, paste(profile_data_path, "/random_profile_", profile_data_name, ".gct", sep = ""), quote = F, col.name = F, row.name = F) # first line for the gct file
  write.table(line2, paste(profile_data_path, "/random_profile_", profile_data_name, ".gct", sep = ""), quote = F, col.name = F, row.name = F, append = T) # 2nd line for the gct file
  write.table(line3, paste(profile_data_path, "/random_profile_", profile_data_name, ".gct", sep = ""), quote = F, col.name = F, row.name = F, append = T) # 2nd line for the gct file


  random_profile <- data[1]
  for (i in 1:num_perm) {
    a <- data[-1][cbind(seq(1:num_rows), sample(1:num_cols, num_rows, replace = T))] # the first column is the decription column, ignore it.
    random_profile <- cbind(random_profile, a)
    if (i %% 100 == 0) {
      print(i)
    }
  }
  write.table(random_profile, paste(profile_data_path, "/random_profile_", profile_data_name, ".gct", sep = ""), sep = "\t", quote = F, col.name = F, append = T)


  print("Random profiles was genereated!")

  # source(msig_library_file)
  selected.models <- select.mods
  OPAM.apply.model.2(
    input.ds = paste(profile_data_path, "/random_profile_", profile_data_name, ".gct", sep = ""),
    models.dir = MOD_file,
    models = selected.models,
    raw.score.outfile = paste(profile_data_path, "/random_raw.score_", profile_data_name, ".gct", sep = ""),
    norm.score.outfile = "", model.score.outfile = "", prob.outfile = "", graphics.off = T
  )

  print("SsGSEA was performed on random profiles!")

  OPAM.apply.model.2(
    input.ds = paste(profile_data_path, profile_data_name, sep = "/"),
    models.dir = MOD_file,
    models = selected.models,
    raw.score.outfile = paste(profile_data_path, "/raw.score_", profile_data_name, ".gct", sep = ""),
    norm.score.outfile = "", model.score.outfile = "", prob.outfile = "", graphics.off = T
  )

  print("SsGSEA was performed on the original profiles!")

  random_result <- read.table(paste(profile_data_path, "/random_raw.score_", profile_data_name, ".gct", sep = ""), header = T, row.names = 1, skip = 2, sep = "\t")
  random_result <- random_result[-1]
  random_result <- t(random_result)
  original_result <- read.table(paste(profile_data_path, "/raw.score_", profile_data_name, ".gct", sep = ""), header = T, row.names = 1, skip = 2, sep = "\t")
  original_result <- original_result[-1]
  original_result <- t(original_result)


  p_result <- original_result
  for (i in 1:dim(original_result)[1]) {
    p_result[i, ] <- colSums(sweep(random_result, 2, original_result[i, ]) >= 0)
  }


  # NEW ADD 2/3
  file.remove(paste(profile_data_path, "/random_profile_", profile_data_name, ".gct", sep = ""))
  file.remove(paste(profile_data_path, "/random_raw.score_", profile_data_name, ".gct", sep = ""))
  file.remove(paste(profile_data_path, "/raw.score_", profile_data_name, ".gct", sep = ""))


  # NEW ADD 3/3
  colnames(p_result) <- paste(colnames(p_result), "pval", sep = "_")
  p_result <- t(apply(p_result, 1, function(x) {
    (x + 1) / (number_perms + 1)
  }))


  write.table(cbind(original_result, p_result), paste(profile_data_path, "/p_result_", profile_data_name, ".txt", sep = ""), sep = "\t", quote = F)

  print("P_values for each subtype was calculated!")

  print("Finished!")
}
