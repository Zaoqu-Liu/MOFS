rm(list = ls())

# Load required libraries
# The libraries loaded here are fundamental for different types of data handling and analysis:
# - tidyverse: For efficient data manipulation and visualization.
# - clusterProfiler: For enrichment analysis involving pathways and gene sets.
library(tidyverse)
library(clusterProfiler)

# Prepare gene sets from GMT files
# These gene sets are essential for pathway analysis, which will help identify biological pathways enriched in different gene sets.
# c2: Canonical pathways, c5: GO (Gene Ontology), h: Hallmark pathways
gene_set_files <- c('c2.cp.v2022.1.Hs.symbols.gmt',
                    'c5.go.v2022.1.Hs.symbols.gmt',
                    'h.all.v2022.1.Hs.symbols.gmt')

gene_sets <- lapply(gene_set_files, function(file) {
  genes <- clusterProfiler::read.gmt(file)  # Read GMT files to obtain gene sets
  genes$term <- as.character(genes$term)  # Convert gene set terms to character type
  return(genes)
}) %>% Reduce(rbind, .)  # Combine all gene sets into a single dataframe

# Load multi-omics data
# The multi-omics data contains expression values for mRNA and protein, allowing us to integrate different types of molecular information.
load('multi-omics-intergration.rda')

# Filter gene sets to include only genes present in the mRNA and protein datasets
# This ensures that we focus only on genes with available data, improving the relevance and reliability of the analysis.
gene_sets_mRNA <- gene_sets[gene_sets$gene %in% rownames(mm$Log2FPKM), ]
gene_sets_protein <- gene_sets[gene_sets$gene %in% rownames(mm$Protein), ]

# Split gene sets by term and filter to retain those with at least 15 genes
# Gene sets with fewer genes may not provide robust insights, so we only keep larger sets for better statistical power.
gene_sets_mRNA_list <- split(gene_sets_mRNA$gene, gene_sets_mRNA$term)
mRNA_gene_count <- sapply(gene_sets_mRNA_list, length)
gene_sets_mRNA_filtered <- gene_sets_mRNA_list[names(mRNA_gene_count)[mRNA_gene_count >= 15]]

gene_sets_protein_list <- split(gene_sets_protein$gene, gene_sets_protein$term)
protein_gene_count <- sapply(gene_sets_protein_list, length)
gene_sets_protein_filtered <- gene_sets_protein_list[names(protein_gene_count)[protein_gene_count >= 15]]

# Single-sample pathway activity analysis using modified Wilcoxon gene set test
# This function calculates pathway activity scores for each sample using a modified Wilcoxon gene set test (ssMwwGST).
# It uses parallel computation to speed up the processing, especially when dealing with large datasets.
ssMwwGST <- function(expression_data, gene_sets, n_cores = 8) {
  library(yaGST)  # For performing gene set testing
  library(doMC)   # For parallel computation to speed up the analysis
  
  # Standardize the gene expression data by calculating Z-scores for each gene
  means <- rowMeans(expression_data)
  sds <- apply(expression_data, 1, sd)
  
  registerDoMC(n_cores)  # Register the number of cores for parallel processing
  results <- foreach(sample_idx = 1:ncol(expression_data)) %dopar% {
    current_sample <- (expression_data[, sample_idx] - means) / sds
    ranked_list <- sort(current_sample, decreasing = TRUE)  # Rank genes for each sample
    
    # Perform the modified Wilcoxon test for each gene set
    pathway_results <- lapply(gene_sets, function(genes) {
      mwwGST(ranked_list, geneSet = genes, minLenGeneSet = 15, alternative = "two.sided", verbose = FALSE)
    })
    pathway_results <- pathway_results[sapply(pathway_results, length) != 0]  # Filter out empty results
    
    # Extract NES and p-values from the test results
    NES <- sapply(pathway_results, function(x) x$log.pu)
    p_values <- sapply(pathway_results, function(x) x$p.value)
    
    list(NES = NES, p_values = p_values)
  }
  
  # Compile the results into matrices for further analysis
  NES_matrix <- sapply(results, function(x) x$NES)
  p_value_matrix <- sapply(results, function(x) x$p_values)
  colnames(NES_matrix) <- colnames(p_value_matrix) <- colnames(expression_data)
  
  # Adjust p-values for multiple testing using the FDR method
  FDR_matrix <- t(apply(p_value_matrix, 1, function(x) p.adjust(x, method = "fdr")))
  
  list(NES = NES_matrix, p_values = p_value_matrix, FDR = FDR_matrix)
}

# Run single-sample pathway activity analysis on mRNA and protein data
# The analysis identifies the activity of different biological pathways for each sample based on the expression data.
res_ssMwwGST_mRNA <- ssMwwGST(mm$Log2FPKM, gene_sets_mRNA_filtered, n_cores = 15)
save(res_ssMwwGST_mRNA, file = 'Transcript-pathway-ssMwwGST.rda')

res_ssMwwGST_protein <- ssMwwGST(mm$Protein, gene_sets_protein_filtered, n_cores = 15)
save(res_ssMwwGST_protein, file = 'Protein-pathway-ssMwwGST.rda')

# Load pre-computed clustering results
# The clustering results are used to understand how samples group together based on their multi-omics profiles.
load('pre-trail-results.rda')
Cluster_info <- allres$coca_res$Cluster
Cluster_info <- Cluster_info[Cluster_info$Sample %in% allres$Core_set, ]
Cluster_info$Cluster <- paste0('C', Cluster_info$Cluster)  # Rename clusters for easier reference

# Differential Enrichment Analysis of Pathways (DEAP)
# This function identifies differentially enriched pathways for each cluster by comparing pathway activity between samples in a given cluster and those not in that cluster.
lzq_DEAPath <- function(Cluster_data, ssMwwGST_results, dea_FDR_threshold = 0.001, dea_gap_threshold = 1.5) {
  colnames(Cluster_data)[1] <- 'Sample'
  unique_clusters <- unique(Cluster_data$Cluster)  # Extract the unique cluster labels
  
  # Extract the NES data from the single-sample pathway analysis results
  NES_data <- as.data.frame(ssMwwGST_results$NES)
  NES_data <- NES_data[, Cluster_data$Sample]  # Keep only the samples that are present in the Cluster_data
  
  # Perform differential pathway analysis for each cluster
  pathway_results <- map(unique_clusters, function(cluster) {
    cluster_samples <- Cluster_data$Sample[Cluster_data$Cluster == cluster]  # Samples in the current cluster
    other_samples <- Cluster_data$Sample[Cluster_data$Cluster != cluster]  # Samples not in the current cluster
    
    res <- apply(NES_data, 1, function(values) {
      # Perform Wilcoxon rank-sum test to compare pathway activity between cluster and other samples
      fit <- wilcox.test(values[cluster_samples], values[other_samples])
      median_diff <- median(values[cluster_samples]) - median(values[other_samples])
      mean_diff <- mean(values[cluster_samples]) - mean(values[other_samples])
      
      # Store the results for each pathway
      data.frame(Cluster = cluster, P = fit$p.value, Statistic = as.numeric(fit$statistic),
                 Median_gap = median_diff, Mean_gap = mean_diff)
    }) %>% Reduce(rbind, .)
    
    # Adjust p-values for multiple testing
    res$FDR <- p.adjust(res$P, method = 'BH')
    res$Pathway <- rownames(NES_data)
    # Clean pathway names for easier interpretation
    res$Pathway_clean <- gsub('_', ' ', tolower(gsub('([A-Z]{1,})_.*', '', res$Pathway)))
    res$Pathway_clean <- Hmisc::capitalize(res$Pathway_clean)
    res <- dplyr::select(res, Cluster, Pathway_clean, FDR, Median_gap, Mean_gap, everything())
    return(res)
  })
  
  # Filter pathways based on threshold for median differences and FDR
  filtered_results <- map(pathway_results, function(res) {
    res[res$Median_gap >= dea_gap_threshold & res$FDR < dea_FDR_threshold, ]
  })
  
  list(all_pathways = pathway_results, significant_pathways = filtered_results, NES = NES_data, Cluster = Cluster_data)
}

# Run DEAP on mRNA and protein data
# This analysis identifies pathways that are significantly enriched in different clusters, helping to understand the biological differences between the clusters.
trans_dea_path <- lzq_DEAPath(Cluster_info, res_ssMwwGST_mRNA, dea_FDR_threshold = 0.001, dea_gap_threshold = 1.5)
pro_dea_path <- lzq_DEAPath(Cluster_info, res_ssMwwGST_protein, dea_FDR_threshold = 0.001, dea_gap_threshold = 1)
## Protein threshold is lowered because the number of protein measurements is significantly lower compared to mRNA

save(trans_dea_path, pro_dea_path, file = 'ssMwwGST-Pathway-DEA-ProTrans-results.rda')
load('ssMwwGST-Pathway-DEA-ProTrans-results.rda')

# Generate heatmap for differentially enriched pathways in mRNA and protein data
# The heatmaps help visualize pathway activity across samples in different clusters, highlighting patterns of enrichment.
create_heatmap <- function(dea_path_results, NES_data, Cluster_data, file_name) {
  pathway_ids <- lapply(dea_path_results$significant_pathways, function(x) x$Pathway_clean)  # Get significant pathway names
  pathway_ids_flat <- unlist(pathway_ids)  # Flatten the list of pathways
  NES_filtered <- NES_data[pathway_ids_flat, ]  # Subset NES data to only include significant pathways
  Cluster_data <- Cluster_data[order(Cluster_data$Cluster), ]  # Order samples by cluster
  NES_filtered <- NES_filtered[, Cluster_data$Sample]  # Reorder NES data by sample
  
  # Create the heatmap using ComplexHeatmap
  library(ComplexHeatmap)
  Heatmap(NES_filtered, show_row_names = FALSE, show_column_names = FALSE, name = 'NES',
          col = colorRamp2(c(-2, 0, 2), c('navy', 'white', 'firebrick3')),
          cluster_columns = FALSE, cluster_rows = FALSE,
          column_split = Cluster_data$Cluster,  # Split columns by cluster to visualize differences between clusters
          row_split = rep(LETTERS[1:3], times = sapply(pathway_ids, length)))
}

# Generate heatmaps for mRNA and protein data
generate_heatmap(trans_dea_path, trans_dea_path$NES, Cluster_info, 'mRNA-pathway-heatmap.pdf')
generate_heatmap(pro_dea_path, pro_dea_path$NES, Cluster_info, 'protein-pathway-heatmap.pdf')

# Differential expression analysis of individual genes
# This function compares the gene expression levels between samples in a specific cluster and those in other clusters.
lzq_DEAgene <- function(data, Cluster_data) {
  unique_clusters <- unique(Cluster_data$Cluster)  # Extract unique cluster labels
  
  # For each cluster, calculate fold change in gene expression compared to other clusters
  results <- lapply(unique_clusters, function(cluster) {
    cluster_samples <- Cluster_data$Sample[Cluster_data$Cluster == cluster]  # Samples in the current cluster
    other_samples <- Cluster_data$Sample[Cluster_data$Cluster != cluster]  # Samples not in the current cluster
    
    res <- apply(data, 1, function(values) {
      fc <- mean(values[cluster_samples]) / mean(values[other_samples])  # Fold change calculation
    }) %>% as.data.frame() %>% tibble::rownames_to_column('ID') %>%
      rename('FC' = '.') %>% filter(FC != 0) %>%  # Remove entries with no fold change
      mutate(logFC = log2(FC)) %>% arrange(desc(logFC))  # Calculate log fold change and sort by descending logFC
    return(res)
  })
  
  return(results)
}

# Load gene expression data for RNA-seq and protein data
# Prepares RNA-seq and protein data for differential expression analysis.
load('RNA_seq_FPKM.rda')
fpkm <- fpkm %>% group_by(geneID) %>% summarise_all(max) %>% column_to_rownames('geneID')
fpkm <- fpkm[, allres$Core_set]  # Keep only core samples
fpkm <- fpkm[apply(fpkm, 1, function(x) sum(x == 0) <= 0.7 * ncol(fpkm)), ]  # Filter genes with too many zero values
ann <- data.table::fread('geneInfo.tab', data.table = FALSE, header = FALSE)
ann <- ann[ann$V2 %in% rownames(fpkm), ]
fpkm <- fpkm[ann$V2[ann$V3 == 'protein_coding'], ] %>% na.omit()  # Keep only protein-coding genes

load('Protein.rda')

# Run differential expression analysis on RNA-seq and protein data
dea_gene_mRNA <- lzq_DEAgene(fpkm, Cluster_info)
dea_gene_protein <- lzq_DEAgene(Protein, Cluster_info)

# Run GSEA analysis on DE genes for RNA-seq and protein data
# Gene Set Enrichment Analysis (GSEA) is used to determine whether a set of genes shows statistically significant, concordant differences between two biological states.
run_GSEA <- function(dea_gene_results, gene_sets) {
  lapply(dea_gene_results, function(gene_data) {
    ranked_genes <- gene_data$logFC  # Rank genes based on their log fold change
    names(ranked_genes) <- gene_data$ID  # Assign gene IDs as names for ranking
    gsea_result <- GSEA(ranked_genes, pvalueCutoff = 1, TERM2GENE = gene_sets)  # Perform GSEA
    
    # Sort GSEA results by normalized enrichment score (NES)
    sorted_gsea <- gsea_result[order(gsea_result$NES, decreasing = TRUE), ]
    sorted_gsea$ID_clean <- gsub('_', ' ', tolower(gsub('([A-Z]{1,})_.*', '', sorted_gsea$ID)))  # Clean pathway names
    sorted_gsea$ID_clean <- Hmisc::capitalize(sorted_gsea$ID_clean)
    dplyr::select(sorted_gsea, ID_clean, everything())
    
    list(GSEA_project = gsea_result, GSEA_matrix = sorted_gsea)
  })
}

# Perform GSEA on RNA-seq and protein DE results
GSEA_results_mRNA <- run_GSEA(dea_gene_mRNA, ll_trans)
GSEA_results_protein <- run_GSEA(dea_gene_protein, ll_pro)

# Save GSEA results
save(GSEA_results_mRNA, GSEA_results_protein, file = 'GSEA-Pathway-ProTrans-results.rda')
