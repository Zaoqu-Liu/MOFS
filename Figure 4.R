# Load necessary libraries and set environment
rm(list = ls())
library(tidyverse)
library(ggsci)
library(ComplexHeatmap)
library(circlize)
library(dendextend)
library(seriation)
library(survival)
library(survminer)
library(patchwork)

# Define color palette for plots
# Here we define a consistent color palette that will be used across different plots for visual consistency.
# The colors are chosen to be distinct to easily differentiate between clusters in the plots.
colors <- c('#119da4', '#FF6666', '#ffc857')

# Load the data, specifically mutation summary for downstream analyses
# The Mutation_summary object is expected to contain a list of dataframes detailing various mutation statistics,
# including data for tumor mutation burden, SNP burden, INDEL burden, and other mutation-related metrics.
load('Mutation_summary.rda')

# MOFS Cluster Labels Visualization
# Create a small visual legend showing the cluster names with their respective colors.
# This helps in understanding the different subtypes (MOFS1, MOFS2, MOFS3) used throughout the analysis.
mofs_cluster_labels_plot <- ggplot() + 
  geom_label(data = data.frame(ID = paste0('MOFS', 1:3)),
             aes(label = ID, x = ID, fill = ID), 
             y = 0.5, color = "white", fontface = 'bold', 
             size = 12/.pt, hjust = 0.5, vjust = 0.5) + 
  scale_fill_manual(values = colors) + 
  theme_void() + 
  theme(legend.position = 'none', plot.margin = unit(c(0, 0.5, 0, 0.5), "cm"))

# Tumor Mutation Burden Analysis
# Compare the distribution of Tumor Mutation Burden (TMB) between different MOFS clusters.
# We use jitter plots to show individual data points, box plots to show median and variability,
# and violin plots to show the distribution density, providing a comprehensive visualization of TMB differences among MOFS clusters.
data <- Mutation_summary$Mutsummary
# Replace 'C' with 'MOFS' in cluster names to make them more interpretable.
data$Cluster <- gsub('C', 'MOFS', data$Cluster)
# Filter out samples with extreme TMB values (greater than 100) to focus on typical TMB ranges.
data <- data[data$TMB < 100, ]

plot_tmb <- ggplot(data, aes(Cluster, log2(TMB + 1))) +
  geom_jitter(shape = 21, size = 2.7, width = 0.2, aes(fill = Cluster, color = Cluster)) +
  geom_boxplot(outlier.colour = NA, aes(fill = Cluster), notch = TRUE, color = 'black', size = 0.6, alpha = 0.65) +
  geom_violin(alpha = 0.5, aes(fill = Cluster), color = NA, trim = TRUE) +
  # Add statistical comparison to check for significant differences between clusters.
  stat_compare_means(label.y = 5.2, label.x = 1, fontface = 'plain', size = 4.5) +
  scale_fill_manual(values = colors) +
  scale_color_manual(values = colors) +
  theme_bw(base_rect_size = 2.5) +
  labs(x = NULL, y = 'Tumor Mutation Burden') +
  theme(axis.text.y = element_text(size = 12, colour = 'black'),
        axis.text.x = element_blank(),
        axis.title.x = element_text(size = 15, colour = 'black', face = 'bold'),
        axis.title.y = element_text(size = 15, colour = 'black', face = 'bold'),
        panel.grid = element_blank(),
        panel.grid.major = element_line(color = "#cacfd2", linetype = "dashed"),
        panel.background = element_rect(fill = '#f3f6f6'),
        plot.title = element_text(hjust = 0.5, size = 15, colour = 'black', face = 'bold'),
        legend.position = 'none',
        axis.ticks.x = element_blank(),
        plot.margin = unit(c(0.1, 0.5, 0, 0.5), "cm"))

# Combine MOFS Cluster Labels with Tumor Mutation Burden Plot for Comprehensive Visualization
# The combined plot helps provide context to which MOFS cluster is represented by which color, making it easier to interpret the results.
mofs_cluster_labels_plot / plot_tmb + plot_layout(heights = c(1, 15))

# Single Nucleotide Polymorphism (SNP) Burden Analysis
# Compare the SNP burden between different MOFS clusters.
# Similar to the TMB analysis, jitter, box, and violin plots are used to provide a detailed comparison.
plot_snp_burden <- ggplot(data, aes(Cluster, log2(SNP + 1))) +
  geom_jitter(shape = 21, size = 2.7, width = 0.2, aes(fill = Cluster, color = Cluster)) +
  geom_boxplot(outlier.colour = NA, aes(fill = Cluster), notch = TRUE, color = 'black', size = 0.6, alpha = 0.65) +
  geom_violin(alpha = 0.5, aes(fill = Cluster), color = NA, trim = TRUE) +
  # Add statistical comparison to check for significant differences between clusters.
  stat_compare_means(label.y = 5.2, label.x = 1, fontface = 'plain', size = 4.5) +
  scale_fill_manual(values = colors) +
  scale_color_manual(values = colors) +
  theme_bw(base_rect_size = 2) +
  labs(x = NULL, y = 'SNP Burden') +
  theme(axis.text.y = element_text(size = 12, colour = 'black'),
        axis.text.x = element_blank(),
        axis.title.x = element_text(size = 15, colour = 'black', face = 'bold'),
        axis.title.y = element_text(size = 15, colour = 'black', face = 'bold'),
        panel.grid = element_blank(),
        panel.grid.major = element_line(color = "#cacfd2", linetype = "dashed"),
        panel.background = element_rect(fill = '#f3f6f6'),
        plot.title = element_text(hjust = 0.5, size = 15, colour = 'black', face = 'bold'),
        legend.position = 'none',
        axis.ticks.x = element_blank(),
        plot.margin = unit(c(0.1, 0.5, 0, 0.5), "cm"))

# Combine Cluster Labels with SNP Burden Plot
mofs_cluster_labels_plot / plot_snp_burden + plot_layout(heights = c(1, 15))

# INDEL Burden Analysis
# Compare the burden of Insertions and Deletions (INDELs) between different MOFS clusters.
# Jitter plots, box plots, and violin plots are used for a detailed visualization, allowing for the identification of potential differences in INDEL burden among clusters.
plot_indel_burden <- ggplot(data, aes(Cluster, log2(INDEL + 1))) +
  geom_jitter(shape = 21, size = 2.7, width = 0.2, aes(fill = Cluster, color = Cluster)) +
  geom_boxplot(outlier.colour = NA, aes(fill = Cluster), notch = FALSE, color = 'black', size = 0.6, alpha = 0.65) +
  geom_violin(alpha = 0.5, aes(fill = Cluster), color = NA, trim = TRUE) +
  # Add statistical comparison to check for significant differences between clusters.
  stat_compare_means(label.y = 0.8, label.x = 1, fontface = 'plain', size = 4.5) +
  scale_fill_manual(values = colors) +
  scale_color_manual(values = colors) +
  theme_bw(base_rect_size = 2) +
  labs(x = NULL, y = 'INDEL Burden') +
  theme(axis.text.y = element_text(size = 12, colour = 'black'),
        axis.text.x = element_blank(),
        axis.title.x = element_text(size = 15, colour = 'black', face = 'bold'),
        axis.title.y = element_text(size = 15, colour = 'black', face = 'bold'),
        panel.grid = element_blank(),
        panel.grid.major = element_line(color = "#cacfd2", linetype = "dashed"),
        panel.background = element_rect(fill = '#f3f6f6'),
        plot.title = element_text(hjust = 0.5, size = 15, colour = 'black', face = 'bold'),
        legend.position = 'none',
        axis.ticks.x = element_blank(),
        plot.margin = unit(c(0.1, 0.5, 0, 0.5), "cm"))

# Combine Cluster Labels with INDEL Burden Plot
mofs_cluster_labels_plot / plot_indel_burden + plot_layout(heights = c(1, 15))

# Custom function to sort mutation matrices for visualization
# This function sorts the input matrix based on mutation status, with non-mutated samples ('N') first and mutated samples ('Mutated') later for better visual distinction.
lzq_sortMutationMatrix <- function(dd){
  dd2 <- t(dd) %>% as.data.frame()
  for (i in rev(colnames(dd2))) {
    dd2 <- dd2[order(dd2[, i]), ]
  }
  dd <- t(dd2) %>% as.data.frame()
  return(dd)
}

# Load mutation data and process it for visualization
# The Mutation_summary contains several dataframes for different mutation analyses.
maf <- Mutation_summary$MAF
mutation_summary <- Mutation_summary$Mutsummary %>% arrange(Cluster) %>%
  mutate(Cluster = gsub('C', 'MOFS', Cluster))

# Extract mutation information for top 54 genes
# We reformat the matrix to indicate whether a gene is 'Mutated' or 'Non-Mutated'.
top_mutations <- Mutation_summary$Mut_top54[, mutation_summary$Sample]
top_mutations[top_mutations == 0] <- 'N'  # Represent non-mutated entries as 'N'
top_mutations[top_mutations != 'N'] <- 'Mutated'  # Represent mutated entries as 'Mutated'

# Calculate mutation frequency for each gene
# This helps identify genes that are frequently mutated in each MOFS cluster.
gene_mutation_frequency <- apply(top_mutations, 1, function(x) { sum(x == 'Mutated') / length(x) }) %>% as.numeric() %>% round(2)
names(gene_mutation_frequency) <- rownames(top_mutations)
gene_mutation_frequency <- sort(gene_mutation_frequency, decreasing = TRUE)
top_mutations <- top_mutations[names(gene_mutation_frequency),]

# Calculate mutation frequencies by cluster
# We calculate the mutation frequency for genes within each specific cluster, which helps in comparison across different MOFS clusters.
mutation_frequency_mofs1 <- apply(top_mutations[, 1:34], 1, function(x) { sum(x == 'Mutated') / length(x) }) %>% as.numeric() %>% round(2)
mutation_frequency_mofs2 <- apply(top_mutations[, 35:67], 1, function(x) { sum(x == 'Mutated') / length(x) }) %>% as.numeric() %>% round(2)
mutation_frequency_mofs3 <- apply(top_mutations[, 68:116], 1, function(x) { sum(x == 'Mutated') / length(x) }) %>% as.numeric() %>% round(2)

# Convert frequencies to percentages for better readability
mutation_frequency_mofs1 <- paste0(mutation_frequency_mofs1 * 100, '%')
mutation_frequency_mofs2 <- paste0(mutation_frequency_mofs2 * 100, '%')
mutation_frequency_mofs3 <- paste0(mutation_frequency_mofs3 * 100, '%')
mutation_frequency_mofs3 <- ifelse(str_length(mutation_frequency_mofs3) == 2, paste0(mutation_frequency_mofs3, '  '), mutation_frequency_mofs3)

# Annotate mutation results summary for heatmap
# We combine multiple statistics, such as differential P-values, Cox P-values, and hazard ratios to give a complete summary.
mutation_results_summary <- Mutation_summary$Mut_top54_summary$results
mutation_results_summary <- mutation_results_summary[match(names(gene_mutation_frequency), mutation_results_summary$Gene),]
mutation_summary_text <- paste0(
  mutation_frequency_mofs3,
  '       ',
  sprintf('%.2f', mutation_results_summary$Diff_P),
  '    ',
  sprintf('%.2f', mutation_results_summary$COX_P),
  '    ',
  sprintf('%.2f', mutation_results_summary$HR),
  ' [',
  sprintf('%.2f', mutation_results_summary$HR95L),
  '-',
  sprintf('%.2f', mutation_results_summary$HR95H),
  ']'
)

# Define significant genes for annotation purposes
# These include significant genes identified by survival analysis and differential mutation analysis.
survival_significant_genes <- c('MCM10', 'LRP2', 'TP53')
differential_significant_genes <- c('USH2A', 'SCN5A', 'PLEC', 'DNAH3')

# Set row colors for significant genes to differentiate them visually in the heatmap
row_colors <- rep('black', nrow(top_mutations))
row_colors[which(rownames(top_mutations) %in% survival_significant_genes)] <- '#0077b6'
row_colors[which(rownames(top_mutations) %in% differential_significant_genes)] <- '#e63946'

# Define row annotations for heatmap
annotation_mofs1 <- rowAnnotation(foo = anno_text(mutation_frequency_mofs1, gp = gpar(fontsize = 8, col = row_colors)))
annotation_mofs2 <- rowAnnotation(foo = anno_text(mutation_frequency_mofs2, gp = gpar(fontsize = 8, col = row_colors)))
annotation_mofs3 <- rowAnnotation(foo = anno_text(mutation_summary_text, gp = gpar(fontsize = 8, col = row_colors)))

# Define top annotations for each cluster
names(colors) <- paste0('MOFS', 1:3)
top_annotation_mofs1 <- HeatmapAnnotation(Cluster = mutation_summary$Cluster[1:34], border = FALSE, show_legend = FALSE, show_annotation_name = FALSE, simple_anno_size = unit(4, 'mm'), col = list(Cluster = colors))
top_annotation_mofs2 <- HeatmapAnnotation(Cluster = mutation_summary$Cluster[35:67], border = FALSE, show_legend = FALSE, show_annotation_name = FALSE, simple_anno_size = unit(4, 'mm'), col = list(Cluster = colors))
top_annotation_mofs3 <- HeatmapAnnotation(Cluster = mutation_summary$Cluster[68:116], border = FALSE, show_legend = FALSE, show_annotation_name = FALSE, simple_anno_size = unit(4, 'mm'), col = list(Cluster = colors))

# Generate heatmaps for each MOFS cluster and then combine them
# Sorting mutation matrices within each cluster allows for better visual separation of samples based on mutation status.
sorted_mutations_mofs1 <- lzq_sortMutationMatrix(top_mutations[, 1:34])
heatmap_mofs1 <- Heatmap(as.matrix(sorted_mutations_mofs1),
                         col = c('Mutated' = '#67AB9F', 'N' = 'WhiteSmoke'),
                         right_annotation = annotation_mofs1,
                         row_names_gp = gpar(fontsize = 9, fontface = 'italic', col = row_colors),
                         top_annotation = top_annotation_mofs1,
                         cluster_rows = FALSE, cluster_columns = TRUE, row_names_side = 'left',
                         show_column_names = FALSE, column_title = NULL,
                         show_heatmap_legend = FALSE)

sorted_mutations_mofs2 <- lzq_sortMutationMatrix(top_mutations[, 35:67])
heatmap_mofs2 <- Heatmap(as.matrix(sorted_mutations_mofs2),
                         col = c('Mutated' = '#67AB9F', 'N' = 'WhiteSmoke'),
                         right_annotation = annotation_mofs2,
                         row_names_gp = gpar(fontsize = 9, fontface = 'italic', col = row_colors),
                         top_annotation = top_annotation_mofs2,
                         cluster_rows = FALSE, cluster_columns = TRUE, row_names_side = 'left',
                         show_column_names = FALSE, show_row_names = FALSE, column_title = NULL,
                         show_heatmap_legend = FALSE)

sorted_mutations_mofs3 <- lzq_sortMutationMatrix(top_mutations[, 68:116])
heatmap_mofs3 <- Heatmap(as.matrix(sorted_mutations_mofs3),
                         col = c('Mutated' = '#67AB9F', 'N' = 'WhiteSmoke'),
                         right_annotation = annotation_mofs3,
                         row_names_gp = gpar(fontsize = 9, fontface = 'italic', col = row_colors),
                         top_annotation = top_annotation_mofs3,
                         cluster_rows = FALSE, cluster_columns = TRUE, row_names_side = 'left',
                         show_column_names = FALSE, show_row_names = FALSE, column_title = NULL,
                         show_heatmap_legend = FALSE)

# Combine heatmaps for all MOFS clusters for comprehensive visualization
combined_heatmap <- heatmap_mofs1 + heatmap_mofs2 + heatmap_mofs3
draw(combined_heatmap)

# Oncogenic Pathways Analysis
# Analyze the oncogenic pathways in different clusters to gain insights into potential pathways driving mutations.
c1_maf <- subsetMaf(maf, mutation_summary$Sample[mutation_summary$Cluster == 'MOFS1'])
c2_maf <- subsetMaf(maf, mutation_summary$Sample[mutation_summary$Cluster == 'MOFS2'])
c3_maf <- subsetMaf(maf, mutation_summary$Sample[mutation_summary$Cluster == 'MOFS3'])

# Generate pathway mutation plots for each MOFS cluster
pathway_mutations_mofs1 <- OncogenicPathways(maf = c1_maf, panelWidths = c(1, 4, 4))
pathway_mutations_mofs2 <- OncogenicPathways(maf = c2_maf, panelWidths = c(1, 4, 4))
pathway_mutations_mofs3 <- OncogenicPathways(maf = c3_maf, panelWidths = c(1, 4, 4))

# Summary
# The analysis explores the mutation rates and patterns in the top 54 genes across MOFS clusters. Detailed heatmaps are generated to visualize the distribution of mutations, and significant genes are highlighted to provide insights into their roles in different clusters. Additionally, oncogenic pathway analysis was conducted to understand potential signaling pathways involved in the tumorigenesis of each cluster, contributing to understanding molecular characteristics and potential therapeutic strategies.

# Clear workspace
rm(list = ls())

# Define color palette for plots
color_palette <- c('#119da4', '#FF6666', '#ffc857')

# Load CNV summary data
load('CNV_summary.rda')
cnv_data <- CNV_summary$CNVsummary

# Preprocess Cluster Information
# Replace cluster labels to a more descriptive label (e.g., MOFS1, MOFS2, MOFS3)
cnv_data$Cluster <- gsub('C', 'MOFS', cnv_data$Cluster)
cnv_data$Cluster <- factor(cnv_data$Cluster, levels = paste0('MOFS', 3:1))

# Plotting CNV Broad Burden Analysis
# This plot visualizes the CNV broad burden across different MOFS clusters.
# The use of jitter, boxplot, and violin plot aims to provide detailed insight into data distribution.
p_broad_burden <- ggplot(cnv_data, aes(x = B_Burden, y = Cluster)) +
  geom_jitter(shape = 21, size = 2.7, width = 0.2, aes(fill = Cluster, color = Cluster)) +
  geom_boxplot(outlier.colour = NA, aes(fill = Cluster), notch = TRUE, color = 'black', size = 0.6, alpha = 0.65) +
  geom_violin(alpha = 0.5, aes(fill = Cluster), color = NA, trim = TRUE) +
  stat_compare_means(
    comparisons = list(c('MOFS1', 'MOFS2'), c('MOFS1', 'MOFS3'), c('MOFS2', 'MOFS3')),
    label = 'p.signif', tip.length = 0, vjust = 0.8, size = 5
  ) +
  scale_fill_manual(values = color_palette) +
  scale_color_manual(values = color_palette) +
  theme_bw(base_rect_size = 2.5) +
  labs(y = NULL, x = 'CNV Broad Burden') +
  theme(
    axis.text.x = element_text(size = 12, colour = 'black'),
    axis.text.y = element_blank(),
    axis.title.y = element_text(size = 15, colour = 'black', face = 'bold'),
    axis.title.x = element_text(size = 15, colour = 'black', face = 'bold'),
    panel.grid = element_blank(),
    panel.grid.major = element_line(color = '#cacfd2', linetype = 'dashed'),
    panel.background = element_rect(fill = '#f3f6f6'),
    plot.title = element_text(hjust = 0.5, size = 15, colour = 'black', face = 'bold'),
    legend.position = 'none',
    axis.ticks.y = element_blank(),
    plot.margin = unit(c(0, 0, 0, 0.5), 'cm')
  )

# Plotting CNV Focal Burden Analysis
# Similar to the broad burden analysis, this plot represents the CNV focal burden across MOFS clusters.
p_focal_burden <- ggplot(cnv_data, aes(x = F_Burden, y = Cluster)) +
  geom_jitter(shape = 21, size = 2.7, width = 0.2, aes(fill = Cluster, color = Cluster)) +
  geom_boxplot(outlier.colour = NA, aes(fill = Cluster), notch = TRUE, color = 'black', size = 0.6, alpha = 0.65) +
  geom_violin(alpha = 0.5, aes(fill = Cluster), color = NA, trim = TRUE) +
  stat_compare_means(
    comparisons = list(c('MOFS1', 'MOFS2'), c('MOFS1', 'MOFS3'), c('MOFS2', 'MOFS3')),
    label = 'p.signif', tip.length = 0, vjust = 0.8, size = 5
  ) +
  scale_fill_manual(values = color_palette) +
  scale_color_manual(values = color_palette) +
  theme_bw(base_rect_size = 2.5) +
  labs(y = NULL, x = 'CNV Focal Burden') +
  theme(
    axis.text.x = element_text(size = 12, colour = 'black'),
    axis.text.y = element_blank(),
    axis.title.y = element_text(size = 15, colour = 'black', face = 'bold'),
    axis.title.x = element_text(size = 15, colour = 'black', face = 'bold'),
    panel.grid = element_blank(),
    panel.grid.major = element_line(color = '#cacfd2', linetype = 'dashed'),
    panel.background = element_rect(fill = '#f3f6f6'),
    plot.title = element_text(hjust = 0.5, size = 15, colour = 'black', face = 'bold'),
    legend.position = 'none',
    axis.ticks.y = element_blank(),
    plot.margin = unit(c(0, 0.5, 0, 0), 'cm')
  )

# Creating Label Plot for MOFS Clusters
# The label plot is used to visually distinguish MOFS clusters in subsequent combined plots.
label_plot <- ggplot() + 
  geom_label(
    data = data.frame(Cluster = factor(paste0('MOFS', 1:3), paste0('MOFS', 3:1))),
    aes(label = Cluster, y = Cluster, fill = Cluster), 
    x = 0.5, color = 'white', fontface = 'bold',
    size = 11 / .pt, hjust = 0.5, vjust = 0.5
  ) + 
  scale_fill_manual(values = rev(color_palette)) + 
  theme_void() + 
  theme(legend.position = 'none', plot.margin = unit(c(0, 0.04, 0, 0.1), 'cm'))

# Combine Plots for Overall Visualization
# Combining the broad and focal burden analyses with cluster labels to create a consolidated figure.
combined_plot <- p_broad_burden + label_plot + p_focal_burden + plot_layout(widths = c(5, 1, 5))

# Plot Heatmap Analysis for CNV Broad Burden
# A heatmap is used to visualize the CNV alterations across different clusters, providing insight into cluster-specific CNV changes.
cnv_data <- cnv_data[order(cnv_data$Cluster), ]
heatmap_annotation <- HeatmapAnnotation(
  Subtype = anno_block(
    gp = gpar(fill = color_palette, lwd = 2), 
    labels = paste0('MOFS', 1:3),
    height = unit(5, 'mm'),
    labels_gp = gpar(cex = 0.8, col = 'white', fontface = 'bold')
  )
)

# Draw Heatmap of CNV Broad Burden
# Using complex heatmap for a more detailed exploration of CNV features within each cluster.
cnv_broad_heatmap <- Heatmap(
  CNV_summary$CNV_broad, name = 'CNV',
  color_space = 'RGB', border = TRUE, border_gp = gpar(col = 'black', lwd = 2.5),
  col = colorRamp2(breaks = c(-1, 0, 1), colors = c('#0077b6', 'white', '#e63946')),
  cluster_rows = FALSE, cluster_columns = FALSE,
  column_split = cnv_data$Cluster,
  row_title = NULL, column_title = NULL,
  top_annotation = heatmap_annotation,
  heatmap_legend_param = list(border = TRUE)
)
draw(cnv_broad_heatmap, heatmap_legend_side = 'left')





