rm(list = ls())

# Load required libraries
library(tidyverse)   # Includes ggplot2, dplyr, and other tools for data manipulation and visualization
library(ggsci)       # Provides color palettes for scientific graphics
library(ComplexHeatmap)  # Useful for creating complex heatmaps
library(dendextend)  # To manipulate dendrograms
require(psych)       # Provides various statistical functions and tools
# Load the dataset which contains results for different combinations of cluster numbers and selection criteria
load('All-possible-selection-clusternum-final.rda')

# Prepare dataset ------------------------------------------------------------------------
# The dataset "res2" contains a variety of features, including metrics for different combinations of selected variables (e.g., mRNA, Protein, Pathway, MRI), along with clustering quality scores (CPI, GAP, and Score).
# The goal is to identify the optimal combination of features that provide the best clustering quality, as indicated by the 'Score' metric.
# Assign appropriate column names to improve readability and sort based on the 'Score' column, starting from the highest.
colnames(res2) <- c("mRNA_number", "Protein_number", "Path_number", "MRI_number", 'Cluster_number', 'CPI', 'GAP', 'Score')
res2 <- res2[order(res2$Score, decreasing = TRUE), ]  # Sort dataset by Score
rownames(res2) <- NULL  # Remove row names

# Filter the dataset to focus on the combination of features that produces the best clustering result
# This filtering ensures we examine and visualize the top-ranking combination of features
res <- res2[res2$mRNA_number == res2$mRNA_number[1] & 
              res2$Protein_number == res2$Protein_number[1] &
              res2$Path_number == res2$Path_number[1] &
              res2$MRI_number == res2$MRI_number[1], ]

# Transform the filtered data into a longer format for visualization
# This allows visualization of CPI and GAP metrics across the cluster numbers
res2 <- pivot_longer(res, cols = 6:7, names_to = 'nn', values_to = 'vv')

# Define colors for plotting
cols <- pal_npg()(10)[1:2]

# Create a line and scatter plot to compare CPI and GAP statistics across different cluster numbers
# CPI (Cluster Prediction Index) and GAP-statistics are commonly used metrics to assess the quality of clustering, indicating how well-separated clusters are.
p_cpi_gap <- ggplot(res2, aes(Cluster_number, vv)) +
  annotate("rect", xmin = 2.8, xmax = 3.2, ymin = 0.495, ymax = 0.63, alpha = .6, fill = '#EEDF7C') +
  geom_line(aes(group = nn), size = 0.8) +
  geom_point(aes(fill = nn), size = 4, shape = 21, stroke = 1.5) +
  scale_y_continuous(name = "Cluster prediction index", 
                     sec.axis = sec_axis(trans = ~., name = "GAP-statistics")) +
  scale_fill_manual(values = cols) +
  labs(x = 'Cluster numbers', y = NULL) +
  theme_bw(base_rect_size = 2) +
  theme(panel.grid = element_blank(),
        panel.grid.major = element_line(color = "#cacfd2", linetype = "dashed"),
        panel.background = element_rect(fill = '#f3f6f6'),
        axis.text.x = element_text(size = 12, colour = 'black'),
        axis.text.y = element_text(size = 12, colour = 'black'),
        axis.title.x = element_text(size = 15, colour = 'black', face = 'bold'),
        axis.title.y.left = element_text(size = 15, colour = cols[1], face = 'bold'),
        axis.title.y.right = element_text(size = 15, colour = cols[2], face = 'bold', angle = 90),
        plot.title = element_text(hjust = 0.5, size = 15, colour = 'black', face = 'bold'),
        legend.position = 'none',
        legend.background = element_blank(),
        legend.title = element_blank(),
        legend.key = element_blank())

# Generate consensus heatmap for clustering visualization ------------------------------------------------
# Consensus clustering is a widely used method for assessing the robustness of clusters by aggregating results from multiple clustering runs.
# Load previously computed clustering results and consensus matrix
load('pre-trail-results.rda')

# Extract clustering results and consensus matrix from the loaded data
consClust <- allres$coca_res$optimal$consensusClass
consMatrix <- allres$coca_res$optimal$consensusMatrix
colnames(consMatrix) <- names(consClust)
rownames(consMatrix) <- names(consClust)
consClust <- sort(consClust)
consMatrix <- consMatrix[names(consClust), names(consClust)]

# Prepare colors for the heatmap annotations
cols <- c('#119da4', '#FF6666', '#ffc857')

# Create a top annotation for the heatmap to indicate the identified subtypes
# The subtypes (e.g., MOFS1, MOFS2, MOFS3) represent different distinct clusters identified in the dataset
Top = HeatmapAnnotation(Subtype = paste0('MOFS', consClust),
                        border = FALSE,
                        show_legend = TRUE,
                        show_annotation_name = FALSE,
                        simple_anno_size = unit(5, 'mm'),
                        col = list(Subtype = c('MOFS1' = cols[1],
                                               'MOFS2' = cols[2],
                                               'MOFS3' = cols[3])))

# Create and color dendrogram based on hierarchical clustering
# The dendrogram visually represents how individual data points cluster together at different levels of similarity
row_dend <- as.dendrogram(hclust(dist(consMatrix)))
row_dend <- color_branches(row_dend, k = 3, col = rev(cols))

# Draw the consensus heatmap with customized options
Heatmap(consMatrix, name = ' ',
        border = TRUE,
        col = colorRampPalette(c("#f3f6f6", "#053061"))(51),
        show_row_names = FALSE,
        show_column_names = FALSE,
        cluster_rows = row_dend,
        cluster_columns = TRUE,
        top_annotation = Top,
        color_space = "RGB",
        show_column_dend = FALSE,
        show_row_dend = TRUE,
        row_dend_width = unit(5, 'mm'),
        row_dend_reorder = TRUE,
        row_dend_gp = gpar(lwd = 2),
        row_gap = unit(0, 'mm'),
        column_gap = unit(0, 'mm'),
        row_title = NULL, column_title = NULL)

# Principal Component Analysis (PCA) for visualization of sample separation ---------------------------------------
# PCA is a dimensionality reduction technique that helps visualize patterns and clusters in high-dimensional data.
# Here, the PCA is used to show how different samples separate based on the features included in the analysis.
jaccarddis_core <- allres$JaccardDis[allres$Core_set, allres$Core_set] %>% as.data.frame()
pc <- prcomp(jaccarddis_core, center = TRUE, scale. = TRUE)
f <- summary(pc)$importance
trg <- predict(pc, jaccarddis_core) %>% as.data.frame()
pp <- trg[, 1:2]
pp <- merge(pp, allres[["coca_res"]][["Cluster"]], by.x = 0, by.y = 1)
pp$Subtype <- paste0('MOFS', pp$Cluster)

# Plot PCA results using ggplot
p_pca <- ggplot(pp, aes(PC1, PC2, fill = Subtype)) +
  scale_fill_manual(values = cols) +
  scale_color_manual(values = cols) +
  stat_ellipse(aes(color = Subtype, fill = Subtype), geom = "polygon",
               alpha = 0.2, type = "norm", level = 0.9) +
  geom_point(size = 3, shape = 21, color = 'black') +
  labs(x = 'PC1 (50.9%)', y = 'PC2 (26.8%)') +
  theme_bw(base_rect_size = 2) +
  theme(panel.grid = element_blank(),
        panel.grid.major = element_line(color = "#cacfd2", linetype = "dashed"),
        panel.background = element_rect(fill = '#f3f6f6'),
        axis.text.x = element_text(size = 12, colour = 'black'),
        axis.text.y = element_text(size = 12, colour = 'black'),
        axis.title.x = element_text(size = 15, colour = 'black', face = 'bold'),
        axis.title.y = element_text(size = 15, colour = 'black', face = 'bold'),
        plot.title = element_text(hjust = 0.5, size = 15, colour = 'black', face = 'bold'),
        legend.position = 'none',
        legend.background = element_blank(),
        legend.title = element_blank(),
        legend.key = element_blank()) +
  annotate('text', x = -4.5, y = -8, label = 'MOFS3', fontface = 'bold.italic') +
  annotate('text', x = 8.2, y = -8, label = 'MOFS1', fontface = 'bold.italic') +
  annotate('text', x = 8.2, y = 10, label = 'MOFS2', fontface = 'bold.italic')

# Heatmap for clustering all samples for a comprehensive view of subtypes ------------------------------------------
# This heatmap aims to provide a detailed visual representation of all samples' clustering, which helps understand subtype distributions across the dataset.
hh <- allres$All_Cluster
cc <- allres[["coca_res"]][["Cluster"]]
cc$Subtype <- paste0('MOFS', cc$Cluster)
cc <- cc[order(cc$Cluster), ]
hh <- hh[, cc$Sample]

Top = HeatmapAnnotation(Subtype = cc$Subtype,
                        border = FALSE,
                        show_legend = TRUE,
                        show_annotation_name = FALSE,
                        simple_anno_size = unit(4, 'mm'),
                        col = list(Subtype = c('MOFS1' = cols[1],
                                               'MOFS2' = cols[2],
                                               'MOFS3' = cols[3])))

Heatmap(hh, name = ' ',
        border = TRUE,
        col = colorRampPalette(c("#f3f6f6", "#053061"))(10),
        show_row_names = TRUE,
        show_column_names = FALSE,
        cluster_rows = TRUE,
        cluster_columns = TRUE,
        top_annotation = Top,
        color_space = "RGB",
        show_column_dend = FALSE,
        show_row_dend = TRUE,
        row_dend_width = unit(5, 'mm'),
        row_dend_reorder = TRUE,
        row_dend_gp = gpar(lwd = 1.5),
        column_split = cc$Subtype,
        row_gap = unit(0, 'mm'), column_gap = unit(0, 'mm'),
        row_title = NULL, column_title = NULL)

# Analysis of Proportion of Ambiguous Clustering (PAC) -----------------------------------------
# PAC measures the proportion of samples whose cluster membership is uncertain.
PAC <- allres$coca_res$PAC
PAC$tt <- ifelse(PAC$PAC == min(PAC$PAC), 'A', 'B')

# Create plot for PAC analysis
p_pac <- ggplot(PAC, aes(factor(K), PAC, group = 1)) +
  geom_line(size = 0.8, color = 'white') +
  geom_point(size = 4, shape = 21, stroke = 1.5, aes(fill = tt), color = 'white') +
  scale_fill_npg() +
  labs(y = 'Proportion of ambiguous clustering', x = 'Cluster number K') +
  annotate('segment', x = 3.4, xend = 2.2, y = 0.195, yend = 0.185,
           size = 1.5, arrow = arrow(), alpha = 0.8, color = '#EEDF7C') +
  theme(axis.text = element_text(size = 12, colour = 'black'),
        axis.title = element_text(size = 14, colour = 'black'),
        plot.title = element_text(size = 15, hjust = 0.5),
        legend.position = 'none')

# Calinski and Harabasz index analysis ---------------------------------------------------------
# The Calinski and Harabasz index evaluates clustering quality, with higher values indicating better-defined clusters.
Calinsky <- data.frame(K = 2:6, Calinsky = allres$Calinsky[1:5])
Calinsky$tt <- ifelse(Calinsky$Calinsky == max(Calinsky$Calinsky), 'A', 'B')

# Create plot for Calinski and Harabasz index
p_calinsky <- ggplot(Calinsky, aes(factor(K), Calinsky, group = 1)) +
  geom_line(size = 0.8, color = 'white') +
  geom_point(size = 4, shape = 21, stroke = 1.5, aes(fill = tt), color = 'white') +
  scale_fill_npg() +
  labs(y = 'Calinski and Harabasz index', x = 'Cluster number K') +
  annotate('segment', x = 3, xend = 2.2, y = 20, yend = 26,
           size = 1.5, arrow = arrow(), alpha = 0.8, color = '#EEDF7C') +
  theme(axis.text = element_text(size = 12, colour = 'black'),
        axis.title = element_text(size = 14, colour = 'black'),
        plot.title = element_text(size = 15, hjust = 0.5),
        legend.position = 'none')

# Silhouette analysis ------------------------------------------------------------------------
# The silhouette analysis assesses the quality of the clustering by measuring how similar each point is to its own cluster compared to other clusters.
aSil <- allres$Silhouette
par(bty = "o", mgp = c(2.5, 0.33, 0), mar = c(5.1, 2.1, 3.1, 2.1) + 0.1, las = 1, tcl = -0.25)
plot(aSil, col = cols, main = 'Silhouette plot', border = NA)
abline(v = 0.4, lty = 2, col = "grey50", lwd = 2)
