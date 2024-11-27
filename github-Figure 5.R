rm(list = ls())

# Load required libraries for data processing and visualization.
library(tidyverse)  # Includes ggplot2 for data visualization and dplyr for data manipulation
library(ggsci)  # Provides aesthetic color palettes
library(ComplexHeatmap)  # For heatmap generation
library(circlize)  # Provides circular visualization functions, often used with ComplexHeatmap
library(dendextend)  # To extend functionality of dendrograms, useful for hierarchical clustering
library(seriation)  # For finding the optimal ordering of data matrices
library(survival)  # Used for survival analysis
library(survminer)  # Provides utilities for survival analysis and generating Kaplan-Meier plots
library(patchwork)  # For arranging multiple ggplot figures

# Load Data
# Load preprocessed data for Tumor Microenvironment (TME) and Immune Checkpoint (ICP) related analysis.
load('TME_data.rda')  # Contains RNA-based TME data
load('ICP_data.rda')  # Contains Immune Checkpoint related data

# Define Custom Color Palette for Clusters
# This color palette will be used throughout the plots to ensure consistency between figures.
cluster_colors <- c('#119da4', '#FF6666', '#ffc857')

# Generate Labels for Clusters
# Creating labels for different MOFS clusters to be displayed in the figures.
cluster_labels_plot <- ggplot() + 
  geom_label(data = data.frame(Cluster = paste0('MOFS', 1:3)), 
             aes(label = Cluster, x = Cluster, fill = Cluster), 
             y = 0.5, color = "white", fontface = 'bold', 
             size = 12/.pt, hjust = 0.5, vjust = 0.5) + 
  scale_fill_manual(values = cluster_colors) + 
  theme_void() + 
  theme(legend.position = 'none', plot.margin = unit(c(0.1, 0.5, 0.1, 0.5), "cm"))

# Display the Label Plot
print(cluster_labels_plot)

# Figure: Tumor Purity Analysis ---------------------------------------------------
# This section aims to visualize the tumor purity data across different MOFS clusters.
# Load RNA TME data, modify cluster labels for consistency, and visualize tumor purity metrics.
tme_data <- RNA_tme$res
# Rename the clusters to represent MOFS categories.
tme_data$Cluster <- gsub('C', 'MOFS', tme_data$Cluster)

# Plot Tumor Purity Analysis for Different MOFS Clusters
tumor_purity_plot <- ggplot(tme_data, aes(x = Cluster, y = Purity)) +
  geom_jitter(shape = 21, size = 2.7, width = 0.2, aes(fill = Cluster, color = Cluster)) +
  geom_boxplot(outlier.colour = NA, aes(fill = Cluster), notch = FALSE, color = 'black', size = 0.6, alpha = 0.65) +
  geom_violin(alpha = 0.5, aes(fill = Cluster), color = NA, trim = TRUE) +
  # Add pairwise comparisons between clusters for statistical significance.
  stat_compare_means(comparisons = list(c('MOFS1', 'MOFS2'), c('MOFS1', 'MOFS3'), c('MOFS2', 'MOFS3')), label = 'p.signif', tip.length = 0, vjust = 0.7, size = 5) +
  scale_fill_manual(values = cluster_colors) +
  scale_color_manual(values = cluster_colors) +
  theme_classic(base_rect_size = 2.5) +
  labs(x = NULL, y = 'Tumor purity') +
  theme(
    axis.text.y = element_text(size = 12, colour = 'black'),
    axis.text.x = element_blank(),
    axis.title.x = element_blank(),
    axis.title.y = element_text(size = 15, colour = 'black', face = 'bold'),
    axis.line.x = element_blank(),
    axis.line.y = element_line(linewidth = 1.2),
    panel.grid = element_blank(),
    panel.grid.major = element_line(color = "#cacfd2", linetype = "dashed"),
    panel.background = element_rect(fill = '#f3f6f6'),
    plot.title = element_text(hjust = 0.5, size = 15, colour = 'black', face = 'bold'),
    legend.position = 'none',
    axis.ticks.x = element_blank(),
    plot.margin = unit(c(0.1, 0.5, 0, 0.5), "cm")
  )

# Display the Tumor Purity Plot
print(tumor_purity_plot)

# Figure: Immune Score Analysis ---------------------------------------------------
# This figure aims to visualize the immune score across different MOFS clusters.
# Immune scores represent the extent of immune cell infiltration in the tumor microenvironment.
immune_score_plot <- ggplot(tme_data, aes(x = Cluster, y = ImmuneScore)) +
  geom_jitter(shape = 21, size = 2.7, width = 0.2, aes(fill = Cluster, color = Cluster)) +
  geom_boxplot(outlier.colour = NA, aes(fill = Cluster), notch = FALSE, color = 'black', size = 0.6, alpha = 0.65) +
  geom_violin(alpha = 0.5, aes(fill = Cluster), color = NA, trim = TRUE) +
  stat_compare_means(comparisons = list(c('MOFS1', 'MOFS2'), c('MOFS1', 'MOFS3'), c('MOFS2', 'MOFS3')), label = 'p.signif', tip.length = 0, vjust = 0.7, size = 5) +
  scale_fill_manual(values = cluster_colors) +
  scale_color_manual(values = cluster_colors) +
  theme_classic(base_rect_size = 2.5) +
  labs(x = NULL, y = 'Immune Score') +
  theme(
    axis.text.y = element_text(size = 12, colour = 'black'),
    axis.text.x = element_blank(),
    axis.title.x = element_blank(),
    axis.title.y = element_text(size = 15, colour = 'black', face = 'bold'),
    axis.line.x = element_blank(),
    axis.line.y = element_line(linewidth = 1.2),
    panel.grid = element_blank(),
    panel.grid.major = element_line(color = "#cacfd2", linetype = "dashed"),
    panel.background = element_rect(fill = '#f3f6f6'),
    plot.title = element_text(hjust = 0.5, size = 15, colour = 'black', face = 'bold'),
    legend.position = 'none',
    axis.ticks.x = element_blank(),
    plot.margin = unit(c(0.1, 0.5, 0, 0.5), "cm")
  )

# Display the Immune Score Plot
print(immune_score_plot)

# Figure: Stromal Score Analysis ---------------------------------------------------
# This plot illustrates the stromal score for the different MOFS clusters. The stromal score indicates the presence of stromal cells in the tumor microenvironment, which can influence the tumor's response to treatment and progression.
stromal_score_plot <- ggplot(tme_data, aes(x = Cluster, y = StromalScore)) +
  geom_jitter(shape = 21, size = 2.7, width = 0.2, aes(fill = Cluster, color = Cluster)) +
  geom_boxplot(outlier.colour = NA, aes(fill = Cluster), notch = FALSE, color = 'black', size = 0.6, alpha = 0.65) +
  geom_violin(alpha = 0.5, aes(fill = Cluster), color = NA, trim = TRUE) +
  stat_compare_means(comparisons = list(c('MOFS1', 'MOFS2'), c('MOFS1', 'MOFS3'), c('MOFS2', 'MOFS3')), label = 'p.signif', tip.length = 0, vjust = 0.7, size = 5) +
  scale_fill_manual(values = cluster_colors) +
  scale_color_manual(values = cluster_colors) +
  theme_classic(base_rect_size = 2.5) +
  labs(x = NULL, y = 'Stromal Score') +
  theme(
    axis.text.y = element_text(size = 12, colour = 'black'),
    axis.text.x = element_blank(),
    axis.title.x = element_blank(),
    axis.title.y = element_text(size = 15, colour = 'black', face = 'bold'),
    axis.line.x = element_blank(),
    axis.line.y = element_line(linewidth = 1.2),
    panel.grid = element_blank(),
    panel.grid.major = element_line(color = "#cacfd2", linetype = "dashed"),
    panel.background = element_rect(fill = '#f3f6f6'),
    plot.title = element_text(hjust = 0.5, size = 15, colour = 'black', face = 'bold'),
    legend.position = 'none',
    axis.ticks.x = element_blank(),
    plot.margin = unit(c(0.1, 0.5, 0, 0.5), "cm")
  )

# Display the Stromal Score Plot
print(stromal_score_plot)

# Combining Plots for Comparative Analysis ----------------------------------------
# The next step is to arrange the different score plots (tumor purity, immune score, stromal score) for comparative analysis across the MOFS clusters.
combined_plot <- tumor_purity_plot + cluster_labels_plot + immune_score_plot + cluster_labels_plot + stromal_score_plot + plot_layout(ncol = 1, heights = c(15, 1, 15, 1, 15))

# Display the Combined Plot
print(combined_plot)

# Explanation:
# - The jitter and box plots are combined to provide an overview of the distribution of data points for each MOFS cluster.
# - Violin plots add information about the data density across different values.
# - The `stat_compare_means()` function is used to highlight statistically significant differences between clusters.
# - The combined plot allows for visual comparison between tumor purity, immune score, and stromal score, illustrating how different MOFS clusters vary in terms of their microenvironmental characteristics.

# Important Considerations:
# The consistent use of colors for the clusters (MOFS1, MOFS2, MOFS3) across all figures allows readers to intuitively link the results and quickly identify trends. The goal of these analyses is to characterize the distinct microenvironmental profiles of the different MOFS clusters, which can have significant implications for understanding tumor heterogeneity and guiding treatment strategies.


# Figure 7 Analysis ---------------------------------------------------------------

# Clear the workspace to ensure no variables or objects interfere with the analysis
rm(list = ls())

# Load necessary data for the analysis, which includes tumor microenvironment (TME) data and immune checkpoint (ICP) data
load('TME.rda')
load('ICP.rda')

# Extract TME-related RNA expression data and rename some columns for better readability
rna_tme_data <- RNA_tme$res
rna_tme_data$Cluster <- gsub('C','MOFS',rna_tme_data$Cluster) # Rename clusters to MOFS
colnames(rna_tme_data) # Check column names of the dataset

# Retain only selected columns: Sample ID, Cluster information, and specific cell types
rna_tme_data <- rna_tme_data[,c(1,2,9:22)]
colnames(rna_tme_data)[c(7,13)] <- c('Dendritic_Cells','Oligodendrocytes') # Rename columns for clarity

# Initialize an empty dataframe for storing statistical test results
test_results <- data.frame()     

# Perform a one-way ANOVA test for each cell type to determine statistical significance between different clusters
for (cell_type in colnames(rna_tme_data)[-c(1,2)]) {
  fit <- oneway.test(get(cell_type) ~ Cluster, data = rna_tme_data) # Perform ANOVA
  p_value <- fit$p.value
  test_results <- rbind(test_results, data.frame(CellType = cell_type, P_Value = p_value))
}

# Adjust p-values using Benjamini-Hochberg method for controlling the false discovery rate
test_results$Adjusted_P_Value <- p.adjust(test_results$P_Value, method = 'BH')

# Load a custom function to annotate significance levels based on adjusted p-values
source('lzq_signif_p.R')
test_results$Significance <- lzq_signif_p(test_results$Adjusted_P_Value)

# Combine TME cell results with ICP differential expression results for comprehensive analysis
combined_results <- rbind(test_results, ICP$ICP_diff)
combined_results$Significance <- ifelse(combined_results$Significance == 'ns', '', combined_results$Significance)

# Prepare a gene list for visualization combining TME and ICP genes
gene_list <- rbind(data.frame(GeneID = test_results$CellType, Category = 'TME Cells', Label = test_results$CellType), ICP$ICP_genelist)
gene_list$Category <- factor(gene_list$Category, unique(gene_list$Category))
gene_list$GeneID <- gsub('\.', '-', gene_list$GeneID) # Replace "." with "-" for consistency

# Merge expression data from TME and ICP datasets for heatmap visualization
expression_data <- merge(rna_tme_data[,-2], t(ICP$ICP_expr_mRNA), by.x = 1, by.y = 0) %>%
  tibble::column_to_rownames('Sample') %>%
  t() %>%
  as.data.frame()
rownames(expression_data) <- gsub('\.', '-', rownames(expression_data))

# Arrange cluster information for ordering columns in the heatmap
cluster_info <- rna_tme_data[,1:2] %>% arrange(Cluster) %>% mutate(Cluster = gsub('C', 'MOFS', Cluster))
expression_data <- expression_data[gene_list$GeneID, cluster_info$Sample] %>% na.omit()

# Update gene list to match the filtered data
gene_list <- gene_list[match(rownames(expression_data), gene_list$GeneID),]

# Standardize expression data by scaling each row to have zero mean and unit variance
expression_data <- t(scale(t(expression_data)))

# Generate a heatmap to visualize expression patterns across different clusters
heatmap_plot <- Heatmap(
  expression_data,
  name = 'Expression Level',
  color_space = "RGB",
  border = TRUE,
  border_gp = gpar(col = 'black', lwd = 3),
  col = colorRamp2(c(-2, -1, 0, 1, 2), c('#0077b6', '#48cae4', 'white', ggplot2::alpha('#e3170a', 0.5), '#e3170a')),
  cluster_rows = FALSE,
  cluster_columns = FALSE,
  show_column_dend = FALSE,
  show_row_dend = FALSE,
  show_row_names = TRUE,
  row_names_side = 'left',
  row_names_gp = gpar(fontsize = 9),
  show_column_names = FALSE,
  row_split = gene_list$Category,
  row_title_gp = gpar(fontsize = 12, col = 'black', fontface = 'bold'),
  column_split = cluster_info$Cluster,
  column_title = NULL,
  top_annotation = HeatmapAnnotation(
    Subtype = anno_block(
      gp = gpar(fill = c('#119da4', '#FF6666', '#ffc857'), lwd = 3),
      labels = paste0('MOFS', 1:3),
      height = unit(5, 'mm'),
      labels_gp = gpar(cex = 0.8, col = "white", fontface = 'bold')
    )
  ),
  right_annotation = rowAnnotation(foo = anno_text(combined_results$Significance, gp = gpar(cex = 1, col = 'black'))),
  gap = unit(2, 'mm'),
  column_gap = unit(2, 'mm'),
  heatmap_legend_param = list(
    labels_gp = gpar(fontsize = 10, fontface = 'bold'),
    border = TRUE,
    legend_height = unit(50, 'mm'),
    title_gp = gpar(fontsize = 10, fontface = "bold")
  )
)

# This analysis aims to visualize the differential expression of tumor microenvironment (TME) cells and immune checkpoint (ICP) genes across different MOFS clusters. By using a heatmap, we can easily identify patterns of gene expression that distinguish these clusters, offering insight into the biological differences between them. The heatmap also highlights the statistically significant differences, helping prioritize targets for further investigation.

# Load necessary libraries for the radar chart
library(fmsb) 
library(RColorBrewer)

# Define colors for different groups
cols <- c('#119da4', '#FF6666', '#ffc857')

# Create a radar chart to visualize immune checkpoint-related data
radarchart(ICP$radar_data,
           axistype = 0, # Use 0 to not display data labels, or 1 to display them
           seg = 7, # Number of segments for the axes
           caxislabels = seq(0, 7, 1), # Labels for the axis, representing scores from 0 to 7
           calcex = 0.8, # Size of the axis labels
           pty = 32, # Type of point, 32 means points are not displayed
           pcol = cols, # Color of the lines for different groups
           plwd = 3, # Line width
           plty = 1, # Line type (solid lines)
           pfcol = ggplot2::alpha(cols, 0.1), # Color fill for different groups with transparency
           cglty = 2, # Type of the grid line (dashed)
           cglwd = 1, # Width of the grid line
           cglcol = 'grey30', # Color of the grid line
           axislabcol = "black", # Color of axis labels and numbers
           title = NULL, # Title of the radar chart (set to NULL if no title)
           maxmin = TRUE, # TRUE means that the first and second rows should contain max and min values
           na.itp = FALSE, # NA values should be treated as zero, rather than interpolated
           centerzero = TRUE, # Set this to TRUE if zero should be treated as the center of the chart
           vlabels = NULL, # Variable labels, set to NULL to use the column names as labels
           vlcex = 1 # Font size of the variable labels
)

# Generate bar plot and box plot for further analysis
library(ggsci)

# Bar plot of response to Anti-PD-1 immunotherapy
my2 <- dd %>%
  count(Cluster, Type) %>%
  group_by(Cluster) %>%
  summarise(Res = Type, nn = n / sum(n))
my2 <- my2[order(my2$Res, decreasing = TRUE),]
my2 <- group_by(my2, Cluster) %>% mutate(label_y = cumsum(nn) - 0.5 * nn)

my2$ll <- ifelse(my2$nn != 0, paste0(round(my2$nn, 2) * 100, '%'), '')
my2$Res <- ifelse(my2$Res == 'NR', 'Non-Responsers', 'Responsers')

p1 <- ggplot(my2, aes(Cluster, nn, fill = Res)) +
  geom_bar(stat = 'identity', color = 'black', width = 0.6) +
  geom_text(aes(y = label_y, label = ll), fontface = 'bold', color = 'white', size = 4.5) +
  scale_fill_manual(values = c('#023e8a', '#c1121f')) + 
  theme_classic() +
  scale_y_continuous(expand = c(0, 0)) +
  labs(x = NULL, y = 'Anti-PD-1 immunotherapy', fill = NULL) +
  scale_x_discrete(expand = c(0, 0.5)) +
  theme(
    axis.text.y = element_text(size = 12, colour = 'black'),
    axis.text.x = element_text(size = 15, colour = 'black', face = 'bold', angle = 60, hjust = 1),
    axis.title.x = element_blank(),
    axis.title.y = element_text(size = 15, colour = 'black', face = 'bold'),
    axis.line.x = element_blank(),
    axis.line.y = element_line(linewidth = 1.2),
    plot.title = element_text(size = 15, hjust = 0.5, colour = 'black', face = 'bold'),
    panel.grid = element_blank(),
    panel.grid.major = element_line(color = "#cacfd2", linetype = "dashed"),
    panel.background = element_rect(fill = '#f3f6f6'),
    legend.position = 'none',
    legend.text = element_text(size = 13, colour = 'black', face = 'bold'),
    legend.title = element_blank(),
    legend.background = element_rect(fill = NA),
    legend.key.width = unit(4, 'mm'),
    legend.key.height = unit(4, 'mm')
  )

# Box plot of subtype activities
p2 <- ggplot(dd2, aes(nn, vv, fill = Type, color = Type)) +
  scale_fill_manual(values = c('#023e8a', '#c1121f')) +
  scale_color_manual(values = c('#023e8a', '#c1121f')) +
  geom_boxplot(outlier.colour = 'grey30', notch = FALSE, outlier.size = 0.6, color = 'black', size = 0.6, alpha = 0.9) +
  scale_x_discrete(expand = c(0, 0.5)) +
  scale_y_continuous(expand = c(0, 0), limits = c(0, 1)) +
  stat_compare_means(label = 'p.signif', size = 7, label.y = 0.92) +
  labs(x = NULL, y = 'Subtype activities') +
  theme(
    axis.text.y = element_text(size = 12, colour = 'black'),
    axis.text.x = element_text(size = 15, colour = 'black', face = 'bold', angle = 60, hjust = 1),
    axis.title.x = element_blank(),
    axis.title.y = element_text(size = 15, colour = 'black', face = 'bold'),
    axis.line.x = element_blank(),
    axis.line.y = element_line(linewidth = 1.2),
    plot.title = element_text(size = 15, hjust = 0.5, colour = 'black', face = 'bold'),
    panel.grid = element_blank(),
    panel.grid.major = element_line(color = "#cacfd2", linetype = "dashed"),
    panel.background = element_rect(fill = '#f3f6f6'),
    legend.position = 'none',
    legend.text = element_text(size = 13, colour = 'black', face = 'bold'),
    legend.title = element_blank(),
    legend.background = element_rect(fill = NA),
    legend.key.width = unit(4, 'mm'),
    legend.key.height = unit(4, 'mm')
  )

# Combine the plots using patchwork
library(patchwork)
p1 + p2 + plot_layout(widths = c(2, 4))











