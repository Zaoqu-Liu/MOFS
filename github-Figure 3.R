rm(list = ls())
library(tidyverse)  # Used for data manipulation and visualization
library(ggsci)      # Provides scientific color schemes
library(ComplexHeatmap) # To create complex heatmaps
library(circlize)   # Helper functions for circular visualizations
library(dendextend) # Manipulate and visualize dendrograms
library(seriation)  # Functions for seriation (finding an ordering)
library(survival)   # Functions for survival analysis
library(survminer)  # Graphical representation of survival data

# Define color scheme used throughout the analysis
color_scheme <- c('#119da4', '#FF6666', '#ffc857')


# Custom survival plot function
# This function generates customized Kaplan-Meier survival curves for the input data.
# Arguments:
# - data: The dataset containing survival information.
# - time/event: Columns for survival time and event status.
# - var: The grouping variable for survival analysis.
# - Additional parameters allow for detailed customization of plot aesthetics.
lzq_survplot <- function(data,
                         time = 'OS.time',
                         event = 'OS',
                         var = 'Cluster',
                         color = c('#119da4','#FF6666','#ffc857'),
                         title = NULL,
                         xlab.lab = 'Time in months',
                         ylab.lab = 'Overall survival',
                         rect.size = 1.2,
                         mytheme = 'bw',
                         surv.median.line = 'none',
                         legend.title = '',
                         legend.labs = c('MOFS1', 'MOFS2', 'MOFS3'),
                         legend.positions = c(0.8,0.9),
                         plot.table = TRUE) {
  require(patchwork)
  require(survival)
  require(tidyverse)
  require(survminer)
  
  # Prepare survival data for plotting
  tmp <- data[, c(time, event, var)]
  colnames(tmp) <- c('time', 'event', 'var')
  fit <- survfit(Surv(time, event) ~ var, tmp)  # Fit survival model
  
  # Custom plot formatting function for labels
  spf <- function(label = label) {
    customize_labels <- function (p, font.title = NULL, font.subtitle = NULL, font.caption = NULL,
                                  font.x = NULL, font.y = NULL, font.xtickslab = NULL, font.ytickslab = NULL) {
      original.p <- p
      if (is.ggplot(original.p)) list.plots <- list(original.p)
      else if (is.list(original.p)) list.plots <- original.p
      else
        stop("Can't handle an object of class ", class(original.p))
      
      .set_font <- function(font) {
        font <- ggpubr:::.parse_font(font)
        ggtext::element_markdown (
          size = font$size,
          face = font$face,
          colour = font$color
        )
      }
      
      for (i in 1:length(list.plots)) {
        p <- list.plots[[i]]
        if (is.ggplot(p)) {
          if (!is.null(font.title))
            p <- p + theme(plot.title = .set_font(font.title))
          if (!is.null(font.subtitle))
            p <- p + theme(plot.subtitle = .set_font(font.subtitle))
          if (!is.null(font.caption))
            p <- p + theme(plot.caption = .set_font(font.caption))
          if (!is.null(font.x))
            p <- p + theme(axis.title.x = .set_font(font.x))
          if (!is.null(font.y))
            p <- p + theme(axis.title.y = .set_font(font.y))
          if (!is.null(font.xtickslab))
            p <- p + theme(axis.text.x = .set_font(font.xtickslab))
          if (!is.null(font.ytickslab))
            p <- p + theme(axis.text.y = .set_font(font.ytickslab))
          list.plots[[i]] <- p
        }
      }
      if (is.ggplot(original.p)) list.plots[[1]]
      else list.plots
    }
    
    # Set legend position depending on whether risk table is included
    u <- c(ifelse(plot.table, 'none', legend.positions[1]), ifelse(plot.table, 'none', legend.positions[2]))
    if (is.character(u)) { u <- u[1] }
    if (mytheme == 'bw') { mytheme2 <- theme_bw(base_rect_size = 2) }
    if (mytheme == 'classic') { mytheme2 <- theme_classic() }
    
    # Create survival plot using survminer::ggsurvplot
    pp <- ggsurvplot(fit,
                     tmp,
                     pval = TRUE,
                     pval.method = TRUE,
                     ylab = NULL,
                     xlab = xlab.lab,
                     size = 1.3,
                     conf.int = FALSE,
                     surv.median.line = surv.median.line,
                     legend.title = legend.title,
                     legend.labs = legend.labs,
                     legend = u,
                     risk.table = plot.table,
                     risk.table.pos = 'out',
                     tables.col = "strata",
                     risk.table.title = "Number at risk",
                     risk.table.height = .3,
                     risk.table.y.text.col = TRUE,
                     risk.table.y.text = TRUE,
                     risk.table.y.title = FALSE,
                     palette = color,
                     font.main = 15,
                     ggtheme = mytheme2)
    
    # Customize labels for the plot
    pp$plot <- customize_labels(pp$plot,
                                font.x = c(15, "bold", "black"),
                                font.y = c(15, "bold", "black"),
                                font.xtickslab = c(12, "plain", "black"),
                                font.ytickslab = c(12, "plain", 'black'))
    
    # Adjust theme for axis lines
    if (plot.table) { xlab.lab <- NULL }
    if (mytheme == 'bw') { axis.line.x <- element_line(); axis.line.y <- element_line() }
    if (mytheme == 'classic') { axis.line.x <- element_blank(); axis.line.y <- element_line(linewidth = rect.size) }
    
    # Customize the main plot theme
    pp$plot <- pp$plot + labs(y = ylab.lab, x = xlab.lab, title = title) +
      theme(plot.title = element_text(face = "bold", colour = "black", size = 16, hjust = 0.5),
            panel.background = element_rect(fill = "#f3f6f6", color = NA),
            panel.grid.minor = element_blank(),
            panel.grid.major = if (surv.median.line == 'none') element_line(color = "#cacfd2", linetype = "dashed") else element_blank(),
            axis.line.x = axis.line.x,
            axis.line.y = axis.line.y,
            legend.position = u,
            legend.text = element_text(face = "bold", colour = "black", size = 13),
            legend.background = element_blank(),
            legend.key = element_blank())
    
    # Customize the risk table, if included
    if (plot.table) {
      pp$table <- customize_labels(pp$table,
                                   font.title = c(15, "bold", "black"),
                                   font.x = c(15, "bold", "black"),
                                   font.y = c(15, "bold", "black"),
                                   font.xtickslab = c(12, "plain", "black"),
                                   font.ytickslab = c(12, "bold")) +
        theme(panel.background = element_rect(fill = "#f3f6f6", color = NA),
              panel.grid.minor = element_blank(),
              panel.grid.major.x = element_blank(),
              panel.grid.major.y = element_line(color = "#cacfd2", linetype = "dashed"),
              axis.line.x = axis.line.x,
              axis.line.y = axis.line.y,
              legend.position = 'none')
      pp.out = pp$plot + pp$table + patchwork::plot_layout(nrow = 2, heights = c(3, 1))
    } else {
      pp.out = pp$plot
    }
    return(pp.out)
  }
  
  print(spf(label))  # Print the plot
  return(spf(label))  # Return the plot
}

# Load classifier results for downstream analysis
load('EnsembleClassifier_results.rda')

# Analysis of Clinical Data Integration and Survival Curves -----------------------
# This section integrates clustering results with clinical information to assess the survival outcome for different clusters.
# Merging the cluster assignments with clinical data helps assess how distinct molecular subtypes relate to survival outcomes, specifically overall survival (OS) and progression-free survival (PFS).

load('ClinicalData.rda')
load('pretrial_results.rda')
merged_data <- merge(allres$coca_res$Cluster, GBMdata$Clin, by = 1)
merged_data <- merged_data[merged_data$Sample %in% allres$Core_set,] # Filter to include only the core set of samples
merged_data$Cluster <- paste0('MOFS', merged_data$Cluster) # Label clusters for easier identification
merged_data$OS <- ifelse(merged_data$OS == 'Yes', 1, 0) # Convert survival status to binary format (1 = yes, 0 = no)
merged_data$Recurrence <- ifelse(merged_data$Recurrence == 'Yes', 1, 0) # Convert recurrence status similarly

# Generate Kaplan-Meier plots for both OS and PFS across the identified clusters
lzq_survplot(data = merged_data, time = 'OS.time', event = 'OS', var = 'Cluster',
             ylab.lab = 'Overall survival',
             title = 'FAHZZU (n = 116)', legend.title = '')

lzq_survplot(data = merged_data, time = 'PFS.time', event = 'Recurrence', var = 'Cluster',
             ylab.lab = 'Progression-free survival',
             title = 'FAHZZU (n = 65)', legend.title = '')

# Kaplan-Meier Survival Analysis for External Datasets ----------------------------
# This analysis focuses on external validation datasets, converting the reported survival times from years to months for comparability.
# The datasets are cross-validated to assess the robustness of subtype classification and its association with patient outcomes.

dataset_id <- 'GSE72951'
dataset_results <- EC_res[[dataset_id]]$res
range(dataset_results$OS.time, na.rm = TRUE)
dataset_results$OS.time <- dataset_results$OS.time * 12  # Convert from years to months
dataset_results$PFS.time <- dataset_results$PFS.time * 12 # Convert similarly

# Generate Kaplan-Meier plots for OS and PFS
lzq_survplot(data = dataset_results, time = 'OS.time', event = 'OS', var = 'Cluster',
             ylab.lab = 'Overall survival',
             title = paste0(dataset_id, ' (n = 110)'), legend.title = '')
ggsave(paste0('Results/Figure4D_', dataset_id, '_Subtype_OS.pdf'), height = 4.6, width = 4.6)

# Additional Analysis for Different Platforms --------------------------------------
# Subsequent analyses are performed on data from various platforms, such as GPL570, GPL6480, GPL96, etc.,
# to determine the association between identified molecular subtypes and overall survival.

names(EC_res)
dataset_id <- 'GPL570'
dataset_results <- EC_res[[dataset_id]]$res
range(dataset_results$OS.time, na.rm = TRUE)
dataset_results$OS.time <- dataset_results$OS.time * 12

lzq_survplot(data = dataset_results, time = 'OS.time', event = 'OS', var = 'Cluster',
             ylab.lab = 'Overall survival',
             title = paste0('GBM_', dataset_id, ' (n = 197)'), legend.title = '')
ggsave(paste0('Results/Figure4E_', dataset_id, '_Subtype_OS.pdf'), height = 4.6, width = 4.6)

# Cox Proportional Hazard Analysis for Module Scoring ------------------------------
# Here, the module scores for three molecular subtypes are computed to assess the effect of the activity level of each cluster on overall survival.
# Each subtype activity is assessed for its impact on survival using Cox proportional hazard models. These visualizations help understand the risk association of each subtype.

load('ClassifierTrainingData.rda')
load('ClinicalData.rda')

# Define a function to compute hazard ratios using Cox proportional hazards model
lzq_getHRdf <- function(data, var = 'C1', survar = 'OS') {
  require(survival)
  require(tidyverse)
  
  # Construct survival formula based on provided survival variable
  survival_formula <- switch(survar,
                             'OS' = as.formula(Surv(OS.time, OS) ~ get(var)),
                             'RFS' = as.formula(Surv(RFS.time, RFS) ~ get(var)),
                             'DFS' = as.formula(Surv(DFS.time, DFS) ~ get(var)),
                             'PFS' = as.formula(Surv(PFS.time, PFS) ~ get(var)),
                             'DSS' = as.formula(Surv(DSS.time, DSS) ~ get(var)))
  
  # Fit Cox proportional hazards model
  cox_model <- coxph(survival_formula, data)
  cox_summary <- summary(cox_model)
  p_value <- cox_summary$coefficient[, 5]
  coefficient <- cox_summary$coefficient[, 1]
  
  # Predict risk and compute confidence intervals
  predicted_risk <- predict(cox_model, type = "risk", se.fit = TRUE)
  hr <- predicted_risk$fit
  high <- hr + 2 * predicted_risk$se.fit
  low <- hr - 2 * predicted_risk$se.fit
  
  # Create data frame for plotting
  result_data <- data.frame(var = data[, var], HR = hr, HR.95H = high, HR.95L = low) %>%
    arrange(var) %>% mutate(x = 1:length(var))
  
  cat(paste0('Cox P value = ', round(p_value, 3)))
  cat(paste0('\nCox coefficient = ', round(coefficient, 3)))
  
  return(list(res = result_data, cox.coef = sprintf('%.3f', coefficient), cox.p = ifelse(p_value < 0.001, '0.001', sprintf('%.3f', p_value))))
}

expression_data <- rna
clinical_data <- ClinicalData
colnames(clinical_data)[1] <- 'SampleID'
clinical_data$OS <- ifelse(clinical_data$OS == 'Yes', 1, 0)  # Convert OS to binary format
clinical_data <- clinical_data[clinical_data$SampleID %in% Clu$Sample,]  # Match clinical data to sample identifiers
expression_data <- expression_data[, match(clinical_data$SampleID, colnames(expression_data))]  # Match gene expression data similarly

fit_results <- lzq_EnsembleClassfiers(expression_data, clinical_data, Need_P_cutoff = TRUE)
save(fit_results, file = 'coxfit_results.rda')

# Generate Hazard Ratio (HR) Curves for MOFS Subtypes ------------------------------
# The resulting HR plots demonstrate the relative risk associated with the activity of each subtype.
# These visualizations help quantify the impact of subtype activity on patient outcomes.

tmp <- lzq_getHRdf(data = fit_results$res, var = 'C1', survar = 'OS')
plot_hr_curve1 <- ggplot(tmp$res, aes(var, HR)) +
  geom_ribbon(aes(ymin = HR.95L, ymax = HR.95H), fill = color_scheme[1], alpha = 0.7) +
  geom_line(linewidth = 1.3, color = '#023e8a') +
  scale_x_continuous(breaks = c(0.1, 0.5, 0.9), labels = c('Low', 'MOFS1 activity', 'High'), expand = c(0, 0)) +
  scale_y_continuous(expand = c(0, 0)) +
  labs(x = NULL, y = 'Relative HR', title = NULL) +
  theme_bw(base_rect_size = 2) +
  theme(panel.grid = element_blank(),
        panel.grid.major = element_line(color = "#cacfd2", linetype = "dashed"),
        panel.background = element_rect(fill = '#f3f6f6'),
        axis.ticks.x = element_blank(),
        axis.text.x = element_text(size = 13, colour = c('#0a9396', '#582f0e', '#ae2012'), face = 'bold'),
        axis.text.y = element_text(size = 12, colour = 'black'),
        axis.title.x = element_blank(),
        axis.title.y = element_text(size = 15, colour = 'black', face = 'bold'),
        plot.title = element_text(hjust = 0.5, size = 15, colour = 'black', face = 'bold'),
        legend.position = 'none',
        legend.background = element_blank(),
        legend.title = element_blank(),
        legend.key = element_blank()) +
  annotate('text', x = 0.4, y = 1.25, fontface = 'bold.italic', hjust = 0, size = 4.8,
           label = paste0('Coefficient = ', tmp$cox.coef, '\nP = ', tmp$cox.p))
plot_hr_curve1

# Comparative Subtype Classification with External Schemes --------------------------
# To validate the consistency of subtype classification, the MOFS subtypes were compared against other classification systems, such as Wang and Phillips schemes.
# Sankey diagrams were generated to illustrate the correspondence between the subtypes across different schemes.
library(ggsankey)
sankey_data <- Clu %>% ggsankey::make_long(Phillips, FAHZZU, Wang)
sankey_data$node <- factor(sankey_data$node, c(unique(Clu$FAHZZU)[c(1, 3, 2)], unique(Clu$Phillips)[c(1, 3, 2)], unique(Clu$Wang)[c(1, 2, 3)]))

# Plot Sankey diagram
ggplot(sankey_data, aes(x = x, next_x = next_x, node = node, next_node = next_node, fill = factor(node), label = node)) +
  ggsankey::geom_sankey(flow.alpha = .6, node.color = "gray30") +
  ggsankey::geom_sankey_label(size = 3, color = "black", fill = "#9dbebb", fontface = 'bold') +
  scale_fill_manual(values = c(color_scheme, rev(color_scheme), color_scheme)) +
  ggsankey::theme_alluvial() +
  labs(x = NULL) +
  ggsankey::theme_sankey(base_size = 18) +
  theme(legend.position = "none", plot.title = element_text(hjust = .5))
