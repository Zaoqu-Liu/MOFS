rm(list = ls())
library(tidyverse)  # Data manipulation and visualization
library(furrr)  # Parallel processing
library(caret)  # Data partitioning and model training
library(glmnet)  # Lasso and Ridge regression
load('radiology.rda')
# Load preliminary results
load('preliminary-results.rda')
multi_omics_data <- multi_omics_data[, preliminary_results$Core_set]
MOFS_labels <- preliminary_results$coca_results$Cluster
MOFS_labels <- MOFS_labels[colnames(multi_omics_data), ]
MOFS_labels <- factor(MOFS_labels, paste0('MOFS', 1:3))

# Merge MOFS labels with multi-omics data
merged_data <- merge(MOFS_labels, t(multi_omics_data), by.x = "Sample", by.y = 0) %>%
  tibble::column_to_rownames('Sample')

# Perform bootstrap analysis to identify important features for each MOFS subtype
bootstrap_analysis <- function(data, target_label, label_name) {
  # Subset data to include only the specified label
  data$MOFS_Label <- ifelse(data$MOFS_Label == target_label, 1, 0)
  
  # Initial feature filtering using univariate logistic regression
  significant_features <- data[, -1] %>% apply(., 2, function(feature) {
    model_fit <- summary(glm(data$MOFS_Label ~ feature, family = binomial))
    return(model_fit$coefficients[2, 4])
  })
  
  # Retain features with p-value < 0.01
  significant_features <- significant_features[significant_features < 0.01]
  filtered_data <- data[, c('MOFS_Label', names(significant_features))]
  
  # Perform bootstrap resampling to validate significant features
  plan(multisession)  # Enable parallel processing
  bootstrap_results <- future_map(colnames(filtered_data[, -1]), function(feature_name) {
    bootstrap_p_values <- map_vec(1:1000, function(iteration) {
      # Create training and validation datasets using 70% of the data
      partition_indices <- createDataPartition(filtered_data$MOFS_Label, p = 0.7, list = FALSE) %>% as.numeric()
      training_data <- filtered_data[partition_indices, ]
      
      # Fit logistic regression model to resampled data
      model_fit <- summary(glm(training_data$MOFS_Label ~ training_data[, feature_name], family = binomial))
      return(model_fit$coefficients[2, 4])
    })
    return(bootstrap_p_values)
  }, .progress = TRUE)
  
  # Summarize bootstrap results to identify stable features
  bootstrap_summary <- Reduce(cbind, bootstrap_results) %>% as.data.frame()
  selected_features <- data.frame(
    Feature = colnames(filtered_data[, -1]),
    SelectedCount = apply(bootstrap_summary, 2, function(p_values) sum(p_values < 0.05))
  )
  stable_features <- selected_features$Feature[selected_features$SelectedCount > 950]
  
  return(stable_features)
}

# Perform bootstrap analysis for each MOFS subtype
features_MOFS1 <- bootstrap_analysis(merged_data, 'MOFS1', 'MOFS1')
features_MOFS2 <- bootstrap_analysis(merged_data, 'MOFS2', 'MOFS2')
features_MOFS3 <- bootstrap_analysis(merged_data, 'MOFS3', 'MOFS3')

# Identify unique features for each MOFS subtype
unique_features_1 <- setdiff(features_MOFS1, c(features_MOFS2, features_MOFS3))
unique_features_2 <- setdiff(features_MOFS2, c(features_MOFS1, features_MOFS3))
unique_features_3 <- setdiff(features_MOFS3, c(features_MOFS1, features_MOFS2))

# Combine all unique features
combined_features <- Reduce(union, list(unique_features_1, unique_features_2, unique_features_3))

# Perform Lasso regression to further refine selected features
lasso_analysis <- function(data, features, target_label) {
  # Subset data for the specified MOFS label and features
  subset_data <- data[, c('MOFS_Label', features)]
  subset_data$MOFS_Label <- ifelse(subset_data$MOFS_Label == target_label, 1, 0)
  
  # Fit Lasso regression model using cross-validation
  set.seed(1234)
  cv_fit <- cv.glmnet(as.matrix(subset_data[, -1]), subset_data$MOFS_Label, family = "binomial")
  lasso_coefficients <- coef(cv_fit, s = "lambda.min")
  selected_features <- which(lasso_coefficients != 0)[-1] - 1
  refined_features <- colnames(subset_data[, -1])[selected_features]
  
  return(refined_features)
}

# Refine features for each MOFS subtype using Lasso regression
refined_features_MOFS1 <- lasso_analysis(merged_data, unique_features_1, 'MOFS1')
refined_features_MOFS2 <- lasso_analysis(merged_data, unique_features_2, 'MOFS2')
refined_features_MOFS3 <- lasso_analysis(merged_data, unique_features_3, 'MOFS3')

# Combine refined features from all MOFS subtypes
final_features <- Reduce(union, list(refined_features_MOFS1, refined_features_MOFS2, refined_features_MOFS3))

# Prepare final dataset for machine learning model training
final_data <- merged_data[, c('MOFS_Label', final_features)]
feature_mapping <- data.frame(
  FeatureLabel = paste0('Feature_', 1:length(final_features)),
  OriginalFeature = colnames(final_data)[-1]
)
colnames(final_data)[-1] <- feature_mapping$FeatureLabel

# Save final dataset and feature mapping for future use
save(feature_mapping, final_data, file = 'MOFS_features_for_AI.rda')


## Import Libraries
rm(list = ls())
library(NeuralNetTools)  # For neural network analysis and visualization
library(ggpol)  # For plotting, especially confusion matrices
library(tidyverse)  # For data manipulation and visualization
library(neuralnet)  # For neural network modeling
library(furrr)  # For parallel processing of functions

## Load Data
# Load MRI feature dataset for classification
load('MRI-feature_for_AI.rda')
data$Cluster <- factor(data$Cluster)  # Convert cluster labels to a factor type
set.seed(123)  # Set seed for reproducibility

# Split Data
# Randomly split data into training (80%) and testing (20%)
train_indices <- sample(x = 1:nrow(data), size = 0.8 * nrow(data))
training_data <- data[train_indices,]
testing_data <- data[-train_indices,]

## Hyperparameter Tuning for Neural Network Model
# Define ranges for tuning parameters
learning_rates <- c(0.1, 0.05, 0.01)
num_neurons_layer1 <- 3:20
num_neurons_layer2 <- 0:20
algorithms <- c('backprop', 'rprop+', 'rprop-')
parameter_combinations <- expand.grid(learning_rates, num_neurons_layer1, num_neurons_layer2, algorithms)

# Plan parallel processing for model training
plan(multisession)

# Train models for different combinations of hyperparameters
results <- future_map(1:nrow(parameter_combinations), function(i) {
  set.seed(123)  # Set seed for consistency in model training
  combination <- parameter_combinations[i,]
  num_neurons <- if (combination[3] == 0) combination[2] else combination[2:3]
  num_neurons <- as.numeric(num_neurons)
  
  # Train neural network model
  model <- tryCatch(
    neuralnet(
      Cluster ~ ., 
      data = training_data,
      hidden = num_neurons,  # Define hidden layers and number of neurons per layer
      err.fct = 'sse',  # Error function; 'sse' for sum of squared errors
      act.fct = 'logistic',  # Activation function; logistic or tanh
      linear.output = FALSE,  # False indicates that output is non-linear
      learningrate = combination[1],  # Learning rate for backpropagation algorithm
      algorithm = combination[4],  # Algorithm to optimize weights
      likelihood = FALSE  # Calculate AIC/BIC if set to TRUE
    ),
    error = function(x) NA
  )
  
  # Evaluate model accuracy on training, testing, and full datasets
  train_accuracy <- tryCatch(sum(paste0('C', apply(predict(model, training_data[,-1]), 1, which.max)) == training_data$Cluster) / nrow(training_data), error = function(x) NA)
  validation_accuracy <- tryCatch(sum(paste0('C', apply(predict(model, testing_data), 1, which.max)) == testing_data$Cluster) / nrow(testing_data), error = function(x) NA)
  overall_accuracy <- tryCatch(sum(paste0('C', apply(predict(model, data), 1, which.max)) == data$Cluster) / nrow(data), error = function(x) NA)
  
  return(c(train_accuracy, validation_accuracy, overall_accuracy))
}, .progress = TRUE)

# Combine results into a single data frame for further analysis
results_summary <- Reduce(rbind, results)
colnames(results_summary) <- c('train_accuracy', 'validation_accuracy', 'overall_accuracy')
results_summary <- cbind(parameter_combinations, results_summary)

## Train Final Neural Network Model for Classification
# Reload data
load('MRI-feature_for_AI.rda')
data$Cluster <- gsub('C', 'MOFS', data$Cluster)  # Rename clusters to MOFS
data$Cluster <- factor(data$Cluster)
set.seed(123)

# Split data into training and testing sets
train_indices <- sample(x = 1:nrow(data), size = 0.8 * nrow(data))
training_data <- data[train_indices,]
testing_data <- data[-train_indices,]

# Train final model with chosen hyperparameters
set.seed(123)
final_model <- tryCatch(
  neuralnet(
    Cluster ~ ., 
    data = training_data,
    hidden = c(13, 6),  # Chosen number of neurons in two hidden layers
    err.fct = 'sse',
    act.fct = 'logistic',
    linear.output = FALSE,
    learningrate = 0.1,
    algorithm = 'rprop+'  # Resilient backpropagation algorithm
  ),
  error = function(x) NA
)

## Model Evaluation and Visualization
# Variable importance analysis using connection weights
olden(final_model) + ggtitle("Variable importance using connection weights") + coord_flip()

# Plotting the generalized weights for each output
neuralnet::gwplot(final_model, selected.response = 'MOFS1', highlight = TRUE)
neuralnet::gwplot(final_model, selected.response = 'MOFS2', highlight = TRUE)
neuralnet::gwplot(final_model, selected.response = 'MOFS3', highlight = TRUE)

# Plot network structure for visual inspection
plotnet(
  mod_in = final_model, 
  nid = FALSE,  # Do not use neural interpretation diagram
  bias = FALSE,
  circle_cex = 5,  # Size of the node
  node_labs = FALSE,
  var_labs = FALSE,
  circle_col = '#0096c7',  # Color for the nodes
  bord_col = NULL,
  max_sp = TRUE  # Stretch each layer
)

## Confusion Matrix for Model Predictions
# Predict on entire dataset
predicted_values <- predict(final_model, data)
predicted_labels <- apply(predicted_values, 1, which.max)  # Get predicted class labels

# Confusion matrix visualization
confusion_matrix <- table(data$Cluster, predicted_labels) %>%
  as.data.frame() %>%
  mutate(predicted_labels = paste0('MOFS', predicted_labels)) %>%
  pivot_wider(names_from = 'predicted_labels', values_from = 'Freq') %>%
  tibble::column_to_rownames('Var1')

## Heatmap of Confusion Matrix
library(circlize)
library(ComplexHeatmap)

cell_function <- function(j, i, x, y, width, height, fill) {
  if (confusion_matrix[i, j] > 4) 
    grid.text(confusion_matrix[i, j], x, y, gp = gpar(fontsize = 14, fontface = 'bold', col = 'black'))
  else 
    grid.text(confusion_matrix[i, j], x, y, gp = gpar(fontsize = 14, fontface = 'plain', col = 'black'))
}

colors <- c('#119da4', '#FF6666', '#ffc857')

TopAnnotation <- HeatmapAnnotation(
  Subtype = paste0('MOFS', 1:3),
  border = FALSE,
  show_legend = FALSE,
  show_annotation_name = FALSE,
  simple_anno_size = unit(3, 'mm'),
  col = list(Subtype = c('MOFS1' = colors[1], 'MOFS2' = colors[2], 'MOFS3' = colors[3]))
)

LeftAnnotation <- HeatmapAnnotation(
  Subtype = paste0('MOFS', 1:3),
  border = FALSE,
  show_legend = FALSE,
  show_annotation_name = FALSE,
  simple_anno_size = unit(3, 'mm'),
  which = 'row',
  col = list(Subtype = c('MOFS1' = colors[1], 'MOFS2' = colors[2], 'MOFS3' = colors[3]))
)

# Define color gradients for heatmap
blue <- "#204F8D"
light_blue <- "#498EB9"
dark_white <- "#B6D1E8"
white <- "#E6EAF7"

Heatmap(
  confusion_matrix,
  name = ' ',
  border = TRUE,
  rect_gp = gpar(col = 'black'),
  top_annotation = TopAnnotation,
  left_annotation = LeftAnnotation,
  col = colorRamp2(c(0, 3, 6, 9), c(white, dark_white, light_blue, blue)),
  cell_fun = cell_function,
  show_row_names = FALSE,
  show_column_names = FALSE,
  cluster_rows = FALSE,
  cluster_columns = FALSE,
  show_heatmap_legend = FALSE,
  column_title = 'Predicted',
  row_title = 'Actual',
  gap = unit(1.2, 'mm'),
  column_gap = unit(1.5, 'mm')
)

## ROC Curves for Model Performance
library(pROC)

# Generate ROC curves for each class
roc_data <- list()
for (class_label in unique(data$Cluster)) {
  roc_curve <- roc(ifelse(data$Cluster == class_label, 1, 0), predicted_values[, as.numeric(gsub('MOFS', '', class_label))])
  roc_data[[class_label]] <- roc_curve
}

# Visualize ROC curves
colors <- c('#119da4', '#FF6666', '#ffc857')

for (i in seq_along(roc_data)) {
  ggplot(data.frame(x = 1 - roc_data[[i]]$specificities, y = roc_data[[i]]$sensitivities), aes(x, y)) +
    geom_line(size = 1.3, color = colors[i]) +
    labs(x = '1 - Specificity', y = 'Sensitivity', color = NULL, title = paste('ROC Curve for', names(roc_data)[i])) +
    theme_bw(base_rect_size = 2.5) +
    geom_abline(slope = 1, color = 'grey70') +
    scale_x_continuous(expand = c(0.01, 0.01)) +
    scale_y_continuous(expand = c(0.01, 0.01)) +
    theme(
      axis.text.y = element_text(size = 12, colour = 'black'),
      axis.text.x = element_text(size = 12, colour = 'black'),
      axis.title.x = element_text(size = 15, colour = 'black', face = 'bold'),
      axis.title.y = element_text(size = 15, colour = 'black', face = 'bold'),
      panel.grid = element_blank(),
      panel.grid.major = element_line(color = "#cacfd2", linetype = "dashed"),
      panel.background = element_rect(fill = '#f3f6f6'),
      plot.title = element_text(hjust = 0.5, size = 15, colour = 'black', face = 'bold'),
      legend.text = element_text(size = 12),
      legend.key = element_blank(),
      legend.background = element_blank(),
      legend.position = c(0.995, 0.03),
      legend.justification = c(1, 0)
    )
}

## Summary and Objective Explanation
# The goal of this code is to build, tune, and evaluate a neural network classifier on MRI features
# of a medical dataset. The model is used to predict cluster membership based on features provided.
# We utilized a neural network with various hyperparameter combinations to determine the best performance.
# Different learning rates, hidden layer configurations, and training algorithms were explored.
# We evaluated the model performance by calculating accuracy on the training, testing, and overall datasets.
# Finally, we used several visualizations such as connection weight importance, confusion matrix heatmaps,
# and ROC curves to assess the classifier's effectiveness and to understand its decision-making process.
