#' @title Calculate Coefficient of Variation for a Numeric Vector
#' @description This function calculates the coefficient of variation (CV) for a given numeric vector.
#' @author Zaoqu Liu; E-mail: liuzaoqu@163.com
#' @param V A numeric vector.
#' @return The coefficient of variation (CV) for the input vector.

CV <- function(V) {
  mean_value <- mean(V, na.rm = TRUE)
  sd_value <- sd(V, na.rm = TRUE)
  return(sd_value / mean_value)
}

#' @title Calculate Coefficient of Variation for a Data Frame
#' @description This function calculates the coefficient of variation (CV) for each column in a data frame.
#' @author Zaoqu Liu; E-mail: liuzaoqu@163.com
#' @param df A data frame with row samples and column features.
#' @return A sorted vector of CV values for each column in the data frame, in decreasing order.

CV.df <- function(df) {
  d <- sapply(as.data.frame(df), CV)
  return(sort(d, decreasing = TRUE))
}

#' @title Calculate Standard Deviation for a Data Frame
#' @description This function calculates the standard deviation (SD) for each column in a data frame.
#' @author Zaoqu Liu; E-mail: liuzaoqu@163.com
#' @param df A data frame with row samples and column features.
#' @return A sorted vector of SD values for each column in the data frame, in decreasing order.

SD.df <- function(df) {
  d <- sapply(as.data.frame(df), \(x)sd(x, na.rm = TRUE))
  return(sort(d, decreasing = TRUE))
}

#' @title Calculate Mean Value for a Data Frame
#' @description This function calculates the mean value for each column in a data frame.
#' @author Zaoqu Liu; E-mail: liuzaoqu@163.com
#' @param df A data frame with row samples and column features.
#' @return A vector of mean values for each column in the data frame.

Mean.df <- function(df) {
  sapply(as.data.frame(df), \(x)mean(x, na.rm = TRUE))
}

#' @title Calculate Median Value for a Data Frame
#' @description This function calculates the median value for each column in a data frame.
#' @author Zaoqu Liu; E-mail: liuzaoqu@163.com
#' @param df A data frame with row samples and column features.
#' @return A vector of median values for each column in the data frame.

Median.df <- function(df) {
  sapply(as.data.frame(df), \(x)median(x, na.rm = TRUE))
}

#' @title Calculate Median Absolute Deviation for a Data Frame
#' @description This function calculates the median absolute deviation (MAD) for each column in a data frame.
#' @author Zaoqu Liu; E-mail: liuzaoqu@163.com
#' @param df A data frame with row samples and column features.
#' @return A sorted vector of MAD values for each column in the data frame, in decreasing order.

MAD.df <- function(df) {
  d <- sapply(as.data.frame(df), \(x)mad(x, na.rm = TRUE))
  return(sort(d, decreasing = TRUE))
}

#' @title Minmax Normalization for a Numeric Vector
#' @description This function performs min-max normalization on a numeric vector.
#' @author Zaoqu Liu; E-mail: liuzaoqu@163.com
#' @param x A numeric vector.
#' @return A normalized numeric vector with values between 0 and 1.

minmax <- function(x) {
  (x - min(x)) / (max(x) - min(x))
}

#' @title Minmax Normalization for a Data Frame
#' @description This function performs min-max normalization on each column in a data frame.
#' @author Zaoqu Liu; E-mail: liuzaoqu@163.com
#' @param data A data frame with row samples and column features.
#' @return A data frame with normalized values between 0 and 1 for each column.

minmax.df <- function(data) {
  as.data.frame(apply(data, 2, minmax))
}

#' @title Select Hypervariable Features
#' @description This function selects hypervariable features from a data frame based on specified variance calculation methods.
#' @author Zaoqu Liu; E-mail: liuzaoqu@163.com
#' @param data A data frame with row features and column samples.
#' @param method Method of calculating variance ("sd", "mad", "cv").
#' @param top.percent The top percent of hypervariable features based on variance.
#' @param top.number The top number of hypervariable features based on variance.
#' @param custom.features A vector of custom features to select.
#' @return A vector of selected hypervariable features.
#' @export
#' @examples
#' \dontrun{
#' df <- data.frame(matrix(rnorm(1000), nrow = 100, ncol = 10))
#' selected_features <- Select.Features(df, method = "mad", top.percent = 10)
#' }
Select.Features <- function(data, method = "mad",
                            top.percent = NULL,
                            top.number = 1000,
                            custom.features = NULL) {
  if (!is.null(custom.features)) {
    ff <- intersect(rownames(data), custom.features)
  } else {
    if (tolower(method) == "sd") {
      vars <- SD.df(t(data))
    }
    if (tolower(method) == "mad") {
      vars <- MAD.df(t(data))
    }
    if (tolower(method) == "cv") {
      vars <- CV.df(t(data))
    }

    if (is.null(top.percent)) {
      if (top.number > length(vars)) {
        stop("The number of genes you specify must be smaller than the number of rows in the matrix!")
      } else {
        ff <- names(vars)[1:top.number]
      }
    } else {
      if (top.percent > 100) {
        stop("The top.percent must be smaller than 100 [1-100]!")
      } else {
        nn <- ceiling(nrow(data) * top.percent / 100)
        ff <- names(vars)[1:nn]
      }
    }
  }
  return(ff)
}
