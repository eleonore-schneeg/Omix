#' Imputation of missing values based on distribution
#'
# Since missing values are associated with proteins with low levels of expression,
# we can substitute the missing values with numbers that are considered “small”
# in each sample. We can define this statistically by drawing from a normal
# distribution with a mean that is down-shifted from the sample mean and a
# standard deviation that is a fraction of the standard deviation of
# the sample distribution.
#' @param df Data frame to be imputed
#' @param width Coefficient shrinking standard deviation with. Default to 0.3
#' @param downshift Coefficient shifting the mean of imputed values. Default to 1.8
#'
#' @return Imputed data frame
#' @export

.impute_distribution <- function(df, width = 0.3, downshift = 1.8) {
  for (i in colnames(df)) {
    non_missing <- !is.na(df[[i]])
    temp.sd <- width * sd(df[[i]][non_missing], na.rm = TRUE) # shrink sd width
    temp.mean <- mean(df[[i]][non_missing], na.rm = TRUE)
    downshift * sd(df[[i]], na.rm = TRUE) # shift mean of imputed values
    n.missing <- sum(is.na(df[[i]]))
    df[[i]][is.na(df[[i]])] <- rnorm(n.missing, mean = temp.mean, sd = temp.sd)
  }
  return(df)
}

#' Imputation of missing values based on 50% minimum value

#' @param df Data frame to be imputed
#'
#' @return Imputed data frame
#'
#' @importFrom matrixStats rowMins
#' @export
#'

.impute_minimum_value <- function(df) {
  for (i in colnames(df)) {
    protein_minimum_imputation <- 0.5 * matrixStats::rowMins(as.matrix(df), na.rm = T)
    missing <- is.na(df[[i]])
    df[[i]][missing] <- protein_minimum_imputation[missing]
  }
  return(df)
}



