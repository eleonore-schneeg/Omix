#' Data filtering function
#'
#' @param matrix data frame containing LOG2 data
#' @param groups A character vector dictating the grouping
#' @param min_sample minimum % of samples that should have non missing
#' protein value, or else the protein is excluded
#' @param at_least_one TRUE means to keep the row if min_count is met for at
#' least one condition, FALSE means min_count must be met across all
#' groups for retention
#'
#' @return Filtered protein matrix
#' @export
#'

filter_protein <- function(matrix,
                           groups,
                           min_sample = 0.5,
                           at_least_one = FALSE) {
  df <- matrix
  cond.names <- list()

  if (!is.null(groups)) {
    for (group in levels(groups)) {
      cond.names[[group]] <- colnames(matrix)[groups == group]
    }

    cond.filter <- sapply(1:length(cond.names), function(i) {
      df2 <- df[, cond.names[[i]]]
      df2 <- as.matrix(df2)
      sums <- rowSums(is.finite(df2))
      sums >= min_sample * length(cond.names[[i]])
    })

    if (at_least_one) {
      KEEP <- apply(cond.filter, 1, any)
    } else {
      KEEP <- apply(cond.filter, 1, all)
    }
    df <- df[KEEP, ]
  }

  if (is.null(groups)) {
    idx <- rowSums(is.finite(as.matrix(df))) >= dim(df)[2] * min_sample
    df <- df[idx, ]
  }

  return(df)
}
