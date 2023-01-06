################################################################################
#' Proteomics normalisation batch correction using median centering method
#'
#' @param matrix
#' @param batch_map
#' @param batch
#'
#' @return
#' @importFrom proBatch center_feature_batch_medians_df quantile_normalize_dm
#' matrix_to_long long_to_matrix
#' @export
#'
#' @examples
median_centering_correction <- function(matrix, batch_map, batch) {
  batch_map <- batch_map[colnames(matrix), ]
  log_matrix <- matrix
  median_normalized_matrix <- proBatch::quantile_normalize_dm(as.matrix(log_matrix))
  median_normalized_plot <- proBatch::matrix_to_long(median_normalized_matrix)
  median_normalized_plot$Batch <- batch_map[[paste(batch)]][match(median_normalized_plot$FullRunName, batch_map$FullRunName)]
  featurelevel_centering <- proBatch::center_feature_batch_medians_df(median_normalized_plot,
    sample_id_col = "FullRunName",
    batch_col = paste(batch)
  )
  median_center <- proBatch::long_to_matrix(featurelevel_centering)
  median_center <- as.data.frame(median_center)
  return(median_center)
}
