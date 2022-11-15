################################################################################
#' Proteomics normalisation batch correction using Combat
#'
#' @param matrix
#' @param batch_map
#'
#' @return
#' @importFrom proBatch correct_with_ComBat_df quantile_normalize_dm
#' matrix_to_long long_to_matrix
#' @export
#'
#' @examples
combat_correction <- function(matrix, batch_map) {
  batch_map <- batch_map[colnames(matrix), ]
  log_matrix <- matrix
  median_normalized_matrix <- proBatch::quantile_normalize_dm(as.matrix(log_matrix))
  median_normalized_omitna <- proBatch::matrix_to_long(na.omit(median_normalized_matrix))
  median_normalized_omitna$Batch <- batch_map[[paste(batch)]][match(median_normalized_omitna$FullRunName, batch_map$FullRunName)]
  comBat_df <- proBatch::correct_with_ComBat_df(median_normalized_omitna,
    sample_id_col = "FullRunName",
    batch_col = paste(batch)
  )
  combat <- proBatch::long_to_matrix(comBat_df)
  combat <- data.frame(combat)
  return(combat)
}
