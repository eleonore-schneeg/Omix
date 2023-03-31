################################################################################
#' Proteomics normalisation batch correction using Combat
#'
#' @param matrix protein abundance matrix
#' @param batch_map generated internally to `batch_correction_protein()`
#' @param batch Batch levels
#'
#' @return Combat corrected data frame
#'
#' @family Pre-processing
#'
#' @importFrom proBatch correct_with_ComBat_df quantile_normalize_dm
#' matrix_to_long long_to_matrix
#' @export

combat_correction <- function(matrix,
                              batch_map,
                              batch) {
  batch_map <- batch_map[colnames(matrix), ]
  log_matrix <- matrix
  median_normalized_matrix <- proBatch::quantile_normalize_dm(
    as.matrix(log_matrix)
  )
  median_normalized_omitna <- proBatch::matrix_to_long(
    na.omit(median_normalized_matrix)
  )
  median_normalized_omitna$Batch <- batch_map[[batch]][match(
    median_normalized_omitna$FullRunName, batch_map$FullRunName
  )]
  median_normalized_omitna$Batch <- as.factor(median_normalized_omitna$Batch)

  if (all(isTRUE(table(batch_map[[batch]]) > 1)) == FALSE) {
    stop(cli::cli_alert_danger(
      paste("One batch is present in a single sample, use median_centering
      parameter instead")
    ))
  }
  combat_df <- proBatch::correct_with_ComBat_df(
    median_normalized_omitna,
    sample_id_col = "FullRunName",
    batch_col = "Batch",
    feature_id_col = "peptide_group_label",
    measure_col = "Intensity"
  )
  combat <- proBatch::long_to_matrix(combat_df)
  combat <- as.data.frame(combat)
  return(combat)
}
