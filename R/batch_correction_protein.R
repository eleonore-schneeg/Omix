################################################################################
#' Calling different proteomics correction methods.
#' Internal to `process_protein()`
#'
#' @param matrix Protein abundance matrix
#' @param batch_map Automatically generated in `process_protein()`
#' @param batch Batch column name
#' @param correction_method `Combat` or `Median_centering`
#'
#' @return Batch corrected matrix
#'
#' @family Pre-processing
#'
#' @export

batch_correction_protein <- function(matrix,
                                     batch_map,
                                     batch,
                                     correction_method) {
  args <- list(
    matrix = matrix,
    batch_map = batch_map,
    batch = batch
  )
  if (correction_method == "Combat") {
    cli::cli_alert_success("Performing combat batch correction")
    matrix <- do.call(combat_correction, args)
  }
  if (correction_method == "Median_centering") {
    cli::cli_alert_success("Performing median centering batch correction")
    matrix <- do.call(median_centering_correction, args)
  }
  return(matrix)
}
