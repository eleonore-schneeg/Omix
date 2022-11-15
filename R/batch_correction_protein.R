################################################################################
#' Calling different proteomics correction methods
#'
#' @param matrix
#' @param batch_map
#' @param batch
#' @param correction_method
#'
#' @return
#' @export
#'
#' @examples
batch_correction_protein <- function(matrix,
                                     batch_map,
                                     batch,
                                     correction_method) {
  args <- list(matrix = matrix,
               batch_map = batch_map,
               batch = batch)
  if (correction_method == "Combat") {
    cli::cli_alert_success("Performing combat batch correction")
    res <- do.call(combat_correction, args)
  }
  if (correction_method == "Median_centering") {
    cli::cli_alert_success("Performing median centering batch correction")
    res <- do.call(median_centering_correction, args)
  }
  return(res)
}
