#' Adds a custom processed transcriptomics slot to multiassay object
#'
#' @param multiassay
#' @param custom_processed_df
#'
#' @return
#' @export
#'
#' @examples
processed_rna <- function(multiassay,
                          custom_processed_df) {
  matrix <- as.data.frame(custom_processed_df)
  map <- MultiAssayExperiment::sampleMap(multiassay)
  map_df <- data.frame(map@listData)
  map_df <- map_df[which(map_df$assay == "rna_raw"), ]
  map_df <- map_df[match(
    colnames(matrix),
    map_df$colname
  ), ][c("primary", "colname")]
  map_df$assay <- "rna_processed"

  rm.rna_processed <- !grepl("rna_processed", names(experiments(multiassay)))
  multiassay <- multiassay[, , rm.rna_processed]
  multiassay <- c(multiassay,
                  rna_processed = matrix,
                  sampleMap = map_df
  )

  metadata(multiassay)$parameters_processing_rna$custom_processed_df <- TRUE

  cli::cli_alert_success("Custom processed transcriptomics added!")
}



#' Adds custom processed dataset to multiassay object
#'
#' @param multiassay
#' @param custom_processed_df
#'
#' @return
#' @export
#'
#' @examples

counts_rna <- function(multiassay,
                       custom_processed_df) {
  matrix <- as.data.frame(custom_processed_df)
  map <- MultiAssayExperiment::sampleMap(multiassay)
  map_df <- data.frame(map@listData)
  map_df <- map_df[which(map_df$assay == "rna_raw"), ]
  map_df <- map_df[match(
    colnames(matrix),
    map_df$colname
  ), ][c("primary", "colname")]
  map_df$assay <- "rna_counts"

  multiassay_temp <- c(multiassay,
                       rna_counts = matrix,
                       sampleMap = map_df
  )

  metadata(multiassay_temp)$parameters_processing_rna$counts_df <- TRUE
  return(multiassay_temp)

  cli::cli_alert_success("Raw transcriptomics counts added!")
}
