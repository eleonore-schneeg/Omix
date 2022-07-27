################################################################################
#' Generate a MultiAssayExperiment object from single-omics data
#'
#' Generate a MultiAssayExperiment object from single-omics data matrix
#' and metadata. Currently supports transcriptomics and proteomics data only.
#'
#' @param rawdata_rna A data-frame containing raw RNA count where rows and
#' columns correspond to genes and samples respectively
#' @param rawdata_protein A data-frame containing raw protein abundance where
#' rows and columns correspond to genes and samples respectively
#' @param map_rna A data-frame of two columns named `primary` and `colname`
#' where primary should contain unique sample name with a link to sample
#' metadata and colname is the column names of the `rawdata_rna` data-frame
#' @param map_protein A data-frame of two columns named `primary` and `colname`
#' where primary should contain unique sample name with a link to sample
#' metadata and colname is the column names of the `rawdata_protein` data-frame
#' @param metadata_rna A data-frame containing `rna` assay specific metadata
#' where rownames are same as the colnames of `rawdata_rna` data.frame
#' @param metadata_protein A data-frame containing `protein` assay specific
#' metadata where rownames are same as the colnames of `rawdata_protein`
#' data.frame
#' @param sample_metadata A data-frame containing sample level metadata for
#' both omics assays
#' @param map_by_column The common column name to link `metadata_rna`
#' and `metadata_protein` to the `sample_metadata`
#' @param rna_id_type The gene ID type used for the rownames of `rawdata_rna`.
#' Possible values are `ensembl_gene_id`, `gene_name`.
#' @param protein_id_type The protein ID type used for the rownames of
#' `rawdata_protein`. Possible values are `uniprot_id`, `gene_name`.
#' @param rna_qc_data Logical whether to add rna QC data.
#' @param rna_qc_data_matrix A data.frame containing RNA qc level data where
#' rownames are same as the colnames of `rawdata_rna` data.frame.
#'
#' @return a MultiAssayExperiment object
#'
#' @importFrom MultiAssayExperiment MultiAssayExperiment listToMap colData
#' @importFrom SummarizedExperiment SummarizedExperiment
#' @importFrom magrittr set_names
#' @importFrom cli cli_alert_danger style_bold cli_alert_success
#'
#' @export
generate_mae <- function(rawdata_rna,
                         rawdata_protein,
                         map_rna,
                         map_protein,
                         metadata_rna,
                         metadata_protein,
                         sample_metadata,
                         map_by_column,
                         rna_id_type,
                         protein_id_type,
                         rna_qc_data = TRUE,
                         rna_qc_data_matrix) {
  if (!all(colnames(rawdata_rna) == rownames(metadata_rna))) {
    stop(cli::cli_alert_danger(
      paste("The columns in the", cli::style_bold("rawdata_rna"),
            "should be same as the rows in", cli::style_bold("metadata_rna"),
            sep = " "
      )
    ))
  }

  if (!all(colnames(rawdata_protein) == rownames(metadata_protein))) {
    stop(cli::cli_alert_danger(
      paste("The columns in the", cli::style_bold("rawdata_rna"),
            "should be same as the rows in", cli::style_bold("metadata_rna"),
            sep = " "
      )
    ))
  }

  if (!all(colnames(map_rna) == c("primary", "colname"))) {
    stop(cli::cli_alert_danger(
      paste(cli::style_bold("primary"), "and", cli::style_bold("colname"),
            "columns are not found in", cli::style_bold("map_rna!"),
            sep = " "
      )
    ))
  }

  if (!all(colnames(map_protein) == c("primary", "colname"))) {
    stop(cli::cli_alert_danger(
      paste(cli::style_bold("primary"), "and", cli::style_bold("colname"),
            "columns are not found in", cli::style_bold("map_protein!"),
            sep = " "
      )
    ))
  }

  if (!all(colnames(rawdata_rna) == rownames(rna_qc_data_matrix))) {
    stop(cli::cli_alert_danger(
      paste("The columns in the", cli::style_bold("rawdata_rna"),
            "should be same as the rows in", cli::style_bold("rna_qc_data_matrix"),
            sep = " "
      )
    ))
  }


  map_l <- list(
    map_rna,
    map_protein
  )

  names(map_l) <- c(
    "rna_raw",
    "protein_raw"
  )

  dfmap <- MultiAssayExperiment::listToMap(map_l)

  rawdata_rna <- rawdata_rna[!duplicated(rownames(rawdata_rna)), ]

  se_rna <- SummarizedExperiment::SummarizedExperiment(
    assays = list("rna_raw" = rawdata_rna),
    colData = metadata_rna,
    rowData = magrittr::set_names(
      data.frame(gene_id = as.character(rownames(rawdata_rna))),
      rna_id_type
    )
  )

  se_rna@metadata$metadata <- metadata_rna

  if (rna_qc_data) {
    se_rna@metadata$rna_qc_data <- rna_qc_data_matrix[colnames(rawdata_rna), ]
  }

  rawdata_protein <- rawdata_protein[!duplicated(rownames(rawdata_protein)), ]

  se_protein <- SummarizedExperiment::SummarizedExperiment(
    assays = list("protein_raw" = rawdata_protein),
    colData = metadata_protein,
    rowData = magrittr::set_names(
      data.frame(gene_id = as.character(rownames(rawdata_protein))),
      protein_id_type
    )
  )

  se_protein@metadata$metadata <- metadata_protein

  MultiAssay <- MultiAssayExperiment::MultiAssayExperiment(
    list(
      "rna_raw" = se_rna,
      "protein_raw" = se_protein
    ),
    colData = sample_metadata,
    sampleMap = dfmap
  )

  cli::cli_alert_success("MultiAssayExperiment object generated")

  return(MultiAssay)
}
