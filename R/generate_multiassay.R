################################################################################
#' Generate a MultiAssayExperiment object from multiple single-omics data
#'
#' Generate a MultiAssayExperiment object from single-omics data matrix
#' and metadata. Currently supports transcriptomics and proteomics data only.
#'
#' @param rawdata_rna A data-frame containing raw RNA count where rows and
#' columns correspond to genes and samples respectively
#' @param rawdata_protein A data-frame containing raw protein abundance where
#' rows and columns correspond to genes and samples respectively
#' @param individual_to_sample Logical whether individual ID and sample names in
#' raw data matrices are the same. If they are different, `mpa_rna` and
#'  `map_protein`
#' dataframes should be provided. Default to `FALSE`.
#' @param map_rna A data-frame of two columns named `primary` and `colname`
#' where primary should contain unique ID name with a link to individual
#' metadata and colname is the column names of the `rawdata_rna` data-frame
#' @param map_protein A data-frame of two columns named `primary` and `colname`
#' where primary should contain unique ID name with a link to individual
#' metadata and colname is the column names of the `rawdata_protein` data-frame
#' @param metadata_rna A data-frame containing `rna` assay specific metadata
#' where rownames are same as the colnames of `rawdata_rna` data.frame
#' @param metadata_protein A data-frame containing `protein` assay specific
#' metadata where rownames are same as the colnames of `rawdata_protein`
#' data.frame
#' @param individual_metadata A data-frame containing sample level metadata
#' @param map_by_column The common column name to link `metadata_rna`
#' and `metadata_protein` to the `individual_metadata`
#' @param rna_qc_data Logical whether to add rna QC data.
#' @param rna_qc_data_matrix A data.frame containing RNA qc level data where
#' rownames are same as the colnames of `rawdata_rna` data.frame.
#' @param organism The organism type for transcripts and proteins.
#' Possible values are `human` and `mouse`. Default to `human`.
#'
#' @return a MultiAssayExperiment object
#'
#' @family Pre-processing
#'
#' @importFrom MultiAssayExperiment MultiAssayExperiment listToMap colData
#' @importFrom SummarizedExperiment SummarizedExperiment
#' @importFrom magrittr set_names
#' @importFrom cli cli_alert_danger style_bold cli_alert_success
#' @importFrom org.Mm.eg.db org.Mm.eg.db
#' @importFrom EnsDb.Hsapiens.v86 EnsDb.Hsapiens.v86
#' @export
generate_multiassay <- function(rawdata_rna,
                                rawdata_protein,
                                individual_to_sample = FALSE,
                                map_rna,
                                map_protein,
                                metadata_rna,
                                metadata_protein,
                                individual_metadata,
                                map_by_column,
                                rna_qc_data = FALSE,
                                rna_qc_data_matrix,
                                organism = "human") {
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
      paste("The columns in the", cli::style_bold("rawdata_protein"),
        "should be same as the rows in", cli::style_bold("metadata_protein"),
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


  individual_metadata <- as.data.frame(individual_metadata)

  if (individual_to_sample == TRUE) {
    map_rna <- data.frame(
      primary = colnames(rawdata_rna),
      colname = colnames(rawdata_rna),
      stringsAsFactors = FALSE
    )

    map_protein <- data.frame(
      primary = colnames(rawdata_protein),
      colname = colnames(rawdata_protein),
      stringsAsFactors = FALSE
    )
  }

  if (individual_to_sample == FALSE) {
    map_l <- list(
      data.frame(map_rna),
      data.frame(map_protein)
    )
  }

  names(map_l) <- c(
    "rna_raw",
    "protein_raw"
  )

  dfmap <- MultiAssayExperiment::listToMap(map_l)

  rawdata_rna <- rawdata_rna[!duplicated(rownames(rawdata_rna)), ]

  rna_id_type <- get_ID_type(rownames(rawdata_rna))

  se_rna <- SummarizedExperiment::SummarizedExperiment(
    assays = list("rna_raw" = rawdata_rna),
    colData = S4Vectors::DataFrame(metadata_rna),
    rowData = magrittr::set_names(
      S4Vectors::DataFrame(gene_id = as.character(rownames(rawdata_rna))),
      rna_id_type
    )
  )

  se_rna@metadata$metadata <- metadata_rna

  if (organism == "human") {
    EnsDb <- EnsDb.Hsapiens.v86::EnsDb.Hsapiens.v86
  }
  if (organism == "mouse") {
    EnsDb <- org.Mm.eg.db::org.Mm.eg.db
  }

  if (rna_id_type == "ensembl_gene_id") {
    # ID conversion using biomart
    cli::cli_alert_success("Ensembl ID conversion to gene symbol")
    cli::cli_alert_success("Retrieval of gene biotype")

    rna_values <- se_rna@elementMetadata@listData[[rna_id_type]]

    rna_df <- AnnotationDbi::select(EnsDb,
      keys = rna_values,
      columns = c("SYMBOL", "GENEBIOTYPE")
    )

    se_rna@elementMetadata@listData$gene_name <- rna_df$SYMBOL[match(
      rna_values, rna_df$GENEID
    )]


    se_rna@elementMetadata@listData$gene_biotype <- rna_df$GENEBIOTYPE[match(
      rna_values, rna_df$GENEID
    )]
  }

  if (rna_id_type == "gene_name") {
    rna_values <- se_rna@elementMetadata@listData[[rna_id_type]]
    cli::cli_alert_success("Retrieval of gene biotype")
    rna_df <- AnnotationDbi::select(EnsDb,
      keys = rna_values,
      columns = c("SYMBOL", "GENEBIOTYPE")
    )

    se_rna@elementMetadata@listData$gene_biotype <- rna_df$GENEBIOTYPE[match(
      rna_values, rna_df$GENEID
    )]
  }

  if (rna_qc_data) {
    se_rna@metadata$rna_qc_data <- rna_qc_data_matrix[colnames(rawdata_rna), ]
  }

  rawdata_protein <- rawdata_protein[!duplicated(rownames(rawdata_protein)), ]
  protein_id_type <- get_ID_type(rownames(rawdata_protein))

  se_protein <- SummarizedExperiment::SummarizedExperiment(
    assays = list("protein_raw" = rawdata_protein),
    colData = S4Vectors::DataFrame(metadata_protein),
    rowData = magrittr::set_names(
      S4Vectors::DataFrame(gene_id = as.character(rownames(rawdata_protein))),
      protein_id_type
    )
  )

  se_protein@metadata$metadata <- metadata_protein

  if (protein_id_type == "uniprot_id") {
    # ID conversion using biomart
    cli::cli_alert_success("UniProt ID conversion to gene name ")
    protein_values <- se_protein@elementMetadata@listData[[protein_id_type]]
    protein_values <- sub("\\-.*", "", protein_values)


    protein_df <- AnnotationDbi::select(EnsDb,
      keys = protein_values,
      keytype = "UNIPROTID",
      columns = c("SYMBOL")
    )

    se_protein@elementMetadata@listData$gene_name <- protein_df$SYMBOL[match(
      protein_values,
      protein_df$UNIPROTID
    )]
  }


  rownames(individual_metadata) <- individual_metadata[, map_by_column]

  multiassay <- MultiAssayExperiment::MultiAssayExperiment(
    list(
      "rna_raw" = se_rna,
      "protein_raw" = se_protein
    ),
    colData = S4Vectors::DataFrame(individual_metadata),
    sampleMap = S4Vectors::DataFrame(dfmap)
  )

  cli::cli_alert_success("RNA raw data loaded")
  print(se_rna)
  cli::cli_alert_success("Protein raw data loaded")
  print(se_protein)
  cli::cli_alert_success("MultiAssayExperiment object generated!")

  return(multiassay)
}
