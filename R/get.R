################################################################################
#' Title
#'
#' @param multiassay
#' @param omic
#' @param slot
#'
#' @return
#' @export
#'
#' @examples
.get_metadata <- function(multiassay,
                          omic,
                          slot=c('rna_processed','protein_processed')) {
  data <- multiassay@ExperimentList@listData[[paste(slot)]]

  if (omic == "protein") {
    raw <- MultiAssayExperiment::getWithColData(multiassay,
                                                i = "protein_raw",
                                                mode = "replace",
                                                verbose = F)
  }
  if (omic == "rna") {
    raw <- MultiAssayExperiment::getWithColData(multiassay,
                                                i = "rna_raw",
                                                mode = "replace",
                                                verbose = F)
  }
  colData <- data.frame(SummarizedExperiment::colData(raw))
  if (omic == "protein") {
    order <- map_protein$primary[match(colnames(raw), map_protein$colname)]
  }

  if (omic == "rna") {
    order <- map_rna$primary[match(colnames(raw), map_rna$colname)]
  }

  colData <- colData[order, ]

  return(colData)
}

#' Title
#'
#' @param id
#' @param omic
#' @param from
#' @param to
#'
#' @return
#' @export
#'
#' @examples
.get_ID_names <- function(id,
                          omic = "protein",
                          from = "uniprot_id",
                          to = "gene_name") {
  if (omic == "protein") {
    elementMetadata <- data.frame(multiassay@ExperimentList@listData[["protein_raw"]]@elementMetadata@listData)
    new_id <- elementMetadata[[paste(to)]][match(id, elementMetadata[[paste(from)]])]
    new_id[is.na(new_id)] <- id[is.na(new_id)]
  }
  if (omic == "rna") {
    elementMetadata <- data.frame(multiassay@ExperimentList@listData[["rna_raw"]]@elementMetadata@listData)
    new_id <- elementMetadata[[paste(to)]][match(id, elementMetadata[[paste(from)]])]
    new_id[is.na(new_id)] <- id[is.na(new_id)]
  }

  return(new_id)
}

.get_background <- function(multiassay,
                            of=c('full','protein','rna')) {
  if(of=='full'){
  background <- union(
    multiassay@ExperimentList@listData[["rna_raw"]]@elementMetadata@listData[["gene_name"]],
    multiassay@ExperimentList@listData[["protein_raw"]]@elementMetadata@listData[["gene_name"]]
  )
  }

  if(of=='rna'){
    background=multiassay@ExperimentList@listData[["rna_raw"]]@elementMetadata@listData[["gene_name"]]

  }
  if(of=='protein'){
  background= multiassay@ExperimentList@listData[["protein_raw"]]@elementMetadata@listData[["gene_name"]]

  }

  return(background)
}


.get_ID_type <- function(character_vector){

  test=character_vector[1]

  if(isTRUE(grepl('^[A-Z0-9-]+$|^C[0-9XY]+orf[0-9]+$',
           test))){
    ID_type='gene_name'
  }

  if(isTRUE(grepl('ENS',test))){
    ID_type='ensembl_gene_id'
  }

  if(isTRUE(!grepl("\\D", test))){
    ID_type='entrez_gene_id'
  }

  if(isTRUE(grepl('[OPQ][0-9][A-Z0-9]{3}[0-9]|[A-NR-Z][0-9]([A-Z][A-Z0-9]{2}[0-9]){1,2}',
        test))){
    ID_type='uniprot_id'
  }

  if(isTRUE(!exists('ID_type'))){
    ID_type='not_identified'
    cli::cli_warn('ID type not identified!')
  }

  return(ID_type)
}


.get_multimodal_object <- function(multiassay,
                                   slots = c(
                                     "rna_processed",
                                     "protein_processed"
                                   ),
                                   intersect_genes=FALSE,
                                   ID_type=c('gene_name','original')
                                   ){

  multi <- multiassay[, , c(slots[1], slots[2])]
  complete <- complete.cases(multi)
  CompleteMulti <- multi[, complete.cases(multi), ]

  if(intersect_genes==TRUE){

    rownames(CompleteMulti@ExperimentList@listData[["rna_processed"]])=
      make.unique(.get_ID_names(rownames(
        CompleteMulti@ExperimentList@listData[["rna_processed"]]),
                  omic = "rna",
                  from = "ensembl_gene_id",
                  to = "gene_name"
    ))

    rownames(CompleteMulti@ExperimentList@listData[["protein_processed"]])=
      make.unique(.get_ID_names(rownames(
      CompleteMulti@ExperimentList@listData[["protein_processed"]]),
                  omic = "protein",
                  from = "uniprot_id",
                  to = "gene_name"
    ))

    CompleteMulti<- intersectRows(CompleteMulti[, ,slots])

  }

  metadata <- data.frame(colData(CompleteMulti))
  multimodal_omics <- lapply(CompleteMulti@ExperimentList@listData, as.matrix)
  multimodal_omics <- lapply(multimodal_omics, "colnames<-", rownames(metadata))

  if(ID_type=='gene_name'){

    rownames(multimodal_omics[[1]])=
      make.unique(.get_ID_names(rownames(
        CompleteMulti@ExperimentList@listData[["rna_processed"]]),
        omic = "rna",
        from = "ensembl_gene_id",
        to = "gene_name"
      ))

    rownames(multimodal_omics[[2]])=
      make.unique(.get_ID_names(rownames(
      CompleteMulti@ExperimentList@listData[["protein_processed"]]),
      omic = "protein",
      from = "uniprot_id",
      to = "gene_name"
    ))
  }


  return(list(multimodal_object=    multimodal_omics,
              metadata=   metadata))

}

