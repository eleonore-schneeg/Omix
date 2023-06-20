#' Integrated pathway enrichment with multiGSEA
#'
#' @param integration Options are pvalue padj weights
#' @param DE_comparison_dataframe if integration == pvalue or padj, provide the differential expression comparison dataframe generated with `single_omic_comparisons()`
#' @param weights if integration == weights, provide weights dataframe generated with `extract_weights()`
#' @param databases Options are: all kegg panther pathbank pharmgkb reactome smpdb wikipathways or a combination thereof. Use NULL for GMT file only
#' @param additional_database_gmt provide GMT file path
#'
#' @return List of results
#'
#' @family Enrichment analysis
#'
#' @importFrom multiGSEA initOmicsDataStructure getMultiOmicsFeatures multiGSEA extractPvalues combinePvalues
#' @importFrom stats p.adjust

#' @export
#'
integrated_enrichment_mutliGSEA <- function(integration = "pvalue",
                                            DE_comparison_dataframe = NULL,
                                            weights = NULL,
                                            databases = 'reactome',
                                            additional_database_gmt = FALSE) {
  # create data structure
  omics_data <- multiGSEA::initOmicsDataStructure(layer = c(
    "transcriptome",
    "proteome"
  ))


  if (integration == "pvalue") {
    scores_group <- DE_comparison_dataframe
    scores_group <- scores_group[!duplicated(scores_group$gene_name), ]
    rownames(scores_group)=scores_group$gene_name

    #rankFeatures calculates the a local statistic ls based on the direction of the fold change and the magnitude of its significance
    omics_data$transcriptome <- multiGSEA::rankFeatures(
      scores_group$log2FoldChange_transcriptomics,
      -log10(scores_group$pvalue_transcriptomics)
    )
    keep=!is.infinite(omics_data$transcriptome)
    omics_data$transcriptome <- omics_data$transcriptome[!is.infinite(omics_data$transcriptome)]
    names(omics_data$transcriptome) <- scores_group$gene_name[keep]

    omics_data$proteome <- multiGSEA::rankFeatures(
      scores_group$log2FoldChange_proteomics,
      -log10(scores_group$pvalue_proteomics)
    )
    keep=!is.infinite(omics_data$proteome)
    omics_data$proteome <- omics_data$proteome[!is.infinite(omics_data$proteome)]
    names(omics_data$proteome) <- scores_group$gene_name[keep]

  }

  if (integration == "padj") {
    scores_group <- DE_comparison_dataframe
    scores_group <- scores_group[!duplicated(scores_group$gene_name), ]
    rownames(scores_group)=scores_group$gene_name


    ## add transcriptome layer
    omics_data$transcriptome <- multiGSEA::rankFeatures(
      scores_group$log2FoldChange_transcriptomics,
      scores_group$padj_transcriptomics
    )
    keep=!is.infinite(omics_data$transcriptome)
    omics_data$transcriptome <- omics_data$transcriptome[!is.infinite(omics_data$transcriptome)]
    names(omics_data$transcriptome) <- scores_group$gene_name[keep]

    omics_data$proteome <- multiGSEA::rankFeatures(
      scores_group$log2FoldChange_proteomics,
      scores_group$padj_proteomics
    )
    keep=!is.infinite(omics_data$proteome)
    omics_data$proteome <- omics_data$proteome[!is.infinite(omics_data$proteome)]
    names(omics_data$proteome) <- scores_group$gene_name[keep]

  }

  if (integration == "weights") {
    ## add transcriptome layer
    rna <- weights$weights_df$rna$Weights
    names(rna) <- weights$weights_df$rna$Feature
    sort <- sort(abs(rna), decreasing = T, index.return = T)
    rna <- rna[sort$ix]
    names(rna) <- names(rna[sort$ix])
    omics_data$transcriptome <- rna

    ## add proteome layer
    prot <- weights$weights_df$protein$Weights
    names(prot) <- weights$weights_df$protein$Feature
    sort <- sort(abs(prot), decreasing = T, index.return = T)
    prot <- prot[sort$ix]
    names(prot) <- names(prot[sort$ix])
    omics_data$proteome <- prot
  }

if(!is.null(databases)){
  layers <- names(omics_data)
  pathways <- multiGSEA::getMultiOmicsFeatures(
    dbs = databases, layer = layers,
    returnTranscriptome = "SYMBOL",
    returnProteome = "SYMBOL"
  )
}

  if (!is.null(additional_database_gmt)) {

    myGO <- ActivePathways::read.GMT(additional_database_gmt)
    myGO <- lapply(myGO , function(x) {
      x$genes
    })

    pathways_GO <- list()
    pathways_GO$transcriptome <- myGO
    pathways_GO$proteome <- myGO
    if(!is.null(databases)){
      pathways <- Map(c, pathways, pathways_GO)
    }
    if(is.null(databases)){
      pathways <- pathways_GO
    }

  }

  enrichment_scores <- multiGSEA::multiGSEA(pathways, omics_data)

  df <- multiGSEA::extractPvalues(
    enrichmentScores = enrichment_scores,
    pathwayNames = names(pathways[[1]])
  )

  df$combined_pval <- multiGSEA::combinePvalues(df, method = "fisher")
  df$combined_padj <- stats::p.adjust(df$combined_pval, method = "BH")

  df <- cbind(data.frame(pathway = names(pathways[[1]])), df)

  return(df)
}
