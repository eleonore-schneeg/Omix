################################################################################
#' Title
#'
#' @param multiassay
#' @param integration
#' @param component
#'
#' @return
#' @export
#'
#' @examples
integrative_results_supervised <- function(multiassay,
                                integration = "DIABLO",
                                component = 1,
                                correlation_threshold=0.5,
                                disease_id='MONDO_0004975') {
  if (integration == "DIABLO" | integration == "sMBPLS") {
    model <- multiassay@metadata$integration[[paste(integration)]]

    matrix <- cbind(model$X[[1]], model$X[[2]])

    loadings=mixOmics::selectVar(model, comp=component)
    Loadings <- c(
      loadings[["mRNA"]][["name"]],
      loadings[["proteins"]][["name"]]
    )

    df1=tibble::rownames_to_column(loadings$mRNA$value, "feature")
    df1$`Omic.layer`='rna'
    df2=tibble::rownames_to_column(loadings$proteins$value, "feature")
    df2$`Omic.layer`='protein'

    data_loadings=rbind(df1,df2)

    if (integration == "sMBPLS"){
    if(model[["loadings"]][["Y"]][1]==1){

      data_loadings_positive=data_loadings[data_loadings$value.var >0,]
      data_loadings_negative=data_loadings[data_loadings$value.var <0,]
    }
    if(model[["loadings"]][["Y"]][1]!=1){
      data_loadings_positive=data_loadings[data_loadings$value.var <0,]
      data_loadings_negative=data_loadings[data_loadings$value.var >0,]
    }
    }

    if (integration == "DIABLO"){
    if (model[["loadings"]][["Y"]][1] > 0) {
      data_loadings_positive <- data_loadings[data_loadings$value.var > 0, ]
      data_loadings_negative <- data_loadings[data_loadings$value.var < 0, ]
    }
    if (model[["loadings"]][["Y"]][1] < 0) {
      data_loadings_positive <- data_loadings[data_loadings$value.var < 0, ]
      data_loadings_negative <- data_loadings[data_loadings$value.var > 0, ]
    }
    }

    cli::cli_alert_success("RETRIEVAL OF PROTECTIVE AND DETRIMENTAL MULTI-OMICS SIGNATURE")
    rna_positive <- data_loadings_positive$feature[data_loadings_positive$Omic.layer == "rna"]
    rna_negative <- data_loadings_negative$feature[data_loadings_negative$Omic.layer == "rna"]
    protein_positive <- data_loadings_positive$feature[data_loadings_positive$Omic.layer == "protein"]
    protein_negative <- data_loadings_negative$feature[data_loadings_negative$Omic.layer == "protein"]

    detrimental <- list(
      rna = rna_positive,
      protein = protein_positive
    )

    protective <- list(
      rna = rna_negative,
      protein = protein_negative
    )

    cli::cli_alert_success("BIPARTITE NETWORKS GENERATION")
    protective_network <- .multiomics_network(multiassay,
      list = protective,
      correlation_threshold = correlation_threshold
    )

    detrimental_network <- .multiomics_network(multiassay,
      list = detrimental,
      correlation_threshold = correlation_threshold
    )


    ### ID MAPPING
    cli::cli_alert_success("ID MAPPING")
    rna_positive <- .get_ID_names(rna_positive,
      omic = "rna",
      from = "ensembl_gene_id",
      to = "gene_name"
    )
    rna_negative <- .get_ID_names(rna_negative,
      omic = "rna",
      from = "ensembl_gene_id",
      to = "gene_name"
    )
    protein_positive <- .get_ID_names(protein_positive,
      omic = "protein",
      from = "uniprot_id",
      to = "gene_name"
    )
    protein_negative <- .get_ID_names(protein_negative,
      omic = "protein",
      from = "uniprot_id",
      to = "gene_name"
    )

    detrimental <- list(
      rna = rna_positive,
      protein = protein_positive
    )
    protective <- list(
      rna = rna_negative,
      protein = protein_negative
    )

    cli::cli_alert_success("DETECTION OF NETWORK COMMUNITIES")
    detrimental_communities <- .communities_network(detrimental_network)

    protective_communities <- .communities_network(protective_network)


    cli::cli_alert_success("CELL TYPE ENRICHMENT OF COMMUNITIES")

    detrimental_cell_type <- cell_type_enrichment(
      multiassay = multiassay,
      communities = detrimental_communities
    )

    protective_cell_type <- cell_type_enrichment(
      multiassay = multiassay,
      communities = protective_communities
    )



    cli::cli_alert_success("FUNCTIONAL ENRICHMENT")
    background=.get_background(multiassay)

    enrichr_detrimental=lapply(detrimental,pathway_analysis_enrichr)
    enrichr_protective=lapply(protective,pathway_analysis_enrichr)

    library(org.Hs.eg.db)
    enrichGO_detrimental=lapply(detrimental,
                                function(x){clusterProfiler::enrichGO(gene = x,
                                             OrgDb         = org.Hs.eg.db,
                                             universe = background,
                                             keyType       = 'SYMBOL',
                                             ont           = "BP",
                                             pAdjustMethod = "BH",
                                             pvalueCutoff  = 0.05,
                                             qvalueCutoff  = 0.05,
                                             minGSSize = 2,
                                             maxGSSize = 500) })

    enrichGO_protective=lapply(protective,function(x){clusterProfiler::enrichGO(gene = x,
                                                      OrgDb         = org.Hs.eg.db,
                                                      universe = background,
                                                      keyType       = 'SYMBOL',
                                                      ont           = "BP",
                                                      pAdjustMethod = "BH",
                                                      pvalueCutoff  = 0.05,
                                                      qvalueCutoff  = 0.05,
                                                      minGSSize = 2,
                                                      maxGSSize = 500) })


   #### Communities comparisons / mix proteomics transcriptomics if highly corr
    detrimental_comp <- clusterProfiler::compareCluster(geneCluster =
                         detrimental_communities,
                         fun="enrichGO",
                         OrgDb='org.Hs.eg.db',
                         keyType="SYMBOL",
                         universe=background,
                         ont = "BP",
                         minGSSize=2,
                         pAdjustMethod = "BH",
                         pvalueCutoff  = 0.05,
                         qvalueCutoff  = 0.05)


    detrimental_comp <- enrichplot::pairwise_termsim(detrimental_comp)

    protective_comp <- clusterProfiler::compareCluster(geneCluster =
                       protective_communities,
                       fun="enrichGO",
                       OrgDb='org.Hs.eg.db',
                       keyType="SYMBOL",
                       universe=background,
                       ont = "BP",
                       minGSSize=2,
                       pAdjustMethod = "BH",
                       pvalueCutoff  = 0.05,
                       qvalueCutoff  = 1)

    protective_comp <- enrichplot::pairwise_termsim(protective_comp)


    ### OPEN TARGETS
    cli::cli_alert_success("TESTING MULTI_OMICS SIGNATURE ASSOCIATION
                           WITH DISEASE OF INTEREST VIA THE OPENTARGETS API")

    opentarget=get_diseaseAssociations_df(disease_id=disease_id,size= 3000)
    opentarget_detrimental=plot_filter_OpenTarget(opentarget_results=
                           opentarget,genes=unlist(detrimental))
    opentarget_protective=plot_filter_OpenTarget(opentarget_results=
                           opentarget,genes=unlist(protective))



    #### RESULTS OBJECT
    integrative_results= list(detrimental=list( multiomics_signature=detrimental,
                                                network=detrimental_network,
                                                network_communities=detrimental_communities,
                                                cell_type_communities=detrimental_cell_type,
                                                enrichment=enrichr_detrimental,
                                                open_targets=    opentarget_detrimental),
                              protective=list(multiomics_signature=protective,
                                              network=protective_network,
                                              network_communities=protective_communities,
                                              cell_type_communities=protective_cell_type,
                                              enrichment=enrichr_protective,
                                              open_targets=    opentarget_protective),
                              comparison=list(detrimental=    detrimental_comp,
                                              protective=    protective_comp ))


return(integrative_results)
  }
}
