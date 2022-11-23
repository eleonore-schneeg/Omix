################################################################################
#' Downstream analysis for integrated multi-omics data
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
                                           correlation_threshold = 0.2,
                                           disease_id = "MONDO_0004975",
                                           enrichment_method='enrichr') {
  cli::cli_alert_success("RETRIEVAL OF PROTECTIVE AND DETRIMENTAL MULTI-OMICS SIGNATURE")
  signature <- extract_multiomic_signature(multiassay,
    integration = integration,
    component = component
  )

  protective <- signature$protective
  detrimental <- signature$detrimental

  cli::cli_alert_success("BIPARTITE NETWORKS GENERATION")
  protective_network <- .multiomics_network(multiassay,
    list = protective,
    correlation_threshold = correlation_threshold
  )

  detrimental_network <- .multiomics_network(multiassay,
    list = detrimental,
    correlation_threshold = correlation_threshold
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

  protective=lapply(protective, function(x) {
    sub("*\\.[0-9]", "", x)
  })

  detrimental=lapply(detrimental, function(x) {
    sub("*\\.[0-9]", "", x)
  })


  background <- .get_background(multiassay, of = "full")

  if(enrichment_method=='enrichr'){
  enrichr_detrimental <- lapply(detrimental, pathway_analysis_enrichr)
  enrichr_protective <- lapply(protective, pathway_analysis_enrichr)
  }


  if(enrichment_method=='enrichGO'){
  library(org.Hs.eg.db)
  enrichGO_detrimental <- lapply(detrimental, function(x) {
      clusterProfiler::enrichGO(
        gene = x,
        OrgDb = org.Hs.eg.db,
        universe = background,
        keyType = "SYMBOL",
        ont = "BP",
        pAdjustMethod = "BH",
        pvalueCutoff = 0.05,
        qvalueCutoff = 0.05,
        minGSSize = 2,
        maxGSSize = 500
      )
    })

  enrichGO_protective <- lapply(protective, function(x) {
    clusterProfiler::enrichGO(
      gene = x,
      OrgDb = org.Hs.eg.db,
      universe = background,
      keyType = "SYMBOL",
      ont = "BP",
      pAdjustMethod = "BH",
      pvalueCutoff = 0.05,
      qvalueCutoff = 0.05,
      minGSSize = 2,
      maxGSSize = 500
    )
  })
}


  # Comparison of communities
  cli::cli_alert_success("FUNCTIONALLY COMPARING MULTIOMICS NETWORKS COMMUNITIES")

  detrimental_communities_x <- lapply(detrimental_communities, function(x) x[length(x) >= 2])
  detrimental_communities_x  <-   detrimental_communities[lapply(detrimental_communities , length) > 0]

  detrimental_comp <- clusterProfiler::compareCluster(
    geneCluster = detrimental_communities_x,
    fun = "enrichGO",
    OrgDb = "org.Hs.eg.db",
    keyType = "SYMBOL",
    universe = background,
    ont = "BP",
    minGSSize = 2,
    pAdjustMethod = "BH",
    pvalueCutoff = 0.05,
    qvalueCutoff = 0.05
  )

  detrimental_comp <- enrichplot::pairwise_termsim(detrimental_comp)

  protective_communities_x <- lapply(protective_communities, function(x) x[length(x) >= 2])
  protective_communities_x  <-   protective_communities[lapply(protective_communities, length) > 0]

  protective_comp <- clusterProfiler::compareCluster(
    geneCluster = protective_communities_x,
    fun = "enrichGO",
    OrgDb = "org.Hs.eg.db",
    keyType = "SYMBOL",
    universe = background,
    ont = "BP",
    minGSSize = 2,
    pAdjustMethod = "BH",
    pvalueCutoff = 0.05,
    qvalueCutoff = 1
  )

  protective_comp <- enrichplot::pairwise_termsim(protective_comp)

  ### Open Targets Platform

  cli::cli_alert_success("TESTING MULTI_OMICS SIGNATURE ASSOCIATION
                           WITH DISEASE OF INTEREST VIA THE OPENTARGETS API")

  opentarget <- get_diseaseAssociations_df(disease_id = disease_id, size = 3000)
  opentarget_detrimental <- plot_filter_OpenTarget(
    opentarget_results = opentarget, genes = unlist(detrimental)
  )
  opentarget_protective <- plot_filter_OpenTarget(
    opentarget_results = opentarget, genes = unlist(protective)
  )



  #### RESULTS OBJECT
  integrative_results <- list(
    detrimental = list(
      multiomics_signature = detrimental,
      network = detrimental_network,
      network_communities = detrimental_communities,
      cell_type_communities = detrimental_cell_type,
      enrichment = enrichr_detrimental,
      open_targets = opentarget_detrimental
    ),
    protective = list(
      multiomics_signature = protective,
      network = protective_network,
      network_communities = protective_communities,
      cell_type_communities = protective_cell_type,
      enrichment = enrichr_protective,
      open_targets = opentarget_protective
    ),
    comparison = list(
      detrimental = detrimental_comp,
      protective = protective_comp
    )
  )


  return(integrative_results)
}
