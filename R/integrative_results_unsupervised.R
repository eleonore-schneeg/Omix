integrative_results_unsupervised <- function(multiassay,
                                             integration = "MOFA",
                                             dependent = "diagnosis",
                                             correlation_threshold = 0.5,
                                             disease_id = "MONDO_0004975") {
  model <- multiassay@metadata$integration[[paste(integration)]]
  covariates <- colnames(samples_metadata(model))
  corr <- correlate_factors_with_covariates(model,
    covariates = covariates,
    plot = "r",
    return_data = T
  )
  plot_factor_cor(model)
  plot_factor(model,
    factors = 4,
    color_by = "gpath",
    add_violin = TRUE,
    dodge = TRUE
  )
  # 2 progression check estimates
  # 3 microglia?
  # 4 pathology

  correlate_factors_with_covariates(model, covariates = covariates, plot = "r")
  plot_variance_explained(model, max_r2 = 5)

  selected_factor <- MOFA_get_relevant_factor(MOFAobject = model, covariate = dependent)

  sign <- MOFA_sign_relevant_factor(model,
    covariate = dependent,
    relevant_factor = selected_factor
  )


  weights <- MOFA2::get_weights(model,
    view = "all",
    factor = "all",
    abs = FALSE,
    scale = TRUE,
    as.data.frame = FALSE
  )

  rna <- data.frame(weights = weights[[1]][, selected_factor])
  protein <- data.frame(weights = weights[[2]][, selected_factor])

  top <- 100
  if (sign == 1) {
    ## high pathology in positive factor (negative corr to outcome)
    rna_positive <- rownames(rna %>% top_n(top))
    rna_negative <- rownames(rna %>% top_n(-top))

    protein_positive <- rownames(protein %>% top_n(top))
    protein_negative <- rownames(protein %>% top_n(-top))
  }

  if (sign == -1) {
    ## high pathology in negative factor (positive corr to outcome)
    rna_positive <- rownames(rna %>% top_n(-top))
    rna_negative <- rownames(rna %>% top_n(top))

    protein_positive <- rownames(protein %>% top_n(-top))
    protein_negative <- rownames(protein %>% top_n(top))
  }

  cli::cli_alert_success("RETRIEVAL OF PROTECTIVE AND DETRIMENTAL MULTI-OMICS SIGNATURE")

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
  background <- .get_background(multiassay)

  enrichr_detrimental <- lapply(detrimental, pathway_analysis_enrichr)
  enrichr_protective <- lapply(protective, pathway_analysis_enrichr)

  library(org.Hs.eg.db)
  enrichGO_detrimental <- lapply(
    detrimental,
    function(x) {
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
    }
  )

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


  #### Communities comparisons / mix proteomics transcriptomics if highly corr
  detrimental_comp <- clusterProfiler::compareCluster(
    geneCluster =
      detrimental_communities,
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

  protective_comp <- clusterProfiler::compareCluster(
    geneCluster =
      protective_communities,
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


  ### OPEN TARGETS
  cli::cli_alert_success("TESTING MULTI_OMICS SIGNATURE ASSOCIATION
                           WITH DISEASE OF INTEREST VIA THE OPENTARGETS API")

  opentarget <- get_diseaseAssociations_df(disease_id = disease_id, size = 3000)
  opentarget_detrimental <- plot_filter_OpenTarget(
    opentarget_results =
      opentarget, genes = unlist(detrimental)
  )
  opentarget_protective <- plot_filter_OpenTarget(
    opentarget_results =
      opentarget, genes = unlist(protective)
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


# selects the top factor that is most strongly associated with outcomes of interest
MOFA_get_relevant_factor <- function(MOFAobject, covariate) {
  corr <- correlate_factors_with_covariates(MOFAobject,
    covariates = covariate,
    plot = "r",
    return_data = T
  )


  relevant_factor <- rownames(corr)[which(abs(corr) == max(abs(corr)))]
  relevant_factor_number <- as.numeric(substr(relevant_factor, 7, 9))

  return(relevant_factor_number)
}

MOFA_sign_relevant_factor <- function(MOFAobject, covariate, relevant_factor) {
  corr <- correlate_factors_with_covariates(MOFAobject,
    covariates = covariate,
    plot = "r",
    return_data = T
  )
  sign <- sign(corr[relevant_factor, ])
  return(sign)
}
