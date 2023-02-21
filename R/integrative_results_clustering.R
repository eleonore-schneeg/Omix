#' Title
#'
#' @param multiassay
#' @param intersect_genes
#' @param ID_type
#' @param integration
#' @param
#' @param annotations_cat
#' @param res.path
#' @param correlation_threshold
#'
#' @return
#' @export
#'
#' @examples
integrative_results_clustering <- function(multiassay,
                                           intersect_genes = FALSE,
                                           slots,
                                           ID_type,
                                           integration,
                                           enrichment_method,
                                           annotations_continuous,
                                           annotations_cat,
                                           res.path,
                                           correlation_threshold = 0.5,
                                           geneset = NULL) {
  multimodal_object <- .get_multimodal_object(
    multiassay = multiassay,
    slots = slots,
    intersect_genes = intersect_genes,
    ID_type = ID_type
  )

  multimodal <- multimodal_object[[1]]
  metadata <- multimodal_object[[2]]
  multimodal <- lapply(multimodal, data.frame)

  cluster <-multiassay@metadata$integration[[paste(integration)]]

  plot_data <- MOVICS::getStdiz(
    data = multimodal,
    halfwidth = c(2, 2),
    centerFlag = c(T, T),
    scaleFlag = c(T, T)
  )


  columns <- c(annotations_continuous, annotations_cat)
  annCol <- metadata[, c(paste(columns)), drop = FALSE]

  # annColors <- list(
  #   amyloid = circlize::colorRamp2(
  #     breaks = c(
  #       min(annCol$amyloid),
  #       median(annCol$amyloid),
  #       max(annCol$amyloid)
  #     ),
  #     colors = c("#0000AA", "#555555", "#AAAA00")
  #   ),
  #   nft = circlize::colorRamp2(
  #     breaks = c(
  #       min(annCol$nft),
  #       median(annCol$nft),
  #       max(annCol$nft)
  #     ),
  #     colors = c("#0000AA", "#555555", "#AAAA00")
  #   ),
  #   gpath = circlize::colorRamp2(
  #     breaks = c(
  #       min(annCol$gpath),
  #       median(annCol$gpath),
  #       max(annCol$gpath)
  #     ),
  #     colors = c("#0000AA", "#555555", "#AAAA00")
  #   ),
  #   progression = c(
  #     "fast_progressing" = "green",
  #     "slow_progressing" = "blue"
  #   ),
  #   cog_diag = c(
  #     "CO" = "green",
  #     "MCI" = "blue",
  #     "AD" = "red"
  #   ),
  #   braaksc = c("3" = "blue", "4" = "red")
  # )

  cli::cli_alert_success("MULTI-OMICS CLUSTERING HEATMAP")
  # heatmap <- MOVICS::getMoHeatmap(
  #   data = plot_data,
  #   row.title = c("mRNA", "Proteins"),
  #   is.binary = c(F, F), # the 4th data is mutation which is binary
  #   legend.name = c("mRNA.counts", "Protein.abundance"),
  #   clust.res = cluster$clust.res, # consensusMOIC results
  #   clust.dend = NULL, # show no dendrogram for samples
  #   show.rownames = c(F, F), # specify for each omics data
  #   show.colnames = FALSE, # show no sample names
  #   show.row.dend = c(F, F),
  #   annCol = annCol, # show no dendrogram for features
  #   annRow = annColors, # no selected features
  #   width = 10, # width of each subheatmap
  #   height = 5, # height of each subheatmap
  #   fig.name = "Multiomics clustering heatmap",
  #   fig.path = res.path
  # )

  annCol <- mutate_if(annCol, is.character, as.factor)

  surv.info <- annCol
  surv.info$Subtype <- cluster$clust.res$clust

  cli::cli_alert_success("RELATING CLUSTERS TO CLINICAL INFORMATION")
  clinical.clust <- MOVICS::compClinvar(
    moic.res = cluster,
    var2comp = surv.info, # data.frame needs to summarize (must has row names of samples)
    strata = "Subtype",
    factorVars = c("var1", "var2", "var3", "var4"),
    doWord = TRUE,
    includeNA = F, # generate .docx file in local path
    tab.name = "Clinical features per clusters"
  )

  surv.info$Cluster <- paste0("CS", surv.info$Subtype)

  cli::cli_alert_success("RUNNING DIFFERENTIAL EXPRESSION ANALYSIS")
  res1 <- clustering_DE_analysis(
    normalized_data = multimodal$rna_processed,
    colData = surv.info,
    dependent = "Cluster",
    levels = c("CS2", "CS1"),
    log2FoldChange = 0.5
  )


  res2 <- clustering_DE_analysis(
    normalized_data = multimodal$protein_processed,
    colData = surv.info,
    dependent = "Cluster",
    levels = c("CS2", "CS1"),
    log2FoldChange = 0.2
  )

  DE_res <- list(
    rna = res1,
    proteins = res2
  )

  Up <- list()
  Down <- list()

  for (i in names(DE_res$rna$sig_feature)) {
    Up[[i]] <- list(
      rna = DE_res$rna$sig_feature[[i]]$up,
      protein = DE_res$protein$sig_feature[[i]]$up
    )
    Down[[i]] <- list(
      rna = DE_res$rna$sig_feature[[i]]$down,
      protein = DE_res$protein$sig_feature[[i]]$down
    )
  }

  cli::cli_alert_success("BIPARTITE NETWORKS GENERATION OF DIFFERENTIALLY EXPRESSED FEATURES IN EACH CLUSTER")

  Down_network <- list()
  Up_network <- list()
  Up_communities <- list()
  Down_communities <- list()
  Up_cell_type <- list()
  Down_cell_type <- list()

  for (i in names(Up)) {
    cli::cli_alert_success(paste("UP REGULATED FEATURES (", i, ") DOWNSTREAM ANALYSIS"))

    cluster <- str_split(names(Up), "-", n = 2, simplify = FALSE)
    cluster <- readr::parse_number(cluster[[1]][1])

    Up_network[[i]] <- .multiomics_network_cluster(
      multiassay = multiassay,
      integration = "iCluster",
      cluster = cluster,
      list = Up[[i]],
      correlation_threshold = correlation_threshold
    )
    cli::cli_alert_success("DETECTION OF NETWORK COMMUNITIES")
    Up_communities[[i]] <- .communities_network(Up_network[[i]])

    cli::cli_alert_success("CELL TYPE ENRICHMENT OF COMMUNITIES")
    Up_cell_type[[i]] <- cell_type_enrichment(
      multiassay = multiassay,
      communities = Up_communities[[i]]
    )
  }

  for (i in names(Down)) {
    cli::cli_alert_success(paste("DOWN REGULATED FEATURES (", i, ") DOWNSTREAM ANALYSIS"))

    cluster <- str_split(names(Down), "-", n = 2, simplify = FALSE)
    cluster <- readr::parse_number(cluster[[1]][1])
    Down_network[[i]] <- .multiomics_network_cluster(
      multiassay = multiassay,
      integration = "iCluster",
      cluster = cluster,
      list = Down[[i]],
      correlation_threshold = correlation_threshold
    )



    cli::cli_alert_success("DETECTION OF NETWORK COMMUNITIES")
    Down_communities[[i]] <- .communities_network(Down_network[[i]])

    cli::cli_alert_success("CELL TYPE ENRICHMENT OF COMMUNITIES")

    Down_cell_type[[i]] <- cell_type_enrichment(
      multiassay = multiassay,
      communities = Down_communities[[i]]
    )
  }


  cli::cli_alert_success("FUNCTIONAL ENRICHMENT OF UP AND DOWN REGULATED FEATURES")
  for (i in names(Up)) {
    Up[[i]] <- lapply(Up[[i]], function(x) {
      sub("*\\.[0-9]", "", x)
    })
  }

  for (i in names(Down)) {
    Down[[i]] <- lapply(Down[[i]], function(x) {
      sub("*\\.[0-9]", "", x)
    })
  }

  background <- .get_background(multiassay = multiassay, of = "full")


  if (enrichment_method == "enrichr") {
    enrichr_Down <- list()
    enrichr_Up <- list()
    for (i in names(Up)) {
      enrichr_Down[[i]] <- lapply(Down[[i]], pathway_analysis_enrichr)
      enrichr_Up[[i]] <- lapply(Up[[i]], pathway_analysis_enrichr)
    }
  }

  if (enrichment_method == "enrichKEGG") {
    df <- DE_res[["rna"]][["CS1-CS2"]]
    df <- df[which(df$de != "Not sig"), ]
    entrez <- mapIds(org.Hs.eg.db, df$gene_name, "ENTREZID", "SYMBOL")
    kk <- clusterProfiler::enrichKEGG(
      gene = entrez,
      organism = "hsa",
      pvalueCutoff = 0.05
    )

    top_kk <- top_n(kk@result, 1, Count)
    top_kk_id <- top_kk$ID
    list <- df$log2FoldChange
    names(list) <- mapIds(org.Hs.eg.db, df$gene_name, "ENTREZID", "SYMBOL")

    top_kk_id <- pathview(
      gene.data = list,
      pathway.id = paste(top_kk_id),
      species = "hsa",
      limit = list(gene = round(max(abs(list)), 1))
    )
  }

  if (!is.null(custom_geneset)) {
    custom_Down <- list()
    custom_Up <- list()

    for (i in names(Up)) {
      custom_Down[[i]] <- lapply(Down[[i]], function(x) {
        enrichment_custom(x, background, geneset, adj = "fdr", verbose = FALSE)
      })

      custom_Up[[i]] <- lapply(Down[[i]], function(x) {
        enrichment_custom(x, background, geneset, adj = "fdr", verbose = FALSE)
      })
    }
  }

  # Comparison of communities
  cli::cli_alert_success("FUNCTIONALLY COMPARING UP AND DOWN REGULATED FEATURES")

  comp_list <- list()
  comp <- list()

  for (i in names(Up)) {
    for (j in names(Up[[i]])) {
      comp_list[[i]][[j]] <- list(
        up = Up[[i]][[j]],
        down = Down[[i]][[j]]
      )

      cli::cli_alert_success(paste("Comparing up vs down regulated", j, "in", i))

      comp[[i]][[j]] <- clusterProfiler::compareCluster(
        geneCluster = comp_list[[i]][[j]],
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
    }
  }

  #### RESULTS OBJECT
  integrative_results <- list(
    clinical = clinical.clust,
    DE = DE_res,
    Down = list(
      multiomics_signature = Down,
      network = Down_network,
      network_communities = Down_communities,
      cell_type_communities = Down_cell_type,
      enrichment = enrichr_Down
    ),
    Up = list(
      multiomics_signature = Up,
      network = Up_network,
      network_communities = Up_communities,
      cell_type_communities = Up_cell_type,
      enrichment = enrichr_Up
    ),
    comparison = comp
  )


  return(integrative_results)
}
