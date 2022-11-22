################################################################################
#' Functional enrichment analysis using enrichR
#'
#' Performs impacted pathway analysis with a list of genes.
#'
#' @param interest_genes vector containing genes of interest
#' Column names should be gene, logFC, pval and padj respectively.
#' @param enrichment_database Name of the database for enrichment. User can
#' specify one or more database names from [enrichR::listEnrichrDbs()].
#' @param is_output If TRUE a folder will be created and results of enrichment
#' analysis will be saved otherwise a R list will be returned. Default FALSE
#' @param output_dir Path for the output directory. Default is current dir.
#'
#' @return enrichment_result a list of data.frames containing enrichment output
#' and a list of plots of top 10 significant genesets.
#'
#' @family Impacted pathway analysis
#'
#' @importFrom enrichR enrichr listEnrichrDbs
#' @importFrom cli cli_alert_info cli_text
#' @importFrom ggplot2 ggplot ggsave
#' @importFrom cowplot theme_cowplot background_grid
#' @importFrom stringr str_wrap
#' @importFrom assertthat assert_that
#' @importFrom dplyr %>% mutate
#' @importFrom purrr map map_chr discard
#'
#' @export
pathway_analysis_enrichr <- function(interest_gene = NULL,
                                     project_name = NULL,
                                     enrichment_database = c(
                                       "GO_Molecular_Function_2021",
                                       "GO_Cellular_Component_2021",
                                       "GO_Biological_Process_2021",
                                       "WikiPathways_2021_Human",
                                       "Reactome_2016",
                                       "KEGG_2021_Human",
                                       "MSigDB_Hallmark_2020",
                                       "BioCarta_2016"
                                     ),
                                     is_output = FALSE,
                                     output_dir = ".",
                                     min_overlap= 3,
                                     plot_n = 20) {
  library(enrichR)
  library(stringr)
  library(dplyr)
  library(ggplot2)

  dbs <- enrichR::listEnrichrDbs()

  res <- enrichR::enrichr(
    genes = as.character(interest_gene),
    databases = enrichment_database
  )

  res <- purrr::discard(res, function(x) {
    dim(x)[1] == 0
  })

  enrichr_res <- lapply(names(res), function(x) {
    dt <- res[[x]] %>%
      mutate(database = x)
  })

  names(enrichr_res) <- names(res)

  enrichr_res <- purrr::map(enrichr_res, ~ .format_res_table_enrichr(.))

  enrichr_res <- purrr::discard(enrichr_res, function(x) {
    dim(x)[1] == 0
  })


  if (length(enrichr_res) == 0) {
    cli::cli_text(
      "{.strong No significant impacted pathways found at FDR <= 0.05! }"
    )
    enrichr_res <- NULL
  } else {
    enrichr_res$plot <- lapply(
      enrichr_res,
      function(dt) .dotplot_enrichr(dt, plot_n)
    )

    project_name <- project_name
    output_dir <- output_dir
    sub_dir <- "enrichr_output"
    output_dir_path <- file.path(output_dir, sub_dir)
    project_dir <- file.path(output_dir_path, paste(project_name, sep = ""))

    if (isTRUE(is_output)) {
      dir.create(output_dir_path, showWarnings = FALSE)
      dir.create(project_dir, showWarnings = FALSE)
      lapply(
        names(enrichr_res)[names(enrichr_res) != "plot"],
        function(dt) {
          write.table(enrichr_res[dt],
            file = paste(project_dir, "/", dt, ".tsv", sep = ""),
            row.names = FALSE,
            col.names = gsub(
              dt, "", colnames(enrichr_res[[dt]])
            ), sep = "\t"
          )
        }
      )

      lapply(
        names(enrichr_res$plot),
        function(p) {
          ggplot2::ggsave(paste(project_dir, "/", p, ".png", sep = ""),
            enrichr_res$plot[[p]],
            device = "png", height = 8,
            width = 10, units = "in", dpi = 300
          )
        }
      )
    } else {
      cli::cli_alert_info("Output is returned as a list!")
    }

    enrichr_res$metadata$project_name <- project_name
    enrichr_res$metadata$enrichment_database <- enrichment_database
  }

  return(enrichr_res)
}


#' Format result table
#' @keywords internal

.format_res_table_enrichr <- function(res,
                                      min_overlap=3) {
  res_table <- res %>%
    as.data.frame() %>%
    dplyr::transmute(
      # geneset = res,
      # geneset = .get_geneset(Term),
      description = gsub("\\(GO:.*|Homo sapiens.R-HSA.*|WP.*", "", Term),
      size = as.numeric(gsub(".*\\/", "", Overlap)),
      overlap = as.numeric(gsub("\\/.*", "", Overlap)),
      odds_ratio = round(Odds.Ratio, 2),
      pval = as.numeric(format(P.value, format = "e", digits = 2)),
      FDR = as.numeric(format(Adjusted.P.value, format = "e", digits = 2)),
      Genes = Genes
    )

  res_table<- res_table[which(res_table$overlap >=min_overlap),]
}
# res_table$geneset <- ifelse(is.na(res_table$geneset),
#   res_table$description,
#   res_table$geneset
# )
# res_table <- res_table %>% dplyr::mutate("-Log10(FDR)" = as.numeric(
#   format(-log10(FDR), format = "e", digits = 2)
# ))
#   res_table$genes <- res$Genes
#
#   text_out <- "thyroid|glomerular|renal|retina|nigra|vessel|estrogen|steroid|androgen|artery|bone|skeletal|muscle|aorta|cartilage|pancreatic|myoblast|embryonic|amyotrophic|neural tube|virus|circadian|ectoderm|stem cell|vitamin|chylomicron|coronary|osteoclast|addiction|tumor|myometrial|prolactin|glioblastoma|sensory|cancer|carcinoma|hepatitis|oocyte|cardiomyocyte|heart|cardiac|eye|kidney|viral|ear|infection|auditory|Allograft|lupus|graft|rett|nose"
#
#   res_table <- res_table %>%
#     dplyr::filter(!grepl(text_out, description, ignore.case = T)) %>%
#     filter(size >= 5 & size <= 300) %>%
#     mutate(FDR = p.adjust(pval, method = "BH"))
#
#   res_table <- res_table[res_table$FDR <= 0.1, ]
#   return(res_table)
# }

# .get_geneset <- function(term) {
#   # geneset <- purrr::map_chr(
#   #   as.character(term),
#   #   ~ strsplit(., "(", fixed = TRUE)[[1]][2]
#   # )
#
#   geneset <- purrr::map_chr(as.character(term),
#                             ~ str_extract(. , "GO:.*|R-HSA.*|WP.*"))
#   geneset <- gsub("\\)|Homo sapiens", "", geneset)
#   geneset <- as.character(geneset)
#   return(geneset)
# }

#' dotplot for ORA. x axis perturbation, y axis description
#' @importFrom stats reorder
#' @keywords internal


.dotplot_enrichr <- function(dt, plot_n = 20) {
  dt <- na.omit(dt)
  dt <- top_n(dt, plot_n, -FDR)
  dt$description <- stringr::str_wrap(dt$description, 40)

  ggplot2::ggplot(dt, aes(
    x = odds_ratio,
    y = stats::reorder(description, odds_ratio)
  )) +
    geom_point(aes(fill = FDR, size = size),
      shape = 21, alpha = 0.7, color = "black"
    ) +
    scale_size(name = "size", range = c(3, 8)) +
    xlab("Total Odds Ratio") +
    ylab("") +
    scale_fill_gradient(
      low = "navy", high = "gold", name = "FDR",
      guide = guide_colorbar(reverse = TRUE),
      limits = c(0, 0.1),
      aesthetics = c("fill")
    ) +
    guides(size = guide_legend(
      override.aes = list(fill = "gold", color = "gold")
    )) +
    cowplot::theme_cowplot() +
    cowplot::background_grid()
}

dot_plot_enrichr_semantics <- function(dt,
                                       semantics = "microglial") {
  semantics_rows <- dt$description[grepl(paste0(semantics), dt$description, fixed = TRUE)]
  dt <- na.omit(dt)
  dt <- dt[which(dt$description %in% semantics_rows), ]
  dt$description <- stringr::str_wrap(dt$description, 40)

  plot <- ggplot2::ggplot(dt, aes(
    x = odds_ratio,
    y = stats::reorder(description, odds_ratio)
  )) +
    geom_point(aes(fill = FDR, size = size),
      shape = 21, alpha = 0.7, color = "black"
    ) +
    scale_size(name = "size", range = c(3, 8)) +
    xlab("Total Odds Ratio") +
    ylab("") +
    scale_fill_gradient(
      low = "navy", high = "gold", name = "FDR",
      guide = guide_colorbar(reverse = TRUE),
      limits = c(0, 0.1),
      aesthetics = c("fill")
    ) +
    guides(size = guide_legend(
      override.aes = list(fill = "gold", color = "gold")
    )) +
    cowplot::theme_cowplot() +
    cowplot::background_grid()

  return(plot)
}



dotplot_enrichment_MOFA <- function(enrichment.results, factor, alpha = 0.1, max.pathways = 25,
                                    text_size = 1.0, dot_size = 5.0) {
  # Sanity checks
  stopifnot(is.numeric(alpha))
  stopifnot(length(factor) == 1)
  if (is.numeric(factor)) factor <- colnames(enrichment.results$pval.adj)[factor]
  if (!factor %in% colnames(enrichment.results$pval)) {
    stop(paste0("No gene set enrichment calculated for factor ", factor))
  }

  # get p-values
  p.values <- enrichment.results$pval.adj

  # Get data
  tmp <- data.frame(
    pvalues = p.values[, factor, drop = TRUE],
    pathway = rownames(p.values)
  )

  # Filter out pathways
  tmp <- tmp[tmp$pvalue <= alpha, , drop = FALSE]
  if (nrow(tmp) == 0) stop("No siginificant pathways at the specified alpha threshold")

  # If there are too many pathways enriched, just keep the 'max_pathways' more significant
  if (nrow(tmp) > max.pathways) tmp <- head(tmp[order(tmp$pvalue), ], n = max.pathways)

  # Convert pvalues to log scale
  tmp$logp <- -log10(tmp$pvalue + 1e-100)

  # order according to significance
  tmp$pathway <- factor(tmp$pathway <- rownames(tmp), levels = tmp$pathway[order(tmp$pvalue, decreasing = TRUE)])
  tmp$start <- 0

  ggplot2::ggplot(tmp, aes_string(x = "logp", y = "pathway")) +
    geom_point(aes(fill = logp),
      size = 5,
      shape = 21, alpha = 0.7, color = "black"
    ) +
    scale_size(name = "size", range = c(3, 8)) +
    xlab("-log10(p value)") +
    ylab("") +
    scale_fill_gradient(
      low = "navy", high = "gold", name = "-log10(p value)",
      guide = guide_colorbar(reverse = TRUE),
      limits = c(0, max(tmp$logp)),
      aesthetics = c("fill")
    ) +
    cowplot::theme_cowplot() +
    cowplot::background_grid()
}




dotplot_enrichment_MOFA_semantics <- function(enrichment.results, semantics, factor, alpha = 0.05, max.pathways = 25,
                                              text_size = 1.0, dot_size = 5.0) {
  factor <- seq(1, dim(enrichment.results$pval.adj)[2], by = 1)
  # get p-values
  p.values <- enrichment.results$pval.adj

  # Get data
  tmp <- data.frame(
    pvalues = p.values[, factor, drop = TRUE],
    pathway = rownames(p.values)
  )

  # order according to significance
  semantics_rows <- rownames(tmp)[grepl(paste0(semantics), rownames(tmp), fixed = TRUE)]
  tmp <- tmp[which(rownames(tmp) %in% semantics_rows), ]

  data_long <- gather(tmp, factor, pvalue, pvalues.Factor1:pvalues.Factor11, factor_key = TRUE)
  data_long$factor <- substr(data_long$factor, 9, 20)
  paths <- unique(data_long$pathway)
  data_long$logpvalue <- -log10(data_long$pvalue)

  semantics_plots <- list()

  for (i in paths) {
    print(i)
    data <- data_long[which(data_long$pathway == i), ]
    plot <- ggplot2::ggplot(data, aes_string(x = "logpvalue", y = "factor")) +
      geom_point(aes(fill = logpvalue),
        size = 5,
        shape = 21, alpha = 0.7, color = "black"
      ) +
      scale_size(name = "size", range = c(3, 8)) +
      xlab("-log10(p value)") +
      xlim(c(0, 1)) +
      ylab("") +
      ggtitle(paste0(i)) +
      scale_fill_gradient(
        low = "navy", high = "gold", name = "-log10(p value)",
        guide = guide_colorbar(reverse = TRUE),
        limits = c(0, max(data$logpvalue)),
        aesthetics = c("fill")
      ) +
      cowplot::theme_cowplot() +
      cowplot::background_grid()

    semantics_plots[[i]] <- plot
  }
  return(semantics_plots)
}


enrichment_custom <- function(genes, reference, genesets, adj = "fdr", verbose = FALSE) {
  tab <- lapply(1:length(genesets), function(i) {
    if (verbose == TRUE) {
      cat("processing term", i, names(genesets)[i], "\n")
    }

    ra <- length(genes) / length(reference)
    geneset_size <- length(genesets[[i]])
    expect <- geneset_size * ra

    reference <- reference[!reference %in% genes]
    RinSet <- sum(reference %in% genesets[[i]])
    RninSet <- length(reference) - RinSet
    GinSet <- sum(genes %in% genesets[[i]])
    enrichment_ratio <- round(GinSet / expect, 2)
    overlap_genes <- paste(intersect(genes, genesets[[i]]), collapse = ",")

    GninSet <- length(genes) - GinSet
    fmat <- matrix(c(GinSet, RinSet, GninSet, RninSet),
      nrow = 2,
      ncol = 2, byrow = F
    )
    colnames(fmat) <- c("inSet", "ninSet")
    rownames(fmat) <- c("genes", "reference")
    fish <- fisher.test(fmat, alternative = "greater")
    pval <- as.numeric(format(fish$p.value, format = "e", digits = 2))
    inSet <- RinSet + GinSet
    pct_overlap <- round((GinSet / inSet) * 100, 2)
    res <- c(GinSet, inSet, geneset_size, pct_overlap, enrichment_ratio, pval, overlap_genes)
    res
  })
  rtab <- do.call("rbind", tab)
  rtab <- data.frame(as.vector(names(genesets)), rtab)
  rtab <- rtab[order(rtab[, 4]), ]
  colnames(rtab) <- c(
    "TermID",
    "genes", "all", "geneset_size",
    "pct_overlap", "enrichment_ratio",
    "pval", "overlap_genes"
  )
  tab.out <- rtab %>%
    dplyr::mutate_at(c(
      "genes", "all", "geneset_size",
      "pct_overlap", "enrichment_ratio",
      "pval"
    ), as.character) %>%
    dplyr::mutate_at(c(
      "genes", "all", "geneset_size",
      "pct_overlap", "enrichment_ratio",
      "pval"
    ), as.numeric) %>%
    dplyr::mutate(padj = p.adjust(pval), method = adj) %>%
    dplyr::mutate(padj = as.numeric(format(padj, format = "e", digits = 2)))
  tab.out <- tab.out %>%
    dplyr::select(-overlap_genes, everything(), overlap_genes)

  return(tab.out)
}


cell_type_enrichment <- function(multiassay,
                                 communities) {

  background_genes <- .get_background(multiassay)
  load(file = "~/RDS_ukdri/multiomics/multiomics/pipeline_17_05/3_UNI_OMIC/CellTypeData_all_ds.rda") # loads ctd

  hits <- communities
  hits <- lapply(hits, function(x) {
    x[x %in% background_genes]
  })
  hits <- lapply(hits, function(x) x[length(x) >= 4])
  hits <- hits[lapply(hits, length) > 0]

  reps <- 100
  annotLevel <- 1


  sct <- ctd[2] # want to get cell type annotation 2

  results <- list()
  results <- lapply(hits, function(x) {
    EWCE::bootstrap_enrichment_test(
      sct_data = sct,
      sctSpecies = "human",
      genelistSpecies = "human",
      hits = x,
      reps = reps,
      bg = background_genes
    )
  })


  resultsPlots <- lapply(results, function(x) x[["results"]])

  plots <- list()
  plots <- lapply(resultsPlots, function(x) {
    EWCE::ewce_plot(
      total_res = x,
      mtc_method = "bonferroni",
      ctd = ctd
    )
  })

  plot_cell_enrichment <- list()
  plot_cell_enrichment <- lapply(plots, function(x) x[["plain"]])
  plot_cell_enrichment <- lapply(names(plot_cell_enrichment), function(x) {
    plot_cell_enrichment[[x]] + ggplot2::ggtitle(paste("Cell type enrichment in community #", x))
  })

  names(resultsPlots) <- names(hits)
  names(plot_cell_enrichment) <- names(hits)
  cell_type_enrichment <- list(results = resultsPlots, plots = plot_cell_enrichment)

  return(cell_type_enrichment)
}
