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
#' @family Enrichment analysis
#'
#' @importFrom enrichR enrichr
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
                                     min_overlap = 3,
                                     plot_n = 20) {
  eval(parse(text = "enrichR:::.onAttach()")) # R CMD check workaround

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


#' Format result table for `pathway_analysis_enrichr()`
#'
#' @param res Internal to `pathway_analysis_enrichr()`
#' @param min_overlap Default to 3
#' @importFrom dplyr transmute
#'
#' @keywords internal

.format_res_table_enrichr <- function(res,
                                      min_overlap = 3) {
  res_table <- res %>%
    as.data.frame() %>%
    dplyr::transmute(
      description = gsub("\\(GO:.*|Homo sapiens.R-HSA.*|WP.*", "", Term),
      size = as.numeric(gsub(".*\\/", "", Overlap)),
      overlap = as.numeric(gsub("\\/.*", "", Overlap)),
      odds_ratio = round(Odds.Ratio, 2),
      pval = as.numeric(format(P.value, format = "e", digits = 2)),
      FDR = as.numeric(format(Adjusted.P.value, format = "e", digits = 2)),
      Genes = Genes
    )

  res_table <- res_table[which(res_table$overlap >= min_overlap), ]
}

#' dotplot for ORA. x axis perturbation, y axis description
#' @importFrom stats reorder
#' @importFrom dplyr top_n
#' @importFrom stringr str_wrap
#' @importFrom cowplot theme_cowplot background_grid
#' @import ggplot2
#' @keywords internal


.dotplot_enrichr <- function(dt, plot_n = 20) {
  dt <- na.omit(dt)
  dt <- dplyr::top_n(dt, plot_n, -FDR)
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

#' Curates dotplot for `pathway_analysis_enrichr()` on a required semantics
#'
#' @param dt functional enrichment results
#' @param semantics vector of biological terms
#'
#' @return functional enrichement plot
#' @keywords internal

dot_plot_enrichr_semantics <- function(dt,
                                        semantics = "microglial") {
  semantics_rows <- dt$description[grepl(paste0(semantics),
                                         dt$description, fixed = TRUE)]
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


#' Functional enrichment based on given reference annotations
#'
#' @param genes Gene list (vector)
#' @param reference Reference background (vector)
#' @param genesets List of genesets
#' @param adj Multiple testing adjustment method. Default to 'fdr'
#' @param verbose Default to FALSE
#'
#' @return Functional enrichment result
#'
#' @family Enrichment analysis
#'
#' @importFrom dplyr mutate_at mutate select everything
#' @export

enrichment_custom <- function(genes,
                              reference,
                              genesets,
                              adj = "fdr",
                              verbose = FALSE) {
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
      ncol = 2,
      byrow = FALSE
    )
    colnames(fmat) <- c("inSet", "ninSet")
    rownames(fmat) <- c("genes", "reference")
    fish <- fisher.test(fmat, alternative = "greater")
    pval <- as.numeric(format(fish$p.value, format = "e", digits = 2))
    inSet <- RinSet + GinSet
    pct_overlap <- round((GinSet / inSet) * 100, 2)
    res <- c(GinSet, inSet, geneset_size, pct_overlap,
             enrichment_ratio, pval, overlap_genes)
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
    dplyr::select(-overlap_genes, dplyr::everything(), overlap_genes)

  return(tab.out)
}


#' Cell type enrichment of communities
#'
#' @param multiassay  Multiassay experiment object generated by Omix
#' @param communities communities generated from `communities_network()`
#'
#' @return List of results
#'
#' @family Enrichment analysis
#'
#' @importFrom EWCE bootstrap_enrichment_test
#' @importFrom ggplot2 ggtitle
#' @export
#'

cell_type_enrichment <- function(multiassay,
                                 communities) {
  background_genes <- get_background(multiassay)

  hits <- communities
  hits <- lapply(hits, function(x) {
    x[x %in% background_genes]
  })
  hits <- lapply(hits, function(x) x[length(x) >= 4])
  hits <- hits[lapply(hits, length) > 0]

  reps <- 100
  annotLevel <- 1

  sct <- ctd
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
    plot_cell_enrichment[[x]] +
      ggplot2::ggtitle(paste("Cell type enrichment in community #", x))
  })

  names(resultsPlots) <- names(hits)
  names(plot_cell_enrichment) <- names(hits)
  cell_type_enrichment <- list(results = resultsPlots,
                               plots = plot_cell_enrichment)

  return(cell_type_enrichment)
}
