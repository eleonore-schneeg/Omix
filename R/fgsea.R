#' Performs functional gene set enrichment analysis
#'
#' @param gene_list vector of genes
#' @param pval threshold for adjusted pvalue
#'
#' @return List of results
#'
#' @family Enrichment analysis
#'
#' @importFrom fgsea fgsea collapsePathways
#' @importFrom dplyr filter arrange desc %>%
#' @importFrom msigdbr msigdbr
#' @importFrom data.table as.data.table
#' @export
#'

GSEA <- function(gene_list, pval) {
  set.seed(54321)

  if (any(duplicated(names(gene_list)))) {
    warning("Duplicates in gene names")
    gene_list <- gene_list[!duplicated(names(gene_list))]
  }
  if (!all(order(gene_list, decreasing = TRUE) == 1:length(gene_list))) {
    warning("Gene list not sorted")
    gene_list <- sort(gene_list, decreasing = TRUE)
  }

  cgp_gene_sets <- msigdbr::msigdbr(species = "human",
                                    category = "C5",
                                    subcategory = "GO:BP")
  msigdbr_list <- split(x = cgp_gene_sets$gene_symbol, f = cgp_gene_sets$gs_name)
  names(msigdbr_list) <- sub(".*GOBP_", "", names(msigdbr_list))
  names(msigdbr_list) <- gsub("_", " ", names(msigdbr_list))
  names(msigdbr_list) <- tolower(names(msigdbr_list))
  myGO <- msigdbr_list

  fgRes <- suppressWarnings({
    fgsea::fgsea(
      pathways = myGO,
      stats = gene_list,
      minSize = 15, ## minimum gene set size
      maxSize = 400, ## maximum gene set size
      nperm = 10000
    )
  }) %>%
    as.data.frame() %>%
    dplyr::filter(padj < !!pval) %>%
    dplyr::arrange(dplyr::desc(NES))
  message(paste("Number of signficant gene sets =", nrow(fgRes)))

  message("Collapsing Pathways -----")
  concise_pathways <- fgsea::collapsePathways(data.table::as.data.table(fgRes),
    pathways = myGO,
    stats = gene_list
  )
  fgRes <- fgRes[fgRes$pathway %in% concise_pathways$mainPathways, ]
  message(paste("Number of gene sets after collapsing =", nrow(fgRes)))

  fgRes$Enrichment <- ifelse(fgRes$NES > 0, "Up-regulated", "Down-regulated")
  filtRes <- rbind(
    head(fgRes, n = 10),
    tail(fgRes, n = 10)
  )

  total_up <- sum(fgRes$Enrichment == "Up-regulated")
  total_down <- sum(fgRes$Enrichment == "Down-regulated")
  header <- paste0("Top 10 (Total pathways: Up=", total_up,
                   ", Down=", total_down, ")")

  colos <- setNames(
    c("firebrick2", "dodgerblue2"),
    c("Up-regulated", "Down-regulated")
  )

  g1 <- ggplot(filtRes, aes(reorder(pathway, NES), NES)) +
    geom_point(aes(fill = Enrichment, size = size), shape = 21) +
    scale_fill_manual(values = colos) +
    scale_size_continuous(range = c(2, 10)) +
    geom_hline(yintercept = 0) +
    coord_flip() +
    labs(
      x = "Pathway", y = "Normalized Enrichment Score",
      title = header
    )

  output <- list("Results" = fgRes, "Plot" = g1)
  return(output)
}
