################################################################################
#' Performs single omic analysis of transcriptomics layer.
#'
#' @param multiassay Multiassay object genereted by Omix
#' @param dependent Dependent variable in the differential expression analysis
#' @param levels Levels for dependent variable. Control group should be put
#' first.
#' @param fiter_protein_coding Logical whether to filter out non protein coding
#' genes
#'
#' @return
#' @importFrom dplyr filter
#' @importFrom DESeq2 DESeq results
#' @export
#'
#' @examples
rna_single_analysis <- function(multiassay,
                                dependent = "diagnosis",
                                levels = NULL,
                                filter_protein_coding=TRUE) {
  dds <- multiassay@metadata[["dds"]]
  dds <- DESeq2::DESeq(dds)
  res <- list()

  if (!is.numeric(dds[[dependent]])) {
    if (!is.null(levels)) {
      dds[[dependent]] <- factor(dds[[dependent]], levels = levels)
    }
    perCombs <- split(
      t(combn(levels(dds[[dependent]]), 2)),
      seq(nrow(t(combn(levels(dds[[dependent]]), 2))))
    )

    names <- list()
    for (i in 1:length(perCombs)) {
      names[[i]] <- c(perCombs[[i]][2], perCombs[[i]][1])
    }

    list_contrast <- lapply(names, function(x) {
      c(dependent, x)
    })
    res <- lapply(list_contrast, function(x) {
      DESeq2::results(dds,
        lfcThreshold = 0,
        alpha = 0.05,
        contrast = x,
        independentFiltering = TRUE
      )
    })
  }

  names2 <- list()
  for (i in 1:length(perCombs)) {
    names2[i] <- paste0(perCombs[[i]][2], "vs", perCombs[[i]][1])
  }

  names(res) <- unlist(names2)

  gene_id_conversion <- as.data.frame(multiassay@ExperimentList@listData[["rna_raw"]]@elementMetadata@listData)

  res <- lapply(res, format_res_deseq,
    n_label = 20,
    gene_id_conversion = gene_id_conversion
  )

  list_DEG_limma <- list()
  list_DEG_limma <- lapply(res, function(x) {
    x$gene_name[which(x$padj <= 0.05 &
      x$gene_type == "protein_coding")]
  })
  names(list_DEG_limma) <- unlist(names2)

  res$plot <- lapply(res, function(x) {
    volcano_plot_deseq(dt = x)
  })

  res$sig_gene <- list_DEG_limma
  metadata(multiassay)$DEG <- res

  return(multiassay)
}
