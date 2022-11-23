#' Title
#'
#' @param colData
#' @param dependent
#' @param covariates
#' @param levels
#' @param log2FoldChange
#'
#' @return
#' @export
#'
#' @examples
clustering_DE_analysis <- function(normalized_data=multimodal$rna_processed,
                                  colData,
                                  dependent = "Cluster",
                                  levels = c('CS2','CS1'),
                                  log2FoldChange = 0.5 ) {
  library(stats)
  library(lme4)
  library(limma)
  library(dplyr)

  abundance <- normalized_data

  dependent_var <- paste(dependent, collapse = " + ")

  if (!is.null(levels)) {
    colData[[dependent]] <- factor(colData[[dependent]], levels = levels)
  }

  for (i in colnames(colData)) {
    assign(i, colData[, i])
  }

  model_formula <- stats::as.formula(
        sprintf(
          "~0+ %s",
          dependent_var
        )
      )


    design <- model.matrix(model_formula)
    colnames(design) <- c(levels(colData[, dependent]))

    perCombs <- split(
      t(combn(levels(colData[, dependent]), 2)),
      seq(nrow(t(combn(levels(colData[, dependent]), 2))))
    )


    expressions <- list()
    for (i in 1:length(perCombs)) {
      expressions[i] <- paste0(perCombs[[i]][2], "-", perCombs[[i]][1])
    }

    x <- unlist(expressions)

    contr.matrix <- limma::makeContrasts(contrasts = x, levels = colnames(design))

    fit <- limma::lmFit(abundance, design)
    vfit <- limma::contrasts.fit(fit, contrasts = contr.matrix)
    efit <- limma::eBayes(vfit)
    de=limma::topTable(efit, n = Inf)

    list_results <- list()
    for (i in 1:length(perCombs)) {
      list_results[[i]] <- assign(paste0(
        "limma.results_",
        colnames(limma::topTable(efit))[i]
      ), limma::topTable(efit, coef = i, n = Inf))
    }

    names(list_results) <-     expressions


  list_results <- lapply(list_results, function(x) {
    cbind(x, Identifier = rownames(x))
  })

  res <- lapply(list_results, format_res_limma, n_label = 10,
                gene_id_conversion =  NULL,
                log2FoldChange = log2FoldChange  )

  list_DEP_limma <- list()
  list_DEP_limma <- lapply(res, function(x) {
    up=x$gene_name[which(x$adj.P.Val <= 0.05 & x$de=='Up')]
    down=x$gene_name[which(x$adj.P.Val <= 0.05 & x$de=='Down')]
    return(list(up=up,
                down=down))
  })

  res$plot <- lapply(res, volcano_plot_limma,
                     log2FoldChange = log2FoldChange  )
  res$sig_feature <- list_DEP_limma

  names(res) <- gsub(pattern = ".", replacement = "vs", x = names(res), fixed = TRUE)
  names(res$sig_feature) <- gsub(
    pattern = ".", replacement = "vs",
    x = names(res$sig_feature), fixed = TRUE
  )
  names(res$plot) <- gsub(
    pattern = ".", replacement = "vs",
    x = names(res$plot), fixed = TRUE
  )

return(res)
}
