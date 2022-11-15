#' Title
#'
#' @param multiassay
#' @param slot
#' @param dependent
#' @param covariates
#' @param reference
#'
#' @return
#' @export
#'
#' @examples
protein_single_analysis <- function(multiassay,
                                    slot = "protein_processed",
                                    dependent = "diagnosis",
                                    covariates = c("pmi", "msex", "age_death"),
                                    levels = c("Control", "Mid", "Late")) {
  library(MultiAssayExperiment)
  library(stats)
  library(lme4)
  library(limma)
  library(dplyr)

  abundance <- multiassay@ExperimentList@listData[[paste(slot)]]

  protein_raw <- MultiAssayExperiment::getWithColData(multiassay, i = "protein_raw", mode = "replace", verbose = F)
  colData <- data.frame(SummarizedExperiment::colData(protein_raw))
  order <- map_protein$primary[match(colnames(abundance), map_protein$colname)]
  colData <- colData[order, ]
  colData <- colData %>% mutate(across(where(is.character), as.factor))


  cov_var <- paste(covariates, collapse = " + ")
  dependent_var <- paste(dependent, collapse = " + ")



  if (!is.null(levels)) {
    colData[[dependent]] <- factor(colData[[dependent]], levels = levels)
  }

  for (i in colnames(colData)) {
    assign(i, colData[, i])
  }


  if (!is.numeric(colData[, dependent])) {
    model_formula <- stats::as.formula(
      sprintf(
        "~0+ %s + %s",
        dependent_var,
        cov_var
      )
    )

    design <- model.matrix(model_formula)

    colnames(design) <- c(levels(colData[, dependent]), covariates)

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

    list_results <- list()
    for (i in 1:length(perCombs)) {
      list_results[[i]] <- assign(paste0(
        "limma.results_",
        colnames(limma::topTable(efit))[i]
      ), limma::topTable(efit, coef = i, n = Inf))
    }

    names(list_results) <- colnames(topTable(efit))[1:length(perCombs)]
  }

  if (is.numeric(colData[, dependent])) {
    model_formula <- stats::as.formula(
      sprintf(
        "~%s + %s",
        dependent_var,
        cov_var
      )
    )

    design <- model.matrix(model_formula)
    fit <- lmFit(abundance, design)
    efit <- eBayes(fit)

    list_results <- list()
    for (i in 1:length(1)) {
      list_results[[i]] <- assign(
        paste0(
          "limma.results_",
          colnames(limma::topTable(efit))[i]
        ),
        limma::topTable(efit, coef = paste(dependent), adjust = "BH", number = Inf)
      )
    }

    names(list_results) <- paste(dependent)
  }

  list_results <- lapply(list_results, function(x) {
    cbind(x, Identifier = rownames(x))
  })
  gene_id_conversion <- as.data.frame(multiassay@ExperimentList@listData[["protein_raw"]]@elementMetadata@listData)

  res <- lapply(list_results, format_res_limma, n_label = 10, gene_id_conversion = gene_id_conversion)

  list_DEP_limma <- list()
  list_DEP_limma <- lapply(res, function(x) {
    x$Identifier[which(x$adj.P.Val <= 0.05)]
  })

  res$plot <- lapply(res, volcano_plot_limma)
  res$sig_protein <- list_DEP_limma

  names(res) <- gsub(pattern = ".", replacement = "vs", x = names(res), fixed = TRUE)
  names(res$sig_protein) <- gsub(
    pattern = ".", replacement = "vs",
    x = names(res$sig_protein), fixed = TRUE
  )
  names(res$plot) <- gsub(
    pattern = ".", replacement = "vs",
    x = names(res$plot), fixed = TRUE
  )

  metadata(multiassay)$DEP <- res

  return(multiassay)
}
