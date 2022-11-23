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
protein_DE_analysis <- function(multiassay,
                                    slot = "protein_processed",
                                    dependent = "diagnosis",
                                    covariates = NA,
                                    levels = c("Control", "Mid", "Late"),
                                    log2FoldChange = 0.5 ) {
  library(MultiAssayExperiment)
  library(stats)
  library(lme4)
  library(limma)
  library(dplyr)

  abundance <- multiassay@ExperimentList@listData[[paste(slot)]]


  if (("parameters_processing_protein" %in% names(multiassay@metadata)) == FALSE) {
    stop(cli::cli_alert_danger(
      paste("parameters_processing_protein not found in metadata, please run",
            cli::style_bold("process_protein"),
            sep = " "
      )
    ))
  }
  denoised <-multiassay@metadata$parameters_processing_protein$denoise
  cov <-multiassay@metadata$parameters_processing_protein$covariates


  if (isTRUE(denoised & isTRUE(cov == covariates))) {
    stop(cli::cli_alert_danger(
      paste("Processed protein data is already for designated covariates, model will be run for unadjusted covariates only")
    ))
  }

  covariates <- setdiff(covariates,cov)
  protein_raw <- MultiAssayExperiment::getWithColData(multiassay, i = "protein_raw", mode = "replace", verbose = F)
  colData <- data.frame(SummarizedExperiment::colData(protein_raw))
  order <- map_protein$primary[match(colnames(abundance), map_protein$colname)]
  colData <- colData[order, ]
  colData <- colData %>% mutate(across(where(is.character), as.factor))

  if(!is.na(covariates)){
  cov_var <- paste(covariates, collapse = " + ")

  if (any(is.na(colData[cov_var]))) {
    stop(cli::cli_alert_danger(
      paste("Missing value in covariates, please use other covariates")
    ))
  }

  }else{
    cov_var <- NA
  }

  dependent_var <- paste(dependent, collapse = " + ")

  if (!is.null(levels)) {
    colData[[dependent]] <- factor(colData[[dependent]], levels = levels)
  }

  for (i in colnames(colData)) {
    assign(i, colData[, i])
  }


  if (!is.numeric(colData[, dependent])) {
    if(!is.na(cov_var)){
    model_formula <- stats::as.formula(
      sprintf(
        "~0+ %s + %s",
        dependent_var,
        cov_var
      )
    )
   }
    if(is.na(cov_var)){
        model_formula <- stats::as.formula(
          sprintf(
            "~0+ %s",
            dependent_var
          )
        )
    }

    design <- model.matrix(model_formula)


    if(!is.na(cov_var)){
    colnames(design) <- c(levels(colData[, dependent]), cov_var)
    }else{
      colnames(design) <- c(levels(colData[, dependent]))
    }

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
    if(!is.na(cov_var)){
    model_formula <- stats::as.formula(
      sprintf(
        "~%s + %s",
        dependent_var,
        cov_var
      )
    )
    }

    if(is.na(cov_var)){
        model_formula <- stats::as.formula(
          sprintf(
            "~0+ %s",
            dependent_var
          )
        )
      }


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

  res <- lapply(list_results, format_res_limma, n_label = 10,
                gene_id_conversion = gene_id_conversion,
                log2FoldChange = log2FoldChange )

  list_DEP_limma <- list()
  list_DEP_limma <- lapply(res, function(x) {
    up=x$gene_name[which(x$adj.P.Val <= 0.05 & x$de=='Up')]
    down=x$gene_name[which(x$adj.P.Val <= 0.05 & x$de=='Down')]
    return(list(up=up,
                down=down))
  })

  res$plot <- lapply(res, volcano_plot_limma,
                     log2FoldChange =log2FoldChange )
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
