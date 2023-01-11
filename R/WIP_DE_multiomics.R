#' Performs DE on multiomics object
#'
#' @param multiassay
#' @param dependent
#' @param covariates
#' @param levels
#'
#' @return
#' @importFrom MultiAssayExperiment MultiAssayExperiment listToMap colData
#'  getWithColData sampleMap
#' @importFrom SummarizedExperiment SummarizedExperiment
#' @importFrom stats as.formula
#' @importFrom magrittr set_names
#' @importFrom cli cli_alert_danger style_bold cli_alert_success
#' @importFrom limma makeContrasts lmFit contrasts.fit eBayes topTable
#'
#' @export
#'
#' @examples
DE_multiomics <- function(multiassay,
                          dependent = "diagnosis",
                          covariates = c("pmi", "msex", "age_death"),
                          levels = c("Control", "Mid", "Late")) {
  count <- multiassay@metadata[["dds"]]@assays@data@listData[["counts"]]
  count <- as.data.frame(count)
  multiassay_temp <- counts_rna(multiassay, custom_processed_df = count)

  multimodal <- .get_multimodal_object(multiassay_temp,
    slots = c(
      "rna_counts",
      "protein_processed"
    )
  )
  # sanity check
  # colnames(multimodal$multimodal_object$rna_counts) == rownames(multimodal$metadata)
  rna <- multimodal$multimodal_object$rna_counts
  protein <- multimodal$multimodal_object$protein_processed
  colData <- multimodal$metadata

  ####
  colData <- colData %>% mutate(across(where(is.character), as.factor))

  if (!is.null(covariates)) {
    cov_var <- paste(covariates, collapse = " + ")

    if (any(is.na(colData[covariates]))) {
      stop(cli::cli_alert_danger(
        paste("Missing value in covariates, please use other covariates")
      ))
    }
  } else {
    cov_var <- NULL
  }

  dependent_var <- paste(dependent, collapse = " + ")

  if (!is.null(levels)) {
    colData[[dependent]] <- factor(colData[[dependent]], levels = levels)
  }

  for (i in colnames(colData)) {
    assign(i, colData[, i])
  }


  ### CATEGORICAL DEPENDENT
  if (!is.numeric(colData[, dependent])) {
    if (!is.na(cov_var)) {
      model_formula <- stats::as.formula(
        sprintf(
          "~0+ %s + %s",
          dependent_var,
          cov_var
        )
      )
    }
    if (is.na(cov_var)) {
      model_formula <- stats::as.formula(
        sprintf(
          "~0+ %s",
          dependent_var
        )
      )
    }

    design <- model.matrix(model_formula)


    if (!is.null(cov_var)) {
      colnames(design) <- c(levels(colData[, dependent]), covariates)
    } else {
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

    contr.matrix <- limma::makeContrasts(contrasts = c(x, covariates), levels = design)

    # RNA
    d0 <- edgeR::DGEList(rna)
    d0 <- edgeR::calcNormFactors(d0)
    y <- limma::voom(d0, design)
    fit <- limma::lmFit(y, design)
    vfit <- limma::contrasts.fit(fit, contrasts = contr.matrix)
    efit <- limma::eBayes(vfit)

    df_rna <- list()
    for (i in 1:length(perCombs)) {
      df_rna[[i]] <- assign(paste0(
        "limma.results_",
        colnames(limma::topTable(efit))[i]
      ), limma::topTable(efit, coef = i, n = Inf))
    }

    names(df_rna) <- colnames(topTable(efit))[1:length(perCombs)]
    names(df_rna) <- gsub(pattern = ".", replacement = "vs", x = names(df_rna), fixed = TRUE)

    # PROTEIN
    fit <- limma::lmFit(protein, design)
    vfit <- limma::contrasts.fit(fit, contrasts = contr.matrix)
    efit <- limma::eBayes(vfit)
    df_protein <- list()
    for (i in 1:length(perCombs)) {
      df_protein[[i]] <- assign(paste0(
        "limma.results_",
        colnames(limma::topTable(efit))[i]
      ), limma::topTable(efit, coef = i, n = Inf))
    }

    names(df_protein) <- colnames(topTable(efit))[1:length(perCombs)]
    names(df_protein) <- gsub(pattern = ".", replacement = "vs", x = names(df_protein), fixed = TRUE)
  }

  ### NUMERIC DEPENDENT
  if (is.numeric(colData[, dependent])) {
    if (!is.na(cov_var)) {
      model_formula <- stats::as.formula(
        sprintf(
          "~%s + %s",
          dependent_var,
          cov_var
        )
      )
    }

    if (is.na(cov_var)) {
      model_formula <- stats::as.formula(
        sprintf(
          "~%s",
          dependent_var
        )
      )
    }


    design <- model.matrix(model_formula, data = colData)
    # RNA
    d0 <- edgeR::DGEList(rna)
    d0 <- edgeR::calcNormFactors(d0)
    y <- limma::voom(d0, design)
    fit <- limma::lmFit(y, design)
    efit <- limma::eBayes(fit)
    df_rna <- limma::topTable(efit, coef = paste(dependent), n = Inf, adjust = "BH")

    # PROTEIN
    fit <- limma::lmFit(protein, design)
    efit <- limma::eBayes(fit)
    df_protein <- limma::topTable(efit, coef = paste(dependent), n = Inf, adjust = "BH")
  }


  DE_list <- list(
    rna = df_rna,
    protein = df_protein
  )

  # DE_list=lapply(DE_list,format_DE)

  return(DE_list)
}


format_volcano <- function(dt,
                           padj,
                           log2FoldChange,
                           n_label) {
  dt$de[dt$adj.P.Val <= padj & dt$logFC > log2FoldChange] <- "Up"
  dt$de[dt$adj.P.Val <= padj & dt$logFC < -log2FoldChange] <- "Down"
  dt$de[dt$adj.P.Val > padj] <- "Not sig"
  dt$de <- factor(dt$de, levels = c("Up", "Down", "Not sig"))
  n_up <- sum(dt$de == "Up", na.rm = TRUE)
  n_down <- sum(dt$de == "Down", na.rm = TRUE)

  dt$label <- NA
  dt$gene_name <- rownames(dt)
  dt$gene_name <- sub("*\\.[0-9]", "", dt$gene_name)
  dt <- dt %>%
    dplyr::select(gene_name, everything())

  if (n_up > 0) {
    top_up <- dt %>%
      dplyr::filter(de == "Up") %>%
      dplyr::top_n(min(n_up, n_label), wt = -padj) %>%
      pull(gene_name)

    dt$label[which(dt$gene_name %in% top_up)] <- "Yes"
  }
  if (n_down > 0) {
    top_down <- dt %>%
      dplyr::filter(de == "Down") %>%
      dplyr::top_n(min(n_down, n_label), wt = -padj) %>%
      pull(gene_name)
    dt$label[which(dt$gene_name %in% top_down)] <- "Yes"
  }

  ### VOLCANO
  volcano <- ggplot2::ggplot(dt, aes(x = logFC, y = -log10(adj.P.Val),
                                     colour = de)) +
    geom_point(aes(x = logFC, y = -log10(adj.P.Val), fill = de, colour = de),
               show.legend = T, alpha = 0.5) +
    theme_classic() +
    scale_colour_manual(
      name = NULL,
      aesthetics = c("colour", "fill"),
      values = c("#DC0000FF", "#3C5488FF", "grey"),
      label = c("Up-regulated", "Down-regulated", "Not significant"),
      breaks = c("Up", "Down", "Not sig")
    ) +
    ggrepel::geom_text_repel(
      data = dt,
      aes(logFC, y = -log10(adj.P.Val), label = ifelse(label == "Yes",
                                                       gene_name, "")),
                                                       max.overlaps = 20
    ) +
    xlab(bquote(Log[2] * " (fold-change)")) +
    ylab(bquote("-" * Log[10] * " (adjusted p-value)")) +
    geom_vline(
      xintercept = c(-log2FoldChange, log2FoldChange),
      linetype = 2, size = 0.2, alpha = 0.5
    ) +
    geom_hline(
      yintercept = -log10(padj),
      linetype = 2, size = 0.2, alpha = 0.5
    )
  volcano
}


format_DE <- function(dt,
                      adj.P.Val = 0.05,
                      log2FoldChange = 0.2) {
  dt$de[dt$adj.P.Val <= adj.P.Val & dt$logFC > log2FoldChange] <- "Up"
  dt$de[dt$adj.P.Val <= adj.P.Val & dt$logFC < -log2FoldChange] <- "Down"
  dt$de[dt$adj.P.Val > adj.P.Val] <- "Not sig"
  dt$de <- factor(dt$de, levels = c("Up", "Down", "Not sig"))
  n_up <- sum(dt$de == "Up", na.rm = TRUE)
  n_down <- sum(dt$de == "Down", na.rm = TRUE)

  dt$label <- NA
  dt$gene_name <- rownames(dt)
  dt$gene_name <- sub("*\\.[0-9]", "", dt$gene_name)
  dt <- dt %>%
    dplyr::select(gene_name, everything())
  return(dt)
}


single_omic_comparison <- function(DEG,
                                   DEP,
                                   filtering = FALSE,
                                   filtering_options = c("both_significant", "either"),
                                   pvalue = c("padj", "pval")) {
  t <- merge(DEG, DEP, by = "gene_name")
  t$direction <- ifelse(t$logFC.x >= 0 & t$logFC.y >= 0 |
    t$logFC.x <= 0 & t$logFC.y <= 0,
  "Concordant", "Discordant"
  )

  if (pvalue == "padj") {
    t$pvalue_category <- ifelse(t$adj.P.Val.x <= 0.05 & t$adj.P.Val.y <= 0.05, "double_platform",
      ifelse(t$adj.P.Val.x <= 0.05 & t$adj.P.Val.y > 0.05, "transcriptomics_only",
        ifelse(t$adj.P.Val.x > 0.05 & t$adj.P.Val.y <= 0.05, "proteomics_only", "never_significant")
      )
    )
  }

  if (pvalue == "pval") {
    t$pvalue_category <- ifelse(t$P.Value.x <= 0.05 & t$P.Value.y <= 0.05, "double_platform",
      ifelse(t$P.Value.x <= 0.05 & t$P.Value.y > 0.05, "transcriptomics_only",
        ifelse(t$P.Value.x > 0.05 & t$P.Value.y <= 0.05, "proteomics_only", "never_significant")
      )
    )
  }

  if (filtering != TRUE) {
    comparison_volcano <- ggpubr::ggscatter(t,
      x = "logFC.x",
      y = "logFC.y", label = "gene_name",
      repel = TRUE, color = "direction",
      font.label = c(8, "plain")
    ) +
      xlab("log2FoldChange Transcriptome") +
      ylab("log2FoldChange Proteome") +
      geom_vline(
        xintercept = 0, linetype = "dotted",
        color = "grey", size = 1
      ) +
      geom_hline(
        yintercept = 0, linetype = "dotted",
        color = "grey", size = 1
      ) +
      theme(
        axis.text.x = element_text(face = "bold", size = 10),
        axis.text.y = element_text(face = "bold", size = 10),
        axis.text = element_text(size = 200)
      ) + theme_classic(base_size = 15)
  }

  if (filtering == TRUE) {
    if (filtering_options == "both_significant") {
      t_double <- t[which(t$pvalue_category == "double_platform"), ]
    }
    if (filtering_options == "either") {
      t_double <- t[which(t$pvalue_category != "never_significant"), ]
    }
    comparison_volcano <- ggpubr::ggscatter(t_double,
      x = "logFC.x",
      y = "logFC.y", label = "gene_name",
      repel = TRUE, color = "direction",
      font.label = c(8, "plain")
    ) +
      xlab("log2FoldChange Transcriptome") +
      ylab("log2FoldChange Proteome") +
      geom_vline(
        xintercept = 0, linetype = "dotted",
        color = "grey", size = 1
      ) +
      geom_hline(
        yintercept = 0, linetype = "dotted",
        color = "grey", size = 1
      ) +
      theme(
        axis.text.x = element_text(face = "bold", size = 10),
        axis.text.y = element_text(face = "bold", size = 10),
        axis.text = element_text(size = 200)
      ) + theme_classic(base_size = 15)
  }
  comparison_volcano
}


single_omic_comparison_df <- function(DEG,
                                      DEP,
                                      filtering = FALSE,
                                      filtering_options = c("both_significant", "either"),
                                      pvalue = c("padj", "pval")) {
  t <- merge(DEG, DEP, by = "gene_name")

  t$direction <- ifelse(t$logFC.x >= 0 & t$logFC.y >= 0 |
    t$logFC.x <= 0 & t$logFC.y <= 0,
  "Concordant", "Discordant"
  )

  if (pvalue == "padj") {
    t$pvalue_category <- ifelse(t$adj.P.Val.x <= 0.05 & t$adj.P.Val.y <= 0.05, "double_platform",
      ifelse(t$adj.P.Val.x <= 0.05 & t$adj.P.Val.y > 0.05, "transcriptomics_only",
        ifelse(t$adj.P.Val.x > 0.05 & t$adj.P.Val.y <= 0.05, "proteomics_only", "never_significant")
      )
    )
  }

  if (pvalue == "pval") {
    t$pvalue_category <- ifelse(t$P.Value.x <= 0.05 & t$P.Value.y <= 0.05, "double_platform",
      ifelse(t$P.Value.x <= 0.05 & t$P.Value.y > 0.05, "transcriptomics_only",
        ifelse(t$P.Value.x > 0.05 & t$P.Value.y <= 0.05, "proteomics_only", "never_significant")
      )
    )
  }

  if (filtering != TRUE) {
    return(t)
  }
  if (filtering == TRUE) {
    if (filtering_options == "both_significant") {
      t_double <- t[which(t$pvalue_category == "double_platform"), ]
    }
    if (filtering_options == "either") {
      t_double <- t[which(t$pvalue_category != "never_significant"), ]
    }

    return(t_double)
  }
}
