################################################################################
#' Formatting differential expression results from limma
#'
#' @param dt Differential expression results dataframe from limma
#' @param gene_id_conversion Default to NULL.
#' Can provide a dataframe with `gene_name`
#' and `uniprot_id` columns to easily convert nomenclature.
#' @param log2FoldChange Log2FC absolute threshold
#' @param n_label maximum number of label on volcano plot
#' @param padj adjusted p value threshold. Default to 0.05
#'
#' @return Formatted data frame
#'
#' @family Helper
#'
#' @importFrom dplyr filter top_n pull
#' @export
#'


format_res_limma <- function(dt,
                             gene_id_conversion = NULL,
                             log2FoldChange = 0,
                             n_label = 20,
                             padj = 0.05,
                             ...) {
  if (!is.null(gene_id_conversion)) {
    dt$gene_name <- gene_id_conversion$gene_name[
      match(dt$Identifier, gene_id_conversion$uniprot_id)]
  }

  if (is.null(gene_id_conversion)) {
    dt$gene_name <- dt$Identifier
  }

  dt$de[dt$adj.P.Val <= padj & dt$logFC > log2FoldChange] <- "Up"
  dt$de[dt$adj.P.Val <= padj & dt$logFC < -log2FoldChange] <- "Down"
  dt$de[dt$adj.P.Val > padj] <- "Not sig"
  dt$de <- factor(dt$de, levels = c("Up", "Down", "Not sig"))
  n_up <- sum(dt$de == "Up", na.rm = TRUE)
  n_down <- sum(dt$de == "Down", na.rm = TRUE)

  dt$label <- NA
  dt$log2FoldChange <- dt$logFC
  dt$padj <- dt$adj.P.Val
  dt$pvalue <- dt$P.Value


  if (n_up > 0) {
    top_up <- dt %>%
      dplyr::filter(de == "Up") %>%
      dplyr::top_n(min(n_up, n_label), wt = -padj) %>%
      dplyr::pull(Identifier)

    dt$label[dt$Identifier %in% top_up] <- "Yes"
  }
  if (n_down > 0) {
    top_down <- dt %>%
      dplyr::filter(de == "Down") %>%
      dplyr::top_n(min(n_down, n_label), wt = -padj) %>%
      dplyr::pull(Identifier)
    dt$label[dt$Identifier %in% top_down] <- "Yes"
  }

  dt <- dt[, c(
    "Identifier", "gene_name", "log2FoldChange", "pvalue", "padj", "de",
    "label"
  )]

  return(dt)
}



#' Formatting differential expression results from DESEQ2
#'
#' @param dt Differential expression results dataframe from DESEQ2
#' @param log2FoldChange  Log2FC absolute threshold
#' @param padj adjusted p value threshold. Default to 0.05
#' @param n_label maximum number of label on volcano plot
#' @param gene_id_conversion Default to NULL. Can provide a dataframe
#' with `gene_name` and `uniprot_id` columns to easily convert nomenclature.
#' @param filter_protein_coding Logical as whether to filter protein coding gene.
#' Default to TRUE.
#'
#' @return Formatted data frame
#'
#' @family Helper
#'
#' @importFrom dplyr mutate filter arrange top_n pull
#' @export
#'

format_res_deseq <- function(dt,
                             log2FoldChange = 0,
                             padj = 0.05,
                             n_label = 10,
                             gene_id_conversion = NULL,
                             filter_protein_coding = TRUE) {
  dt <- dt %>%
    as.data.frame() %>%
    dplyr::mutate(ensembl_name = rownames(.)) %>%
    dplyr::filter(!is.na(padj), !is.na(log2FoldChange)) %>%
    dplyr::arrange(padj)

  if (!is.null(gene_id_conversion)) {
    dt$gene_name <- gene_id_conversion$gene_name[match(
      dt$ensembl_name, gene_id_conversion$ensembl_gene_id)]
    dt$gene_type <- gene_id_conversion$gene_biotype[match(
      dt$ensembl_name, gene_id_conversion$ensembl_gene_id)]
  }


  dt$de[dt$padj <= padj & dt$log2FoldChange > log2FoldChange] <- "Up"
  dt$de[dt$padj <= padj & dt$log2FoldChange < -log2FoldChange] <- "Down"
  dt$de <- factor(dt$de, levels = c("Up", "Down", "Not sig"))
  n_up <- sum(dt$de == "Up" & dt$gene_type == "protein_coding", na.rm = TRUE)
  n_down <- sum(dt$de == "Down" & dt$gene_type == "protein_coding", na.rm = TRUE)

  dt$label <- NA
  if (filter_protein_coding == TRUE) {
    dt <- dt[which(dt$gene_type == "protein_coding"), ]
  }

  if (n_up > 0) {
    top_up <- dt %>%
      dplyr::filter(de == "Up" & gene_type == "protein_coding") %>%
      dplyr::top_n(min(n_up, n_label), wt = -padj) %>%
      dplyr::pull(gene_name)

    dt$label[dt$gene_name %in% top_up] <- "Yes"
  }
  if (n_down > 0) {
    top_down <- dt %>%
      dplyr::filter(de == "Down" & gene_type == "protein_coding") %>%
      dplyr::top_n(min(n_down, n_label), wt = -padj) %>%
      dplyr::pull(gene_name)
    dt$label[dt$gene_name %in% top_down] <- "Yes"
  }

  top_protein_coding_gene <- dt %>%
    dplyr::filter(label == "Yes") %>%
    dplyr::pull(gene_name)

  attr(dt, "top_protein_coding_gene") <- top_protein_coding_gene

  n_up <- dt %>%
    dplyr::filter(de == "Up") %>%
    dplyr::pull(gene_name) %>%
    length()

  n_down <- dt %>%
    dplyr::filter(de == "Down") %>%
    dplyr::pull(gene_name) %>%
    length()


  res_summary <- setNames(c(n_up, n_down), c("up", "down"))
  dt$gene_name <- sub("*\\.[0-9]", "", dt$gene_name)
  dt <- dt[, c(
    "ensembl_name", "gene_name", "log2FoldChange", "pvalue", "padj", "de",
    "label", "gene_type"
  )]

  attr(dt, "summary") <- res_summary
  return(dt)
}




#' Builds a volcano plot from limma results
#'
#' @param dt Formatted Differential expression results from `format_res_limma()`
#' @param log2FoldChange  Log2FC absolute threshold
#' @param padj adjusted p value threshold. Default to 0.05
#'
#' @return Volcano plot
#'
#' @family Plotting
#'
#' @import ggplot2
#' @importFrom ggrepel geom_text_repel
#' @export

volcano_plot_limma <- function(dt,
                               log2FoldChange = 0,
                               padj = 0.05) {
  ggplot2::ggplot(dt, ggplot2::aes(x = log2FoldChange, y = -log10(padj),
                          label = gene_name, colour = de, repel = TRUE)) +
    geom_point(aes(x = log2FoldChange, y = -log10(padj),
                   fill = de, colour = de), show.legend = T, alpha = 0.5) +
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
      aes(log2FoldChange, y = -log10(padj),
          label = ifelse(label == "Yes", gene_name, "")),
      max.overlaps = 2000
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
}


#' Builds a volcano plot from DESEQ2 results
#'
#' @param dt Formatted Differential expression results from `format_res_deseq()`
#' @param log2FoldChange  Log2FC absolute threshold
#' @param padj adjusted p value threshold. Default to 0.05
#'
#' @return Volcano plot
#'
#' @family Plotting
#'
#' @import ggplot2
#' @importFrom ggrepel geom_text_repel
#' @export
#'

volcano_plot_deseq <- function(dt,
                               log2FoldChange = 0,
                               padj = 0.05) {
  ggplot2::ggplot(dt, aes(x = log2FoldChange, y = -log10(padj),
                          label = gene_name, colour = de, repel = TRUE)) +
    geom_point(aes(x = log2FoldChange, y = -log10(padj),
                   fill = de, colour = de), show.legend = T, alpha = 0.5) +
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
      aes(log2FoldChange, y = -log10(padj),
          label = ifelse(label == "Yes", gene_name, "")),
      max.overlaps = 2000
    ) +
    xlab(bquote(Log[2] * " (fold-change)")) +
    ylab(bquote("-" * Log[10] * " (adjusted p-value)")) +
    geom_vline(
      xintercept = c(-log2FoldChange, log2FoldChange),
      linetype = 2, size = 0.2, alpha = 0.5
    ) +
    geom_hline(
      yintercept = -log10(0.05),
      linetype = 2, size = 0.2, alpha = 0.5
    )
}

#' Builds an interactive volcano plot using plotly
#'
#' @param data Formatted Differential expression results from
#' `format_res_deseq()` or `format_res_limma()`
#' @param log2FoldChange  Log2FC absolute threshold. Default to 0.25
#' @param padj adjusted p value threshold. Default to 0.05
#'
#' @return Interactive volcano plot
#'
#' @family Plotting
#'
#' @import ggplot2
#' @importFrom plotly ggplotly
#' @export
#'
volcano_interactive <- function(data,
                                log2FoldChange = 0.25,
                                padj = 0.05) {
  volcano <- ggplot2::ggplot(data, aes(x = log2FoldChange, y = -log10(padj),
                                       label = gene_name, colour = de,
                                       text = paste("Gene:", gene_name, "<br>",
    "Log2FC:", round(log2FoldChange, 3), "<br>",
    "Adjust pval;", round(padj, 3)
  ))) +
    geom_point() +
    geom_text(
      data = data,
      aes(log2FoldChange - 0.05, y = -log10(padj) + 0.1,
          label = ifelse(label == "Yes", gene_name, ""))
    ) +
    geom_vline(
      xintercept = c(-log2FoldChange, log2FoldChange),
      linetype = 2, size = 0.2, alpha = 0.5
    ) +
    geom_hline(
      yintercept = -log10(padj),
      linetype = 2, size = 0.2, alpha = 0.5
    ) +
    scale_color_manual(values = c(
      "Up" = "#DC0000FF",
      "Down" = "#3C5488FF",
      "Not sig" = "grey"
    ))

  plot <- plotly::ggplotly(volcano, tooltip = c("text"))
  plot
}



#' Builds an interactive comparison plot using plotly
#'
#' @param data Formatted Comparative Differential expression results
#' from `single_omic_comparisons()`
#'
#' @return Interactive plot
#'
#' @family Plotting
#'
#' @import ggplot2
#' @importFrom plotly ggplotly
#' @export
#'

volcano_interactive_comparison <- function(data) {
  volcano <- ggplot2::ggplot(data, aes(x = log2FoldChange.x,
                                       y = log2FoldChange.y,
                                       label = gene_name,
                                       colour = direction,
                                       text = paste(
    "Gene:", gene_name, "<br>",
    "Log2FC transcriptomics:", round(log2FoldChange.x, 3), "<br>",
    "Log2FC proteomics;", round(log2FoldChange.y, 3), "<br>",
    "pval transcriptomics:", round(pvalue.x, 3), "<br>",
    "pval proteomics:", round(pvalue.y, 3), "<br>",
    "padj transcriptomics:", round(padj.x, 3), "<br>",
    "padj proteomics:", round(padj.y, 3), "<br>"
  ))) +
    xlab("log2FoldChange Transcriptome") +
    ylab("log2FoldChange Proteome") +
    geom_point(size = 0.4) +
    geom_text(
      data = data,
      aes(log2FoldChange.x - 0.05, y = log2FoldChange.y + 0.1,
          label = ifelse(pvalue_category == "double_platform", gene_name, ""))
    ) +
    geom_vline(
      xintercept = 0,
      linetype = 2, size = 0.2, alpha = 0.5
    ) +
    geom_hline(
      yintercept = 0,
      linetype = 2, size = 0.2, alpha = 0.5
    ) +
    scale_color_manual(values = c(
      "Concordant" = "#3C5488FF",
      "Discordant" = "#DC0000FF"
    ))

  plot <- plotly::ggplotly(volcano, tooltip = c("text"))
  plot
}
