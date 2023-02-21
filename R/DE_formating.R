################################################################################
#' Title
#'
#' @param dt
#' @param gene_id_conversion
#' @param log2FoldChange
#' @param n_label
#' @param padj
#' @param ylim
#'
#' @return
#' @importFrom dplyr filter arrange pull top_n
#' @export
#'
#' @examples

format_res_limma <- function(dt,
                             gene_id_conversion = NULL,
                             log2FoldChange = 0,
                             n_label = 20,
                             padj = 0.05,
                             ylim = c(0, 10)) {

  if (!is.null(gene_id_conversion)) {
    dt$gene_name <- gene_id_conversion$gene_name[match(dt$Identifier, gene_id_conversion$uniprot_id)]
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


  if (n_up > 0) {
    top_up <- dt %>%
      dplyr::filter(de == "Up") %>%
      dplyr::top_n(min(n_up, n_label), wt = -padj) %>%
      pull(Identifier)

    dt$label[dt$Identifier %in% top_up] <- "Yes"
  }
  if (n_down > 0) {
    top_down <- dt %>%
      dplyr::filter(de == "Down") %>%
      dplyr::top_n(min(n_down, n_label), wt = -padj) %>%
      pull(Identifier)
    dt$label[dt$Identifier %in% top_down] <- "Yes"
  }
  return(dt)
}



#' Title
#'
#' @param dt
#' @param log2FoldChange
#' @param padj
#' @param n_label
#' @param gene_id_conversion
#' @param filter_protein_coding
#'
#' @return
#' @importFrom dplyr filter arrange pull
#' @export
#'
#' @examples
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
    dt$gene_name <- gene_id_conversion$gene_name[match(dt$ensembl_name, gene_id_conversion$ensembl_gene_id)]
    dt$gene_type <- gene_id_conversion$gene_biotype[match(dt$ensembl_name, gene_id_conversion$ensembl_gene_id)]
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

  attr(dt, "summary") <- res_summary
  return(dt)
}




#' Title
#'
#' @param dt
#' @param log2FoldChange
#' @param n_label
#' @param padj
#' @param ylim
#'
#' @return
#' @export
#'
#' @examples
volcano_plot_limma <- function(dt,
                               log2FoldChange = 0,
                               n_label = 10,
                               padj = 0.05,
                               ylim = c(0, 10)) {
  library(ggplot2)

  ggplot2::ggplot(dt, aes(x = logFC, y = -log10(adj.P.Val), label = gene_name, colour = de, repel = TRUE)) +
    geom_point(aes(x = logFC, y = -log10(adj.P.Val), fill = de, colour = de), show.legend = T, alpha = 0.5) +
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
      aes(logFC, y = -log10(adj.P.Val), label = ifelse(label == "Yes", gene_name, "")), max.overlaps = 2000
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


#' Title
#'
#' @param dt
#' @param log2FoldChange
#' @param n_label
#' @param padj
#' @param ylim
#'
#' @return
#' @export
#'
#' @examples
volcano_plot_deseq <- function(dt,
                               log2FoldChange = 0,
                               n_label = 10,
                               padj = 0.05,
                               ylim = c(0, 10)) {
  library(ggplot2)

  ggplot2::ggplot(dt, aes(x = log2FoldChange, y = -log10(padj), label = gene_name, colour = de, repel = TRUE)) +
    geom_point(aes(x = log2FoldChange, y = -log10(padj), fill = de, colour = de), show.legend = T, alpha = 0.5) +
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
      aes(log2FoldChange, y = -log10(padj), label = ifelse(label == "Yes", gene_name, "")), max.overlaps = 2000
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

volcano_interactive <- function(data, log2FoldChange = 0.25) {
  volcano <- ggplot2::ggplot(data, aes(x = log2FoldChange, y = -log10(padj), label = gene_name, colour = de, text = paste(
    "Gene:", gene_name, "<br>",
    "Log2FC:", log2FoldChange, "<br>",
    "Adjust pval;", padj
  ))) +
    geom_point() +
    geom_text(
      data = data,
      aes(log2FoldChange - 0.05, y = -log10(padj) + 0.6, label = ifelse(label == "Yes", gene_name, ""))
    ) +
    geom_vline(
      xintercept = c(-log2FoldChange, log2FoldChange),
      linetype = 2, size = 0.2, alpha = 0.5
    ) +
    geom_hline(
      yintercept = -log10(0.05),
      linetype = 2, size = 0.2, alpha = 0.5
    ) +
    scale_color_manual(values = c(
      "Up" = "#DC0000FF",
      "Down" = "#3C5488FF",
      "Not sig" = "grey"
    ))

  plot <- ggplotly(volcano, tooltip = c("text"))
  plot
}
