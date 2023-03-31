#' Extract_weigths from MOFA model for a given factor of interest at a
#' defined absolute threshold
#'
#' @param model MOFA model from integration in
#' `multiassay@metadata$integration$MOFA`
#' @param factor Factor to extract weight from
#' @param threshold Absolute threshold to filter weights
#' @param sense_check_variable default to `NULL`. High weights should
#' coincide with stronger correlation if the sense_check_variable is an
#' important driver of variation in the designated factor. Will be used to
#' generate a plot if set to a covariate.
#'
#' @return List with positive and negative weights above and below the
#' specified threshold, distirbution plots of the weights, and dataframes
#' of the weights for each omic layer.
#'
#' @family Helper
#'
#' @export
#'
extract_weigths <- function(model,
                            factor = 1,
                            threshold = 0.3,
                            sense_check_variable = NULL) {
  weights_rna_1 <- MOFA2::get_weights(model,
    view = "mRNA",
    factor = factor,
    abs = FALSE,
    scale = TRUE,
    as.data.frame = FALSE
  )


  weights_prot_1 <- MOFA2::get_weights(model,
    view = "proteins",
    factor = factor,
    abs = FALSE,
    scale = TRUE,
    as.data.frame = FALSE
  )

  rna_1 <- data.frame(weights = weights_rna_1)
  protein_1 <- data.frame(weights = weights_prot_1)

  rna_1$Row.names <- sub("\\_.*", "", rownames(rna_1))
  protein_1$Row.names <- sub("\\_.*", "", rownames(protein_1))

  colnames(rna_1) <- c("Weights", "Feature")
  colnames(protein_1) <- c("Weights", "Feature")

  rna_1$model_feature <- rownames(rna_1)
  protein_1$model_feature <- rownames(protein_1)

  protein_positive <- protein_1 %>%
    arrange(desc(Weights), desc(Feature)) %>%
    filter(Weights >= threshold) %>%
    pull(Feature)

  rna_positive <- rna_1 %>%
    arrange(desc(Weights), desc(Feature)) %>%
    filter(Weights >= threshold) %>%
    pull(Feature)

  protein_negative <- protein_1 %>%
    arrange(desc(Weights), desc(Feature)) %>%
    filter(Weights <= -threshold) %>%
    pull(Feature)

  rna_negative <- rna_1 %>%
    arrange(desc(Weights), desc(Feature)) %>%
    filter(Weights <= -threshold) %>%
    pull(Feature)

  Weights_up <- list(
    rna = unique(rna_positive),
    protein = unique(protein_positive)
  )
  Weights_down <- list(
    rna = unique(rna_negative),
    protein = unique(protein_negative)
  )

  rna_1 <- rna_1 %>%
    arrange(desc(Weights), desc(Feature))

  protein_1 <- protein_1 %>%
    arrange(desc(Weights), desc(Feature))

  ## Check distribution
  distrib_rna <- ggplot2::ggplot(rna_1, aes(x = Weights)) +
    geom_density() +
    theme_classic() +
    geom_vline(
      xintercept = threshold,
      linetype = "dotted",
      color = "blue",
      linewidth = 0.8
    ) +
    geom_vline(
      xintercept = -threshold, linetype = "dotted",
      color = "blue", linewidth = 0.8
    ) +
    ylab("Density of RNA weights")

  distrib_protein <- ggplot2::ggplot(protein_1, aes(x = Weights)) +
    geom_density() +
    theme_classic() +
    geom_vline(
      xintercept = threshold,
      linetype = "dotted",
      color = "blue",
      linewidth = 0.8
    ) +
    geom_vline(
      xintercept = -threshold, linetype = "dotted",
      color = "blue", size = 0.8
    ) +
    ylab("Density of Proteins weights")


  ggplot2::ggplot(protein_1, aes(x = correlation)) +
    geom_density() +
    theme_classic() +
    geom_vline(
      xintercept = threshold,
      linetype = "dotted",
      color = "blue",
      linewidth = 0.8
    ) +
    geom_vline(
      xintercept = -threshold, linetype = "dotted",
      color = "blue", linewidth = 0.8
    ) +
    ylab("Density of Proteins weights")


  #### Check weight/correlation to sense check variable
  if (!is.null(sense_check_variable)) {
    metadata <- samples_metadata(model)

    t_p <- t(data.frame(model@data$proteins))
    cor_p <- stats::cor(metadata[, sense_check_variable], t_p,
                        use = "pairwise.complete.obs")
    protein_1$Correlation <- cor_p[match(protein_1$model_feature,
                                         colnames(cor_p))]

    t_r <- t(data.frame(model@data$mRNA))
    cor_r <- stats::cor(metadata[, sense_check_variable], t_r,
                        use = "pairwise.complete.obs")
    rna_1$Correlation <- cor_r[match(rna_1$model_feature,
                                     colnames(cor_r))]

    cor_weights_rna <- ggpubr::ggscatter(
      rna_1,
      "Weights",
      "Correlation"
    ) + ggplot2::geom_vline(
      xintercept = c(-threshold, threshold),
      linetype = "dotted",
      color = "blue",
      linewidth = 0.8
    ) +
      ggplot2::geom_hline(
        yintercept = c(-threshold, threshold),
        linetype = "dotted",
        color = "blue",
        linewidth = 0.8
      ) +
      viridis::scale_color_viridis() +
      ggplot2::ylab(paste(
        "Correlation of each transcript with",
        sense_check_variable
      ))

    cor_weights_protein <- ggpubr::ggscatter(
      protein_1,
      "Weights",
      "Correlation"
    ) + ggplot2::geom_vline(
      xintercept = c(-threshold, threshold),
      linetype = "dotted",
      color = "blue",
      linewidth = 0.8
    ) +
      ggplot2::geom_hline(
        yintercept = c(-threshold, threshold),
        linetype = "dotted",
        color = "blue",
        linewidth = 0.8
      ) +
      viridis::scale_color_viridis() +
      ggplot2::ylab(paste(
        "Correlation of each protein with",
        sense_check_variable
      ))
  }

  if (!is.null(sense_check_variable)) {
    return(result = list(
      weights = list(
        ranked_weights_positive = Weights_up,
        ranked_weights_negative = Weights_down
      ),
      distribution_plot = list(
        rna = distrib_rna,
        protein = distrib_protein
      ),
      weights_cor_plot = list(
        rna = cor_weights_rna,
        protein = cor_weights_protein
      ),
      weights_df = list(
        rna = rna_1,
        protein = protein_1
      )
    ))
  }

  if (is.null(sense_check_variable)) {
    return(result = list(
      weights = list(
        ranked_weights_positive = Weights_up,
        ranked_weights_negative = Weights_down
      ),
      distribution_plot = list(
        rna = distrib_rna,
        protein = distrib_protein
      ),
      weights_df = list(
        rna = rna_1,
        protein = protein_1
      )
    ))
  }
}
