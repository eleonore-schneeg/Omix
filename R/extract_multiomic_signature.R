#' Title
#'
#' @param multiassay
#' @param integration
#' @param component
#'
#' @return
#' @export
#'
#' @examples
extract_multiomic_signature <- function(multiassay,
                                        integration = "DIABLO",
                                        component = 1) {
  model <- multiassay@metadata$integration[[paste(integration)]]
  matrix <- cbind(model$X[[1]], model$X[[2]])

  loadings <- mixOmics::selectVar(model, comp = component)
  Loadings <- c(
    loadings[["mRNA"]][["name"]],
    loadings[["proteins"]][["name"]]
  )

  df1 <- tibble::rownames_to_column(loadings$mRNA$value, "feature")
  df1$`Omic.layer` <- "rna"
  df2 <- tibble::rownames_to_column(loadings$proteins$value, "feature")
  df2$`Omic.layer` <- "protein"

  data_loadings <- rbind(df1, df2)
  data_loadings_positive <- data_loadings[data_loadings$value.var > 0, ]
  data_loadings_negative <- data_loadings[data_loadings$value.var < 0, ]

  rna_positive <- data_loadings_positive$feature[data_loadings_positive$Omic.layer == "rna"]
  rna_negative <- data_loadings_negative$feature[data_loadings_negative$Omic.layer == "rna"]
  protein_positive <- data_loadings_positive$feature[data_loadings_positive$Omic.layer == "protein"]
  protein_negative <- data_loadings_negative$feature[data_loadings_negative$Omic.layer == "protein"]

  if (integration == "sMBPLS") {
    if (model[["loadings"]][["Y"]][1] == 1) {
      detrimental <- list(
        rna = rna_positive,
        protein = protein_positive
      )

      protective <- list(
        rna = rna_negative,
        protein = protein_negative
      )
    }
    if (model[["loadings"]][["Y"]][1] != 1) {
      detrimental <- list(
        rna = rna_negative,
        protein = protein_negative
      )

      protective <- list(
        rna = rna_positive,
        protein = protein_positive
      )
    }
  }
  if (integration == "DIABLO") {
    if (model[["loadings"]][["Y"]][1] > 0) {
      detrimental <- list(
        rna = rna_negative,
        protein = protein_negative
      )

      protective <- list(
        rna = rna_positive,
        protein = protein_positive
      )
    }
    if (model[["loadings"]][["Y"]][1] < 0) {
      detrimental <- list(
        rna = rna_negative,
        protein = protein_negative
      )

      protective <- list(
        rna = rna_positive,
        protein = protein_positive
      )
    }
  }

  list <- list(
    protective = protective,
    detrimental = detrimental
  )


  return(list)
}
