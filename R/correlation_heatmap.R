#' Correlation heatmap between Factors and covariates of interest
#'
#' @param integrated_object MOFA object
#' @param covariates covariates to be used for correlation
#'
#' @return Correlation heatmao
#'
#' @family Plotting 
#'
#' @importFrom stats cor
#' @importFrom psych corr.test
#' @importFrom corrplot corrplot
#' @import RColorBrewer
#' 
#' @export
correlation_heatmap <- function(integrated_object,
                                covariates = c(
                                  "AD", "Braak", "MTG",
                                  "SOM", "amyloid", "AT8", "PHF1"
                                )) {
  pr <- integrated_object@expectations[["Z"]][["group1"]]
  metadata <- integrated_object@samples_metadata[, covariates]

  metadata <- sapply(metadata, unclass)
  dim <- ncol(integrated_object@expectations$Z$group1)
  mat <- stats::cor(cbind(pr, metadata), method = c("pearson"),
                    use = "pairwise.complete.obs")
  pmat <- psych::corr.test(cbind(pr, metadata), method = c("pearson"))$p
  corrplot <- corrplot::corrplot(mat[(dim + 1):ncol(pmat), 1:dim],
    insig = "label_sig", p.mat = pmat[(dim + 1):ncol(pmat), 1:dim],
    sig.level = c(0.001, 0.01, 0.05), pch.cex = 1.3,
    col = colorRampPalette(c(
      "#2166AC", "#4393C3", "#92C5DE",
      "#D1E5F0", "#FDDBC7", "#F4A582", "#D6604D", "#B2182B"
    ))(50)
  )
  corrplot
  return(corrplot)
}



#' Correlation heatmap between Factors and covariates of interest
#'
#' @param integrated_object integrated object
#' @param covariates covariates to be used for correlation
#' @param metadata metadata
#'
#' @return Correlation heatmap
#'
#' @family Plotting 
#'
#' @importFrom stats cor
#' @importFrom psych corr.test
#' @importFrom corrplot corrplot
#'
#' @export
correlation_heatmap_supervised <- function(integrated_object,
                                           metadata,
                                           covariates = c("AD", "Braak", "amyloid", "pTau", "PHF1")) {
  pr <- integrated_object$variates$Y
  metadata <- metadata[, covariates]

  metadata <- sapply(metadata, as.numeric)
  dim <- ncol(pr)
  mat <- stats::cor(cbind(pr, metadata),
                    method = c("pearson"),
                    use = "pairwise.complete.obs"
  )
  pmat <- psych::corr.test(cbind(pr, metadata), method = c("pearson"))$p
  corrplot <- corrplot::corrplot(mat[(dim + 1):ncol(pmat), 1:dim],
                                 insig = "label_sig", p.mat = pmat[(dim + 1):ncol(pmat), 1:dim],
                                 sig.level = c(0.001, 0.01, 0.05), pch.cex = 1.3,
                                 col = colorRampPalette(c(
                                   "#2166AC", "#4393C3", "#92C5DE",
                                   "#D1E5F0", "#FDDBC7", "#F4A582", "#D6604D", "#B2182B"
                                 ))(50)
  )
  corrplot
  return(corrplot)
}




#' Correlation heatmap between Factors and covariates of interest
#'
#' @param integrated_object integrated object
#' @param covariates covariates to be used for correlation
#' @param metadata metadata
#'
#' @return Correlation heatmap
#'
#' @family Plotting 
#'
#' @importFrom stats cor
#' @importFrom psych corr.test
#' @importFrom corrplot corrplot
#'
#' @export
correlation_heatmap_clusters<- function(integrated_object,
                                        metadata,
                                        covariates = c("AD", "Braak", "amyloid", "pTau", "PHF1")) {
  clusters <- data.frame(cluster=as.factor(integrated_object$clusters))
  pr <- model.matrix(~ 0 + cluster, clusters)
  metadata <-  metadata[,covariates]
  metadata <- sapply(metadata, as.numeric)
  dim <- ncol(pr)
  mat <- stats::cor(cbind(pr, metadata),
                    method = c("pearson"),
                    use = "pairwise.complete.obs"
  )
  pmat <- psych::corr.test(cbind(pr, metadata), method = c("pearson"))$p
  corrplot <- corrplot::corrplot(mat[(dim + 1):ncol(pmat), 1:dim],
                                 insig = "label_sig", p.mat = pmat[(dim + 1):ncol(pmat), 1:dim],
                                 sig.level = c(0.001, 0.01, 0.05), pch.cex = 1.3,
                                 col = colorRampPalette(c(
                                   "#2166AC", "#4393C3", "#92C5DE",
                                   "#D1E5F0", "#FDDBC7", "#F4A582", "#D6604D", "#B2182B"
                                 ))(50)
  )
  corrplot
  return(corrplot)
}

