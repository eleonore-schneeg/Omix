#' Pseudotime_inference based on integrated embeddings and clusters
#'
#' @param model MOFA model from integration in `multiassay@metadata$integration$MOFA`
#' @param clusters cluster generated from MOFA2::cluster_samples or other
#' clustering method
#' @param time_factor MOFA factor corresponding to time
#' @param second_factor MOFA factor used as other axis of variation in
#' the embedding
#' @param lineage lineage
#' @param start.clus Starting cluser
#' @param end.clus Ending cluser
#'
#' @return A list object with `model` and `plot` slots
#'
#' @family Pseudotime inference
#'
#' @importFrom slingshot slingshot as.SlingshotDataSet
#' @importFrom viridis scale_color_viridis
#' @importFrom MOFA2 plot_factors samples_metadata
#' @import ggplot2
#' @export

#'

pseudotime_inference <- function(model,
                                 clusters,
                                 time_factor = 8,
                                 second_factor = 1,
                                 lineage = "Lineage1",
                                 start.clus = 1,
                                 end.clus=2) {
  gg <- MOFA2::plot_factors(model,
    factors = c(time_factor, second_factor),
    scale = FALSE
  )

  embeddings <- gg$data[, c("x", "y")]

  rownames(embeddings) <- gg$data[, c("sample")]
  embeddings <- embeddings[names(clusters$cluster), ]

  sds <- slingshot::slingshot(embeddings, clusters$cluster,
                              start.clus =start.clus,
                              end.clus=end.clus)
  df <- slingshot::as.SlingshotDataSet(sds)

  pseudotime <- as.data.frame(sds@assays@data@listData[["pseudotime"]])

  MOFA2::samples_metadata(model)$pseudotime <- pseudotime[paste(lineage)][match(
    names(clusters$cluster),
    rownames(pseudotime)
  ), ]

  embeddings <- data.frame(embeddings)
  embeddings$pseudotime <- pseudotime[, paste(lineage)]

  df2 <- data.frame(df@curves[[paste(lineage)]]$s)

  plot <- ggpubr::ggscatter(embeddings,
    x = "x", y = "y",
    color = "pseudotime"
  ) +
    viridis::scale_color_viridis() +
    ggplot2::geom_point(data = df2, ggplot2::aes(x = x, y = y)) +
    ggplot2::xlab(paste("Factor", time_factor, " ~ Pseudotime")) +
    ggplot2::ylab(paste("Factor", second_factor))

  return(list(
    model = model,
    plot = plot
  ))
}
