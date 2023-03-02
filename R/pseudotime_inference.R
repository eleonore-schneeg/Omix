#' pseudotime_inference
#'
#' @param model
#' @param clusters
#' @param time_factor
#' @param second_factor
#'
#' @return
#' @export
#'
#' @examples
pseudotime_inference <- function(model,
                                 clusters,
                                 time_factor = 8,
                                 second_factor = 1,
                                 pseudotime_var = "pseudotime", ...) {
  gg <- MOFA2::plot_factors(model,
    factors = c(time_factor, second_factor),
    color_by = pseudotime_var,
    scale = F
  )

  embeddings <- gg$data[, c("x", "y")]

  rownames(embeddings) <- gg$data[, c("sample")]
  embeddings <- embeddings[names(clusters$cluster), ]

  sds <- slingshot::slingshot(embeddings, clusters$cluster,...)
  df <- slingshot::as.SlingshotDataSet(sds)

  pseudotime <- as.data.frame(sds@assays@data@listData[["pseudotime"]])

  MOFA2::samples_metadata(model)$pseudotime <- pseudotime$Lineage1[match(
    names(clusters$cluster),
    rownames(pseudotime)
  )]

  embeddings <- data.frame(embeddings)
  embeddings$pseudotime <- pseudotime$Lineage1

  df2 <- data.frame(df@curves$Lineage1$s)

  plot <- ggpubr::ggscatter(embeddings,
    x = "x", y = "y",
    color = paste(pseudotime_var)
  ) +
    viridis::scale_color_viridis() +
    geom_point(data = df2, aes(x = x, y = y)) +
    xlab(paste("Factor", time_factor, " ~ Pseudotime")) +
    ylab(paste("Factor", second_factor, " ~ Neuropathology"))

  return(list(
    model = model,
    plot = plot
  ))
}
