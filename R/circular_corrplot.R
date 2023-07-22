#' Circular heatmap of module-trait correlations and module-cell type enrichments
#'
#' @param moduleTraitCor  moduleTraitCor from `multiomics_modules()`
#' @param moduleColors moduleColors from `multiomics_modules()`
#' @param cell_type cell type plots from `cell_type_enrichment()`
#' @param names modulesLabels from `multiomics_modules()`
#'
#' @import circlize
#' @return Circular correlation plot
#' @family Plotting
#' @export
#'
circular_corrplot <- function(moduleTraitCor,
                              moduleColors,
                              cell_type,
                              names = NULL) {
  ## mat corr
  mat1 <- t(moduleTraitCor)
  col_fun1 <- circlize::colorRamp2(c(-1, 0, 1), c("red", "white", "blue"))

  # mat 0

  # mat modules
  mat0 <- data.frame(rownames(mat1))
  rownames(mat0) <- rownames(mat1)
  mat0$rownames.mat1. <- seq(1:length(mat0$rownames.mat1.))

  if (!is.null(names)) {
    rownames(mat0) <- names
  }

  col_fun0 <- list()

  for (i in moduleColors) {
    col_fun0[i] <- i
  }
  names(col_fun0) <- substring(rownames(mat1), 3)


  ## mat 2
  df <- do.call("rbind", cell_type)
  df <- df[, 1]

  df <- do.call("rbind", df)

  df$color <- sub("\\..*", "", rownames(df))

  df2 <- df[, c("q", "CellType", "color")]
  df3 <- stats::reshape(df2, idvar = "color", timevar = "CellType", direction = "wide")
  rownames(df3) <- df3$color
  df3$color <- NULL
  df3 <- abs(df3)
  rownames(df3) <- paste0("ME", rownames(df3))
  colnames(df3) <- stringr::str_remove(colnames(df3), "q.")
  col_fun2 <- circlize::colorRamp2(c(0, 0.05, 1), c("black", "yellow", "white"))

  merge_=merge(data.frame(mat1),df3,all.y=T)

  ### PLOT
  circlize::circos.clear()
  circlize::circos.par(gap.after = c(40))

  # mat 0
  # circos.par(start.degree = 30, gap.degree = 1, points.overflow.warning = FALSE)
  circlize::circos.heatmap(mat0,
    col = col_fun0,
    track.height = 0.05,
    bg.border = "gray50",
    rownames.side = "inside"
  )
  # mat 1
  circlize::circos.heatmap(mat1, col = col_fun1, track.height = 0.2, cell.border = "black", cluster = TRUE)
  circlize::circos.track(track.index = get.current.track.index(), panel.fun = function(x, y) {
    if (CELL_META$sector.numeric.index == 1) { # the last sector
      cn <- rev(colnames(mat1))
      n <- length(cn)
      circlize::circos.text(rep(CELL_META$cell.xlim[2], n) + convert_x(1, "mm"), # x coordinate
        (1:n), # y coordinate
        cn, # label
        cex = 0.5, adj = c(0, 1), facing = "inside"
      )
    }
  }, bg.border = NA)

  circlize::circos.heatmap(df3, col = col_fun2, track.height = 0.2, cell.border = "grey")
  circlize::circos.track(track.index = circlize::get.current.track.index(), panel.fun = function(x, y) {
    if (CELL_META$sector.numeric.index == 1) { # the last sector
      cn <- rev(colnames(df3))
      n <- length(cn)
      circlize::circos.text(rep(CELL_META$cell.xlim[2], n) + convert_x(1, "mm"), # x coordinate
        (1:n), # y coordinate
        cn, # label
        cex = 0.5, adj = c(0, 1), facing = "inside"
      )
    }
  }, bg.border = NA)


  lgd <- ComplexHeatmap::Legend(
    title = "Pearson's correlation",
    col_fun = col_fun1, at = c(-1, 0, 1)
  )
  lgd2 <- ComplexHeatmap::Legend(
    title = "BH-corrected q-values",
    col_fun = col_fun2, at = c(0, 0.05, 1)
  )
  lgd_list_vertical2 <- ComplexHeatmap::packLegend(lgd, lgd2)

  ComplexHeatmap::draw(lgd_list_vertical2,
    x = unit(50, "mm"),
    y = unit(1, "mm"), just = c("right", "bottom")
  )
}
