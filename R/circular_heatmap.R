#' Title
#'
#' @param network_communities The name of the network communities
#' @param cell_type_communities Cell type enrichment of the community
#'
#' @export
#'
circular_heatmap <- function(network_communities,
                             cell_type_communities,
                             data =
                               "bundle_WGCNA_rna/rna-networkConstruction-auto_new.RData") {
  # compute PC1 eigen gene of each community

  int_results_clustering$Up$network_communities
  int_results_clustering$Up$cell_type_communities$`CS1-CS2`$results
  # correlated with covariates
  # cell type result


  lnames <- load(paste(data))


  nGenes <- ncol(datExpr)
  library(ComplexHeatmap)
  library(circlize) # >= 0.4.10

  nSamples <- nrow(datExpr)
  # Recalculate MEs with color labels
  MEs0 <- moduleEigengenes(datExpr, moduleColors)$eigengenes
  MEs <- orderMEs(MEs0)
  moduleTraitCor <- cor(MEs, datTraits, use = "p")
  rownames(moduleTraitCor) <- substring(rownames(moduleTraitCor), 3)
  moduleTraitPvalue <- corPvalueStudent(moduleTraitCor, nSamples)



  ## mat corr
  colnames(moduleTraitCor)
  mat1 <- moduleTraitCor[, c("braaksc", "amyloid", "gpath", "nft", "activated_microglia", "cogng_path_slope")]
  col_fun1 <- colorRamp2(c(-0.4, 0, 0.4), c("blue", "white", "red"))


  # mat 0

  # mat modules
  mat0 <- data.frame(rownames(mat1))
  rownames(mat0) <- rownames(mat1)
  mat0$rownames.mat1. <- seq(1:length(mat0$rownames.mat1.))



  ## annotate
  df0 <- do.call("rbind", enrichment)
  df0 <- df0[, 3]

  # pick top
  top <- lapply(df0, function(x) {
    top_n(x, 1, odds_ratio)
  })
  top2 <- lapply(top, function(x) {
    x[1, ]
  })
  top3 <- lapply(top2, function(x) {
    x$description
  })
  top4 <- unlist(top3)
  mat0$process <- top4


  col_fun0 <- list()
  for (i in modNames) {
    col_fun0[i] <- i
  }
  names(col_fun0) <- substring(rownames(mat1), 3)


  ## mat 2
  df <- do.call("rbind", cell_type_enrichment$plots)
  df <- df[, 1]

  df <- do.call("rbind", df)

  df$color <- sub("\\..*", "", rownames(df))

  df2 <- df[, c("q", "CellType", "color")]
  df3 <- reshape(df2, idvar = "color", timevar = "CellType", direction = "wide")
  rownames(df3) <- df3$color
  df3$color <- NULL
  df3 <- abs(df3)
  rownames(df3) <- paste0("ME", rownames(df3))
  colnames(df3) <- str_remove(colnames(df3), "q.")
  col_fun2 <- colorRamp2(c(0, 0.05, 1), c("black", "yellow", "white"))


  ### PLOT
  circos.clear()
  circos.par(gap.after = c(40))

  # mat 0
  # circos.par(start.degree = 30, gap.degree = 1, points.overflow.warning = FALSE)
  circos.heatmap(mat0,
    col = col_fun0,
    track.height = 0.05,
    bg.border = "gray50",
    rownames.side = "inside"
  )
  # mat 1
  circos.heatmap(mat1, col = col_fun1, track.height = 0.2, cell.border = "black", cluster = TRUE)
  circos.track(track.index = get.current.track.index(), panel.fun = function(x, y) {
    if (CELL_META$sector.numeric.index == 1) { # the last sector
      cn <- rev(colnames(mat1))
      n <- length(cn)
      circos.text(rep(CELL_META$cell.xlim[2], n) + convert_x(1, "mm"), # x coordinate
        (1:n), # y coordinate
        cn, # label
        cex = 0.5, adj = c(0, 1), facing = "inside"
      )
    }
  }, bg.border = NA)

  circos.heatmap(df3, col = col_fun2, track.height = 0.2, cell.border = "grey")
  circos.track(track.index = get.current.track.index(), panel.fun = function(x, y) {
    if (CELL_META$sector.numeric.index == 1) { # the last sector
      cn <- rev(colnames(df3))
      n <- length(cn)
      circos.text(rep(CELL_META$cell.xlim[2], n) + convert_x(1, "mm"), # x coordinate
        (1:n), # y coordinate
        cn, # label
        cex = 0.5, adj = c(0, 1), facing = "inside"
      )
    }
  }, bg.border = NA)


  lgd <- ComplexHeatmap::Legend(title = "Pearson's correlation", col_fun = col_fun1, at = c(-1, 0, 1))
  lgd2 <- ComplexHeatmap::Legend(title = "BH-corrected q-values", col_fun = col_fun2, at = c(0, 0.05, 1))
  lgd_list_vertical2 <- ComplexHeatmap::packLegend(lgd, lgd2)

  ComplexHeatmap::draw(lgd_list_vertical2, x = unit(50, "mm"), y = unit(1, "mm"), just = c("right", "bottom"))
}
