#' Generates a plot relating modules eigenvalue along pseudotime
#'
#' @param metadata
#' @param modules
#' @param variable
#' @param pseudotime_var
#'
#' @return
#' @export
#'
#' @examples
pseudotime_ordering_modules <- function(metadata,
                                        modules,
                                        variable,
                                        pseudotime_var = "pseudotime") {
  range01 <- function(x) {
    (x - min(x)) / (max(x) - min(x))
  }
  myPalette <- colorRampPalette(rev(RColorBrewer::brewer.pal(11, "Spectral")))

  pseudotime_ordering_modules <- list()

  for (i in modules) {
    metadata[, i] <- range01(metadata[, i])
    plot <- ggplot2::ggplot(metadata, aes_string(pseudotime_var, paste(i))) +
      ggplot2::geom_point(aes_string(colour = variable)) +
      ggplot2::geom_smooth(se = F) +
      ylim(0, 1) +
      theme_classic() +
      ggplot2::ylab(paste("Scaled changes (EigenValue)", i)) +
      viridis::scale_color_viridis()

    pseudotime_ordering_modules[[i]] <- plot
  }

  names(pseudotime_ordering_modules) <- modules

  return(pseudotime_ordering_modules)
}
