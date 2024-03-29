% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/integration_visualisations.R
\name{plot_optimal_cluster}
\alias{plot_optimal_cluster}
\title{plot_optimal_cluster}
\usage{
plot_optimal_cluster(optk1, optk2, try.N.clust)
}
\arguments{
\item{optk1}{Cluster Prediction Index}

\item{optk2}{Gap statistics}

\item{try.N.clust}{A integer vector to indicate possible choices of number of clusters.}
}
\value{
plot
}
\description{
plot optimal number of clusters from \code{getClustNum()}
Function adapted from package \code{MOVICS}  (see ref)
}
\references{
Lu, X., et al. (2020). MOVICS: an R package for multi-omics integration and visualization in cancer subtyping. Bioinformatics, 36(22-23), 5539–5541.
}
\seealso{
Other Plotting: 
\code{\link{circular_corrplot}()},
\code{\link{communities_network}()},
\code{\link{community_graph}()},
\code{\link{correlation_heatmap_clusters}()},
\code{\link{correlation_heatmap_supervised}()},
\code{\link{correlation_heatmap}()},
\code{\link{getStdiz}()},
\code{\link{interactive_network}()},
\code{\link{multiomics_heatmap}()},
\code{\link{plot_OpenTarget}()},
\code{\link{plot_communities}()},
\code{\link{plot_components}()},
\code{\link{plot_loadings}()},
\code{\link{visualise_pathway}()},
\code{\link{volcano_interactive_comparison}()},
\code{\link{volcano_interactive}()},
\code{\link{volcano_plot_deseq}()},
\code{\link{volcano_plot_limma}()}
}
\concept{Plotting}
