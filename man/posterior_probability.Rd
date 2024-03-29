% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/posterior_probability.R
\name{posterior_probability}
\alias{posterior_probability}
\title{Retrieves posterior probabilities from iCluster integration}
\usage{
posterior_probability(
  integrated_object,
  multimodal,
  metadata,
  threshold = 0.9,
  sense_check_variable = NULL
)
}
\arguments{
\item{integrated_object}{metadata in multiassay@metadata$integration$iCluster}

\item{multimodal}{metadata in multiassay@metadata$mulitmodal_object$omics}

\item{metadata}{metadata in multiassay@metadata$mulitmodal_object$metadata}

\item{threshold}{Absolute threshold to filter weights}

\item{sense_check_variable}{default to \code{NULL}. High weights should
coincide with stronger correlation if the sense_check_variable is an
important driver of variation in the designated factor. Will be used to
generate a plot if set to a covariate.}
}
\value{
List with positive and negative weights above and below the
specified threshold, distirbution plots of the weights, and dataframes
of the weights for each omic layer.
}
\description{
Retrieves posterior probabilities from iCluster integration
}
\seealso{
Other Multi-omic integration downstream analysis: 
\code{\link{OpenTarget_dataframe}()},
\code{\link{Transcription_Factor_enrichment_supervised}()},
\code{\link{Transcription_Factor_enrichment}()},
\code{\link{circos_TF}()},
\code{\link{clustering_DE_analysis}()},
\code{\link{extract_loadings}()},
\code{\link{extract_multiomic_signature}()},
\code{\link{extract_weigths}()},
\code{\link{integrative_results_clustering}()},
\code{\link{integrative_results_sparse_supervised}()},
\code{\link{integrative_results_supervised}()},
\code{\link{integrative_results_unsupervised}()},
\code{\link{multiomics_modules}()},
\code{\link{multiomics_network_matrix}()},
\code{\link{multiomics_network}()},
\code{\link{plot_module_trajectory}()},
\code{\link{pseudotime_inference}()},
\code{\link{sign_signature}()}
}
\concept{Multi-omic integration downstream analysis}
