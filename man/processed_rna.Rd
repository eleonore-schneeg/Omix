% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/custom_rna.R
\name{processed_rna}
\alias{processed_rna}
\title{Adds a custom processed transcriptomics slot to multiassay object}
\usage{
processed_rna(multiassay, custom_processed_df)
}
\arguments{
\item{multiassay}{Multiassay experiment object generated by Omix}

\item{custom_processed_df}{custom processed data frame}
}
\value{
Adds custom processed dataset to Multiassay experiment
object generated by Omix
}
\description{
Adds a custom processed transcriptomics slot to multiassay object
}
\seealso{
Other Helper: 
\code{\link{counts_rna}()},
\code{\link{format_res_deseq}()},
\code{\link{format_res_limma}()},
\code{\link{get_background}()},
\code{\link{get_metadata}()},
\code{\link{processed_proteomics}()}
}
\concept{Helper}
