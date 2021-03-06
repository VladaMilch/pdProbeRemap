% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/3_annotate.T.alignment.R
\name{annotate.T.alignment}
\alias{annotate.T.alignment}
\title{annotate Transcriptomic Alignment data.frame}
\usage{
annotate.T.alignment(TAlignment, Annotation)
}
\arguments{
\item{TAlignment}{transcriptomic alignment \code{data.frame}}

\item{Annotation}{data.frame with regerence genome annotation}
}
\value{
a \code{"data.frame"}.
}
\description{
maps transcript IDs to gene IDs, strand, ...
}
\details{
This function maps transcript IDs from the Transcriptomic alignment \code{data.frame}
to the ones from the Reference Annotation \code{data.frame}, and extracts gene ids.
}
\seealso{
\code{filter.T.alignment}, \code{read.gtf? from another package -- how to parce reference genome annotations}
}
\author{
Vladislava Milchevskaya \email{milchv@gmail.com}
}
