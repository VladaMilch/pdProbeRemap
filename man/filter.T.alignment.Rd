% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/4_filter.T.alignment.R
\name{filter.T.alignment}
\alias{filter.T.alignment}
\title{filter out non-unique transcriptomic mappings}
\usage{
filter.T.alignment(AT, level)
}
\arguments{
\item{AT}{annoated transcriptomic alignment, \code{"data.frame"}}

\item{level}{character string, takes values \code{"gene"} or \code{"transcript"}}
}
\value{
a \code{"data.frame"}
}
\description{
filter out non-unique mappings in Transcriptomic alignment
}
\details{
Annotated transcriptomic alignment \code{data.frame} processed,
probe IDs with more than one gene or transcript (see \code{"level"}) in the data.frame called
non-unique mappings and filtered out, all rows with these probes excluded.
}
\author{
Vladislava Milchevskaya \email{milchv@gmail.com}
}
