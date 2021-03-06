% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/0_DataClasses.R
\docType{class}
\name{Alignments}
\alias{Alignments}
\alias{Alignments-class}
\title{The Alignments class}
\description{
This class contains an example. This line goes into the description
}
\details{
This line and the next ones go into the details.
This line thus appears in the details as well.
}
\note{
You can still add notes
}
\section{Slots}{

 \describe{
   \item{\code{path_AlignmentToGenome_Sam}:}{character indicating path to the genomic alignment, \code{".sam"}}
   \item{\code{path_AlignmentToTranscriptome_Bam}:}{character indicating path to the transcriptomic alignment, \code{".bam"}}
   \item{\code{path_AlignmentsDir}:}{character indicating path to the directory with both alignments, \code{".bam"}}
   \item{\code{samGenome}:}{list of 2, contains header list and the genomic alignment \code{"data.frame"}}
   \item{\code{samTranscriptome}:}{list of 2, contains header list and the transcriptomic alignment \code{"data.frame"}}
 }
}

\author{
Vladislava Milchevskaya
}
