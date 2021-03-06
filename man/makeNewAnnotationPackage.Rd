% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/94_makeNewAnnotationPackage.R
\name{makeNewAnnotationPackage}
\alias{makeNewAnnotationPackage}
\title{What function does (short)}
\usage{
makeNewAnnotationPackage(alignment.dir, Annotation, outputDir, level,
  min_probe_number, quiet = FALSE, pkgNameSUFFIX, fastaDir)
}
\arguments{
\item{alignment.dir}{path to the directory that contains alignment files}

\item{Annotation}{regerence genome annotation \code{data.frame}}

\item{outputDir}{output directory}

\item{level}{character string, can take values "gene" or "transcript", specifying if new probe sets should be formed per gene or per transcript}

\item{min_probe_number}{numeric, minimal number or probes in a probe set}

\item{quiet}{logical}

\item{pkgNameSUFFIX}{character, suffix to the package name, usually starts with dor: '.hereIsYourSuffix'}
}
\value{
a \code{"output type"}.
}
\description{
Short description string here
}
\details{
Long description string
}
\examples{
data("Alignments_class_example")
data("Annotation_example")
data("seed_example")
alignmentDir <- 
  unlist(strsplit(dir(system.file("extdata",package="pdProbeRemap"), 
                      pattern="example_drosophila.Aligned.out.sam",
                      full.names=TRUE), split = "example_drosophila.Aligned.out.sam"))[1]

outputDir_example = system.file("data", package = "pdProbeRemap")

\dontrun{
makeNewAnnotationPackage(alignment.dir = alignmentDir, 
                         Annotation = Annotation_example, 
                         outputDir = outputDir_example, 
                         level = "gene", 
                         min_probe_number = 1, 
                         pkgNameSUFFIX = ".example")
                         }
}
\seealso{
\code{function or class name}
}
\author{
Vladislava Milchevskaya \email{milchv@gmail.com}
}
