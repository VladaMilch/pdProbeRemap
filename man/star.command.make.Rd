% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/10_2_star_command.R
\name{star.command.make}
\alias{star.command.make}
\title{generate command for STAR aligner}
\usage{
star.command.make(path_to_STAR, sequence_fasta, genomeDir, organism, outfile,
  runThreadN = 6)
}
\arguments{
\item{path_to_STAR}{path to STAR aligner}

\item{sequence_fasta}{path to FASTA file with probe sequences (most likely, you will have to run alignment on a cluster, and copy the file from you local machine to the cluster)}

\item{genomeDir}{path to reference genome annotation, build for STAR}

\item{organism}{character, "Human", "Mouse" or "Dmel". Note! STAR parameters differ for Drosophila!}

\item{outfile}{outFileNamePrefix, add a dot in the end.}

\item{runThreadN}{outBAMsortingThreadN in STAR, 6 by default}
}
\value{
returns a character, STAR command
}
\description{
make command for STAR aligner
}
\details{
Provides a set of STAR parameters to make alignment for the probes, and returns a bash command as a character string. 
Needed mostly because alignments are run on clusters, but not local machines.
}
\examples{
path_to_STAR = "/path/to/STAR"
sequence_fasta = "/path/to/sequence.fasta"
genomeDir = "/path/to/reference/genome/"
organism = "Human"
outfile = "output_file."
star.command.make(path_to_STAR, sequence_fasta, genomeDir, organism, outfile)
}
\author{
Vladislava Milchevskaya \email{milchv@gmail.com}
}
