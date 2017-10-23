#' The Alignments class
#'
#' This class contains an example. This line goes into the description
#'
#' This line and the next ones go into the details.
#' This line thus appears in the details as well.
#'
#'@section Slots:
#'  \describe{
#'    \item{\code{path_AlignmentToGenome_Sam}:}{character indicating path to the genomic alignment, \code{".sam"}}
#'    \item{\code{path_AlignmentToTranscriptome_Bam}:}{character indicating path to the transcriptomic alignment, \code{".bam"}}
#'    \item{\code{path_AlignmentsDir}:}{character indicating path to the directory with both alignments, \code{".bam"}}
#'    \item{\code{samGenome}:}{list of 2, contains header list and the genomic alignment \code{"data.frame"}}
#'    \item{\code{samTranscriptome}:}{list of 2, contains header list and the transcriptomic alignment \code{"data.frame"}}
#'  }
#'
#' @note You can still add notes
#' @name Alignments
#' @rdname Alignments
#' @aliases Alignments-class
#' @exportClass Alignments
#' @author Vladislava Milchevskaya

setClass("Alignments",
         representation
         (
         path_AlignmentToGenome_Sam = "character",
         path_AlignmentToTranscriptome_Bam = "character",
         path_AlignmentsDir = "character",
         samGenome = "list",
         samTranscriptome = "list"
         )
)

setClass("AnnotatedAlignments",
         representation
         (
           Genomic = "data.frame",
           Transcriptomic = "data.frame",
           Annotation = "data.frame"
         )
)
