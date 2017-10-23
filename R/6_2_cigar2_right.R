#' Length of the last M-region
#'
#' If a probe maps exon-exon border region of a transcript,
#' it does not appear in the Genomic Alignment as a mapping with
#' CIGAR = "25M", but as mapping with xxM in the beginning of
#' the CIGAR string, and xxM in the end,
#' which represents length of the regions mapped to the left and right exons correspondingly.
#' The \code{cigar2_right} function gets length of the right exon mapping,
#' if CIGAR string has M of the right-most position.
#'
#' @title cigar - to - length of the right exon mapping
#' @param cigar_string CIGAR character string
#' @return returns \code{integer}.
#' @seealso \code{cigar2_left}, \code{cigar2length}
#' @author Vladislava Milchevskaya \email{milchv@gmail.com}

# right  exon #
cigar2_right <- function(cigar_string)
{
  M_values = as.numeric(unlist(stringr::str_extract_all(
    stringr::str_extract_all(cigar_string, "[0-9]+[M]$"), "[0-9]+")))
  return(M_values)
}
