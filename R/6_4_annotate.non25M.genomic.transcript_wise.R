#' Map genomic alignment coordinated to transcript coordinates
#'
#' Mappings in the genomic alignment \code{data.frame} are devided into two parts:
#' those that have \code{CIGAR == "25M"} and those that have  \code{CIGAR != "25M"}.
#' The first ones are mapped to exon coordinated from the Reference Annotation,
#' and the second ones are suspects to be exon-exon border mappings,
#' which is why they are mapped to transcript coordinates.
#'
#' @title Map non-25M genomic mapppings to transcripts coordinates
#' @param Alignment.genomic.non25M \code{data.frame} of raw genomic alignment, subsetted to non-25M mappings (\code{CIGAR != "25M"})
#' @param Annotation \code{data.frame}, reference genome annotation, must contatin \code{feature} column, and feature names must contain \code{"transcript"} or \code{"RNA"}
#' @return returns a \code{data.frame}
#' @seealso \code{annotate.G.alignment}
#' @author Vladislava Milchevskaya \email{milchv@gmail.com}
annotate.non25M.genomic.transcrits_wise <- 
  function( Alignment.genomic.non25M, Annotation)
{
  Annotation = Annotation[which(grepl("transcript|RNA",Annotation$feature)), ,
                          drop = FALSE]
  if (nrow(Annotation) == 0)
  {
    stop(message =
           "Error: Input full or transcriptomic annotation data.frame! 
         Check the column names (must be feature) 
         and feature names (must contain transcript or RNA) \n")
  }
  Annotation = Annotation[!duplicated(Annotation), ]
  stopifnot(all(Alignment.genomic.non25M$CIGAR != "25M"))

  Alignment.genomic.non25M = 
    Alignment.genomic.non25M[!duplicated(Alignment.genomic.non25M), ]
  AGquery = 
    cbind(Alignment.genomic.non25M, 
          sapply(Alignment.genomic.non25M$CIGAR, FUN = cigar2length) + 
            Alignment.genomic.non25M$POS - 1)

  colnames(AGquery)[ncol(AGquery)] <- "endPOS"
  colnames(AGquery)[   c( which(colnames(AGquery) == "RNAME"),
                          which(colnames(AGquery) == "POS"),
                          which(colnames(AGquery) == "endPOS")
  )]                      <-  c("chr", "start.q", "end.q")

  referenceDF_input = Annotation[ ,c("id", "seqid", "feature", 
                                     "start", "end", "strand", 
                                     "transcript_id", "gene_id")]
  AGinclusion.transcript = 
    findInclusion( AGquery, referenceDF  = referenceDF_input)
  return(AGinclusion.transcript)
}
