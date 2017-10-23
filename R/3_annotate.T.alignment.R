#' maps transcript IDs to gene IDs, strand, ...
#'
#' This function maps transcript IDs from the Transcriptomic alignment \code{data.frame}
#' to the ones from the Reference Annotation \code{data.frame}, and extracts gene ids.
#'
#' @title annotate Transcriptomic Alignment data.frame
#' @param TAlignment transcriptomic alignment \code{data.frame}
#' @param Annotation data.frame with regerence genome annotation
#' @return a \code{"data.frame"}.
#' @seealso \code{filter.T.alignment}, \code{read.gtf? from another package -- how to parce reference genome annotations}
#' @author Vladislava Milchevskaya \email{milchv@gmail.com}
###################################
#
###################################
annotate.T.alignment <- function(TAlignment, Annotation)
{
  #logger <- create.logger()
  #logfile(logger) <- file.path('alignment_processing.log')
  #level(logger) <- "INFO"

  Annotation = Annotation[which(grepl("transcript|RNA",Annotation$feature)), ,
                          drop = FALSE]
  if (nrow(Annotation) == 0)
  {
    stop(message =
            "Error: Input full or transcriptomic annotation data.frame! 
         Check the column names (must be feature) 
         and feature names (must contain transcript or RNA) \n")
  }

  if (  !all (TAlignment$RNAME %in% Annotation$transcript_id)  ) # if same annotation used for alignment (with STAR) and for processing now
  {
    #warn(logger, message = "Not all transcript IDs from the Transcriptomic ALignment appear in the Annotation... \n")
    warning( 
      message = "Not all transcript IDs from the Transcriptomic Alignment 
      appear in the Annotation... for DROSOPHILA: pseudogenes? \n")
  }

  prev_nrow = nrow(TAlignment)
  TAlignment = subset(TAlignment, CIGAR == "25M")
  curr_nrow = nrow(TAlignment)

  if(curr_nrow / prev_nrow < 0.9)
  {
    #warn(logger, message = paste0("N of mappings in raw Alignment.transcriptome: ", prev_nrow, "   after filtering for 25M cigar:  ", curr_nrow ))
    warning( message = "Ideally, almost all of the probe sequences, 
             it they align to the Transcriptome, 
             should have CIGAR string 25M. 
             \nIf not, check alignment files and bam-to-sam convertion!\n")
  }

  AT = merge(TAlignment, Annotation[,c("gene_id", 
                                       "transcript_id", 
                                       "strand", 
                                       "feature")], 
             by.x = "RNAME", 
             by.y = "transcript_id" )
  colnames(AT)[1:2] <- c("transcript_id", "probe_id")
  ATT = AT[,-which(names(AT) %in% c("RNEXT","PNEXT", "TLEN")), 
           drop = FALSE] # relevant only for paired-end alignment (not the case with probes)
  return(ATT)
}
