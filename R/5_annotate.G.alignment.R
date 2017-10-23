# !see documentation below! #
############################################################
########## genome alignment annotation #############
############################################################
# probes filtering: only CIGAR == 25M
# use exon level annotation

findInclusion <- function(queryDF, referenceDF)
{
  stopifnot(  c("chr","start.q","end.q") %in% colnames(queryDF))
  stopifnot(  c("seqid", "start", "end", "strand") %in% colnames(referenceDF))

  referenceDF = subset(referenceDF, seqid %in% intersect(queryDF$chr, 
                                                         referenceDF$seqid))
  queryDF = subset(queryDF, 
                   chr %in% intersect(queryDF$chr, referenceDF$seqid))

  gr0 = with(referenceDF, 
             GenomicRanges::GRanges(seqid, 
                                    IRanges::IRanges(start=start, end=end), 
                                    strand = strand))
  gr1 = with(queryDF, 
             GenomicRanges::GRanges(chr, 
                                    IRanges::IRanges(start=start.q, end=end.q)))
  hits = IRanges::findOverlaps(gr1, gr0)
  resultOverlap = cbind(queryDF[hits@from,], #queryHits
                        referenceDF[hits@to,]) #subjectHits
  resultInclusion = 
    resultOverlap[which(resultOverlap$start.q >= resultOverlap$start & 
                          resultOverlap$end.q <= resultOverlap$end),]
  return(resultInclusion)
}

annotate.G.alignment_exon <- function(GAlignment, Annotation)
{
  starttime = Sys.time()
  if (length(which(names(GAlignment) %in% c("RNEXT","PNEXT", "TLEN"))!=0))
  {GAlignment = 
    GAlignment[,-which(names(GAlignment) %in% c("RNEXT","PNEXT", "TLEN"))]} # relevant only for paired-end alignment (not the case with probes)
  Annotation = Annotation[which(grepl("exon", Annotation$feature)),] # dont restrict to mRNA! to be able to do a more stringent filtering later, and forbid alignment to non-coding RNAs
  if (nrow(Annotation) == 0) 
  {
    stop(
      message = "Error: Input full or annotation data.frame! 
      (needs 'exon' present in features column)")
  }

  message(paste0("\nNumber of probes in the alignment before processing: ", 
                 length(unique(GAlignment$QNAME))))

  GAlignment25 = subset(GAlignment, GAlignment$CIGAR == "25M")
  message(paste0("\nNumber of 25M probes: ", 
                 length(unique(GAlignment25$QNAME))))

  AGquery = cbind(GAlignment25, GAlignment25$POS + 25 - 1)
  colnames(AGquery)[ncol(AGquery)] <- "endPOS"
  colnames(AGquery)[   c( which(colnames(AGquery) == "RNAME"),
                          which(colnames(AGquery) == "POS"),
                          which(colnames(AGquery) == "endPOS")
  )]                      <-  c("chr", "start.q", "end.q")

  AGinclusion.transcript = findInclusion( AGquery, Annotation)
  colnames(AGinclusion.transcript)[
    which(colnames(AGinclusion.transcript) == "QNAME")] <- "probe_id"
  endtime = Sys.time()
  simpleMessage(
    paste0("Calculation time for the genome alignment annotation: ", 
           endtime-starttime, "\n"))
  return(AGinclusion.transcript)
}


#' annotate probe genomic alignment with \code{CIGAR = 25M} only
#'
#' Raw genomic alignment \code{data.frame} subsetted for mappings with \code{CIGAR = 25M} only,
#' these are mapped to the genomic coordinates of exons from the Reference Annotation \code{data.frame}.
#' A probes that falls in the region of annotated exon, is mapped to all transcripts corresponding to the exon.
#'
#' @title annotates probes genomic alignment
#' @param GAlignment \code{data.frame} of the raw genomic alignment (output of \code{processAlignments} function)
#' @param Annotation reference genome annotation,  \code{data.frame}
#' @return a \code{"data.frame"}.
#' @seealso \code{processAlignments}
#' @author Vladislava Milchevskaya \email{milchv@gmail.com}

annotate.G.alignment = annotate.G.alignment_exon
