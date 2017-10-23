#' finds exon-exon border probes and corresponding transcripts
#'
#' Identifies which probes fall into a transcript region in the reference genome annotation,
#' and map simultaneously to the left end of one exon and to the left end of another exon of the same transcript,
#' based on their CIGAR string and coordinates in the raw genomic alignment.
#' Only mapping that have \code{CIGAR != "25M"} are considered.
#'
#' @title What function does (short)
#' @param Als object of the class \code{Alignment}
#' @param Annotation \code{data.frame}, reference genome annotation
#' @param level character string, equals \code{"gene"} or \code{"transcript"}
#' @return a \code{"output type"}.
#' @seealso \code{function or class name}
#' @author Vladislava Milchevskaya \email{milchv@gmail.com}
getExonExonBorder <- function(Als, Annotation, level)
{
  Annotation = subset(Annotation, strand %in% c("+", "-"))
  stopifnot(nrow(Annotation) > 0)

  stopifnot(level %in% c("gene", "transcript"))

  options(stringsAsFactors = FALSE)
  AGraw.non25M = Als@samGenome$x[which(Als@samGenome$x$CIGAR != "25M"), c("QNAME", "CIGAR", "POS", "RNAME")]
  if (nrow(AGraw.non25M) == 0) {return( data.frame() )}
  
  AGraw.non25M.annotated = annotate.non25M.genomic.transcrits_wise(AGraw.non25M, Annotation = Annotation) # function 6_4 # to take only probes that fell in the transcr. region
  aa = length(unique(AGraw.non25M.annotated$QNAME)); bb = length(unique(AGraw.non25M$QNAME))
  message( paste0(aa, " probes out of all ", bb, " non-25M probes fell into transcripts regions..") )

  m_left = sapply(AGraw.non25M.annotated$CIGAR, cigar2_left)
  m_right = sapply(AGraw.non25M.annotated$CIGAR, cigar2_right)
  summa = m_left + m_right

  AGraw.non25M.annotated_RL = cbind(AGraw.non25M.annotated, m_left = sapply(AGraw.non25M.annotated$CIGAR, cigar2_left), m_right = sapply(AGraw.non25M.annotated$CIGAR, cigar2_right), summa)
  AGraw.non25M.annotated_RL_1 = subset(AGraw.non25M.annotated_RL, summa == 25)
  AGraw.non25M.annotated_RL_1 = AGraw.non25M.annotated_RL_1[!duplicated(AGraw.non25M.annotated_RL_1), ]

  stopifnot(all(AGraw.non25M.annotated_RL_1$m_left > 0))
  stopifnot(all(AGraw.non25M.annotated_RL_1$m_right > 0))

  end.left.M = AGraw.non25M.annotated_RL_1$start.q + AGraw.non25M.annotated_RL_1$m_left - 1
  start.right.M = AGraw.non25M.annotated_RL_1$end.q - AGraw.non25M.annotated_RL_1$m_right + 1

  AGraw.non25M.annotated_RL_1_1 = as.data.frame(cbind(AGraw.non25M.annotated_RL_1, end.left.M, start.right.M))
  ANN_exon = subset(Annotation, feature == "exon")
  stopifnot(nrow(ANN_exon) > 0)

  AGquery_left = AGraw.non25M.annotated_RL_1_1[ , c("QNAME",  "chr","start.q", "end.left.M")]
  colnames(AGquery_left)[ncol(AGquery_left)] <- "end.q"
  AGinclusion.LEFT = findInclusion( AGquery_left, referenceDF  = ANN_exon)[ ,c("QNAME", "start.q", "end.q", "id", "feature",  "transcript_id" ,"start", "end", "gene_id", "chr")]

  AGquery_right = AGraw.non25M.annotated_RL_1_1[ , c("QNAME",  "chr","end.q", "start.right.M")]
  colnames(AGquery_right)[ncol(AGquery_right)] <- "start.q"
  AGinclusion.RIGHT = findInclusion( AGquery_right, referenceDF  = ANN_exon)[ ,c("QNAME", "start.q", "end.q", "id", "feature",  "transcript_id" ,"start", "end", "gene_id", "chr")]

  AG_bothSides = merge(AGinclusion.LEFT[ , c("QNAME", "transcript_id", "chr", "start.q", "end.q", "start", "end")],
                       AGinclusion.RIGHT[ , c("QNAME", "transcript_id", "gene_id", "start.q", "end.q", "start", "end")],
                       by.x = c("QNAME", "transcript_id"),
                       by.y = c("QNAME", "transcript_id"),
                       suffixes = c(".left", ".right"))
  if (any(duplicated(AG_bothSides)))
  {
    AG_bothSides = AG_bothSides[duplicated(AG_bothSides), ]
  }
  AG_bothSides1 = AG_bothSides[order(as.numeric(AG_bothSides$QNAME)), ]
  AG_bothSides2 = subset(AG_bothSides1, start.q.right == start.right)
  AG_bothSides3 = subset(AG_bothSides2, end.q.left == end.left)

  ExonExonBorderProbes = AG_bothSides3[ ,c("QNAME", "transcript_id","gene_id", "chr") ]
  ExonExonBorderProbes = ExonExonBorderProbes[!duplicated(ExonExonBorderProbes), ]
  colnames(ExonExonBorderProbes)[which(colnames(ExonExonBorderProbes) == "QNAME")]  <- "probe_id"

  if (level == "gene")
  {
    ExonExonBorderProbes_G = ExonExonBorderProbes[ ,c("probe_id","gene_id", "chr") ]
    ExonExonBorderProbes_G = ExonExonBorderProbes_G[!duplicated(ExonExonBorderProbes_G), ]
    return(ExonExonBorderProbes_G)
  }
  if (level == "transcript")
  {
    return(ExonExonBorderProbes)
  }
}

