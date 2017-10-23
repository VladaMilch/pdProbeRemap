
#' Short description string here
#'
#' Long description string
#'
#' @title What function does (short)
#' @param Alignments object of class \code{Alignments}, 
#' prodused by \code{processAlignments} function
#' @param Annotation \code{data.frame} reference genome annotation
#' @param level \code{character} string, equals  \code{"gene"} 
#' or \code{"transcript"}
#' @return returns \code{data.frame}.
#' @seealso \code{function or class name}
#' @author Vladislava Milchevskaya \email{milchv@gmail.com}

alignments2probesetsRaw  <- function(Alignments, 
                                     Annotation, 
                                     level = "Input_gene_or_transcript")
{
  stopifnot(level %in% c("gene", "transcript"))
  ##### check annotation #####
  Annotation = subset(Annotation, strand %in% c("+", "-"))
  stopifnot(nrow(Annotation) > 0)

  ##### annotate T alignment ##### (needed for the gene level, and ?)
  AT = annotate.T.alignment(
    TAlignment = Alignments@samTranscriptome$x, 
    Annotation = Annotation)
  stopifnot(nrow(AT) > 0)
  # pseudogenes in warning - shall i do something with it? #

  ##### exclude probes with non-unique mappings on the corresponding level ###
  # filter AT, subset AG for the good AT probes #
  AT.filtered = filter.T.alignment(AT, level = level) 
  #first AT, so that to subset AG properly based on AT
  stopifnot(nrow(AT.filtered) > 0)

  good_transcriptomic_probes = unique(AT.filtered$probe_id) 
  #good transcriptomic probes
  
  bad_transcriptomic_probes = setdiff(AT$probe_id, AT.filtered$probe_id)
  message( paste0("Good probes in transcriptomic alignment: ",  
                  length(good_transcriptomic_probes), 
                  ", bad probes: ", 
                  length(bad_transcriptomic_probes)))

  Alignment.genome.subsetted = 
    subset(Alignments@samGenome$x, 
           QNAME %in% good_transcriptomic_probes) 
  # only good transcriptomic probes
  
  Alignment.genome.subsetted = 
    data.frame(cbind(Alignment.genome.subsetted, 
                     instanceID = c(1:nrow(Alignment.genome.subsetted)))) 
  # instance id to track mappings, not probes

  AG.25M = subset(Alignment.genome.subsetted, CIGAR == "25M")
  AG.ne25M = subset(Alignment.genome.subsetted, CIGAR != "25M")
  AG.25M.annotated = annotate.G.alignment(AG.25M, Annotation = Annotation)
  AG.exonexonborder.annotated = 
    getExonExonBorder(Als = Alignments, 
                      Annotation = Annotation, 
                      level = "transcript") 
  # to make sure transcript ids are kept 
  # (gene-probe correspondance will remain the same)
  
  
  instances.25M.annotated = unique(AG.25M.annotated$instanceID)
  instances.25M.NONannotated = setdiff(AG.25M$instanceID, 
                                       AG.25M.annotated$instanceID)
  #instances.non25M = AG.ne25M$instanceID

  probes.A25 = unique(AG.25M.annotated$probe_id)
  probes.neAnnot_25 = 
    unique(subset(Alignment.genome.subsetted, 
                  instanceID %in%  instances.25M.NONannotated)$QNAME)
  #probes.ne25 = unique(subset(Alignment.genome.subsetted, 
  #               instanceID %in%  instances.non25M)$QNAME) 
  # those that ARE in good transcr ptobes

  #good.G.probes = setdiff( union(probes.A25, probes.ne25), probes.neAnnot_25)
  #stopifnot(length(union( union(probes.A25, probes.ne25), 
  # probes.neAnnot_25)) == length(unique(Alignment.genome.subsetted$QNAME)))
  bad_genomic_probes = probes.neAnnot_25 # probes that have non-annotated 25M

  ATT = subset(AT.filtered, !(probe_id %in% bad_genomic_probes))
  AGG = subset(AG.25M.annotated, !(probe_id %in% bad_genomic_probes))
  ATG = merge(ATT[ ,c("probe_id", "transcript_id", "gene_id", "strand")], 
              AGG[ ,c("probe_id", "transcript_id", "gene_id", "strand")])
  if(nrow(AG.exonexonborder.annotated) > 0)
  {
     AGG.exex = subset(AG.exonexonborder.annotated, 
                       !(probe_id %in% 
                           union(bad_genomic_probes, 
                                 bad_transcriptomic_probes)))
     ATTG = rbind(ATG[ , c("probe_id", "transcript_id", "gene_id")],
                  AGG.exex[ ,c("probe_id", "transcript_id", "gene_id")])
  }else
  {
    ATTG = ATG[ , c("probe_id", "transcript_id", "gene_id")]
  }
  #probes_not_in_merged = setdiff(ATT$probe_id, ATG$probe_id) 
  # these are good probes (and mostly exon-exon border probes) 
  # AND probes fallen in the transcript, but may be not all in exons 
  # (these need to be filtered out)

  ATTG = ATTG[!duplicated(ATTG), ]

  ATTG.filtered = filter.ATTG.alignment.df(ATTG, level = level)
  return(ATTG.filtered)
}

# internal function #
filter.ATTG.alignment.df <- function(AG, level = level) 
  #non-unique mappings filtered out
{
  if(level == "gene")
  {aaa1 = AG[,c("probe_id", "gene_id")]}

  if(level == "transcript")
  {aaa1 = AG[,c("probe_id", "transcript_id")]}
  aaa2 = aaa1[!duplicated(aaa1),]

  ad2 = aaa2[duplicated(aaa2[,c("probe_id")]),"probe_id"]
  non_unique = (unique(ad2))
  AG.filtered = subset(AG, !(probe_id %in% non_unique))
  return(AG.filtered)
}

