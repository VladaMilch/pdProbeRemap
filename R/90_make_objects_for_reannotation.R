#' convert ProbeSets to a format needed for package building
#'
#' Converts ProbeSets \code{data.frame} into a list of two \code{data.frame}s:
#' \code{ProbeSetDataAnnotation} and \code{ProbeSetDataFrame}.
#' The \code{ProbeSetDataAnnotation} contains 4 colunms:
#' \code{fsetname}, 
#' \code{fsetid}, 
#' \code{strand}, 
#' \code{chrom}, where \code{fsetname} is an Ensembl gene/transcript ID.
#' The \code{ProbeSetDataFrame} contains two colunms: 
#' \code{fsetid} and \code{fid},
#' where \code{fid} is a probe ID from the original array annotation.
#'
#' @title Converts ProbeSets to package building input
#' @param Probesets \code{data.frame}, must contain columns 
#' \code{"transcript_id"}, \code{"gene_id"} and \code{"probe_id"}
#' @param level \code{character}, must be "gene" or "transcript"
#' @param Annotation reference genome annotation \code{data.frame}
#' @param shift numeric, 
#' a constant that is being added to all main probe set ids, 
#' in order to make them larger than control probe set ids 
#' (and thus avoid same probe set ids when unwanted)
#' @return a \code{list} with two \code{data.frame}s
#' @author Vladislava Milchevskaya \email{milchv@gmail.com}

# generating objects to feed to probe set making scripts #
# probe sets raw --> probe sets + feature annotation
make_objects_for_reannotation <-
  function(
    Probesets, #gene set object
    level = "gene",
    Annotation, 
    shift)
  {
    stopifnot(class(Probesets) == "data.frame")
    stopifnot(class(Annotation) == "data.frame")
    stopifnot(level %in% c("gene", "transcript"))
    stopifnot("feature" %in% colnames(Annotation))
    stopifnot(c("probe_id", "transcript_id", "gene_id") %in%  
                colnames(Probesets))
    stopifnot(c("transcript_id", "strand", "seqid", "id")  %in% 
                colnames(Annotation))
    stopifnot(class(shift) == "numeric")
    
    Annotation = Annotation[grep("transcript|RNA|rna|Transcript|Rna", 
                                 Annotation$feature), ]
    Annotation = subset(Annotation, strand %in% c("+", "-"))

    if (level == "transcript")
    {
      AnnotationSUBS = subset(Annotation, transcript_id %in% 
                                as.vector(Probesets$transcript_id))
      AnnotationSUBS = AnnotationSUBS[,c("transcript_id", "strand", "seqid")]
      AnnotationSUBS = AnnotationSUBS[!duplicated(AnnotationSUBS),]
      rrr1 = merge(Probesets, AnnotationSUBS)
      rrr1 = rrr1[,c("transcript_id", "strand", "seqid", "probe_id")]
      rrr1<-subset(rrr1, !duplicated(rrr1))
      head(rrr1)

    #stopifnot(colnames(Probesets[[1]])[ncol(Probesets[[1]])] == 
      #"transcript_id")
      fsetname = unique(as.vector(Probesets$transcript_id)) # for transcript_id
      fsetid <-  
        as.numeric(Annotation[match(fsetname, 
                                    Annotation$transcript_id),"id"]) + shift
      fset.table <- as.data.frame(cbind(fsetid, fsetname))

      psann <- subset(rrr1[,c("transcript_id", "strand", "seqid")], 
                      !duplicated(rrr1[,c("transcript_id", 
                                          "strand", 
                                          "seqid")]))
      psann2 = merge(fset.table, psann, 
                     by.y = "transcript_id", 
                     by.x = "fsetname")
      fset.table$fsetid[match(rrr1$transcript_id, 
                              fset.table$fsetname)]
      p.df <- as.data.frame(
        cbind(rrr1$probe_id, 
              fset.table$fsetid[match(rrr1$transcript_id, 
                                      fset.table$fsetname)]))
      colnames(p.df) <- c("fid","fsetid")
      #
    }
    if (level == "gene")
    {
      AnnotationSUBS = subset(Annotation, gene_id %in% 
                                as.vector(Probesets$gene_id))
      AnnotationSUBS = AnnotationSUBS[,c("gene_id", "strand", "seqid")]
      AnnotationSUBS = AnnotationSUBS[!duplicated(AnnotationSUBS),]
      rrr1 = merge(Probesets, AnnotationSUBS)
      rrr1 = rrr1[,c("gene_id", "strand", "seqid", "probe_id")]
      rrr1<-subset(rrr1, !duplicated(rrr1))
      head(rrr1)

      #stopifnot(colnames(Probesets[[1]])[ncol(Probesets[[1]])] == "gene_id")
      fsetname = unique(as.vector(Probesets$gene_id)) # for gene_id
      fsetid <-  
        as.numeric(Annotation[match(fsetname, Annotation$gene_id),"id"]) + shift
      fset.table <- as.data.frame(cbind(fsetid, fsetname))

      psann <- subset(rrr1[,c("gene_id", "strand", "seqid")], 
                      !duplicated(rrr1[,c("gene_id", "strand", "seqid")]))
      psann2 = merge(fset.table, psann, 
                     by.y = "gene_id", 
                     by.x = "fsetname")
      fset.table$fsetid[match(rrr1$gene_id, 
                              fset.table$fsetname)]
      p.df <- as.data.frame(
        cbind(rrr1$probe_id, 
              fset.table$fsetid[match(rrr1$gene_id, 
                                      fset.table$fsetname)]))
      colnames(p.df) <- c("fid","fsetid")
    }
    colnames(psann2)[4] <- "chrom"
    #head(psann2)
    #head(p.df)

    PSDA <- psann2
    PSDF <- p.df
    res = list(PSDA, PSDF)
    
    # # check if all fine #
    # random_row_index <- sample(1:nrow(Probesets), 1)
    # random_fid <- Probesets[random_row_index, "probe_id"]
    # random_fsetid <- subset(PSDF, fid == random_fid)$fsetid
    # random_fsetname <- subset(PSDA, fsetid == random_fsetid)$fsetname
    # stopifnot(
    #   random_fsetname == Probesets[random_row_index, c("gene_id")] | 
    #     random_fsetname == Probesets[random_row_index, c("transcript_id")])
    # message("Alll fine!")
    # # all fine! #
    
    names(res) <- c( "ProbeSetDataAnnotation", "ProbeSetDataFrame")
    return(res)
  }

#### helper check functions ####

check_ProbeSetDataAnnotation <- function(ProbeSetDataAnnotation)
{
  check1 <- is.data.frame(ProbeSetDataAnnotation)
  if(!check1) stop("ERROR: the ProbeSetDataAnnotation must be a data.frame")
  
  check2 <- all(c("fsetid", "fsetname", "strand", "chrom") %in% 
                  colnames(ProbeSetDataAnnotation))
  if(!check2) 
    stop("ERROR: the ProbeSetDataAnnotation must contain these columns: 
         \n fsetid, fsetname, strand, chrom \n")
  
  #check3 <- 
  # all(apply(ProbeSetDataAnnotation[,c("fsetid", "strand", "chrom")], 
  # 2, is.numeric))
  #if(!check3) stop("ERROR: the fsetid,  
  # strand and chrom columns must be numeric \n")
  
}

check_ProbeSetDataFrame <- function( ProbeSetDataFrame )
{
  check1 = (is.data.frame(ProbeSetDataFrame) || is.matrix(ProbeSetDataFrame))
  if (!check1) 
    stop( "ERROR: the probe-to-probeset correspondence object 
          must be a matrix or a data.frame. \n" )
  
  check2 <- setequal(colnames(ProbeSetDataFrame), c("fid", "fsetid"))  
  if (!check2) 
    stop( "ERROR: wrong columns / column names if the ProbeSetDataFrame. \n" )
  
  check3 <- all(apply(ProbeSetDataFrame, 2, is.numeric))
  if (!check2) 
    stop( "ERROR: both fid and fsetid must be numeric. \n" )
}
