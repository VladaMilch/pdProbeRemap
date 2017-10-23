#' Reads probe alignment files and creates an \code{Alignments} object
#'
#' This funtion reads genomic alignment file (\code{Aligned.out.sam} postfix 
#' required)
#' and transcriptomic alignment file (\code{Aligned.toTranscriptome.out.bam} 
#' postfix requires)
#' and generated an object of class \code{"Alignments"}
#'
#' @title Read probe alignments from a directory
#' @param alignment.dir the path to a directory that contains genomic 
#' and transcriptomic alignments
#' @return a \code{"Alignments"}.
#' @author Vladislava Milchevskaya \email{milchv@gmail.com}
#' @examples 
#' sam_file <- 
#'    dir(system.file("extdata",package="pdProbeRemap"), 
#'    pattern="example_drosophila.Aligned.out.sam",full.names=TRUE)
#' alignment_dir <- 
#'    strsplit(sam_file,  split = "example_drosophila.Aligned.out.sam")[[1]]
#' Als <- processAlignments(alignment_dir)
#' 
#' @include 0_DataClasses.R
#' @export
#' @importFrom reshape2 melt
#' @importFrom Rsamtools asSam
### 2nd function ####
processAlignments <- function(alignment.dir = NA)
{
  #if (is.na(alignment.dir) & (!(is.na(outputDir))))
  #{alignment.dir = file.path(outputDir, "Alignments")}
  alignment.to.genome.path = list.files(path = alignment.dir,
                                        pattern = "Aligned.out.sam",
                                        ignore.case = TRUE,
                                        full.names = TRUE)
  alignment.to.transcriptome.path = 
    list.files(path = alignment.dir,
               pattern = "Aligned.toTranscriptome.out.bam",
               ignore.case = TRUE,
               full.names = TRUE)
  alignment.to.transcriptome.name = 
    list.files(path = alignment.dir,
               pattern = "Aligned.toTranscriptome.out.bam",
               ignore.case = TRUE,
               full.names = FALSE)

  if (!all(length(alignment.to.genome.path) == 1 ,
           length(alignment.to.transcriptome.path) == 1,
           length(alignment.to.transcriptome.name) == 1,
           length(alignment.dir) == 1))
  {
    message = 
      paste('Wrong alignment directory path or not all alignment files found.')
    stop(message = message)
  }

  alignment.to.transcriptome.name.woBAM = 
    unlist(strsplit(alignment.to.transcriptome.name, split = ".bam"))

  Als <- new("Alignments")
  Als@path_AlignmentToGenome_Sam <- alignment.to.genome.path
  Als@path_AlignmentToTranscriptome_Bam <- alignment.to.transcriptome.path
  Als@path_AlignmentsDir <- alignment.dir

  samGenome <- read.sam(alignment.to.genome.path)
  samFileTranscriptome <- 
    Rsamtools::asSam(
      alignment.to.transcriptome.path, 
      destination = file.path(alignment.dir, 
                              alignment.to.transcriptome.name.woBAM), 
      overwrite = TRUE)
  samTranscriptome <- read.sam(samFileTranscriptome)

  # check that the alignments are OK #
  g25 = reshape2::melt(table(samGenome$x$CIGAR == "25M"))
  if ( all(c("TRUE", "FALSE") %in% g25$Var1))
  {
    if (g25[which(g25$Var1 == TRUE),"value"] <= 10*g25[which(g25$Var1 == FALSE),
                                                       "value"])
    {
      warning("Most of the probe mappings in the Genomic Alignment 
              should have CIGAR = 25M. 
              Now >10% of the mappings have different CIGAR string! \n")}
  }
  t25 = reshape2::melt(table(samTranscriptome$x$CIGAR == "25M"))
  if (nrow(t25) > 1)
  {
    if(t25[which(g25$Var1 == TRUE),"value"] <= 10*t25[which(g25$Var1 == FALSE),
                                                      "value"])
    {
      warning("Most of the probe mappings in the 
              Genomic Alignment should have CIGAR = 25M. 
              Now >10% of the mappings have different CIGAR string! \n")
    }
  }

  if (nrow(t25) == 1 & t25[1,1] == FALSE)
  {stop("Not even a single mapping in Transcriptome ALignment 
        with CIGAR = 25M. Corrupted input file? \n")}

  Als@samGenome <- samGenome
  Als@samTranscriptome <- samTranscriptome
  return(Als)
}

