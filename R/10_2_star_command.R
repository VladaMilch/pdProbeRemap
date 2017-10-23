
#' make command for STAR aligner
#'
#' Provides a set of STAR parameters to make alignment for the probes, and returns a bash command as a character string. 
#' Needed mostly because alignments are run on clusters, but not local machines.
#'
#' @title generate command for STAR aligner
#' @param path_to_STAR path to STAR aligner
#' @param sequence_fasta path to FASTA file with probe sequences (most likely, you will have to run alignment on a cluster, and copy the file from you local machine to the cluster)
#' @param genomeDir path to reference genome annotation, build for STAR
#' @param organism character, "Human", "Mouse" or "Dmel". Note! STAR parameters differ for Drosophila!
#' @param outfile outFileNamePrefix, add a dot in the end.
#' @param runThreadN outBAMsortingThreadN in STAR, 6 by default
#' @return returns a character, STAR command
#' @examples
#' path_to_STAR = "/path/to/STAR"
#' sequence_fasta = "/path/to/sequence.fasta"
#' genomeDir = "/path/to/reference/genome/"
#' organism = "Human"
#' outfile = "output_file."
#' star.command.make(path_to_STAR, sequence_fasta, genomeDir, organism, outfile)
#' @author Vladislava Milchevskaya \email{milchv@gmail.com}
#' @export
star.command.make = function(path_to_STAR,
                             sequence_fasta,
                             genomeDir,
                             organism,
                             outfile,
                             runThreadN = 6)
{
  if(!(organism == "Human" || organism == "Mouse" || organism == "Dmel"))
  {
    warning(message = "organism not Human, Mouse or Dmel")
  }
  readFilesIn = sequence_fasta
  star.parameters = list(
    "genomeDir" = genomeDir,
    "readFilesIn" = readFilesIn,
    #     "readFilesCommand" = "zcat -fc",
    "runThreadN"= runThreadN,
    "genomeLoad"="NoSharedMemory",
    "limitBAMsortRAM"="128000000000", # required if shared memory
    "outFilterMultimapNmax"="20",
    "outFilterMultimapScoreRange"="1",
    "outFilterMismatchNmax"="999", # maximum number of mismatches per pair, 999 = switches off
    "outFilterMismatchNoverLmax"="0.04", # relativ maximum number of mismatches per pair, e.g. for 2 x 100 (paired end) leading to: 0.04 * 200 = 8
    "outFilterMatchNmin"="16",
    "outFilterMatchNminOverLread"="0.66",
    "outFilterScoreMinOverLread"="0.66",
    "outFilterType"="BySJout", # reduces the number of "spurious" junctions
    "alignIntronMin"="20", # minimum intron length
    "alignIntronMax"="1000000", # maximal intron length, the majority should be 500kb-750kb
    #     "alignMatesGapMax"="1000000", # maximum genomic distance between mates
    "alignSJoverhangMin"="8", # minimum overhang for unannotated junctions
    "alignSJDBoverhangMin"="1", # minimum overhang for annotated junctions
    "alignSoftClipAtReferenceEnds"="No", #  option which prevents soft clipping of alignments at the reference (chromosome) ends, for compatibility with Cufflinks/Cuffmerge
    "chimSegmentMin"="20", #  20 for 2 x 75 = a chimeric alignment with 130b on one chromosome and 20b on the other will be output, while 135 + 15 won't be.
    "quantMode"="TranscriptomeSAM",
    #    "quantTranscriptomeBAMcompression"="1", #set the transcriptome compression level
    "outFileNamePrefix"=outfile, # add a dot at the end
    "outSAMattributes"="Standard",
    "outSAMmode"="Full",
    "outSAMstrandField"="intronMotif",
    "outReadsUnmapped"="Fastx",
    #"outSAMtype"="BAM SortedByCoordinate",
    "outBAMcompression"="1",
    "outBAMsortingThreadN"=as.character(runThreadN), 
    "outStd"="BAM_SortedByCoordinate"
  )
  #  }
  if (organism == "Dmel")
  {
    star.parameters[["alignIntronMax"]]<- "500000"
    star.parameters[["outFilterMismatchNmax"]]<- "1"
    star.parameters[["outFilterScoreMinOverLread"]]<- "0"
    star.parameters[["outFilterMatchNminOverLread"]]<- "0"
    star.parameters[["alignSJoverhangMin"]]<- "8"
    star.parameters[["alignSJDBoverhangMin"]]<- "1"
  }
  mm = reshape2::melt(star.parameters)
  colnames(mm) = c("value", "attribute")
  mm_vector = apply(mm, 1, FUN = function(z){paste0(" --", z["attribute"], " ", z["value"])})
  star.params.string = paste(mm_vector, collapse = '')
  res = paste(path_to_STAR, star.params.string, paste0(" > ", outfile, "bam") )
  return(res)
}
