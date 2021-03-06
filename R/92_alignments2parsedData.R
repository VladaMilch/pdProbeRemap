#' Short description string here
#'
#' Long description string
#'
#' @title What function does (short)
#' @param alignment.dir path to the directory that contains alignment files
#' @param Annotation regerence genome annotation \code{data.frame}
#' @param outputDir output directory
#' @param level character string, can take values "gene" or "transcript", 
#' specifying if new probe sets should be formed per gene or per transcript
#' @param min_probe_number minimal number or probes in a probe set
#' @param package_seed seed generated by function 
#' \code{makeOriginalPackageObject}
#' @return a \code{"output type"}.
#' @seealso \code{function or class name}
#' @author Vladislava Milchevskaya \email{milchv@gmail.com}
alignments2parsedData <- function(alignment.dir, 
                                  Annotation,
                                  outputDir,
                                  level,
                                  min_probe_number, 
                                  package_seed)
{
  stopifnot(dir.exists(alignment.dir))
  
    stopifnot(class(Annotation) == "data.frame")
  stopifnot("exon" %in% Annotation$feature)
  stopifnot(length(grep("transcript|RNA|rna|Rna", Annotation$feature)) > 0)
  
  stopifnot(dir.exists(outputDir))
  
  stopifnot(level %in% c("gene", "transcript"))
  
  stopifnot(class(min_probe_number) == "numeric")
  stopifnot(min_probe_number >=1 )
  
  stopifnot(  setequal(names(package_seed),   
                       c("object", "parsedDataOLD", "AnnotationType")))
  stopifnot(class(package_seed[[1]]) %in% 
              c("AffyGenePDInfoPkgSeed", "AffyExpressionPDInfoPkgSeed"))
  stopifnot(package_seed[[3]] %in% c("pgf", "cdf"))
  
  
  ParsedData <- package_seed$parsedDataOLD
  Als <- processAlignments(alignment.dir = alignment.dir)
  PS_mus = alignments2probesetsRaw(Alignments = Als, 
                                   Annotation = Annotation, level = level)
  PS_mus_cleaned <- cleanProbeSets(ProbeSets = PS_mus,
                                   level = level,
                                   min_probe_number = min_probe_number)
  
  PS_controls <- addControlProbes(ProbeSetsClean = PS_mus_cleaned, 
                                  outputDir = outputDir, 
                                  alignmentType = package_seed[[3]])
  
  shift = max(as.numeric(PS_controls$fsetid))
  PSD <- make_objects_for_reannotation(Probesets = PS_mus_cleaned, 
                                       level = level, 
                                       Annotation = Annotation, 
                                       shift = shift)
  
  save(PSD, file = file.path(outputDir, paste0("ProbeSetData", ".RData")))
  
  check_ProbeSetDataAnnotation(PSD$ProbeSetDataAnnotation)
  check_ProbeSetDataFrame(PSD$ProbeSetDataFrame)
  
  if (package_seed[[3]] == "cdf")
  {NewParsedData <- 
    build_new_ParsedDataOld_object_CDF(level = level, 
                                       ParsedDataOld = ParsedData, 
                                       PS_controls = PS_controls, 
                                       Annotation = Annotation, 
                                       ProbeSetData = PSD)
  }
  
  if (package_seed[[3]] == "pgf")
  # {
  #   package_seed[[3]] <- "pgf2cdf"
  #   NewParsedData <- 
  #     build_new_ParsedDataOld_object_CDF(level = level, 
  #                                        ParsedDataOld = ParsedData, 
  #                                        PS_controls = PS_controls, 
  #                                        Annotation = Annotation, 
  #                                        ProbeSetData = PSD)
  # }
  {NewParsedData <-
    build_new_ParsedDataOld_object_PGF(level = level,
                                       ParsedDataOld = ParsedData,
                                       PS_controls = PS_controls,
                                       Annotation = Annotation,
                                       ProbeSetData = PSD,
                                       outputDir = outputDir)
  }
  
  # debugging msg start
  cat(paste0("From inside alignments2parsedData:  ", names(NewParsedData)))
  # debugging msg end
  return(NewParsedData)
}

