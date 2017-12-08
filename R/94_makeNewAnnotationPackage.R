#' Short description string here
#'
#' Long description string
#'
#' @title What function does (short)
#' @param alignment.dir path to the directory that contains alignment files
#' @param Annotation regerence genome annotation \code{data.frame}
#' @param outputDir output directory
#' @param level character string, can take values "gene" or "transcript", specifying if new probe sets should be formed per gene or per transcript
#' @param min_probe_number numeric, minimal number or probes in a probe set
#' @param pkgNameSUFFIX character, suffix to the package name, usually starts with dor: '.hereIsYourSuffix'
#' @param quiet logical
#' @param fastaDir directory used for generation of the FASTA files
#' @return a \code{"output type"}.
#' @seealso \code{function or class name}
#' @author Vladislava Milchevskaya \email{milchv@gmail.com}
#' @examples
#' data("Alignments_class_example")
#' data("Annotation_example")
#' data("seed_example")
#' alignmentDir <- 
#'   unlist(strsplit(dir(system.file("extdata",package="pdProbeRemap"), 
#'                       pattern="example_drosophila.Aligned.out.sam",
#'                       full.names=TRUE), split = "example_drosophila.Aligned.out.sam"))[1]
#' 
#' outputDir_example = system.file("data", package = "pdProbeRemap")
#' 
#' \dontrun{
#' makeNewAnnotationPackage(alignment.dir = alignmentDir, 
#'                          Annotation = Annotation_example, 
#'                          outputDir = outputDir_example, 
#'                          level = "gene", 
#'                          min_probe_number = 1, 
#'                          pkgNameSUFFIX = ".example")
#'                          }
#' @export
makeNewAnnotationPackage <- function(alignment.dir, 
                                     Annotation,
                                     outputDir,
                                     level,
                                     min_probe_number, 
                                     quiet=FALSE, 
                                     pkgNameSUFFIX, 
                                     fastaDir)
{
    if(!file.exists( file.path(fastaDir, "tmpProbesTableDir/seed3.RData")))
    {stop("Please provide the directory used to generate fasta files in the previous step (fastaDir). \n")}
    
    load(file.path(fastaDir, "tmpProbesTableDir/seed3.RData"))
    new_pd <- makeNewAnnotationPackage_inner(alignment.dir = alignment.dir, 
                                   Annotation = Annotation,
                                   package_seed = seed3,
                                   outputDir = outputDir,
                                   level = level,
                                   min_probe_number = min_probe_number, 
                                   quiet=quiet, 
                                   pkgNameSUFFIX = pkgNameSUFFIX)
    # debugging msg start
    cat(paste0(names(new_pd)))
    # debugging msg end

    invisible(new_pd)
}


# inner function #
makeNewAnnotationPackage_inner <- function(alignment.dir, 
                                           Annotation,
                                           package_seed,
                                           outputDir,
                                           level,
                                           min_probe_number, 
                                           quiet=FALSE, 
                                           pkgNameSUFFIX
)
{
    stopifnot(class(pkgNameSUFFIX) == "character")
    stopifnot(class(quiet) == "logical")
  
    # checks for input parameters for alignments2parsedData inside the function  
    NEW_ParsedData <- alignments2parsedData(alignment.dir = alignment.dir, 
                                            Annotation = Annotation,
                                            outputDir = outputDir,
                                            level = level,
                                            min_probe_number = min_probe_number, 
                                            package_seed = package_seed)
    
    make_package_from_parsedData(parsedData = NEW_ParsedData, 
                                 object = package_seed[[1]], #seed
                                 destDir = outputDir,
                                 quiet = quiet, 
                                 pkgNameSUFFIX = pkgNameSUFFIX,
                                 annotation_type = package_seed[[3]])
  
    return(NEW_ParsedData)
}