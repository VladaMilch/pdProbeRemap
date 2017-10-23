# transfered from /Users/milchevs/ownCloud/ProbeRemap/Functions/make.original.package.object.function.R #
make.original.package.object.fromPGF <- function(
  pgf, clf, prob, mps, trans,
  author, email, 
  organism, species, outputDir
)
{
  if(! all( sapply(c(pgf, clf, prob, mps, trans), length) == 1))
    stop(message = "\nError: incorrect input for orifinal annotation. \n")

  object <- new("AffyGenePDInfoPkgSeed",
                pgfFile = pgf, clfFile = clf, coreMps=mps, transFile=trans,
                probeFile = prob, author = author,
                email = email,
                biocViews = "AnnotationData",
                organism = organism, species = species)
  batch_size=10000
  quiet=FALSE
  unlink=FALSE
  geneArray <- object@geneArray
  stopifnot(is.logical(geneArray))
  object@chipName <- chipName(object)
  

  if (geneArray){
    msg <- "Parsing annotation files for Affymetrix Gene ST Array"
  }else{
    msg <- "Parsing annotation files for Affymetrix Exon ST Array"
  }

  msgBar()
  message(msg)
  message("PGF.........: ", basename(object@pgfFile))
  message("CLF.........: ", basename(object@clfFile))
  message("Probeset....: ", basename(object@probeFile))
  message("Transcript..: ", basename(object@transFile))
  message("Core MPS....: ", basename(object@coreMps))
  if (!geneArray){
    message("Full MPS....: ", basename(object@fullMps))
    message("Extended MPS: ", basename(object@extendedMps))
  }
  msgBar()

  PDolD <- combinePgfClfProbesetsMps_andWriteAllProbesFile(object@pgfFile,
                                     object@clfFile,
                                     object@probeFile,
                                     object@coreMps,
                                     object@fullMps,
                                     object@extendedMps,
                                     verbose=!quiet,
                                     geneArray=geneArray, 
                                     outputDir = outputDir)
  res <- list(object, PDolD)
  names(res) <- c("object", "parsedDataOLD")
  return(res)

}



## INPUT ###
make.original.package.object.fromCDF <- function(
  cdf, cel, tab,
  author, email, 
  organism, species, outputDir
)
{
  if(! all( sapply(c(cdf, cel, tab), length) == 1))
    stop(message = "\nError: incorrect input for orifinal annotation. \n")

  object <- new("AffyExpressionPDInfoPkgSeed",
                cdfFile = cdf, celFile = cel,
                tabSeqFile = tab,  author = author,
                email = email,
                biocViews = "AnnotationData",
                organism = organism, species = species)
  batch_size=10000
  quiet=FALSE
  unlink=FALSE

  object@chipName <- chipName(object)
  
  msgBar()
  message("Parsing annotation files for Affymetrix Expression array")
  message("CDF...............: ", basename(object@cdfFile))
  message("CEL...............: ", basename(object@celFile))
  message("Sequence TAB-Delim: ", basename(object@tabSeqFile))
  msgBar()

  parsedData <- parseCdfCelProbe_andWriteAllProbesFile(object@cdfFile,
                                                       object@celFile,
                                                       object@tabSeqFile,
                                                       verbose=!quiet, 
                                                       outputDir = outputDir)
  hasMM <- nrow(parsedData[["mmFeatures"]]) > 0
  #############################################################################

  res <- list(object, parsedData)
  names(res) <- c("object", "parsedDataOLD")
  return(res)

}


##########################
#' read in original array annotation
#'
#' read in original array annotation from Affymetrix (CDF-based or PGF-based),
#' and produce an object required for an annotation package,
#' compatible with \code{oligo} or \code{affy}.
#'
#' @title make original annotation object
#' @param originalAnnotationDir directory that contains annotation files: 
#' \code{".pgf"}, \code{".clf"},  \code{".probeset.csv"},  \code{".mps"} and  
#' \code{".transcript.csv"} files, if annotation is PGF-based;
#' \code{".cdf"} -- Chip Definition file, \code{".CEL"} -- example CEL file, 
#' \code{"probe_tab"} -- tabulated .txt probes file
#' @param celFilePath path to CEL file, if not present in 
#' \code{originalAnnotationDir} and annotation is CDF-based
#' @param author \code{character}, author
#' @param pkgNameSUFFIX \code{character}, suffix to package name
#' @param email \code{character}, email
#' @param organism \code{character} (ex: "Human")
#' @param species \code{character} (ex: "Homo Sapiens")
#' @param outputDir path to user-spesified output directory
#' @return a \code{"output type"}.
#' @seealso \code{function or class name}
#' @author Vladislava Milchevskaya \email{milchv@gmail.com}
#' @examples 
#' alignmentDir = "/Users/milchevs/Documents/Biology/PROJECTS/Paul_Bertone/ReANNot/MoGene2/"
#' load("/Users/milchevs/Documents/Biology/PROJECTS/Paul_Bertone/ReANNot/Annotations/AnnotationDataFrame.Mus_musculus.GRCm38.82.RData")
#' dir1 = "/Users/milchevs/ownCloud/ProbeRemapSteps14/inst/extdata/mogene2_affyLibFiles_pgfBased/"
#' outDir = "/Users/milchevs/ownCloud/ProbeRemapSteps14/OutPutDir/"
#' o1 = makeOriginalPackageObject(originalAnnotationDir = dir1, organism = "Mouse", species = "Mus musculus", outputDir = outDir )
#' @export
makeOriginalPackageObject <-function(originalAnnotationDir,
                                        #pgf = "InputFilePath", 
                                        #clf  = "InputFilePath", 
                                        #prob  = "InputFilePath", 
                                        #mps  = "InputFilePath", 
                                        #trans  = "InputFilePath",
                                        #cdf = "InputFilePath",
                                        celFilePath = "InputFilePath",
                                        #tab = "InputFilePath",
                                        author = "Unknown", email = "un@known", 
                                        pkgNameSUFFIX = '.originalPackage',
                                        organism, species, outputDir)
{
  if(!(organism %in% c("Dmel", "Human", "Mouse")))
  {stop(message = "Organism can only take values 
        'Human', 'Mouse' or 'Dmel' (for Drosophila). Pick one of these.")}
  
  if(!(species %in% 
       c("Drosophila melanogaster", "Homo sapiens", "Mus musculus")))
  {stop(message = "species can only take values 
        'Homo sapiens', 
        'Mus musculus' or 
        'Drosophila melanogaster'. Pick one of these.")}
  
  if(organism == "Dmel" | species == "Drosophila melanogaster")
  {stopifnot(all(c(organism == "Dmel", species == "Drosophila melanogaster")))}
  
  if(organism == "Human" | species == "Homo sapiens")
  {stopifnot(all(c(organism == "Human", species == "Homo sapiens")))}
  
  if(organism == "Mouse" | species == "Mus musculus")
  {stopifnot(all(c(organism == "Mouse", species == "Mus musculus")))}
  
  if(!(dir.exists(originalAnnotationDir)))
  {stop("Wrong path to original Annotation directory! 
        Input existing directory path.")}
  
  if(!(dir.exists(outputDir)))
  {stop("Wrong path to Output directory! 
        Input existing directory path.")}
  
  baseDir <- originalAnnotationDir
  (pgf <- list.files(baseDir, pattern = ".pgf", ignore.case = TRUE,
                     full.names = TRUE))
  (cdf <- list.files(baseDir, pattern = ".cdf", ignore.case = TRUE,
                     full.names = TRUE))
  #   if (length(pgf) > 1)
  #   {stop(message = "More than 1 PGF file in originalAnnotationDir. \n") }
  if (!(identical(pgf, character(0))))
  {
    if( length(pgf) == 1 & pgf  != "InputFilePath")
    {
      (pgf <- list.files(baseDir, pattern = ".pgf", 
                         ignore.case = TRUE,
                         full.names = TRUE))
      (clf <- list.files(baseDir, pattern = ".clf",
                         ignore.case = TRUE,
                         full.names = TRUE))
      (prob <- list.files(baseDir, pattern = ".probeset.csv",
                          ignore.case = TRUE,
                          full.names = TRUE))
      mps <- list.files(baseDir, pattern = "mps$", 
                        ignore.case = TRUE,
                        full.names = TRUE)
      trans <- list.files(baseDir, pattern="transcript.file.csv|transcript.csv",
                          ignore.case = TRUE,full.names=TRUE)

      stopifnot(  length( c(pgf, clf, prob, mps, trans) ) == 5)


      res = make.original.package.object.fromPGF(
        pgf = pgf, clf = clf, prob = prob, mps = mps, trans = trans,
        author=author, 
        email=email, 
        organism = organism, 
        species = species, 
        outputDir = outputDir)
      
      simpleMessage( "PGF-based Affymetrix chip annotation parced. \n")
      AnnotationType = "pgf"
      res_wtype = c(res, AnnotationType)
      names(res_wtype)[3] <- "AnnotationType"
      return(res_wtype)
    }
  }
  if (!(identical(cdf, character(0))))  
    {
        if(length(cdf)==1 & cdf != "InputFilePath") 
        {
            (cdf <- list.files(baseDir, pattern = "cdf", ignore.case = TRUE,
                               full.names = TRUE))
            (cel <- list.files( baseDir, pattern = ".CEL", ignore.case = TRUE,
                        full.names = TRUE))
            (tab <- list.files(baseDir, pattern = "probe_tab",
                                ignore.case = TRUE,
                                full.names = TRUE))
            if (celFilePath == "InputFilePath" & 
                !(identical(cel, character(0))))
            {
            if(length(cel) > 1)
            {stop(message = "Only one CEL file must be provided!")}
            }
            if (celFilePath != "InputFilePath")
            {
                if(!(identical(cel, character(0))))
                {
                    warning("Two CEL files available: 
                        in the annotation directory 
                        and another one provided as input to this function. 
                        The one provided as input to the function 
                        will be used! \n")
                }
                cel = celFilePath
            }
            stopifnot(  length( c(cdf, cel, tab) ) == 3)

            res = make.original.package.object.fromCDF(
                cdf = cdf, 
                cel = cel, 
                tab = tab,
                author=author, 
                email=email, 
                organism = organism, 
                species = species, 
                outputDir = outputDir)
            simpleMessage( "CDF-based Affymetric chip annotation parced. \n")
            AnnotationType = "cdf"
            res_wtype = c(res, AnnotationType)
            names(res_wtype)[3] <- "AnnotationType"
            return(res_wtype)  
        }
    }
    if ((identical(cdf, character(0))) & (identical(pgf, character(0))) )
    { simpleMessage("Error: Original Annotation Input Files are missing. \n")}
}

#PDolD <- 
#  make.original.package.object(originalAnnotationDir, originalPackageDir)

# head(PDolD$featureSet) ##
# head(PDolD$pmFeatures) ##
# head(PDolD$mmFeatures) ##
# head(PDolD$geometry)
# head(PDolD$pmSequence)
# head(PDolD$mmSequence)
# (PDolD$chrom_dict)
# (PDolD$level_dict)
# (PDolD$type_dict)
# head(PDolD$core) ##
#
#
