#' Fasta files from a pd.annotation object
#'
#' Fasta files with pobe sequences generated from a parsed data object
#' (probe sequences and ids from original annotation).
#'
#' @title write probes in FASTA 
#' @param originalAnnotationDir directory that contains annotation files: 
#' \code{".pgf"}, \code{".clf"},  \code{".probeset.csv"},  
#' \code{".mps"} and  \code{".mps"} files, if annotation is PGF-based;
#' \code{".cdf"} -- Chip Definition file, \code{".CEL"} -- example CEL file, 
#' \code{"probe_tab"} -- tabulated .txt probes file
#' @param celFilePath path to CEL file, if not present in 
#' \code{originalAnnotationDir} and annotation is CDF-based
#' @param organism \code{character} (ex: "Human")
#' @param species \code{character} (ex: "Homo Sapiens")
#' @param author \code{character}, name of the annotation package author
#' @param email \code{character} e-mail of the annotation package author 
#' (must contain "@")
#' @param outputDir path to user-spesified output directory
#' @param reverse \code{logical}, if \code{reverse = TRUE}, probe sequences 
#' are written in reverse order (not reverse complementary)
#' @return writes \code{"FASTA"} file in directory 
#' \code{outpurDir/Probe_Sequences_Fasta/} (created the directory if needed)
#' @author Vladislava Milchevskaya \email{milchv@gmail.com}
#' @import BiocGenerics
#' @import GenomicRanges
#' @import IRanges
#' @importFrom seqinr write.fasta
#' @importFrom oligo cleanPlatformName
#' @examples 
#' dir1 = "/Users/milchevs/ownCloud/ProbeRemapSteps14/inst/extdata/mogene2_affyLibFiles_pgfBased/"
#' outDir = "/Users/milchevs/ownCloud/ProbeRemapSteps14/OutPutDir/"
#' dir.create(outDir)
#' annotation2fasta(originalAnnotationDir = dir1, 
#'                 organism = "Mouse", species = "Mus musculus", 
#'                 author = "Vladislava Milchevskaya",
#'                 email = "milchv@gmail.com",
#'                 outputDir = outDir, 
#'                 reverse = FALSE)
#' @export 

annotation2fasta <- function(originalAnnotationDir,
                            celFilePath = "InputFilePath", author, email, 
                            organism, species, outputDir, reverse = FALSE)
{
    stopifnot(class(author) == "character")
    stopifnot(class(email) == "character")
    stopifnot(grep("@", email) == 1)

    if(class(reverse) != "logical")
    {stop("Parameter reverse can only take values TRUE or FALSE.")}

    fastaDir <- "Probe_Sequences_Fasta"
    mainDir = outputDir
    subDir = fastaDir
    if(dir.exists(file.path(mainDir, subDir)))
    {
    ll = list.files(path = file.path(mainDir, subDir), pattern = ".FASTA", 
                    ignore.case = TRUE)
        if(length(ll) > 0)
        {
            warning(paste0("There are FASTA files in directory ", 
                        file.path(mainDir, subDir)), immediate. = TRUE)
        }
    }

    if(!(organism %in% c("Dmel", "Human", "Mouse")))
    {stop(message = "Organism can only take values 'Human', 
        'Mouse' or 'Dmel' (for Drosophila). Pick one of these.")}
    
    if(!(species %in% c("Drosophila melanogaster", 
                        "Homo sapiens", 
                        "Mus musculus")))
    {stop(message = "species can only take values 'Homo sapiens', 
        'Mus musculus' or 'Drosophila melanogaster'. Pick one of these.")}
    
    if(organism == "Dmel" | species == "Drosophila melanogaster")
    {
        stopifnot(
            all(
                c(
                    organism == "Dmel", 
                    species == "Drosophila melanogaster")))
    }
    
    if(organism == "Human" | species == "Homo sapiens")
    {stopifnot(all(c(organism == "Human", species == "Homo sapiens")))}
    
    if(organism == "Mouse" | species == "Mus musculus")
    {stopifnot(all(c(organism == "Mouse", species == "Mus musculus")))}
    
    if(!(dir.exists(originalAnnotationDir)))
    {stop("Wrong path to original Annotation directory! 
        Input existing directory path.")}
    
    if(!(dir.exists(outputDir)))
    {stop("Wrong path to Output directory! Input existing directory path.")}
    
    seed3 <- makeOriginalPackageObject(
        originalAnnotationDir = originalAnnotationDir,
        celFilePath = celFilePath, author = author, email = email, 
        organism = organism, species = species, outputDir = outputDir)
    
    stopifnot(dir.exists(file.path(outputDir, "tmpProbesTableDir")))
    save(seed3, file = file.path(outputDir, "tmpProbesTableDir", "seed3.RData"))
    
    object <- seed3[["object"]]
    ParcedDataOld <- seed3[["parsedDataOLD"]]
    annotationType <- seed3[[3]]
    
    species <- object@species
    organism <- object@organism
    chip <- pdInfoBuilder::chipName(object)
    platformName <- cleanPlatformName(chip)
    
    if(seed3[[3]] == "cdf")
    {
        probesCDF <- read.csv(file = file.path(outputDir, 
                                             "tmpProbesTableDir",
                                             "ProbesCDF.csv"), 
                            header = TRUE)
        probesCDF = probesCDF[!is.na(probesCDF$sequence), ]
        PP = probesCDF[ ,c("fid", "sequence")]
    }
    if(seed3[[3]] == "pgf")
    {
        probesPGF <- read.csv(file = file.path(outputDir, 
                                             "tmpProbesTableDir",
                                             "ProbesPGF.csv"), 
                            header = TRUE)
        probesPGF = probesPGF[!is.na(probesPGF$sequence),  ]
        PP = probesPGF[ ,c("fid", "sequence")]
    }
    
    ifelse(!dir.exists(file.path(mainDir, subDir)), 
        dir.create(file.path(mainDir, subDir)), FALSE)
    message(paste0("Probe_Sequences_Fasta directory created in ", 
                file.path(outputDir)))
    
    if(any(duplicated(PP)))
    { 
        PP = PP[!duplicated(PP), ]
    }
    

    #### write fasta ### 
    
    if (!reverse)
    {
        fasta_name = paste0(chip, "_probe_sequences_to_align.fasta")
        fasta_path_full = file.path(outputDir, fastaDir, fasta_name)
        if(file.exists(fasta_path_full))
        {warning(paste0("FASTA file with name \n", 
                    fasta_path_full,"\nexists! Can not write FASTA file."))}
    
        if(!file.exists(fasta_path_full))
        {
            simpleMessage("\nWriting FASTA file with probe sequences...")
            NewProbeIds = as.numeric(PP$fid)
            seqinr::write.fasta(as.list(PP$sequence), NewProbeIds, nbchar = 60, 
                            file.out = fasta_path_full,
                            open = "w")
            message("OK\n")
        }
    fp = fasta_path_full
    }
    
      
    if(reverse)
    {
        fasta_reverse_name = 
            paste0(chip, "_probe_sequences_to_align_REVERSE.fasta")
        fasta_reverse_path_full = 
            file.path(outputDir, fastaDir, fasta_reverse_name)
        if(file.exists(fasta_reverse_name))
        {warning(paste0("FASTA file with name \n", 
                    fasta_reverse_name, 
                    "\nexists! File was NOT overwritten."))} 
    
        if(!file.exists(fasta_reverse_name))
        {
            simpleMessage(
            "Writing FASTA file with REVERSE 
            (not complementary) probe sequences..."
            )
            NewProbeIds_reverseprobes = as.numeric(PP$fid)
            reverse_probes = as.vector(reverse(DNAStringSet(PP$sequence)))
            seqinr::write.fasta(as.list(reverse_probes), 
                                NewProbeIds_reverseprobes, 
                                nbchar = 60, 
                                file.out = fasta_reverse_path_full, 
                                open = "w")
            message("OK\n") 
        }

        fp = fasta_reverse_path_full
    }
    
    simpleMessage(paste0("\nFASTA file with probe sequences:    \n",  fp, "\n"))
    
    #BeforeAl = c(object_ParcedData_annotationType, 
    # fasta_path_full, fasta_reverse_path_full)
    #return(BeforeAl)
    return(fp)
    
}

############# ##########
# parsed2fasta <- function(object_ParcedData_annotationType, 
#                          outputDir, reverse = FALSE)
# {
#   message(Sys.time())
#   object <- object_ParcedData_annotationType[["object"]]
#   ParcedDataOld <- object_ParcedData_annotationType[["parsedDataOLD"]]
#   annotationType <- object_ParcedData_annotationType[[3]]
#   
#   species <- object@species
#   organism <- object@organism
#   chip <- ame(object)
#   platformName <- cleanPlatformName(chip)
#   message(Sys.time())
#   
#   ######### create fasta files ########
#   ifelse(!dir.exists(outputDir), dir.create(outputDir), FALSE)
#   fastaDir <- "Probe_Sequences_Fasta"
#   mainDir = outputDir
#   subDir = fastaDir
#   ifelse(!dir.exists(file.path(mainDir, subDir)), 
#   dir.create(file.path(mainDir, subDir)), FALSE)
#   message(paste0("Probe_Sequences_Fasta directory created in ", outputDir))
#   message(Sys.time())
#   
#   ParcedDataOld$pmSequence <- 
#       ParcedDataOld$pmSequence[!duplicated(ParcedDataOld$pmSequence),]
#   #dim(ParcedDataOld$pmSequence  )
#   #### write fasta ### 
#   
#   if (!reverse)
#   {
#     fasta_name = paste0(chip, "_probe_sequences_to_align.fasta")
#     fasta_path_full = file.path(outputDir, fastaDir, fasta_name)
#     message(Sys.time())
#     
#     simpleMessage("\nWriting FASTA file with probe sequences...")
#     NewProbeIds = as.numeric(ParcedDataOld$pmSequence$fid)
#     seqinr::write.fasta(as.list(ParcedDataOld$pmSequence$sequence), 
#             NewProbeIds, nbchar = 60, 
#                 file.out = fasta_path_full,
#                 open = "w")
#     message("OK\n")
#     message(Sys.time())
#     
#     fp = fasta_path_full
#   }
#   
#   if(reverse)
#   {
#     
#     fasta_reverse_name = 
# paste0(chip, "_probe_sequences_to_align_REVERSE.fasta")
#     fasta_reverse_path_full = 
#        file.path(outputDir, fastaDir, fasta_reverse_name)
#     message(Sys.time())
#     
#     simpleMessage("Writing FASTA file with REVERSE 
#         (not complementary) probe sequences...")
#     NewProbeIds_reverseprobes = as.numeric(ParcedDataOld$pmSequence$fid)
#     reverse_probes = 
#          as.vector(reverse(DNAStringSet(ParcedDataOld$pmSequence$sequence)))
#     seqinr::write.fasta(as.list(reverse_probes), 
#             NewProbeIds_reverseprobes, nbchar = 60, 
#                 file.out = fasta_reverse_path_full, 
#                 open = "w")
#     message("OK\n")
#     message(Sys.time())
#     
#     fp = fasta_reverse_path_full
#   }
#   
#   simpleMessage(paste0("\nFASTA file with probe sequences:    \n",  fp, "\n"))
#   
#   #BeforeAl = c(object_ParcedData_annotationType, 
#        fasta_path_full, fasta_reverse_path_full)
#   #return(BeforeAl)
#   return(fp)
# } 

