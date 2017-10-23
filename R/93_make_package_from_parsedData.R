#### make the package out of the modified parsedDataNEW object ####
#load("/Users/milchevs/Documents/Biology/PROJECTS/
# Paul_Bertone/PROBE_REANNOTATION/HUGene_ParsedDataListNEWconverted.RData")
#source("/Users/milchevs/Documents/Biology/PROJECTS/
# Paul_Bertone/PROBE_REANNOTATION/make_my_package_main.R")

#rm(parsedData)

make_package_from_parsedData <-   function( parsedData, 
                                            object, #seed
                                            destDir,
                                            quiet=FALSE, 
                                            pkgNameSUFFIX = ".VLADA3",
                                            annotation_type = "PGF_or_CDF")
{
    if (!(annotation_type == "PGF" || 
        annotation_type == "CDF" || 
        annotation_type == "pgf" || 
        annotation_type == "cdf"))
    {stop(message="Error: input manually annotation type! PGF or CDF. \n")}
    if (annotation_type == "PGF" || annotation_type == "pgf")
    {  
        make_package_from_parsedData_fromPGF(      parsedData = parsedData, 
                                               object = object, #seed
                                               destDir = destDir,
                                               quiet = quiet, 
                                               pkgNameSUFFIX = pkgNameSUFFIX, 
                                               newPackageDir = destDir)
    }
    if (annotation_type == "CDF"|| annotation_type == "cdf")
    {
      make_package_from_parsedData_fromCDF(      parsedData = parsedData, 
                                               object = object, #seed
                                               destDir = destDir,
                                               quiet = quiet, 
                                               pkgNameSUFFIX = pkgNameSUFFIX, 
                                               newPackageDir = destDir)
    }
    simpleMessage("\n Done! \n")
}
###################################################
############ from PGF-like annotation #############



make_package_from_parsedData_fromPGF <- 
    function( parsedData, 
            object, #seed
            destDir,
            quiet, 
            pkgNameSUFFIX, 
            newPackageDir)
{
    mainDir = getwd()
    if (file.exists(paste0(mainDir, "/",newPackageDir))) {
      message("newPackageDir exists in the current directory, 
              and is a directory")
    } else if (file.exists(newPackageDir)) {
      message("newPackageDir exists in mainDir but is a file")
      # you will probably want to handle this separately
    } else {
      message("newPackageDir does not exist in mainDir - creating")
      dir.create(file.path(mainDir, newPackageDir))
    }
    geneArray <- object@geneArray
    stopifnot(is.logical(geneArray))
    if (geneArray){
      msg <- "Building annotation package for Affymetrix Gene ST Array"
    }else{
      msg <- "Building annotation package for Affymetrix Exon ST Array"
    }
    #######################################################################
    ## Part i) get array info (chipName, pkgName, dbname)
    #######################################################################
    chip <- chipName(object)
    pkgName <- paste0(cleanPlatformName(chip), pkgNameSUFFIX)
    extdataDir <- file.path(destDir, pkgName, "inst", "extdata")
    dbFileName <- paste(pkgName, "sqlite", sep=".")
    dbFilePath <- file.path(extdataDir, dbFileName)
    #
    print(chip)
    print(pkgName)
    print(extdataDir)
    print(dbFileName)    
    print(dbFilePath) 
    
    #######################################################################
    ## Part iii) Create package from template
    #######################################################################
    pdInfoClass <- ifelse(geneArray, "AffyGenePDInfo", "AffyExonPDInfo")
    syms <- list(MANUF=object@manufacturer,
                 VERSION=object@version,
                 GENOMEBUILD=object@genomebuild,
                 AUTHOR=object@author,
                 AUTHOREMAIL=object@email,
                 LIC=object@license,
                 DBFILE=dbFileName,
                 CHIPNAME=chip,
                 PKGNAME=pkgName,
                 PDINFONAME=pkgName,
                 PDINFOCLASS=pdInfoClass,
                 GEOMETRY=parsedData[["geometry"]])
    templateDir <- system.file("pd.PKG.template",
                               package="pdInfoBuilder")
    createPackage(pkgname=pkgName, destinationDir=destDir,
                  originDir=templateDir, symbolValues=syms,
                  quiet=quiet, unlink = TRUE)
    dir.create(extdataDir, recursive=TRUE)
    
    #######################################################################
    ## Part iv) Create SQLite database
    ## FIX ME: Fix ordering of the tables
    #######################################################################
    conn <- dbConnect(dbDriver("SQLite"), dbname=dbFilePath)
    increaseDbPerformance(conn)
    
    ## Adding new tables
    dbCreateTable(conn,
                  "chrom_dict",
                  chromDictTable[["col2type"]],
                  chromDictTable[["col2key"]])
    dbCreateTable(conn,
                  "level_dict",
                  levelDictTable[["col2type"]],
                  levelDictTable[["col2key"]])
    dbCreateTable(conn,
                  "type_dict",
                  typeDictTable[["col2type"]],
                  typeDictTable[["col2key"]])
    dbCreateTable(conn,
                  "core_mps",
                  mpsSchema[["col2type"]],
                  mpsSchema[["col2key"]])
    if (!geneArray){
      dbCreateTable(conn,
                    "full_mps",
                    mpsSchema[["col2type"]],
                    mpsSchema[["col2key"]])
      dbCreateTable(conn,
                    "extended_mps",
                    mpsSchema[["col2type"]],
                    mpsSchema[["col2key"]])
    }
    ## end adding
    
    dbCreateTable(conn,
                  "featureSet",
                  exonTranscriptionFeatureSetSchema[["col2type"]],
                  exonTranscriptionFeatureSetSchema[["col2key"]])
    
    if (geneArray){
      dbCreateTable(conn, "pmfeature",
                    genePmFeatureSchema[["col2type"]],
                    genePmFeatureSchema[["col2key"]])
    }else{
      dbCreateTable(conn,
                    "pmfeature",
                    exonTranscriptionPmFeatureSchema[["col2type"]],
                    exonTranscriptionPmFeatureSchema[["col2key"]])
    }
    containsMm <- nrow(parsedData[["mmFeatures"]]) > 0
    if (containsMm)
      dbCreateTable(conn,
                    "mmfeature",
                    exonTranscriptionMmFeatureSchema[["col2type"]],
                    exonTranscriptionMmFeatureSchema[["col2key"]])
    
    ## Inserting data in new tables
    dbInsertDataFrame(conn, "chrom_dict", parsedData[["chrom_dict"]],
                      chromDictTable[["col2type"]], !quiet)
    dbInsertDataFrame(conn, "level_dict", parsedData[["level_dict"]],
                      levelDictTable[["col2type"]], !quiet)
    dbInsertDataFrame(conn, "type_dict", parsedData[["type_dict"]],
                      typeDictTable[["col2type"]], !quiet)
    dbInsertDataFrame(conn, "core_mps", parsedData[["core"]],
                      mpsSchema[["col2type"]], !quiet)
    if (!geneArray){
      dbInsertDataFrame(conn, "full_mps", parsedData[["full"]],
                        mpsSchema[["col2type"]], !quiet)
      dbInsertDataFrame(conn, "extended_mps", parsedData[["extended"]],
                        mpsSchema[["col2type"]], !quiet)
    }
    ## end inserting
    
    dbInsertDataFrame(conn, "featureSet", parsedData[["featureSet"]],
                      exonTranscriptionFeatureSetSchema[["col2type"]], !quiet)
    ###
    
    ###
    
    if (geneArray){
      dbInsertDataFrame(conn, "pmfeature", parsedData[["pmFeatures"]],
                        genePmFeatureSchema[["col2type"]], !quiet)
    }else{
      dbInsertDataFrame(conn, "pmfeature", parsedData[["pmFeatures"]],
                        exonTranscriptionPmFeatureSchema[["col2type"]], !quiet)
    }
    if (containsMm)
      dbInsertDataFrame(conn, "mmfeature", parsedData[["mmFeatures"]],
                        exonTranscriptionMmFeatureSchema[["col2type"]], !quiet)
    
    dbCreateTableInfo(conn, !quiet)
    
    ## Create indices
    if (geneArray){
      dbCreateIndex(conn, 
                    "idx_pmfsetid", "pmfeature", "fsetid", 
                    FALSE, verbose=!quiet)
      dbCreateIndex(conn, "idx_pmfid", "pmfeature", 
                    "fid", FALSE, verbose=!quiet)
    }else{
      dbCreateIndicesPm(conn, !quiet)
    }
    dbCreateIndicesFs(conn, !quiet)
    dbCreateIndex(conn, "idx_core_meta_fsetid", 
                  "core_mps", "meta_fsetid", FALSE, verbose=!quiet)
    dbCreateIndex(conn, "idx_core_fsetid", 
                  "core_mps", "fsetid", FALSE, verbose=!quiet)
    if (!geneArray){
      dbCreateIndex(conn, "idx_full_meta_fsetid", 
                    "full_mps", "meta_fsetid", 
                    FALSE, verbose=!quiet)
      dbCreateIndex(conn, "idx_full_fsetid", "full_mps", 
                    "fsetid", FALSE, verbose=!quiet)
      dbCreateIndex(conn, "idx_extended_meta_fsetid", 
                    "extended_mps", "meta_fsetid", FALSE, verbose=!quiet)
      dbCreateIndex(conn, "idx_extended_fsetid", "extended_mps", 
                    "fsetid", FALSE, verbose=!quiet)
    }
    
    if (containsMm){
      dbCreateIndex(conn, "idx_mmfsetid", "mmfeature", 
                    "fsetid", FALSE, verbose=!quiet)
      dbCreateIndex(conn, "idx_mmfid", "mmfeature", 
                    "fid", FALSE, verbose=!quiet)
    }
    
    dbGetQuery(conn, "VACUUM")
    dbDisconnect(conn)
    
    #######################################################################
    ## Part v) Save sequence DataFrames
    ## FIX ME: Fix ordering of the tables to match xxFeature tables
    #######################################################################
    datadir <- file.path(destDir, pkgName, "data")
    dir.create(datadir)
    pmSequence <- parsedData[["pmSequence"]]
    pmSeqFile <- file.path(datadir, "pmSequence.rda")
    if (!quiet) message("Saving DataFrame object for PM.")
    save(pmSequence, file=pmSeqFile, compress='xz')
    if (containsMm){
      mmSequence <- parsedData[["mmSequence"]]
      mmSeqFile <- file.path(datadir, "mmSequence.rda")
      if (!quiet) message("Saving DataFrame object for MM.")
      save(mmSequence, file=mmSeqFile, compress='xz')
    }
    
    
    #######################################################################
    ## Part vi) Save NetAffx Annotation to extdata
    #######################################################################
    if (!quiet) message("Saving NetAffx Annotation... ", appendLF=FALSE)
    netaffxProbeset <- annot2fdata(object@probeFile)
    save(netaffxProbeset, file=file.path(extdataDir,
                                         'netaffxProbeset.rda'), compress='xz')
    netaffxTranscript <- annot2fdata(object@transFile)
    save(netaffxTranscript, file=file.path(extdataDir,
                                           'netaffxTranscript.rda'), 
         compress='xz')
    if (!quiet) msgOK()
    
    if (!quiet) message("Done.")
    #})
    
}


### END: now the package is made ###

#make_package_from_parsedData(tmp)


make_package_from_parsedData_fromCDF <- function( parsedData, 
                                                  object, #seed
                                                  destDir,
                                                  quiet, 
                                                  pkgNameSUFFIX, 
                                                  newPackageDir)
{
    msgBar()
    message("Building annotation package for Affymetrix Expression array\n")
    message("CDF...............: ", basename(object@cdfFile), "\n")
    message("CEL...............: ", basename(object@celFile), "\n")
    message("Sequence TAB-Delim: ", basename(object@tabSeqFile), "\n")
    msgBar()
    
    #######################################################################
    ## Part i) get array info (chipName, pkgName, dbname)
    #######################################################################
    chip <- chipName(object)
    pkgName <- paste0(cleanPlatformName(chip), pkgNameSUFFIX)
    extdataDir <- file.path(destDir, pkgName, "inst", "extdata")
    dbFileName <- paste(pkgName, "sqlite", sep=".")
    dbFilePath <- file.path(extdataDir, dbFileName)
    
    #######################################################################
    ## Part ii) parse data obrained as an argument
    #######################################################################
    hasMM <- nrow(parsedData[["mmFeatures"]]) > 0
    #######################################################################
    ## Part iii) Create package from template
    #######################################################################
    syms <- list(MANUF=object@manufacturer,
                 VERSION=object@version,
                 GENOMEBUILD=object@genomebuild,
                 AUTHOR=object@author,
                 AUTHOREMAIL=object@email,
                 LIC=object@license,
                 DBFILE=dbFileName,
                 CHIPNAME=chip,
                 PKGNAME=pkgName,
                 PDINFONAME=pkgName,
                 PDINFOCLASS="AffyExpressionPDInfo",
                 GEOMETRY=parsedData[["geometry"]])
    templateDir <- system.file("pd.PKG.template",
                               package="pdInfoBuilder")
    createPackage(pkgname=pkgName, destinationDir=destDir,
                originDir=templateDir, symbolValues=syms,
                quiet=quiet, unlink = TRUE)
  dir.create(extdataDir, recursive=TRUE)
  #######################################################################
  ## Part iv) Create SQLite database
  ## FIX ME: Fix ordering of the tables
  #######################################################################
  conn <- dbConnect(dbDriver("SQLite"), dbname=dbFilePath)
  increaseDbPerformance(conn)
  dbCreateTable(conn,
                "featureSet",
                affyHTExpressionFeatureSetSchema[["col2type"]],
                affyHTExpressionFeatureSetSchema[["col2key"]])
  
  dbCreateTable(conn,
                "pmfeature",
                affyHTExpressionPmFeatureSchema[["col2type"]],
                affyHTExpressionPmFeatureSchema[["col2key"]])
  
    if (hasMM)
      dbCreateTable(conn,
                    "mmfeature",
                    affyHTExpressionMmFeatureSchema[["col2type"]],
                    affyHTExpressionMmFeatureSchema[["col2key"]])
  
    dbInsertDataFrame(conn, "featureSet", parsedData[["featureSet"]],
                      affyHTExpressionFeatureSetSchema[["col2type"]], !quiet)
    dbInsertDataFrame(conn, "pmfeature", parsedData[["pmFeatures"]],
                      affyHTExpressionPmFeatureSchema[["col2type"]], !quiet)
    
    if (hasMM)
    dbInsertDataFrame(conn, "mmfeature", parsedData[["mmFeatures"]],
                      affyHTExpressionMmFeatureSchema[["col2type"]], !quiet)
    
    dbCreateTableInfo(conn, !quiet)
    
    ## Create indices
    dbCreateIndicesPm(conn, !quiet)
    dbCreateIndicesFs(conn, !quiet)
    
    dbGetQuery(conn, "VACUUM")
    dbDisconnect(conn)
    
    #######################################################################
    ## Part v) Save sequence DataFrames
    ## FIX ME: Fix ordering of the tables to match xxFeature tables
    #######################################################################
    datadir <- file.path(destDir, pkgName, "data")
    dir.create(datadir)
    pmSequence <- parsedData[["pmSequence"]]
    pmSeqFile <- file.path(datadir, "pmSequence.rda")
    if (!quiet) message("Saving DataFrame object for PM.\n")
    save(pmSequence, file=pmSeqFile, compress='xz')
    if (hasMM){
      mmSequence <- parsedData[["mmSequence"]]
      mmSeqFile <- file.path(datadir, "mmSequence.rda")
      if (!quiet) message("Saving DataFrame object for MM.\n")
      save(mmSequence, file=mmSeqFile, compress='xz')
    }
    if (!quiet) message("Done.\n")

}