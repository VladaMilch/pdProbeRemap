### original code from pdInfoBUilder package ###

checkFields <- function(needed, found){
    if (!all(needed %in% found))
    stop("Looking for fields '", 
        paste(needed, sep="", collapse="', '"), "', 
        but found '",
        paste(found, sep="", collapse="', '"), "'.")
    TRUE
}

#######################################################################
## SECTION A - Db Schema
#######################################################################
affyHTExpressionFeatureSetSchema <- list(col2type=c(
    fsetid="INTEGER",
    strand="INTEGER",
    man_fsetid="INTEGER"),
    col2key=c(
    fsetid="PRIMARY KEY"
    ))

affyHTExpressionPmFeatureSchema <- list(col2type=c(
    fid="INTEGER",
    fsetid="INTEGER",
    x="INTEGER",
    y="INTEGER"
),
col2key=c(
    fid="PRIMARY KEY"
))

affyHTExpressionMmFeatureSchema <- list(col2type=c(
    fid="INTEGER",
    fsetid="INTEGER",
    x="INTEGER",
    y="INTEGER",
    fidpm="INTEGER"
),
col2key=c(
    fid="PRIMARY KEY"
))

affyHTExpressionBgFeatureSchema <- list(col2type=c(
    fid="INTEGER",
    fsetid="INTEGER",
    x="INTEGER",
    y="INTEGER"
),
col2key=c(
    fid="PRIMARY KEY"
))


parseCdfCelProbe <- function(cdfFile, celFile, probeFile, verbose=TRUE){
    if (verbose) msgParsingFile(cdfFile)
    cdf <- readCdf(cdfFile)
    if (verbose) msgOK()
    
    if (verbose) msgParsingFile(celFile)
    cel <- readCelHeader(celFile)
    if (verbose) msgOK()
    geometry <- c(cel[["rows"]], cel[["cols"]])
    rm(cel)
    
    if (verbose) msgParsingFile(probeFile)
    cols <- c("probe.x", "probe.y", "probe.sequence")
    probeSeq <- read.delim(probeFile, stringsAsFactors=FALSE)
    names(probeSeq) <- tolower(names(probeSeq))
    ok <- checkFields(cols, names(probeSeq))
    probeSeq <- probeSeq[, cols]
    rm(cols, ok)
    names(probeSeq) <- c("x", "y", "sequence")
    if (verbose) msgOK()
    
    strands <- sapply(cdf, "[[", "unitdirection")
    strands <- ifelse(tolower(strands) == "sense",
    as.integer(SENSE),
    as.integer(ANTISENSE))
    if (verbose) simpleMessage("Getting information for featureSet table... ")
    featureSet <- data.frame(fsetid=1:length(strands),
         man_fsetid=names(strands),
         strand=strands,
         stringsAsFactors=FALSE)
    rm(strands)
    if (verbose) msgOK()
    
    extractFromGroups <- function(x){
    ## x is a list and has "groups" as component
    ngroups <- length(x[["groups"]])
    natoms <- 
        sapply(x[["groups"]], 
               "[[", "natoms") * sapply(x[["groups"]], "[[", "ncellsperatom")
    mfsetid <- unlist(mapply('rep', names(x[["groups"]]), each=natoms))
    mfsetid <- as.character(mfsetid)
    probes <- lapply(x[["groups"]],
     function(y){
         data.frame(x=y[["x"]],
        y=y[["y"]],
        isPm=!(y[["tbase"]]==y[["pbase"]]),
        atom=y[["atom"]])
     })
    probes <- do.call("rbind", probes)
    probes[["man_fsetid"]] <- mfsetid
    return(probes)
    }
    
    xy2i <- function(x, y, geom)
    as.integer(geom[1]*y+x+1)
    
    if (verbose) 
        simpleMessage("Getting information for pm/mm feature tables... ")
    allProbes <- lapply(cdf, extractFromGroups)
    allProbes <- do.call("rbind", allProbes)
    allProbes[["fid"]] <- xy2i(allProbes[["x"]], allProbes[["y"]], geometry)
    allProbes[["fsetid"]] <- featureSet[match(allProbes[["man_fsetid"]],
            featureSet[["man_fsetid"]]),
        "fsetid"]
    if (verbose) msgOK()
    if (verbose) simpleMessage("Combining probe information 
                                with sequence information... ")
    allProbes <- 
        merge(allProbes, probeSeq,
            by.x=c("x", "y"),
            by.y=c("x", "y"),
            all.x=TRUE)
    rm(probeSeq)
    if (verbose) msgOK()
    
    if (verbose) simpleMessage("Getting PM probes and sequences... ")
    geometry <- paste(geometry, collapse=";")
    cols <- c("fid", "fsetid", "x", "y", "atom")
    cols2 <- c("fid", "sequence")
    pmidx <- which(allProbes[["isPm"]])
    pmFeatures <- allProbes[pmidx, cols]
    pmSequence <- allProbes[pmidx, cols2]
    pmSequence <- pmSequence[order(pmSequence[["fid"]]),]
    if (verbose) msgOK()
    
    if (any(naseq <- is.na(pmSequence[["sequence"]])))
    warning("Probe sequences were not found for all PM probes. ",
    "These probes will be removed from the pmSequence object.")
    pmSequence <- pmSequence[!naseq,]
    
    pmSequence <- DataFrame(fid=pmSequence[["fid"]],
        sequence=DNAStringSet(pmSequence[["sequence"]]))
    
    mmFeatures <- allProbes[-pmidx, cols]
    mmSequence <- allProbes[-pmidx, cols2]
    if (any(naseq <- is.na(mmSequence[["sequence"]])))
    warning("Probe sequences were not found for all MM probes. ",
    "These probes will be removed from the mmSequence object.")
    mmSequence <- mmSequence[!naseq,]
    
    rm(pmidx, allProbes, naseq)
    
    cols1 <- c("fsetid", "atom")
    cols2 <- c("fsetid", "fid", "atom")
    matchpm <- merge(mmFeatures[, cols2], pmFeatures[, cols2],
     by.x=cols1, by.y=cols1)[, c("fid.x", "fid.y")]
    rm(cols1, cols2)
    names(matchpm) <- c("fid", "fidpm")
    
    mmFeatures <- merge(mmFeatures, matchpm, by.x="fid", by.y="fid")
    mmFeatures[["atom"]] <- NULL
    pmFeatures[["atom"]] <- NULL
    if (verbose) message("Done parsing.")
    return(list(featureSet=featureSet,
    pmSequence=pmSequence,
    pmFeatures=pmFeatures,
    mmSequence=mmSequence,
    mmFeatures=mmFeatures,
    geometry=geometry))
}

#######################################################################
## SECTION D - Package Maker
##     This shouldn't be extremely hard.
##     The idea is to: i) get array info (name, pkgname, dbname)
##     ii) parse data; iii) create pkg from template;
##     iv) dump the database
#######################################################################

setMethod("makePdInfoPackage", "AffyExpressionPDInfoPkgSeed",
    function(object, destDir=".", batch_size=10000, quiet=FALSE, unlink=FALSE) {
    
    msgBar()
    message("Building annotation package for Affymetrix Expression array")
    message("CDF...............: ", basename(object@cdfFile))
    message("CEL...............: ", basename(object@celFile))
    message("Sequence TAB-Delim: ", basename(object@tabSeqFile))
    msgBar()
    
    #######################################################################
    ## Part i) get array info (chipName, pkgName, dbname)
    #######################################################################
    chip <- chipName(object)
    pkgName <- cleanPlatformName(chip)
    extdataDir <- file.path(destDir, pkgName, "inst", "extdata")
    dbFileName <- paste(pkgName, "sqlite", sep=".")
    dbFilePath <- file.path(extdataDir, dbFileName)
    
    #######################################################################
    ## Part ii) parse data. This should return a list of data.frames.
    ##    The names of the elements in the list are table names.
    #######################################################################
    parsedData <- parseCdfCelProbe(object@cdfFile,
            object@celFile,
            object@tabSeqFile,
            verbose=!quiet)
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
        quiet=quiet)
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
    if (!quiet) message("Saving DataFrame object for PM.")
    save(pmSequence, file=pmSeqFile, compress='xz')
    if (hasMM){
    mmSequence <- parsedData[["mmSequence"]]
    mmSeqFile <- file.path(datadir, "mmSequence.rda")
    if (!quiet) message("Saving DataFrame object for MM.")
    save(mmSequence, file=mmSeqFile, compress='xz')
    }
    if (!quiet) message("Done.\n")
    })
