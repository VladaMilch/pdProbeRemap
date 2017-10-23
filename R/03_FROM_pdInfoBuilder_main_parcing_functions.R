#there are functions from the original package, just run this code as a source
#source("./Functions/pdInfoBuilderFROM/pdInfoBuilder_R_utils.R")

simpleMessage <- function(...)
    message(..., appendLF=FALSE)

#r

getLevelSchema <- function()
    data.frame(level=as.integer(1:5),
     level_id=c("core", "extended", "full", "free", "ambiguous"),
     stringsAsFactors=FALSE)

getTypeSchema <- function()
    data.frame(type=as.integer(1:17),
     type_id=c("main", "main->junctions", "main->psrs",
         "main->rescue", "control->affx", "control->chip",
         "control->bgp->antigenomic",
         "control->bgp->genomic", "normgene->exon",
         "normgene->intron", "rescue->FLmRNA->unmapped",
         "control->affx->bac_spike", "oligo_spike_in",
         "r1_bac_spike_at", "control->affx->polya_spike",
         "control->affx->ercc",
         "control->affx->ercc->step"),
     stringsAsFactors=FALSE)

msgParsingFile <- function(fname)
    simpleMessage("Parsing file: ", basename(fname))

msgOK <- function() message("OK\n")

msgBar <- function(){
    n <- options()[["width"]]
    message(paste(c(rep("=", n), "\n"), collapse=""))
}


## Define schema

#######################################################################
## SECTION A - Db Schema
#######################################################################
chromDictTable <- list(col2type=c(
    chrom="INTEGER",
    chrom_id="TEXT"),
    col2key=c(
    chrom="PRIMARY KEY"
    ))

levelDictTable <- list(col2type=c(
    level="INTEGER",
    level_id="TEXT"),
    col2key=c(
    level="PRIMARY KEY"
    ))

typeDictTable <- list(col2type=c(
    type="INTEGER",
    type_id="TEXT"),
    col2key=c(
    type="PRIMARY KEY"
    ))

exonTranscriptionFeatureSetSchema <- list(col2type=c(
    fsetid="INTEGER",
    strand="INTEGER",
    start="INTEGER",
    stop="INTEGER",
    transcript_cluster_id="INTEGER",
    exon_id="INTEGER",
    crosshyb_type="INTEGER",
    level="INTEGER",
    chrom="INTEGER",
    type="INTEGER"),
    col2key=c(
    fsetid="PRIMARY KEY",
    chrom="REFERENCES chrom_dict(chrom_id)",
    level="REFERENCES level_dict(level_id)",
    type ="REFERENCES type_dict(type_id)"
    ))

## merge both pmfeature schemes???
## apparently fid in gene arrays can't be key?
exonTranscriptionPmFeatureSchema <- list(col2type=c(
    fid="INTEGER",
    fsetid="INTEGER",
    atom="INTEGER",
    x="INTEGER",
    y="INTEGER"
),
col2key=c(
    fid="PRIMARY KEY"
))

genePmFeatureSchema <- list(col2type=c(
    fid="INTEGER",
    fsetid="INTEGER",
    atom="INTEGER",
    x="INTEGER",
    y="INTEGER"
),
col2key=c(
    fsetid="REFERENCES featureSet(fsetid)"
))


exonTranscriptionMmFeatureSchema <- list(col2type=c(
    fid="INTEGER",
    fsetid="INTEGER",
    atom="INTEGER",
    x="INTEGER",
    y="INTEGER"
),
col2key=c(
    fsetid="REFERENCES featureSet(fsetid)"
))

exonTranscriptionBgFeatureSchema <- list(col2type=c(
    fid="INTEGER",
    fsetid="INTEGER",
    fs_type="TEXT",
    f_type="TEXT",
    x="INTEGER",
    y="INTEGER"
),
col2key=c(
    fid="PRIMARY KEY"
))

geneBgFeatureSchema <- list(col2type=c(
    fid="INTEGER",
    fsetid="INTEGER",
    fs_type="TEXT",
    f_type="TEXT",
    x="INTEGER",
    y="INTEGER"
),
col2key=c(
    fsetid="REFERENCES featureSet(fsetid)"
))

mpsSchema <- list(col2type=c(
    meta_fsetid="INTEGER",
    transcript_cluster_id="INTEGER",
    fsetid="INTEGER"
),
col2key=c(
    fsetid="REFERENCES featureSet(fsetid)"
))

#######################################################################
## SECTION B - Utils - Already done.
#######################################################################

#######################################################################
## SECTION C - Parser specific for Affy Exon array
##     This will take PGF/CLF and process (in memory)
##     the data, generating data.frames for each table
##     (featureSet, pmfeature, mmfeature*, bgfeature**)
##     to be created in the db.
#######################################################################

.idxToIdx <- function(pgfFrom, pgfTo, ids) {
    starts <- pgfFrom[ids]
    ends <- pgfFrom[ids+1] - 1
    ends[is.na(ends)] <- length(pgfTo)
    mapply(":", starts, ends, SIMPLIFY=FALSE)
}

.combineIdx <- function(psToAtomIdx, atomToProbeIdx) {
    probesPerAtom <- with(atomToProbeIdx,
        sapply(split(atomIdx, atomIdx), length))
    cbind(probesetIdx=rep(psToAtomIdx$probesetIdx, probesPerAtom),
    atomToProbeIdx)
}

probesetIdxToAtomIdx <- function(pgf, probesetIdx) {
    atoms <- .idxToIdx(pgf[["probesetStartAtom"]],
         pgf[["atomId"]], probesetIdx)
    data.frame(probesetIdx=rep(probesetIdx, sapply(atoms, length)),
     atomIdx=unlist(atoms))
}

atomIdxToProbeIdx <- function(pgf, atomIdx) {
    probes <- .idxToIdx(pgf[["atomStartProbe"]],
        pgf[["probeId"]], atomIdx)
    data.frame(atomIdx=rep(atomIdx, sapply(probes, length)),
     probeIdx=unlist(probes))
}

probesetIdxToTripletIdx <- function(pgf, probesetIdx) {
    df1 <- probesetIdxToAtomIdx(pgf, probesetIdx)
    df2 <- atomIdxToProbeIdx(pgf, df1$atomIdx)
    .combineIdx(df1, df2)
}

parseProbesetCSV <- function(probeFile, verbose=TRUE){
    ## Variables
    SENSE <- as.integer(0)
    ANTISENSE <- as.integer(1)

    ###################################################
    ## TABLES TO ADD
    ###################################################
    if (verbose) simpleMessage("Creating dictionaries... ")
    ## chromosome dictionary moved to after reading
    ##    the CSV file

    ## level_schema
    level_schema <- getLevelSchema()
    ##     level_schema <- data.frame(level=as.integer(1:5),
    ##        level_id=c("core", "extended", 
    ## "full", "free", "ambiguous"),
    ##        stringsAsFactors=FALSE)
    ##

    ## type_schema
    type_schema <- getTypeSchema()
    ##     type_schema <- data.frame(type=as.integer(1:8),
    ##         type_id=c("main", "control->affx",
    ##             "control->chip",
    ##             "control->bgp->antigenomic",
    ##             "control->bgp->genomic",
    ##             "normgene->exon",
    ##             "normgene->intron",
    ##             "rescue->FLmRNA->unmapped"),
    ##         stringsAsFactors=FALSE)

    if (verbose) msgOK()


    ## the "probesets" df is to be the featureSet table
    if (verbose) msgParsingFile(probeFile)
    probesets <- utils::read.csv(probeFile, comment.char="#",
        stringsAsFactors=FALSE, na.strings="---")
    if (verbose) msgOK()
    cols <- c("probeset_id", "seqname", "strand", "start", "stop",
    "transcript_cluster_id", "exon_id",
    "crosshyb_type", "level", "probeset_type")
    probesets <- probesets[, cols]
    cols[1] <- "fsetid"
    names(probesets) <- cols
    rm(cols)

    chromosome_schema <- createChrDict(probesets[["seqname"]])

    probesets[["chrom"]] <- 
        match(probesets[["seqname"]], chromosome_schema[["chrom_id"]])
    probesets[["seqname"]] <- NULL
    probesets[["strand"]] <- 
        ifelse(probesets[["strand"]] == "-", ANTISENSE, SENSE)
    probesets[["level"]] <- 
        match(tolower(probesets[["level"]]), level_schema[["level_id"]])
    probesets[["type"]] <- 
        match(probesets[["probeset_type"]], type_schema[["type_id"]])
    probesets[["probeset_type"]] <- NULL

    list(
        probesets=probesets, 
        level=level_schema, 
        chromosome=chromosome_schema, 
        type=type_schema)
}

## Add tables: core, extended, full
mpsParser <- function(mpsFile, verbose=TRUE){
    if (verbose) msgParsingFile(mpsFile)
    mps <- 
        read.delim(
            mpsFile, 
            comment.char="#", 
            stringsAsFactors=FALSE, 
            header=TRUE)
    cols <- c("probeset_id", "transcript_cluster_id", "probeset_list")
    mps <- mps[cols]
    psids <- lapply(strsplit(mps[["probeset_list"]], " "), as.integer)
    psids <- cbind(mps[rep(1:nrow(mps),
         sapply(psids, length)), 1:2],
     item=unlist(psids))
    names(psids) <- c("meta_fsetid", "transcript_cluster_id", "fsetid")
    if (verbose) msgOK()
    return(psids)
}

parsePgfClf <- function(pgfFile, clfFile, verbose=TRUE){
    if (verbose) msgParsingFile(pgfFile)
    pgf <- readPgf(pgfFile)
    if (verbose) msgOK()
    if (verbose) msgParsingFile(clfFile)
    clf <- readClf(clfFile)
    if (verbose) msgOK()
    if (verbose) simpleMessage("Creating initial table for probes... ")
    geom <- paste(clf[["dims"]], collapse=";")
    triplet <- probesetIdxToTripletIdx(pgf, 1:length(pgf[["probesetId"]]))
    fid <- pgf[["probeId"]][triplet[["probeIdx"]]]
    i <- match(fid, pgf[["probeId"]])
    ii <- match(fid, clf[["id"]])
    probes.table <- data.frame(fid=fid,
         man_fsetid=pgf[["probesetName"]][triplet[["probesetIdx"]]],
         fsetid=pgf[["probesetId"]][triplet[["probesetIdx"]]],
         pstype=pgf[["probesetType"]][triplet[["probesetIdx"]]],
         atom=pgf[["atomId"]][triplet[["atomIdx"]]],
         x=clf[["x"]][ii],
         y=clf[["y"]][ii],
         ptype=pgf[["probeType"]][i],
         sequence=pgf[["probeSequence"]][i],
         stringsAsFactors=FALSE)
    rm(i, ii, triplet, fid, pgf, clf)
    if (verbose) msgOK()
    return(list(probes.table=probes.table, geometry=geom))
}

combinePgfClfProbesetsMps <- function(pgfFile, clfFile, probeFile,
            coreMps, fullMps, extendedMps,
            geneArray=FALSE, WT=TRUE,
            verbose=TRUE){
    ## WT = Whole Transcript arrays (Gene/Exon ST)
    tmp <- parsePgfClf(pgfFile=pgfFile, clfFile=clfFile, verbose=verbose)
    probes.table <- tmp[["probes.table"]]
    geom <- tmp[["geometry"]]
    rm(tmp)

    if (WT){
    probesetInfo <- parseProbesetCSV(probeFile, verbose=verbose)

    ## levels table
    ## id
    ## desc
    level_dict <- probesetInfo[["level"]]

    ## chromosome table
    ## id
    ## chrom_id
    chrom_dict <- probesetInfo[["chromosome"]]

    ## types table
    ## id
    ## type_id
    type_dict <- probesetInfo[["type"]]


    ## featureSet table - Fields
    ## probeset_id
    ## strand
    ## start
    ## stop
    ## transcript_cluster_id
    ## exon_id
    ## crosshyb_type
    ## level
    ## chrom
    ## type
    featureSet <- probesetInfo[["probesets"]]
    missFeatureSet <- setdiff(unique(probes.table[["fsetid"]]),
        unique(featureSet[["fsetid"]]))
    if (length(missFeatureSet) > 0){
    missFS = data.frame(fsetid=missFeatureSet)
    cols <- names(featureSet)
    cols <- cols[cols != "fsetid"]
    for (i in cols)
    missFS[[i]] <- NA
    missFS <- missFS[, names(featureSet)]
    featureSet <- rbind(featureSet, missFS)
    rm(missFS, cols, i)
    }
    rm(missFeatureSet)
    }else{
    featureSet <- unique(probes.table[, c('fsetid', 'man_fsetid', 'pstype')])
    type_dict <- getTypeSchema()
    featureSet[['type']] <- match(tolower(featureSet[['pstype']]),
            type_dict[['type_id']])
    if (any(is.na(featureSet[['type']]))){
    found <- paste('     ', sort(unique(featureSet$pstype)), collapse='\n')
    expct <- paste('     ', sort(unique(type_dict$type_id)), collapse='\n')
    txt <- 
        paste(
            'The type_dict template is incomplete.\nTemplate contains:\n', 
            expct, 
            '\n', 
            'Data contains:\n', 
            found, 
            sep='')
    stop(txt)
    }
    featureSet[['pstype']] <- NULL
    }

    ## pmfeature table - Fields
    ##    fid
    ##    fsetid
    ##    chr (NA)
    ##    location (NA)
    ##    x
    ##    y
    ## IMPORTANT:
    ##    ignoring strand
    ##    keeping atom to match with MM's
    pmFeatures <- subset(probes.table,
         substr(probes.table[["ptype"]], 1, 2) == "pm",
         select=c("fid", "fsetid", "atom", "x", "y", "sequence"))

    pmSequence <- pmFeatures[, c("fid", "sequence")]
    pmFeatures[["sequence"]] <- NULL
    pmSequence <- pmSequence[order(pmSequence[["fid"]]),]
    pmSequence <- DataFrame(fid=pmSequence[["fid"]],
        sequence=DNAStringSet(pmSequence[["sequence"]]))

    ## mmfeature table - Fields
    ##    fid
    ##    fid of matching pm
    ##    x
    ##    y
    ## IMPORTANT:
    ##    ignoring strand
    ##    keeping atom to match with MM's
    ##    ADD sequence for MM
    mmFeatures <- subset(probes.table, substr(probes.table$ptype, 1, 2) =="mm",
        select=c("fid", "fsetid", "atom", "x", "y", "sequence"))
    if (nrow(mmFeatures) > 0){
    mmSequence <- mmFeatures[, c("fid", "sequence")]
    mmFeatures[["sequence"]] <- NULL
    mmSequence <- mmSequence[order(mmSequence[["fid"]]),]
    mmSequence <- DataFrame(fid=mmSequence[["fid"]],
        sequence=DNAStringSet(mmSequence[["sequence"]]))
    }else{
    mmFeatures <- data.frame()
    mmSequence <- data.frame()
    }

    ## IMPORTANT: for the moment, bgfeature will 
    ##    contain everything (that is PM) but 'main'
    ## bgfeature table - Fields
    ##    fid
    ##    x
    ##    y
    ##    fs_type: featureSet type: genomic/antigenomic
    ##    f_type: pm/mm at/st
    ## old code:
    ## subset using cols
    ## cols <- c("fid", "fsetid", "pstype", "ptype", "x", "y", "sequence")
    rm(probes.table)

    if (WT){
    core <- mpsParser(coreMps, verbose=verbose)
    if (!geneArray){
    extended <- mpsParser(extendedMps, verbose=verbose)
    full <- mpsParser(fullMps, verbose=verbose)
    }

    ## Here we should have the following tables available:
    ##    featureSet: fsetid, type
    ##    pmfeature: fid, fsetid, atom, x, y
    ##    bgfeature: fid, fsetid, fs_type, f_type, x, y    - NOT ANYMORE
    ##    pmSequence: fid, sequence
    ##    bgSequence: fid, sequence    - NOT ANYMORE
    ##    core, extended, full: meta_fsetid, trancript_cluster_id, fsetid
    ##    mmfeatures/mmSequence

    out <- list(featureSet=featureSet, pmFeatures=pmFeatures,
    mmFeatures=mmFeatures, geometry=geom,
    pmSequence=pmSequence, mmSequence=mmSequence,
    chrom_dict=chrom_dict, level_dict=level_dict,
    type_dict=type_dict, core=core)
    if (!geneArray){
    out[["extended"]] <- extended
    out[["full"]] <- full
    }
    }else{
    out <- list(featureSet=featureSet, pmFeatures=pmFeatures,
    mmFeatures=mmFeatures, geometry=geom,
    pmSequence=pmSequence, mmSequence=mmSequence,
    type_dict=type_dict)
    }
    return(out)
}

combinePgfClfProbesetsMps0 <- function(pgfFile, clfFile, probeFile,
             coreMps, fullMps, extendedMps,
             geneArray=FALSE, verbose=TRUE){
    tmp <- parsePgfClf(pgfFile=pgfFile, clfFile=clfFile, verbose=verbose)
    probes.table <- tmp[["probes.table"]]
    geom <- tmp[["geometry"]]
    rm(tmp)
    probesetInfo <- parseProbesetCSV(probeFile, verbose=verbose)

    ## levels table
    ## id
    ## desc
    level_dict <- probesetInfo[["level"]]

    ## chromosome table
    ## id
    ## chrom_id
    chrom_dict <- probesetInfo[["chromosome"]]

    ## types table
    ## id
    ## type_id
    type_dict <- probesetInfo[["type"]]


    ## featureSet table - Fields
    ## probeset_id
    ## strand
    ## start
    ## stop
    ## transcript_cluster_id
    ## exon_id
    ## crosshyb_type
    ## level
    ## chrom
    ## type
    featureSet <- probesetInfo[["probesets"]]
    missFeatureSet <- setdiff(unique(probes.table[["fsetid"]]),
        unique(featureSet[["fsetid"]]))
    if (length(missFeatureSet) > 0){
    missFS = data.frame(fsetid=missFeatureSet)
    cols <- names(featureSet)
    cols <- cols[cols != "fsetid"]
    for (i in cols)
    missFS[[i]] <- NA
    missFS <- missFS[, names(featureSet)]
    featureSet <- rbind(featureSet, missFS)
    rm(missFS, cols, i)
    }
    rm(missFeatureSet)

    ## pmfeature table - Fields
    ##    fid
    ##    fsetid
    ##    chr (NA)
    ##    location (NA)
    ##    x
    ##    y
    ## IMPORTANT:
    ##    ignoring strand
    ##    keeping atom to match with MM's
    pmFeatures <- subset(probes.table,
         substr(probes.table[["ptype"]], 1, 2) == "pm",
         select=c("fid", "fsetid", "atom", "x", "y", "sequence"))

    pmSequence <- pmFeatures[, c("fid", "sequence")]
    pmFeatures[["sequence"]] <- NULL
    pmSequence <- pmSequence[order(pmSequence[["fid"]]),]
    pmSequence <- DataFrame(fid=pmSequence[["fid"]],
        sequence=DNAStringSet(pmSequence[["sequence"]]))

    ## mmfeature table - Fields
    ##    fid
    ##    fid of matching pm
    ##    x
    ##    y
    ## IMPORTANT:
    ##    ignoring strand
    ##    keeping atom to match with MM's
    ##    ADD sequence for MM
    mmFeatures <- subset(probes.table, substr(probes.table$ptype, 1, 2) =="mm",
         select=c("fid", "fsetid", "atom", "x", "y", "sequence"))
    if (nrow(mmFeatures) > 0){
    mmSequence <- mmFeatures[, c("fid", "sequence")]
    mmFeatures[["sequence"]] <- NULL
    mmSequence <- mmSequence[order(mmSequence[["fid"]]),]
    mmSequence <- DataFrame(fid=mmSequence[["fid"]],
        sequence=DNAStringSet(mmSequence[["sequence"]]))
    }else{
    mmFeatures <- data.frame()
    mmSequence <- data.frame()
    }

    ## IMPORTANT: for the moment, 
    ##    bgfeature will contain everything (that is PM) but 'main'
    ## bgfeature table - Fields
    ##    fid
    ##    x
    ##    y
    ##    fs_type: featureSet type: genomic/antigenomic
    ##    f_type: pm/mm at/st
    ## old code:
    ## subset using cols
    ## cols <- c("fid", "fsetid", "pstype", "ptype", "x", "y", "sequence")
    rm(probes.table)

    core <- mpsParser(coreMps, verbose=verbose)
    if (!geneArray){
    extended <- mpsParser(extendedMps, verbose=verbose)
    full <- mpsParser(fullMps, verbose=verbose)
    }

    ## Here we should have the following tables available:
    ##    featureSet: fsetid, type
    ##    pmfeature: fid, fsetid, atom, x, y
    ##    bgfeature: fid, fsetid, fs_type, f_type, x, y    - NOT ANYMORE
    ##    pmSequence: fid, sequence
    ##    bgSequence: fid, sequence    - NOT ANYMORE
    ##    core, extended, full: meta_fsetid, trancript_cluster_id, fsetid
    ##    mmfeatures/mmSequence

    out <- list(featureSet=featureSet, pmFeatures=pmFeatures,
    mmFeatures=mmFeatures, geometry=geom,
    pmSequence=pmSequence, mmSequence=mmSequence,
    chrom_dict=chrom_dict, level_dict=level_dict,
    type_dict=type_dict, core=core)
    if (!geneArray){
    out[["extended"]] <- extended
    out[["full"]] <- full
    }

    return(out)
}



#Define something else for those functions
#######################################################################
################## ORIGINAL MAKE PDINFO PACKAGE METHOD/FUNCTION SCRIPT (PART 1)

## Define schema

#######################################################################
## SECTION A - Db Schema
#######################################################################
chromDictTable <- list(col2type=c(
    chrom="INTEGER",
    chrom_id="TEXT"),
    col2key=c(
    chrom="PRIMARY KEY"
    ))

levelDictTable <- list(col2type=c(
    level="INTEGER",
    level_id="TEXT"),
    col2key=c(
    level="PRIMARY KEY"
    ))

typeDictTable <- list(col2type=c(
    type="INTEGER",
    type_id="TEXT"),
    col2key=c(
    type="PRIMARY KEY"
    ))

exonTranscriptionFeatureSetSchema <- list(col2type=c(
    fsetid="INTEGER",
    strand="INTEGER",
    start="INTEGER",
    stop="INTEGER",
    transcript_cluster_id="INTEGER",
    exon_id="INTEGER",
    crosshyb_type="INTEGER",
    level="INTEGER",
    chrom="INTEGER",
    type="INTEGER"),
    col2key=c(
    fsetid="PRIMARY KEY",
    chrom="REFERENCES chrom_dict(chrom_id)",
    level="REFERENCES level_dict(level_id)",
    type ="REFERENCES type_dict(type_id)"
    ))

## merge both pmfeature schemes???
## apparently fid in gene arrays can't be key?
exonTranscriptionPmFeatureSchema <- list(col2type=c(
    fid="INTEGER",
    fsetid="INTEGER",
    atom="INTEGER",
    x="INTEGER",
    y="INTEGER"
),
col2key=c(
    fid="PRIMARY KEY"
))

genePmFeatureSchema <- list(col2type=c(
    fid="INTEGER",
    fsetid="INTEGER",
    atom="INTEGER",
    x="INTEGER",
    y="INTEGER"
),
col2key=c(
    fsetid="REFERENCES featureSet(fsetid)"
))


exonTranscriptionMmFeatureSchema <- list(col2type=c(
    fid="INTEGER",
    fsetid="INTEGER",
    atom="INTEGER",
    x="INTEGER",
    y="INTEGER"
),
col2key=c(
    fsetid="REFERENCES featureSet(fsetid)"
))

exonTranscriptionBgFeatureSchema <- list(col2type=c(
    fid="INTEGER",
    fsetid="INTEGER",
    fs_type="TEXT",
    f_type="TEXT",
    x="INTEGER",
    y="INTEGER"
),
col2key=c(
    fid="PRIMARY KEY"
))

geneBgFeatureSchema <- list(col2type=c(
    fid="INTEGER",
    fsetid="INTEGER",
    fs_type="TEXT",
    f_type="TEXT",
    x="INTEGER",
    y="INTEGER"
),
col2key=c(
    fsetid="REFERENCES featureSet(fsetid)"
))

mpsSchema <- list(col2type=c(
    meta_fsetid="INTEGER",
    transcript_cluster_id="INTEGER",
    fsetid="INTEGER"
),
col2key=c(
    fsetid="REFERENCES featureSet(fsetid)"
))

#######################################################################
## SECTION B - Utils - Already done.
#######################################################################

#######################################################################
## SECTION C - Parser specific for Affy Exon array
##     This will take PGF/CLF and process (in memory)
##     the data, generating data.frames for each table
##     (featureSet, pmfeature, mmfeature*, bgfeature**)
##     to be created in the db.
#######################################################################

###########

.idxToIdx <- function(pgfFrom, pgfTo, ids) {
    starts <- pgfFrom[ids]
    ends <- pgfFrom[ids+1] - 1
    ends[is.na(ends)] <- length(pgfTo)
    mapply(":", starts, ends, SIMPLIFY=FALSE)
}

.combineIdx <- function(psToAtomIdx, atomToProbeIdx) {
    probesPerAtom <- with(atomToProbeIdx,
        sapply(split(atomIdx, atomIdx), length))
    cbind(probesetIdx=rep(psToAtomIdx$probesetIdx, probesPerAtom),
    atomToProbeIdx)
}

probesetIdxToAtomIdx <- function(pgf, probesetIdx) {
    atoms <- .idxToIdx(pgf[["probesetStartAtom"]],
         pgf[["atomId"]], probesetIdx)
    data.frame(probesetIdx=rep(probesetIdx, sapply(atoms, length)),
     atomIdx=unlist(atoms))
}

atomIdxToProbeIdx <- function(pgf, atomIdx) {
    probes <- .idxToIdx(pgf[["atomStartProbe"]],
        pgf[["probeId"]], atomIdx)
    data.frame(atomIdx=rep(atomIdx, sapply(probes, length)),
     probeIdx=unlist(probes))
}

probesetIdxToTripletIdx <- function(pgf, probesetIdx) {
    df1 <- probesetIdxToAtomIdx(pgf, probesetIdx)
    df2 <- atomIdxToProbeIdx(pgf, df1$atomIdx)
    .combineIdx(df1, df2)
}

parseProbesetCSV <- function(probeFile, verbose=TRUE){
    ## Variables
    SENSE <- as.integer(0)
    ANTISENSE <- as.integer(1)

    ###################################################
    ## TABLES TO ADD
    ###################################################
    if (verbose) simpleMessage("Creating dictionaries... ")
    ## chromosome dictionary moved to after reading
    ##    the CSV file

    ## level_schema
    level_schema <- getLevelSchema()
    ##     level_schema <- data.frame(level=as.integer(1:5),
    ##        level_id=c("core", "extended", "full", "free", "ambiguous"),
    ##        stringsAsFactors=FALSE)
    ##

    ## type_schema
    type_schema <- getTypeSchema()
    ##     type_schema <- data.frame(type=as.integer(1:8),
    ##         type_id=c("main", "control->affx",
    ##             "control->chip",
    ##             "control->bgp->antigenomic",
    ##             "control->bgp->genomic",
    ##             "normgene->exon",
    ##             "normgene->intron",
    ##             "rescue->FLmRNA->unmapped"),
    ##         stringsAsFactors=FALSE)

    if (verbose) msgOK()


    ## the "probesets" df is to be the featureSet table
    if (verbose) msgParsingFile(probeFile)
    probesets <- utils::read.csv(probeFile, comment.char="#",
        stringsAsFactors=FALSE, na.strings="---")
    if (verbose) msgOK()
    cols <- c("probeset_id", "seqname", "strand", "start", "stop",
    "transcript_cluster_id", "exon_id",
    "crosshyb_type", "level", "probeset_type")
    probesets <- probesets[, cols]
    cols[1] <- "fsetid"
    names(probesets) <- cols
    rm(cols)

    chromosome_schema <- createChrDict(probesets[["seqname"]])

    probesets[["chrom"]] <- 
        match(probesets[["seqname"]], chromosome_schema[["chrom_id"]])
    probesets[["seqname"]] <- NULL
    probesets[["strand"]] <- 
        ifelse(probesets[["strand"]] == "-", ANTISENSE, SENSE)
    probesets[["level"]] <- 
        match(tolower(probesets[["level"]]), level_schema[["level_id"]])
    probesets[["type"]] <- 
        match(probesets[["probeset_type"]], type_schema[["type_id"]])
    probesets[["probeset_type"]] <- NULL

    list(   probesets=probesets, 
            level=level_schema, 
            chromosome=chromosome_schema, 
            type=type_schema)
}

## Add tables: core, extended, full
mpsParser <- function(mpsFile, verbose=TRUE){
    if (verbose) msgParsingFile(mpsFile)
    mps <- 
        read.delim(
            mpsFile, comment.char="#", stringsAsFactors=FALSE, header=TRUE)
    cols <- c("probeset_id", "transcript_cluster_id", "probeset_list")
    mps <- mps[cols]
    psids <- lapply(strsplit(mps[["probeset_list"]], " "), as.integer)
    psids <- 
        cbind(mps[rep(1:nrow(mps),sapply(psids, length)), 1:2],
        item=unlist(psids))
    names(psids) <- c("meta_fsetid", "transcript_cluster_id", "fsetid")
    if (verbose) msgOK()
    return(psids)
}

parsePgfClf <- function(pgfFile, clfFile, verbose=TRUE){
    if (verbose) msgParsingFile(pgfFile)
    pgf <- readPgf(pgfFile)
    if (verbose) msgOK()
    if (verbose) msgParsingFile(clfFile)
    clf <- readClf(clfFile)
    if (verbose) msgOK()
    if (verbose) simpleMessage("Creating initial table for probes... ")
    geom <- paste(clf[["dims"]], collapse=";")
    triplet <- probesetIdxToTripletIdx(pgf, 1:length(pgf[["probesetId"]]))
    fid <- pgf[["probeId"]][triplet[["probeIdx"]]]
    i <- match(fid, pgf[["probeId"]])
    ii <- match(fid, clf[["id"]])
    probes.table <- 
        data.frame(fid=fid,
            man_fsetid=pgf[["probesetName"]][triplet[["probesetIdx"]]],
            fsetid=pgf[["probesetId"]][triplet[["probesetIdx"]]],
            pstype=pgf[["probesetType"]][triplet[["probesetIdx"]]],
            atom=pgf[["atomId"]][triplet[["atomIdx"]]],
            x=clf[["x"]][ii],
            y=clf[["y"]][ii],
            ptype=pgf[["probeType"]][i],
            sequence=pgf[["probeSequence"]][i],
            stringsAsFactors=FALSE)
    rm(i, ii, triplet, fid, pgf, clf)
    if (verbose) msgOK()
    return(list(probes.table=probes.table, geometry=geom))
}

combinePgfClfProbesetsMps <- 
    function(pgfFile, clfFile, probeFile,
                coreMps, fullMps, extendedMps,
                geneArray=FALSE, WT=TRUE,
                verbose=TRUE)
{
    ## WT = Whole Transcript arrays (Gene/Exon ST)
    tmp <- parsePgfClf(pgfFile=pgfFile, clfFile=clfFile, verbose=verbose)
    probes.table <- tmp[["probes.table"]]
    geom <- tmp[["geometry"]]
    rm(tmp)

    if (WT){
    probesetInfo <- parseProbesetCSV(probeFile, verbose=verbose)

    ## levels table
    ## id
    ## desc
    level_dict <- probesetInfo[["level"]]

    ## chromosome table
    ## id
    ## chrom_id
    chrom_dict <- probesetInfo[["chromosome"]]

    ## types table
    ## id
    ## type_id
    type_dict <- probesetInfo[["type"]]


    ## featureSet table - Fields
    ## probeset_id
    ## strand
    ## start
    ## stop
    ## transcript_cluster_id
    ## exon_id
    ## crosshyb_type
    ## level
    ## chrom
    ## type
    featureSet <- probesetInfo[["probesets"]]
    missFeatureSet <- setdiff(unique(probes.table[["fsetid"]]),
        unique(featureSet[["fsetid"]]))
    if (length(missFeatureSet) > 0){
    missFS = data.frame(fsetid=missFeatureSet)
    cols <- names(featureSet)
    cols <- cols[cols != "fsetid"]
    for (i in cols)
    missFS[[i]] <- NA
    missFS <- missFS[, names(featureSet)]
    featureSet <- rbind(featureSet, missFS)
    rm(missFS, cols, i)
    }
    rm(missFeatureSet)
    }else{
    featureSet <- unique(probes.table[, c('fsetid', 'man_fsetid', 'pstype')])
    type_dict <- getTypeSchema()
    featureSet[['type']] <- 
        match(tolower(featureSet[['pstype']]),
                type_dict[['type_id']])
    if (any(is.na(featureSet[['type']]))){
    found <- paste('     ', sort(unique(featureSet$pstype)), collapse='\n')
    expct <- paste('     ', sort(unique(type_dict$type_id)), collapse='\n')
    txt <- 
        paste('The type_dict template is incomplete.\nTemplate contains:\n', 
                expct, '\n', 'Data contains:\n', found, sep='')
    stop(txt)
    }
    featureSet[['pstype']] <- NULL
    }

    ## pmfeature table - Fields
    ##    fid
    ##    fsetid
    ##    chr (NA)
    ##    location (NA)
    ##    x
    ##    y
    ## IMPORTANT:
    ##    ignoring strand
    ##    keeping atom to match with MM's
    pmFeatures <- subset(probes.table,
            substr(probes.table[["ptype"]], 1, 2) == "pm",
            select=c("fid", "fsetid", "atom", "x", "y", "sequence"))

    pmSequence <- pmFeatures[, c("fid", "sequence")]
    pmFeatures[["sequence"]] <- NULL
    pmSequence <- pmSequence[order(pmSequence[["fid"]]),]
    pmSequence <- DataFrame(fid=pmSequence[["fid"]],
        sequence=DNAStringSet(pmSequence[["sequence"]]))

    ## mmfeature table - Fields
    ##    fid
    ##    fid of matching pm
    ##    x
    ##    y
    ## IMPORTANT:
    ##    ignoring strand
    ##    keeping atom to match with MM's
    ##    ADD sequence for MM
    mmFeatures <- subset(probes.table, substr(probes.table$ptype, 1, 2) =="mm",
            select=c("fid", "fsetid", "atom", "x", "y", "sequence"))
    if (nrow(mmFeatures) > 0){
    mmSequence <- mmFeatures[, c("fid", "sequence")]
    mmFeatures[["sequence"]] <- NULL
    mmSequence <- mmSequence[order(mmSequence[["fid"]]),]
    mmSequence <- DataFrame(fid=mmSequence[["fid"]],
        sequence=DNAStringSet(mmSequence[["sequence"]]))
    }else{
    mmFeatures <- data.frame()
    mmSequence <- data.frame()
    }

    ## IMPORTANT: for the moment, bgfeature will contain everything 
    ## (that is PM) but 'main'
    ## bgfeature table - Fields
    ##    fid
    ##    x
    ##    y
    ##    fs_type: featureSet type: genomic/antigenomic
    ##    f_type: pm/mm at/st
    ## old code:
    ## subset using cols
    ## cols <- c("fid", "fsetid", "pstype", "ptype", "x", "y", "sequence")
    rm(probes.table)

    if (WT){
    core <- mpsParser(coreMps, verbose=verbose)
    if (!geneArray){
    extended <- mpsParser(extendedMps, verbose=verbose)
    full <- mpsParser(fullMps, verbose=verbose)
    }

    ## Here we should have the following tables available:
    ##    featureSet: fsetid, type
    ##    pmfeature: fid, fsetid, atom, x, y
    ##    bgfeature: fid, fsetid, fs_type, f_type, x, y    - NOT ANYMORE
    ##    pmSequence: fid, sequence
    ##    bgSequence: fid, sequence    - NOT ANYMORE
    ##    core, extended, full: meta_fsetid, trancript_cluster_id, fsetid
    ##    mmfeatures/mmSequence

    out <- list(featureSet=featureSet, pmFeatures=pmFeatures,
    mmFeatures=mmFeatures, geometry=geom,
    pmSequence=pmSequence, mmSequence=mmSequence,
    chrom_dict=chrom_dict, level_dict=level_dict,
    type_dict=type_dict, core=core)
    if (!geneArray){
    out[["extended"]] <- extended
    out[["full"]] <- full
    }
    }else{
    out <- list(featureSet=featureSet, pmFeatures=pmFeatures,
    mmFeatures=mmFeatures, geometry=geom,
    pmSequence=pmSequence, mmSequence=mmSequence,
    type_dict=type_dict)
    }
    return(out)
}

combinePgfClfProbesetsMps0 <- function(pgfFile, clfFile, probeFile,
            coreMps, fullMps, extendedMps,
            geneArray=FALSE, verbose=TRUE){
    tmp <- parsePgfClf(pgfFile=pgfFile, clfFile=clfFile, verbose=verbose)
    probes.table <- tmp[["probes.table"]]
    geom <- tmp[["geometry"]]
    rm(tmp)
    probesetInfo <- parseProbesetCSV(probeFile, verbose=verbose)

    ## levels table
    ## id
    ## desc
    level_dict <- probesetInfo[["level"]]

    ## chromosome table
    ## id
    ## chrom_id
    chrom_dict <- probesetInfo[["chromosome"]]

    ## types table
    ## id
    ## type_id
    type_dict <- probesetInfo[["type"]]


    ## featureSet table - Fields
    ## probeset_id
    ## strand
    ## start
    ## stop
    ## transcript_cluster_id
    ## exon_id
    ## crosshyb_type
    ## level
    ## chrom
    ## type
    featureSet <- probesetInfo[["probesets"]]
    missFeatureSet <- setdiff(unique(probes.table[["fsetid"]]),
        unique(featureSet[["fsetid"]]))
    if (length(missFeatureSet) > 0){
    missFS = data.frame(fsetid=missFeatureSet)
    cols <- names(featureSet)
    cols <- cols[cols != "fsetid"]
    for (i in cols)
    missFS[[i]] <- NA
    missFS <- missFS[, names(featureSet)]
    featureSet <- rbind(featureSet, missFS)
    rm(missFS, cols, i)
    }
    rm(missFeatureSet)

    ## pmfeature table - Fields
    ##    fid
    ##    fsetid
    ##    chr (NA)
    ##    location (NA)
    ##    x
    ##    y
    ## IMPORTANT:
    ##    ignoring strand
    ##    keeping atom to match with MM's
    pmFeatures <- subset(probes.table,
        substr(probes.table[["ptype"]], 1, 2) == "pm",
        select=c("fid", "fsetid", "atom", "x", "y", "sequence"))

    pmSequence <- pmFeatures[, c("fid", "sequence")]
    pmFeatures[["sequence"]] <- NULL
    pmSequence <- pmSequence[order(pmSequence[["fid"]]),]
    pmSequence <- DataFrame(fid=pmSequence[["fid"]],
        sequence=DNAStringSet(pmSequence[["sequence"]]))

    ## mmfeature table - Fields
    ##    fid
    ##    fid of matching pm
    ##    x
    ##    y
    ## IMPORTANT:
    ##    ignoring strand
    ##    keeping atom to match with MM's
    ##    ADD sequence for MM
    mmFeatures <- subset(probes.table, substr(probes.table$ptype, 1, 2) =="mm",
        select=c("fid", "fsetid", "atom", "x", "y", "sequence"))
    if (nrow(mmFeatures) > 0){
    mmSequence <- mmFeatures[, c("fid", "sequence")]
    mmFeatures[["sequence"]] <- NULL
    mmSequence <- mmSequence[order(mmSequence[["fid"]]),]
    mmSequence <- DataFrame(fid=mmSequence[["fid"]],
        sequence=DNAStringSet(mmSequence[["sequence"]]))
    }else{
    mmFeatures <- data.frame()
    mmSequence <- data.frame()
    }

    ## IMPORTANT: for the moment, bgfeature will contain everything 
    ## (that is PM) but 'main'
    ## bgfeature table - Fields
    ##    fid
    ##    x
    ##    y
    ##    fs_type: featureSet type: genomic/antigenomic
    ##    f_type: pm/mm at/st
    ## old code:
    ## subset using cols
    ## cols <- c("fid", "fsetid", "pstype", "ptype", "x", "y", "sequence")
    rm(probes.table)

    core <- mpsParser(coreMps, verbose=verbose)
    if (!geneArray){
    extended <- mpsParser(extendedMps, verbose=verbose)
    full <- mpsParser(fullMps, verbose=verbose)
    }

    ## Here we should have the following tables available:
    ##    featureSet: fsetid, type
    ##    pmfeature: fid, fsetid, atom, x, y
    ##    bgfeature: fid, fsetid, fs_type, f_type, x, y    - NOT ANYMORE
    ##    pmSequence: fid, sequence
    ##    bgSequence: fid, sequence    - NOT ANYMORE
    ##    core, extended, full: meta_fsetid, trancript_cluster_id, fsetid
    ##    mmfeatures/mmSequence

    out <- list(featureSet=featureSet, pmFeatures=pmFeatures,
    mmFeatures=mmFeatures, geometry=geom,
    pmSequence=pmSequence, mmSequence=mmSequence,
    chrom_dict=chrom_dict, level_dict=level_dict,
    type_dict=type_dict, core=core)
    if (!geneArray){
    out[["extended"]] <- extended
    out[["full"]] <- full
    }

    return(out)
}
#########
# all commented stuff from original function
#######################################################################
## SECTION D - Package Maker
##     This shouldn't be extremely hard.
##     The idea is to: i) get array info (name, pkgname, dbname)
##     ii) parse data; iii) create pkg from template;
##     iv) dump the database
#######################################################################

# setMethod("makePdInfoPackage", "AffySTPDInfoPkgSeed",
#     function(object, destDir=".", batch_size=10000, 
#     quiet=FALSE, unlink=FALSE) {
#     geneArray <- object@geneArray
#     stopifnot(is.logical(geneArray))
#     if (geneArray){
#     msg <- "Building annotation package for Affymetrix Gene ST Array"
#     }else{
#     msg <- "Building annotation package for Affymetrix Exon ST Array"
#     }
#
#     msgBar()
#     message(msg)
#     message("PGF.........: ", basename(object@pgfFile))
#     message("CLF.........: ", basename(object@clfFile))
#     message("Probeset....: ", basename(object@probeFile))
#     message("Transcript..: ", basename(object@transFile))
#     message("Core MPS....: ", basename(object@coreMps))
#     if (!geneArray){
#     message("Full MPS....: ", basename(object@fullMps))
#     message("Extended MPS: ", basename(object@extendedMps))
#     }
#     msgBar()
#
#     #######################################################################
#     ## Part i) get array info (chipName, pkgName, dbname)
#     #######################################################################
#     chip <- chipName(object)
#     pkgName <- cleanPlatformName(chip)
#     extdataDir <- file.path(destDir, pkgName, "inst", "extdata")
#     dbFileName <- paste(pkgName, "sqlite", sep=".")
#     dbFilePath <- file.path(extdataDir, dbFileName)
#
#     #######################################################################
#     ## Part ii) parse data. This should return a list of data.frames.
#     ##    The names of the elements in the list are table names.
#     #######################################################################
#     parsedData <- combinePgfClfProbesetsMps(object@pgfFile,
#                 object@clfFile,
#                 object@probeFile,
#                 object@coreMps,
#                 object@fullMps,
#                 object@extendedMps,
#                 verbose=!quiet,
#                 geneArray=geneArray)
############
