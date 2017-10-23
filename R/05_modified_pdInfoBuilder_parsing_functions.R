parseCdfCelProbe_andWriteAllProbesFile <- function(cdfFile, celFile, probeFile, verbose=TRUE, outputDir)
{
  if(!(dir.exists(outputDir)))
  {
    stop("Wrong path to Output directory! Input existing directory path")
  }
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
    natoms <- sapply(x[["groups"]], "[[", "natoms") * sapply(x[["groups"]], "[[", "ncellsperatom")
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
  
  if (verbose) simpleMessage("Getting information for pm/mm feature tables... ")
  allProbes <- lapply(cdf, extractFromGroups)
  allProbes <- do.call("rbind", allProbes)
  allProbes[["fid"]] <- xy2i(allProbes[["x"]], allProbes[["y"]], geometry)
  allProbes[["fsetid"]] <- featureSet[match(allProbes[["man_fsetid"]],
                                            featureSet[["man_fsetid"]]),
                                      "fsetid"]
  if (verbose) msgOK()
  if (verbose) simpleMessage("Combining probe information with sequence information... ")
  allProbes <- merge(allProbes, probeSeq,
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
  
  ### probes writing ####
  if(dir.exists(file.path(outputDir, "tmpProbesTableDir")))
  {
    warning(message = "The temporary directory exists?? Must not be the case on users machine!")
  }
  dir.create(file.path(outputDir, "tmpProbesTableDir"))
  
  write.csv(allProbes[!is.na(allProbes$sequence), ], 
              file = file.path(outputDir, "tmpProbesTableDir", "ProbesCDF.csv"), 
              row.names = FALSE, quote = FALSE)
  ### probes writing end ####
  
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

combinePgfClfProbesetsMps_andWriteAllProbesFile <- function(pgfFile, clfFile, probeFile,
                                      coreMps, fullMps, extendedMps,
                                      geneArray=FALSE, WT=TRUE,
                                      verbose=TRUE, outputDir)
{
  if(!(dir.exists(outputDir)))
  {
    stop("Wrong path to Output directory! Input existing directory path")
  }
  ## WT = Whole Transcript arrays (Gene/Exon ST)
  tmp <- parsePgfClf(pgfFile=pgfFile, clfFile=clfFile, verbose=verbose)
  probes.table <- tmp[["probes.table"]]
  geom <- tmp[["geometry"]]
  rm(tmp)
  
  ### probes writing ####
  if(dir.exists(file.path(outputDir, "tmpProbesTableDir")))
  {
    warning(message = "The temporary directory exists?? Must not be the case on users machine!")
  }
  dir.create(file.path(outputDir, "tmpProbesTableDir"))
  
  write.csv(probes.table, 
              file = file.path(outputDir, "tmpProbesTableDir", "ProbesPGF.csv"), 
              row.names = FALSE, quote = FALSE)
  ### probes writing end ####
  
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
      found <- paste('   ', sort(unique(featureSet$pstype)), collapse='\n')
      expct <- paste('   ', sort(unique(type_dict$type_id)), collapse='\n')
      txt <- paste('The type_dict template is incomplete.\nTemplate contains:\n', expct, '\n', 'Data contains:\n', found, sep='')
      stop(txt)
    }
    featureSet[['pstype']] <- NULL
  }
  
  ## pmfeature table - Fields
  ##  fid
  ##  fsetid
  ##  chr (NA)
  ##  location (NA)
  ##  x
  ##  y
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
  ##  fid
  ##  fid of matching pm
  ##  x
  ##  y
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
  
  ## IMPORTANT: for the moment, bgfeature will contain everything (that is PM) but 'main'
  ## bgfeature table - Fields
  ##  fid
  ##  x
  ##  y
  ##  fs_type: featureSet type: genomic/antigenomic
  ##  f_type: pm/mm at/st
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
    ##  featureSet: fsetid, type
    ##  pmfeature: fid, fsetid, atom, x, y
    ##  bgfeature: fid, fsetid, fs_type, f_type, x, y  - NOT ANYMORE
    ##  pmSequence: fid, sequence
    ##  bgSequence: fid, sequence  - NOT ANYMORE
    ##  core, extended, full: meta_fsetid, trancript_cluster_id, fsetid
    ##  mmfeatures/mmSequence
    
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