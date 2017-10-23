#' build_new_ParsedDataOld_object_PGF, build_new_ParsedDataOld_object_CDF
#'
#' 
#' @title build new ParsedDataOld object PGF or CDF
#' @return list with the same structure as original ParsedData object, containing NewParsedData
#' @author Vladislava Milchevskaya \email{milchv@gmail.com}
#' @importFrom dplyr arrange
#' @param level "gene" or "transcript"
#' @param ParsedDataOld old arced data 
#' @param PS_controls ids of control probe sets
#' @param Annotation reference genome annotation data.frame
#' @param ProbeSetData an object with new probe groupings
#' @param outputDir output directory, the one that was used to store Fasta files

# needed ofjects:
#   featureSet
#   pmFeatures
#   mmFeatures
#   geometry
#   pmSequence
#   mmSequence
#   chrom_dict
#   level_dict
#   type_dict
#   core

build_new_ParsedDataOld_object_PGF <- function(level,
                                               ParsedDataOld, 
                                               PS_controls, 
                                               Annotation,
                                               ProbeSetData, 
                                               outputDir 
                                               )
{
  options(stringsAsFactors = FALSE)
  #######################################
  ############## newLeveDict ############
  newLeveDict <- ParsedDataOld$level_dict
  ############## newTypeDict ############
  newTypeDict <- ParsedDataOld$type_dict
  ############## newGeometry ############
  newGeometry <- ParsedDataOld$geometry
  
  # bugs come from this package that converts an annotation to a data frame
  if("" %in% colnames(Annotation)) 
  {
    Annotation = Annotation[, setdiff(colnames(Annotation), "")]
  }
  
  PS_controls = subset(PS_controls, 
                       !(fid %in% ProbeSetData$ProbeSetDataFrame$fid))
  
  ########################################
  ############## newChromDict ############
  # from old feature set find out chrom ids for control probes # they have none!
  type_id_main <- 
    ParsedDataOld$type_dict[which(ParsedDataOld$type_dict$type_id == "main"), 
                            "type"]
  nonMain_from_feature_set <- 
    subset(
      subset(ParsedDataOld$featureSet, 
             type != type_id_main | is.na(type)), 
      fsetid %in% PS_controls$fsetid)
  unique_chrom <- unique(ProbeSetData$ProbeSetDataAnnotation$chrom)
  num_chrom <- suppressWarnings(as.numeric(unique_chrom))
  num_chrom[which(is.na(num_chrom))] <- rep(100000, 
                                            length(which(is.na(num_chrom))))
  unique_chrom[order(num_chrom)]
  
  chrom_df <- data.frame(unique_chrom, num_chrom)
  chrom_df = dplyr::arrange(chrom_df, 
                            num_chrom, 
                            nchar(unique_chrom), 
                            unique_chrom)
  newChromDict = data.frame(chrom = c(1:length(chrom_df$unique_chrom)), 
                            chrom_id = chrom_df$unique_chrom)
  
  ##########################################
  ############   newFeatureSet   ###########
  Annotation_tr = subset(Annotation[grep("RNA|Rna|rna|transcript" ,
                                         Annotation$feature), ], 
                         strand %in% c("+", "-"))

  ProbeSetData$ProbeSetDataAnnotation$strand <- 
    ifelse((ProbeSetData$ProbeSetDataAnnotation$strand == "-" | 
              ProbeSetData$ProbeSetDataAnnotation$strand == 1), 1, 0)
  mm <- merge(newChromDict, ProbeSetData$ProbeSetDataAnnotation, 
              by.x = "chrom_id", 
              by.y = "chrom")
  mmm <- data.frame(mm, 
                    type = rep(type_id_main, nrow(mm)),
                    level = rep(NA, nrow(mm)), 
                    exon_id = rep(0, nrow(mm)), 
                    crosshyb_type = rep(1, nrow(mm)), 
                    transcript_cluster_id = mm$fsetid
  )

  ############## get start stop coordinates for FeatureSet table ############
  if( level == "transcript")
  {
    start_stop_coordinates_per_transcript <- 
      Annotation_tr[ ,c("start", "end", "transcript_id")]
    if(any(duplicated(start_stop_coordinates_per_transcript)))
    {
      start_stop_coordinates_per_transcript = 
        start_stop_coordinates_per_transcript[
          !duplicated(start_stop_coordinates_per_transcript), ]
    }
    mmm1 <- merge(mmm, start_stop_coordinates_per_transcript, 
                  by.x = "fsetname", 
                  by.y = "transcript_id")
    
  }
  
  if( level == "gene")
  {
    a4 = dplyr::arrange(
      Annotation_tr, 
      seqid, 
      gene_id, 
      start)[ ,c("gene_id", "start")]
    a5 = dplyr::arrange(Annotation_tr, 
                        seqid, 
                        gene_id, 
                        desc(end))[ ,c("gene_id", "end")]
    
    aa5 = a5[!duplicated(a5$gene_id), ]
    aa4 = a4[!duplicated(a4$gene_id), ]
    start_stop_coordinates_per_gene <- merge(aa4, aa5)
    
    # small tests that start/ stop coord of genes picked right #
    stopifnot(nrow(aa5) == nrow(start_stop_coordinates_per_gene))
    
    example_gene_id <- 
      as.vector(sample(start_stop_coordinates_per_gene$gene_id, 1))
    stopifnot(
      min(
        subset(
          start_stop_coordinates_per_gene, 
          gene_id == example_gene_id)$start) == 
        subset(start_stop_coordinates_per_gene, 
               gene_id == example_gene_id)$start)
    stopifnot(
      max(
        subset(
          start_stop_coordinates_per_gene, 
          gene_id == example_gene_id)$end) == 
        subset(start_stop_coordinates_per_gene, 
               gene_id == example_gene_id)$end)
    message("start stop coordinated of genes correct /n")
    rm(a4, a5, aa4, aa5)
    
    mmm1 <- merge(mmm, 
                  start_stop_coordinates_per_gene, 
                  by.x = "fsetname", 
                  by.y = "gene_id")
  }
  
  if(any(duplicated(mmm1)))
  {
    mmm1 = mmm1[!duplicated(mmm1), ]
  }

  write.csv(mmm1[, c("fsetname", "fsetid")], 
            file = file.path(outputDir, "Chip_probeset_annotation.csv"), 
            row.names = FALSE, 
            quote = FALSE)
  
  colnames(mmm1)[which(colnames(mmm1) == "end")] <- "stop"
  mmm1_1 <- mmm1[ ,c(colnames(ParsedDataOld$featureSet))]
  mmm1_2 <- rbind(nonMain_from_feature_set, mmm1_1)
  if(any(duplicated(mmm1_2$fsetid)))
  {
    mmm1_2 = mmm1_2[!duplicated(mmm1_2$fsetid), ]
  }
  
  # mmm1_2$fsetid <- as.character(mmm1_2$fsetid)
  # mmm1_2$transcript_cluster_id<- as.character(mmm1_2$transcript_cluster_id)
  # 
  # mmm1_2$strand <- as.numeric(mmm1_2$strand)
  # mmm1_2$start <- as.numeric(mmm1_2$start)
  # mmm1_2$stop <- as.character(mmm1_2$stop)
  # mmm1_2$exon_id <- as.numeric(mmm1_2$exon_id)
  # mmm1_2$crosshyb_type <- as.numeric(mmm1_2$crosshyb_type)
  # mmm1_2$type <- as.numeric(mmm1_2$type)
  rownames(mmm1_2) <- 1:nrow(mmm1_2)
  newFeatureSet <- as.data.frame(mmm1_2)

  ####################################
  ############   newCore   ###########
  #table(is.na(ParsedDataOld$core$transcript_cluster_id)) # all FALSE
  nonMain_core <- subset(ParsedDataOld$core, 
                         fsetid %in% nonMain_from_feature_set$fsetid)
  MainCore <- data.frame(meta_fsetid = mmm1$fsetid, 
                         transcript_cluster_id = mmm1$fsetid, 
                         fsetid = mmm1$fsetid )
  newCore <- as.data.frame(rbind(nonMain_core, MainCore))
  if(any(duplicated(newCore$fsetid)))
  {
    newCore = newCore[!duplicated(newCore$fsetid), ]
  }
  
  ##########################################
  ############   newPmFeatures   ###########
  nonMain_pmFeatures <- subset(ParsedDataOld$pmFeatures, 
                               fsetid %in% nonMain_from_feature_set$fsetid)
  ProbeSetData$ProbeSetDataFrame <- 
    as.data.frame(apply(ProbeSetData$ProbeSetDataFrame, 2, as.integer))
  ttt = merge(ProbeSetData$ProbeSetDataFrame, 
              ParsedDataOld$pmFeatures[ ,c("fid", "atom", "x", "y")])
  ttt = dplyr::arrange(ttt, fsetid)
  
  ttt_all <- rbind(nonMain_pmFeatures, ttt)
  ttt_all = dplyr::arrange(ttt_all, atom, fsetid)

  newPmFeatures <- ttt_all
  
  ##########################################
  ############   newPmSequence   ###########
  ppp <- BiocGenerics::subset(ParsedDataOld$pmSequence, 
                              fid %in% newPmFeatures$fid)
  newPmSequence <- ppp
  
  ##########################################
  ############   newMmFeatures   ###########
  if(length(ParsedDataOld$mmFeatures) > 0)
  {
    if (nrow(ParsedDataOld$mmFeatures)>0)
    {
      newMmFeatures <- subset(ParsedDataOld$mmFeatures, 
                              !(fid %in% newPmFeatures$fid))
    }
  }else{newMmFeatures <- ParsedDataOld$mmFeatures}
  ###########################################
  ############   newMmSequences   ###########
  if(length(ParsedDataOld$mmSequences) > 0)
  {
    if (nrow(ParsedDataOld$mmSequences)>0)
    {
      newMmSequence <- BiocGenerics::subset(ParsedDataOld$mmSequence, 
                                            !(fid %in% newPmSequences$fid))
    }
  }else{newMmSequence <- ParsedDataOld$mmSequence}
  ###########################################
  NewParsedData <- list(newFeatureSet, 
                        newPmFeatures, 
                        newMmFeatures, 
                        newGeometry, 
                        newPmSequence, 
                        newMmSequence,
                        newChromDict, 
                        newLeveDict, 
                        newTypeDict, 
                        newCore)
  names(NewParsedData) <- names(ParsedDataOld)
  
  
  
  return(NewParsedData)
}


build_new_ParsedDataOld_object_CDF  <- function(level,
                                               ParsedDataOld, 
                                               PS_controls,
                                               Annotation, 
                                               ProbeSetData)
  
  
{
  PS_controls = subset(PS_controls, !(fid %in% ProbeSetData$ProbeSetDataFrame$fid))
  
  modify_featureSet_CDFlike <- function(old_featureSet, 
                                        newProbeSetDataAnnotation, 
                                        newProbeSetDataFrame, 
                                        ControlPS)
  {
    options(stringsAsFactors = FALSE)
    
    stopifnot(class(ControlPS) == "data.frame")
    stopifnot(setequal(colnames(ControlPS), 
                       c("x", "y", "isPm", 
                         "atom", "man_fsetid", 
                         "fid", "fsetid", "sequence")))
    stopifnot(setequal(colnames(old_featureSet), 
                       c("fsetid", "man_fsetid", "strand")))
    stopifnot(length(unique(newProbeSetDataFrame$fsetid)) == 
                nrow(newProbeSetDataAnnotation))
    
    FS_control <- subset(old_featureSet, 
                         man_fsetid %in% ControlPS$man_fsetid)
    
    n_of_probesets = nrow(newProbeSetDataAnnotation)
    FS_treatnemt <- 
      data.frame(as.numeric(as.vector(newProbeSetDataAnnotation$fsetid)),
                 as.vector(newProbeSetDataAnnotation$fsetname),
                 as.vector(newProbeSetDataAnnotation$strand))
    colnames(FS_treatnemt) <- colnames(old_featureSet)
    
    FS_all <- rbind(FS_control, FS_treatnemt)
    if(any(duplicated(FS_all)))
    {   FS_all = FS_all[!duplicated(FS_all), ]  }
    
    return(FS_all)
  }

  modify_pmFeatures_CDFlike <- function(old_pmFeatures, ControlPS, ProbeSetData)
  {
    #old_pmFeatures <- ParsedData$pmFeatures
    oldpm3 <- old_pmFeatures[ ,c("fid", "x", "y")]
    mm3_control <- merge(ControlPS[ ,c("fid", "fsetid","x", "y")], oldpm3)
    stopifnot(setequal(ControlPS$fid, mm3_control$fid))
    
    mm3_treatment <- merge(ProbeSetData$ProbeSetDataFrame, oldpm3)
    if(!(setequal(unique(mm3_treatment$fid), 
                       unique(ProbeSetData$ProbeSetDataFrame$fid))))
    {
      aa = length(unique(mm3_treatment$fid))
      bb = length(unique(ProbeSetData$ProbeSetDataFrame$fid))
      cc = length(intersect(mm3_treatment$fid, 
                            ProbeSetData$ProbeSetDataFrame$fid))
      warning(message = paste0("setequal(unique(mm3_treatment$fid), unique(ProbeSetData$ProbeSetDataFrame$fid)) is not TRUE: 
                               they are ", aa, " and ", bb, ", with intersection of ", cc ))
      # these does not have to be equal sets,
      # as NEW probes can become PM 
      # (especially for the chips that have same amount of MM and PM).
      # Nevertheless, these sets should not differ much.
    }
    
    mm3_all = rbind(mm3_treatment[ ,c("fid", "fsetid", "x", "y" )], 
                    mm3_control[ ,c("fid", "fsetid", "x", "y" )])
    
    if(any(duplicated(mm3_all)))
    {
      mm3_all = mm3_all[!duplicated(mm3_all), ]
    }
    if(any(duplicated(mm3_all$fid)))
    {warning("Duplicated fid-s in pmFeatures data.frame! \n")}
    
    mm3_all = mm3_all[order(as.numeric(mm3_all$x*max(mm3_all$y) + mm3_all$y)), ]
    return(mm3_all)
  }
  
  modify_mmFeatures_CDFlike <- function(old_mmFeatures, newPmFeatures)
  {
    mm = merge(old_mmFeatures, 
               newPmFeatures, 
               by.x = "fidpm", by.y = "fid", 
               suffixes = c("_OLD", "_NEW"))
    mm5 = mm[ ,c("fid", "fsetid_NEW", "x_OLD", "y_OLD", "fidpm")]
    colnames(mm5) <- colnames(old_mmFeatures)

    mm5 = subset(mm5, !(fid %in% newPmFeatures$fid))
    
    # test #
    rand_idx <- sample(1:nrow(mm5), 1)
    stopifnot(
      subset(
        newPmFeatures, 
        fid == mm5[rand_idx, ]$fidpm)$fsetid == mm5[rand_idx, ]$fsetid)
    
    stopifnot( 
      subset(old_mmFeatures, 
             fidpm == mm5[rand_idx, ]$fidpm)$x == mm5[rand_idx, ]$x)
    stopifnot( 
      subset(old_mmFeatures, 
             fidpm == mm5[rand_idx, ]$fidpm)$y == mm5[rand_idx, ]$y)
    # all fine! #
    return(mm5)
  }

  
  newFeatureSet = 
    modify_featureSet_CDFlike(
      old_featureSet = ParsedDataOld$featureSet, 
      newProbeSetDataAnnotation = ProbeSetData$ProbeSetDataAnnotation, 
      newProbeSetDataFrame = ProbeSetData$ProbeSetDataFrame, 
      ControlPS = PS_controls)
  
  newPmFeatures <- 
    modify_pmFeatures_CDFlike(
      old_pmFeatures = ParsedDataOld$pmFeatures, 
      ControlPS = PS_controls, 
      ProbeSetData = ProbeSetData)
  
  newPmSequence <- 
    BiocGenerics::subset(ParsedDataOld$pmSequence, fid %in% newPmFeatures$fid)
                            
  newMmFeatures <- 
    modify_mmFeatures_CDFlike(
      old_mmFeatures = ParsedDataOld$mmFeatures, 
      newPmFeatures = newPmFeatures)
  
  newMmSequence <- 
    BiocGenerics::subset(
      ParsedDataOld$mmSequence, !(fid %in% newPmFeatures$fid))
  
  NewParsedData <- list(newFeatureSet, 
                        newPmSequence, 
                        newPmFeatures, 
                        newMmSequence, 
                        newMmFeatures, 
                        ParsedDataOld$geometry)
  
  names(NewParsedData) <- names(ParsedDataOld)
  
  # debugging msg start
  cat(paste0("From inside build_new_ParsedDataOld_object_CDF:  ", names(NewParsedData)))
  # debugging msg end
  
  return(NewParsedData)
}