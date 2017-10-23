#' Filter out probe sets that have too few probes.
#'
#' Filter out probe sets that have too few probes.
#'
#' @title filter out small probe sets
#' @param ProbeSets \code{data.frame}
#' @param level \code{character}, must equal "gene" or "transcript"
#' @param min_probe_number \code{numeric}, probe sets with >= probes will remain
#' @return a \code{"output type"}.
#' @author Vladislava Milchevskaya \email{milchv@gmail.com}
# clean probe sets according to probe sets filtering parameters #
cleanProbeSets <- function(ProbeSets, level, min_probe_number = 3)
{
  stopifnot(level %in% c("gene", "transcript"))
  stopifnot(class(ProbeSets) == "data.frame")
  stopifnot(setequal(colnames(ProbeSets), 
                     c("transcript_id", "gene_id", "probe_id")))

  if(level == "transcript")
  {
    mm1 = ProbeSets
    mm2 = mm1[!duplicated(mm1[ ,c("probe_id","transcript_id")]), ,drop = FALSE]
    mm2s = split(mm2, mm2$transcript_id )
    mm2sl = lapply(mm2s, FUN = function(zz) length(unique(zz$probe_id)))
    mm2sl_melted = reshape2::melt(mm2sl)
    m3 = subset(mm2sl_melted, value >= min_probe_number)
    ProbeSets_filtered = subset(ProbeSets, transcript_id %in% unique(m3$L1))
  }
  if(level == "gene")
  {
    mm1 = ProbeSets
    mm2 = mm1[!duplicated(mm1[ ,c("probe_id","gene_id")]), ,drop = FALSE]
    mm2s = split(mm2, mm2$gene_id )
    mm2sl = lapply(mm2s, FUN = function(zz) length(unique(zz$probe_id)))
    mm2sl_melted = reshape2::melt(mm2sl)
    m3 = subset(mm2sl_melted, value >= min_probe_number)
    ProbeSets_filtered = subset(ProbeSets, gene_id %in% unique(m3$L1))
    ProbeSets_filtered$transcript_id <- rep(NA, nrow(ProbeSets_filtered))
  }
  ProbeSets_filtered = ProbeSets_filtered[!duplicated(ProbeSets_filtered), ]
  ProbeSets_filtered = ProbeSets_filtered[order(ProbeSets_filtered$probe_id), ]
  if(any(duplicated(ProbeSets_filtered$probe_id)))
  {
    message("Probes have multipple mappings! 
            Check level value: must agree in all function. 
            Preferably use level = 'gene' everywhere.")
  }
  return(ProbeSets_filtered)
}

# tmp = cleanProbeSets(PS_dro, level = "transcript")
# head(tmp)
# dim(tmp)
#
