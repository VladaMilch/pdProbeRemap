#' filter out non-unique mappings in Transcriptomic alignment
#'
#' Annotated transcriptomic alignment \code{data.frame} processed,
#' probe IDs with more than one gene or transcript (see \code{"level"}) in the data.frame called
#' non-unique mappings and filtered out, all rows with these probes excluded.
#'
#' @title filter out non-unique transcriptomic mappings
#' @param AT annoated transcriptomic alignment, \code{"data.frame"}
#' @param level character string, takes values \code{"gene"} or \code{"transcript"}
#' @return a \code{"data.frame"}
#' @author Vladislava Milchevskaya \email{milchv@gmail.com}

## filteres out probes that have non-unique mappings ##

filter.T.alignment <- function(AT, level)
{
  stopifnot( c("probe_id", "gene_id", "transcript_id") %in% colnames(AT) )

  if(level == "gene")
  {aaa1 = AT[,c("probe_id", "gene_id")]}
  if(level == "transcript")
  {aaa1 = AT[,c("probe_id", "transcript_id")]}

  aaa2 = aaa1[!duplicated(aaa1),]
  ad2 = aaa2[duplicated(aaa2[,c("probe_id")]),]
  non_unique = (unique(ad2$probe_id))
  message(paste0("Non-unique mappings detected: ", length(non_unique), " \n"))
  aaa3 = subset(aaa2, !(probe_id %in% non_unique))

  AT.filtered = subset(AT, !(probe_id %in% non_unique))
  stopifnot(class(AT.filtered) == "data.frame")

  if(nrow(AT.filtered) == 0)
  {stop (message = "No probes left after filtering. \n")}

  return(AT.filtered)
}
