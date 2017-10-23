#' Converts CIGAR string of the mappint to length of the mapping
#'
#' Converts CIGAR string to length of the mapping fragment; copied from Java script to some alignment processing tool..
#'
#' @title CIGAR-to-length
#' @param cigar_string character, CIGAR string from the alignment \code{data.frame}
#' @return returns \code{integer}.
#' @examples
#' cigar_character = "20M1S2N5M"
#' mapping_length = cigar2length(cigar_character)
#' mapping_length
#' @author Vladislava Milchevskaya \email{milchv@gmail.com}
#' @importFrom stringr str_extract_all
#' @export
cigar2length <- function(cigar_string)
{
  cigar_allowed = c(as.character(c(0:9)), 
                    "M", "I", "D", "N", "S", "H", "P", "X", "=")
  cigar_by_symbol = unlist(strsplit(cigar_string, split = ""))
  if (!all(cigar_by_symbol %in% cigar_allowed))
    stop(
      message = "Wrong symbols in CIGAR string (check alognment SAM files). \n")

  M_values = as.numeric(unlist(stringr::str_extract_all(
    stringr::str_extract_all(cigar_string, "[0-9]+[M]"), "[0-9]+")))
  D_values = as.numeric(unlist(stringr::str_extract_all(
    stringr::str_extract_all(cigar_string, "[0-9]+[D]"), "[0-9]+")))
  N_values = as.numeric(unlist(stringr::str_extract_all(
    stringr::str_extract_all(cigar_string, "[0-9]+[N]"), "[0-9]+")))
  EQ_values = as.numeric(unlist(stringr::str_extract_all(
    stringr::str_extract_all(cigar_string, "[0-9]+[=]"), "[0-9]+")))
  X_values = as.numeric(unlist(stringr::str_extract_all(
    stringr::str_extract_all(cigar_string, "[0-9]+[X]"), "[0-9]+")))
  P_values = as.numeric(unlist(stringr::str_extract_all(
    stringr::str_extract_all(cigar_string, "[0-9]+[P]"), "[0-9]+")))

  res_length = 
    sum(M_values) + 
    sum(D_values) + 
    sum(N_values) + 
    sum(EQ_values) + 
    sum(X_values) + 
    sum(P_values)
  
  return(res_length)
}
