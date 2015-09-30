#' dateF
#' 
#' date formatted for console messages
dateF <- function() {
  return(format(Sys.time(), "[ %a %b %d %H:%M:%S %Y ]"))
} 

#' showDateRetrieved
#' 
#' Show the dates when the databases (RAP-DB to MSU IDs, OryzabaseGeneListEn and
#' Rice ID IRGSP Build5 database) were last downloaded to generate
#' \code{R/sysdata.rda}.
#' @return Returns a named character vector of date retrieved for each database.
showDateRetrieved <- function() {
  return(c(
    RapMsuRefSeq = attr(RapMsuRefSeq, 'dateRetrieved'),
    RAPMSU = attr(RAPMSU, 'dateRetrieved'),
    GeneListWithSynonyms = attr(GeneListWithSynonyms, 'dateRetrieved')
  ))
}