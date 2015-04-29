#' dateF
#' 
#' date formatted for console messages
dateF <- function() {
  return(format(Sys.time(), "[%a %b %d %H:%M:%S %Y]"))
} 

#' showDateRetrieved
#' 
#' Show the dates when the databases (RAP-DB to MSU IDs, OryzabaseGeneListEn and
#' Rice ID IRGSP Build5 database) were last downloaded to generate
#' \code{R/sysdata.rda}.
#' @return Returns a named character vector of date retrieved for each database.
showDateRetrieved <- function() {
  getDatestamps <- function(sysData){
    env <- new.env()
    load(sysData, env)
    c(
      RapMsuRefSeq = attr(env[['RapMsuRefSeq']], 'dateRetrieved'),
      RAPMSU = attr(env[['RAPMSU']], 'dateRetrieved'),
      GeneListByID.frame = attr(env[['GeneListByID.frame']], 'dateRetrieved')
    )
  }
  getDatestamps('R/sysdata.rda')
}