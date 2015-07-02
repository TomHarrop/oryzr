#' SearchByGeneSymbol
#' 
#' Search Rap-DB and oryzabase information by gene symbol and return identifiers
#' for potential matches.
#' 
#' Symbols are matched against the OryzabaseGeneListEn database, downloaded from
#' \href{http://www.shigen.nig.ac.jp/rice/oryzabase/download/gene}{Oryzabase}. 
#' Rap-DB identifiers from the OryzabaseGeneListEn database are converted to MSU
#' identifiers using the table of RAP-DB and MSU IDs 
#' \href{http://rapdb.dna.affrc.go.jp/download/irgsp1.html}{available from 
#' Rap-DB}. \emph{n.b.} gene symbols are not necessarily unique so the 
#' \code{data.table} returned by \code{SearchByGeneSymbol} must be manually 
#' verified.
#' 
#' @param GeneSymbols Character vector of gene symbols
#'   
#' @return Returns a \code{data.frame} with MSU IDs as \code{rownames}, and
#'   columns RapID (RAP-DB gene identifier), symbols (CSGNL recommended symbols)
#'   and names (CSGNL recommended names).
#'   
#' @import data.table
#' @export
#' 
#' @examples
#' SearchByGeneSymbol(c('FON1'))
#' SearchByGeneSymbol(c('SPL7'))

SearchByGeneSymbol <- function(GeneSymbols) {
  # Get RapDB ID from OryzabaseGeneListEn
  symbolMatches <- as.data.table(do.call(
    rbind,
    lapply(GeneSymbols, function(x)
      GeneListByID.frame[grep(x, GeneListByID.frame$symbol,
                              ignore.case = TRUE),]
    )))
  # Match RapDB ID to MSU ID
  RapMsuMatches <- as.data.table(do.call(
    rbind,
    lapply(unique(symbolMatches$RAP_id), function(x)
      RAPMSU[grep(x, RAPMSU$Rap_ID, ignore.case = TRUE),]
    )))
  # Trim transcript numbers
  RapMsuMatches[,MSU_ID := gsub("\\..*", "", MSU_ID)]
  
  # Return gene information
  return(
  oryzr::LocToGeneName(RapMsuMatches$MSU_ID, plotLabels = FALSE)
  )
  
}

