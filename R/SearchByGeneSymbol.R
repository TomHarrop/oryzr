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
#' @param \code{GeneSymbols} Character vector of gene symbols
#' @param shortLabels Boolean (default \code{FALSE}). Tries to return a sensible
#'   label \emph{e.g.} for plotting. See help for \code{LocToRefSeq} for more
#'   information.
#' @param return.synonyms Boolean (default \code{FALSE}). Return a \strong{long}
#'   data.table including a column with the synonyms for each gene. See help for
#'   \code{LocToRefSeq} for more information.
#' 
#' @return Returns a \code{data.frame} with MSU IDs as \code{rownames}, and 
#'   columns RapID (RAP-DB gene identifier), symbols (CSGNL recommended symbols)
#'   and names (CSGNL recommended names), unless either of \code{shortLabels} or
#'   \code{return.synonyms} is \code{TRUE}. See help for \code{LocToRefSeq} for
#'   more information.
#'   
#' @import data.table
#' @export
#' 
#' @examples
#' SearchByGeneSymbol(c('FON1'))
#' SearchByGeneSymbol(c('SPL7'))

SearchByGeneSymbol <- function(GeneSymbols, shortLabels = FALSE, return.synonyms = FALSE) {
  # Get RapDB ID from OryzabaseGeneListEn
  getMatch <- function(GeneSymbol) {
    GeneListWithSynonyms[grepl(toupper(GeneSymbol),
                               toupper(CGSNL.Gene.Symbol)) |
                           grepl(toupper(GeneSymbol),
                                 toupper(symbol_synonyms))]
  }
  symbolMatches <- do.call(rbind, lapply(GeneSymbols, getMatch))
  # Match RapDB ID to MSU ID
  rapMsuMatches <- symbolMatches[, RAPMSU[Rap_ID == unique(RAP_id), MSU_ID]]
  # Return gene information
  return(
    oryzr::LocToGeneName(rapMsuMatches, shortLabels = shortLabels,
                         return.synonyms = return.synonyms)
  )
}

