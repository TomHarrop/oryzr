#' LocToGeneName
#' 
#' Convert MSU identifiers to CSGNL gene names and symbols and generate a label
#' for plotting if requested.
#' 
#' Conversion from MSU to Rap-DB identifiers uses the table of RAP-DB and MSU 
#' IDs \href{http://rapdb.dna.affrc.go.jp/download/irgsp1.html}{available from 
#' Rap-DB}. CSGNL symbols and gene names are retrieved from the 
#' OryzabaseGeneListEn database, downloaded from 
#' \href{http://www.shigen.nig.ac.jp/rice/oryzabase/download/gene}{Oryzabase}.
#' 
#' @param LOCs Character vector of MSU locus identifiers
#' @param shortLabels Boolean (default \code{FALSE}). Tries to return a 
#'   sensible label \emph{e.g.} for plotting: either the LOC_ID for unnamed 
#'   genes or the name for named genes, with synonyms in brackets.
#' @param return.synonyms Boolean (default \code{FALSE}). Return a 
#'   \strong{long} data.table including a column with the synonyms for each 
#'   gene.
#' @return Returns a \code{data.table} with columns MsuID (MSUv7 gene 
#'   identifier), RapID (RAP-DB gene identifier), symbols (CSGNL recommended 
#'   symbol), names (CSGNL recommended name), MsuAnnotation 
#'   (\href{http://rice.plantbiology.msu.edu/}{TIGR} annotations), 
#'   OgroObjective and OgroRef (manual annotation from 
#'   \href{http://qtaro.abr.affrc.go.jp/ogro}{OGRO} database) and optionally 
#'   labels for plotting. If \code{return.synonyms} is \code{TRUE}, returns a 
#'   \strong{long \code{data.table}} additionally containing columns
#'   symbol_synonyms and name_synonyms.
#'   
#' @export
#' 
#' @import data.table
#'   
#' @examples 
#' LocToGeneName(c('LOC_Os01g03340', 'LOC_Os01g09410', 'LOC_Os01g17396', 'LOC_Os01g18800', 'LOC_Os01g44069', 'LOC_Os01g45470', 'LOC_Os01g53880', 'LOC_Os01g64730', 'LOC_Os01g71820', 'LOC_Os02g07110', 'LOC_Os02g26680', 'LOC_Os02g26700', 'LOC_Os02g41904', 'LOC_Os02g42880', 'LOC_Os02g44630', 'LOC_Os03g15530', 'LOC_Os03g22270', 'LOC_Os03g61160', 'LOC_Os04g58710', 'LOC_Os06g06750', 'LOC_Os07g32170', 'LOC_Os08g24930', 'LOC_Os08g31980', 'LOC_Os09g07154', 'LOC_Os09g25330', 'LOC_Os09g31438', 'LOC_Os10g07229', 'LOC_Os10g18870', 'LOC_Os10g26340', 'LOC_Os10g42210', 'LOC_Os10g42220'))
#' LocToGeneName(c('LOC_Os03g03034', 'LOC_Os06g36560', 'LOC_Os01g50910', "LOC_Os03g60430"))
#' LocToGeneName(c("LOC_Os01g03840", "LOC_Os01g10110", "LOC_Os01g10504", "LOC_Os01g40630", "LOC_Os01g61480", "LOC_Os01g66030", "LOC_Os01g66290", "LOC_Os02g15950", "LOC_Os02g45770", "LOC_Os02g52340", "LOC_Os03g11600", "LOC_Os03g11614", "LOC_Os03g51690", "LOC_Os03g54170", "LOC_Os03g57240", "LOC_Os03g60430", "LOC_Os04g32510", "LOC_Os04g49150", "LOC_Os04g51000", "LOC_Os05g41760", "LOC_Os06g11330", "LOC_Os06g45460", "LOC_Os06g49840", "LOC_Os06g50340", "LOC_Os07g13170", "LOC_Os07g15770", "LOC_Os07g41370", "LOC_Os07g42410", "LOC_Os07g47330", "LOC_Os08g06480", "LOC_Os08g07740", "LOC_Os08g39890", "LOC_Os08g41950", "LOC_Os09g26999", "LOC_Os09g32948", "LOC_Os10g33780", "LOC_Os11g38270", "LOC_Os12g10540"))

LocToGeneName <- function(LOCs, shortLabels = FALSE, return.synonyms = FALSE) {
  LocToRap <- data.table(MsuID = unique(LOCs), key = "MsuID")
  # get RAP IDs
  LocToRap <- rap.msu.ogro[LocToRap, .(
    MsuID = MSU_ID,
    RapID = locus_id,
    OgroObjective = ogro.objective,
    OgroRef = ogro.doi
  )]
  LocToRap[RapID == "None", RapID := NA]
  setkey(LocToRap, "RapID")
  
  # get labels
  LocToLabels <- GeneListWithSynonyms[LocToRap, .(
    MsuID,
    RapID,
    symbols = CGSNL.Gene.Symbol,
    names = CGSNL.Gene.Name,
    symbol_synonyms,
    name_synonyms,
    OgroObjective,
    OgroRef)]
  # replace blanks with NA
  rmNaDt = function(DT) {
    for (j in seq_len(ncol(DT)))
      set(DT,which(DT[[j]] == ""),j,NA)
  }
  rmNaDt(LocToLabels)
  
  # 1. loc, no symbol, no synonyms
  LocToLabels[is.na(symbols) & is.na(symbol_synonyms), labels := MsuID,
              by = MsuID]
  
  # 2. symbol, no synonyms
  LocToLabels[(!is.na(symbols)) & is.na(symbol_synonyms), labels := symbols,
              by = MsuID]
  
  # 3. loc, no symbol, synonyms
  LocToLabels[is.na(symbols) & !is.na(symbol_synonyms),
              labels := paste0(MsuID, " (", paste(unique(symbol_synonyms), collapse = "/"), ")"),
              by = MsuID]
  
  # 4. symbol and synonyms
  LocToLabels[(!is.na(symbols)) & !is.na(symbol_synonyms),
              labels := paste0(symbols, " (", paste(unique(symbol_synonyms), collapse = "/"), ")"),
              by = MsuID]
  
  setkey(LocToLabels, MsuID)
  
  # get MSU annotation
  LocToLabels <- msu.annotation.collapsed[LocToLabels, .(
    MsuID = locus,
    MsuAnnotation = annotation,
    RapID,symbols,names,symbol_synonyms,name_synonyms,OgroObjective,OgroRef,
    labels
  )]
  setkey(LocToLabels, MsuID)
  
  # don't return synonyms or labels. In this case return a data.frame.
  if (!shortLabels & !return.synonyms) {
    warning("RAP-DB data changed. LocToGeneName now returns a data.table")
    return(unique(LocToLabels[,.(MsuID,
                                 RapID,
                                 symbols,
                                 names,
                                 MsuAnnotation,
                                 OgroObjective,
                                 OgroRef)],
           by = c("MsuID", "symbols", "names", "MsuAnnotation")))
  }
  
  # return short labels, don't return synonyms
  if (shortLabels & !return.synonyms){
    return(unique(LocToLabels[,.(MsuID,
                                 RapID,
                                 symbols,
                                 names,
                                 MsuAnnotation,
                                 OgroObjective,
                                 OgroRef,
                                 labels)],
           by = c("MsuID", "symbols", "names", "MsuAnnotation")))
  }
  
  # don't return short labels, return synonyms
  if (!shortLabels & return.synonyms) {
    warning("Synonyms requested, returning **LONG** data.table")
    return(LocToLabels[, .(MsuID,
                           RapID,
                           symbols,
                           symbol_synonyms,
                           names,
                           name_synonyms,
                           MsuAnnotation,
                           OgroObjective,
                           OgroRef)])
  }
  
  # short labels and synonyms
  if (shortLabels & return.synonyms) {
    warning("Synonyms requested, returning **LONG** data.table")
  }
  # return the whole table
  return(LocToLabels[, .(MsuID,
                         RapID,
                         symbols,
                         symbol_synonyms,
                         names,
                         name_synonyms,
                         MsuAnnotation,
                         OgroObjective,
                         OgroRef,
                         labels)])
} 
