#' RapToMsu
#' 
#' Convert Rap-DB identifiers to MSU identifiers.
#' 
#' Conversion uses the table of RAP-DB and MSU IDs
#' \href{http://rapdb.dna.affrc.go.jp/download/irgsp1.html}{available from 
#' Rap-DB}.
#' 
#' @param rap.id Character vector of RAP-DB locus identifiers
#' @return Returns a \code{data.frame} with two colums: rap.id and msu.id.
#'   
#' @export
#' 
#' @import data.table
#'   
#' @examples 
#' RapToMsu(c("Os02g0278700", "Os04g0611800", "Os06g0570600", "Os06g0569900", "Os06g0110000"))
#' RapToMsu(c("Os10g0138100", "Os03g0790900", "Os07g0281700", "Os02g0635200", "Os02g0635000", "Os01g0645400", "Os05g0528600", "Os01g0732600", "Os01g0224700", "Os12g0512000", "Os07g0437000", "Os04g0128900", "Os03g0162000", "Os05g0150500", "Os02g0759700", "Os04g0395600", "Os11g0462900", "Os11g0515500", "Os01g0178500") 

RapToMsu <- function(rap.id) {
  rap.id.table <- data.table(rap.id, key = "rap.id")
  setkey(RAPMSU, Rap_ID)
  matched.genes <- RAPMSU[rap.id.table, .(rap.id = Rap_ID, msu.id = MSU_ID)]
  matched.genes[msu.id == "None", msu.id := NA]
  data.frame(matched.genes)
  }