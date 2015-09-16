library(data.table)

# RAP-MSU ----------------------------------------

# download RAP-MSU from rap-DB and save it as Rdata
temp <- tempfile()
download.file("http://rapdb.dna.affrc.go.jp/download/archive/RAP-MSU.txt.gz", 
    temp)
RAPMSU <- data.table::as.data.table(read.delim(temp, sep = "\t", header = FALSE,
                                               stringsAsFactors = FALSE,
                                               col.names = c("Rap_ID", "MSU_Transcripts_dirty")))
# tidy MSU column
splitDirtyTransctripts <- function(ids) {
  splitIds <- unlist(strsplit(ids, ",", fixed = TRUE))
  unique(gsub("\\.[[:digit:]]+$", "", splitIds))
}
RAPMSU <- RAPMSU[, .(MSU_ID = splitDirtyTransctripts(MSU_Transcripts_dirty)),
                 by = Rap_ID]
setkey(RAPMSU, MSU_ID)

# timestamp
attr(RAPMSU, "dateRetrieved") <- date()

# GeneListById ----------------------------------------

# get a cookie from shigen and extract the jsessionid
temp <- tempfile()
download.file("http://www.shigen.nig.ac.jp/rice/oryzabase/gene/download", 
    temp, method = "wget")
rl <- readLines(temp)
jsessionid <- gsub(".*(jsessionid=[0-9A-Z]+).*", "\\1", rl[grep("jsessionid", 
    rl)][1])

# use the jsessionid to make a new url and download the file
fileUrl <- paste0("http://www.shigen.nig.ac.jp/rice/oryzabase/gene/download;", 
    jsessionid, "?classtag=GENE_EN_LIST")
temp <- tempfile()
download.file(fileUrl, temp, method = "wget")
GeneList <- read.delim(temp, stringsAsFactors = FALSE, na.strings = "_")

# mung the geneList
GeneList[is.na(GeneList)] <- ""
GeneList.table <- data.table::data.table(GeneList, key = "RAP.ID")
GeneListByID <- GeneList.table[!"", .(
  RAP.ID = gsub(" ", "", RAP.ID),
  CGSNL.Gene.Symbol,
  Gene.symbol.synonym.s.,
  CGSNL.Gene.Name,
  Gene.name.synonym.s.,
  Gramene.ID)]
data.table::setnames(GeneListByID, old = "RAP.ID", new = "RAP_id")
data.table::setkey(GeneListByID, "RAP_id")

# tidy data: one line per symbol
# remove symbols that == gene name here
splitSymbol <- function(symbols, CGSNL.Gene.Symbol) {
  splitSymbols <- unlist(strsplit(symbols, ",|;|/"))
  # remove whitespace
  splitSymbols <- sapply(splitSymbols, gsub, pattern = "^[[:space:]]+|[[:space:]]+$", replacement = "")
  # remove "Os"
  splitSymbols <- sapply(splitSymbols, gsub, pattern = "^Os[[:space:]]*", replacement = "")
  # remove symbols that are similar to CGSNL.Gene.Symbol
  splitSymbols <- splitSymbols[
    !toupper(gsub("[^[:alnum:]]", "", splitSymbols)) %in%
      toupper(gsub("[^[:alnum:]]", "", CGSNL.Gene.Symbol))]
  # small variations in formatting
  if(length(splitSymbols) > 1) {
    alnum <- sapply(splitSymbols, gsub, pattern = "[^[:alnum:]]", replacement = "")
  keep <- sapply(2:length(alnum), function(i)
    !(toupper(alnum[i]) %in% toupper(alnum[1:(i-1)])))
  return(splitSymbols[c(TRUE, keep)])
  } else {return(splitSymbols)}
}
splitName <- function(geneNames) {
  splitNames <- unlist(strsplit(geneNames, ",|;|/"))
  # remove whitespace
  splitNames <- sapply(splitNames, gsub, pattern = "^[[:space:]]+|[[:space:]]+$", replacement = "")
  # small variations in formatting
  if(length(splitNames) > 1) {
  alnum <- sapply(splitNames, gsub, pattern = "[^[:alnum:]]", replacement = "")
  keep <- sapply(2:length(alnum), function(i)
    !(toupper(alnum[i]) %in% toupper(alnum[1:(i-1)])))
  return(splitNames[c(TRUE, keep)])
  } else {return(unique(splitNames))}
}
  
symbol_synonyms.table <- GeneListByID[, .(
  symbol_synonyms = splitSymbol(Gene.symbol.synonym.s., CGSNL.Gene.Symbol)), by = RAP_id]
name_synonyms.table <- GeneListByID[, .(
  name_synonyms = splitName(Gene.name.synonym.s.)), by = RAP_id]
synonyms <- merge(symbol_synonyms.table, name_synonyms.table, allow.cartesian = TRUE)

# combine symbols with other information
GeneListWithSynonyms <- synonyms[GeneListByID[, .(RAP_id, CGSNL.Gene.Symbol,
                                                  CGSNL.Gene.Name, Gramene.ID)],
                                 allow.cartesian = TRUE]
data.table::setcolorder(GeneListWithSynonyms,
                        c('RAP_id', 'CGSNL.Gene.Symbol', 'symbol_synonyms',
                          'CGSNL.Gene.Name', 'name_synonyms', 'Gramene.ID'))
data.table::setkey(GeneListWithSynonyms, NULL)
GeneListWithSynonyms <- unique(GeneListWithSynonyms)
data.table::setkey(GeneListWithSynonyms, 'RAP_id')

# timestamp
attr(GeneListWithSynonyms, "dateRetrieved") <- date()

# RapMsuRefSeq ----------------------------------------

# download data.table of rice id build table from oryzabase
temp <- tempfile()
download.file("http://www.shigen.nig.ac.jp/rice/oryzabase/tool/riceIdChecker/build5/download", 
    method = "auto", temp)
RapMsuRefSeq <- data.table::data.table(read.table(unzip(temp, unzip(temp, 
    list = TRUE)[, 1], exdir = "data-raw", overwrite = TRUE), header = TRUE, 
    fill = TRUE, sep = "\t", stringsAsFactors = FALSE))
data.table::setkey(RapMsuRefSeq, tigrId)

# timestamp
attr(RapMsuRefSeq, "dateRetrieved") <- date()

# Save data ----------------------------------------

devtools::use_data(RAPMSU, GeneListWithSynonyms, RapMsuRefSeq, internal = TRUE, 
    overwrite = TRUE) 
