library(data.table)

###########
# RAP-MSU #
###########

# download RAP-MSU from rap-DB and save it as Rdata
temp <- tempfile()
rap.msu.url <-
  "http://rapdb.dna.affrc.go.jp/download/archive/RAP-MSU_2016-08-05.txt.gz"
download.file(rap.msu.url, temp)
RAPMSU <- data.table::as.data.table(
  read.delim(temp, sep = "\t", header = FALSE, stringsAsFactors = FALSE,
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
attr(RAPMSU, "dateRetrieved") <- Sys.time()

################
# GeneListById #
################

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
  splitSymbols <- unlist(strsplit(symbols, ",|;|/|[[:space:]]|\\(|\\)"))
  #splitSymbols <- unlist(strsplit(symbols, "[^[:alnum:]]"))
  # remove whitespace
  splitSymbols <- sapply(splitSymbols, gsub,
                         pattern = "^[[:space:]]+|[[:space:]]+$",
                         replacement = "")
  # remove "Os"
  splitSymbols <- sapply(splitSymbols, gsub, pattern = "^Os[[:space:]]*",
                         replacement = "")
  # remove blanks
  splitSymbols <- splitSymbols[sapply(splitSymbols, nchar) != 0]
  # remove symbols that are similar to CGSNL.Gene.Symbol
  splitSymbols <- splitSymbols[
    !toupper(gsub("[^[:alnum:]]", "", splitSymbols)) %in%
      toupper(gsub("[^[:alnum:]]", "", CGSNL.Gene.Symbol))]
  # small variations in formatting
  if(length(splitSymbols) > 1) {
    alnum <- sapply(splitSymbols, gsub, pattern = "[^[:alnum:]]",
                    replacement = "")
  keep <- sapply(2:length(alnum), function(i)
    !(toupper(alnum[i]) %in% toupper(alnum[1:(i-1)])))
  return(splitSymbols[c(TRUE, keep)])
  } else {return(splitSymbols)}
}
splitName <- function(geneNames) {
  splitNames <- unlist(strsplit(geneNames, ",|;|/"))
  # remove whitespace
  splitNames <- sapply(splitNames, gsub,
                       pattern = "^[[:space:]]+|[[:space:]]+$",
                       replacement = "")
  # small variations in formatting
  if(length(splitNames) > 1) {
  alnum <- sapply(splitNames, gsub, pattern = "[^[:alnum:]]", replacement = "")
  keep <- sapply(2:length(alnum), function(i)
    !(toupper(alnum[i]) %in% toupper(alnum[1:(i-1)])))
  return(splitNames[c(TRUE, keep)])
  } else {return(unique(splitNames))}
}
  
symbol_synonyms.table <- GeneListByID[, .(
  symbol_synonyms = splitSymbol(Gene.symbol.synonym.s., CGSNL.Gene.Symbol)),
  by = RAP_id]
name_synonyms.table <- GeneListByID[, .(
  name_synonyms = splitName(Gene.name.synonym.s.)), by = RAP_id]
synonyms <- merge(symbol_synonyms.table, name_synonyms.table, 
                  allow.cartesian = TRUE)

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
attr(GeneListWithSynonyms, "dateRetrieved") <- Sys.time()

####################
# TIGR annotations #
####################

temp <- tempfile()
download.file("ftp://ftp.plantbiology.msu.edu/pub/data/Eukaryotic_Projects/o_sativa/annotation_dbs/pseudomolecules/version_7.0/all.dir/all.locus_brief_info.7.0",
              method = "auto", temp)
msu.annotation <- data.table(read.delim(file = temp, sep = "\t",header = TRUE,
                      fill = TRUE), key = 'locus')  
# get 1 annotation per LOC ID
msu.annotation.collapsed <- unique(msu.annotation[, .(
  annotation = paste(unique(annotation), collapse = ", ")),
  by = locus])

# timestamp
attr(msu.annotation.collapsed, "dateRetrieved") <- Sys.time()

################
# RapMsuRefSeq #
################

# download data.table of rice id build table from oryzabase
temp <- tempfile()
download.file(
  "http://shigen.nig.ac.jp/rice/oryzabase/tool/riceIdChecker/download", 
  method = "auto", temp)
RapMsuRefSeq <- data.table::data.table(read.table(unzip(temp, unzip(temp, 
    list = TRUE)[, 1], exdir = "data-raw", overwrite = TRUE), header = TRUE, 
    fill = TRUE, sep = "\t", stringsAsFactors = FALSE))
data.table::setkey(RapMsuRefSeq, tigrId)

# timestamp
attr(RapMsuRefSeq, "dateRetrieved") <- Sys.time()

########
# OGRO #
########

temp <- tempfile()
download.file("http://qtaro.abr.affrc.go.jp/ogro/table/export?format=csv", 
              temp)
ogro <- data.table(
  (read.delim(temp, sep = ",", header = TRUE, stringsAsFactors = FALSE)),
  key = "locus_id")

# remove html quotation marks from objective field
ogro[, objective := gsub("&quot;", "", objective, fixed = TRUE)]
# ignore na ids, collapse doi and character
ogro.processed <- ogro[!is.na(locus_id), .(
  ogro.objective = paste(unique(objective), collapse = " "),
  ogro.doi = paste(unique(doi), collapse = ", ")),
  by = locus_id]

setkey(RAPMSU, Rap_ID)
rap.msu.ogro <- ogro.processed[RAPMSU]
setkey(rap.msu.ogro, MSU_ID)

# timestamp
attr(rap.msu.ogro, "dateRetrieved") <- Sys.time()

###############################
# MSUv7 GFF for GenomicRanges #
###############################

# urls
signon.url <- "https://signon.jgi.doe.gov/signon/create"
annot.url <- 
  "http://genome.jgi.doe.gov/Osativa/download/_JAMO/53112ab649607a1be00559b0/Osativa_204_v7.0.gene_exons.gff3.gz"

# check dependencies
if (!requireNamespace("httr", quietly = TRUE)) {
  stop("httr needed to download GTF from JGI. Please install it.",
       call. = FALSE)
}
if (!requireNamespace("GenomicRanges", quietly = TRUE)) {
  stop("GenomicRanges needed to download GTF from JGI. Please install it.",
       call. = FALSE)
}
if (!requireNamespace("rtracklayer", quietly = TRUE)) {
  stop("rtracklayer needed to download GTF from JGI. Please install it.",
       call. = FALSE)
}

# get login/password
message(
  paste("You need to provide your JGI login and password to download the",
        "MSUv7 GFF from Phytozome",
        sep = "\n")
  )
while(is.na(jgi.email)) {
  jgi.email <- readline("Email address: ")
  jgi.email <- ifelse(grepl("[^[:blank:]+]@[^[:blank:]+]", jgi.email),
                      jgi.email,
                      NA)}
while(is.na(jgi.password)) {
  jgi.password <- readline("Password: ")
  jgi.password <- ifelse(grepl("[[:blank:]]", jgi.password),
                         NA,
                         jgi.password)}

# set a handle
jgi.handle <- httr::handle(url = "http://jgi.doe.gov")

# login
cookie.response <- httr::POST(
  url = signon.url,
  body = list(login = jgi.email,
              password = jgi.password),
  encode = "form",
  handle = jgi.handle)

# get the file
file.response <- httr::GET(
  url = annot.url,
  handle = jgi.handle)

# write to a temporary file
tmp <- tempfile(fileext = ".gff3.gz")
writeBin(httr::content(file.response), tmp)

# gffread
gff <- rtracklayer::import.gff(tmp, format = 'gff3',
                               genome = 'Osativa_204_v7',
                               feature.type="exon")
gff.exons <- GenomicRanges::split(
  x = gff,
  f = GenomicRanges::elementMetadata(gff)$Parent)[[1]]

# timestamp
attr(gff.exons, "dateRetrieved") <- Sys.time()

#############
# Save data #
#############

devtools::use_data(RAPMSU, GeneListWithSynonyms, RapMsuRefSeq,
                   msu.annotation.collapsed, rap.msu.ogro, gff.exons,
                   internal = TRUE, 
                   overwrite = TRUE) 
