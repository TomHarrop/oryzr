###############
### RAP-MSU ###
###############

# download RAP-MSU from rap-DB and save it as Rdata
temp <- tempfile()
download.file('http://rapdb.dna.affrc.go.jp/download/archive/RAP-MSU.txt.gz', temp)
RAPMSU <- read.delim(temp, sep = '\t', header = FALSE, stringsAsFactors = FALSE, col.names = c('Rap_ID', 'MSU_ID'))
attr(RAPMSU, 'dateRetrieved') <- date()

####################
### GeneListById ###
####################

phytozome <- biomaRt::useMart('phytozome_mart', dataset = 'phytozome')
gramene <- biomaRt::useMart('ENSEMBL_MART_PLANT', dataset = 'osativa_eg_gene')

# get a cookie from shigen and extract the jsessionid
temp <- tempfile()
download.file('http://www.shigen.nig.ac.jp/rice/oryzabase/gene/download', temp, method = 'wget')
rl <- readLines(temp)
jsessionid <- gsub('.*(jsessionid=[0-9A-Z]+).*', "\\1", rl[grep('jsessionid', rl)][1])

# use the jsessionid to make a new url and download the file
fileUrl <- paste0("http://www.shigen.nig.ac.jp/rice/oryzabase/gene/download;", jsessionid, "?classtag=GENE_EN_LIST")
temp <- tempfile()
download.file(fileUrl, temp, method = 'wget')
GeneList <- read.delim(temp, stringsAsFactors = FALSE, na.strings = "_")

# mung the geneList
GeneList[is.na(GeneList)] <- ''
GeneListByID <- data.table::data.table(GeneList[GeneList$RAP_id != '',
                                                ][,c(7, 12, 13, 2:5, 18)
                                                  ],
                                       key = 'RAP_id')
GeneListByID[,RAP_ID := toupper(gsub(' ', '', RAP_id))]
GeneListByID[,RAP_id := NULL][,RAP_id := RAP_ID][,RAP_ID := NULL]
data.table::setkey(GeneListByID, RAP_id)

# concatenate the recommended gene symbols / names

pasteSymbols <- GeneListByID[,paste(recommended_gene_symbol, collapse = ','),
                             by = RAP_id
                             ][, recommended_gene_symbol := V1
                               ][, V1 := NULL]
pasteNames <- GeneListByID[,paste(recommended_gene_name, collapse = ','),
                           by = RAP_id][, recommended_gene_name := V1
                                        ][, V1 := NULL]
GeneListByID <- GeneListByID[, recommended_gene_symbol := NULL
                             ][, recommended_gene_name := NULL
                               ][pasteSymbols][pasteNames]
GeneListByID <- unique(GeneListByID)

# convert back to data.frame

GeneListByID.frame <- data.frame(GeneListByID)
rownames(GeneListByID.frame) <- GeneListByID.frame$RAP_id

attr(GeneListByID.frame, 'dateRetrieved') <- date()

#################
### Save data ###
#################

devtools::use_data(RAPMSU, GeneListByID.frame, internal = TRUE, overwrite = TRUE)
