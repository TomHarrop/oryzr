#' LocToRefSeq
#' 
#' Convert MSU identifiers to Rap-DB identifiers and map them to RefSeq mRNA 
#' identifiers.
#' 
#' Conversion from MSU to Rap-DB identifiers uses the table of RAP-DB and MSU 
#' IDs \href{http://rapdb.dna.affrc.go.jp/download/irgsp1.html}{available from 
#' Rap-DB}. CSGNL symbols and gene names are retrieved from the 
#' OryzabaseGeneListEn database, downloaded from 
#' \href{http://www.shigen.nig.ac.jp/rice/oryzabase/download/gene}{Oryzabase}. 
#' Conversion to refSeq mRNA IDs uses the Rice ID IRGSP Build5 database 
#' downloaded from 
#' \href{http://www.shigen.nig.ac.jp/rice/oryzabase/download/riceId}{Oryzabase}.
#' If biomaRt is installed, it can be used for secondary searches for unmatched 
#' MSU IDs.
#' 
#' @param LOCs Character vector of MSU locus identifiers
#' @param useBiomart Boolean (default TRUE). If biomaRt is installed, use it for
#'   secondary searches for genes not matched in the Oryzabase data.
#' @return Returns a list of two \code{data.frame}s, matched for MSU loci that 
#'   were successfully matched to refSeq IDs and notMatched for ones that 
#'   weren't. The \code{data.frame}s each have four columns: tigrId (MSU ID),
#'   refseqRnaNucleotideAccessionNo (refSeq mRNA ID), rapLocus (Rap-DB ID) and
#'   symbols (CSGNL recommended gene symbol(s)).
#'   
#' @import data.table
#' @export
#' @examples 
#' LocToRefSeq(c('LOC_Os01g03340', 'LOC_Os01g09410', 'LOC_Os01g17396', 'LOC_Os01g18800', 'LOC_Os01g44069', 'LOC_Os01g45470', 'LOC_Os01g53880', 'LOC_Os01g64730', 'LOC_Os01g71820', 'LOC_Os02g07110', 'LOC_Os02g26680', 'LOC_Os02g26700', 'LOC_Os02g41904', 'LOC_Os02g42880', 'LOC_Os02g44630', 'LOC_Os03g15530', 'LOC_Os03g22270', 'LOC_Os03g61160', 'LOC_Os04g58710', 'LOC_Os06g06750', 'LOC_Os07g32170', 'LOC_Os08g24930', 'LOC_Os08g31980', 'LOC_Os09g07154', 'LOC_Os09g25330', 'LOC_Os09g31438', 'LOC_Os10g07229', 'LOC_Os10g18870', 'LOC_Os10g26340', 'LOC_Os10g42210', 'LOC_Os10g42220'))
#' LocToRefSeq(c('LOC_Os03g03034', 'LOC_Os06g36560', 'LOC_Os01g50910'))

LocToRefSeq <- function(LOCs, useBiomart = TRUE) {
    
    # search RapMsuRefSeq by LOC
    message(paste(dateF(), "Searching Oryzabase database for", length(LOCs), 
        "record(s)"))
    # data.table search
    if (requireNamespace("OpenRepGrid", quietly = TRUE)) {
        # sapply over LOCs with grep, unlist results and select columns to
        # output
        matchedRecords <- RapMsuRefSeq[unique(unlist(OpenRepGrid::sapply_pb(LOCs, 
            function(x) grep(x, RapMsuRefSeq$tigrId, ignore.case = TRUE)))), 
            ][, list(rapLocus, tigrId, refseqRnaNucleotideAccessionNo)]
    } else {
        # sapply over LOCs with grep, unlist results and select columns to
        # output
        matchedRecords <- RapMsuRefSeq[unique(unlist(sapply(LOCs, function(x) grep(x, 
            RapMsuRefSeq$tigrId, ignore.case = TRUE)))), ][, list(rapLocus, 
            tigrId, refseqRnaNucleotideAccessionNo)]
    }
    # remove trailing transcript numbers from tigrId and refSeq ID
    matchedRecords[, `:=`(tigrId, gsub("\\..*$", "", tigrId))][, `:=`(refseqRnaNucleotideAccessionNo, 
        gsub("\\..*$", "", refseqRnaNucleotideAccessionNo))]
    data.table::setkey(matchedRecords, tigrId, rapLocus)
    
    # search again in RAPMSU for unmatched LOC IDs
    unmatchedLocs <- matchedRecords[rapLocus == ""]$tigrId
    if (length(unmatchedLocs) > 0) {
        message(paste(dateF(), "Searching RAP-DB ID database for", length(unmatchedLocs), 
            "gene(s) not matched in Oryzabase"))
        if (requireNamespace("OpenRepGrid", quietly = TRUE)) {
            rematchedLocs <- RAPMSU[unique(unlist(OpenRepGrid::sapply_pb(unmatchedLocs, 
                function(x) grep(x, RAPMSU$MSU_ID, ignore.case = TRUE)))), 
                ]
        } else {
            rematchedLocs <- RAPMSU[unique(unlist(sapply(unmatchedLocs, 
                function(x) grep(x, RAPMSU$MSU_ID, ignore.case = TRUE)))), 
                ]
        }
        # remove trailng Tx numbers, remove unmatched locs
        rematchedLocs <- data.table(rematchedLocs)
        rematchedLocs <- rematchedLocs[, `:=`(tigrId, gsub("\\..*$", "", 
            MSU_ID))][, `:=`(rapLocus, toupper(Rap_ID))][, `:=`(Rap_ID, 
            NULL)][, `:=`(MSU_ID, NULL)][!rapLocus == "NONE"]
        data.table::setkey(rematchedLocs, tigrId)
        
        # re-join tables, remove NAs and duplicates
        rematchedRecords <- rbind(matchedRecords, rematchedLocs, fill = TRUE)
        rematchedRecords[is.na(refseqRnaNucleotideAccessionNo)]$refseqRnaNucleotideAccessionNo <- ""
        rematchedRecords <- unique(rematchedRecords, by = c("rapLocus", 
            "tigrId", "refseqRnaNucleotideAccessionNo"))
    } else {
        rematchedRecords <- unique(matchedRecords, by = c("rapLocus", "tigrId", 
            "refseqRnaNucleotideAccessionNo"))
    }
    # uppercase rapLocus and setkey for match with ensembl during biomaRt
    # search
    rematchedRecords[, `:=`(rapLocus, toupper(rapLocus))]
    data.table::setkey(rematchedRecords, rapLocus)
    
    unmatchedRecords <- rematchedRecords[refseqRnaNucleotideAccessionNo == 
        ""]
    # only search rapLocus ids that really don't already have a refSeqID
    unmatchedIds <- unmatchedRecords[!rapLocus %in% rematchedRecords[!refseqRnaNucleotideAccessionNo == 
        ""]$rapLocus]$rapLocus
    
    # biomaRt search for unmatched records
    if (useBiomart & length(unmatchedIds > 0)) {
        if (!requireNamespace("biomaRt", quietly = TRUE)) {
            stop("biomaRt is not installed. Either install it or use useBiomart == FALSE")
        }
        
        message("Preparing biomaRt object")
        ensembl <- biomaRt::useMart(biomart = "ENSEMBL_MART_PLANT", dataset = "osativa_eg_gene")
        
        # search BioMart
        message(paste(dateF(), "Running biomaRt query for", length(unmatchedIds), 
            "record(s) with no RefSeqID in Oryzabase"))
        refSeq <- data.table::data.table(biomaRt::getBM(attributes = c("ensembl_gene_id", 
            "refseq_mrna"), filters = "ensembl_gene_id", values = unique(unmatchedIds), 
            mart = ensembl), key = "ensembl_gene_id")
        # only process records that returned a hit
        bmMatch <- refSeq[!refseq_mrna == ""][, list(refseq_mrna = paste(refseq_mrna, 
            collapse = "/")), by = ensembl_gene_id]
        
        # join BM hits to rematchedRecords and remove NAs
        joinedResults <- bmMatch[rematchedRecords]
        joinedResults[is.na(refseq_mrna)]$refseq_mrna <- ""
        
    } else {
        joinedResults <- rematchedRecords
        joinedResults$refseq_mrna <- ""
        data.table::setnames(joinedResults, "rapLocus", "ensembl_gene_id")
    }
    message(paste(dateF(), "Merging results"))
    # fast join by tigrId to concatenate duplicates (tigrId with >1
    # rapLocus or refSeqID)
    
    # group refSeq IDs, don't group blanks; then group rapLocus IDs, don't
    # group blanks
    LocToRefSeq.table <- joinedResults[, list(refseqRnaNucleotideAccessionNo = list(unique(c(refseq_mrna[!refseq_mrna == 
        ""], refseqRnaNucleotideAccessionNo[!refseqRnaNucleotideAccessionNo == 
        ""]))), rapLocus = list(unique(ensembl_gene_id[!ensembl_gene_id == 
        ""]))), by = tigrId]
    
    # add short gene names
    message(paste(dateF(), "Adding short gene names from Oryzabase"))
    LocToRefSeq.table[, `:=`(symbols, LocToGeneName(LocToRefSeq.table$tigrId)$symbols)]
    # remove NAs
    LocToRefSeq.table$symbols[is.na(LocToRefSeq.table$symbols)] <- ""
    
    # swap commas for slashes
    swapFunction <- function(x) {
        sapply(x, paste, collapse = "/")
    }
    # rapLocus and RefSeqID are lists of lists, so can use paste to
    # concatenate and convert to character vector simultaneously, but
    # symbols is already a character vector, so use gsub instead
    LocToRefSeq.table[, `:=`(rapLocus, swapFunction(rapLocus))][, `:=`(refseqRnaNucleotideAccessionNo, 
        swapFunction(refseqRnaNucleotideAccessionNo))][, `:=`(symbols, 
        gsub(",", "/", symbols))]
    setkey(LocToRefSeq.table, tigrId)
    # split into found / not found and return
    
    matched <- LocToRefSeq.table[!refseqRnaNucleotideAccessionNo == ""]
    notMatched <- LocToRefSeq.table[refseqRnaNucleotideAccessionNo == ""]
    
    message(paste(dateF(), "Done. Found refSeqIDs for", dim(matched)[1], 
        "records"))
    
    return(list(matched = matched, notMatched = notMatched))
} 
