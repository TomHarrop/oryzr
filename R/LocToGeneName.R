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
#' @param plotLabels Boolean (default \code{TRUE}). Should labels for plotting be 
#'   generated? Labels will be of the format LOC_ID\\nGene name(s) unless 
#'   \code{shortLabels == TRUE}.
#' @param shortLabels Boolean (default \code{FALSE}). Return labels in format 
#'   LOC_ID\\nGene symbols(s). Has no effect if \code{plotLabels == FALSE}.
#' @param simpleLabels Boolean (default \code{FALSE}). Return \emph{really} short 
#'   labels (CSGNL symbol or MSU ID \strong{only}). Useful \emph{e.g}. in 
#'   heatmaps with lots of genes. Has no effect if \code{plotLabels == FALSE}
#'   and is overriden by \code{shortLabels == TRUE}.
#' @return Returns a \code{data.frame} with columns RapID (RAP-DB gene identifier), 
#'   symbols (CSGNL recommended symbols) and names (CSGNL recommended names) and
#'   optionally labels for plotting.
#'   
#' @export
#' @examples 
#' LocToGeneName(c('LOC_Os01g03340', 'LOC_Os01g09410', 'LOC_Os01g17396', 'LOC_Os01g18800', 'LOC_Os01g44069', 'LOC_Os01g45470', 'LOC_Os01g53880', 'LOC_Os01g64730', 'LOC_Os01g71820', 'LOC_Os02g07110', 'LOC_Os02g26680', 'LOC_Os02g26700', 'LOC_Os02g41904', 'LOC_Os02g42880', 'LOC_Os02g44630', 'LOC_Os03g15530', 'LOC_Os03g22270', 'LOC_Os03g61160', 'LOC_Os04g58710', 'LOC_Os06g06750', 'LOC_Os07g32170', 'LOC_Os08g24930', 'LOC_Os08g31980', 'LOC_Os09g07154', 'LOC_Os09g25330', 'LOC_Os09g31438', 'LOC_Os10g07229', 'LOC_Os10g18870', 'LOC_Os10g26340', 'LOC_Os10g42210', 'LOC_Os10g42220'))
#' LocToGeneName(c('LOC_Os03g03034', 'LOC_Os06g36560', 'LOC_Os01g50910'))

LocToGeneName <- function(LOCs, plotLabels = TRUE, shortLabels = FALSE, 
    simpleLabels = FALSE) {
    
    # this needs to be fine tuned. Currently, if one MSU ID matches more
    # than one RAP-DB id, only the first one is used.
    if (requireNamespace("OpenRepGrid", quietly = TRUE)) {
        message("Matching MSU IDs to Rap-DB IDs. This could take some time.")
        Rap.vector <- OpenRepGrid::sapply_pb(LOCs, function(x) RAPMSU[grepl(x, 
            RAPMSU$MSU_ID, ignore.case = TRUE), ][1, ]$Rap_ID)
    } else {
        message("Matching MSU IDs to Rap-DB IDs. This could take some time. If you would like to see a progress bar please install OpenRepGrid from CRAN.")
        Rap.vector <- sapply(LOCs, function(x) RAPMSU[grepl(x, RAPMSU$MSU_ID, 
            ignore.case = TRUE), ][1, ]$Rap_ID)
    }
    Rap.vector[Rap.vector == "None"] <- NA
    
    # merge with GeneListByID
    symbols <- sapply(Rap.vector, function(x) GeneListByID.frame[toupper(x), 
        ]$recommended_gene_symbol)
    symbols[symbols == ""] <- NA
    names <- sapply(Rap.vector, function(x) GeneListByID.frame[toupper(x), 
        ]$recommended_gene_name)
    names[names == ""] <- NA
    LocToLabels <- data.frame(row.names = names(Rap.vector), RapID = Rap.vector, 
        symbols = symbols, names = names, stringsAsFactors = FALSE)
    
    # deal with genes that only have synonyms on CSGNL (which would return
    # an empty field with a comma, normally)
    LocToLabels$symbols[!grepl("[A-Za-z]", LocToLabels$symbols)] <- NA
    LocToLabels$names[!grepl("[A-Za-z]", LocToLabels$names)] <- NA
    
    # generate an optional labels field for graphing
    if (plotLabels & !shortLabels & !simpleLabels) {
        labels <- is.na(LocToLabels$symbols) & is.na(LocToLabels$names)
        LocToLabels$labels <- labels
        for (i in 1:dim(LocToLabels)[1]) {
            if (LocToLabels[i, ]$labels) {
                LocToLabels[i, ]$labels <- rownames(LocToLabels)[i]
            } else {
                LocToLabels[i, ]$labels <- paste(gsub("[/,]", "/\n", LocToLabels[i, 
                  ]$names), paste("(", LocToLabels[i, ]$symbols, ")", sep = ""), 
                  rownames(LocToLabels)[i], sep = "\n")
            }
        }
    }
    if (shortLabels & plotLabels) {
        labels <- is.na(LocToLabels$symbols) & is.na(LocToLabels$names)
        LocToLabels$labels <- labels
        for (i in 1:dim(LocToLabels)[1]) {
            if (LocToLabels[i, ]$labels) {
                LocToLabels[i, ]$labels <- rownames(LocToLabels)[i]
            } else {
                LocToLabels[i, ]$labels <- paste(LocToLabels[i, ]$symbols, 
                  rownames(LocToLabels)[i], sep = "\n")
            }
        }
    }
    if (plotLabels & simpleLabels & !shortLabels) {
        LocToLabels$labels <- LocToLabels$symbols
        LocToLabels$labels[is.na(LocToLabels$labels)] <- rownames(LocToLabels)[is.na(LocToLabels$labels)]
    }
    
    return(LocToLabels)
    
} 
