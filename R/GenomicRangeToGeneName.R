#' GenomicRangeToGeneName
#' 
#' Return information for all of the genes in a genomic region.
#' 
#' @param seqnames A \code{chr} containing the sequence (\emph{i.e.} 
#'   chromosome) name. Must match the chromosome name in the MSUv7 GFF, 
#'   \emph{e.g.} chromosome 1 == \code{'Chr1'}
#' @param start A \code{numeric} containing the starting coordinate on 
#'   \code{seqnames}
#' @param end A \code{numeric} containing the end coordinat on \code{seqnames}
#' @return Returns a \code{data.table} of annotated genes found in the MSUv7 
#'   GFF between \code{start} and \code{end} on chromosome \code{seqnames}.
#'   
#' @seealso \link{\code{LocToGeneName}} for a description of the returned
#'   \code{data.table}.
#'   
#' @import data.table
#'   
#' @export

# seqnames <- "Chr1"
# start <- 34365
# end <- 55475

GenomicRangeToGeneName <- function(seqnames, start, end){
  my_range <- GenomicRanges::GRanges(
    seqnames = seqnames,
    ranges = IRanges::IRanges(
      start = start,
      end = end
    ),
    strand = "*",
    seqinfo = GenomicRanges::seqinfo(gff.exons))
  my_exons <- as.data.table(
    GenomicRanges::findOverlaps(gff.exons, my_range))
  my_hits <- as.data.table(
    gff.exons[my_exons[, queryHits]]
  )
  my_hits[, parent_transcript := unlist(Parent)]
  my_hits[, gene := gsub("\\..+$", "", parent_transcript)]
  my_hits[, LocToGeneName(unique(gene))]
}
