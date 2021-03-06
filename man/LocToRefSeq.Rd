% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/LocToRefSeq.R
\name{LocToRefSeq}
\alias{LocToRefSeq}
\title{LocToRefSeq}
\usage{
LocToRefSeq(LOCs, useBiomart = TRUE)
}
\arguments{
\item{LOCs}{Character vector of MSU locus identifiers}

\item{useBiomart}{Boolean (default \code{TRUE}). If biomaRt is installed, use it for
secondary searches for genes not matched in the Oryzabase data.}
}
\value{
Returns a list of two \code{data.frame}s, matched for MSU loci that 
  were successfully matched to refSeq IDs and notMatched for ones that 
  weren't. The \code{data.frame}s each have four columns: tigrId (MSU ID), 
  refseq_mrna (refSeq mRNA ID), ensembl_gene_id (Rap-DB ID) and 
  symbols (CSGNL recommended gene symbol(s)).
}
\description{
Convert MSU identifiers to Rap-DB identifiers and map them to RefSeq mRNA 
identifiers. Currently only conversions \strong{from} MSU ID are implemented.
}
\details{
Conversion from MSU to Rap-DB identifiers uses the table of RAP-DB and MSU 
IDs \href{http://rapdb.dna.affrc.go.jp/download/irgsp1.html}{available from 
Rap-DB}. CSGNL symbols and gene names are retrieved from the 
OryzabaseGeneListEn database, downloaded from 
\href{http://www.shigen.nig.ac.jp/rice/oryzabase/download/gene}{Oryzabase}. 
Conversion to refSeq mRNA IDs uses the Rice ID IRGSP Build5 database 
downloaded from 
\href{http://www.shigen.nig.ac.jp/rice/oryzabase/download/riceId}{Oryzabase}.
If biomaRt is installed, it can be used for secondary searches for unmatched 
MSU IDs.
}
\examples{
LocToRefSeq(c('LOC_Os01g03340', 'LOC_Os01g09410', 'LOC_Os01g17396', 'LOC_Os01g18800', 'LOC_Os01g44069', 'LOC_Os01g45470', 'LOC_Os01g53880', 'LOC_Os01g64730', 'LOC_Os01g71820', 'LOC_Os02g07110', 'LOC_Os02g26680', 'LOC_Os02g26700', 'LOC_Os02g41904', 'LOC_Os02g42880', 'LOC_Os02g44630', 'LOC_Os03g15530', 'LOC_Os03g22270', 'LOC_Os03g61160', 'LOC_Os04g58710', 'LOC_Os06g06750', 'LOC_Os07g32170', 'LOC_Os08g24930', 'LOC_Os08g31980', 'LOC_Os09g07154', 'LOC_Os09g25330', 'LOC_Os09g31438', 'LOC_Os10g07229', 'LOC_Os10g18870', 'LOC_Os10g26340', 'LOC_Os10g42210', 'LOC_Os10g42220'))
LocToRefSeq(c('LOC_Os03g03034', 'LOC_Os06g36560', 'LOC_Os01g50910'))
}
