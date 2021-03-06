oryzr
=====

Work with *Oryza sativa* gene identifiers
-----------------------------------------

This R package is for converting between different types of gene identifiers for *Oryza sativa* L. ssp. *japonica* cv. Nipponbare genes (MSU, Rap-DB, CSGNL, RefSeq mRNA).

#### Functions

-   LocToGeneName: Convert MSU identifiers to CSGNL gene names and symbols and generate a label for plotting if requested.
-   LocToRefSeq: Convert MSU identifiers to Rap-DB identifiers and map them to RefSeq mRNA identifiers.
-   SearchByGeneSymbol: Search Rap-DB and oryzabase information by gene symbol and return identifiers for potential matches.
-   RapToMsu: Convert Rap-DB identifiers to MSU identifiers.
-   GenomicRangeToGeneName: Return information for all of the genes in a genomic region.

#### Installation

-   Requires [`devtools`](https://github.com/hadley/devtools): `install.packages("devtools")`
-   Install `oryzr` from GitHub: `devtools::install_github("TomHarrop/oryzr")`
