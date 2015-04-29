oryzr
=====

Work with *Oryza sativa* gene identifiers
-----------------------------------------

This R package is for converting between different types of gene identifiers for *Oryza sativa* L. ssp. *japonica* cv. Nipponbare genes (MSU, Rap-DB, CSGNL, RefSeq mRNA). Currently, only conversions from MSU ID are implemented.

#### Functions

-   LocToGeneName: Convert MSU identifiers to CSGNL gene names and symbols and generate a label for plotting if requested.
-   LocToRefSeq: Convert MSU identifiers to Rap-DB identifiers and map them to RefSeq mRNA identifiers.

#### Installation

-   Make sure you have [Hadley Wickham's `devtools`](https://github.com/hadley/devtools) installed: `install.packages("devtools")`
-   Install `oryzr` from GitHub: `devtools::install_github('Tom14-/oryzr')`
