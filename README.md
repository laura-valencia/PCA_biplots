This script contains two function used to produce PCA biplots from the output of a DESeq2 (R) run on RNA-seq count data.
The functions include 1) the obtention of a dds object to extract counts in an R-readable way, and 2) parsing a DESeq2 output file with differentially expressed genes, next to the previously obtained counts, to obtain a table with the top genes contributing to variation.
The function2 also produces the PCA biplot.

