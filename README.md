# R Analysis for scRNAseq from Young/Old mice

## Data

The data comes from 6 old and 6 young mice. The mice genotypes and tissue area can be 
found in the metadata file metadata_scRNAseq.csv.

The raw reads have been processed with CellRanger version 5. The resulting 
filtered matrices are in the "data" folder.

## Src
 In source I add the code used to do the analysis. 
 Most of the code and the explanations for the use of different method are extracted from the [Sanger institute scRNA-seq workflow](https://scrnaseq-course.cog.sanger.ac.uk/website/index.html) and the Bioconductor book, [Orchestrating Single-Cell Analysis with Bioconductor](https://bioconductor.org/books/release/OSCA/)
 
 - import.R imports the raw counts from 10x and the metadata into a sce object. 
 - The first cell and gene quality control is in QC.Rmd
 
 
## Docs
 The html documents generated with the source code and the plots are separated 
 in this folder.