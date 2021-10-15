# R Analysis for scRNAseq from Young/Old mice

## Data

The data comes from 6 old and 6 young mice. The mice genotypes and tissue area can be 
found in the metadata file metadata_scRNAseq.csv.

The raw reads have been processed with CellRanger version 5. The resulting 
filtered matrices are in the "data" folder.

## Src
 In source I add the code used to do the analysis. 
 Most of the code and the explanations for the use of different method are extracted from the [Sanger institute scRNA-seq workflow](https://scrnaseq-course.cog.sanger.ac.uk/website/index.html) and the Bioconductor book, [Orchestrating Single-Cell Analysis with Bioconductor](https://bioconductor.org/books/release/OSCA/)
 
-   The first cell and gene quality control is in QC

-   The normalisation by deconvolution is in normalise

-   Feature selection and dimensional reduction in feature_selection_dimred

-   Clustering at different resolutions, clustering_01

-   First rough annotation in annotation_01

-   Cluster QC in clusterQC_k5, now includes doublets detection

-   Stricter cell and gene QC; dimensional reduction and feature selection with the cleaned data

-   Clustering with different resolutions clustering_02

-   Annotation annotation_02

-   Differential gene expression between WT and KO for each celltype DE_WT_KO_celltype and each cluster DE_WT_KO_k20

-   Differential abundance of cells between clusters between WT and KO DA_WT_KO

-   Cell cycle in cell_cycle

-   Subset only the oligodendroglia in oligos

Interactive [Shinyapp](https://annawilliams.shinyapps.io/shinyApp_jpriller) with the analysis until this stage.

 
 
## Docs
 The html documents generated with the source code and the plots are separated 
 in this folder. (no echo when knitting, to hide the code for a nicer view)