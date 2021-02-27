# R Analysis for scRNAseq from Young/Old mice

## Data

The data comes from 6 old and 6 young mice. The mice genotypes and tissue area can be 
found in the metadata file metadata_scRNAseq.csv.

The raw reads have been processed with CellRanger version 5. The resulting 
filtered matrices are in the "data" folder.

## Src
 In source I add the code used to do the analysis. 
 - import.R imports the raw counts from 10x and the metadata into a sce object. 
