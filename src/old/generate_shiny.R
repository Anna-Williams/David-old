library(SingleCellExperiment)
library(ShinyCell)
library(here)
library(DropletUtils)

age <- "old" 

sce <- readRDS(here("processed", age, "sce_anno.RDS"))
#[provisional until I correct all the objects with the import.r modification]
 barcodes_mod <- paste0("1_", sce$Barcode)
 colnames(sce)<- make.unique(barcodes_mod)
# # test with seurat
 seu <- as.Seurat(sce)


conf <- createConfig(seu)
# #Delete some unnecessary metadata
# conf <- delMeta(conf, meta.to.del = c( "subsets_mt_sum", "subsets_mt_detected",
#                        "outlier", "ratio_detected_sum", "outlier_ratio",
#                        "discard", "sizeFactor", "total")
#                 )
# # Change name of some metadata
# conf <- modMetaName(conf,
#                     meta.to.mod = c("sum", "detected", "subsets_mt_percent", "celltype"),
#                     new.name = c("umi counts", "detected genes", "% mt genes", "cell type"))
# 

makeShinyApp(seu, conf, gene.mapping = TRUE, shiny.title = "ShinyCell app")
