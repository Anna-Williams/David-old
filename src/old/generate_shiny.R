library(SingleCellExperiment)
library(ShinyCell)
library(here)

age <- "old" 

sce <- readRDS(here("processed", age, "sce_anno_02.RDS"))

# remove unwanted dimensional reductions
 reducedDim(sce, "PCA_all") <- NULL
 reducedDim(sce, "PCA_coldata") <- NULL


conf <- createConfig(sce)
#Delete some unnecessary metadata
conf <- delMeta(conf, meta.to.del = c( "subsets_mt_sum", "subsets_mt_detected",
                       "outlier", "ratio_detected_sum", "outlier_ratio",
                       "discard", "sizeFactor", "total", "subsets_ribo_sum",
                       "subsets_ribo_detected", "filter_out", "original_sample_name",
                       "dbl_exploratory")
                )
# Change name of some metadata
conf <- modMetaName(conf,
                    meta.to.mod = c("sum", "detected", "subsets_mt_percent", "subsets_ribo_percent",
                                    "celltype"),
                    new.name = c("umi counts", "detected genes", "% mt genes", "% ribo genes",
                                 "cell type"))
# 

makeShinyApp(sce, conf, gene.mapping = TRUE, 
             shiny.title = "ShinyCell app",
             default.dimred = c("TSNE1", "TSNE2"))


