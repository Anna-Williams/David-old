## Second cell and gene QC

# set up ----
# load packages
library(here)
library(scran) # for normalisation

age <- "old"

sce <- readRDS(here("processed/old/sce_clusterQC.RDS"))

## Cell QC ---

# set stricter thresholds
discard_lib_high <- sce$sum > 20000
discard_expr_high <- sce$detected > 50000
discard_expr_low <- sce$sum < 500
discard_lib_low <- sce$detected < 300
discard_mt <- sce$subsets_mt_percent > 15

sce$discard <-  discard_expr_low | discard_lib_low | discard_mt

summary_discard <- data.frame(
  lib_size_high = c(sum(discard_lib_high), 20000),
  expression_high = c(sum(discard_expr_high), 50000),
  lib_size_low = c(sum(discard_lib_low), 300),
  expression_low = c(sum(discard_expr_low), 500),
  mt_pct = c(sum(discard_mt), 15),
  total = c(sum(sce$discard), NA)
)

# very few cells are being deleted by high expression thresholds, doublets have been removed succesfully
write.csv(summary_discard,
          here("outs/old/reason_for_discard_QC_02.csv"))

# subset the cells
sce <- sce[, sce$discard == FALSE]
## Gene QC ---

# at least 10 cells should express the gene
keep_feature <- rowSums(counts(sce) > 0) > 10
sce <- sce[keep_feature,]


## Save object with QC ---
saveRDS(sce, here("processed/old/sce_QC_02.RDS"))

## Normalise by deconvolution ---

# For reproducibility
set.seed(100)
# Quick clustering to pool samples together and deal with 0 counts
quick_clusters <- quickCluster(sce)
# Calculate size factors
sce <-
  computeSumFactors(sce, cluster = quick_clusters, min.mean = 0.1)
# Check that there are not negative size factors
summary(sizeFactors(sce))
# Apply size factors and log transform them
sce <- logNormCounts(sce)
# save object
saveRDS(sce, here("processed", age,  "sce_norm_02.RDS"))


