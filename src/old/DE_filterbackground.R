
library(scater) # for aggregate counts
library(edgeR) #for De
library(here) # reproducible paths
library(DropletUtils) # ambient RNA
library(scuttle) # modify gene names
library(scran) # DE

age <- "old"
project <- "old"
sce <- readRDS(here("processed", age, "sce_anno_02.RDS"))

# Calculate ambient files
if (!file.exists(here("processed", project, "ambient.RDS"))) {
  # import data
  # to calculate the ambient the raw data should be used
  samples <- dir(here("data/raw/"))
  # create empty vector to hold ambient
  ambient <- vector("list", length(samples))
  names(ambient) <- samples
  # generate profile for each sample
  for(sample in samples){
    matrix <- here("data/raw/", sample, "raw_feature_bc_matrix")
    sce_sample <- read10xCounts(matrix, sample, version = "auto", col.names = TRUE)
    # calculate ambient counts, with good turning FALSE because we will use the counts for further estimations and round false because we already have integers. 
    ambient[[sample]] <- ambientProfileEmpty(counts(sce_sample), good.turing = FALSE, round = FALSE)
  }
  
  # clean up ambient output
  ambient <-  do.call(cbind, ambient)
  
  # use symbol nomenclature
  rownames(ambient) <- uniquifyFeatureNames(
    ID=rowData(sce_sample)$ID,
    names=rowData(sce_sample)$Symbol
  )
  
  saveRDS(ambient, (here("processed", project, "ambient.RDS")))
} else {
  ambient <- readRDS((here("processed", project, "ambient.RDS")))
}


#### add Ambient to the k20 DE
de_results <- readRDS(here("processed",age, "DE_k20_de_results.RDS"))


# create summed object
summed <- aggregateAcrossCells(sce, id=colData(sce)[,c("clusters_named", "Sample")])

# Removing all pseudo-bulk samples with 'insufficient' cells.
summed_filt <- summed[, summed$ncells >= 10]


# maximum proportion for each gene that could be due to ambient
# I need here to sort out the ambient vecotr, so I can do all clusters at once, repeating n (often 6) times
# the ambient profile. the ambient profile needs to be calculated on groups of cells that are similar: We'll use 
# the clusters we run the DEG on, like this it can be added to each of these. It also needs to run per sample (obviously you don't calculate the ambient with another sample), but we can take for each gene the avareage of the proportion obtained in each sample. For this again we will need to sort out the fact I'm doing this for all clusters at once, as we only want to average for each cluster. 

# create an aggregate ambient, duplicating the information to match the format of the summed sce
n=0
aggregate_ambient <- matrix(nrow = nrow(ambient), ncol=ncol(summed_filt))
for( sample in summed_filt$Sample){
  n = n+1
  aggregate_ambient[,n] <-  ambient[,sample]
}
# add the gene names again
rownames(aggregate_ambient) <- rownames(ambient)

# subset for same genes as in the summed sce
aggregate_ambient <- aggregate_ambient[row.names(summed_filt),]

# upper bound of gene contribution from ambient contamination
max_ambient <- ambientContribMaximum(counts(summed_filt), ambient= aggregate_ambient, mode = "proportion")
microglia_genes <- readLines(here("data", "microglia_genes.csv"))
microglia_ambient <- ambientContribNegative(counts(summed_filt), ambient = aggregate_ambient, features =  microglia_genes, mode = "proportion")


# mean of result obtained between all samples -> needs to be done for each cluster

# create a list to hold the updated DFs ( i tried with lapply but not there yet)

# add to the list of dataframes with results for each cluster the ambient info

de_result_ambient <- lapply(names(de_results), function(cluster) {
  de_results_clu <- de_results[[cluster]]
  # calculate mean of all samples for that cluster
  cluster_max_ambient <- rowMeans(max_ambient[, summed_filt$clusters_named == cluster], na.rm = TRUE)
  # for the ambient calculated from microglia genes exlude the immune and microglia clusters
  if(cluster %in% c("Microglia", "Immune")){
    cluster_microglia_ambient <- NA
  }else{
    cluster_microglia_ambient <- rowMeans(microglia_ambient[, summed_filt$clusters_named == cluster], na.rm = TRUE)
  }
  # combine with previous dataframe
  cbind(de_results_clu, maxAmbient = cluster_max_ambient, microgliaAmbient = cluster_microglia_ambient, minAmbient = pmin(cluster_max_ambient, cluster_microglia_ambient))
})
names(de_result_ambient) <- names(de_results)

#save R objects
dir.create(here("processed",project), showWarnings = FALSE)
saveRDS(de_result_ambient, here("processed",project, "DE_results_ambient_edgeR_k20.RDS"))


# save results DE tables, different files for each cluster
lapply(names(de_result_ambient), function(cluster){
  de_results_clu <- de_result_ambient[[cluster]]
  write.csv(de_results_clu, here("outs", age, "DE_k20", "de_results", paste0("de_results_", cluster, ".csv")), quote = FALSE)
}
)

### Repeat Add ambient to CELLTYPE comparisons

# create summed object
summed <- aggregateAcrossCells(sce, id=colData(sce)[,c("celltype", "Sample")])

# Removing all pseudo-bulk samples with 'insufficient' cells.
summed_filt <- summed[, summed$ncells >= 10]


# maximum proportion for each gene that could be due to ambient
# I need here to sort out the ambient vecotr, so I can do all clusters at once, repeating n (often 6) times
# the ambient profile. the ambient profile needs to be calculated on groups of cells that are similar: We'll use 
# the clusters we run the DEG on, like this it can be added to each of these. It also needs to run per sample (obviously you don't calculate the ambient with another sample), but we can take for each gene the avareage of the proportion obtained in each sample. For this again we will need to sort out the fact I'm doing this for all clusters at once, as we only want to average for each cluster. 

# create an aggregate ambient, duplicating the information to match the format of the summed sce
n=0
aggregate_ambient <- matrix(nrow = nrow(ambient), ncol=ncol(summed_filt))
for( sample in summed_filt$Sample){
  n = n+1
  aggregate_ambient[,n] <-  ambient[,sample]
}
# add the gene names again
rownames(aggregate_ambient) <- rownames(ambient)

# subset for same genes as in the summed sce
aggregate_ambient <- aggregate_ambient[row.names(summed_filt),]

# upper bound of gene contribution from ambient contamination
max_ambient <- ambientContribMaximum(counts(summed_filt), ambient= aggregate_ambient, mode = "proportion")
microglia_genes <- readLines(here("data", "microglia_genes.csv"))
microglia_ambient <- ambientContribNegative(counts(summed_filt), ambient = aggregate_ambient, features =  microglia_genes, mode = "proportion")


# mean of result obtained between all samples -> needs to be done for each cluster

# create a list to hold the updated DFs ( i tried with lapply but not there yet)

# add to the list of dataframes with results for each cluster the ambient info
#de_results <- readRDS(here("processed",age, "DE_k20_de_results.RDS"))

# I did not save the celltype one, repeat the analysis 
summed_filt$genotype <- factor(summed_filt$genotype, levels = c("WT", "KO"))
de_results <- pseudoBulkDGE(summed_filt, 
                            label=summed_filt$celltype,
                            design= ~factor(chip) + genotype,
                            coef="genotypeKO",
                            condition=summed_filt$genotype 
)

# Here add the list to the results
de_result_ambient <- lapply(names(de_results), function(cluster) {
  de_results_clu <- de_results[[cluster]]
  # calculate mean of all samples for that cluster
  cluster_max_ambient <- rowMeans(max_ambient[, summed_filt$celltype == cluster], na.rm = TRUE)
  # for the ambient calculated from microglia genes exlude the immune and microglia clusters
  if(cluster %in% c("Microglia", "Immune")){
    cluster_microglia_ambient <- NA
  }else{
    cluster_microglia_ambient <- rowMeans(microglia_ambient[, summed_filt$celltype == cluster], na.rm = TRUE)
  }
  # combine with previous dataframe
  cbind(de_results_clu, maxAmbient = cluster_max_ambient, microgliaAmbient = cluster_microglia_ambient, minAmbient = pmin(cluster_max_ambient, cluster_microglia_ambient))
})
names(de_result_ambient) <- names(de_results)

#save R objects
dir.create(here("processed",project), showWarnings = FALSE)
saveRDS(de_result_ambient, here("processed",project, "DE_results_ambient_edgeR_celltype.RDS"))


# save results DE tables, different files for each cluster
lapply(names(de_result_ambient), function(cluster){
  de_results_clu <- de_result_ambient[[cluster]]
  write.csv(de_results_clu, here("outs", age, "DE_celltype", "de_results", paste0("de_results_", cluster, ".csv")), quote = FALSE)
}
)
