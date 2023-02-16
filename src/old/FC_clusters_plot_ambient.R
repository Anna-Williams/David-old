# SCRIPT FOR PLOT WITH FC ON Y AXIS, GROUPS IN X AXIS

## set up ---
library(ggplot2) # manipulate dfs and ggplots
library(tibble)
library(forcats)
library(dplyr)
library(here) # reproducible paths

# functions

plotFCbyCluster <- function(de_results, # scran pseudoBulkDGE output
                              ambient = FALSE,
                              ambient_threshold = 0.1,
                              FDR_threshold = 0.05,
                              logFC_threshold = 0) {
  
  # plots each gene organised on the X axis by clusters and representing the 
  # FC on the y axis and the colour if they are significant
  
  
  # transform all the DFrame to a standart df format
  de_results_dfs <- lapply(
    de_results,
    function(x) {
      na.omit(as.data.frame(x)) %>%
        rownames_to_column("gene")
    }
  )
  # merge all the dfs from the list in a single big dataframe
  de_results_df <- bind_rows(de_results_dfs, .id = "cluster") %>%
    # add column that indicates if the result is significant
    mutate(Significant = ifelse(FDR < FDR_threshold,
                                yes = ifelse(logFC > logFC_threshold, "upregulated", "downregulated"),
                                no = "not-significant")) %>% 
    mutate(logFC_colour = ifelse(Significant == "not-significant", NA, logFC))  %>%
    # order the levels to display in correct order in plot
    mutate(cluster = forcats::fct_relevel(cluster, unique(cluster)))  %>% 
    # sort so the significant values are plotted on top of the non significants
    mutate(Significant = forcats::fct_relevel(Significant, c("upregulated", "downregulated", "not-significant"))) %>%
    arrange(desc(Significant))
  
  if (ambient == TRUE) {
    # filter the genes that are more than 10% ambient (or other threshold)
    de_results_df <- de_results_df %>%
      filter(minAmbient < ambient_threshold)
  }

    ## plot ---
    de_results_df %>%
      ggplot(mapping = aes(x = cluster, y = logFC, color = logFC_colour)) +
      geom_jitter() +
      scale_colour_gradient2(high = scales::muted("red"), mid = "white", low = scales::muted("blue")) +
      theme_minimal() +
      theme(
        axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1),
        legend.title = element_blank()
      ) +
      ylab("Log Fold Change") +
      xlab(element_blank())
}

plotDEbyCluster <- function(de_results, # scran pseudoBulkDGE output
                            ambient = FALSE,
                            ambient_threshold = 0.1,
                            FDR_threshold = 0.05,
                            logFC_threshold = 0,
                            data = FALSE) {
  
  # plots each gene organised on the X axis by clusters and representing the 
  # FC on the y axis 
  # and the FDR with the colour gradient: 
  # smaller FDR have more intense colour (red if up and blue if down), non-significant genes are grey. 

  # transform all the DFrame to a standart df format
  de_results_dfs <- lapply(
    de_results,
    function(x) {
      na.omit(as.data.frame(x)) %>%
        rownames_to_column("gene")
    }
  )
  # merge all the dfs from the list in a single big dataframe
  de_results_df <- bind_rows(de_results_dfs, .id = "cluster") %>%
    # add column that indicates if the result is significant
    mutate(Significant = ifelse(FDR < FDR_threshold,
      yes = ifelse(logFC > logFC_threshold, "upregulated", "downregulated"),
      no = "not-significant"
    )) %>%
    mutate(FDR_colour = ifelse(Significant == "not-significant", NA,
      ifelse(Significant == "downregulated", 0 - (FDR_threshold - FDR), FDR_threshold - FDR)
    ))  %>%
    # order the levels to display in correct order in plot
    mutate(cluster = forcats::fct_relevel(cluster, unique(cluster)))  %>% 
    # sort so the significant values are plotted on top of the non significants
    mutate(Significant = forcats::fct_relevel(Significant, c("upregulated", "downregulated", "not-significant"))) %>%
    arrange(desc(Significant))

  if (ambient == TRUE) {
    # filter the genes that are more than 10% ambient (or other threshold)
    de_results_df <- de_results_df %>%
      filter(minAmbient < ambient_threshold)
  }
  
  ## plot ---
    # set up legend colours
    max_col <- max(na.omit(de_results_df$FDR_colour))
    min_col <- min(na.omit(de_results_df$FDR_colour))
    ## plot ---
    de_results_df %>%
      ggplot(mapping = aes(x = cluster, y = logFC, color = FDR_colour)) +
      geom_jitter() +
      scale_colour_gradientn(
        colours = pals::coolwarm(),
        breaks = c(max_col, 0, min_col),
        labels = c(round(FDR_threshold - max_col, 2), FDR_threshold, round(FDR_threshold + min_col, 2)),
        name = "adj. p-value"
      ) +
      theme_minimal() +
      theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1)) +
      ylab("Log Fold Change") +
      xlab(element_blank())

}

filterDE <- function(de_results,
                     ambient = FALSE,
                     ambient_threshold = 0.1,
                     FDR_threshold = 0.05,
                     logFC_threshold = 0) {
  
  # filters the scran pseudoBulkDGE output, returning a combined df with all the clusters
  # appended
  
  
  # transform all the DFrame to a standart df format
  de_results_dfs <- lapply(
    de_results,
    function(x) {
      na.omit(as.data.frame(x)) %>%
        rownames_to_column("gene")
    }
  )
  # merge all the dfs from the list in a single big dataframe
  de_results_df <- bind_rows(de_results_dfs, .id = "cluster") %>%
    # add column that indicates if the result is significant
    mutate(Significant = ifelse(FDR < FDR_threshold,
                                yes = ifelse(logFC > logFC_threshold, "upregulated", "downregulated"),
                                no = "not-significant")) %>% 
    filter(Significant != "not-significant") %>% 
    mutate(cluster = forcats::fct_relevel(cluster, unique(cluster)))
  
  
  if (ambient == TRUE) {
    # filter the genes that are more than 10% ambient (or other threshold)
    de_results_df <- de_results_df %>%
      filter(minAmbient < ambient_threshold)
  }
}

getAmbientGenes <- function(de_results,
                              ambient_threshold = 0.1,
                              FDR_threshold = 0.05) {
  
  # To know which genes are filtered out with a certain ambient threshold: extracts from the scran pseudoBulkDGE output, 
  # returning a combined df with all the clusters appended
  
  # transform all the DFrame to a standart df format
  de_results_dfs <- lapply(
    de_results,
    function(x) {
      na.omit(as.data.frame(x)) %>%
        rownames_to_column("gene")
    }
  )
  # merge all the dfs from the list in a single big dataframe
  de_results_df <- bind_rows(de_results_dfs, .id = "cluster") %>%
    filter(minAmbient > ambient_threshold) %>%
    # add column that indicates if the result is significant
    mutate(Significant = ifelse(FDR < FDR_threshold,
                                yes = ifelse(logFC > 0, "upregulated", "downregulated"),
                                no = "not-significant"
    )) %>%
    mutate(cluster = forcats::fct_relevel(cluster, names(de_results))) %>% 
    # sort so the significant values are plotted on top of the non significants
    arrange(desc(Significant))
}


# project
project <- "old"

de_results <- readRDS(here("processed",project, "DE_results_ambient_edgeR_k20.RDS"))

# the output from edgeR for sc is
# a list with the DE results, each element named as one of the clusters and
# contain a DFrame with the DE for that cluster

# sort the order from the list # use levels from sce clsuters
de_results <-
  de_results[c("BAMs", "Lymphocytes","DCs","Astrocyte_1", "Astrocyte_2", "Astrocyte_3", "OligoAstro",  "OPCs","Oligo_1", "Oligo_2",  "Endothelial","fEndothelia","Mural_cells",  "ChP_epithelial", "iNeurons & NRPs", "mNeurons")]


plotDEbyCluster(de_results, ambient=TRUE, ambient_threshold = 0.25) 
dir.create(here("outs", project, "DE_k20", "plots"))
ggsave(here("outs", project, "DE_k20", "plots", "DE_clusters.pdf"),
       height = 7, width = 10
)

## For the celltypes ###
#load
de_results <- readRDS(here("processed",project, "DE_results_ambient_edgeR_celltype.RDS"))
#sort
de_results <-  de_results[c( "BAMs", "Lymphocytes","DCs","Astrocyte", "OligoAstro",  "OPCs","Oligo",  "Endothelial","Mural_cells", "ChP_epithelial", "Neuron")]
#plot
plotDEbyCluster(de_results, ambient=TRUE, ambient_threshold = 0.25) 
dir.create(here("outs", project, "DE_celltype", "plots"))
ggsave(here("outs", project, "DE_celltype", "plots", "DE_celltypes.pdf"),
       height = 7, width = 10
)

## For the Main celltypes ###
#load
de_results <- readRDS(here("processed",project, "DE_results_ambient_edgeR_celltype.RDS"))
#sort
de_results <-  de_results[c( "Astrocyte", "OPCs","Oligo",  "Endothelial","Mural_cells")]
#plot
plotDEbyCluster(de_results, ambient=TRUE, ambient_threshold = 0.25) 
dir.create(here("outs", project, "DE_celltype", "plots"))
ggsave(here("outs", project, "DE_celltype", "plots", "DE_main_celltypes.pdf"),
       height = 7, width = 10
)