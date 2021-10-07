library(here) #reproducible paths
library(scater) #feature plots
library(patchwork) # agregate plots
library(pals) # for palettes with large n #kelly()22, #polychrome()#36, cols25()

age <- "old"

# remove the black and white from the pallete, still 20 colours left
kelly_col <- unname(kelly()[-c(1,2)])

if (!file.exists(here("processed", age, "sce_anno_02.RDS"))) {
  sce <- readRDS(here("processed", age, "sce_clusters_02.RDS"))
}else{
  sce <- readRDS(here("processed", age, "sce_anno_02.RDS"))
}


plotTSNE(sce, colour_by = "clusters_named", text_by = "clusters_named", text_size = 4) + scale_color_manual(values = kelly_col)


## plot 

# import data
de_results <- readRDS(here("processed", age, "DE_k20_de_results.RDS"))
# transform to df format
de_results_dfs <- lapply(de_results, 
                         function(x) {
                           na.omit(as.data.frame(x))
                         }
)
# merge all in a single dataframe
de_results_df <- bind_rows(de_results_dfs, .id = "cluster") %>% 
  select(cluster, logFC, FDR) %>% 
  mutate(Significant = ifelse(FDR < 0.05, 
                              yes = ifelse(logFC > 0, "upregulated", "downregulated"),
                              no = "not-significant")) %>% 
  mutate(Significant = forcats::fct_relevel(Significant, c("upregulated", "downregulated", "not-significant"))) %>% 
  mutate(cluster = forcats::fct_relevel(cluster, c("Astrocyte_1", "Astrocyte_2", "Astrocyte_3", "OligoAstro", "Oligo_1", "Oligo_2", "OPCs", "mNeurons", "iNeurons & NRPs", "Lymphocytes", "BAMs", "DCs",  "Endothelial", "fEndothelia", "Mural_cells", "ChP_epithelial"))) %>% 
  arrange(desc(Significant))

ggplot(data = de_results_df, mapping = aes(x=cluster, y=logFC, color=Significant)) + 
  geom_jitter(size = 0.5) + 
  scale_colour_manual(values = c("darkred", "darkgreen", "#999999" )) +
  theme_minimal() + 
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1),
        legend.title = element_blank()) +
  ylab("logFC fire-mice to wt") +
  xlab(element_blank())
