#SCRIPT FOR PLOT WITH FC ON Y AXIS, GROUPS IN X AXIS

## set up ---

library(here) #reproducible paths
library(scater) #feature plots
library(patchwork) # agregate plots
library(pals) # for palettes with large n #kelly()22, #polychrome()#36, cols25()

# project
age <- "old"

# colours
# remove the black and white from the pallete, still 20 colours left
kelly_col <- unname(kelly()[-c(1,2)])

## import and transform data ---
# this is the output from edgeR
# a list with the DE results, each element named as one of the clusters and 
# contain a DFrame with the DE for that cluster
de_results <- readRDS(here("processed", age, "DE_k20_de_results.RDS"))

# transform all the DFrame to a standart df format
de_results_dfs <- lapply(de_results, 
                         function(x) {
                           na.omit(as.data.frame(x))
                         }
)
# merge all the dfs from the list in a single big dataframe
de_results_df <- bind_rows(de_results_dfs, .id = "cluster") %>% 
  # only keep relevant columns
  select(cluster, logFC, FDR) %>% 
  # add column that indicates if the result is significant
  mutate(Significant = ifelse(FDR < 0.05, 
                              yes = ifelse(logFC > 0, "upregulated", "downregulated"),
                              no = "not-significant")) %>% 
  # order the levels to display in correct order in plot
  mutate(Significant = forcats::fct_relevel(Significant, c("upregulated", "downregulated", "not-significant"))) %>% 
  mutate(cluster = forcats::fct_relevel(cluster, c("Astrocyte_1", "Astrocyte_2", "Astrocyte_3", "OligoAstro", "Oligo_1", "Oligo_2", "OPCs", "mNeurons", "iNeurons & NRPs", "Lymphocytes", "BAMs", "DCs",  "Endothelial", "fEndothelia", "Mural_cells", "ChP_epithelial"))) %>% 
  # sort so the significant values are plotted on top of the non significants
  arrange(desc(Significant))

## plot ---
ggplot(data = de_results_df, mapping = aes(x=cluster, y=logFC, color=Significant)) + 
  geom_jitter(size = 0.5) + 
  scale_colour_manual(values = c("darkred", "darkgreen", "#999999" )) +
  theme_minimal() + 
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1),
        legend.title = element_blank()) +
  ylab("logFC fire-mice to wt") +
  xlab(element_blank())
