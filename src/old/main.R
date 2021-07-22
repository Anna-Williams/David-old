# main script that calls other scripts in correct order
library(here)

age <- "old"

# Create sce object
source(here("src", age, "import.R"))
# Gene and Cell QC,
rmarkdown::render(here("src", age, "QC_01.Rmd"))
# normalisation
rmarkdown::render(here("src", age, "normalise_01.Rmd"))
#feature selection and dimensional reduction
rmarkdown::render(here("src", age, "feature_selection_dimred_01.Rmd"))
# clustering
rmarkdown::render(here("src", age, "clustering_01.Rmd"))
# annotation
rmarkdown::render(here("src", age, "annotation_01.Rmd"))
# cluster QC
rmarkdown::render(here("src", age, "cluster_QC_k5.Rmd"))
# stricter thresholds in cell and gene QC
source(here("src", age, "QC_normalise_02.R"))
# redo a feature selection and dimensional reduction
rmarkdown::render(here("src", age, "feature_selection_dimred_02.Rmd"))
# redo clustering 
rmarkdown::render(here("src", age, "clustering_02.Rmd"))
# redo annotation
rmarkdown::render(here("src", age, "annotation_02.Rmd"))
# shiny app
source(here("src", age, generate_shiny.R))