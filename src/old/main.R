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
rmarkdown::render(here("src", age, "annotation.Rmd"))
# cluster QC
rmarkdown::render(here("src", age, "cluster_QC_k5.Rmd"))