### plots ###
# libraries 
library(here) # reproducible paths
library(scater) # plot reduced dims
library(dplyr) #manipulate df
library(pals) # for palettes with large n #kelly()22, #polychrome()#36, cols25()
library("RColorBrewer") # for continues palette

# remove the black and white from the pallete, still 20 colours left
kelly_col <- unname(kelly()[-c(1,2)])

# set up
age <- "old"
sce <- readRDS(here("processed", age, "sce_anno_02.RDS"))

# define colors
myPalette <- colorRampPalette(c("grey85","#FFF7EC","#FEE8C8","#FDD49E","#FDBB84", 
                                         "#FC8D59","#EF6548","#D7301F","#B30000","#7F0000"))


## feature plots with serpina3n and c4b
# scales
sc_serp <- scale_colour_gradientn(colours = myPalette(100), limits=c(0, 5), name = "Serpina3n")
sc_c4b <- scale_colour_gradientn(colours = myPalette(100), limits=c(0, 6), name = "C4b")
plotTSNE(sce, colour_by = "Serpina3n" ) + sc_serp
plotTSNE(sce, colour_by = "C4b" ) + sc_c4b

# divide the two objects
sce_ko <- sce[,sce$genotype == "KO"]
sce_ctr <- sce[,sce$genotype == "WT"]

# plot serpina3n side by side
pdf(here("outs", age, "plots", "Serpina3n.pdf"), width = 14)
gridExtra::grid.arrange(
  plotTSNE(sce_ctr, colour_by = "Serpina3n") + sc_serp +
    ggtitle("Control"), 
  plotTSNE(sce_ko, colour_by = "Serpina3n") + sc_serp +
    ggtitle("Fire mice"),
  ncol = 2
)
dev.off()
# plot c4b side by side
pdf(here("outs", age, "plots", "C4b.pdf"), width = 14)
gridExtra::grid.arrange(
  plotTSNE(sce_ctr, colour_by = "C4b", point_alpha = 1) + sc_c4b +
    ggtitle("Control"), 
  plotTSNE(sce_ko, colour_by = "C4b", point_alpha = 1) + sc_c4b +
    ggtitle("Fire mice"),
  ncol = 2
)
dev.off()
