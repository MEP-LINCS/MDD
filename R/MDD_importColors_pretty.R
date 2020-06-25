# This script imports RNAseq data into the environment.
# Sample annotations are also imported.
# Complex Heatmap annotations are also created from the sample annotations.
library(tidyverse)
colFile            <- "../misc/MDD_color_codes.csv"

###############################################################################

if(!grepl("R$", getwd())) {
  setwd("R")
}

col <- read.csv(colFile)
col <- list(
  Ligand = dplyr::slice(col, 1:8),
  Time = dplyr::slice(col, 10:15),
  # secondLigand = dplyr::slice(col, 17:18),
  # replicate = dplyr::slice(col, 19:24),
  collection = dplyr::slice(col, 26:28)
)

col <-
  lapply(col, function(x) {
    x <- setNames(as.character(x[, 2]), x[, 1])
  })
names(col$Ligand)[1] <- "CTRL"
names(col$Ligand)[6:8] <- sprintf("%s+EGF", names(col$Ligand)[6:8])

order <- c("CTRL", "PBS", "HGF", "OSM", "EGF", "BMP2+EGF", "IFNG+EGF", "TGFB+EGF")
col$Ligand <- col$Ligand[order]
