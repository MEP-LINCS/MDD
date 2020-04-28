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
  ligand = dplyr::slice(col, 1:8),
  experimentalTimePoint = dplyr::slice(col, 10:15),
  # secondLigand = dplyr::slice(col, 17:18),
  # replicate = dplyr::slice(col, 19:24),
  collection = dplyr::slice(col, 26:28)
)

col <-
  lapply(col, function(x) {
    x <- setNames(as.character(x[, 2]), x[, 1])
  })

order <- c("ctrl", "PBS", "BMP2", "IFNG", "TGFB", "HGF", "OSM", "EGF")
col$ligand <- col$ligand[order]
