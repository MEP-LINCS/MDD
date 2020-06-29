library(RColorBrewer)
library(tidyverse)
library(plotly)
library(apeglm)
library(DESeq2)
library(circlize)
library(ComplexHeatmap)
library(cmapR)
library(cowplot)
library(fgsea)
library(reshape2)

load("../RNAseq/MDD_RNAseq_shrunkenResults.Rdata")
outDir <- "../RNAseq/misc/"
hallmark <- parse.gmt("../misc/h.all.v6.2.symbols.gmt") %>% 
  lapply(., function(x) x[[3]])

# theme_set(theme_cowplot())

makePalette <- function(data, squish_range = c(0.005, .995)) {
  # The purpose of this function is to make a ColorRamp2 palette by
  # "squishing" extreme values, then using the squished values and the palette
  # from Mark to produce a palette for ComplexHeatmap
  pal_choice <- "RdBu"
  
  squish <- function(mat, lower_thresh = .005, upper_thresh = .995) {
    # The purpose of this function is to "squish" a matrix.
    # Values in the matrix that are more extreme than input quantiles
    # are "squished" down to the quantile value.
    range <- quantile(unlist(mat), c(lower_thresh, upper_thresh))
    
    if (range[1] > -1) {
      range[1] <- -1
    }
    
    if (range[2] < 1) {
      range[2] <- 1
    }
    
    mat <- apply(mat, c(1,2), function(x) {
      if (x > range[2]) {
        x <- range[2]
      } else if (x < range[1]) {
        x <- range[1]
      }
      return(x)
    })
    return(mat)
  }
  
  pal <- brewer.pal(n=11,pal_choice)
  rc1 <- colorRampPalette(colors = c(pal[1], pal[2]), space="Lab")(11)
  for(i in 2:10){
    tmp <- colorRampPalette(colors = c(pal[i], pal[i+1]), space = "Lab")(10)
    rc1 <- c(rc1,tmp)
  }
  
  if (!is_null(squish_range)) {
    data <- squish(data, 
                   lower_thresh = squish_range[1], 
                   upper_thresh = squish_range[2])
  }
  
  rb1 <- seq(min(data), 0, length.out=50+1)
  rb2 <- seq(0, max(data), length.out=50+1)[-1]
  rampbreaks <- c(rb1, rb2)
  palette <- colorRamp2(breaks = rampbreaks,
                        colors = rev(rc1))
}

# Getting shrunken Log2FC estimates
# lfc <- lapply(resultsNames(dds)[-1], function(x) {
#   x <- lfcShrink(dds, coef = x, type = "apeglm")
# })
# 
# names(lfc) <- str_extract(resultsNames(dds)[-1], "[:alnum:]+_[248]{2}")
# Creating matrix of shrunken Log2FC estimates

##############################################################################
# fGSEA analysis with list of log2FC
L2FClist <- lapply(res_RNA_df, function(x) {
  # x <- data.frame(x)
  x <- x %>%
    # rownames_to_column("ensembl_gene_id") %>% 
    # left_join(annoTable) %>% 
    filter(!is.na(padj)) %>% 
    distinct(hgnc_symbol, .keep_all = TRUE)
  x <- setNames(x$log2FoldChange, x$hgnc_symbol)
  x
})

fGSEA <- lapply(L2FClist, function(x) {
  set.seed(5)
  x <- fgsea(hallmark, x, 10000)
})

fGSEASummaryNES <- sapply(fGSEA, function(x) {
  x <- x %>% arrange(pathway)
  x <- setNames(x[, 5], x[, 1])
  x
  }) %>% 
  t() %>% 
  data.frame() %>% 
  rownames_to_column("experimentalCondition") %>% 
  melt()

fGSEASummaryFDR <- sapply(fGSEA, function(x) {
  x <- x %>% arrange(pathway)
  x <- setNames(x[, 3] < 0.05, x[, 1])
  x
  }) %>% 
  t() %>% 
  data.frame() %>% 
  rownames_to_column("experimentalCondition") %>% 
  melt(id.vars = "experimentalCondition") %>% 
  dplyr::rename(significant = value)

remove <- 
  fGSEASummaryFDR %>% 
  group_by(variable) %>% 
  summarize(remove = sum(significant) < 1) %>% 
  filter(remove) %>% 
  pull(variable) %>% 
  as.character()

fGSEASummary <- 
  fGSEASummaryNES %>% 
  left_join(fGSEASummaryFDR) %>% 
  filter(!(variable %in% remove))  %>% 
  mutate(variable = str_remove(variable, "HALLMARK_"))


write_csv(dcast(fGSEASummary, variable ~ experimentalCondition),
          path = sprintf("%s/MDD_RNAseq_HallmarkNES.csv", outDir))
