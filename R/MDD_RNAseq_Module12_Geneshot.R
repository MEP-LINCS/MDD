library(DESeq2)
library(biomaRt)
library(tidyverse)
library(plotly)
library(cowplot)

if(!grepl("R$", getwd())) {setwd("R")}

# gSfile <- "../misc/MDD_EMT_geneshot_autorif.tsv"
# gSfile <- "../misc/MDD_cellcycle_geneshot_autorif.tsv"
gSfile <- "../misc/MDD_immune_interferon_geneshot_autorif.tsv"
moduleDir <- "../misc/Ligand_tables"
resFile <- "../RNAseq/MCF10A_RNAseq_IFNG48_Results.Rdata"

outDir <- "../plots/Module12_Geneshot"
outPathWidget <- "/Users/derrickd/workspace/MDD/plots/Module12_Geneshot/MDD_RNAseq_geneshot_scatterplot.html"

###############################################################################
# Importing geneShot results for gSfile
geneShot <- read.table(gSfile,
                       sep =  "\t",
                       header = TRUE) %>% 
  dplyr::rename(hgnc_symbol = Gene)

# importing module features
myFiles        <- list.files(moduleDir, full.names = TRUE)
names(myFiles) <- str_replace(str_extract(myFiles, "cluster_[0-9]+"), 
                              "cluster", "module")
fT <- lapply(myFiles, read.csv,
             header = FALSE, stringsAsFactors = FALSE)

fT <-
  lapply(fT, function(x) {
    colnames(x) <- c("feature", sprintf("V%s", seq(2, 10, 1)), "type", "symbol")
    x
  }
  )

moduleRNAseqFeatures <- lapply(fT, function(x) {
  x <- x %>%
    filter(type == "RNA") %>%
    distinct(symbol) %>%
    pull(symbol)
})
moduleRNAseqFeatures <- moduleRNAseqFeatures[c(1, 5:12, 2:4)]

module12DF <-
  data.frame(hgnc_symbol = moduleRNAseqFeatures$module_12,
             stringsAsFactors = FALSE) %>% 
  left_join(geneShot)

# Importing RNAseq results
load(resFile)

###############################################################################

a <-
  IFNG48_RNA %>% 
  filter(!is.na(padj)) %>% 
  right_join(module12DF) %>% 
  filter(!is.na(padj)) %>% 
  mutate(padj = case_when(padj == 0 ~ 1.726570e-301,
                          padj != 0 ~ padj)) %>% 
  ggplot(aes(x = log2FoldChange, 
             y = -log10(padj),
             col = Rank,
             text = hgnc_symbol)) +
  geom_point(size = 1.5, alpha = .75) +
  xlab("IFNG_48 log2FC") +
  ylab("IFNG_48 -log10(padj)")+
  ggtitle("Module 12 Genes",
          "geneshot: 'immune', 'interferon'") +
  geom_hline(yintercept = 2) +
  geom_vline(xintercept = 1.5) +
  scale_color_gradient(low = "pink", high = "red")

b <-
  IFNG48_RNA %>% 
  filter(!is.na(padj)) %>% 
  right_join(module12DF) %>% 
  filter(!is.na(padj)) %>% 
  mutate(padj = case_when(padj == 0 ~ 1.726570e-301,
                          padj != 0 ~ padj)) %>% 
  ggplot(aes(x = log2FoldChange, 
             y = -log10(padj),
             col = Rank,
             text = hgnc_symbol)) +
  geom_point(size = 1.5, alpha = .75) +
  xlab("IFNG_48 log2FC") +
  ylab("IFNG_48 -log10(padj)") +
  ggtitle("Module 12 Genes",
          "geneshot: 'immune', 'interferon'") +
  xlim(c(1.5, 19)) +
  ylim(c(2, 305)) +
  scale_color_gradient(low = "pink", high = "red")

c <-
  IFNG48_RNA %>% 
  filter(!is.na(padj)) %>% 
  right_join(module12DF) %>% 
  filter(!is.na(padj)) %>% 
  mutate(padj = case_when(padj == 0 ~ 1.726570e-301,
                          padj != 0 ~ padj)) %>% 
  ggplot(aes(x = log2FoldChange, 
             y = -log10(padj),
             col = Rank,
             label = hgnc_symbol)) +
  geom_text(size = 1.5) +
  xlab("IFNG_48 log2FC") +
  ylab("IFNG_48 -log10(padj)")+
  ggtitle("Module 12 Genes",
          "geneshot: 'immune', 'interferon'") +
  xlim(c(3, 19)) +
  ylim(c(2, 305)) +
  scale_color_gradient(low = "pink", high = "red")


if (!dir.exists(outDir)) {dir.create(outDir)}

pdf(sprintf("%s/MDD_RNAseq_Module12_geneshot_scatterplot.pdf", outDir), 
    height = 5, width = 5)
a
b
dev.off()

b_plotly <- ggplotly(b)

htmlwidgets::saveWidget(b_plotly, file = outPathWidget)

pdf(sprintf("%s/MDD_RNAseq_Module12_geneshot_scatterplot_labels.pdf", outDir), 
    height = 7, width = 7)
c
dev.off()

