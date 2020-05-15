library(biomaRt)
library(DESeq2)
library(tidyverse)
library(ComplexHeatmap)
library(circlize)
library(knitr)
library(kableExtra)
library(cowplot)
library(cmapR)
theme_set(theme_cowplot())

OSM_Reactome <- "/Users/derrickd/workspace/MDD/misc/OSM_REACTOME_GENESET.csv"

source("MDD_RNAseq_heatmapFunctions.R")
source("MDD_import_RNAseq_ensembl.R")
source("MDD_importColors_pretty.R")

###############################################################################
# Changing ensembl gene IDs to hgnc_symbols
RNAseqL3Z <- 
  RNAseqL3Z %>% 
  rownames_to_column("ensembl_gene_id") %>% 
  left_join(at.RNA) %>% 
  distinct(hgnc_symbol, .keep_all = TRUE) %>% 
  filter(hgnc_symbol != "") %>% 
  dplyr::select(-ensembl_gene_id, -gene_biotype) %>% 
  column_to_rownames("hgnc_symbol")

sa.RNA.L3$experimentalCondition[sa.RNA.L3$experimentalCondition == "ctrl_0"] <- "CTRL_0"
sa.RNA.L3$Ligand[sa.RNA.L3$experimentalCondition == "ctrl_0"] <- "CTRL"

MDD_RNAseq_Zscore_Medians <-
  RNAseqL3Z %>% 
  rownames_to_column("hgnc_symbol") %>% 
  pivot_longer(-hgnc_symbol, values_to = "z", names_to = "specimenID") %>% 
  left_join(sa.RNA.L3) %>% 
  group_by(hgnc_symbol, experimentalCondition) %>% 
  summarize(medianZ = median(z, na.rm = TRUE)) %>% 
  mutate(medianZ = replace_na(medianZ, 0)) %>% 
  ungroup %>% 
  pivot_wider(names_from = experimentalCondition, values_from = medianZ) %>% 
  column_to_rownames("hgnc_symbol")

# getting Peroxisome, UV-Response, and EMT Hallmark Gene Sets
which <- c("HALLMARK_PEROXISOME", "HALLMARK_EPITHELIAL_MESENCHYMAL_TRANSITION",
           "HALLMARK_UV_RESPONSE_UP", "HALLMARK_UV_RESPONSE_DN", "HALLMARK_MYC_TARGETS_V1")

# Reading in gene sets
hallmarks <- parse.gmt("../misc/h.all.v6.2.symbols.gmt")
hallmarks <- hallmarks[which]

sa.RNA.L4 <-
  sa.RNA.L3 %>% 
  distinct(experimentalCondition, .keep_all = TRUE) %>% 
  mutate(experimentalCondition = fct_inorder(as.factor(experimentalCondition))) %>% 
  dplyr::select(experimentalCondition, Ligand, Time) %>% 
  mutate(Ligand = fct_inorder(as.factor(Ligand)),
         Time   = fct_inorder(as.factor(Time)))

MDD_RNAseq_Zscore_Medians <- MDD_RNAseq_Zscore_Medians[, as.character(sa.RNA.L4$experimentalCondition)]

ha.RNA.L4 <- HeatmapAnnotation(
  df = sa.RNA.L4 %>% dplyr::select(-experimentalCondition), col = col
)

exprFilt <-
  RNAseqL3 %>% 
  data.frame %>% 
  rownames_to_column("ensembl_gene_id") %>% 
  pivot_longer(-ensembl_gene_id, names_to = "specimenID", values_to = "log2fpkm") %>% 
  left_join(sa.RNA.L3, by = "specimenID") %>% 
  group_by(experimentalCondition, ensembl_gene_id) %>% 
  summarize(expressed = sum(log2fpkm >= 0.5) >= 2) %>% 
  ungroup %>% 
  left_join(at.RNA) %>% 
  group_by(hgnc_symbol) %>% 
  summarize(expressed = any(expressed, na.rm = TRUE)) %>% 
  mutate(expressed = replace_na(expressed, FALSE))

MDD_RNAseq_Zscore_Medians <- MDD_RNAseq_Zscore_Medians[exprFilt %>% 
                                                         filter(expressed) %>% 
                                                         pull(hgnc_symbol), ]

plots <- lapply(hallmarks, function(x) {
  notShown <- setdiff(x$entry, rownames(MDD_RNAseq_Zscore_Medians))
  print(notShown)
  shown <- base::intersect(x$entry, rownames(MDD_RNAseq_Zscore_Medians))
  print(length(shown))
  # pdf(sprintf("%s/MDD_RNAseq_%s_Heatmap.pdf", outDir, x$head), height = 10, width = 8)
  HeatmapL4C(MDD_RNAseq_Zscore_Medians[shown, ], 
             showRow = TRUE)
  # dev.off()

})


outDir <- "../plots/RNAseq_heatmaps/Selected_Hallmark_Gene_Sets_Zscore"

if (!dir.exists(outDir)) {
  dir.create(outDir, recursive = TRUE)
}

pdf(sprintf("%s/MDD_RNAseq_%s_medianZHeatmap.pdf", outDir, hallmarks$HALLMARK_PEROXISOME$head), height = 12, width = 8)
plots$HALLMARK_PEROXISOME
dev.off()

pdf(sprintf("%s/MDD_RNAseq_%s_medianZHeatmap.pdf", outDir, hallmarks$HALLMARK_EPITHELIAL_MESENCHYMAL_TRANSITION$head), height = 20, width = 8)
plots$HALLMARK_EPITHELIAL_MESENCHYMAL_TRANSITION
dev.off()

pdf(sprintf("%s/MDD_RNAseq_%s_medianZHeatmap.pdf", outDir, hallmarks$HALLMARK_UV_RESPONSE_UP$head), height = 15, width = 8)
plots$HALLMARK_UV_RESPONSE_UP
dev.off()

pdf(sprintf("%s/MDD_RNAseq_%s_medianZHeatmap.pdf", outDir, hallmarks$HALLMARK_UV_RESPONSE_DN$head), height = 13, width = 8)
plots$HALLMARK_UV_RESPONSE_DN
dev.off()

pdf(sprintf("%s/MDD_RNAseq_%s_medianZHeatmap.pdf", outDir, hallmarks$HALLMARK_MYC_TARGETS_V1$head), height = 22, width = 8)
plots$HALLMARK_MYC_TARGETS_V1
dev.off()

###############################################################################

outDir <- "../plots/RNAseq_heatmaps/Reactome_genesets_Zscore/"

if (!dir.exists(outDir)) {
  dir.create(outDir, recursive = TRUE)
}

osm_reactome <- read.csv(OSM_Reactome, stringsAsFactors = FALSE, header = FALSE)
colnames(osm_reactome) <- c("Reactome", "hgnc_symbol", "uniprot_id")

osm_reactome_filt <-
  osm_reactome %>% 
  filter(hgnc_symbol %in% rownames(MDD_RNAseq_Zscore_Medians))

plotsReactome <- lapply(unique(osm_reactome_filt$Reactome), function(x) {
  
  whichGenes <- 
    osm_reactome_filt %>% 
    filter(Reactome == x) %>% 
    pull(hgnc_symbol)
  
  HeatmapL4C(MDD_RNAseq_Zscore_Medians[whichGenes, ], 
             column_title = x,
             showRow = TRUE)
  
})

plotsReactome[[6]] <-
  HeatmapL4C(MDD_RNAseq_Zscore_Medians[osm_reactome_filt %>% 
                                         pull(hgnc_symbol), ],
             split = osm_reactome_filt %>% pull(Reactome),
             row_title_gp = gpar(fontsize = 7.5),
           column_title = "All",
           showRow = TRUE)

plotsReactome[[7]] <-
  HeatmapL4C(MDD_RNAseq_Zscore_Medians[osm_reactome_filt %>% 
                                         pull(hgnc_symbol) %>% 
                                         unique, ], 
             column_title = "All unique",
             showRow = TRUE)


pdf(sprintf("%s/MDD_RNAseq_%s_Heatmap.pdf", outDir, "Reactome"), height = 6, width = 7)
plotsReactome[[1]]
plotsReactome[[2]]
plotsReactome[[3]]
plotsReactome[[4]]
plotsReactome[[5]]
dev.off()

pdf(sprintf("%s/MDD_RNAseq_%s_Heatmap.pdf", outDir, "ReactomeCombined"), height = 10, width = 7)
plotsReactome[[7]]
plotsReactome[[6]]
dev.off()



