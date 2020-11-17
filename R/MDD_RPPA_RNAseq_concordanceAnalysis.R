# The purpose of this script is to load DE proteins from RPPA_limma.R,
# filter out certain antibodies (e.g. phosphoantibodies, those that map to
# multiple HGNC symbols, etc.), join with RNAseq results, and compare 
# DE genes/proteins.

library(tidyverse)
library(eulerr)
library(ComplexHeatmap)

RPPA_limmaFile  <- "MDD_RPPA_limma.R"
colScript  <- "MDD_importColors_pretty.R"
RNA_resultsFile <- "../RNAseq/Data/MDD_RNAseq_shrunkenResults_df.Rdata"

RNAseq_logFC_threshold <- 1.5
RNAseq_pval_threshold  <- 0.01

outDir <- "../plots/RPPA_RNAseq_concordanceHeatmaps/lfc15_fixedNames"

###############################################################################
if(!grepl("R$", getwd())) { setwd("R") }

source(RPPA_limmaFile)
source(colScript)
###############################################################################

# # Getting color codes for plots
# col <- read.csv(color_annoFile, stringsAsFactors = FALSE)
# col <- list(
#   ligand = dplyr::slice(col, 1:8),
#   experimentalTimePoint = dplyr::slice(col, 10:15),
#   secondLigand = dplyr::slice(col, 17:18),
#   collection = dplyr::slice(col, 26:28)
# )
# 
# col <-
#   lapply(col, function(x) {
#     x <- setNames(x[, 2], x[, 1])
#   })

###############################################################################
# Filtering antibody annotation table to non-phospho, non-multi-protein abs
aTrppaTempFull <- 
  read.csv(antibodyFile, header = TRUE, stringsAsFactors = FALSE) %>% 
  mutate(antibody = str_remove(MDD, "-[:alnum:]-[:alnum:]$")) %>% 
  dplyr::rename(hgnc_symbol = Symbols) %>% 
  dplyr::select(antibody, hgnc_symbol, Sites)%>% 
  # filter(is.na(Sites)) %>% 
  mutate(symbolFilt = !str_detect(hgnc_symbol, " "))
  # filter(symbolFilt) %>% 
  # dplyr::select(-symbolFilt)

aTrppaFilt <- 
  read.csv(antibodyFile, header = TRUE, stringsAsFactors = FALSE) %>% 
  mutate(antibody = str_remove(MDD, "-[:alnum:]-[:alnum:]$")) %>% 
  dplyr::rename(hgnc_symbol = Symbols) %>% 
  dplyr::select(antibody, hgnc_symbol, Sites)%>% 
  filter(is.na(Sites)) %>% 
  mutate(symbolFilt = !str_detect(hgnc_symbol, " ")) %>% 
  filter(symbolFilt) %>% 
  dplyr::select(-symbolFilt)

res <-
  res %>% 
  data.frame %>% 
  rownames_to_column("antibody") %>% 
  right_join(aTrppaFilt) %>%
  dplyr::select(hgnc_symbol, contains("experimentalCondition"))

# Loading RNASeq data
load(RNA_resultsFile)

setdiff(res$hgnc_symbol, res_RNA_df$PBS_24$hgnc_symbol)

res$hgnc_symbol[res$hgnc_symbol == "AIM1"]  <- "CRYBG1"
res$hgnc_symbol[res$hgnc_symbol == "CD29"]  <- "ITGB1"
res$hgnc_symbol[res$hgnc_symbol == "CDC2"]  <- "CDK1"
res$hgnc_symbol[res$hgnc_symbol == "PAR"]   <- "F2R"
res$hgnc_symbol[res$hgnc_symbol == "PDGFR"] <- "PDGFRB"
res$hgnc_symbol[res$hgnc_symbol == "PDHK1"] <- "PDK1"
res$hgnc_symbol[res$hgnc_symbol == "RIP"]   <- "RIPK1"
res$hgnc_symbol[res$hgnc_symbol == "C12ORF5"] <- "TIGAR"

setdiff(res$hgnc_symbol, res_RNA_df$PBS_24$hgnc_symbol)

### Filtering RNAseq results to those with corresponding RPPA antibodies
overlap_genes <- base::intersect(res_RNA_df$PBS_24$hgnc_symbol,
                                 res$hgnc_symbol)

res_RNA <-
  Reduce(bind_rows, res_RNA_df) %>% 
  filter(hgnc_symbol %in% overlap_genes)

res_temp <-
  res_RNA %>% 
  mutate(RNA_significant = case_when(padj <= RNAseq_pval_threshold & log2FoldChange >= RNAseq_logFC_threshold ~ "1",
                                     padj <= RNAseq_pval_threshold & log2FoldChange <= -RNAseq_logFC_threshold ~ "-1",
                                     padj > .01 ~ "0",
                                     is.na(padj) ~ "0",
                                     log2FoldChange < RNAseq_logFC_threshold & log2FoldChange > -RNAseq_logFC_threshold ~ "0"))

### Heatmaps: RNAseq/RPPA concordance

RNA_matched <- 
  res_temp %>% 
  dplyr::select(hgnc_symbol, experimentalCondition, RNA_significant) %>% 
  mutate(RNA_significant = as.numeric(RNA_significant)) %>% 
  pivot_wider(names_from = experimentalCondition, values_from = RNA_significant) %>% 
  data.frame %>% 
  column_to_rownames("hgnc_symbol")

colnames(res) <- str_remove(colnames(res), "experimentalCondition")
res <- data.frame(res) %>% 
  column_to_rownames("hgnc_symbol")
  
hmapForRowOrder <- Heatmap(RNA_matched, cluster_columns = FALSE, 
                           show_row_names = FALSE, column_title = "RNAseq",
                           name = "RNAseq", show_column_names = FALSE) 
RNA_matched <- RNA_matched[row_order(hmapForRowOrder), ]
RPPA_matched <- res[rownames(RNA_matched), colnames(RNA_matched)]


# if(!dir.exists(dir)) {
  # dir.create(dir)
# }

combined.meta.forHA <-
  RPPA.meta %>% 
  dplyr::select(experimentalCondition, ligand, experimentalTimePoint) %>% 
  dplyr::rename(Ligand = ligand,
                Time = experimentalTimePoint) %>% 
  filter(Time %in% c(24, 48)) %>% 
  distinct(experimentalCondition, .keep_all = TRUE) %>% 
  mutate(experimentalCondition = as.character(experimentalCondition)) %>% 
  mutate(Ligand = fct_recode(Ligand, 
                             "CTRL" = "ctrl",
                             "BMP2+EGF" = "BMP2",
                             "IFNG+EGF" = "IFNG",
                             "TGFB+EGF" = "TGFB")) %>% 
         # experimentalCondition = fct_recode(experimentalCondition,
                                            # "CTRL_0" = "ctrl_0")) %>% 
  mutate(Ligand = fct_relevel(Ligand, "CTRL", "PBS",
                              "HGF", "OSM", "EGF",
                              "BMP2+EGF", "IFNG+EGF", "TGFB+EGF")) %>% 
  arrange(Ligand, Time) %>% 
  mutate(experimentalCondition = fct_inorder(experimentalCondition))

columnOrder <- 
  combined.meta.forHA %>% 
  pull(experimentalCondition) %>% 
  as.character

names(col_MDD)[c(1,2)] <- c("Ligand", "Time")
ha <- HeatmapAnnotation(df = dplyr::select(combined.meta.forHA, -experimentalCondition),
                        col = col_MDD)

ha1 <- HeatmapAnnotation(df = dplyr::select(combined.meta.forHA, -experimentalCondition),
                         col = col_MDD, show_legend = FALSE)
names(ha1) <- c("", " ")

ind <- match(c("JAK2", "CD274", "EGFR", "MYC",
          "CCNB1", "DUSP4", "RB1", "RPS6", "IRS1"),
          rownames(RNA_matched))

ra <- rowAnnotation(hgnc_symbol = anno_mark(at = ind, labels = rownames(RNA_matched)[ind]))

RNA_matched_factors <-
  RNA_matched %>%
  rownames_to_column("gene") %>% 
  pivot_longer(-gene) %>% 
  mutate(value = as.factor(value)) %>% 
  mutate(value = fct_recode(value,
                            "Upregulated" = "1",
                            "Downregulated" = "-1",
                            "No Change" = "0")) %>% 
  pivot_wider() %>% 
  column_to_rownames("gene")
RNA_matched_factors <- RNA_matched_factors[, columnOrder]

RPPA_matched_factors <-
  RPPA_matched %>%
  rownames_to_column("gene") %>% 
  pivot_longer(-gene) %>% 
  mutate(value = as.factor(value)) %>% 
  mutate(value = fct_recode(value,
                            "Upregulated" = "1",
                            "Downregulated" = "-1",
                            "No Change" = "0")) %>% 
  pivot_wider() %>% 
  column_to_rownames("gene")
RPPA_matched_factors <- RPPA_matched_factors[, columnOrder]

if (!dir.exists(outDir)) {
  dir.create(outDir)
}

pdf(sprintf("%s/MDD_RPPA_RNAseq_significanceHeatmap.pdf", outDir), height = 9, width = 15)
Heatmap(RNA_matched_factors, cluster_columns = FALSE, 
        cluster_rows = FALSE,
        show_row_names = FALSE, column_title = "RNAseq", 
        name = "RNAseq", bottom_annotation = ha1,
        col = c("Upregulated" = "Red",
                "No Change" = "gray95",
                "Downregulated" = "Blue"),
        heatmap_legend_param = c(labels = fct_inorder(as.factor(unique(levels(RNA_matched_factors$PBS_48)))),
                                 title = "RNAseq",
                                 legend_gp = gpar(fill = c("blue", "white", "red"))),
        show_column_names = FALSE) + 
Heatmap(RPPA_matched_factors, cluster_columns = FALSE, 
        show_row_names = FALSE, column_title = "RPPA", 
        name = "RPPA",  bottom_annotation = ha, right_annotation = ra,
        col = c("Upregulated" = "Red",
                "No Change" = "gray95",
                "Downregulated" = "Blue"),
        heatmap_legend_param = c(labels = unique(levels(RNA_matched_factors$PBS_48)),
                                 title = "RPPA",
                                 legend_gp = gpar(fill = c("blue", "white", "red"))),
        show_column_names = FALSE)
dev.off()


### Reshaping data

# The differential expression analyses are reshaped into a list with `assay_significance` names
# (e.g. "RPPA_Up" or "RNA_noChange"), so that they may be compared in a set analysis.

RNA_toJoin <-
  RNA_matched %>% 
  rownames_to_column("hgnc_symbol") %>% 
  pivot_longer(-hgnc_symbol, 
               names_to = "experimentalCondition", 
               values_to = "RNAseq")

RPPA_toJoin <-
  RPPA_matched %>% 
  rownames_to_column("hgnc_symbol") %>% 
  pivot_longer(-hgnc_symbol, 
               names_to = "experimentalCondition", 
               values_to = "RPPA")

joined <- 
  RNA_toJoin %>% 
  left_join(RPPA_toJoin) %>% 
  pivot_longer(-c(hgnc_symbol, experimentalCondition), 
               names_to = "assay", values_to = "significant") %>% 
  mutate(symbol_condition = paste(hgnc_symbol, experimentalCondition, sep = "_"),
         significant = case_when(significant == 1  ~ "Up",
                                 significant == -1 ~ "Down",
                                 significant == 0  ~ "NoChange"))

DEFeatureList <-
  joined %>% 
  group_by(assay, significant) %>% 
  summarize(sc = list(symbol_condition))
DEFeatureList <- setNames(DEFeatureList$sc, 
                          nm = paste(DEFeatureList$assay, 
                                     DEFeatureList$significant, 
                                     sep = "_"))
intersect(DEFeatureList$RNAseq_Down, DEFeatureList$RPPA_Up)
# This euler diagram shows overlaps between the "up" and "down" sets of genes. 
# The "noChange" groups are not shown.

pdf(sprintf("%s/MDD_RPPA_RNAseq_significanceEuler.pdf", outDir), height = 6, width = 6)
euler(DEFeatureList[!grepl("NoChange", names(DEFeatureList))]) %>% 
  plot(quantities = TRUE, 
       fills = c("steelblue4", "firebrick4", "steelblue1", "firebrick2"))
dev.off()

intersect(DEFeatureList$RNAseq_NoChange,
          DEFeatureList$RPPA_NoChange) %>% 
  length

c(DEFeatureList$RNAseq_NoChange,
          DEFeatureList$RPPA_NoChange) %>% 
  unique %>% 
  length

c(DEFeatureList$RNAseq_Up,
  DEFeatureList$RPPA_Up) %>% 
  unique %>% 
  length
  
c(DEFeatureList$RNAseq_Down,
  DEFeatureList$RPPA_Down) %>% 
  unique %>% 
  length



concordant_up <-
  base::intersect(DEFeatureList$RNAseq_Up,
                  DEFeatureList$RPPA_Up) %>% 
  data.frame(hgnc_symbol = str_split_fixed(., "_", 2)[, 1],
             experimentalCondition = str_split_fixed(., "_", 2)[, 2],
             stringsAsFactors = FALSE) %>% 
  mutate(direction = "Up")

concordant_down <-
  base::intersect(DEFeatureList$RNAseq_Down,
                  DEFeatureList$RPPA_Down) %>% 
  data.frame(hgnc_symbol = str_split_fixed(., "_", 2)[, 1],
             experimentalCondition = str_split_fixed(., "_", 2)[, 2],
             stringsAsFactors = FALSE) %>% 
  mutate(direction = "Down")

concordant_up
concordant_down

# res_RNA_working <-
#   res_RNA %>% 
#   mutate(toMatch = paste(hgnc_symbol, experimentalCondition ,sep = "_"))

concordant_TW <-
  rbind(concordant_up, concordant_down)

concordant_TW %>% 
  arrange(desc(direction), experimentalCondition, hgnc_symbol) %>% 
  dplyr::select(-.) %>% 
  write_csv(sprintf("%s/MDD_RNAseq_RPPA_concordanceTable.csv", outDir))

