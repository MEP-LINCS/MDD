# The purpose of this script is to load DE proteins from RPPA_limma.R,
# filter out certain antibodies (e.g. phosphoantibodies, those that map to
# multiple HGNC symbols, etc.), join with RNAseq results, and compare 
# DE genes/proteins.

library(tidyverse)
library(eulerr)

RPPA_limmaFile  <- "MDD_RPPA_limma.R"
color_annoFile  <- "../misc/MDD_color_codes.csv"
RNA_resultsFile <- "../RNAseq/Data/MDD_RNAseq_shrunkenResults_df.Rdata"

RNAseq_logFC_threshold <- 1.5
RNAseq_pval_threshold  <- 0.01

###############################################################################
source(RPPA_limmaFile)
###############################################################################

# Getting color codes for plots
col <- read.csv(color_annoFile, stringsAsFactors = FALSE)
col <- list(
  ligand = dplyr::slice(col, 1:8),
  experimentalTimePoint = dplyr::slice(col, 10:15),
  secondLigand = dplyr::slice(col, 17:18),
  collection = dplyr::slice(col, 26:28)
)

col <-
  lapply(col, function(x) {
    x <- setNames(x[, 2], x[, 1])
  })

###############################################################################
# Filtering antibody annotation table to non-phospho, non-multi-protein abs
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
  
RPPA_matched <- res[rownames(RNA_matched), colnames(RNA_matched)]

Heatmap(RPPA_matched, cluster_columns = FALSE, show_row_names = FALSE, column_title = "RPPA") + 
  Heatmap(RNA_matched, cluster_columns = FALSE, show_row_names = FALSE, column_title = "RNAseq")

dir <- "../plots/RPPA_RNAseq_concordanceHeatmaps/lfc15"

if(!dir.exists(dir)) {
  dir.create(dir)
}

combined.meta.forHA <-
  RPPA.meta %>% 
  dplyr::select(experimentalCondition, ligand, experimentalTimePoint) %>% 
  dplyr::rename(Ligand = ligand,
                Time = experimentalTimePoint) %>% 
  filter(Time %in% c(24, 48)) %>% 
  distinct(experimentalCondition, .keep_all = TRUE)

names(col)[c(1,2)] <- c("Ligand", "Time")

ha <- HeatmapAnnotation(df = dplyr::select(combined.meta.forHA, -experimentalCondition),
                        col = col)

ind <- match(c("JAK2", "CD274", "EGFR", "MYC",
          "CCNB1", "DUSP4", "RB1", "RPS6", "IRS1"),
          rownames(RNA_matched))

ra <- rowAnnotation(hgnc_symbol = anno_mark(at = ind, labels = rownames(RNA_matched)[ind]))

pdf(sprintf("%s/MDD_RPPA_RNAseq_significanceHeatmap.pdf", dir), height = 9, width = 14)
Heatmap(RNA_matched, cluster_columns = FALSE, 
        show_row_names = FALSE, column_title = "RNAseq", 
        name = "RNAseq", bottom_annotation = ha,
        show_column_names = FALSE) + 
Heatmap(RPPA_matched, cluster_columns = FALSE, 
          show_row_names = FALSE, column_title = "RPPA", 
          name = "RPPA",  bottom_annotation = ha, right_annotation = ra,
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

# This euler diagram shows overlaps between the "up" and "down" sets of genes. 
# The "noChange" groups are not shown.



pdf(sprintf("%s/MDD_RPPA_RNAseq_significanceEuler.pdf", dir), height = 6, width = 6)
euler(DEFeatureList[!grepl("NoChange", names(DEFeatureList))]) %>% 
  plot(quantities = TRUE, 
       fills = c("steelblue4", "firebrick4", "steelblue1", "firebrick2"))
dev.off()
 
