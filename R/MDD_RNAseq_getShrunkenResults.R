# The purpose of this script is to run differential expression analyses
# for each condition vs ctrl_0.

library(tidyverse)
library(DESeq2)
library(biomaRt)

MDD_DESeqDataSetFile <- "../RNAseq/Data/MDD_RNAseq_DESeqDataSet.Rdata"
resFile <- "../RNAseq/Data/MDD_RNAseq_DESeq2Results.Rdata"
hgnc_referenceFile <- "/Users/derrickd/workspace/MDD/misc/HGNC_non_alt_loci_set.txt"

outDir <- "../RNAseq/Data"

###############################################################################

formatResults <- function(x, x_nm) {
  x <- 
    data.frame(x) %>% 
    rownames_to_column("ensembl_gene_id") %>% 
    left_join(annoTable) %>% 
    mutate(identifier = case_when(!is.na(hgnc_symbol) ~ hgnc_symbol,
                                  is.na(hgnc_symbol)  ~ ensembl_gene_id)) %>% 
    arrange(padj) %>% 
    distinct(identifier, .keep_all = TRUE) %>% 
    mutate(experimentalCondition = x_nm) %>% 
    dplyr::select(experimentalCondition, ensembl_gene_id, hgnc_symbol, everything()) %>% 
    dplyr::select(-identifier)
  
  x
}

###############################################################################

if( str_extract(getwd(), "[:alnum:]+$") != "R" ) {setwd("R")}

load(MDD_DESeqDataSetFile)

#Get annotations to convert from ensemble to HGNC
hgnc_ref <- read_tsv(hgnc_referenceFile) # for symbols missing from biomaRt

mart <- useMart(biomart = "ENSEMBL_MART_ENSEMBL",
                dataset = "hsapiens_gene_ensembl",
                host = "uswest.ensembl.org")
annoTable <- getBM(attributes = c("ensembl_gene_id",
                                  "hgnc_symbol"),
                   mart = mart) %>% 
  left_join(dplyr::select(hgnc_ref, 
                          ensembl_gene_id, symbol)) %>% 
  mutate(hgnc_symbol = case_when(hgnc_symbol != "" ~ hgnc_symbol,
                                 hgnc_symbol == ""  ~ symbol)) %>% 
  dplyr::select(-symbol)


if (!file.exists(resFile)) {
  res <- list()
  
  for (i in 2:15) {
    res[[i - 1]] <- lfcShrink(dds, coef = i, type = "apeglm")
  }
  
  names(res) <- str_extract(resultsNames(dds)[2:15], "[:alnum:]+_[248]{2}")
  
  save(res, file = resFile)
} else {
  load(resFile)
}


res <- mapply(formatResults, 
              x = res, x_nm = names(res),
              SIMPLIFY = FALSE)

resLong <-
  Reduce(bind_rows, res) %>% 
  mutate(baseMean = round(baseMean, 4),
         log2FoldChange = round(log2FoldChange, 4),
         lfcSE = round(lfcSE, 4),
         pvalue = signif(pvalue, 4),
         padj   = signif(padj, 4))

resLong$ensembl_gene_id %>% unique %>% length
  
allNA <- 
  resLong %>% 
  group_by(ensembl_gene_id) %>% 
  summarize(anyNA = anyNA(padj),
              allNA = !any(!is.na(padj))) %>% 
  ungroup %>% 
  filter(allNA) %>% 
  pull(ensembl_gene_id)

resLongFilt <-
  resLong %>% 
  filter(!(ensembl_gene_id %in% allNA))

write_csv(resLongFilt, path = "MDD_RNAseq_DESeq2_Long.csv")

# save(res_RNA, res_RNA_df, 
#      file = sprintf("%s/MDD_RNAseq_shrunkenResults.Rdata", outDir))
# save(res_RNA_df, 
#      file = sprintf("%s/MDD_RNAseq_shrunkenResults_df.Rdata", outDir))
# 
# load(sprintf("%s/MDD_RNAseq_shrunkenResults.Rdata", outDir))
# 
# sig_genes <- lapply(res_RNA_df, function(X) {
#   X <- X %>%
#     filter(abs(log2FoldChange) >= 1.5,
#            padj <= 0.01) %>% 
#     arrange(log2FoldChange)
#   }
# )
# 
# sig_genes_long <- Reduce(bind_rows, sig_genes) %>% 
#   dplyr::select(experimentalCondition, hgnc_symbol, ensembl_gene_id, log2FoldChange, padj, everything())
# 
# write_csv(sig_genes_long, path = "../RNAseq/Data/MDD_RNAseq_DEGenesLong.csv")
#  