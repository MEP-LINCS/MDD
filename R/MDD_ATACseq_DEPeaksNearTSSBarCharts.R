# The purpsoe of this script is to create quick plots examining new peaks in TSS sites in ATACseq data

library(tidyverse)
sigDir    <- "../ATACseq/Data/signatures/ctrl_0/fdr_05/anno/"
annoFile  <- "../ATACseq/Metadata/MDD_ATACseq_sampleMetadata.csv"
colScript <- "MDD_importColors_default.R"

outDir <- '../plots/ATACseq_DEPeaks_barCharts'

###############################################################################

if(!grepl("R$", getwd())) {
  setwd("R")
}

source(colScript)

###############################################################################

sampleAnno <- 
  read.csv(annoFile, stringsAsFactors = FALSE) %>% 
  distinct(experimentalCondition, .keep_all = TRUE) %>% 
  dplyr::select(experimentalCondition, ligand, experimentalTimePoint)

sigs <- lapply(list.files(sigDir, full.names = TRUE), 
               read.csv, 
               header = TRUE, stringsAsFactors = FALSE)
names(sigs) <- str_extract(list.files(sigDir), "[:alnum:]+[_][248]+")

sigs <-
  mapply(function(signature, condition) {
    signature <- 
      signature %>% 
      dplyr::select(seqnames, start, end, Conc_ctrl_0, Fold, FDR, annotation, geneId, distanceToTSS) %>% 
      mutate(experimentalCondition = condition) %>% 
      left_join(sampleAnno)
  }, 
  signature = sigs, condition = names(sigs), SIMPLIFY = FALSE)

sigLong <- Reduce(bind_rows, sigs)

if (!dir.exists(outDir)) {
  dir.create(outDir, recursive = TRUE)
}

pdf(sprintf("%s/MDD_ATACseq_DEPeaks_nearTSS.pdf", outDir), height = 8, width = 8)
sigLong %>% 
  filter(abs(distanceToTSS) < 1000) %>% 
  filter(Conc_ctrl_0 < 2) %>% 
  filter(Fold > 1) %>% 
  group_by(experimentalCondition, ligand, experimentalTimePoint) %>%
  summarize(n_peaks = n()) %>% 
  ggplot(aes(experimentalCondition, n_peaks, fill = ligand)) +
  geom_col() +
  scale_fill_manual(values = as.character(col$ligand)) +
  theme(axis.text.x = element_text(angle = 90, hjust = 1)) +
  ggtitle(label = "Newly accessible peaks near TSS",
          subtitle = "\nFDR < 0.05; Fold > 1\n|DistanceToTSS| < 1000 bp\nConc_ctrl_0 < 2")

sigLong %>% 
  filter(abs(distanceToTSS) < 1000) %>% 
  filter(Fold > 1) %>% 
  group_by(experimentalCondition, ligand, experimentalTimePoint) %>%
  summarize(n_peaks = n()) %>% 
  ggplot(aes(experimentalCondition, n_peaks, fill = ligand)) +
  geom_col() +
  scale_fill_manual(values = as.character(col$ligand)) +
  theme(axis.text.x = element_text(angle = 90, hjust = 1)) +
  ggtitle(label = "Increased accessibility in peaks near TSS",
          subtitle = "\nFDR < 0.05; Fold > 1\n|DistanceToTSS| < 1000 bp")

sigLong %>% 
  filter(abs(distanceToTSS) < 1000) %>% 
  filter(abs(Fold) > 1) %>% 
  group_by(experimentalCondition, ligand, experimentalTimePoint) %>%
  summarize(n_peaks = n()) %>% 
  ggplot(aes(experimentalCondition, n_peaks, fill = ligand)) +
  geom_col() +
  scale_fill_manual(values = as.character(col$ligand)) +
  theme(axis.text.x = element_text(angle = 90, hjust = 1)) +
  ggtitle(label = "DE peaks near TSS",
          subtitle = "\nFDR < 0.05; |Fold| > 1\n|DistanceToTSS| < 1000 bp")

sigLong %>% 
  filter(abs(Fold) > 1) %>% 
  group_by(experimentalCondition, ligand, experimentalTimePoint) %>%
  summarize(n_peaks = n()) %>% 
  ggplot(aes(experimentalCondition, n_peaks, fill = ligand)) +
  geom_col() +
  scale_fill_manual(values = as.character(col$ligand)) +
  theme(axis.text.x = element_text(angle = 90, hjust = 1)) +
  ggtitle("DE peaks\nFDR < 0.05; |Fold| > 1")
dev.off()

