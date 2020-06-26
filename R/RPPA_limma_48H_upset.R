# The purpose of this script is to perform upset analysis on RPPA data
# from the MCF10A common project. Analysis will be performed on 48-hour
# conditions.

# Starting with RPPA limma analysis
library(UpSetR)
source("MDD_RPPA_limma.R")

RNA_resultsFile <- "../RNAseq/misc/MDD_RNAseq_shrunkenResults_df.Rdata"

outDirPlots <- "../plots/RPPA_limma_upset"
outDirData <- "../RPPA/Data/limma_48H"

###############################################################################

# Filtering results to 48-hour time points only
colnames(res) <- str_remove(colnames(res), "experimentalCondition")

conditionsOfInterest <-
  RPPA.meta %>% 
  filter(experimentalTimePoint == 48) %>% 
  pull(experimentalCondition) %>% 
  as.character() %>% 
  unique

res48_up <- 
  res[, conditionsOfInterest] %>% 
  data.frame %>% 
  rownames_to_column("antibody")

res48_up[res48_up == -1] <- 0

res48_down <- 
  res[, conditionsOfInterest] %>% 
  data.frame %>% 
  rownames_to_column("antibody")

res48_down[res48_down == 1] <- 0
res48_down[res48_down == -1] <- 1

# Writing plots to file

if (!dir.exists(outDirPlots)) {
  dir.create(outDirPlots)
}

pdf(sprintf("%s/MDD_RPPA_48H_Up_Upset.pdf", outDirPlots), height = 8, width = 8, useDingbats = FALSE)
upset(res48_up,
      sets = rev(conditionsOfInterest),
      nintersects = NA,
      text.scale = 2,
      order.by = "freq",
      keep.order = TRUE)
upset(res48_up,
      sets = rev(conditionsOfInterest),
      nintersects = 15,
      text.scale = (20/7),
      order.by = "freq",
      keep.order = TRUE)
dev.off()

pdf(sprintf("%s/MDD_RPPA_48H_Down_Upset.pdf", outDirPlots), height = 8, width = 8, useDingbats = FALSE)
upset(res48_down, 
      sets = rev(conditionsOfInterest), 
      nintersects = NA, 
      text.scale = 2,
      order.by = "freq",
      keep.order = TRUE)
upset(res48_down, 
      sets = rev(conditionsOfInterest), 
      nintersects = 10, 
      text.scale = (24/7),
      order.by = "freq",
      keep.order = TRUE)
dev.off()


# Writing intersections to file

if (!dir.exists(outDirData)) {
  dir.create(outDirData)
}

res48_toWrite <- 
  res[, conditionsOfInterest] %>% 
  data.frame %>% 
  rownames_to_column("antibody")

write.csv(res48_toWrite,
          sprintf("%s/MDD_RPPA_48H_sigProteins_lfc05_p001.csv", outDirData),
          row.names = FALSE,
          quote = FALSE)


