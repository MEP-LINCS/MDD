#MDD GCP analysis

library(tidyverse)
library(scales)
library(RColorBrewer)
library(classInt)
#library(CePa)
#library(cmapR)

diverge.color <- function(data,pal_choice="RdGy",centeredOn=0){
  nHalf=50
  Min <- min(data,na.rm=TRUE)
  Max <- max(data,na.rm=TRUE)
  Thresh <- centeredOn
  pal<-brewer.pal(n=11,pal_choice)
  rc1<-colorRampPalette(colors=c(pal[1],pal[2]),space="Lab")(10)
  for(i in 2:10){
    tmp<-colorRampPalette(colors=c(pal[i],pal[i+1]),space="Lab")(10)
    rc1<-c(rc1,tmp)
  }
  rb1 <- seq(Min, Thresh, length.out=nHalf+1)
  rb2 <- seq(Thresh, Max, length.out=nHalf+1)[-1]
  rampbreaks <- c(rb1, rb2)
  cuts <- classIntervals(data, style="fixed",fixedBreaks=rampbreaks)
  return(list(cuts,rc1))
}


# foo <- read.gct("GCP_Data/GCP_MCF10a_EGF48hr_BatchNormalized.gct")
# foo <- read.gctx.meta("GCP_Data/GCP_MCF10a_raw.gct")
ProbeMetadata <- read_delim("GCP_Data/ProbeMetadata.txt", 
                            "\t", escape_double = FALSE, trim_ws = TRUE)

GCP_Data <- read_tsv("GCP_Data/GCP_MCF10a_log2_noProbeMetadata.txt", skip = 2) %>%
  t() %>%
  as.tibble()

colnames(GCP_Data) <- GCP_Data[1,]
GCP_Data <- GCP_Data %>%
  slice(-1) %>%
  mutate_at(vars(matches("^H3K")), as.numeric) %>%
  filter(!is.na(id))


w1<- GCP_Data %>%
  filter(pert_batch_internal_compound_enumerator == "W1")
w2<- GCP_Data %>%
  filter(pert_batch_internal_compound_enumerator == "W2")
w3<- GCP_Data %>%
  filter(pert_batch_internal_compound_enumerator == "W3")
w4<- GCP_Data %>%
  filter(pert_batch_internal_compound_enumerator == "W4")
w5<- GCP_Data %>%
  filter(pert_batch_internal_compound_enumerator == "W5")

EGF_values <- GCP_Data %>%
  select(matches("^H3K|pert_batch_internal_compound_enumerator|pert_iname|pert_time$")) %>%
  filter(pert_iname == "EGF") %>%
  gather(probe_id,value,-pert_batch_internal_compound_enumerator,-pert_iname,-pert_time)

#identify the NA values
EGF_NAs <- EGF_values %>%
  filter(is.na(value)) %>%
  select(pert_iname, pert_time, probe_id, pert_batch_internal_replicate) 

#calculate the medians from the replicates of the NA values
EGF_imputed <- EGF_NAs %>%
  select(-pert_batch_internal_replicate) %>%
  left_join(EGF_values, by = c("pert_iname", "pert_time", "probe_id")) %>%
  distinct() %>%
  group_by(pert_iname, pert_time, probe_id) %>%
  summarise(value = median(value, na.rm = TRUE)) %>%
  ungroup() %>%
  right_join(EGF_NAs, by = c("pert_iname", "pert_time", "probe_id"))

#remove the NA values and add in the imputed values
EGF_values <- EGF_values %>%
  filter(!is.na(value)) %>%
  bind_rows(EGF_imputed) %>%
  rename(EGF_value = value) %>%
  select(-pert_iname)

#Add the EGF values and normalize to them
EGF_norm<- GCP_Data %>%
  filter(!pert_time == "0") %>%
  select(matches("^H3K|replicate|pert_iname|pert_time$")) %>%
  gather(probe_id,value,-pert_batch_internal_replicate,-pert_iname,-pert_time) %>%
  left_join(EGF_values, by = c("pert_time", "probe_id", "pert_batch_internal_replicate")) %>%
  mutate(value_EGF_norm = value-EGF_value) %>%
  select(-value, -EGF_value) %>%
  mutate(pert_iname = toupper(pert_iname),
         condition = paste(pert_iname, pert_time, pert_batch_internal_replicate, sep = "_"),
         probe_id = factor(probe_id, levels = unique(probe_id))) %>%
  spread(key = probe_id, value = value_EGF_norm) %>%
  mutate(condition =  factor(condition, levels = paste(rep(c("PBS","BMP2", "IFNG", "TGFB", "HGF", "OSM", "EGF"), each = 20), rep(c(4,8, 24, 48), each = 5), 1:5, sep = "_"))) %>%
  arrange(condition)

GCP_metadata <- read_csv("MCF10a_GCP_data_with_OHSU_metadata_vocabulary.csv") %>%
  select(-contains("GMa")) %>%
  slice(-1:-23)


lowProb <- .02
upProb <- .98
m <- EGF_norm %>%
  select(-matches("pert|condition")) %>%
  unlist() %>%
  matrix(nrow = nrow(EGF_norm), byrow = FALSE) %>%
  squish(quantile(.,probs=c(lowProb, upProb), na.rm = TRUE))

rownames(m) <- unique(EGF_norm$condition, incomparables = NA)
colnames(m) <- EGF_norm %>%
  select(-matches("pert|condition")) %>%
  colnames()

m <- t(m)
ccols <- diverge.color(unlist(m), pal_choice = "RdBu", centeredOn = 0)

#data with columns ordered by samples
hm <- heatmap(m,
              Colv=NA,
              #ColSideColors = colAnBar,
              scale="none",
              col=ccols[[2]][length(ccols[[2]]):1],
              breaks=ccols[[1]]$brks,
              margins = c(5,10),
              cexRow = .8,
              cexCol = .8)#data with columns ordered by samples

hm <- heatmap(m,
              Colv=NA,
              Rowv = NA,
              #ColSideColors = colAnBar,
              scale="none",
              col=ccols[[2]][length(ccols[[2]]):1],
              breaks=ccols[[1]]$brks,
              margins = c(5,10),
              cexRow = .8,
              cexCol = .8)

