#RPPA line graphs of MAP pathway members
library(tidyverse)

ligand_cols <- c("CTRL" = "#7A4A2A",
                 "PBS" = "#8dd3c7",
                 "HGF" = "#80b1d3",
                 "OSM" = "#fdb462",
                 "EGF" = "#fb8072",
                 "BMP2+EGF" = "#b3de69",
                 "IFNG+EGF" = "#bebada",
                 "TGFB+EGF" = "#ffd92f")

get_assay_values <- function(assay){
  #browser()

    level3File <- "../RPPA/Data/MDD_RPPA_limma_logFCTable.csv"
 
  MDDannoFile <- "../Metadata/MDD_sample_annotations.csv"
  

    mat <- read_csv(level3File,
                    col_types = cols(
                      antibody = col_character(),
                      experimentalCondition = col_character(),
                      coefficient = col_double(),
                      pvalue = col_double(),
                      adjusted_pvalue = col_double()
                    ))
    
    meta <- NULL

  return(list(mat = mat, meta = meta))
}


plot_line_graphs <- function(df, md, ligand_cols, fc_thresh = 1, assay_name){
  
  windsorize_probs <- get_win_probs(df, .01, .99)
  
  df_med <- df %>%
    prep_data(md = md) %>%
    gather("specimenName", "value", -feature) %>%
    inner_join(md, by = "specimenName") %>%
    select(ligand, secondLigand, experimentalTimePoint, replicate, collection, feature, value) %>%
    filter(!ligand == "EGF") %>%
    rename(Time = experimentalTimePoint,
           Ligand = ligand,
           EGF = secondLigand,
           Replicate = replicate,
           Collection = collection) %>%
    mutate( Ligand = str_replace(Ligand, "BMP2", "BMP2+EGF"),
            Ligand = str_replace(Ligand, "IFNG", "IFNG+EGF"),
            Ligand = str_replace(Ligand, "TGFB", "TGFB+EGF"),
            Ligand = factor(Ligand, levels = names(ligand_cols))) %>%
    group_by(feature, Collection, Time, Ligand) %>%
    summarise(value = median(value)) %>%
    ungroup() 
  
  fc_set <- df_med %>%
    group_by(Time, Ligand) %>%
    filter(abs(value) >= fc_thresh) %>%
    ungroup() %>%
    select(feature) %>%
    distinct() %>%
    inner_join(df_med, by = "feature") %>%
    mutate(Ligand = factor(Ligand, levels = names(ligand_cols)))
  
  title_suffix <- ""
  p <- ggplot(fc_set, aes(x=Time, y=value, colour=Ligand))+
    geom_line()+
    labs(title=paste0("Line Graphs",title_suffix),
         x = "Time (hours)",
         y="Intensity (AU)") +
    scale_x_continuous(breaks = c(0,8,24,48)) +
    scale_color_manual(values = ligand_cols) +
    theme(axis.text.x = element_text(angle=90),
          axis.text.y=element_blank(),
          axis.ticks.y=element_blank(),
          strip.text = element_text(size=10, hjust=0),
          panel.background = element_rect(fill = NA),
          panel.grid.major = element_line(colour = "gray95"))+
    facet_wrap(~feature, scales = "free_y")
  p
  
}
fc_thresh = 0

MAP_members <- read_tsv("~/Downloads/Participating Molecules [R-HSA-450282].tsv")
RPPA_data <-get_assay_values("RPPA")$mat

all_data <- RPPA_data %>%
  mutate(Time = str_remove(experimentalCondition, ".*_"),
         Time = as.integer(Time),
         Ligand = str_remove(experimentalCondition, "_.*"), 
         Ligand = str_replace(Ligand, "BMP2", "BMP2+EGF"),
         Ligand = str_replace(Ligand, "IFNG", "IFNG+EGF"),
         Ligand = str_replace(Ligand, "TGFB", "TGFB+EGF"),
         Ligand = factor(Ligand, levels = names(ligand_cols)))

all_data_with_T0 <- all_data %>%
  filter(Time == 1) %>%
  distinct() %>%
  mutate(coefficient = 0,
         Time = 0) %>%
    bind_rows(all_data)
  
small_data <- all_data %>%
  filter(antibody %in% head(antibody, n = 60))

p <- ggplot(all_data_with_T0, aes(x=Time, y=coefficient, colour=Ligand))+
  geom_line()+
  geom_point(size = .2)+
  labs(title=paste0("Line Graphs, T0 DE values"),
       x = "Time (hours)",
       y="abundance") +
  scale_x_continuous(breaks = c(0,8,24,48)) +
  scale_color_manual(values = ligand_cols) +
  theme(text = element_text(size = 6, colour = "black"),
        # axis.text.y = element_text(color = "black"),
        # axis.text.x = element_text(color = "black"),
        axis.line = element_line(size = .3),
        #axis.text = element_text(angle=0, colour = "black"),
        axis.text.y=element_blank(),
        axis.ticks.y=element_blank(),
        legend.key = element_blank(),
        legend.key.height = unit(.2, 'cm'),
        #strip.text = element_text(size=12, hjust=0, color = "black"),
        strip.background = element_blank(),
        panel.background = element_rect(fill = NA),
        panel.grid.major = element_blank(),
        legend.spacing.y = unit(.01, 'cm')
        #panel.spacing.y = unit(1.1, "lines")
  ) +
  facet_wrap(~antibody, scales = "free_y",nrow = 12)
p

pdf("Grant_figure_MDD_RPPA_T0_DE_line_graphs", width = 24, height = 12)
print(p)
dev.off()


MAP_pathway <- read_csv("../RPPA/Metadata/MDD_RPPA_antibodyAnnotations.csv") %>%
  filter(str_detect(Pathway, "MAP"),
         !Protein == "P90RSKPT359S363") %>%
  mutate(antibody = str_remove(MDD, "-R-V")) %>%
  left_join(all_data, by = c("antibody" = "antibody")) %>%
  mutate(Time = str_remove(experimentalCondition, ".*_"),
         Time = as.integer(Time),
         Ligand = str_remove(experimentalCondition, "_.*"),
         Ligand = str_replace(Ligand, "BMP2", "BMP2+EGF"),
         Ligand = str_replace(Ligand, "IFNG", "IFNG+EGF"),
         Ligand = str_replace(Ligand, "TGFB", "TGFB+EGF"),
         Ligand = factor(Ligand, levels = names(ligand_cols)))
MAP_pathway_with_T0 <-  MAP_pathway %>%
  filter(Time == 1) %>%
  distinct() %>%
  mutate(coefficient = 0,
         Time = 0) %>%
  bind_rows(MAP_pathway)
p <- ggplot(MAP_pathway_with_T0, aes(x=Time, y=coefficient, colour=Ligand))+
  geom_line()+
  geom_point(size = .6)+
  labs(title=paste0("Line Graphs, T0 DE values"),
       x = "Time (hours)",
       y="abundance") +
  scale_x_continuous(breaks = c(0,8,24,48)) +
  scale_color_manual(values = ligand_cols) +
  theme(text = element_text(size = 6, colour = "black"),
        axis.text.y = element_text(color = "black"),
        axis.text.x = element_text(color = "black"),
        axis.line = element_line(size = .3),
        #axis.text = element_text(angle=0, colour = "black"),
        #axis.text.y=element_blank(),
        #axis.ticks.y=element_blank(),
        legend.key = element_blank(),
        legend.key.height = unit(.2, 'cm'),
        #strip.text = element_text(size=12, hjust=0, color = "black"),
        strip.background = element_blank(),
        panel.background = element_rect(fill = NA),
        panel.grid.major = element_blank(),
        legend.spacing.y = unit(.01, 'cm')
        #panel.spacing.y = unit(1.1, "lines")
  ) +
  facet_wrap(~Protein, scales = "free_y",nrow = 2)
p
pdf("Grant_figure_MDD_RPPA_T0_DE_MAP_pathway_line_graphs",height = 3)
print(p)
dev.off()

