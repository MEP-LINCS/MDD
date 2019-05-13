#create circos figures of MDD correlations
library(tidyverse)
library(circlize)
library(ComplexHeatmap)
library(gridBase)

parseCorMatrix <- function(x){
  alpha_factor <- 255
  from_condition <- x$Condition
  to_condition <- names(x)[c(-1,-22,-23)]
  #remove metadata from correlations
  cors <- x[c(-1,-22,-23)] %>%
    unlist()
  #remove correlations within the ligand conditions
  cors[str_detect(names(cors), x[["Ligand"]])] <- 0
  #Remove negative correlations
  cors[cors<0] <- 0
  #scale alpha based on correlations
  alpha <- (cors * alpha_factor) %% 256
  alpha <- as.integer(alpha) %>%
    as.character.hexmode()
  col <- paste0(col_ligands[x$Ligand],alpha)
  df <- tibble(From_Condition = x$Condition,
               To_Condition =  names(x)[c(-1,-22,-23)],
               From_sector = str_remove(from_condition, "_.*"),
               To_sector = str_remove(to_condition, "_.*"),
               From_timepoint = str_remove(from_condition, ".*_"),
               To_timepoint = str_remove(to_condition, ".*_"),
               Col = col,
               Cors = as.numeric(cors)) %>%
    mutate(From_timepoint = str_replace(From_timepoint,"^4$",".5"),
           From_timepoint = str_replace(From_timepoint,"^8$","1.5"),
           From_timepoint = str_replace(From_timepoint,"^24$","2.5"),
           From_timepoint = str_replace(From_timepoint,"^48$","3.5"),
           From_timepoint = as.numeric(From_timepoint),
           To_timepoint = str_replace(To_timepoint,"^4$",".5"),
           To_timepoint = str_replace(To_timepoint,"^8$","1.5"),
           To_timepoint = str_replace(To_timepoint,"^24$","2.5"),
           To_timepoint = str_replace(To_timepoint,"^48$","3.5"),
           To_timepoint = as.numeric(To_timepoint),
           From_timepoint_st = From_timepoint - .45,
           From_timepoint_end = From_timepoint + .45,
           To_timepoint_st = To_timepoint - .45,
           To_timepoint_end = To_timepoint + .45)
}

circlize_plot <- function() {
  circos.initializeWithIdeogram(plotType = NULL)
  
  bed = generateRandomBed(nr = 300)
  bed = generateRandomBed(nr = 300, nc = 2)
  circos.genomicTrackPlotRegion(bed,
                                panel.fun = function(region, value, ...) {
                                  circos.genomicPoints(region, value, cex = 0.5, pch = 16, col = 2:3, ...)
                                })
  
  bed = generateRandomBed(nr = 500, nc = 2)
  circos.genomicTrackPlotRegion(bed,
                                panel.fun = function(region, value, ...) {
                                  circos.genomicLines(region, value, col = 4:5, ...)
                                })
  
  bed1 = generateRandomBed(nr = 100)
  bed1 = bed1[sample(nrow(bed1), 20), ]
  bed2 = generateRandomBed(nr = 100)
  bed2 = bed2[sample(nrow(bed2), 20), ]
  
  circos.genomicLink(bed1, bed2, col = col_fun(bed1[[4]]))
  
  circos.clear()
}
addLink <- function(x){
  cor_thresh <- .45
  if(x["Cors"] > cor_thresh) {
    lwd <- max(0, as.numeric(x["Cors"])-cor_thresh)*30
    #cat("lwd ",lwd, "col ", x["Col"])
    circos.link(x["From_sector"], as.numeric(x["From_timepoint"]), x["To_sector"], as.numeric(x["To_timepoint"]), h = 0.4, lwd = lwd, col =x["Col"])
  }
}

prepare_data <- function(x){
  df <- x %>%
    data.frame() %>%
    rownames_to_column(var ="Condition") %>%
    mutate(Time = str_remove(Condition, ".*_"),
           Time = str_replace_all(Time, "^4$","1"),
           Time = str_replace_all(Time, "^8$","2"),
           Time = str_replace_all(Time, "^24$","3"),
           Time = str_replace_all(Time, "^48$","4"),
           Ligand = str_remove(Condition, "_.*")
           )
}

#return the default ggplot colors
gg_color_hue <- function(n) {
  hues = seq(15, 375, length = n + 1)
  hcl(h = hues, l = 65, c = 100)[1:n]
}

draw_one_track <- function(title_text){
  circos.par("track.height" = 0.2,
             "points.overflow.warning" = FALSE)
  circos.initialize(factors = as.factor(rppa$Ligand), xlim =c(0,4), sector.width = 20)

  circos.track(ylim = c(0, 1), panel.fun = function(x, y) {
    xlim = CELL_META$xlim
    ylim = CELL_META$ylim
    breaks = seq(xlim[1], xlim[2], by = 1)
    n_breaks = length(breaks)
    circos.rect(breaks[-n_breaks], rep(ylim[1], n_breaks - 1),
                breaks[-1], rep(ylim[2], n_breaks - 1),
                col = track_col_names,
                border = NA)
    circos.text(CELL_META$xcenter,
                CELL_META$cell.ylim[2] + uy(3, "mm"),
                CELL_META$sector.index,
                col = c("black","red","black", "red", "red"))
  }) 
  
  title(title_text,
        outer = TRUE,
        cex.main = 1,
        line = -1)
  
  legend(-1.1,
         -.9,
         "TGFB, BMP2 and IFNG\n are paired with EGF",
         bty = "n",
         cex = .6)
}

######end of functions

load("Data/MCF10A_RPPA_GCP_CYCIF_condition_correlation.Rdata")

rppa <- prepare_data(cormats[["RPPA"]])
#cycIF <- prepare_data(cormats[["CYCIF"]])
#GCP <- prepare_data(cormats[["GCP"]])

#Assign ggplot colors to the ligands based on their order in the RPPA dataset
col_ligands <- gg_color_hue(5)
names(col_ligands) <- unique(rppa$Ligand)
#Assign colors to conditions in the track
track_col_names <- c("grey50", "grey25","red3","royalblue3")

pdf("CircosPlots.pdf",width = 7, height = 7)
# for(ligand in unique(rppa$Ligand)){
#   res <- draw_one_track(title_text = "RPPA Correlations")
#   df <- rppa %>%
#     filter(Ligand == ligand)
#   for(i in 1:nrow(df)){
#     foo <- parseCorMatrix(df[i,])
#     res <- apply(foo, 1, addLink)
#   }
#   circos.clear()
# }
res <- draw_one_track(title_text = "RPPA Correlations")
df <- rppa
for(i in 1:nrow(df)){
  foo <- parseCorMatrix(df[i,])
  res <- apply(foo, 1, addLink)
}
circos.clear()

# for(ligand in unique(cycIF$Ligand)){
#   res <- draw_one_track(title_text = "cycIF Correlations")
#   df <- rppa %>%
#     filter(Ligand == ligand)
#   for(i in 1:nrow(df)){
#     foo <- parseCorMatrix(df[i,])
#     res <- apply(foo, 1, addLink)
#   }
#   circos.clear()
# }
res <- draw_one_track(title_text = "cycIF Correlations")
df <- cycIF
for(i in 1:nrow(df)){
  foo <- parseCorMatrix(df[i,])
  res <- apply(foo, 1, addLink)
}
circos.clear()

res <- draw_one_track(title_text = "GCP Correlations")
df <- GCP
for(i in 1:nrow(df)){
  foo <- parseCorMatrix(df[i,])
  res <- apply(foo, 1, addLink)
}
circos.clear()
dev.off()
