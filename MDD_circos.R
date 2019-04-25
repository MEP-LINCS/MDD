#create circos figures of MDD correlations
library(tidyverse)
library(circlize)
library(ComplexHeatmap)
library(gridBase)

load("MCF10A_RPPA_GCP_CYCIF_condition_correlation.Rdata")

prepare_data <- function(x){
  df <- x %>%
    data.frame() %>%
    rownames_to_column(var ="Condition") %>%
    mutate(Time = str_remove(Condition, ".*_"),
           Time = str_replace_all(Time, "^4$","1"),
           Time = str_replace_all(Time, "^8$","2"),
           Time = str_replace_all(Time, "^24$","3"),
           Time = str_replace_all(Time, "^48$","4"),
           Ligand = str_remove(Condition, "_.*"))
}

col_names <- c("grey50", "grey25","red3","royalblue3")

rppa <- prepare_data(cormats[["RPPA"]])
cycIF <- prepare_data(cormats[["CYCIF"]])
GCP <- prepare_data(cormats[["GCP"]])

circos.par("track.height" = 0.2)
circos.initialize(factors = as.factor(rppa$Ligand), xlim =c(0,4), sector.width = 20)

circos.track(ylim = c(0, 1), panel.fun = function(x, y) {
  xlim = CELL_META$xlim
  ylim = CELL_META$ylim
  breaks = seq(xlim[1], xlim[2], by = 1)
  n_breaks = length(breaks)
  circos.rect(breaks[-n_breaks], rep(ylim[1], n_breaks - 1),
              breaks[-1], rep(ylim[2], n_breaks - 1),
              col = col_names,
              border = NA)
  circos.text(CELL_META$xcenter, CELL_META$cell.ylim[2] + uy(5, "mm"),
              CELL_META$sector.index)
}) 
parseCorMatrix <- function(x){
  from_condition <- x$Condition
  to_condition <- names(x)[c(-1,-22,-23)]
  alpha <- x[c(-1,-22,-23)] %>%
    unlist()
  alpha <- (alpha * 16) %% 256
  alpha <- as.integer(alpha) %>%
    as.character.hexmode()
  col <- paste0("#000000",alpha)
  df <- tibble(From_Condition = x$Condition,
                    To_Condition =  names(x)[c(-1,-22,-23)],
                    From_sector = str_remove(from_condition, "_.*"),
                    To_sector = str_remove(to_condition, "_.*"),
                    From_timepoint = str_remove(from_condition, ".*_"),
                    To_timepoint = str_remove(to_condition, ".*_"),
                    Col = col) %>%
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
           From_timepoint_st = From_timepoint - .5,
           From_timepoint_end = From_timepoint + .5,
           To_timepoint_st = To_timepoint - .5,
           To_timepoint_end = To_timepoint + .5)
}

addLink <- function(x){
  # circos.link(x["From_sector"], c(as.numeric(x["From_timepoint_st"]), as.numeric(x["From_timepoint_end"])), x["To_sector"], c(as.numeric(x["To_timepoint_st"]), as.numeric(x["To_timepoint_end"])), h = 0.4, lwd = .01, col =x["Col"])
  circos.link(x["From_sector"], as.numeric(x["From_timepoint"]), x["To_sector"], as.numeric(x["To_timepoint"]), h = 0.4, lwd = 1, col =x["Col"])
}

plot.new()
circle_size = unit(1, "snpc")
pushViewport(viewport(x = 0.5, y = 1, width = circle_size, height = circle_size,
                      just = c("center", "top")))
par(omi = gridOMI(), new = TRUE)
circlize_plot()
upViewport()

draw(lgd_list, y = unit(1, "npc") - circle_size, just = "top")

title("RPPA Correlations")

for(i in 1:20){
  foo <- parseCorMatrix(rppa[i,])
  res <- apply(foo, 1, addLink)
}

circos.initialize(factors = as.factor(rppa$Ligand), xlim =c(0,4), sector.width = 20)
# create legend
lgd = Legend(at = c(4,8,24,48), type = "grid", 
             legend_gp = gpar(col = 4:5, lwd = 2),
             title_position = "topleft",
             title = "Time")
lgd_list <- packLegend(lgd)

circos.track(ylim = c(0, 1), panel.fun = function(x, y) {
  xlim = CELL_META$xlim
  ylim = CELL_META$ylim
  breaks = seq(xlim[1], xlim[2], by = 1)
  n_breaks = length(breaks)
  circos.rect(breaks[-n_breaks], rep(ylim[1], n_breaks - 1),
              breaks[-1], rep(ylim[2], n_breaks - 1),
              col = col_names,
              border = NA)
  circos.text(CELL_META$xcenter, CELL_META$cell.ylim[2] + uy(5, "mm"),
              CELL_META$sector.index)
}) 

title("cycIF Correlations")

for(i in 1:20){
  foo <- parseCorMatrix(cycIF[i,])
  res <- apply(foo, 1, addLink)
}

circos.initialize(factors = as.factor(rppa$Ligand), xlim =c(0,4), sector.width = 20)

circos.track(ylim = c(0, 1), panel.fun = function(x, y) {
  xlim = CELL_META$xlim
  ylim = CELL_META$ylim
  breaks = seq(xlim[1], xlim[2], by = 1)
  n_breaks = length(breaks)
  circos.rect(breaks[-n_breaks], rep(ylim[1], n_breaks - 1),
              breaks[-1], rep(ylim[2], n_breaks - 1),
              col = col_names,
              border = NA)
  circos.text(CELL_META$xcenter, CELL_META$cell.ylim[2] + uy(5, "mm"),
              CELL_META$sector.index)
}) 

draw(lgd_list, x = unit(4, "mm"), y = unit(4, "mm"), just = c("left", "bottom"))
title("GCP Correlations")

for(i in 1:20){
  foo <- parseCorMatrix(GCP[i,])
  res <- apply(foo, 1, addLink)
}
