#import required libraries
library(mavric)
library(vsn)
library(readr)
library(dplyr)
library(tidyverse)
library(rJava)
library(FactoMineR)
library(cluster)
library(DESeq2)
library(doParallel)
library(foreach)
library(SummarizedExperiment)
library(glmnet)
library(gtools)
library(msir)
library(scatterplot3d)
library(stats)
library(venneuler)
library(plot3D)
library(ggplot2)
library(factoextra)
library(reshape2)

#custom functions

#winsorize assay datasets
winsorize<-function(dataframe){
  lower_bound<-as.numeric(quantile(unlist(dataframe), probs=c(0.01)))
  upper_bound<-as.numeric(quantile(unlist(dataframe), probs=c(0.99)))
  
  current_windsored<-dataframe
  current_windsored[current_windsored<lower_bound]<-lower_bound
  current_windsored[current_windsored>upper_bound]<-upper_bound
  
  return(current_windsored)
}

#median center datasets
med_center_replicate<-function(dataframe, design_mat){
  design_mat_formatted<-design_mat %>%
    rownames_to_column('specimenID')
  
  dataframe_centered<-t(dataframe) %>%
    as.data.frame() %>%
    rownames_to_column('specimenID') %>%
    left_join(design_mat_formatted) %>%
    column_to_rownames('specimenID') %>%
    select(-experimentalTimePoint, -ligand) %>%
    group_by(replicate) %>%
    mutate_all(funs(. - median(.))) %>%
    ungroup() %>%
    mutate(specimenID=colnames(dataframe)) %>%
    select(-replicate) %>%
    column_to_rownames('specimenID') %>%
    t %>%
    as.data.frame()
  
  return(dataframe_centered)
}

#venn diagram (adapted from Mavric R package)
plotVars_color<-function (results, annotation, plotzero = TRUE, incvals = TRUE, 
                          plotve = TRUE, maxStress = 0.1) 
{
  if (!is.data.frame(annotation)) 
    annotation <- data.frame(annotation)
  results$pcs <- lapply(results$pcs, function(e) e[as.numeric(rownames(e)) > 
                                                     0, , drop = F])
  vars <- round(unlist(sapply(1:length(results$pcs), function(p) {
    if (nrow(results$pcs[[p]]) == 0) return(0)
    return(sum(results$pcavar[as.numeric(rownames(results$pcs[[p]]))] * 
                 sapply(as.numeric(rownames(results$pcs[[p]])), function(i) totVarEx(results$pcasp, 
                                                                                     annotation, i, p))))
  })))
  vn <- paste(colnames(annotation), paste(vars, "%", sep = ""), 
              sep = "\n")
  pcex <- sort(unique(unlist(lapply(results$pcs, rownames))))
  pcm <- cbind(matrix(unlist(apply(matrix(unlist(sapply(1:length(results$pcs), 
                                                        function(li) {
                                                          rv <- rep(0, length(pcex))
                                                          if (nrow(results$pcs[[li]]) > 0) {
                                                            e <- rownames(results$pcs[[li]])
                                                            for (i in 1:length(e)) rv[match(e[i], pcex)] <- results$pcavar[as.numeric(e[i])] * 
                                                                totVarEx(results$pcasp, annotation, as.numeric(e[i]), 
                                                                         li)
                                                          }
                                                          return(rv)
                                                        })), nrow = length(pcex), dimnames = list(pcex, vn)), 
                                   1, function(i) {
                                     nz <- head(order(i, decreasing = T), n = sum(i != 
                                                                                    0))
                                     if (length(nz) == 1) 
                                       return(rbind(i[nz], vn[nz]))
                                     return(rbind(c(diff(t(embed(i[nz], 2))), i[nz[length(nz)]]), 
                                                  sapply(1:length(nz), function(n) paste(vn[nz[1:n]], 
                                                                                         collapse = "&"))))
                                   })), nrow = 2), matrix(sapply(unlist(sapply(1:length(results$pcs), 
                                                                               function(li) if (nrow(results$pcs[[li]]) == 0) 
                                                                                 return(c(0, vn[li])))), function(i) i), nrow = 2))
  unacc <- 100 - sum(results$pcavar)
  unex <- 100 - (unacc + sum(as.numeric(pcm[1, ])))
  other.var <- c(unex, unacc)
  names(other.var) <- paste(c("Unexplained", "Discarded"), 
                            paste(as.character(round(c(unex, unacc))), "%", sep = ""), 
                            sep = "\n")
  pl.vars <- sort(c(tapply(as.numeric(unlist(pcm[1, ])), unlist(lapply(lapply(strsplit(unlist(pcm[2, 
  ]), "&"), sort), paste, collapse = "&")), sum), other.var), 
  decreasing = incvals)
  if (!plotzero) 
    pl.vars <- pl.vars[pl.vars > 0]
  if (!plotve) 
    return(pl.vars)
  pl.vars.c <- ceiling(pl.vars)
  repeat {
    pl.tmp <- abs(jitter(pl.vars.c))
    if (plotzero) 
      pl.tmp[pl.vars.c == 0] = 0
    pl.ve <- venneuler::venneuler(pl.tmp)
    if (pl.ve$stress < maxStress) 
      break
  }
  if (plotve) 
    fig<-venneuler(pl.tmp, col = c('red','green','purple','black','grey'))
    return(fig)
}

#variance calculations (from Mavric R package)
totVarEx <- function(pcsp, ann, csp, grp) {
  if(class(ann[,grp]) == "factor") return(discVarEx(pcsp, ann, csp, grp))
  return(contVarEx(ann[,grp,drop=F], pcsp[,csp]))
}

discVarEx <- function(pcsp, ann, csp, grp) {
  return(1 - sum(unlist(sapply(levels(ann[, grp]), function(l) {
    cmean <- mean(pcsp[ann[, grp] == l, csp]);
    sapply(pcsp[ann[, grp] == l, csp], function(i) sum((i - cmean)^2))
  })))/sum(scale(pcsp[, csp], scale = F)^2))
}

contVarEx <- function(x, y) apply(x, 2, function(j) quickcor(j, y)^2)

quickcor<-function (x, y) 
  .Call(stats:::C_cor, x, y, 4, FALSE)

#import level 3 datasets
rna_level3<-read_csv('../RNAseq/Data/MDD_RNAseq_Level3.csv')
rppa_level3<-read_csv('../RPPA/Data/MDD_RPPA_Level3.csv')
atac_level3<-read_csv('../ATACseq/Data/MDD_ATACseq_Level3.csv') 
cyclicIF_level3<-read_csv('../cycIF/Data/MDD_cycIF_Level3.csv')
l1000_level3<-read_csv('../L1000/Data/MDD_L1000_Level3.csv')
gcp_level3<-read_csv('../GCP/Data/MDD_GCP_Level3.csv') 

#import metadata
rna_meta<-read_csv('../rnaSeq/Metadata/MDD_RNAseq_sampleMetadata.csv')
atac_meta<-read_csv('../ATACseq/Metadata/MDD_ATACseq_sampleMetadata.csv')
l1000_meta<-read_csv('../L1000/Metadata/MDD_l1000_sampleMetadata.csv')
sample_meta<-read_csv('../Metadata/MDD_sample_annotations.csv')

#format data and metadata
rna_mat<-rna_level3[c(-1)]
rownames(rna_mat)<-t(rna_level3[,1])
filtered_rna_meta<-as.data.frame(rna_meta) %>%
  filter(RNAseq_QCpass==T) 

rppa_mat<-rppa_level3[c(-1)]
rownames(rppa_mat)<-t(rppa_level3[,1])
filtered_rppa_meta<-as.data.frame(sample_meta) %>%
  filter(RPPA_QCpass==T) 

atac_mat<-atac_level3[c(-1)] 
rownames(atac_mat)<-t(atac_level3[,1])
filtered_atac_meta<-as.data.frame(atac_meta) %>%
  filter(ATACseq_QCpass==T) 

cyclicIF_mat<-cyclicIF_level3 %>%
  filter(grepl('_int_med', feature)) %>%
  column_to_rownames('feature')
filtered_cyclicIF_meta<-as.data.frame(sample_meta) %>%
  filter(studySite=='HMS') 

l1000_mat<-l1000_level3[c(-1)]
rownames(l1000_mat)<-t(l1000_level3[,1])
filtered_l1000_meta<-as.data.frame(l1000_meta) %>%
  filter(outlierSample==F & collection=='C2') 

L1000_c2<- l1000_mat %>%
  select(.dots = filtered_l1000_meta$specimenID)
colnames(L1000_c2)<-filtered_l1000_meta$specimenID

#remove missing data from GCP dataset
gcp_missing_removed<-data.frame(t(apply(gcp_level3[,c(2:89)],1,f)))
gcp_mat<-gcp_missing_removed
rownames(gcp_mat)<-t(gcp_level3[,1])
colnames(gcp_mat)<-colnames(gcp_level3[,-1])
filtered_gcp_meta<-as.data.frame(sample_meta) %>%
  filter(GCP_QCpass==T) %>%
  filter(specimenID %in% colnames(gcp_mat)) 

#winsorize all data matrices at 99%
rna_winsorized<-winsorize(rna_mat)
rppa_winsorized<-winsorize(rppa_mat)
atac_winsorized<-winsorize(atac_mat)
cyclicIF_winsorized<-winsorize(cyclicIF_mat)
l1000_winsorized<-winsorize(L1000_c2)
gcp_winsorized<-winsorize(gcp_mat)

#create design matrices for variance estimation
rna_design_mat<-filtered_rna_meta %>%
  dplyr::select(ligand, experimentalTimePoint, replicate) %>%
  mutate_if(is.character, as.factor) 
rownames(rna_design_mat)<-t(filtered_rna_meta[,1])

rppa_design_mat<-filtered_rppa_meta %>%
  dplyr::select(ligand, experimentalTimePoint, replicate) %>%
  mutate_if(is.character, as.factor)
rownames(rppa_design_mat)<-t(filtered_rppa_meta[,1])

atac_design_mat<-filtered_atac_meta %>%
  select(ligand, experimentalTimePoint, replicate) %>%
  mutate_if(is.character, as.factor) 
rownames(atac_design_mat)<-t(filtered_atac_meta[,1])

cyclicIF_design_mat<-filtered_cyclicIF_meta %>%
  select(ligand, experimentalTimePoint, replicate) %>%
  mutate_if(is.character, as.factor) 
rownames(cyclicIF_design_mat)<-t(filtered_cyclicIF_meta[,1])

l1000_design_mat<-filtered_l1000_meta %>%
  select(ligand, experimentalTimePoint, replicate) %>%
  mutate_if(is.character, as.factor) 
rownames(l1000_design_mat)<-t(filtered_l1000_meta[,1])

gcp_design_mat<-filtered_gcp_meta %>%
  select(ligand, experimentalTimePoint, replicate) %>%
  mutate_if(is.character, as.factor) 
rownames(gcp_design_mat)<-t(filtered_gcp_meta[,1])

#median center all dataframes within replicate
rna_med_centered<-med_center_replicate(rna_winsorized, rna_design_mat)
rppa_med_centered<-med_center_replicate(rppa_winsorized, rppa_design_mat)
atac_med_centered<-med_center_replicate(atac_winsorized, atac_design_mat)
cyclicIF_med_centered<-med_center_replicate(cyclicIF_winsorized, cyclicIF_design_mat)
l1000_med_centered<-med_center_replicate(l1000_winsorized, l1000_design_mat)
gcp_med_centered<-med_center_replicate(gcp_winsorized, gcp_design_mat)

#create data matrices for mavric analysis
rna<-as.matrix(rna_med_centered)
dimnames(rna)<-list(rownames(rna), colnames(rna))

rppa<-as.matrix(rppa_med_centered)
dimnames(rppa)<-list(rownames(rppa), colnames(rppa))

atac<-as.matrix(atac_med_centered)
dimnames(atac)<-list(rownames(atac), colnames(atac))

cyclicIF<-as.matrix(cyclicIF_med_centered)
dimnames(cyclicIF)<-list(rownames(cyclicIF), colnames(cyclicIF))

l1000<-as.matrix(l1000_med_centered)
dimnames(l1000)<-list(rownames(l1000), colnames(l1000))

gcp<-as.matrix(gcp_med_centered)
dimnames(gcp)<-list(rownames(gcp), colnames(gcp))

#Run mavric on all assay matrices
contrastlist=list(list('experimentalTimePoint'))

evc_rna<-mavric::estVC(data=rna, annotation = rna_design_mat, contrastlist = contrastlist, sigcor = F, ntop = Inf, alpha = .05, effsize =.02,
                       discthresh = 1, prenorm = T, autosel = F, corMethod = c('pearson'), resampling = c('permutation'), scaleVar = F, minVar = .35)
rna_var_explained<-as.data.frame(plotVars(evc_rna, rna_design_mat)) %>%
  rownames_to_column('variance_type') %>%
  mutate(variance_type=str_remove_all(variance_type,'[0-9]+')) %>%
  mutate(variance_type=str_remove_all(variance_type,'%')) %>%
  mutate(variance_type=str_remove_all(variance_type,' ')) 
colnames(rna_var_explained)<-c('variance_type','RNAseq')

evc_rppa<-mavric::estVC(data=rppa, annotation = rppa_design_mat, contrastlist = contrastlist, sigcor = F, ntop = Inf, alpha = .05, effsize =.02,
                       discthresh = 1, prenorm = T, autosel = F, corMethod = c('pearson'), resampling = c('permutation'), scaleVar = F, minVar = .35)
rppa_var_explained<-as.data.frame(plotVars(evc_rppa, rppa_design_mat)) %>%
  rownames_to_column('variance_type') %>%
  mutate(variance_type=str_remove_all(variance_type,'[0-9]+')) %>%
  mutate(variance_type=str_remove_all(variance_type,'%')) %>%
  mutate(variance_type=str_remove_all(variance_type,' ')) 
colnames(rppa_var_explained)<-c('variance_type','RPPA')

evc_atac<-mavric::estVC(data=atac, annotation = atac_design_mat, contrastlist = contrastlist, sigcor = F, ntop = Inf, alpha = .05, effsize =.02,
                       discthresh = 1, prenorm = T, autosel = F, corMethod = c('pearson'), resampling = c('permutation'), scaleVar = F, minVar = .35)
atac_var_explained<-as.data.frame(plotVars(evc_atac, atac_design_mat)) %>%
  rownames_to_column('variance_type') %>%
  mutate(variance_type=str_remove_all(variance_type,'[0-9]+')) %>%
  mutate(variance_type=str_remove_all(variance_type,'%')) %>%
  mutate(variance_type=str_remove_all(variance_type,' ')) 
colnames(atac_var_explained)<-c('variance_type','ATACseq')

evc_cyclicIF<-mavric::estVC(data=cyclicIF, annotation = cyclicIF_design_mat, contrastlist = contrastlist, sigcor = F, ntop = Inf, alpha = .05, effsize =.02,
                       discthresh = 1, prenorm = T, autosel = F, corMethod = c('pearson'), resampling = c('permutation'), scaleVar = F, minVar = .35)
cyclicIF_var_explained<-as.data.frame(plotVars(evc_cyclicIF, cyclicIF_design_mat)) %>%
  rownames_to_column('variance_type') %>%
  mutate(variance_type=str_remove_all(variance_type,'[0-9]+')) %>%
  mutate(variance_type=str_remove_all(variance_type,'%')) %>%
  mutate(variance_type=str_remove_all(variance_type,' ')) 
colnames(cyclicIF_var_explained)<-c('variance_type','CyclicIF')

evc_l1000<-mavric::estVC(data=l1000, annotation = l1000_design_mat, contrastlist = contrastlist, sigcor = F, ntop = Inf, alpha = .05, effsize =.02,
                       discthresh = 1, prenorm = T, autosel = F, corMethod = c('pearson'), resampling = c('permutation'), scaleVar = F, minVar = .35)
l1000_var_explained<-as.data.frame(plotVars(evc_l1000, l1000_design_mat)) %>%
  rownames_to_column('variance_type') %>%
  mutate(variance_type=str_remove_all(variance_type,'[0-9]+')) %>%
  mutate(variance_type=str_remove_all(variance_type,'%')) %>%
  mutate(variance_type=str_remove_all(variance_type,' ')) 
colnames(l1000_var_explained)<-c('variance_type','L1000')

evc_gcp<-mavric::estVC(data=gcp, annotation = gcp_design_mat, contrastlist = contrastlist, sigcor = F, ntop = Inf, alpha = .05, effsize =.02,
                       discthresh = 1, prenorm = T, autosel = F, corMethod = c('pearson'), resampling = c('permutation'), scaleVar = F, minVar = .35)
gcp_var_explained<-as.data.frame(plotVars(evc_gcp, gcp_design_mat)) %>%
  rownames_to_column('variance_type') %>%
  mutate(variance_type=str_remove_all(variance_type,'[0-9]+')) %>%
  mutate(variance_type=str_remove_all(variance_type,'%')) %>%
  mutate(variance_type=str_remove_all(variance_type,' ')) 
colnames(gcp_var_explained)<-c('variance_type','GCP')


#compile results
total_variation_comparison<-list(atac_var_explained, rna_var_explained, cyclicIF_var_explained,
                                 gcp_var_explained,l1000_var_explained, rppa_var_explained) %>% 
  Reduce(function(dtf1,dtf2) full_join(dtf1, dtf2, by = 'variance_type'), .)
total_variation_comparison[is.na(total_variation_comparison)]<-0

#reassign low variance covariate combinations
total_variation_comparison_formatted<-total_variation_comparison %>%
  mutate(variance_type=gsub('&','/',variance_type)) %>%
  mutate(variance_type=gsub("\n", "", variance_type))

variation_reduced<-total_variation_comparison_formatted[rowSums(total_variation_comparison_formatted[2:7] >= 1) > 0, ]
variation_low<-total_variation_comparison_formatted[rowSums(total_variation_comparison_formatted[2:7]) < 1, ]


#reassign variances
variation_reduced[2,2:7]<-variation_reduced[2,2:7] + variation_low[2,2:7]/3 + variation_low[1,2:7]/2 
variation_reduced[3,2:7]<-variation_reduced[3,2:7] + variation_low[2,2:7]/3 + variation_low[1,2:7]/2 + variation_low[3,2:7]/2
variation_reduced[5,2:7]<-variation_reduced[5,2:7] + variation_low[2,2:7]/3 + variation_low[3,2:7]/2


comparison_data.m<-melt(variation_reduced, id.vars='variance_type')
comparison_data.m$variable<-factor(comparison_data.m$variable, levels=c('ATACseq', 'RNAseq', 'GCP', 'L1000', 'CyclicIF', 'RPPA'))


cov_cols <- c("Discarded" = "#808080",
                 "Unexplained" = "#D3D3D3",
                 "experimentalTimePoint/ligand" = '#0000CC',
                 "replicate" = "#F8766D",
                 "experimentalTimePoint" = "#7CAE00",
                 "ligand" = "#00BFC4")



ggplot(comparison_data.m, aes(x=variable,y=value, fill=factor(comparison_data.m$variance_type, levels=c('Discarded','Unexplained', 
                                                                                                        'experimentalTimePoint/ligand', 
                                                                                                        'replicate','experimentalTimePoint','ligand')))) +
  geom_bar(stat='Identity') + ylab('Variance Explained') + xlab('Assay') + labs(fill='Covariate') +coord_flip() + guides(fill = guide_legend(reverse = TRUE)) +
  scale_fill_manual(values=cov_cols) + theme_bw() + theme(panel.border = element_blank(), panel.grid.major = element_blank(),
                                                          panel.grid.minor = element_blank(), axis.line = element_blank())


#plot rppa venn diagram of variance contributions
rppa_venn<-plotVars_color(results = evc_rppa, annotation = rppa_design_mat, plotzero = TRUE, incvals = TRUE, 
                          plotve = TRUE, maxStress = 0.1)


plot(rppa_venn, col=c('#D3D3D3','#00BFC4', '#808080', '#7CAE00', '#F8766D'))

#plot rppa variance by pc (adapted from Mavric R package)
#get correlations by pc
results=evc_rppa
annotation=rppa_design_mat
plotzero = TRUE
incvals = FALSE
plotve = TRUE
maxStress = .1

#pull out pcs
results$pcs <- lapply(results$pcs, function(e) e[as.numeric(rownames(e)) > 0, , drop = F])

#make into columns
total_var_pcs<-as.data.frame(results$pcavar) %>%
  dplyr::select(totalvar='results$pcavar')
pcex <- sort(unique(unlist(lapply(results$pcs, rownames))))
pcex<-as.numeric(pcex)
total_var_pcs<-as.data.frame(total_var_pcs[c(1:length(pcex)),]) 
pc_order<-as.data.frame(pcex) %>%
  dplyr::select(PC=(pcex))

rv <- rep(0, length(pcex))
e <- rownames(results$pcs[[1]])
for(i in 1:length(e)) rv[match(e[i], pcex)] <- results$pcavar[as.numeric(e[i])]*totVarEx(results$pcasp, annotation, as.numeric(e[i]), 1)

ligand<-as.data.frame(rv) %>%
  dplyr::select(ligand=rv)

rv_2 <- rep(0, length(pcex))
e <- rownames(results$pcs[[2]])
for(i in 1:length(e)) rv_2[match(e[i], pcex)] <- results$pcavar[as.numeric(e[i])]*totVarEx(results$pcasp, annotation, as.numeric(e[i]), 2)

time <- as.data.frame(rv_2) %>%
  dplyr::select(time=rv_2)

rv_3 <- rep(0, length(pcex))
e <- rownames(results$pcs[[3]])
for(i in 1:length(e)) rv_3[match(e[i], pcex)] <- results$pcavar[as.numeric(e[i])]*totVarEx(results$pcasp, annotation, as.numeric(e[i]), 3)

replicate <- as.data.frame(rv_3) %>%
  dplyr::select(replicate=rv_3)

#bind columns of variances together
by_cov<-bind_cols(time,ligand, replicate,pc_order) %>%
  arrange((PC))

by_cov_tot<-bind_cols(by_cov,total_var_pcs) %>%
  dplyr::select(total_variance=`total_var_pcs[c(1:length(pcex)), ]`, replicate, time, ligand, PC) %>%
  mutate(Unexplained=total_variance - replicate - time - ligand ) 

#set unexplained variances to unexplained
by_cov_tot$Unexplained[by_cov_tot$Unexplained<0]<-0

by_cov.m <- melt(by_cov_tot,id.vars = "PC")

cov_cols <- c(
              "Unexplained" = "#D3D3D3",
              "replicate" = "#F8766D",
              "time" = "#7CAE00",
              "ligand" = "#00BFC4")

#plot RPPA variance contributions by PC
rppa_pc<-ggplot(aes(y=value, x=PC, fill = variable), data = subset(by_cov.m, variable!='total_variance' & PC<8)) +
  geom_bar(stat = 'identity') + ylab('Percent Variance Explained') + labs(title='RPPA Variance Explained by PC') +
  coord_flip() + theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
                        panel.background = element_blank(), axis.line = element_line(colour = "grey")) + scale_fill_manual(values=cov_cols)



