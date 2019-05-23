#Author: Mark Dane, copyright 2015-2019
library(tidyverse)
library(readxl)
#library(scales)

localMinima <- function(x, probs=c(.2,.8)){
  #Finds the local minima between the probs quantiles
  #x numeric vector
  #probs interval limits on where to search for the minima
  h <- hist(x,breaks=200, plot=FALSE)
  if(length(h$mids)<2) return(max(x))
  f <- approxfun(h$mids, h$counts)
  o <- optimise(f, interval=quantile(x, probs, na.rm=TRUE))
  if(length(o)>2) stop()
  return(o$minimum)
}

gateOnlocalMinima <- function(x, probs=c(.2,.8), ...){
  if(anyNA(x)) {
    thresh <- NA
  } else {  thresh <- localMinima(x, probs, ...)}
  cluster <- rep.int(as.integer(1),times=length(x))
  cluster[x>thresh] <- as.integer(2)
  return(cluster)
}

wellAN <- function(nrRows,nrCols){
  if(nrRows>702)stop("Too many rows to convert. Well alphanumerics are limited to 2 letters")
  Well=unlist(c(lapply(1:nrRows, function(x){
    paste0(paste0(LETTERS[(x-1)%/%26],LETTERS[(1+(x-1)%%26)]), lapply(1:nrCols,function(x){
      if(x<10) x<-paste0("0",x)
      return(as.character(x))
    }))
  })
  ))
  return(Well)
}
  
readMetadata <-function(xlsFile){
  sheetList<-sapply(excel_sheets(xlsFile), read_excel, path = xlsFile, simplify = FALSE)
  nrRows<-dim(sheetList[[1]])[1]
  nrCols<-dim(sheetList[[1]])[2] - 1
  nrWells=nrRows*nrCols
  sheetDF<-lapply(sheetList,function(df,nrCols){
    #create a dataframe from all rows and columns of each sheet
    dfM<-matrix(t(df[,2:(nrCols+1)]),byrow=TRUE)
  },nrCol = nrCols) %>%
    as_tibble() %>%
    mutate(WellIndex=1:nrWells,
           Well = wellAN(nrRows, nrCols))
  return(sheetDF)
}

#reduce the numeric values to 4 significant digits
shorten <- function(x){
  if(class(x)=="numeric") x <- signif(x,4)
  return(x)
}

fileNames <- dir("../IF/Data/",pattern = "Cyto|Main", full.names = TRUE)
barcodes <- str_extract(fileNames, "LI8[[:alnum:]]*") %>%
  unique()

#Get and clean the cell level raw data
rd <- lapply(barcodes,function(barcode) {
  fns <- fileNames[str_detect(fileNames, barcode)]
  dfl <- lapply(fns, function(fn){
    df <- read_tsv(fn) %>%
      as_tibble(.name_repair = "universal") %>%
      rename_all(.funs = str_remove_all, pattern = "[.]")
  })  
  
  #merge cyto and nuclei data with the same barcode
  df <- inner_join(dfl[[1]],dfl[[2]], by=c("ParentObjectIDMO" = "ObjectID")) %>%
    select(-contains("Parent"))
  df$barcode <- barcode
  return(df)
}) %>% bind_rows() %>%
  filter(Area>150) %>%
  rename(WellIndex = Well,
         MeanIntensity_DAPI = MeanIntensity377,
         TotalIntensity_DAPI = TotalIntensity377,
         MeanIntensity_KRT5 = MeanIntensity485,
         TotalIntensity_KRT5 = TotalIntensity485,
         MeanIntensity_CellMask = MeanIntensity560,
         TotalIntensity_CellMask = TotalIntensity560,
         MeanIntensity_EdU = MeanIntensity650,
         TotalIntensity_EdU = TotalIntensity650)

#Get the well level metadata
md <- lapply(barcodes, function(barcode){
  dt <- readMetadata(paste0("../IF/Metadata/",barcode,".xlsx"))
  dt$barcode <- barcode
  if(barcode %in% c(paste0("LI80200",4:9), paste0("LI8020",10:12))) dt$collection <- "C1"
  if(barcode %in% c(paste0("LI80230",1:9), paste0("LI8023",10:15))) dt$collection <- "C2"
  if(barcode %in% c(paste0("LI80200",4:6), paste0("LI80230",1:3))) dt$replicate <- "A"
  if(barcode %in% c(paste0("LI80200",7:9), paste0("LI80230",4:6))) dt$replicate <- "B"
  if(barcode %in% c(paste0("LI8020",10:12), paste0("LI80230",7:9))) dt$replicate <- "C"
  if(barcode %in% paste0("LI8023",10:11)) dt$replicate <- "D"
  if(barcode %in% paste0("LI8023",12:15)) dt$replicate <- "E"
  dt$Ligand[dt$IncubationTime==0] <- "ctrl"
  return(dt)
}) %>% bind_rows() %>%
  mutate(ligand = Ligand,
         experimentalTimePoint = IncubationTime,
         ligand=str_replace(ligand,"IFNg","IFNG"),
         ligand=str_replace(ligand,"TGFb","TGFB"),
         ligand=str_replace(ligand,"pbs","PBS"),
         specimenName=paste(ligand, experimentalTimePoint, collection, replicate, sep="_")) %>%
  select(specimenName, barcode, WellIndex)
#Merge in the MDD experimental metadata
mdd <- read_csv("../Metadata/MDD_sample_annotations.csv") %>%
  select(specimenName, specimenID) %>%
  right_join(md, by=c("specimenName"))

#Merge the data and metadata and gate some features
l1 <- mdd %>%
  inner_join(rd,by=c("barcode","WellIndex")) %>%
  group_by(barcode) %>%   #Gate at barcode level
  mutate(CellCycleState = gateOnlocalMinima(TotalIntensity_DAPI),
         EdUPositive = gateOnlocalMinima(MeanIntensity_EdU)) %>%
  ungroup()

#for (j in colnames(l1)) data.table::set(cell, j = j, value = shorten(l1[[j]]))
write_csv(l1, "../IF/Data/MDD_IF_Level1.csv")

#Summarize to the field level
parameterNames <- grep("Intensity|Elongation|Perimeter|Area|Proportion|Count|Thresh",colnames(l1),value = TRUE)

field <- l1 %>%
  group_by(barcode, WellIndex, Position, specimenName, specimenID) %>%
  mutate(FieldCellCount=n()) %>%
  mutate(DNA2nProportion = sum(CellCycleState==1)/n()) %>%
  mutate(EdUPositiveProportion = sum(EdUPositive==2)/n()) %>%
  summarise_at(c(parameterNames,"FieldCellCount","DNA2nProportion","EdUPositiveProportion","CellCycleState"),median) %>%
  mutate(Field=Position) %>%
  ungroup() %>%
  select(-Position)


#Summarize the fields to get a well level
parameterNames <- grep("Intensity|Elongation|Perimeter|Area|Proportion|Count|Thresh",colnames(field),value = TRUE)

well <- field %>%
  group_by(barcode, WellIndex, specimenID, specimenName) %>%
  mutate(WellCellCount = median(FieldCellCount)) %>%
  summarise_at(c(parameterNames,"FieldCellCount","DNA2nProportion","EdUPositiveProportion", "WellCellCount","CellCycleState"),median) %>%
  ungroup() %>%
  mutate(replicate = str_remove(specimenName, ".*_"),
         experimentalTimePoint = str_extract(specimenName,"0|24|48"),
         collection = str_extract(specimenName, "C[12]"))
write_csv(well, "../IF/Data/MDD_IF_Level3.csv")

#add EGF timecourse normalized values
well_EGFNorm <- well %>%
  filter(str_detect(specimenName, "EGF")) %>%
  select(specimenName, WellCellCount, DNA2nProportion, MeanIntensity_KRT5, EdUPositiveProportion, ElongationFactor, Perimeter, Area) %>%
  rename(WellCellCountEGF = WellCellCount,
         DNA2nProportionEGF = DNA2nProportion,
         MeanIntensity_KRT5EGF = MeanIntensity_KRT5,
         EdUPositiveProportionEGF = EdUPositiveProportion,
         ElongationFactorEGF = ElongationFactor,
         AreaEGF = Area,
         PerimeterEGF = Perimeter) %>%
  group_by(specimenName) %>%
  summarise_all(median) %>%
  ungroup() %>%
  mutate(replicate = str_remove(specimenName, ".*_"),
         experimentalTimePoint = str_extract(specimenName,"24|48"),
         collection = str_extract(specimenName, "C[12]")) %>%
  select(-specimenName) %>%
  left_join(well, by = c("replicate", "collection", "experimentalTimePoint")) %>%
  mutate(WellCellCount = WellCellCount/WellCellCountEGF,
         DNA2nProportion = DNA2nProportion/DNA2nProportionEGF,
         EdUPositiveProportion = EdUPositiveProportion/EdUPositiveProportionEGF,
         MeanIntensity_KRT5 = MeanIntensity_KRT5/MeanIntensity_KRT5EGF
         ,
         ElongationFactor = ElongationFactor/ElongationFactorEGF,
         Area = Area/AreaEGF,
         Perimeter = Perimeter/PerimeterEGF) %>%
  select(specimenID, WellCellCount, DNA2nProportion, EdUPositiveProportion, MeanIntensity_KRT5, ElongationFactor, Area, Perimeter) %>%
  group_by(specimenID) %>%
  summarise_all(median) %>%
  ungroup %>%
  gather(key = "feature", value = "value", -specimenID) %>%
  spread(key = specimenID, value = value)

write_csv(well_EGFNorm, "../IF/Data/MDD_IF_Level4.csv")
