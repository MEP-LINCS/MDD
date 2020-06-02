#import libraries
library(readr)
library(dplyr)
library(ggplot2)
library(tidyr)
library(stringr)
library(LPCM)


#function for processing image data and assigning metadata
image_process<-function(image_name, ctrl=F, time, collection){
  image_dat<-read.csv(image_name, stringsAsFactors = F)
  cell_count_dat<-image_dat %>%
    select(Count_Nuclei, ImageNumber, Metadata_Plate, Metadata_Run, Metadata_Well)
  if(ctrl==T){
    cell_count_dat$ligand<-'ctrl'
    cell_count_dat$time<-0
  } else{
    cell_count_dat$ligand[cell_count_dat$Metadata_Well=='A1' | cell_count_dat$Metadata_Well=='B1'] = 'PBS'
    cell_count_dat$ligand[cell_count_dat$Metadata_Well=='A2'] = 'EGF' 
    cell_count_dat$ligand[cell_count_dat$Metadata_Well=='A3']="HGF"
    cell_count_dat$ligand[cell_count_dat$Metadata_Well=='A4']='OSM'
    cell_count_dat$ligand[cell_count_dat$Metadata_Well=='B2']='BMP2'
    cell_count_dat$ligand[cell_count_dat$Metadata_Well=='B3']='IFNG'
    cell_count_dat$ligand[cell_count_dat$Metadata_Well=='B4']="TGFB"
    cell_count_dat$time<-time
  }
  
  cell_count_tot<-cell_count_dat
  cell_count_tot$size<-(1344*1024)
  cell_count_tot$null<-.5/sqrt(cell_count_tot$Count_Nuclei/cell_count_tot$size)
  cell_count_tot$collection<-collection
  
  
  return(cell_count_tot)
}

#function for calculating distance from cell to all other cells in a set
Euclidean<-function(data, point1){
  distances<-numeric()
  #remove seed1 from data frame
  data_comp<-anti_join(data, point1, by = c("AreaShape_Center_X", "AreaShape_Center_Y"))
  #calculate and return distances
  for(i in 1:nrow(data_comp)){
    distances[i]<-dist(rbind(data_comp[i,],point1))
  }
  return(distances)
}

#function to find nearest neighbor distances
nearest_neighbor_dist<-function(nuc_data){
  coord_data<-nuc_data[c(7:8)]
  ln<-nrow(coord_data)
  top_distances<-data.frame(firstneighbor=numeric(),secondneighbor=numeric(),thirdneighbor=numeric(),fourthneighbor=numeric())
  
  #search for nearest neighbor within 400 pixel square
  for(i in 1:nrow(coord_data)){
    point<-coord_data[i,]
    high_bounds_x<-as.numeric(point[1] + 200)
    low_bounds_x<-as.numeric(point[1] - 200)
    high_bounds_y<-as.numeric(point[2] + 200)
    low_bounds_y<-as.numeric(point[2] - 200)
    
    search_data<-subset(coord_data, (AreaShape_Center_X > low_bounds_x) & AreaShape_Center_X < high_bounds_x & 
                          (AreaShape_Center_Y > low_bounds_y & AreaShape_Center_Y < high_bounds_y))
    
    #if four neighbors not found, expand square
    if(nrow(search_data)<4){
      high_bounds_x<-as.numeric(point[1] + 400)
      low_bounds_x<-as.numeric(point[1] - 400)
      high_bounds_y<-as.numeric(point[2] + 400)
      low_bounds_y<-as.numeric(point[2] - 400)
      
      search_data<-subset(coord_data, (AreaShape_Center_X > low_bounds_x) & AreaShape_Center_X < high_bounds_x & 
                            (AreaShape_Center_Y > low_bounds_y & AreaShape_Center_Y < high_bounds_y))
      
      if(nrow(search_data)<4){
        high_bounds_x<-as.numeric(point[1] + 600)
        low_bounds_x<-as.numeric(point[1] - 600)
        high_bounds_y<-as.numeric(point[2] + 600)
        low_bounds_y<-as.numeric(point[2] - 600)
        
        search_data<-subset(coord_data, (AreaShape_Center_X > low_bounds_x) & AreaShape_Center_X < high_bounds_x & 
                              (AreaShape_Center_Y > low_bounds_y & AreaShape_Center_Y < high_bounds_y))
        if(nrow(search_data)<4){
          
          search_data<-coord_data
        }}}
    
    #find distance to nearest neighbors
    d<-Euclidean(search_data,point1=point)
    sorted_d<-sort(d)
    
    #remove cells on edges of image
    if(point[1,1]<100 | point[1,2]<100 
       | (1344-point[1,1])<100 | (1024-point[1,2] < 100)){
      top_distances[i,1]<-NA
      top_distances[i,2]<-NA
      top_distances[i,3]<-NA
      top_distances[i,4]<-NA
      
    } else{
      top<-sorted_d[1:4]
      
      #if no neighbor found, set to maximum distance
      top[is.na(top)]<-1344
      top_distances[i,1:4]<-top[1:4]
    }
  }
  
  return(top_distances)
}

#function to calculate number of neighbors in a given boundary
Number_of_Neighbors<-function(nuc_data, distance_bound){
  coord_data<-nuc_data[c(7:8)]
  
  ln<-nrow(coord_data)
  neighbor_counts<-data.frame(neighbor_n=numeric(),boundary_distance=numeric())
  
  #define boundary box for neighbor search
  for(i in 1:nrow(coord_data)){
    point<-coord_data[i,]
    high_bounds_x<-as.numeric(point[1] + distance_bound)
    low_bounds_x<-as.numeric(point[1] - distance_bound)
    high_bounds_y<-as.numeric(point[2] + distance_bound)
    low_bounds_y<-as.numeric(point[2] - distance_bound)
    
    #find the cells within the neighborhood
    search_data<-subset(coord_data, (AreaShape_Center_X > low_bounds_x) & AreaShape_Center_X < high_bounds_x & 
                          (AreaShape_Center_Y > low_bounds_y & AreaShape_Center_Y < high_bounds_y))
    d<-Euclidean(search_data,point1=point)
    
    #first check if cell is on the edge, if so remove it
    if(point[1,1]<distance_bound | point[1,2]<distance_bound 
       | (1344-point[1,1])<distance_bound | (1024-point[1,2] < distance_bound)){
      neighbor_counts[i,1]<-NA
      neighbor_counts[i,2]<-NA
    } 
    else{
      #if no neighbors found, set to 0
      if(anyNA(d)){
        neighbor_counts[i,1]<-0
        neighbor_counts[i,2]<-distance_bound
      }else{
        cells_in_bounds<-length(d[d<distance_bound])
        neighbor_counts[i,1]<-cells_in_bounds
        neighbor_counts[i,2]<-distance_bound
      }
    }
  }
  return(neighbor_counts)
}

#function to use mean-shift clustering to calculate cluster sizes
cluster_sizes<-function(cell_data, h){
  cluster<-data.frame()
  collections<-unique(cell_data$collection)
  
  for(i in 1:length(collections)){
    current_coll<-collections[i]
    current_coll_sub<-cell_data[cell_data$collection==current_coll,]
    run_list<-unique(current_coll_sub$Metadata_Run)
    
    for(j in 1:length(run_list)){
      current_run<-run_list[j]
      current_run_sub<-current_coll_sub[current_coll_sub$Metadata_Run==current_run,]
      image_list<-unique(current_run_sub$ImageNumber)
      
      for(k in 1:length(image_list)){
        current_image<-image_list[k]
        current_image_sub<-current_run_sub[current_run_sub$ImageNumber==current_image,]
        ligand<-current_image_sub$ligand[1]
        time<-current_image_sub$time[1]
        cell_count<-nrow(current_image_sub)
        
        #for each image, check if more then one cell
        if(cell_count>1){
          coord_data<-as.matrix(current_image_sub[c(7:8)])
          
          #perform mean-shift clustering
          ms<-ms(coord_data, scaled=0, plot=F, h=h)
          
          labels<-as.data.frame(ms$cluster.label) %>%
            rename(Cluster_Label=`ms$cluster.label`)
          
          #record the cluster labels
          label_counts<-labels %>%
            group_by(Cluster_Label) %>%
            summarize(n()) %>%
            rename(label_count='n()') %>%
            mutate(collection=current_coll, Metadata_Run=current_run, Image_Cell_Count=cell_count,
                   ImageNumber=current_image, time=time, ligand=ligand) 
        }else{
          labels<-as.data.frame(T)
          label_counts<-labels %>%
            mutate(collection=current_coll, Metadata_Run=current_run, Image_Cell_Count=cell_count,
                   ImageNumber=current_image, time=time, ligand=ligand) %>%
            rename(SingleCell='T')
          
        }
        cluster<-bind_rows(cluster, label_counts)
      }
    }
  }
  return(cluster)
}

#read in metadata
sample_meta<-read_csv('../Metadata/MDD_sample_annotations.csv')

#names of files for nuclear, cytoplasmic, and image csvs
nuc_name=c('../Immunofluorescence/Level0/LI802004_Nuclei.csv', '../Immunofluorescence/Level0/LI802005_Nuclei.csv','../Immunofluorescence/Level0/LI802006_Nuclei.csv',
           '../Immunofluorescence/Level0/LI802007_Nuclei.csv','../Immunofluorescence/Level0/LI802008_Nuclei.csv','../Immunofluorescence/Level0/LI802009_Nuclei.csv',
           '../Immunofluorescence/Level0/LI802010_Nuclei.csv','../Immunofluorescence/Level0/LI802011_Nuclei.csv','../Immunofluorescence/Level0/LI802012_Nuclei.csv',
           '../Immunofluorescence/Level0/LI802301_Nuclei.csv','../Immunofluorescence/Level0/LI802302_Nuclei.csv','../Immunofluorescence/Level0/LI802303_Nuclei.csv',
           '../Immunofluorescence/Level0/LI802304_Nuclei.csv','../Immunofluorescence/Level0/LI802305_Nuclei.csv','../Immunofluorescence/Level0/LI802306_Nuclei.csv',
           '../Immunofluorescence/Level0/LI802307_Nuclei.csv','../Immunofluorescence/Level0/LI802308_Nuclei.csv','../Immunofluorescence/Level0/LI802309_Nuclei.csv',
           '../Immunofluorescence/Level0/LI802310_Nuclei.csv','../Immunofluorescence/Level0/LI802311_Nuclei.csv','../Immunofluorescence/Level0/LI802312_Nuclei.csv',
           '../Immunofluorescence/Level0/LI802313_Nuclei.csv','../Immunofluorescence/Level0/LI802314_Nuclei.csv','../Immunofluorescence/Level0/LI802315_Nuclei.csv')
cyto_name = c('../Immunofluorescence/Level0/LI802004_Cytoplasm.csv', '../Immunofluorescence/Level0/LI802005_Cytoplasm.csv','../Immunofluorescence/Level0/LI802006_Cytoplasm.csv',
              '../Immunofluorescence/Level0/LI802007_Cytoplasm.csv','../Immunofluorescence/Level0/LI802008_Cytoplasm.csv','../Immunofluorescence/Level0/LI802009_Cytoplasm.csv',
              '../Immunofluorescence/Level0/LI802010_Cytoplasm.csv','../Immunofluorescence/Level0/LI802011_Cytoplasm.csv','../Immunofluorescence/Level0/LI802012_Cytoplasm.csv',
              '../Immunofluorescence/Level0/LI802301_Cytoplasm.csv','../Immunofluorescence/Level0/LI802302_Cytoplasm.csv','../Immunofluorescence/Level0/LI802303_Cytoplasm.csv',
              '../Immunofluorescence/Level0/LI802304_Cytoplasm.csv','../Immunofluorescence/Level0/LI802305_Cytoplasm.csv','../Immunofluorescence/Level0/LI802306_Cytoplasm.csv',
              '../Immunofluorescence/Level0/LI802307_Cytoplasm.csv','../Immunofluorescence/Level0/LI802308_Cytoplasm.csv','../Immunofluorescence/Level0/LI802309_Cytoplasm.csv',
              '../Immunofluorescence/Level0/LI802310_Cytoplasm.csv','../Immunofluorescence/Level0/LI802311_Cytoplasm.csv','../Immunofluorescence/Level0/LI802312_Cytoplasm.csv',
              '../Immunofluorescence/Level0/LI802313_Cytoplasm.csv','../Immunofluorescence/Level0/LI802314_Cytoplasm.csv','../Immunofluorescence/Level0/LI802315_Cytoplasm.csv')
image_name = c('../Immunofluorescence/Level0/LI802004_Image.csv', '../Immunofluorescence/Level0/LI802005_Image.csv','../Immunofluorescence/Level0/LI802006_Image.csv',
               '../Immunofluorescence/Level0/LI802007_Image.csv','../Immunofluorescence/Level0/LI802008_Image.csv','../Immunofluorescence/Level0/LI802009_Image.csv',
               '../Immunofluorescence/Level0/LI802010_Image.csv','../Immunofluorescence/Level0/LI802011_Image.csv','../Immunofluorescence/Level0/LI802012_Image.csv',
               '../Immunofluorescence/Level0/LI802301_Image.csv','../Immunofluorescence/Level0/LI802302_Image.csv','../Immunofluorescence/Level0/LI802303_Image.csv',
               '../Immunofluorescence/Level0/LI802304_Image.csv','../Immunofluorescence/Level0/LI802305_Image.csv','../Immunofluorescence/Level0/LI802306_Image.csv',
               '../Immunofluorescence/Level0/LI802307_Image.csv','../Immunofluorescence/Level0/LI802308_Image.csv','../Immunofluorescence/Level0/LI802309_Image.csv',
               '../Immunofluorescence/Level0/LI802310_Image.csv','../Immunofluorescence/Level0/LI802311_Image.csv','../Immunofluorescence/Level0/LI802312_Image.csv',
               '../Immunofluorescence/Level0/LI802313_Image.csv','../Immunofluorescence/Level0/LI802314_Image.csv','../Immunofluorescence/Level0/LI802315_Image.csv')

#read in and process all IF files
total_cell_data<-data.frame()
total_image_data<-data.frame()

#Filter files and assign metadata
for(i in 1:length(nuc_name)){
  n_name<-nuc_name[i]
  c_name<-cyto_name[i]
  i_name<-image_name[i]
  
  plate<-gsub('.*LI80','',n_name)
  run<-as.numeric(str_extract(plate, '([0-9]{4})'))
  
  nucl_data<-read_csv(n_name, col_types = cols())
  cyto_data<-read_csv(c_name, col_types = cols())
  
  nuclei_data<-nucl_data %>%
    select(-contains('Metadata_Channel'),-Metadata_FileLocation,  -Metadata_Frame, -Metadata_Series,
           -Children_Cytoplasm_Count, -Number_Object_Number,  -contains('_Z'), -AreaShape_Center_Z, 
           -Location_Center_X, -Location_Center_Y, -contains('Edge'), -contains('MAD'), -contains('MassDisp'))
  
  cytoplasm_data<-cyto_data %>%
    select(-contains('Metadata_Channel'),-Metadata_FileLocation,  -Metadata_Frame, -Metadata_Series,
           -contains('Location'), -Number_Object_Number, -Parent_Nuclei)
  
  total_data<-full_join(nuclei_data,cytoplasm_data, by=c("ImageNumber", "ObjectNumber",
                                                         "Metadata_Plate", "Metadata_Well", "Metadata_Run"))
  
  #scale intensities properly by 65535
  total_data_int_rm<- total_data %>%
    select(-contains('Intensity'))
  total_data_int <- total_data %>%
    select(contains('Intensity'))
  total_data_adjusted_int<-total_data_int*65535
  total_data<-bind_cols(total_data_int_rm,total_data_adjusted_int)
  
  #set metadata
  if(run==2004 | run==2007 | run==2010 | run==2301 | run==2304 | run==2307 | run==2310 | run==2313){
    ctrl=T
    total_data$ligand<-'ctrl'
    total_data$time<-0
    
  }else if(run==2005 | run==2008 | run==2011 | run==2302 | run==2305 | run==2308 | run==2311 | run==2312){
    ctrl=F
    time=24
    total_data$time<-24
    total_data$ligand[total_data$Metadata_Well == 'A1'] <- 'PBS'
    total_data$ligand[total_data$Metadata_Well == 'A2'] <- 'EGF'
    total_data$ligand[total_data$Metadata_Well == 'A3'] <- 'HGF'
    total_data$ligand[total_data$Metadata_Well == 'A4'] <- 'OSM'
    total_data$ligand[total_data$Metadata_Well == 'B1'] <- 'PBS'
    total_data$ligand[total_data$Metadata_Well == 'B2'] <- 'BMP2'
    total_data$ligand[total_data$Metadata_Well == 'B3'] <- 'IFNG'
    total_data$ligand[total_data$Metadata_Well == 'B4'] <- 'TGFB'
  } else{
    ctrl=F
    time=48
    total_data$time<-48
    total_data$ligand[total_data$Metadata_Well == 'A1'] <- 'PBS'
    total_data$ligand[total_data$Metadata_Well == 'A2'] <- 'EGF'
    total_data$ligand[total_data$Metadata_Well == 'A3'] <- 'HGF'
    total_data$ligand[total_data$Metadata_Well == 'A4'] <- 'OSM'
    total_data$ligand[total_data$Metadata_Well == 'B1'] <- 'PBS'
    total_data$ligand[total_data$Metadata_Well == 'B2'] <- 'BMP2'
    total_data$ligand[total_data$Metadata_Well == 'B3'] <- 'IFNG'
    total_data$ligand[total_data$Metadata_Well == 'B4'] <- 'TGFB'
  }
  
  if(run<2013){
    collection='C1'
    total_data$collection<-collection
  }else{
    collection='C2'
    total_data$collection<-collection
  }
  
  image_data<-image_process(i_name, ctrl=ctrl, time=time, collection=collection)
  
  total_image_data<-bind_rows(image_data, total_image_data)
  total_cell_data<-bind_rows(total_data, total_cell_data)
}

#calculate image means for intensity / area measurements
image_int_means<-total_cell_data %>%
  select(-ObjectNumber, -AreaShape_Center_X, -AreaShape_Center_Y, -contains('Location_')) %>%
  group_by(ImageNumber, Metadata_Plate, Metadata_Run, Metadata_Well, time,
           collection, ligand) %>%
  summarize_if(is.numeric, mean)

###find distances to nearest neighbors
total_cell_distances<-total_cell_data %>%
  group_by(ImageNumber, Metadata_Plate, Metadata_Run, Metadata_Well, time,
           collection, ligand) %>%
  do(data.frame(., e=nearest_neighbor_dist(.)))

#bind to cell data
IF_cells_distance<-(full_join(total_cell_data, total_cell_distances)) %>%
  rename(FirstNeighbor_Dist=e.firstneighbor,SecondNeighbor_Dist=e.secondneighbor,
         ThirdNeighbor_Dist=e.thirdneighbor,FourthNeighbor_Dist=e.fourthneighbor) 

#calculate image means for all measures
distance_images_total<-IF_cells_distance %>%
  select(-ObjectNumber, -AreaShape_Center_X, -AreaShape_Center_Y, -contains('Location_')) %>%
  group_by(ImageNumber, Metadata_Plate, Metadata_Run, Metadata_Well, time,
           collection, ligand) %>%
  summarize_if(is.numeric, funs(mean, .args = list(na.rm=T)))

#join to image data
distance_images_bind<-full_join(distance_images_total, total_image_data, by = c("ImageNumber", "Metadata_Plate", "Metadata_Run", "Metadata_Well", "time", "collection", "ligand")) 

#calculate normalized distance metrics
distance_images_bind$Normalized_First_Neighbor_Dist<-distance_images_bind$FirstNeighbor_Dist/distance_images_bind$null
distance_images_bind$Normalized_Second_Neighbor_Dist<-distance_images_bind$SecondNeighbor_Dist/distance_images_bind$null
distance_images_bind$Normalized_Third_Neighbor_Dist<-distance_images_bind$ThirdNeighbor_Dist/distance_images_bind$null
distance_images_bind$Normalized_Fourth_Neighbor_Dist<-distance_images_bind$FourthNeighbor_Dist/distance_images_bind$null

###find nearest neighbors
total_neighborhood_numbers<-total_cell_data %>%
  group_by(ImageNumber, Metadata_Plate, Metadata_Run, Metadata_Well, time,
           collection, ligand) %>%
  do(data.frame(., e=Number_of_Neighbors(., distance_bound=100)))

#rename neighbor measure
IF_cells_neighbors<-(full_join(total_cell_data, total_neighborhood_numbers)) %>%
  rename(number_neighbors=e.neighbor_n)

#calculate image means for all measures
neighborhood_images_total<-IF_cells_neighbors %>%
  select(-ObjectNumber, -AreaShape_Center_X, -AreaShape_Center_Y, contains('Location_')) %>%
  group_by(ImageNumber, Metadata_Plate, Metadata_Run, Metadata_Well, time,
           collection, ligand) %>%
  summarize_if(is.numeric, funs(mean, .args = list(na.rm=T)))

#join to image data, 
neighborhood_images_bind<-full_join(neighborhood_images_total, total_image_data, by = c("ImageNumber", "Metadata_Plate", "Metadata_Run", "Metadata_Well", "time", "collection", "ligand")) 


###mean shift clustering analysis
total_cluster<-cluster_sizes(cell_data=total_cell_data, h=30)

#summarize to image
image_cluster_summary<-full_join(total_cluster, total_image_data, Joining, by = c("collection", "Metadata_Run", "ImageNumber", "time", "ligand")) %>%
  group_by(collection, Metadata_Run, Metadata_Well, Metadata_Plate, ImageNumber, time, ligand) %>%
  summarize_if(is.numeric, funs(mean, .args = list(na.rm=T)))

#calculate proportion of cells in cluster > n=8
image_cluster_proportions<-total_cluster %>%
  group_by(collection, Metadata_Run, ImageNumber, time, ligand) %>%
  filter(label_count>8) %>%
  mutate(cells_in_large_cluster=sum(label_count)) 

cluster_bind<-full_join(image_cluster_proportions, total_cluster)
cluster_bind[is.na(cluster_bind$cells_in_large_cluster),]$cells_in_large_cluster<-0

#calculate image means for all measures
image_props<-cluster_bind %>%
  group_by(collection, Metadata_Run, ImageNumber, time, ligand) %>%
  mutate(proportion_in_large_cluster=max(cells_in_large_cluster)/mean(Image_Cell_Count)) %>%
  summarise(proportion_in_large_cluster=mean(proportion_in_large_cluster, na.rm = T))


#bind all cell data together
total_distance_data<-full_join(IF_cells_distance, IF_cells_neighbors) 

###rename to conform with SID metadata
total_distance_data$Metadata_Well[total_distance_data$Metadata_Well=='A1'] = 1
total_distance_data$Metadata_Well[total_distance_data$Metadata_Well=='A2'] = 2
total_distance_data$Metadata_Well[total_distance_data$Metadata_Well=='A3']=3
total_distance_data$Metadata_Well[total_distance_data$Metadata_Well=='A4']=4
total_distance_data$Metadata_Well[total_distance_data$Metadata_Well=='B1']=5
total_distance_data$Metadata_Well[total_distance_data$Metadata_Well=='B2']=6
total_distance_data$Metadata_Well[total_distance_data$Metadata_Well=='B3']=7
total_distance_data$Metadata_Well[total_distance_data$Metadata_Well=='B4']=8

total_distance_data$replicate[total_distance_data$Metadata_Run>2003 & total_distance_data$Metadata_Run<2007]='A'
total_distance_data$replicate[total_distance_data$Metadata_Run>2006 & total_distance_data$Metadata_Run<2010]='B'
total_distance_data$replicate[total_distance_data$Metadata_Run>2009 & total_distance_data$Metadata_Run<2013]='C'
total_distance_data$replicate[total_distance_data$Metadata_Run>2300 & total_distance_data$Metadata_Run<2304]='A'
total_distance_data$replicate[total_distance_data$Metadata_Run>2303 & total_distance_data$Metadata_Run<2307]='B'
total_distance_data$replicate[total_distance_data$Metadata_Run>2306 & total_distance_data$Metadata_Run<2310]='C'
total_distance_data$replicate[total_distance_data$Metadata_Run>2309 & total_distance_data$Metadata_Run<2313]='D'
total_distance_data$replicate[total_distance_data$Metadata_Run>2312 & total_distance_data$Metadata_Run<2316]='C'

###prepare final cell file
final_cell_file<-total_distance_data %>%
  ungroup() %>%
  rename(WellIndex=Metadata_Well) %>%
  mutate(specimenName=paste0(ligand,'_',time,'_',collection,'_',replicate)) %>%
  mutate(Metadata_Run=paste0('LI80',Metadata_Run)) %>%
  rename(barcode=Metadata_Run, Number_Neighbors=number_neighbors, experimentalTimePoint=time) %>%
  select(-Metadata_Plate) %>%
  mutate_if(is.numeric, funs(signif, .args = list(digits=4)))

#set specimen ID based on SID metadata
j<-as.vector(sample_meta[match(final_cell_file$specimenName, sample_meta$specimenName),1])
final_cell_file$specimenID<-j$specimenID

#filter cell file for desired columns
final_cell_file_filtered<-final_cell_file[,c(87,88,1,3,4,2,77:79,86,23:38,47:76,5:9,11:22,39:46,80:84)]

#output cell csv
write_csv(final_cell_file_filtered, '../Immunofluorescence/MDD_IF_Cell_File.csv')

###prep image file
#bind all image data together
total_density_image_data<-full_join(distance_images_bind[,c(1:7,70:80)], neighborhood_images_bind[,c(1:7,78:82)]) %>%
  full_join(image_cluster_summary) %>%
  full_join(image_props) %>%
  full_join(image_int_means)

#prep image file
total_density_image_data$Metadata_Well[total_density_image_data$Metadata_Well=='A1'] = 1
total_density_image_data$Metadata_Well[total_density_image_data$Metadata_Well=='A2'] = 2
total_density_image_data$Metadata_Well[total_density_image_data$Metadata_Well=='A3']=3
total_density_image_data$Metadata_Well[total_density_image_data$Metadata_Well=='A4']=4
total_density_image_data$Metadata_Well[total_density_image_data$Metadata_Well=='B1']=5
total_density_image_data$Metadata_Well[total_density_image_data$Metadata_Well=='B2']=6
total_density_image_data$Metadata_Well[total_density_image_data$Metadata_Well=='B3']=7
total_density_image_data$Metadata_Well[total_density_image_data$Metadata_Well=='B4']=8

total_density_image_data$replicate[total_density_image_data$Metadata_Run>2003 & total_density_image_data$Metadata_Run<2007]='A'
total_density_image_data$replicate[total_density_image_data$Metadata_Run>2006 & total_density_image_data$Metadata_Run<2010]='B'
total_density_image_data$replicate[total_density_image_data$Metadata_Run>2009 & total_density_image_data$Metadata_Run<2013]='C'
total_density_image_data$replicate[total_density_image_data$Metadata_Run>2300 & total_density_image_data$Metadata_Run<2304]='A'
total_density_image_data$replicate[total_density_image_data$Metadata_Run>2303 & total_density_image_data$Metadata_Run<2307]='B'
total_density_image_data$replicate[total_density_image_data$Metadata_Run>2306 & total_density_image_data$Metadata_Run<2310]='C'
total_density_image_data$replicate[total_density_image_data$Metadata_Run>2309 & total_density_image_data$Metadata_Run<2313]='D'
total_density_image_data$replicate[total_density_image_data$Metadata_Run>2312 & total_density_image_data$Metadata_Run<2316]='C'

final_image_file<-total_density_image_data %>%
  ungroup() %>%
  rename(WellIndex=Metadata_Well) %>%
  mutate(specimenName=paste0(ligand,'_',time,'_',collection,'_',replicate)) %>%
  mutate(Metadata_Run=paste0('LI80',Metadata_Run)) %>%
  rename(barcode=Metadata_Run, Mean_Cells_per_Cluster=label_count,Number_Neighbors=number_neighbors,
         Expected_Random_Neighbor_Dist=null, experimentalTimePoint=time, Proportion_in_large_cluster=proportion_in_large_cluster) %>%
  select(-Metadata_Plate, -contains('Location'),-size) %>%
  mutate_if(is.numeric, funs(signif, .args = list(digits=4)))

final_image_file$specimenID<-NA
j<-as.vector(sample_meta[match(final_image_file$specimenName, sample_meta$specimenName),1])
final_image_file$specimenID<-j$specimenID

final_image_file_filtered<-final_image_file[,c(86,87,2,3,1,4:6,85,21,39:84,23:25,27:38,13:17,20,22)]

write_csv(final_image_file_filtered, '../Immunofluorescence/MDD_IF_Image_File.csv')

###prep well file
final_well_file<-final_image_file_filtered %>%
  ungroup() %>%
  select(-ImageNumber) %>%
  group_by(specimenName, specimenID, barcode, WellIndex, experimentalTimePoint, collection, ligand, replicate) %>%
  mutate(Well_Cell_Count=sum(Image_Cell_Count)) %>%
  summarize_if(is.numeric, funs(mean, .args = list(na.rm=T))) %>%
  mutate_if(is.numeric, funs(signif, .args = list(digits=4)))

final_well_file_filtered<-final_well_file[c(1:8,78, 9:77)]

write_csv(final_well_file_filtered, '../Immunofluorescence/MDD_IF_Well_File.csv')

###prep specimenID file
SID_file<-final_well_file_filtered %>%
  ungroup() %>%
  select(-WellIndex, -barcode) %>%
  group_by(specimenName, specimenID,experimentalTimePoint, collection, ligand, replicate) %>%
  summarize_if(is.numeric, funs(mean, .args = list(na.rm=T))) %>%
  mutate_if(is.numeric, funs(signif, .args = list(digits=4)))

write_csv(SID_file, '../Immunofluorescence/MDD_IF_SID_File.csv')

###prep collection file
final_collection_file<-final_well_file_filtered %>%
  ungroup() %>%
  select(-WellIndex, -specimenName, -specimenID, -barcode, -replicate) %>%
  group_by(ligand, experimentalTimePoint, collection) %>%
  summarize_if(is.numeric, funs(mean, .args = list(na.rm=T))) %>%
  mutate(conditionName=paste0(ligand,'_',experimentalTimePoint,'_',collection)) %>%
  mutate_if(is.numeric, funs(signif, .args = list(digits=4)))

write_csv(final_collection_file, '../Immunofluorescence/MDD_IF_Collection_File.csv')

