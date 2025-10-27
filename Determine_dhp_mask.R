#Determine mask for dhp
source('hemispheric_PAI_functions.R')
library(terra)
library(hemispheR)
library(ggplot2)
library(patchwork)
library(ggpubr)


file_locations<-c("D:/Tern images")
im_details<-read.csv(paste0(file_locations,'/Tern image details.csv'))
dhp_files<-im_details[im_details$Method=='dhp',]
dhp_folders<-unique(dhp_files$file.dirs)

mask_list<-readRDS('DHP_mask_list.rds')

i=1

#choose a directory to examine:
focus_dir<-dhp_files[which(dhp_files$file.dirs==names(mask_list)[i]),]
im=1
#im<-which(focus_dir$Image.code=='DSC_6155')

dhp<-focus_dir$file.path[im]
dhp_rast<-read_B_and_gamma_correct(dhp,rgb=F)
plot(dhp_rast)

#inspect mask
hem_mask<-mask_list[[i]]
plot_with_circ(dhp_rast,hem_mask)
plot(dhp_rast)


######## code to determine the crop extent of DHCP
width_height<-c(ncol(dhp_rast),nrow(dhp_rast))
long_side<-which(width_height==max(width_height))
short_side<-which(width_height==min(width_height))

portrait<-NA
if(length(short_side)==1){
  crop_start<-(width_height[long_side]-width_height[short_side])/2
  if(long_side==1){
    xmin=crop_start
    xmax=width_height[1]-crop_start
    ymin=0
    ymax=width_height[2]
    portrait=2
  }else{
    xmin=0
    xmax=width_height[1]
    ymin=crop_start
    ymax=width_height[2]-crop_start
    portrait=1
  }
}else{
  xmin=0
  xmax=width_height[1]
  ymin=0
  ymax=width_height[2]
  portrait=0
}

crop_ext<-ext(xmin,xmax,ymin,ymax)

#add the mask
hem_mask<-find_mask(crop_ext,portrait)
plot_with_circ(dhp_rast,hem_mask)

# if it is a full frame image, unhash and use this instead!!
#hem_mask<-find_mask_FF(dhp_rast,portrait)
#plot_with_circ(dhp_rast,hem_mask,circ=F)

mask_list[[i]]<-hem_mask
#names(mask_list)[i]<-dhp_folders[i] #<- check this line works as intended. no new files can be added for this to function correctly.
saveRDS(mask_list,'DHP_mask_list.rds')
