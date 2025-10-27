########### downward facing hemispheric photos
source('hemispheric_PAI_functions.R')
library(terra)
library(hemispheR)
library(ggplot2)
library(patchwork)
library(ggpubr)
library(dplyr)
library(gtools)
library(autothresholdr)
file_locations<-c("D:/Tern images")
im_deets<-read.csv(paste0(file_locations,'/Tern image details.csv'))
imlist<-im_deets[im_deets$Method=='dcp'&im_deets$Image.orientation.2 =='downward',]
imlist<-split(imlist,imlist$file.path)

img_outputs<-"D:/Tern binarised images/DCP downward/"

canopy<-list()
for(i in 1:length(imlist)){
  start<-Sys.time()
  dcp_files<-imlist[[i]]
  dcp<-dcp_files$file.path
  Code<-tools::file_path_sans_ext(basename(dcp))
  
  # read in img generate plot 1:
  rast_im<-read_B_and_gamma_correct(dcp,rgb=T)
  #plot(rast_im)
  
  if(nlyr(rast_im)==3){
    rgbim<-rgb_plot(rast_im)
  }else{
    rgbim<-bw_plot(rast_im,title = '')
  }

  # create a GLA index image
  img.values <- -((-1 * terra::values(rast_im[[which(names(rast_im)=='Red')]]) + 2 * terra::values(rast_im[[which(names(rast_im)=='Green')]]) - 1 * terra::values(rast_im[[which(names(rast_im)=='Blue')]])) / (1 * terra::values(rast_im[[which(names(rast_im)=='Red')]]) + 2 * terra::values(rast_im[[which(names(rast_im)=='Green')]]) + 1 * terra::values(rast_im[[which(names(rast_im)=='Blue')]])))
  
  # create the target raster and assign the values
  rast_im<-rast_im[[which(names(rast_im)=='Blue')]]
  
# code to just stretch to min-max
  range01 <- function(x){
   y=na.omit(x)
   x=((x-min(y))/(max(y)-min(y)))
   return(x)
  }
  
  terra::values(rast_im)<-range01(img.values)
    
  rast_im<-terra::stretch(rast_im,minv=0,maxv=255,minq=0.01,maxq=0.99)

  # Binarise the image
  img.mat <- matrix(terra::values(rast_im),nrow = nrow(rast_im),byrow=T)
  th <- unlist(autothresholdr::auto_thresh(round(img.mat), method='Otsu',ignore_na=TRUE))[1]

  
  img.bw <- terra::classify(rast_im, rbind(c(-Inf,th,0),c( th,Inf,1)))
  base::names(img.bw) <- base::names(rast_im)
  mk<- is.na(img.bw)
  img.bw[mk]<-NA
  # Make a histogram to inspect the classification
  histo<-hist_plot(rast_im,th[1])
  rm(rast_im)
  # Make a plot for inspection
  
  # Determine the gap fraction and gap size distribution
  # alic looked like taken at 57 deg so applying this change to extinction coef (k)
  if(dcp_files$Site=='alic'){
    kidx=0.91
  }else{
    kidx=0.5
  }
  
  # Determine the gap fraction and gap size distribution
  canopy[[i]] <- cbind(as.data.frame(getGap(img.bw,path_id=Code,export.image=F,k=kidx)),dcp_files)
  
  # save final plot
  imbw<-getGap(img.bw,path_id=Code,export.image=T)
  final_im<-rgbim+histo+imbw+plot_layout(ncol=2)
  ggsave(final_im,filename=paste0(img_outputs,Code,'.tiff'),width = 30,height=30,units='cm',dpi=72)
  
  end<-Sys.time()
  print(i)
  print(end-start)
  

  # save it every 50 pics just in case of crash
  if(i%%50==0){
    saveRDS(canopy,'D:/TERN_dcp_down_LAI.rds')
  }
}

saveRDS(canopy,'D:/TERN_dcp_down_LAI.rds')
write.csv(do.call(smartbind,canopy),'D:/TERN_dcp_down_LAI_final.csv',row.names=F,na='')
