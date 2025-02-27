########### downward facing hemispheric photos
source('hemispheric_PAI_functions.R')
library(terra)
library(hemispheR)
library(ggplot2)
library(patchwork)
library(ggpubr)
library(dplyr)
file_locations<-c("D:/Tern images")
im_deets<-read.csv(paste0(file_locations,'/Tern image details.csv'))
imlist<-im_deets[im_deets$Method=='dcp'&im_deets$Image.orientation.2 =='57_degrees',]
imlist<-split(imlist,imlist$file.path)

img_outputs<-"D:/Tern binarised images/DCP 57deg/"
canopy<-readRDS('D:/TERN_dcp_57_LAI.rds')

for(i in 1:length(imlist)){
  start<-Sys.time()
  dcp_files<-imlist[[i]]
  dcp<-dcp_files$file.path
  Code<-tools::file_path_sans_ext(basename(dcp))
  
  # read in img generate plot 1:
  rast_im<-read_B_and_gamma_correct(dcp,rgb=T)
  if(nlyr(rast_im)==3){
    rgbim<-rgb_plot(rast_im)
  }else{
    rgbim<-bw_plot(rast_im,title = '')
  }
  #plot(rast_im)
  
  # create a GLA index image
  #img.values <- -((-1 * terra::values(rast_im[[which(names(rast_im)=='Red')]]) + 2 * terra::values(rast_im[[which(names(rast_im)=='Green')]]) - 1 * terra::values(rast_im[[which(names(rast_im)=='Blue')]])) / (1 * terra::values(rast_im[[which(names(rast_im)=='Red')]]) + 2 * terra::values(rast_im[[which(names(rast_im)=='Green')]]) + 1 * terra::values(rast_im[[which(names(rast_im)=='Blue')]])))
  # create the target raster and assign the values
  rast_im<-rast_im[[which(names(rast_im)=='Blue')]]
  #range01 <- function(x){((x-min(x))/(max(x)-min(x)))*1000}
  #terra::values(rast_im)<-range01(img.values)
  
  # Binarise the image
  img.mat <- matrix(terra::values(rast_im),nrow = nrow(rast_im),byrow=T)
  if(is.na(dcp_files$thresh_QC)){
    th <- autothresholdr::auto_thresh(round(img.mat), method='Otsu',ignore_na=TRUE)
    th <- unlist(th)
  }else{
    th<-dcp_files$thresh_QC
  }
  img.bw <- terra::classify(rast_im, rbind(c(-Inf,th,0),c( th,Inf,1)))
  base::names(img.bw) <- base::names(rast_im)
  mk<- is.na(img.bw)
  img.bw[mk]<-NA
  # Make a histogram to inspect the classification
  histo<-hist_plot(rast_im,th[1])
  rm(rast_im)
  # Make a plot for inspection
  
  
  
  # Determine the gap fraction and gap size distribution
  canopy[[i]] <- cbind(as.data.frame(getGap(img.bw,path_id=Code,export.image=F,k=0.91)),dcp_files)  # change k to reflect image angle assuming spherical LAD
  imbw<-getGap(img.bw,path_id=Code,export.image=T)
  
  final_im<-rgbim+histo+imbw+plot_layout(ncol=2)
  ggsave(final_im,filename=paste0(img_outputs,Code,'.tiff'),width = 30,height=30,units='cm',dpi=72)
  end<-Sys.time()
  print(i)
  print(end-start)
  
  # save it every 50 pics just in case of crash
  if(i%%50==0){
    saveRDS(canopy,'D:/TERN_dcp_57_LAI.rds')
  }
}
saveRDS(canopy,'D:/TERN_dcp_57_LAI.rds')
write.csv(do.call(rbind,canopy),'D:/TERN_dcp_57_LAI.csv',row.names=F,na='')