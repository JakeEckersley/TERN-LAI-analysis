source('hemispheric_PAI_functions.R')
library(terra)
library(hemispheR)
library(ggplot2)
library(patchwork)
library(ggpubr)
library(dplyr)
file_locations<-c("D:/Tern images")
im_deets<-read.csv(paste0(file_locations,'/Tern image details.csv'))
imlist<-im_deets[im_deets$Method=='dcp'&im_deets$Image.orientation.2 =='upward',]
imlist<-split(imlist,imlist$file.path)


img_outputs<-"D:/Tern binarised images/DCP upward/"
canopy<-list()

############# check whether th is working ok here.
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
    rgbim<-bw_plot(rast_im,title = Code)
  }
  #plot(rast_im)
  
  # select the blue band
  rast_im<-rast_im[[which(names(rast_im)=='Blue')]]
  
  # Binarise the image
  img.mat <- matrix(terra::values(rast_im),nrow = nrow(rast_im),byrow=T)
  th <- autothresholdr::auto_thresh(round(img.mat), method='Otsu',ignore_na=TRUE)
  th <- unlist(th)[1]
  
  img.bw <- terra::classify(rast_im, rbind(c(-Inf,th,0),c( th,Inf,1)))
  base::names(img.bw) <- base::names(rast_im)
  mk<- is.na(img.bw)
  img.bw[mk]<-NA
  
  # Make a histogram to inspect the classification
  histo<-hist_plot(rast_im,th[1])
  rm(rast_im)
  
  
  # Make a plot for inspection
  canopy[[i]] <- cbind(as.data.frame(getGap(img.bw,path_id=Code,export.image=F)),dcp_files)

  imbw<-getGap(img.bw,path_id=Code,export.image=T)
  final_im<-rgbim+histo+imbw+plot_layout(ncol=2)
  ggsave(final_im,filename=paste0(img_outputs,Code,'.tiff'),width = 30,height=30,units='cm',dpi=72)
  
  end<-Sys.time()
  print(i)
  print(end-start)
  
  # save it every 50 pics just in case of crash
  if(i%%50==0){
    saveRDS(canopy,'D:/TERN_dcp_LAI.rds')
  }
}
saveRDS(canopy,'D:/TERN_dcp_LAI.rds')
write.csv(do.call(rbind,canopy),'D:/TERN_dcp_LAI.csv',row.names=F,na='')
