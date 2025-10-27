########### downward facing dcp photos
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
#imlist<-split(imlist,imlist$file.path)
QC<-read.csv("C:/Users/22064705/OneDrive - UWA/My Documents/TERN LAI/250803_final_check/TERN_dcp_down_LAI_QC.csv")

img_outputs<-"D:/Tern binarised images/DCP downward/"
#canopy<-readRDS("D:/TERN_dcp_down_LAI.rds")

col_idx<-which(colnames(QC)=='Le')
if(length(col_idx)==1){colnames(QC)[col_idx]<-"PAIe"}
col_idx<-which(colnames(QC)=='L')
if(length(col_idx)==1){colnames(QC)[col_idx]<-"PAI"}

for(i in 1:nrow(QC)){
  
  print(i)
  if(!is.na(QC$QC_1[i])){
    next
  }
  
  if(is.na(QC$thresh_QC[i])){ #check if there is a number in thresh_QC
    
    #create a new one if not.
    QC$thresh_type[i]<-'Auto'
    QC$img_thd[i]<-NA
    
  }else{
    if(QC$thresh_QC[i]==1){
      QC$thresh_QC[i]<-NA
      QC$thresh_type[i]='Auto'
      QC$QC_1[i]<-1
      print(paste0('skipped ',i,'LAI ok'))
      next
    }else{
      if(QC$thresh_QC[i]==0){
        QC$thresh_QC[i]<-NA
        QC$thresh_type[i]='manual'
        QC$img_thd[i]=NA
        QC$QC_1[i]<-1
        print(paste0('LAI ',i,'is 0'))
        next
      }
      
      if(QC$thresh_QC[i]==2){
        QC$QC_1[i]<-0
        QC$thresh_type[i]=NA
        QC$thresh_QC[i]<-NA
        QC$img_thd[i]=NA
        print(paste0('skipped ',i,' dodgy image'))
        next
      }
      
      if(QC$thresh_QC[i]>=3){
        QC$img_thd[i]=QC$thresh_QC[i]
        QC$thresh_QC[i]<-NA
        QC$thresh_type[i]='manual'
        QC$QC_1[i]<-1
        print(paste0('LAI ',i,'recalculated'))
      }
    }
  }
  
  start<-Sys.time()
  dcp<-QC$file.path[i]
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
  
  # this is for repeat analysis. check if the quality control step is done and adjust the threshold if needed.  
  if(is.na(QC$img_thd[i])){
    img_thd <- unlist(autothresholdr::auto_thresh(round(img.mat), method='Otsu',ignore_na=TRUE))[1]
    QC$img_thd[i] <- img_thd
  }else{
    img_thd <- QC$img_thd[i]
  }
  
  
  img.bw <- terra::classify(rast_im, rbind(c(-Inf,img_thd,0),c( img_thd,Inf,1)))
  base::names(img.bw) <- base::names(rast_im)
  mk<- is.na(img.bw)
  img.bw[mk]<-NA
  # Make a histogram to inspect the classification
  histo<-hist_plot(rast_im,img_thd[1])
  rm(rast_im)
  
  # Determine the extinc coef. Alic looked like 57 deg
  if(QC$Site[i]=='alic'){
    kidx=0.91
  }else{
    kidx=0.5
  }
  
  # Determine the gap fraction and gap size distribution
  if(!is.na(img_thd)){
    outputs <- as.data.frame(getGap(img.bw,path_id=Code,export.image=F,k=kidx,th=img_thd))
  }else{
    outputs <- data.frame(path_id=Code,id='Blue',FC=NA,CC=NA,CP=NA,PAIe=0,PAI=0, CI=NA,k=kidx,imgchannel='Blue',gap_method='Macfarlane', gap_thd=NA, img_method='No plant material',img_thd=NA)
  }
  
  # update the QC spreadsheet
  cns<-colnames(outputs)
  for(k in 1:length(cns)){
    QC[i,cns[k]]<-outputs[1,cns[k]]
  }
  
  # Make a plot for inspection
  imbw<-getGap(img.bw,path_id=Code,export.image=T,th=img_thd)
  final_im<-rgbim+histo+imbw+plot_layout(ncol=2)
  ggsave(final_im,filename=paste0(img_outputs,Code,'.tiff'),width = 30,height=30,units='cm',dpi=72)
  
  end<-Sys.time()
  print(i)
  print(end-start)
  
  # save it every 50 pics just in case of crash
  if(i%%50==0){
    write.csv(QC,'D:/TERN_dcpDown_LAI_QC_complete.csv',na='',row.names = F)
  }
}

col_idx<-which(colnames(QC)=='PAIe')
if(length(col_idx)==1){colnames(QC)[col_idx]<-"GLAIe"}
col_idx<-which(colnames(QC)=='PAI')
if(length(col_idx)==1){colnames(QC)[col_idx]<-"GLAI"}

write.csv(QC,'D:/TERN_dcpDown_LAI_QC_complete.csv',na='',row.names = F)
