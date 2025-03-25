########### 57degree facing dcp photos
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
QC<-read.csv("C:/Users/22064705/OneDrive - UWA/My Documents/TERN LAI/250803_final_check/TERN_dcp_57_LAI_QC.csv")

img_outputs<-"D:/Tern binarised images/DCP 57/"

for(i in 1:nrow(QC)){
  
  print(i)
  if(is.na(QC$Duplicate[i])){
    QC$Duplicate[i]<-1
  }
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
  
  if(nlyr(rast_im)==3){
    rgbim<-rgb_plot(rast_im)
  }else{
    rgbim<-bw_plot(rast_im,title = '')
  }
  
  # select the blue band
  rast_im<-rast_im[[which(names(rast_im)=='Blue')]]
  
  # Binarise the image
  img.mat <- matrix(terra::values(rast_im),nrow = nrow(rast_im),byrow=T)
  
  # this is for repeat analysis.
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
  
  # Determine the gap fraction and gap size distribution
  if(!is.na(img_thd)){
    outputs <- as.data.frame(getGap(img.bw,path_id=Code,export.image=F,k=0.91,th=img_thd))
  }else{
    outputs <- data.frame(path_id=Code,id='Blue',FC=NA,CC=NA,CP=NA,PAIe=0,PAI=0, CI=NA,k=0.91,imgchannel='Blue',gap_method='Macfarlane', gap_thd=NA, img_method='No plant material',img_thd=NA)
  }
  
  # update QC csv
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
    write.csv(QC,'D:/TERN_dcp_57_LAI_QC_complete.csv',na='',row.names = F)
  }
}

write.csv(QC,'D:/TERN_dcp_57_LAI_QC_complete.csv',na='',row.names = F)
