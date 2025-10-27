# take the 're-do' ones out of dcp_QC
# change whroo to 0 in dhp LAI set....
# add a comment that there are too many images for the two belinda raised.


# check how many there are in the 12/03/2015 dataset. look to be too few? add the error bars...

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
library(mgc)
file_locations<-c("J:/Tern images")
im_deets<-read.csv(paste0(file_locations,"/re-do subset.csv"))
dcpUpdate<-QC[-which(QC$file.path%in%im_deets$file.path),]

QC<-read.csv("E:/TERN LAI final doc/TERN_dcp_LAI_QC_complete.csv")
dcpUpdate<-QC[-which(QC$file.path%in%im_deets$file.path),]
QC<-QC[which(QC$file.path%in%im_deets$file.path),]
#QC$Duplicate<-1
img_outputs<-"E:/Tern binarised images/DCP 57deg/"

for(i in 1:nrow(QC)){

  print(i)
  
  if(!is.na(QC$QC[i])){
    if(QC$k[i]==0.5){
      if(QC$QC[i]==0){next}
      print('re-do dodgy k')
    }else{
      next
    }
  }
  
  if(QC$QC[i]!=1){
    print('skip - failed QC')
    next
  }
  
  start<-Sys.time()
  dcp<-QC$file.path[i]
  Code<-tools::file_path_sans_ext(basename(dcp))
  rast_im<-read_B_and_gamma_correct(dcp,rgb=T,stretch_im = T)
  
  
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
    write.csv(QC,'E:/TERN_dcp_57_LAI_QC_updated250926.csv',na='',row.names = F)
  }
  
  terra::tmpFiles(remove = TRUE)
  
}
QC2<-read.csv("E:/TERN LAI final doc/TERN_dcp_LAI_QC_complete.csv")
nrow(QC2)

lai_data_57 <- read.csv("E:/TERN_dcp_57_LAI_QC_updated250926.csv")

colnames(QC)[!colnames(QC)%in%colnames(lai_data_57)]
colnames(lai_data_57)[!colnames(lai_data_57)%in%colnames(QC)]

QC<-rbind(QC,lai_data_57)

QC$img_thd[!is.na(QC$th)]
QC[is.infinite(QC$PAIe),'Comments']<-'No plant material'
QC[is.infinite(QC$PAIe),c('PAIe','PAI')]<-0
QC$CanopyProp[QC$CanopyProp==501]<-0

write.csv(QC,'E:/TERN_dcp_57_LAI_QC_updated250926.csv',na='',row.names = F)

#final step once dcp has finished - read in 57 and remove any overlapping file paths in the upward dcp
write.csv(dcpUpdate,'E:/TERN_dcp_LAI_QC_updated250926.csv',na='',row.names = F)



# fix whroo
dhp <- read.csv("E:/TERN LAI final doc/TERN_dhp_LAI_QC_complete.csv")
dhp$Date <- as.Date(as.character(dhp$Date), format = "%Y%m%d")

dhp$QC_1[dhp$Site=="vicd_whroo"&dhp$Plot=="core1ha"&dhp$Date>"2017-01-01"]<-0
write.csv(dhp,'E:/TERN_dhp_LAI_QC_updated250926.csv',na='',row.names = F)
