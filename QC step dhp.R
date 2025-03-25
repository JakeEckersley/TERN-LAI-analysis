source('hemispheric_PAI_functions.R')
library(terra)
library(hemispheR)
library(ggplot2)
library(patchwork)
library(ggpubr)
library(gtools)
#### notes::
# may need to add the stretch here.
# add Michael's new images last.

file_locations<-c("D:/Tern images")
im_deets<-read.csv(paste0(file_locations,'/Tern image details.csv'))
imlist<-im_deets[im_deets$Method=='dhp'&im_deets$Image.orientation.2 =='upward',]
dhp_folders<-unique(imlist$file.dirs)
mask_list<-readRDS('DHP_mask_list.rds')
#unique(im_deets$Site)
#list of camera lens models I need
unique(paste(imlist$Camera,imlist$Site))
unique(imlist$Camera)

#make into a list
#imlist<-split(imlist,imlist$file.path)
img_outputs<-"D:/Tern binarised images/DHP upward/"


#QC check
QC<-read.csv("C:/Users/22064705/OneDrive - UWA/My Documents/TERN LAI/250803_final_check/TERN_dhp_LAI_QC_final.csv")
idx<-which(!imlist$file.path%in%QC$file.path&imlist$Camera!='unknown-camera')
if(length(idx)>0){
  QC<-smartbind(QC,imlist[idx,])
}

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
    QC$th[i]<-NA
    
  }else{
    if(QC$thresh_QC[i]==1){
      QC$thresh_QC[i]<-NA
      QC$thresh_type[i]='Auto'
      QC$QC_1[i]<-1
      print(paste0('skipped ',i,'LAI ok'))
      next
    }
    
    if(QC$thresh_QC[i]==0){
      QC$thresh_QC[i]<-NA
      QC$thresh_type[i]='manual'
      QC$th[i]=NA
  
      # change LAI to 0
      QC[i,which(colnames(QC)%in%c('PAIe_57','PAI_LX_57','PAIe_70','PAI_LX_70'))]<-0
      QC[i,which(colnames(QC)%in%c('Psat_ring_57','satP_azm_57','CI_57','Psat_ring_70','satP_azm_70','CI_70'))]<-NA
      QC$QC_1[i]<-1
      print(paste0('LAI ',i,'is 0'))
      next
    }
    
    if(QC$thresh_QC[i]==2){
      QC$QC_1[i]<-0
      QC$thresh_type[i]=NA
      QC$thresh_QC[i]<-NA
      QC$th[i]=NA
      print(paste0('skipped ',i,' dodgy image'))
      next
    }
    
    if(QC$thresh_QC[i]>=3){
      QC$th[i]=QC$thresh_QC[i]
      QC$thresh_QC[i]<-NA
      QC$thresh_type[i]='manual'
      QC$QC_1[i]<-1
      print(paste0('LAI ',i,'recalculated'))
    }
  }
  
  
  start<-Sys.time()
  dhp<-QC$file.path[i]
  Code<-tools::file_path_sans_ext(basename(dhp))

  # read in img :
  rast_im<-read_B_and_gamma_correct(dhp,rgb=T)
  
  #plot(rast_im)
  
  # apply mask
  hem_mask<-mask_list[[which(names(mask_list)==QC$file.dirs[i])]]
  xy <- terra::xyFromCell(rast_im,1:terra::ncell(rast_im))
  
  hem_mask<-check_mask_rotation(rast_im,hem_mask)
  
  # check which ones are warra and re-do the print image cos they fucked rn....
  circular.mask = (xy[,1] - hem_mask$center_x)^2 + (xy[,2] - hem_mask$center_y)^2 <= hem_mask$radius^2
  terra::values(rast_im)[!circular.mask] <- NA
  #plot(rast_im)
  
  # create a list of lens details
  lens_details<-'equidistant'
  if(QC$Camera[i]=='Canon-EOS-500D'){lens_details<-'orthographic'} #Canon-EOS-500D assumed orthographic based on email chain with jason and Ian
  if(QC$Camera[i]=='Canon-EOS-7D'){lens_details<-'Sigma-4.5'} # we are unsure on this lens likely this one.
  if(QC$Camera[i]=='Nikon-D5600'){lens_details<-'Sigma-4.5'}
  if(QC$Camera[i]=='Nikon-D5500'){lens_details<-'Sigma-4.5'}
  if(QC$Camera[i]=='Nikon-D600'){lens_details<-'orthographic'} # this is a guess; probably actually the sony alpha as image dims + mask are are exactly the same
  if(QC$Camera[i]=='Nikon-D700'){lens_details<-"Nikkor-10.5"}
  if(QC$Camera[i]=='Sony-Alpha-A600'){lens_details<-'orthographic'}
  ##################################################################
  
  #generate plot 1
  if(nlyr(rast_im)==3){
    rgbim<-rgb_plot(rast_im)
  }else{
    rgbim<-bw_plot(rast_im,title=Code)
  }
  
  # select the blue band
  rast_im<-rast_im[[which(names(rast_im)=='Blue')]]
  
  # Binarise the masked image
  img.mat <- matrix(terra::values(rast_im),nrow = nrow(rast_im),byrow=T)
  

  # this is for repeat analysis. check if the quality control step is done and adjust the threshold if needed.  
  if(is.na(QC$th[i])){
    th <- unlist(autothresholdr::auto_thresh(round(img.mat), method='Otsu',ignore_na=TRUE))[1]
    QC$th[i] <- th
  }else{
    th <- QC$th[i]
  }
  
  img.bw <- terra::classify(rast_im, rbind(c(-Inf,th,0),c( th,Inf,1)))
  base::names(img.bw) <- base::names(rast_im)
  mk<- is.na(img.bw)
  img.bw[mk]<-NA
    
  # Make a histogram to inspect the classification
  histo<-hist_plot(rast_im,th)
  rm(rast_im)
    
  # Make a plot for inspection
  img.bw <- as.data.frame(img.bw, xy= TRUE)
  bw_img<-bw_plot(img.bw,title=Code)
    
  # Segment the image
  seg.img.bw<-segment(img.bw,hem_mask,lens=lens_details)
  rm(img.bw)
  PAI_Nils<-getPAINils_fromImg(seg.img.bw,Code)
  PAI_57<-getPAINils_fromImg(adj_seg(seg.img.bw,hem_mask,lens=lens_details),Code)
    
  PAI_57$CI<-(PAI_57$PAIe/PAI_57$PAI_LX)
  PAI_Nils$CI<-(PAI_Nils$PAIe/PAI_Nils$PAI_LX)
  colnames(PAI_57)<-paste0(colnames(PAI_57),'_57')
  colnames(PAI_Nils)<-paste0(colnames(PAI_Nils),'_70')
    
  outputs<-cbind(PAI_57[,c(1,2,5,6,8)],PAI_Nils[,c(1,2,5,6,8)])
  cns<-colnames(outputs)
  for(k in 1:length(cns)){
    QC[i,cns[k]]<-outputs[1,cns[k]]
  }
  
  final_im<-rgbim+histo+bw_img+plot_layout(ncol=2)
  ggsave(final_im,filename=paste0(img_outputs,Code,'.tiff'),width = 30,height=30,units='cm',dpi=76) 
  
  end<-Sys.time()
  print(end-start)
  
  # save it every 50 pics just in case of crash
  if(i%%50==0){
    write.csv(QC,'D:/TERN_dhp_LAI_QC_complete.csv',na='',row.names = F)
  }
}
colnames(QC)
write.csv(QC,'D:/TERN_dhp_LAI_QC_complete.csv',na='',row.names = F)
