source('hemispheric_PAI_functions.R')
library(terra)
library(hemispheR)
library(ggplot2)
library(patchwork)
library(ggpubr)
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
imlist<-split(imlist,imlist$file.path)
img_outputs<-"D:/Tern binarised images/DHP upward/"
#i=1596



#QC check
QC<-read.csv("C:/Users/22064705/OneDrive - UWA/My Documents/TERN LAI/TERN_dhp_LAI.csv")

canopy<-readRDS('D:/TERN_dhp_LAI.rds')
for(i in 3151:length(imlist)){
  start<-Sys.time()
  dhp_files<-imlist[[i]]
  dhp<-dhp_files$file.path
  Code<-tools::file_path_sans_ext(basename(dhp))
  
  if(grepl('unknown',dhp_files$Camera)){next}

  # read in img :
  rast_im<-read_B_and_gamma_correct(dhp,rgb=T)
  
  #plot(rast_im)
  
  # apply mask
  hem_mask<-mask_list[[which(names(mask_list)==dhp_files$file.dirs)]]
  xy <- terra::xyFromCell(rast_im,1:terra::ncell(rast_im))
  
  hem_mask<-check_mask_rotation(rast_im,hem_mask)
  
  # check which ones are warra and re-do the print image cos they fucked rn....
  circular.mask = (xy[,1] - hem_mask$center_x)^2 + (xy[,2] - hem_mask$center_y)^2 <= hem_mask$radius^2
  terra::values(rast_im)[!circular.mask] <- NA
  

  # create a list of lens details
  lens_details<-'equidistant'
  if(dhp_files$Camera=='Canon-EOS-500D'){lens_details<-'orthographic'} #Canon-EOS-500D assumed orthographic based on email chain with jason and Ian
  if(dhp_files$Camera=='Canon-EOS-7D'){lens_details<-'Sigma-4.5'} # we are unsure on this lens likely this one.
  if(dhp_files$Camera=='Nikon-D5600'){lens_details<-'Sigma-4.5'}
  if(dhp_files$Camera=='Nikon-D5500'){lens_details<-'Sigma-4.5'}
  if(dhp_files$Camera=='Nikon-D600'){lens_details<-'orthographic'} # this is a guess; probably actually the sony alpha as image dims + mask are are exactly the same
  if(dhp_files$Camera=='Nikon-D700'){lens_details<-"Nikkor-10.5"}
  if(dhp_files$Camera=='Sony-Alpha-A600'){lens_details<-'orthographic'}
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

  skip_im=F
  # this is for repeat analysis. check if the quality control step is done and adjust the threshold if needed.  
  QCidx<-which(QC$file.path==dhp)
  if(length(QCidx)==1){
    if(!is.na(QC[QCidx,]$th)){
      th<-QC[QCidx,]$th
      if(!is.na(QC[QCidx,]$thresh_QC)){
        if(QC[QCidx,]$thresh_QC==0){
          th<-NA
          dhp_files$QC_1<-1
          dhp_files$thresh_type<-'Manual'
          skip_im=T
        }
        if(QC[QCidx,]$thresh_QC==1){
          canopy[[i]]$thresh_QC<-NA
          canopy[[i]]$thresh_type='Auto'
          canopy[[i]]$QC_1<-1
          print(paste0('skipped ',i,'LAI ok'))
          next
        }
        if(QC[QCidx,]$thresh_QC==2){
          canopy[[i]]$QC_1<-0
          canopy[[i]]$thresh_type=NA
          canopy[[i]]$thresh_QC<-NA
          print(paste0('skipped ',i,' dodgy image'))
          next
        }
        if(QC[QCidx,]$thresh_QC>2){
          dhp_files$QC_1<-1
          dhp_files$thresh_type<-'Manual'
          th<-QC[QCidx,]$thresh_QC
        }
      }
    }else{
      print(paste(i,'no manual threshold - creating new one'))
      th <- autothresholdr::auto_thresh(round(img.mat), method='Otsu',ignore_na=TRUE)
      th <- unlist(th)[1]
      dhp_files$thresh_type<-'Auto'
      skip_im=T
    }
  }else{
    th <- autothresholdr::auto_thresh(round(img.mat), method='Otsu',ignore_na=TRUE)
    th <- unlist(th)[1]
    dhp_files$thresh_type<-'Auto'
  }
  
  if(!is.na(th)){
    img.bw <- terra::classify(rast_im, rbind(c(-Inf,th,0),c( th,Inf,1)))
    base::names(img.bw) <- base::names(rast_im)
    mk<- is.na(img.bw)
    img.bw[mk]<-NA
    
    # Make a histogram to inspect the classification
    histo<-hist_plot(rast_im,th[1])
    rm(rast_im)
  
    # Make a plot for inspection
    img.bw <- as.data.frame(img.bw, xy= TRUE)
    bw_img<-bw_plot(img.bw,title=Code)
  
    # Segment the image
    seg.img.bw<-segment(img.bw,hem_mask,lens=lens_details)
    rm(img.bw)
    PAI_Nils<-getPAINils_fromImg(seg.img.bw,Code)
    PAI_57<-getPAINils_fromImg(adj_seg(seg.img.bw,hem_mask,lens=lens_details),Code)
  
    PAI_57$CC<-(PAI_57$PAIe/PAI_57$PAI_LX)
    PAI_Nils$CC<-(PAI_Nils$PAIe/PAI_Nils$PAI_LX)
    colnames(PAI_57)<-paste0(colnames(PAI_57),'_57')
    colnames(PAI_Nils)<-paste0(colnames(PAI_Nils),'_70')
    
    outputs<-cbind(data.frame(path_id=Code,id='Blue'),PAI_57[,c(1,2,5,6,8)],PAI_Nils[,c(1,2,5,6,8)],data.frame(th=th))
  }else{
    df1<-data.frame(PAIe_57=0,PAI_LX_57=0,Psat_ring_57=NA,satP_azm_57=NA,CC_57=NA)
    df2<-data.frame(PAIe_70=0,PAI_LX_70=0,Psat_ring_70=NA,satP_azm_70=NA,CC_70=NA)
    outputs<-cbind(data.frame(path_id=Code,id='Blue'),df1,df2,data.frame(th=th))
  }
  canopy[[i]] <- cbind(outputs,dhp_files)

  final_im<-rgbim+histo+bw_img+plot_layout(ncol=2)
  
  if(dhp_files$Site=="wrra"){
    skip_im=F
  }
  if(skip_im==F){
    ggsave(final_im,filename=paste0(img_outputs,Code,'.tiff'),width = 30,height=30,units='cm',dpi=76) 
  }
  end<-Sys.time()
  print(i)
  print(end-start)
  
  # save it every 50 pics just in case of crash
  if(i%%50==0){
    saveRDS(canopy,'D:/TERN_dhp_LAI.rds')
  }
}

saveRDS(canopy,'D:/TERN_dhp_LAI.rds')

# read in the previous one before doing this - will have columns for QC
# change CC to CI throughout to be consistent with DHP
write.csv(do.call(rbind,canopy),'D:/TERN_dhp_LAI.csv',row.names=F,na='')
