source('hemispheric_PAI_functions_251020update.R')
library(terra)
library(hemispheR)
library(ggplot2)
library(patchwork)
library(ggpubr)
library(gtools)

######################################
# Read in data
######################################
QC <- read.csv("E:/TERN LAI final doc/TERN_dhp_LAI_QC_complete.csv")
mask_data <- readRDS('DHP_mask_list.rds')
ImgDirQC="J:/BinarisedImagesBase/DHP upward/"
dir.create(ImgDirQC,showWarnings=F)

#####################################
# Fix some dicky stuff
#####################################
colnames(QC)[which(colnames(QC)=="QC_1")]<-"QC"
QC$file.path<-gsub('D:/',"J:/",QC$file.path)
QC$file.dirs<-gsub('D:/',"J:/",QC$file.dirs)
names(mask_data)<-gsub('D:/',"J:/",names(mask_data))

QC<-QC[,-which(colnames(QC)%in%c("thresh_QC","FC","CC","CP","PAIe","PAI","CI","k","imgchannel","gap_method","gap_thd","img_method","img_thd"))]
colnames(QC)[which(colnames(QC)=="th")]<-"AutoThd"
QC$ManualThd <- NA
QC$ManualThd[QC$thresh_type=='Manual']<-QC$AutoThd[QC$thresh_type=='Manual']
QC$AutoThd[QC$thresh_type=='Manual']<-NA
QC$Duplicate[QC$Duplicate!=0]<-1

######################################
# Main function
######################################
DHP_analysis<-function(ImgDat, # single row of data
                       mask_list=mask_data,
                       img_outputs=ImgDirQC,
                       PrintImg=T){
  
  start<-Sys.time()
  dhp<-ImgDat$file.path
  Code<-tools::file_path_sans_ext(basename(dhp))
  
  ###################################
  # read in img and check QC
  ###################################
  rast_im<-read_B_and_gamma_correct(dhp,rgb=T,stretch_im = T)
  
  qc_fail<-ifelse(ImgDat$QC==1,F,T)
  if(is.na(qc_fail)){
    ImgDat$QC=NA
    ImgDat$Comments='QC needed'
  }
  
  ################################
  # apply mask
  ################################
  hem_mask<-mask_list[[which(names(mask_list)==ImgDat$file.dirs)]]
  xy <- terra::xyFromCell(rast_im,1:terra::ncell(rast_im))
  
  hem_mask<-check_mask_rotation(rast_im,hem_mask)
  
  # check which ones are warra and re-do the print image cos they fucked rn....
  circular.mask = (xy[,1] - hem_mask$center_x)^2 + (xy[,2] - hem_mask$center_y)^2 <= hem_mask$radius^2
  terra::values(rast_im)[!circular.mask] <- NA
  #plot(rast_im)
  
  #######################################
  # create a list of lens details
  #######################################
  lens_details<-'equidistant'
  if(ImgDat$Camera=='Canon-EOS-500D'){lens_details<-'orthographic'} #Canon-EOS-500D assumed orthographic based on email chain with jason and Ian
  if(ImgDat$Camera=='Canon-EOS-7D'){lens_details<-'Sigma-4.5'} # we are unsure on this lens likely this one.
  if(ImgDat$Camera=='Nikon-D5600'){lens_details<-'Sigma-4.5'}
  if(ImgDat$Camera=='Nikon-D5500'){lens_details<-'Sigma-4.5'}
  if(ImgDat$Camera=='Nikon-D600'){lens_details<-'orthographic'} # this is a guess; probably actually the sony alpha as image dims + mask are are exactly the same
  if(ImgDat$Camera=='Nikon-D700'){lens_details<-"Nikkor-10.5"}
  if(ImgDat$Camera=='Sony-Alpha-A600'){lens_details<-'orthographic'}

  ##################################
  # check the image isn't weird
  ##################################
  if(!nlyr(rast_im)%in%c(1,3)){
    message('odd value for ',Code)
    outputs <- data.frame(path_id=Code,PAIe=NA,img_thd=NA)
    ImgDat <- cbind(ImgDat[c('Site','Image.type','Plot','Date')],outputs,ImgDat[,c('Image.orientation.2','Image.code','Method','Camera','Download.code','file.path','QC','Comments','thresh_type','Duplicate')])
    ImgDat$Comments <- 'Odd number of layers in image'
    ImgDat$QC <- 0
    return(list(ImgDat,NULL))
  }
  
  ##################################
  # select the blue band
  ##################################
  rast_im<-rast_im[[which(names(rast_im)=='Blue')]]
  img.mat <- matrix(terra::values(rast_im),nrow = nrow(rast_im),byrow=T)
  
  #################################
  # Set binary threshold
  #################################
  if(is.na(ImgDat$AutoThd)){ 
    img_thd <- unlist(autothresholdr::auto_thresh(round(img.mat), method='Otsu',ignore_na=TRUE))[1]
    ImgDat$AutoThd <- img_thd
    if(!is.na(ImgDat$ManualThd)){
      img_thd <- ImgDat$ManualThd
    }else{
      ImgDat$Comments<-'QC needed'
      ImgDat$QC<-NA
      ImgDat$thresh_type<-'Auto'
    }
  }else{
    img_thd <- ImgDat$AutoThd
  }
  
  #################################
  # Binarise image
  #################################
  img.bw <- terra::classify(rast_im, rbind(c(-Inf,img_thd,0),c( img_thd,Inf,1)))
  base::names(img.bw) <- base::names(rast_im)
  mk<- is.na(img.bw)
  img.bw[mk]<-NA
  
  #################################
  # Make evaluation plot
  #################################
  if(PrintImg){  
    png(paste0(img_outputs,Code,'.png'), width = 1600, height = 1200, res = 150)
    
    # Set 2x2 layout (use 3 panels)
    par(mfrow = c(2, 2))
    par(mar = c(4, 4, 2, 1))  # reduce margins a bit
    
    ## 1. RGB image 
    if(nlyr(rast_im)==3){
      plotRGB(rast_im, r = 1, g = 2, b = 3, main = if(is.na(qc_fail)){
          'QC needed'
        }else if(qc_fail==1){
          'QC pass'
        }else if(qc_fail==0){
          'QC fail'
        }
        )
    }else{
      plot(rast_im, main = if(is.na(qc_fail)){
        'QC needed'
      }else if(qc_fail==1){
        'QC pass'
      }else if(qc_fail==0){
        'QC fail'
      }
      )
    }
    
    ## --- 2. Histogram of all pixel values (top right)
    vals <- values(rast_im, na.rm = TRUE)
    vals_all <- as.vector(vals)
    hist(vals_all,
         breaks = 50,
         col = "grey80",
         border = "white",
         main = paste0('thresh = ',ImgDat$AutoThd),
         xlab = "Value")
    abline(v = ImgDat$AutoThd, col = "blue", lwd = 2, lty = 2)
    if (!is.na(ImgDat$ManualThd)) {
      abline(v = img_thd, col = "red", lwd = 2, lty = 2)
      
      legend("topright",
             legend = c("AutoThd", "ManualThd"),
             col = c("blue", "red"),
             lty = 2,
             lwd = 2,
             bty = "n")
    } else {
      legend("topright",
             legend = c("AutoThd"),
             col = c("blue"),
             lty = 2,
             lwd = 2,
             bty = "n")
    }
    ## --- 3. Binary blue band (bottom left)
    plot(img.bw, col = c("black", "white"), legend = FALSE, main = paste0('im = ',Code))
    plot.new()
    text(0.5, 0.5, " ")
    
    dev.off()
  }
  
  
  ########################################
  # Segment the image
  ########################################
  img.bw <- as.data.frame(img.bw, xy= TRUE)
  seg.img.bw<-segment(img.bw,hem_mask,lens=lens_details)
  PAI_Nils<-getPAINils_fromImg(seg.img.bw,Code)
  
  #re-segment for hinge
  seg_57<-adj_seg(seg.img.bw,hem_mask,lens=lens_details)
  PAI_57<-getPAINils_fromImg(seg_57,Code)
  
  PAI_57$CI<-(PAI_57$PAIe/PAI_57$PAI_LX)
  PAI_Nils$CI<-(PAI_Nils$PAIe/PAI_Nils$PAI_LX)
  colnames(PAI_57)<-paste0(colnames(PAI_57),'_57')
  colnames(PAI_Nils)<-paste0(colnames(PAI_Nils),'_70')
  
  ####################################################################
  # Compute proportion of gaps per ring and azimuth_group
  ####################################################################
  prop_df <- aggregate(Gap ~ ring + azimuth_group, data = seg.img.bw, FUN = function(x) mean(x == 1))
  names(prop_df )[3] <- "GapFrac"
  prop_df <-suppressWarnings(cbind(prop_df,ImgDat[,c("path_id","Site","Image.type","Plot","Date","Image.orientation.2","Image.code","Method","Camera","thresh_type","QC","Duplicate","Comments")]))
  
  prop_df_2 <- aggregate(Gap ~ ring + azimuth_group, data = seg_57, FUN = function(x) mean(x == 1))
  names(prop_df_2 )[3] <- "GapFrac"
  prop_df_2 <-suppressWarnings(cbind(prop_df_2,ImgDat[,c("path_id","Site","Image.type","Plot","Date","Image.orientation.2","Image.code","Method","Camera","thresh_type","QC","Duplicate","Comments")]))
  
  outputs<-cbind(PAI_57[,c(1,2,5,6,8)],PAI_Nils[,c(1,2,5,6,8)])
  cbind(ImgDat[,c("path_id","Site","Image.type","Plot","Date")],outputs,ImgDat[,c("Image.orientation.2","Image.code","Method","Camera","thresh_type","AutoThd","ManualThd","QC","Duplicate","Comments")])
  
  #############################################
  # Tidy up
  #############################################
  rm(img.bw,vals,rast_im)
  
  end<-Sys.time()
  message(basename(ImgDat$file.path),'time = ',end-start)
  
  terra::tmpFiles(remove = TRUE)
  
  #############################################
  # Return results
  #############################################
  # this is the default - don't include anything in the site level outputs until the image has been inspected.
  if(is.na(ImgDat$QC)){
    qc_fail = 0
  }
  
  if(qc_fail&ImgDat$Duplicate){
    return(list(ImgDat,prop_df,prop_df_2))
  }else{
    return(list(ImgDat,NULL,NULL))
  }
}



##############################
# Parallel processing
##############################

split_QC <- split(QC, list(QC$Site,QC$Plot,QC$Date))
split_QC<-split_QC[unlist(lapply(split_QC,function(x)nrow(x)>0))]


library(parallel)
library(future.apply)
library(progressr)
#https://henrikbengtsson.github.io/future-tutorial-user2022/reporting-on-progress-updates.html



# # the following parallelises the application of the dhp function. the same can be achieved by running the lapply function in the 'sequential processing' section

# define the speedy processing function with progress reporter..
FastApply <- function(X) {
  # define progressor of total number of rows across all chunks
  total_rows <- sum(sapply(X, nrow))
  p <- progressr::progressor(steps = total_rows)
  
  
  future_lapply(
    X,
    function(df_chunk) {
      results <- vector("list", nrow(df_chunk))
      
      for (i in seq_len(nrow(df_chunk))) {
        results[[i]] <- DHP_analysis(df_chunk[i, ])
        p(message = sprintf("Row %d of %d", i, nrow(df_chunk)))
      }
      return(results)
      
    },
    future.seed = TRUE,
    future.packages = c('terra','hemispheR','ggplot2','patchwork','ggpubr','gtools')
  )
}



handlers(global = TRUE)   # enable global progress reporting


plan(multisession, workers = 20)

  LAI_chunks <- FastApply(split_QC)

plan(sequential)

# end parallel processing
# LAI_chunks[[1]][[1]][[1]]

LAI_df <- do.call(smartbind,
  lapply(LAI_chunks,function(chunk){
    chunk<-lapply(chunk,function(row){
      row<-row[[1]]
    })
    chunk<-do.call(smartbind,chunk)
  })
)

Gaps_70 <- lapply(LAI_chunks,function(chunk){
  chunk<-lapply(chunk,function(row){
    row<-row[[2]]
  })
  chunk<-do.call(smartbind,chunk)
})

Gaps_57 <- lapply(LAI_chunks,function(chunk){
  chunk<-lapply(chunk,function(row){
    row<-row[[3]]
  })
  chunk<-do.call(smartbind,chunk)
})

write.csv(LAI_df,'J:/TERN_dhp_LAI_QC_complete_251026update.csv',na='',row.names = F)
saveRDS(Gaps_70,'J:/TERN_dhp_RawGapFrac_70.rds')
saveRDS(Gaps_57,'J:/TERN_dhp_RawGapFrac_57.rds')
saveRDS(LAI_chunks,'J:/TERN_dhp_RawFastApplyList.rds')

##########################################
### Sequential processing
##########################################
# LAI_chunks <- lapply(split_QC,function(df_chunk) {
#     results <- vector("list", nrow(df_chunk))
#     for (i in seq_len(nrow(df_chunk))) {
#       results[[i]] <- DHP_analysis(df_chunk[i, ])
#       p(message = sprintf("Row %d of %d", i, nrow(df_chunk)))
#     }
#     return(results)
#   })


