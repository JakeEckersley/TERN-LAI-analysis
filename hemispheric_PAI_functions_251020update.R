#determine appropriate mask for DHCP images
read_B_and_gamma_correct<-function(dhp=dcp,gamma=2.2,rgb=T,stretch_im=F){
  imgB <- suppressWarnings(terra::rast(dhp))# remove the spatial data reference warning...
  if(nlyr(imgB)==1){
    img.valuesB<-terra::values(imgB)
    minmB<-min(img.valuesB,na.rm=TRUE)
    maxmB<-max(img.valuesB,na.rm=TRUE)
    img.valuesB=(maxmB-minmB)*((img.valuesB/(maxmB-minmB))^gamma)
    terra::values(imgB)<-img.valuesB
    names(imgB)<-'Blue'
  }else{
    if(nlyr(imgB)>1&rgb==F){
      imgB<-imgB[[3]]
      img.valuesB<-terra::values(imgB)
      minmB<-min(img.valuesB,na.rm=TRUE)
      maxmB<-max(img.valuesB,na.rm=TRUE)
      img.valuesB=round((maxmB-minmB)*((img.valuesB/(maxmB-minmB))^gamma))
      terra::values(imgB)<-img.valuesB
      names(imgB)<-'Blue'
    }
    if(nlyr(imgB)>1&rgb==T){
      for(i in 1:nlyr(imgB)){
        img.valuesB<-terra::values(imgB[[i]])
        minmB<-min(img.valuesB,na.rm=TRUE)
        maxmB<-max(img.valuesB,na.rm=TRUE)
        img.valuesB=round((maxmB-minmB)*((img.valuesB/(maxmB-minmB))^gamma))
        terra::values(imgB[[i]])<-img.valuesB
      }
      names(imgB)<-c('Red','Green','Blue')
    }
  }
  
  
  if(stretch_im==T){
    imgB<-terra::stretch(imgB,minv=0,maxv=255,minq=0.01,maxq=0.99)
  }
  return(imgB)
}


# adjust_threshold <- function(T_orig, dhp, gamma = 2.2, qmin = 0.01, qmax = 0.99) {
#   imgB <- suppressWarnings(terra::rast(dhp))# remove the spatial data reference warning..
#   
#   if(nlyr(imgB)==1){
#     img_values<-terra::values(imgB)
#   }else{
#     img_values<-terra::values(imgB[[3]])
#   }
#   
# 
#   minv <- min(img_values, na.rm=TRUE)
#   maxv <- max(img_values, na.rm=TRUE)
#   
#   # Apply gamma correction
#   T_gamma <- ((T_orig - minv) / (maxv - minv))^gamma * (maxv - minv) + minv
#   
#   # Apply stretch
#   qlo <- quantile(img_values, qmin, na.rm=TRUE)
#   qhi <- quantile(img_values, qmax, na.rm=TRUE)
#   
#   T_stretch <- 255 * (T_gamma - qlo) / (qhi - qlo)
#   return(T_stretch)
# }



find_mask_FF<-function(dhp_rast,portrait){

  mask<-data.frame(
    center_x=xmax(dhp_rast)/2,
    center_y=ymax(dhp_rast)/2,
    radius=round(sqrt((xmax(dhp_rast)/2)^2+(ymax(dhp_rast)/2)^2))-2,
    portrait = portrait # 0 for square, 1 for portrait, 2 for landscape
  )
  crop_ext<-ext(dhp_rast)
  mask<-cbind(mask,as.data.frame(t(unlist(as.list(crop_ext)))))
  return(mask)
}

find_mask<-function(crop_ext,portrait){
  xrad<-(xmax(crop_ext) - xmin(crop_ext))/2
  yrad<-(ymax(crop_ext) - ymin(crop_ext))/2
  if(yrad!=xrad){
    if(which(c(xrad,yrad)==min(c(xrad,yrad)))==1){
      rad<-xrad
    }else{
      rad<-yrad
    }
  }else{
    rad<-yrad
  }


  mask<-data.frame(
    center_x = floor((xmin(crop_ext) + xmax(crop_ext)) / 2),
    center_y = floor((ymin(crop_ext) + ymax(crop_ext)) / 2),
    radius = floor(rad),
    portrait = portrait # 0 for square, 1 for portrait, 2 for landscape
  )
  
  mask<-cbind(mask,as.data.frame(t(unlist(as.list(crop_ext)))))
  return(mask)
}

plot_with_circ<-function(dhp_rast,mask,crop_ext=NULL,circ=T){
  if(is.null(crop_ext)){crop_ext<-ext(dhp_rast)}
  if(nlyr(dhp_rast)==1){
    plot(terra::crop(dhp_rast,crop_ext), colNA = "black")
  }else{
    plotRGB(terra::crop(dhp_rast,crop_ext), colNA = "black")
  }
  if(circ==T){
    dd <- dismo::circles(cbind(mask$center_x,mask$center_y), mask$radius, lonlat=F)
    plot(dd,  add=TRUE, border='red', lwd=3)
    graphics::segments(mask$center_x,mask$center_y,mask$center_x+mask$radius,mask$center_y,col='red',lwd=3)
  }
  
  if(circ==F){
    graphics::segments(mask$center_x,mask$center_y,terra::ext(dhp_rast)[2],terra::ext(dhp_rast)[4],col='red',lwd=3)
  }
}


# next bit - plotting and image safety
check_mask_rotation<-function(dhp_rast,mask){
  if(mask$portrait==0){
    return(mask)
  }else{
    width_height<-c(ncol(dhp_rast),nrow(dhp_rast))
    long_side<-which(width_height==max(width_height))
    short_side<-which(width_height==min(width_height))
    if(long_side==1&mask$portrait==1){
      colnames(mask)[c(1,2)]<-colnames(mask)[c(2,1)]
      return(mask)
    }
    if(long_side==2&mask$portrait==2){
      colnames(mask)[c(1,2)]<-colnames(mask)[c(2,1)]
      return(mask)
    }
    return(mask)
  }
}

# rgb_plot<-function(dhp_rast=rast_im,failed_QC=F){
#   if("INT1U"%in%datatype(dhp_rast)){
#     maxcolval<-255
#   }else{
#     maxcolval <- global(dhp_rast, fun = "max", na.rm = TRUE)$max[1]
#   }
#   RGB(dhp_rast)<-c(1,2,3)
#   imdf <- as.data.frame(dhp_rast, xy= TRUE)
#   if(nlyr(dhp_rast)==3){
#     gg3.1<-ggplot(data = imdf, aes(x = x, y = y))+
#       geom_raster(fill = rgb(r = imdf$Red,
#                              g = imdf$Green,
#                              b = imdf$Blue,
#                              maxColorValue = max(maxcolval)),show.legend = FALSE) +
#       scale_fill_identity()+theme_bw()+
#       theme(line = element_blank(),
#             axis.title.x=element_blank(),
#             axis.text.x=element_blank(),
#             axis.ticks.x=element_blank(),
#             axis.title.y=element_blank(),
#             axis.text.y=element_blank(),
#             axis.ticks.y=element_blank(),
#             legend.position = "none")+
#       scale_x_continuous(minor_breaks = seq(min(imdf$x), max(imdf$x), by = 10))
#     
#     if (failed_QC) {
#       gg3.1 <- gg3.1 + ggtitle("failed QC")
#     }
#     
#     return(gg3.1)
#   }else{
#     gg3.1<-ggplot(data = imdf, aes(x = x, y = y,fill=Blue))+
#       geom_raster(show.legend = FALSE) +
#       scale_fill_grey()+theme_bw()+
#       theme(line = element_blank(),
#             axis.title.x=element_blank(),
#             axis.text.x=element_blank(),
#             axis.ticks.x=element_blank(),
#             axis.title.y=element_blank(),
#             axis.text.y=element_blank(),
#             axis.ticks.y=element_blank(),
#             legend.position = "none")+
#       scale_x_continuous(minor_breaks = seq(min(imdf$x), max(imdf$x), by = 10))
# 
#     if (failed_QC) {
#       gg3.1 <- gg3.1 + ggtitle("failed QC")
#     }
#     return(gg3.1)
#   }
# }
# 
# hist_plot<-function(dhp_rast,thresh){
#   dhp_rast<-na.omit(as.data.frame(values(dhp_rast)))
#   # basic histogram
#   plt.3 <- ggplot(data=dhp_rast,aes(x=Blue)) + 
#     geom_histogram(bins = 255)+theme_bw()+
#     geom_vline(xintercept = thresh)+ggtitle(paste0('Thresh = ',thresh))
#   return(plt.3)
# }
# 
# bw_plot<-function(img.bw,title){
#   img.bw<-as.data.frame(img.bw,xy=T)
#   gg3.2<-ggplot(data = img.bw, aes(x = x, y = y,fill=Blue))+
#     geom_raster()+
#     scale_fill_gradient(low = "black", high = "white")+
#     theme_bw()+
#     theme(line = element_blank(),
#           axis.title.x=element_blank(),
#           axis.text.x=element_blank(),
#           axis.ticks.x=element_blank(),
#           axis.title.y=element_blank(),
#           axis.text.y=element_blank(),
#           axis.ticks.y=element_blank(),
#           legend.position = "none")+ggtitle(paste0('im = ',title))
#   return(gg3.2)
# }

getPAIdcp<-function(img.bw,ImgID=Code,k=0.5,ClumpMethod="Both",GapThdMcf=0.013,GapListOut=T){
  
  if(k > 1){
    if(message==TRUE){
      warning(paste0("The extinction coefficient (k) was set to ",
                     k,
                     ". It ranges typically between 0.4 and 0.9"))
    }
  }
  
  #### determine gap size distribution:
  vals <- matrix(terra::values(img.bw,format='matrix'),nrow = nrow(img.bw),byrow=T)
  y <- mgc::ConnCompLabel(vals)
  yr <- terra::rast(nrows = nrow(img.bw), ncols=ncol(img.bw),vals=y)
  ext <- terra::ext(img.bw)
  terra::set.ext(yr, ext)
  terra::set.names(yr, base::names(img.bw))
  
  # make frequency table
  GapList <- data.frame(table(terra::values(yr)))
  colnames(GapList)[1:2]<-c('GapID','GapSize')
  GapList$band <- names(yr)
  n_pixels <- sum(GapList$GapSize)
  GapList$n_pixels <- n_pixels
  GapList$Code = ImgID
  
  #############################################################
  # classify gaps using the macfarlane method.
  #############################################################
  if(ClumpMethod=='Mcf'|ClumpMethod=='Both'){
      
    ThMcf <- round(n_pixels*GapThdMcf)
    GapList$GapClass <- ifelse(
      GapList$GapID == "0", "Canopy",
      ifelse(GapList$GapSize >= ThMcf, "LargeGap", "SmallGap")
    )
    
    GapPropMcf <- data.frame(
      Code=ImgID,
      Canopy=sum(GapList$GapSize[GapList$GapClass=='Canopy']/n_pixels),
      SmallGap=sum(GapList$GapSize[GapList$GapClass=='SmallGap']/n_pixels),
      LargeGap=sum(GapList$GapSize[GapList$GapClass=='LargeGap']/n_pixels)
    )
  }
  
  #################################################
  # classify gaps using the Alverni method
  #################################################
  if(ClumpMethod=='Alv'|ClumpMethod=='Both'){
    OnlyGaps <- GapList$GapSize[GapList$GapID != "0"]
    GapThdAlv <- round(mean(OnlyGaps) + sd(OnlyGaps) / sqrt(length(OnlyGaps)),0)   # Alivernini threshold
    GapList$GapClassAlv <- ifelse(
      GapList$GapID == "0", "Canopy",
      ifelse(GapList$GapSize >= GapThdAlv, "LargeGap", "SmallGap")
    )
    GapPropAlv <- data.frame(
      Code=ImgID,
      Canopy=sum(GapList$GapSize[GapList$GapClassAlv=='Canopy']/n_pixels),
      SmallGap=sum(GapList$GapSize[GapList$GapClassAlv=='SmallGap']/n_pixels),
      LargeGap=sum(GapList$GapSize[GapList$GapClassAlv=='LargeGap']/n_pixels),
      GapThdAlv = GapThdAlv/n_pixels
    )
  }
  
  #################################################
  # Generate outputs
  #################################################
  if(ClumpMethod=='Mcf'|ClumpMethod=='Both'){
      
    #Canopy attributes
    CrownCov = 1 - GapPropMcf$LargeGap
    GapFrac = 1 - GapPropMcf$Canopy
    CrownPor = 1 - GapPropMcf$Canopy/CrownCov
    
    # PAI calcs
    PAI_Mcf = -mean(CrownCov)*log(mean(CrownPor))/k
    PAIe = -log(mean(GapFrac))/k
    PAI_LX = -mean(log(GapFrac))/k
    PAI_McfLX = -mean(CrownCov)*mean(log(CrownPor))/k
    
    # Clumping Calcs
    ClumpMcf = ((1-mean(CrownPor))*log(mean(GapFrac)))/(log(mean(CrownPor))*mean(GapPropMcf$Canopy))
    ClumpEffMcf = PAIe/PAI_Mcf
    ClumpEffLX = PAIe/PAI_LX
    ClumpEffMcfLX = PAIe/PAI_McfLX
    
    #################################################
    # Generate outputs
    #################################################
    OutputMcf<-data.frame(
      PAIe = PAIe,
      PAI_Mcf = PAI_Mcf,
      PAI_LX = PAI_LX,
      PAI_McfLX = PAI_McfLX,
      ClumpMcf = ClumpMcf,
      ClumpEffMcf = ClumpEffMcf,
      ClumpEffLX = ClumpEffLX,
      ClumpEffMcfLX = ClumpEffMcfLX,
      CanopyFrac = mean(GapPropMcf$Canopy),
      SmallGapMcf = mean(GapPropMcf$SmallGap),
      LargeGapMcf = mean(GapPropMcf$LargeGap),
      GapFrac = mean(GapFrac),
      CrownPorMcf = mean(CrownPor),
      CrownCovMcf = mean(CrownCov)
    )

    # if(terra::nlyr(img.bw)>1){
    #   UncertMcf <- data.frame(
    #     nSamples = terra::nlyr(img.bw),
    #     CanopyFracSD = sd(GapPropMcf$Canopy),
    #     SmallGapMcfSD = sd(GapPropMcf$SmallGap),
    #     LargeGapMcfSD = sd(GapPropMcf$LargeGap),
    #     GapFracSD = sd(GapFrac),
    #     CrownPorMcfSD = sd(CrownPor),
    #     CrownCovMcfSD = sd(CrownCov)
    #   )
    #   OutputMcf <- cbind(OutputMcf,UncertMcf)
    # }
  }
  
  
  
  if(ClumpMethod=='Alv'|ClumpMethod=='Both'){
    
    #Canopy attributes
    CrownCov = 1 - GapPropAlv$LargeGap
    GapFrac = 1 - GapPropAlv$Canopy
    CrownPor = 1 - GapPropAlv$Canopy/CrownCov
    
    # PAI calcs
    PAI_Alv = -mean(CrownCov)*log(mean(CrownPor))/k
    PAIe = -log(mean(GapFrac))/k
    PAI_LX = -mean(log(GapFrac))/k
    PAI_AlvLX = -mean(CrownCov)*mean(log(CrownPor))/k
    
    # Clumping Calcs
    ClumpAlv = ((1-mean(CrownPor))*log(mean(GapFrac)))/(log(mean(CrownPor))*mean(GapPropAlv$Canopy))
    ClumpEffAlv = PAIe/PAI_Alv
    ClumpEffLX = PAIe/PAI_LX
    ClumpEffAlvLX = PAIe/PAI_AlvLX
    
    #################################################
    # Generate outputs
    #################################################
    OutputAlv<-data.frame(
      PAIe = PAIe,
      PAI_Alv = PAI_Alv,
      PAI_LX = PAI_LX,
      PAI_AlvLX = PAI_AlvLX,
      ClumpAlv = ClumpAlv,
      ClumpEffLX = ClumpEffLX,
      ClumpEffAlvLX=ClumpEffAlvLX,
      CanopyFracAlv = mean(GapPropAlv$Canopy),
      SmallGapAlv = mean(GapPropAlv$SmallGap),
      LargeGapAlv = mean(GapPropAlv$LargeGap),
      GapFrac = mean(GapFrac),
      CrownPorAlv = mean(CrownPor),
      CrownCovAlv = mean(CrownCov)
    )
    
    # if(terra::nlyr(img.bw)>1){
    #   UncertAlv <- data.frame(
    #     nSamples = terra::nlyr(img.bw),
    #     CanopyFracSD = sd(GapPropAlv$Canopy),
    #     SmallGapSD = sd(GapPropAlv$SmallGap),
    #     LargeGapSD = sd(GapPropAlv$LargeGap),
    #     GapFracSD = sd(GapFrac),
    #     CrownPorSD = sd(CrownPor),
    #     CrownCovSD = sd(CrownCov)
    #   )
    #   OutputAlv <- cbind(OutputAlv,UncertAlv)
    # }
    
  }

  
  if(ClumpMethod=='Mcf'){
    output <- OutputMcf
  }else if(ClumpMethod=='Alv'){
    output <- OutputAlv
  }else if(ClumpMethod=='Both'){
    output <- cbind(OutputMcf,OutputAlv[,-which(colnames(OutputAlv)%in%colnames(OutputMcf))])
  }
  
  if(terra::nlyr(img.bw)==1){
    output <- cbind(data.frame(Code=ImgID),output)
    output <- output[,-which(grepl("LX",colnames(output)))]
  }

  if(ClumpMethod=='Mcf'|ClumpMethod=='Both'){
    output$GapThMcf <- GapThdMcf
  }
  if(ClumpMethod=='ALv'|ClumpMethod=='Both'){
    output$GapThAlv <- mean(GapPropAlv$GapThdAlv)
  }
  if(GapListOut==T){
    return(list(output,GapList)) # this is updated. fix the rest of the script...
  }else{
    return(output) # this is updated. fix the rest of the script...
  }

}


getPAIdcp.df<-function(GapList=LAI_gaps[[1]],k=0.5,ClumpMethod="All",GapThdMcf=0.013,GapListOut=F){
  
  if(k > 1){
    if(message==TRUE){
      warning(paste0("The extinction coefficient (k) was set to ",
                     k,
                     ". It ranges typically between 0.4 and 0.9"))
    }
  }
  
  nSamples = sum(GapList$GapClass=='Canopy')
  n_pixels <- mean(GapList$n_pixels)
  ImgID <- GapList$Code[1]
  
  #############################################################
  # classify gaps using the macfarlane method.
  #############################################################
  if(ClumpMethod=='Mcf'|ClumpMethod=='All'){
    
    ThMcf <- round(n_pixels*GapThdMcf)
    GapList$GapClass <- ifelse(
      GapList$GapID == "0", "Canopy",
      ifelse(GapList$GapSize >= ThMcf, "LargeGap", "SmallGap")
    )
    
    GapPropMcf <- data.frame(
      Code=ImgID,
      Canopy=sum(GapList$GapSize[GapList$GapClass=='Canopy']/n_pixels)/nSamples,
      SmallGap=sum(GapList$GapSize[GapList$GapClass=='SmallGap']/n_pixels)/nSamples,
      LargeGap=sum(GapList$GapSize[GapList$GapClass=='LargeGap']/n_pixels)/nSamples
    )
  }
  
  #################################################
  # classify gaps using the Alverni method
  #################################################
  if(ClumpMethod=='Alv'|ClumpMethod=='All'){
    OnlyGaps <- GapList$GapSize[GapList$GapID != "0"]
    GapThdAlv <- round(mean(OnlyGaps) + sd(OnlyGaps) / sqrt(length(OnlyGaps)),0)   # Alivernini threshold
    GapList$GapClassAlv <- ifelse(
      GapList$GapID == "0", "Canopy",
      ifelse(GapList$GapSize >= GapThdAlv, "LargeGap", "SmallGap")
    )
    GapPropAlv <- data.frame(
      Code=ImgID,
      Canopy=sum(GapList$GapSize[GapList$GapClassAlv=='Canopy']/n_pixels)/nSamples,
      SmallGap=sum(GapList$GapSize[GapList$GapClassAlv=='SmallGap']/n_pixels)/nSamples,
      LargeGap=sum(GapList$GapSize[GapList$GapClassAlv=='LargeGap']/n_pixels)/nSamples,
      GapThdAlv = GapThdAlv/n_pixels
    )
  }
  
  if(ClumpMethod=='LX'|ClumpMethod=='All'){
    GapPropLX <- list()
    for(i in unique(GapList$ImgIdx)){
      GapPropLX[[i]] <- data.frame(
        ImgIdx=i,
        Canopy=sum(GapList$GapSize[GapList$GapID=='0'&GapList$ImgIdx==i]/n_pixels),
        Gap=sum(GapList$GapSize[GapList$GapID!='0'&GapList$ImgIdx==i]/n_pixels)
      )
    }
    GapPropLX<-do.call(rbind,GapPropLX)
  }
  
  #################################################
  # Generate outputs
  #################################################
  if(ClumpMethod=='Mcf'|ClumpMethod=='All'){
    
    #Canopy attributes
    CrownCov = 1 - GapPropMcf$LargeGap
    GapFrac = 1 - GapPropMcf$Canopy
    CrownPor = 1 - GapPropMcf$Canopy/CrownCov
    
    # PAI calcs
    PAI_Mcf = -mean(CrownCov)*log(mean(CrownPor))/k
    PAIe = -log(mean(GapFrac))/k
    
    # Clumping Calcs
    ClumpMcf = ((1-mean(CrownPor))*log(mean(GapFrac)))/(log(mean(CrownPor))*mean(GapPropMcf$Canopy))
    ClumpEffMcf = PAIe/PAI_Mcf

    
    #################################################
    # Generate outputs
    #################################################
    OutputMcf<-data.frame(
      PAIe = PAIe,
      PAI_Mcf = PAI_Mcf,
      ClumpMcf = ClumpMcf,
      ClumpEffMcf = ClumpEffMcf,
      CanopyFrac = mean(GapPropMcf$Canopy),
      SmallGapMcf = mean(GapPropMcf$SmallGap),
      LargeGapMcf = mean(GapPropMcf$LargeGap),
      GapFrac = mean(GapFrac),
      CrownPorMcf = mean(CrownPor),
      CrownCovMcf = mean(CrownCov)
    )
    
    UncertMcf <- data.frame(
      nSamples = sum(GapList$GapClass=='Canopy')
    )
    OutputMcf <- cbind(OutputMcf,UncertMcf)

  }
  
  
  
  if(ClumpMethod=='Alv'|ClumpMethod=='All'){
    
    #Canopy attributes
    CrownCov = 1 - GapPropAlv$LargeGap
    GapFrac = 1 - GapPropAlv$Canopy
    CrownPor = 1 - GapPropAlv$Canopy/CrownCov
    
    # PAI calcs
    PAI_Alv = -mean(CrownCov)*log(mean(CrownPor))/k
    PAIe = -log(mean(GapFrac))/k

    
    # Clumping Calcs
    ClumpAlv = ((1-mean(CrownPor))*log(mean(GapFrac)))/(log(mean(CrownPor))*mean(GapPropAlv$Canopy))
    ClumpEffAlv = PAIe/PAI_Alv

    #################################################
    # Generate outputs
    #################################################
    OutputAlv<-data.frame(
      PAIe = PAIe,
      PAI_Alv = PAI_Alv,
      ClumpAlv = ClumpAlv,
      CanopyFracAlv = mean(GapPropAlv$Canopy),
      SmallGapAlv = mean(GapPropAlv$SmallGap),
      LargeGapAlv = mean(GapPropAlv$LargeGap),
      GapFrac = mean(GapFrac),
      CrownPorAlv = mean(CrownPor),
      CrownCovAlv = mean(CrownCov)
    )
    
    UncertAlv <- data.frame(
      nSamples = sum(GapList$GapClass=='Canopy')
    )
    OutputAlv <- cbind(OutputAlv,UncertAlv)
    
  }
  
  if(ClumpMethod=='LX'|ClumpMethod=='All'){
    
    #Canopy attributes
    GapFrac = 1 - GapPropLX$Canopy

    PAIe = -log(mean(GapFrac))/k
    PAI_LX = -mean(log(GapFrac))/k
    ClumpLX = PAIe/PAI_LX
  
    
    OutputLX<-data.frame(
      PAIe = PAIe,
      PAI_LX = PAI_LX,
      ClumpLX = ClumpLX
    )
    
    if(ClumpMethod=='All'){
      OutputLX<-cbind(OutputLX,data.frame(
        PAI_McfLX = PAI_Mcf/ClumpEffLX,
        ClumpEffMcfLX = PAIe/(PAI_Mcf/ClumpEffLX)
        )
      )
    }
    
    if(ClumpMethod=='All'){
      OutputLX<-cbind(OutputLX,data.frame(
        PAI_AlvLX = PAI_Alv/ClumpEffLX,
        ClumpEffAlvLX = PAIe/(PAI_Alv/ClumpEffLX)
      )
      )
    }
    
    UncertLX <- data.frame(
      nSamples = sum(GapList$GapClass=='Canopy'),
      CanopyFracSD = sd(GapPropLX$Canopy),
      GapFracSD = sd(GapFrac)
    )
    OutputLX<-cbind(OutputLX,UncertLX)
  }
  
  
  
  if(ClumpMethod=='Mcf'){
    output <- OutputMcf
  }else if(ClumpMethod=='Alv'){
    output <- OutputAlv
  }else if(ClumpMethod=='LX'){
    output <- OutputLX
  }else if(ClumpMethod=='All'){
    output <- cbind(OutputMcf,OutputAlv[,-which(colnames(OutputAlv)%in%colnames(OutputMcf))],OutputLX[,-which(colnames(OutputLX)%in%colnames(OutputMcf))])
  }
  
  if(ClumpMethod=='Mcf'|ClumpMethod=='Both'){
    output$GapThMcf <- GapThdMcf
  }
  if(ClumpMethod=='ALv'|ClumpMethod=='Both'){
    output$GapThAlv <- mean(GapPropAlv$GapThdAlv)
  }
  
  if(GapListOut){
    return(list(output,GapList)) # this is updated. fix the rest of the script...
  }else{
    return(output) # this is updated. fix the rest of the script...
  }
  
}


segment<-function(img=img.bw,
                  mask_obj=mask,
                  nseg = 360/15,  # Number of azimuth segments
                  startVZA=0, # start viewing zenith angle
                  endVZA=70, # end viewing zenith angle
                  maxVZA=90, # maximum view angle of the camera
                  nrings=7, # n-1 for the number of zenith segments
                  add2img=F, # should the output be a dataframe or added to the image? default is false 
                  reduce_outputs=T,  # just keep the necessary outputs (discard intermediate steps)
                  x_col='x',  # x column name
                  y_col='y',  # y column name
                  lens='equidistant'
){
  colnames(img)[3]<-'Gap'
  xc=mask_obj$center_x # x coordinate of mask centre
  yc=mask_obj$center_y # y coordinate of mask centre
  rc=mask_obj$radius # radius of mask

  
  # determine whether it's a classified image or a dataframe
  if(class(img)[1]=="SpatRaster"){
    imgdf<-terra::as.data.frame(img, xy=TRUE)
    base::names(imgdf)<-c('x','y',names(img))
  }
  if(class(img)[1]=='data.frame'){
    imgdf<-img
    colnames(imgdf)[which(colnames(imgdf)==x_col)]<-'x'
    colnames(imgdf)[which(colnames(imgdf)==y_col)]<-'y'
  }
  
  ######################################################
  # lens specific bit - taken from F. Chiannuci code.
  ######################################################
  if(lens=='equidistant'){
    VZAbins<-seq(startVZA,endVZA,(endVZA-startVZA)/nrings)
    rset<-round(rc*VZAbins/maxVZA)
  }
  
  #From Pekin & Macfarlane 2009 (doi:10.3390/rs1041298). Note: it gives relative distance y based on relative distance theta:
  if (lens=="FC-E8") {
    VZAbins<-seq(startVZA,endVZA,(endVZA-startVZA)/nrings)
    rset<-round(rc*(1.06*VZAbins/maxVZA+0.00498*(VZAbins/maxVZA)^2-0.0639*(VZAbins/maxVZA)^3))
  }
  
  #From Pekin & Macfarlane 2009. Note: it gives relative distance y based on relative distance theta:
  if (lens=="Nikkor-10.5") {
    VZAbins<-seq(startVZA,endVZA,(endVZA-startVZA)/nrings)
    rset<-round(rc*(1.13*VZAbins/maxVZA+0.00798*(VZAbins/maxVZA)^2-0.138*(VZAbins/maxVZA)^3))
  }
  
  #From Pekin & Macfarlane 2009. Note: it gives relative distance y based on relative distance theta:
  if (lens=="Sigma-4.5") {
    VZAbins<-seq(startVZA,endVZA,(endVZA-startVZA)/nrings)
    rset<-round(rc*(1.12*VZAbins/maxVZA+0.00598*(VZAbins/maxVZA)^2-0.178*(VZAbins/maxVZA)^3))
  }
  
  #From Hemisfer Lens file. Note: it give relative distance (y) based on theta in radians (x)
  if (lens=="FC-E9") {
    VZAbins<-seq(startVZA,endVZA,(endVZA-startVZA)/nrings)
    rset<-round(rc*(0.6427*(VZAbins*pi/180)+0.0346*(VZAbins*pi/180)^2-0.024491*(VZAbins*pi/180)^3))
  }
  
  #From Hemisfer - Note: it give relative distance (y) based on theta in radians (x)
  if (lens=="Sigma-8") {
    VZAbins<-seq(startVZA,endVZA,(endVZA-startVZA)/nrings)
    rset<-round(rc*(0.75276*(VZAbins*pi/180)-0.073937*(VZAbins*pi/180)^2))
  }
  
  #From Hemisfer - Note: it give relative distance (y) based on theta in radians (x)
  if (lens=="Nikkor-8") {
    VZAbins<-seq(startVZA,endVZA,(endVZA-startVZA)/nrings)
    rset<-round(rc*(0.9192*(VZAbins*pi/180)-0.1792*(VZAbins*pi/180)^2-0.000443*(VZAbins*pi/180)^3))
  }
  
  #From Hemisfer - Note: it give relative distance (y) based on theta in radians (x)
  if (lens=="Nikkor-OP10") {
    VZAbins<-seq(startVZA,endVZA,(endVZA-startVZA)/nrings)
    rset<-round(rc*(1.0168*(VZAbins*pi/180)-0.0573*(VZAbins*pi/180)^2-0.1176*(VZAbins*pi/180)^3))
  }
  
  #From Hemisfer - Note: it give relative distance (y) based on theta in radians (x)
  if (lens=="Soligor") {
    VZAbins<-seq(startVZA,endVZA,(endVZA-startVZA)/nrings)
    rset<-round(rc*(0.677923*(VZAbins*pi/180)-0.029481*(VZAbins*pi/180)^2-0.022084*(VZAbins*pi/180)^3+0.041495*(VZAbins*pi/180)^4-0.016644*(VZAbins*pi/180)^5))
  }
  
  #From Hemisfer - Note: it give relative distance (y) based on theta in radians (x)
  if (lens=="Raynox-CF185") {
    VZAbins<-seq(startVZA,endVZA,(endVZA-startVZA)/nrings)
    rset<-round(rc*(0.5982*(VZAbins*pi/180)+0.024459*(VZAbins*pi/180)^2))
  }
  
  #From Hemisfer - Note: it give relative distance (y) based on theta in radians (x)
  if (lens=="CanonEF-8-15") {
    VZAbins<-seq(startVZA,endVZA,(endVZA-startVZA)/nrings)
    rset<-round(rc*(0.88554*(VZAbins*pi/180)-0.29909*(VZAbins*pi/180)^2)+0.089523*(VZAbins*pi/180)^3)
  }
  
  #From Hemisfer - Note: it give relative distance (y) based on theta in radians (x)
  if (lens=="Kodak-SP360") {
    VZAbins<-seq(startVZA,endVZA,(endVZA-startVZA)/nrings)
    rset<-round(rc*(0.42444*(VZAbins*pi/180)+0.23456*(VZAbins*pi/180)^2)-0.063332*(VZAbins*pi/180)^3)
  }
  
  #From Paul Bourke (http://www.paulbourke.net/dome/fisheyecorrect/) - Note: it give relative distance (y) based on theta in radians (x)
  if (lens=="Laowa-4") {
    VZAbins<-seq(startVZA,endVZA,(endVZA-startVZA)/nrings)
    rset<-round(rc*(0.6241*(VZAbins*pi/180)-0.0252*(VZAbins*pi/180)^2)+0.024*(VZAbins*pi/180)^3-0.0183*(VZAbins*pi/180)^4)
  }
  
  #From Paul Bourke (http://www.paulbourke.net/dome/fisheyecorrect/) - Note: it give relative distance (y) based on theta in radians (x)
  if (lens=="Meike-6.5") {
    VZAbins<-seq(startVZA,endVZA,(endVZA-startVZA)/nrings)
    rset<-round(rc*(0.5165*(VZAbins*pi/180)+0.1106*(VZAbins*pi/180)^2)+0.0617*(VZAbins*pi/180)^3-0.0601*(VZAbins*pi/180)^4)
  }
  
  #From Paul Bourke (http://www.paulbourke.net/dome/fisheyecorrect/) - Note: it give relative distance (y) based on theta in radians (x)
  if (lens=="Meike-3.5") {
    VZAbins<-seq(startVZA,endVZA,(endVZA-startVZA)/nrings)
    rset<-round(rc*(0.6475*(VZAbins*pi/180)-0.002*(VZAbins*pi/180)^2)-0.0331*(VZAbins*pi/180)^3-0.00010171*(VZAbins*pi/180)^4)
  }
  
  #From Paul Bourke (http://www.paulbourke.net/dome/fisheyecorrect/) - Note: it give relative distance (y) based on theta in radians (x)
  if (lens=="Entaniya-M12-220") {
    VZAbins<-seq(startVZA,endVZA,(endVZA-startVZA)/nrings)
    rset<-round(rc*(0.5449*(VZAbins*pi/180)-0.0491*(VZAbins*pi/180)^2)+0.0513*(VZAbins*pi/180)^3-0.0166*(VZAbins*pi/180)^4)
  }
  
  #From Paul Bourke (http://www.paulbourke.net/dome/fisheyecorrect/) - Note: it give relative distance (y) based on theta in radians (x)
  if (lens=="Entaniya-M12-250") {
    VZAbins<-seq(startVZA,endVZA,(endVZA-startVZA)/nrings)
    rset<-round(rc*(0.5146*(VZAbins*pi/180)-0.0453*(VZAbins*pi/180)^2)+0.0374*(VZAbins*pi/180)^3-0.013*(VZAbins*pi/180)^4)
  }
  
  #From Paul Bourke (http://www.paulbourke.net/dome/fisheyecorrect/) - Note: it give relative distance (y) based on theta in radians (x)
  if (lens=="Entaniya-HAL-250") {
    VZAbins<-seq(startVZA,endVZA,(endVZA-startVZA)/nrings)
    rset<-round(rc*(0.5187*(VZAbins*pi/180)-0.0351*(VZAbins*pi/180)^2)+0.0333*(VZAbins*pi/180)^3-0.0137*(VZAbins*pi/180)^4)
  }
  
  #From Paul Bourke (http://www.paulbourke.net/dome/fisheyecorrect/) - Note: it give relative distance (y) based on theta in radians (x)
  if (lens=="Entaniya-HAL-200") {
    VZAbins<-seq(startVZA,endVZA,(endVZA-startVZA)/nrings)
    rset<-round(rc*(0.6093*(VZAbins*pi/180)-0.0187*(VZAbins*pi/180)^2)+0.0235*(VZAbins*pi/180)^3-0.0142*(VZAbins*pi/180)^4)
  }
  
  #From Paul Bourke (http://www.paulbourke.net/dome/fisheyecorrect/) - Note: it give relative distance (y) based on theta in radians (x)
  if (lens=="Entaniya-M12-280") {
    VZAbins<-seq(startVZA,endVZA,(endVZA-startVZA)/nrings)
    rset<-round(rc*(0.5229*(VZAbins*pi/180)-0.043*(VZAbins*pi/180)^2)+0.0253*(VZAbins*pi/180)^3-0.0109*(VZAbins*pi/180)^4)
  }
  
  #From Paul Bourke (http://www.paulbourke.net/dome/fisheyecorrect/) - Note: it give relative distance (y) based on theta in radians (x)
  if (lens=="Canon-RF-5.2") {
    VZAbins<-seq(startVZA,endVZA,(endVZA-startVZA)/nrings)
    rset<-round(rc*0.6031*(VZAbins*pi/180))
  }
  
  #From Paul Bourke (http://www.paulbourke.net/dome/fisheyecorrect/) - Note: it give relative distance (y) based on theta in radians (x)
  if (lens=="iZugar-MKX22") {
    VZAbins<-seq(startVZA,endVZA,(endVZA-startVZA)/nrings)
    rset<-round(rc*(0.6622*(VZAbins*pi/180)-0.0163*(VZAbins*pi/180)^2)+0.0029*(VZAbins*pi/180)^3-0.0169*(VZAbins*pi/180)^4)
  }
  
  #From Paul Bourke (http://www.paulbourke.net/dome/fisheyecorrect/) - Note: it give relative distance (y) based on theta in radians (x)
  if (lens=="iZugar-MKX13") {
    VZAbins<-seq(startVZA,endVZA,(endVZA-startVZA)/nrings)
    rset<-round(rc*(0.6208*(VZAbins*pi/180)-0.0101*(VZAbins*pi/180)^2)+0.0134*(VZAbins*pi/180)^3-0.0047*(VZAbins*pi/180)^4)
  }
  
  #From Paul Bourke (http://www.paulbourke.net/dome/fisheyecorrect/) - Note: it give relative distance (y) based on theta in radians (x)
  if (lens=="iZugar-MKX200") {
    VZAbins<-seq(startVZA,endVZA,(endVZA-startVZA)/nrings)
    rset<-round(rc*(0.5898*(VZAbins*pi/180)-0.0133*(VZAbins*pi/180)^2)+0.0193*(VZAbins*pi/180)^3-0.0099*(VZAbins*pi/180)^4)
  }
  
  #From Paul Bourke (http://www.paulbourke.net/dome/fisheyecorrect/) - Note: it give relative distance (y) based on theta in radians (x)
  if (lens=="iZugar-MKX19") {
    VZAbins<-seq(startVZA,endVZA,(endVZA-startVZA)/nrings)
    rset<-round(rc*(0.6355*(VZAbins*pi/180)-0.0259*(VZAbins*pi/180)^2)+0.0292*(VZAbins*pi/180)^3-0.0152*(VZAbins*pi/180)^4)
  }
  
  #From Paul Bourke (http://www.paulbourke.net/dome/fisheyecorrect/) - Note: it give relative distance (y) based on theta in radians (x)
  if (lens=="DZO-VRCA") {
    VZAbins<-seq(startVZA,endVZA,(endVZA-startVZA)/nrings)
    rset<-round(rc*(0.6984*(VZAbins*pi/180)+0.0173*(VZAbins*pi/180)^2)-0.0421*(VZAbins*pi/180)^3-0.008*(VZAbins*pi/180)^4)
  }
  
  #From Paul Bourke (http://www.paulbourke.net/dome/fisheyecorrect/) - Note: it give relative distance (y) based on theta in radians (x)
  if (lens=="Sunex-DSL315") {
    VZAbins<-seq(startVZA,endVZA,(endVZA-startVZA)/nrings)
    rset<-round(rc*(0.748*(VZAbins*pi/180)-0.0272*(VZAbins*pi/180)^2)-0.0032*(VZAbins*pi/180)^3-0.0198*(VZAbins*pi/180)^4)
  }
  
  #From Paul Bourke (http://www.paulbourke.net/dome/fisheyecorrect/) - Note: it give relative distance (y) based on theta in radians (x)
  if (lens=="Sunex-DSL239") {
    VZAbins<-seq(startVZA,endVZA,(endVZA-startVZA)/nrings)
    rset<-round(rc*(0.7129*(VZAbins*pi/180)-0.0239*(VZAbins*pi/180)^2)+0.0069*(VZAbins*pi/180)^3-0.0172*(VZAbins*pi/180)^4)
  }
  
  #From Paul Bourke (http://www.paulbourke.net/dome/fisheyecorrect/) - Note: it give relative distance (y) based on theta in radians (x)
  if (lens=="Sunex-DSL415") {
    VZAbins<-seq(startVZA,endVZA,(endVZA-startVZA)/nrings)
    rset<-round(rc*(0.587*(VZAbins*pi/180)-0.0231*(VZAbins*pi/180)^2)+0.0351*(VZAbins*pi/180)^3-0.0124*(VZAbins*pi/180)^4)
  }
  
  
  #From Paul Bourke (http://www.paulbourke.net/dome/fisheyecorrect/) - Note: it give relative distance (y) based on theta in radians (x)
  if (lens=="Sunex-DSLR01") {
    VZAbins<-seq(startVZA,endVZA,(endVZA-startVZA)/nrings)
    rset<-round(rc*(0.7811*(VZAbins*pi/180)-0.0363*(VZAbins*pi/180)^2)-0.0144*(VZAbins*pi/180)^3-0.0154*(VZAbins*pi/180)^4)
  }
  
  
  #From Paul Bourke (http://www.paulbourke.net/dome/fisheyecorrect/) - Note: it give relative distance (y) based on theta in radians (x)
  if (lens=="Lensagon-BF10M14522S118"){
    VZAbins<-seq(startVZA,endVZA,(endVZA-startVZA)/nrings)
    rset<-round(rc*(0.6391*(VZAbins*pi/180)-0.0252*(VZAbins*pi/180)^2)+0.0381*(VZAbins*pi/180)^3-0.0214*(VZAbins*pi/180)^4)
  }
  
  #From Paul Bourke (http://www.paulbourke.net/dome/fisheyecorrect/) - Note: it give relative distance (y) based on theta in radians (x)
  if (lens=="Lensagon-BF16M220D"){
    VZAbins<-seq(startVZA,endVZA,(endVZA-startVZA)/nrings)
    rset<-round(rc*(.5574*(VZAbins*pi/180)-0.0371*(VZAbins*pi/180)^2)+0.0412*(VZAbins*pi/180)^3-0.0164*(VZAbins*pi/180)^4)
  }
  
  
  #From Paul Bourke (http://www.paulbourke.net/dome/fisheyecorrect/) - Note: it give relative distance (y) based on theta in radians (x)
  if (lens=="Omnitech-ORIFL190-3"){
    VZAbins<-seq(startVZA,endVZA,(endVZA-startVZA)/nrings)
    rset<-round(rc*(0.7418*(VZAbins*pi/180)-0.0309*(VZAbins*pi/180)^2)-0.0022*(VZAbins*pi/180)^3-0.0176*(VZAbins*pi/180)^4)
  }
  
  #From Paul Bourke (http://www.paulbourke.net/dome/fisheyecorrect/) - Note: it give relative distance (y) based on theta in radians (x)
  if (lens=="Aico-ACHIR01028B10M"){
    VZAbins<-seq(startVZA,endVZA,(endVZA-startVZA)/nrings)
    rset<-round(rc*(0.4896*(VZAbins*pi/180)-0.054*(VZAbins*pi/180)^2)+0.0651*(VZAbins*pi/180)^3-0.0208*(VZAbins*pi/180)^4)
  }
  
  #From Paul Bourke (http://www.paulbourke.net/dome/fisheyecorrect/) - Note: it give relative distance (y) based on theta in radians (x)
  if (lens=="Aico-ACHIR01420B9M"){
    VZAbins<-seq(startVZA,endVZA,(endVZA-startVZA)/nrings)
    rset<-round(rc*(0.6194*(VZAbins*pi/180)-0.0094*(VZAbins*pi/180)^2)+0.0126*(VZAbins*pi/180)^3-0.0041*(VZAbins*pi/180)^4)
  }
  
  # From Paul Bourke (http://www.paulbourke.net/dome/fisheyecorrect/) - Note: it give relative distance (y) based on theta in radians (x)
  # if (lens=="Raynox-DCR-CF187"){
  #   VZAbins<-seq(startVZA,endVZA,(endVZA-startVZA)/nrings)
  #   rset<-round(rc*(0.9442*(VZAbins*pi/180)+0.0138*(VZAbins*pi/180)^2)-0.167*(VZAbins*pi/180)^3+0.0213*(VZAbins*pi/180)^4)
  # }
  
  #From Paul Bourke (http://www.paulbourke.net/dome/fisheyecorrect/) - Note: it give relative distance (y) based on theta in radians (x)
  if (lens=="SMTEC-SL-190"){
    VZAbins<-seq(startVZA,endVZA,(endVZA-startVZA)/nrings)
    rset<-round(rc*(0.7091*(VZAbins*pi/180)-0.0123*(VZAbins*pi/180)^2)-0.0013*(VZAbins*pi/180)^3-0.0126*(VZAbins*pi/180)^4)
  }
  
  #From Paul Bourke (http://www.paulbourke.net/dome/fisheyecorrect/) - Note: it give relative distance (y) based on theta in radians (x)
  if (lens=="Evetar-E3307"){
    VZAbins<-seq(startVZA,endVZA,(endVZA-startVZA)/nrings)
    rset<-round(rc*(0.6929*(VZAbins*pi/180)+0.003*(VZAbins*pi/180)^2)+0.0541*(VZAbins*pi/180)^3-0.0548*(VZAbins*pi/180)^4)
  }
  
  #From Paul Bourke (http://www.paulbourke.net/dome/fisheyecorrect/) - Note: it give relative distance (y) based on theta in radians (x)
  if (lens=="Evetar-E3267A"){
    VZAbins<-seq(startVZA,endVZA,(endVZA-startVZA)/nrings)
    rset<-round(rc*(.5171*(VZAbins*pi/180)-0.055*(VZAbins*pi/180)^2)+0.054*(VZAbins*pi/180)^3-0.019*(VZAbins*pi/180)^4)
  }
  
  #From Paul Bourke (http://www.paulbourke.net/dome/fisheyecorrect/) - Note: it give relative distance (y) based on theta in radians (x)
  if (lens=="Evetar-E3279"){
    VZAbins<-seq(startVZA,endVZA,(endVZA-startVZA)/nrings)
    rset<-round(rc*(0.4914*(VZAbins*pi/180)-0.0616*(VZAbins*pi/180)^2)+0.1315*(VZAbins*pi/180)^3-0.0323*(VZAbins*pi/180)^4)
  }
  
  #From Paul Bourke (http://www.paulbourke.net/dome/fisheyecorrect/) - Note: it give relative distance (y) based on theta in radians (x)
  if (lens=="ArduCam-M25156H18"){
    VZAbins<-seq(startVZA,endVZA,(endVZA-startVZA)/nrings)
    rset<-round(rc*(0.6826*(VZAbins*pi/180)-0.0284*(VZAbins*pi/180)^2)+0.0292*(VZAbins*pi/180)^3-0.0185*(VZAbins*pi/180)^4)
  }
  
  #From Paul Bourke (http://www.paulbourke.net/dome/fisheyecorrect/) - Note: it give relative distance (y) based on theta in radians (x)
  if (lens=="Bosch-Flexidome-7000"){
    VZAbins<-seq(startVZA,endVZA,(endVZA-startVZA)/nrings)
    rset<-round(rc*(.7823*(VZAbins*pi/180)-0.0102*(VZAbins*pi/180)^2)-0.0325*(VZAbins*pi/180)^3-0.0127*(VZAbins*pi/180)^4)
  }
  
  if (lens=="equisolid") {
    VZAbins<-seq(startVZA,endVZA,(endVZA-startVZA)/nrings)
    maxVZAmm<-2*sin(maxVZA/2*pi/180)
    rset<-round(rc*2*sin(VZAbins/2*pi/180)/maxVZAmm)
  }
  
  if (lens=="stereographic") {
    VZAbins<-seq(startVZA,endVZA,(endVZA-startVZA)/nrings)
    maxVZAmm<-2*tan(maxVZA/2*pi/180)
    rset<-round(rc*2*tan(VZAbins/2*pi/180)/maxVZAmm)
  }
  
  if (lens=="orthographic") {
    VZAbins<-seq(startVZA,endVZA,(endVZA-startVZA)/nrings)
    maxVZAmm<-2*sin(maxVZA*pi/180)
    rset<-round(rc*2*sin(VZAbins*pi/180)/maxVZAmm)
  }
  ######################################################
  
  
  
  
  
  
  
  
  #create radial distance from centre
  imgdf$dx<-imgdf$x-xc
  imgdf$dy<-imgdf$y-yc
  imgdf$r<-round(sqrt(imgdf$dx^2+imgdf$dy^2)) # radial distance of each pixel to the centre
  
  
  # Initialize a new column 'azimuth' in the data frame 'imgdf' and set all its values to NA
  imgdf$azimuth = NA
  
  # Calculate angles for points in the data frame based on differences in x and y coordinates
  # For points in the first quadrant (dx > 0 and dy >= 0)
  imgdf$azimuth[imgdf$dx > 0 & imgdf$dy >= 0] <- atan(imgdf$dy[imgdf$dx > 0 & imgdf$dy >= 0] / imgdf$dx[imgdf$dx > 0 & imgdf$dy >= 0])
  
  # For points in the second quadrant (dx > 0 and dy < 0)
  imgdf$azimuth[imgdf$dx > 0 & imgdf$dy < 0] <- atan(imgdf$dy[imgdf$dx > 0 & imgdf$dy < 0] / imgdf$dx[imgdf$dx > 0 & imgdf$dy < 0]) + 2 * pi
  
  # For points in the third quadrant (dx < 0)
  imgdf$azimuth[imgdf$dx < 0] <- atan(imgdf$dy[imgdf$dx < 0] / imgdf$dx[imgdf$dx < 0]) + pi
  
  # For points on the positive y-axis (dx == 0 and dy > 0)
  imgdf$azimuth[imgdf$dx == 0 & imgdf$dy > 0] <- pi / 2
  
  # For points on the negative y-axis (dx == 0 and dy < 0)
  imgdf$azimuth[imgdf$dx == 0 & imgdf$dy < 0] <- 3 * pi / 2
  
  # Convert angles from radians to degrees
  imgdf$azimuth <- imgdf$azimuth * 180 / pi
  
  # Create a new column 'ring' by cutting the 'r' column into segments
  imgdf$ring <- cut(imgdf$r, rset, include.lowest = TRUE,
                    labels = seq(startVZA + ((endVZA - startVZA) / 2 / nrings),
                                 endVZA - ((endVZA - startVZA) / 2 / nrings),
                                 (endVZA - startVZA) / nrings))
  
  if(add2img==F){
    # Drop rows with missing values in the 'ring' column
    imgdf <- imgdf[!is.na(imgdf$ring), ]
    
    # Convert 'ring' to numeric
    imgdf$ring <- as.numeric(as.character(imgdf$ring))
    
    # Create a new column 'azimuth_group' by cutting the 'azimuth' column into segments
    imgdf$azimuth_group <- cut(imgdf$azimuth, seq(0, 2 * pi, 2 * pi / nseg) * 180 / pi,
                               include.lowest = TRUE,
                               labels = rep(seq(2 * pi / nseg, 2 * pi, 2 * pi / nseg)) * 180 / pi)
    
    # Convert 'alpha.to' to numeric
    imgdf$azimuth_group <- as.numeric(as.character(imgdf$azimuth_group))
    
    #reduce outputs
    if(reduce_outputs==T){
      idx<-which(colnames(imgdf)%in%c('Gap','r','azimuth','ring','azimuth_group'))
      imgdf<-imgdf[,idx]
    }
    
    return(imgdf)
  }
  if(add2img==T){
    # Convert 'ring' to numeric
    imgdf$ring <- as.numeric(as.character(imgdf$ring))
    
    # Create a new column 'azimuth_group' by cutting the 'azimuth' column into segments
    imgdf$azimuth_group <- cut(imgdf$azimuth, seq(0, 2 * pi, 2 * pi / nseg) * 180 / pi,
                               include.lowest = TRUE,
                               labels = rep(seq(2 * pi / nseg, 2 * pi, 2 * pi / nseg)) * 180 / pi)
    
    # Convert 'alpha.to' to numeric
    imgdf$azimuth_group <- as.numeric(as.character(imgdf$azimuth_group))
    
    # Convert 'azimuth_group' to a raster layer
    idx_1<-which(colnames(imgdf)%in%c('azimuth_group','ring','x','y'))
    img$r <- imgdf$r # this works do for all check where the rgb vals have gone wtf...
    img$ring<-imgdf$ring
    img$azimuth<-imgdf$azimuth
    img$azimuth_group<-imgdf$azimuth_group
    return(img)
  }
}


adj_seg<-function(Seg_img,
                  mask_obj=mask,
                  nseg = 360/15,  # Number of azimuth segments
                  startVZA=55, # start viewing zenith angle
                  endVZA=60, # end viewing zenith angle
                  maxVZA=90, # maximum view angle of the camera
                  nrings=1,
                  lens='equidistant',
                  reduce_outputs=T){ # n-1 for the number of zenith segments){
  
  
  xc=mask_obj$center_x # x coordinate of mask centre
  yc=mask_obj$center_y # y coordinate of mask centre
  rc=mask_obj$radius # radius of mask
  imgdf<-Seg_img

  
  ######################################################
  # lens specific bit - taken from F. Chiannuci code.
  ######################################################
  if(lens=='equidistant'){
    VZAbins<-seq(startVZA,endVZA,(endVZA-startVZA)/nrings)
    rset<-round(rc*VZAbins/maxVZA)
  }
  
  #From Pekin & Macfarlane 2009 (doi:10.3390/rs1041298). Note: it gives relative distance y based on relative distance theta:
  if (lens=="FC-E8") {
    VZAbins<-seq(startVZA,endVZA,(endVZA-startVZA)/nrings)
    rset<-round(rc*(1.06*VZAbins/maxVZA+0.00498*(VZAbins/maxVZA)^2-0.0639*(VZAbins/maxVZA)^3))
  }
  
  #From Pekin & Macfarlane 2009. Note: it gives relative distance y based on relative distance theta:
  if (lens=="Nikkor-10.5") {
    VZAbins<-seq(startVZA,endVZA,(endVZA-startVZA)/nrings)
    rset<-round(rc*(1.13*VZAbins/maxVZA+0.00798*(VZAbins/maxVZA)^2-0.138*(VZAbins/maxVZA)^3))
  }
  
  #From Pekin & Macfarlane 2009. Note: it gives relative distance y based on relative distance theta:
  if (lens=="Sigma-4.5") {
    VZAbins<-seq(startVZA,endVZA,(endVZA-startVZA)/nrings)
    rset<-round(rc*(1.12*VZAbins/maxVZA+0.00598*(VZAbins/maxVZA)^2-0.178*(VZAbins/maxVZA)^3))
  }
  
  #From Hemisfer Lens file. Note: it give relative distance (y) based on theta in radians (x)
  if (lens=="FC-E9") {
    VZAbins<-seq(startVZA,endVZA,(endVZA-startVZA)/nrings)
    rset<-round(rc*(0.6427*(VZAbins*pi/180)+0.0346*(VZAbins*pi/180)^2-0.024491*(VZAbins*pi/180)^3))
  }
  
  #From Hemisfer - Note: it give relative distance (y) based on theta in radians (x)
  if (lens=="Sigma-8") {
    VZAbins<-seq(startVZA,endVZA,(endVZA-startVZA)/nrings)
    rset<-round(rc*(0.75276*(VZAbins*pi/180)-0.073937*(VZAbins*pi/180)^2))
  }
  
  #From Hemisfer - Note: it give relative distance (y) based on theta in radians (x)
  if (lens=="Nikkor-8") {
    VZAbins<-seq(startVZA,endVZA,(endVZA-startVZA)/nrings)
    rset<-round(rc*(0.9192*(VZAbins*pi/180)-0.1792*(VZAbins*pi/180)^2-0.000443*(VZAbins*pi/180)^3))
  }
  
  #From Hemisfer - Note: it give relative distance (y) based on theta in radians (x)
  if (lens=="Nikkor-OP10") {
    VZAbins<-seq(startVZA,endVZA,(endVZA-startVZA)/nrings)
    rset<-round(rc*(1.0168*(VZAbins*pi/180)-0.0573*(VZAbins*pi/180)^2-0.1176*(VZAbins*pi/180)^3))
  }
  
  #From Hemisfer - Note: it give relative distance (y) based on theta in radians (x)
  if (lens=="Soligor") {
    VZAbins<-seq(startVZA,endVZA,(endVZA-startVZA)/nrings)
    rset<-round(rc*(0.677923*(VZAbins*pi/180)-0.029481*(VZAbins*pi/180)^2-0.022084*(VZAbins*pi/180)^3+0.041495*(VZAbins*pi/180)^4-0.016644*(VZAbins*pi/180)^5))
  }
  
  #From Hemisfer - Note: it give relative distance (y) based on theta in radians (x)
  if (lens=="Raynox-CF185") {
    VZAbins<-seq(startVZA,endVZA,(endVZA-startVZA)/nrings)
    rset<-round(rc*(0.5982*(VZAbins*pi/180)+0.024459*(VZAbins*pi/180)^2))
  }
  
  #From Hemisfer - Note: it give relative distance (y) based on theta in radians (x)
  if (lens=="CanonEF-8-15") {
    VZAbins<-seq(startVZA,endVZA,(endVZA-startVZA)/nrings)
    rset<-round(rc*(0.88554*(VZAbins*pi/180)-0.29909*(VZAbins*pi/180)^2)+0.089523*(VZAbins*pi/180)^3)
  }
  
  #From Hemisfer - Note: it give relative distance (y) based on theta in radians (x)
  if (lens=="Kodak-SP360") {
    VZAbins<-seq(startVZA,endVZA,(endVZA-startVZA)/nrings)
    rset<-round(rc*(0.42444*(VZAbins*pi/180)+0.23456*(VZAbins*pi/180)^2)-0.063332*(VZAbins*pi/180)^3)
  }
  
  #From Paul Bourke (http://www.paulbourke.net/dome/fisheyecorrect/) - Note: it give relative distance (y) based on theta in radians (x)
  if (lens=="Laowa-4") {
    VZAbins<-seq(startVZA,endVZA,(endVZA-startVZA)/nrings)
    rset<-round(rc*(0.6241*(VZAbins*pi/180)-0.0252*(VZAbins*pi/180)^2)+0.024*(VZAbins*pi/180)^3-0.0183*(VZAbins*pi/180)^4)
  }
  
  #From Paul Bourke (http://www.paulbourke.net/dome/fisheyecorrect/) - Note: it give relative distance (y) based on theta in radians (x)
  if (lens=="Meike-6.5") {
    VZAbins<-seq(startVZA,endVZA,(endVZA-startVZA)/nrings)
    rset<-round(rc*(0.5165*(VZAbins*pi/180)+0.1106*(VZAbins*pi/180)^2)+0.0617*(VZAbins*pi/180)^3-0.0601*(VZAbins*pi/180)^4)
  }
  
  #From Paul Bourke (http://www.paulbourke.net/dome/fisheyecorrect/) - Note: it give relative distance (y) based on theta in radians (x)
  if (lens=="Meike-3.5") {
    VZAbins<-seq(startVZA,endVZA,(endVZA-startVZA)/nrings)
    rset<-round(rc*(0.6475*(VZAbins*pi/180)-0.002*(VZAbins*pi/180)^2)-0.0331*(VZAbins*pi/180)^3-0.00010171*(VZAbins*pi/180)^4)
  }
  
  #From Paul Bourke (http://www.paulbourke.net/dome/fisheyecorrect/) - Note: it give relative distance (y) based on theta in radians (x)
  if (lens=="Entaniya-M12-220") {
    VZAbins<-seq(startVZA,endVZA,(endVZA-startVZA)/nrings)
    rset<-round(rc*(0.5449*(VZAbins*pi/180)-0.0491*(VZAbins*pi/180)^2)+0.0513*(VZAbins*pi/180)^3-0.0166*(VZAbins*pi/180)^4)
  }
  
  #From Paul Bourke (http://www.paulbourke.net/dome/fisheyecorrect/) - Note: it give relative distance (y) based on theta in radians (x)
  if (lens=="Entaniya-M12-250") {
    VZAbins<-seq(startVZA,endVZA,(endVZA-startVZA)/nrings)
    rset<-round(rc*(0.5146*(VZAbins*pi/180)-0.0453*(VZAbins*pi/180)^2)+0.0374*(VZAbins*pi/180)^3-0.013*(VZAbins*pi/180)^4)
  }
  
  #From Paul Bourke (http://www.paulbourke.net/dome/fisheyecorrect/) - Note: it give relative distance (y) based on theta in radians (x)
  if (lens=="Entaniya-HAL-250") {
    VZAbins<-seq(startVZA,endVZA,(endVZA-startVZA)/nrings)
    rset<-round(rc*(0.5187*(VZAbins*pi/180)-0.0351*(VZAbins*pi/180)^2)+0.0333*(VZAbins*pi/180)^3-0.0137*(VZAbins*pi/180)^4)
  }
  
  #From Paul Bourke (http://www.paulbourke.net/dome/fisheyecorrect/) - Note: it give relative distance (y) based on theta in radians (x)
  if (lens=="Entaniya-HAL-200") {
    VZAbins<-seq(startVZA,endVZA,(endVZA-startVZA)/nrings)
    rset<-round(rc*(0.6093*(VZAbins*pi/180)-0.0187*(VZAbins*pi/180)^2)+0.0235*(VZAbins*pi/180)^3-0.0142*(VZAbins*pi/180)^4)
  }
  
  #From Paul Bourke (http://www.paulbourke.net/dome/fisheyecorrect/) - Note: it give relative distance (y) based on theta in radians (x)
  if (lens=="Entaniya-M12-280") {
    VZAbins<-seq(startVZA,endVZA,(endVZA-startVZA)/nrings)
    rset<-round(rc*(0.5229*(VZAbins*pi/180)-0.043*(VZAbins*pi/180)^2)+0.0253*(VZAbins*pi/180)^3-0.0109*(VZAbins*pi/180)^4)
  }
  
  #From Paul Bourke (http://www.paulbourke.net/dome/fisheyecorrect/) - Note: it give relative distance (y) based on theta in radians (x)
  if (lens=="Canon-RF-5.2") {
    VZAbins<-seq(startVZA,endVZA,(endVZA-startVZA)/nrings)
    rset<-round(rc*0.6031*(VZAbins*pi/180))
  }
  
  #From Paul Bourke (http://www.paulbourke.net/dome/fisheyecorrect/) - Note: it give relative distance (y) based on theta in radians (x)
  if (lens=="iZugar-MKX22") {
    VZAbins<-seq(startVZA,endVZA,(endVZA-startVZA)/nrings)
    rset<-round(rc*(0.6622*(VZAbins*pi/180)-0.0163*(VZAbins*pi/180)^2)+0.0029*(VZAbins*pi/180)^3-0.0169*(VZAbins*pi/180)^4)
  }
  
  #From Paul Bourke (http://www.paulbourke.net/dome/fisheyecorrect/) - Note: it give relative distance (y) based on theta in radians (x)
  if (lens=="iZugar-MKX13") {
    VZAbins<-seq(startVZA,endVZA,(endVZA-startVZA)/nrings)
    rset<-round(rc*(0.6208*(VZAbins*pi/180)-0.0101*(VZAbins*pi/180)^2)+0.0134*(VZAbins*pi/180)^3-0.0047*(VZAbins*pi/180)^4)
  }
  
  #From Paul Bourke (http://www.paulbourke.net/dome/fisheyecorrect/) - Note: it give relative distance (y) based on theta in radians (x)
  if (lens=="iZugar-MKX200") {
    VZAbins<-seq(startVZA,endVZA,(endVZA-startVZA)/nrings)
    rset<-round(rc*(0.5898*(VZAbins*pi/180)-0.0133*(VZAbins*pi/180)^2)+0.0193*(VZAbins*pi/180)^3-0.0099*(VZAbins*pi/180)^4)
  }
  
  #From Paul Bourke (http://www.paulbourke.net/dome/fisheyecorrect/) - Note: it give relative distance (y) based on theta in radians (x)
  if (lens=="iZugar-MKX19") {
    VZAbins<-seq(startVZA,endVZA,(endVZA-startVZA)/nrings)
    rset<-round(rc*(0.6355*(VZAbins*pi/180)-0.0259*(VZAbins*pi/180)^2)+0.0292*(VZAbins*pi/180)^3-0.0152*(VZAbins*pi/180)^4)
  }
  
  #From Paul Bourke (http://www.paulbourke.net/dome/fisheyecorrect/) - Note: it give relative distance (y) based on theta in radians (x)
  if (lens=="DZO-VRCA") {
    VZAbins<-seq(startVZA,endVZA,(endVZA-startVZA)/nrings)
    rset<-round(rc*(0.6984*(VZAbins*pi/180)+0.0173*(VZAbins*pi/180)^2)-0.0421*(VZAbins*pi/180)^3-0.008*(VZAbins*pi/180)^4)
  }
  
  #From Paul Bourke (http://www.paulbourke.net/dome/fisheyecorrect/) - Note: it give relative distance (y) based on theta in radians (x)
  if (lens=="Sunex-DSL315") {
    VZAbins<-seq(startVZA,endVZA,(endVZA-startVZA)/nrings)
    rset<-round(rc*(0.748*(VZAbins*pi/180)-0.0272*(VZAbins*pi/180)^2)-0.0032*(VZAbins*pi/180)^3-0.0198*(VZAbins*pi/180)^4)
  }
  
  #From Paul Bourke (http://www.paulbourke.net/dome/fisheyecorrect/) - Note: it give relative distance (y) based on theta in radians (x)
  if (lens=="Sunex-DSL239") {
    VZAbins<-seq(startVZA,endVZA,(endVZA-startVZA)/nrings)
    rset<-round(rc*(0.7129*(VZAbins*pi/180)-0.0239*(VZAbins*pi/180)^2)+0.0069*(VZAbins*pi/180)^3-0.0172*(VZAbins*pi/180)^4)
  }
  
  #From Paul Bourke (http://www.paulbourke.net/dome/fisheyecorrect/) - Note: it give relative distance (y) based on theta in radians (x)
  if (lens=="Sunex-DSL415") {
    VZAbins<-seq(startVZA,endVZA,(endVZA-startVZA)/nrings)
    rset<-round(rc*(0.587*(VZAbins*pi/180)-0.0231*(VZAbins*pi/180)^2)+0.0351*(VZAbins*pi/180)^3-0.0124*(VZAbins*pi/180)^4)
  }
  
  
  #From Paul Bourke (http://www.paulbourke.net/dome/fisheyecorrect/) - Note: it give relative distance (y) based on theta in radians (x)
  if (lens=="Sunex-DSLR01") {
    VZAbins<-seq(startVZA,endVZA,(endVZA-startVZA)/nrings)
    rset<-round(rc*(0.7811*(VZAbins*pi/180)-0.0363*(VZAbins*pi/180)^2)-0.0144*(VZAbins*pi/180)^3-0.0154*(VZAbins*pi/180)^4)
  }
  
  
  #From Paul Bourke (http://www.paulbourke.net/dome/fisheyecorrect/) - Note: it give relative distance (y) based on theta in radians (x)
  if (lens=="Lensagon-BF10M14522S118"){
    VZAbins<-seq(startVZA,endVZA,(endVZA-startVZA)/nrings)
    rset<-round(rc*(0.6391*(VZAbins*pi/180)-0.0252*(VZAbins*pi/180)^2)+0.0381*(VZAbins*pi/180)^3-0.0214*(VZAbins*pi/180)^4)
  }
  
  #From Paul Bourke (http://www.paulbourke.net/dome/fisheyecorrect/) - Note: it give relative distance (y) based on theta in radians (x)
  if (lens=="Lensagon-BF16M220D"){
    VZAbins<-seq(startVZA,endVZA,(endVZA-startVZA)/nrings)
    rset<-round(rc*(.5574*(VZAbins*pi/180)-0.0371*(VZAbins*pi/180)^2)+0.0412*(VZAbins*pi/180)^3-0.0164*(VZAbins*pi/180)^4)
  }
  
  
  #From Paul Bourke (http://www.paulbourke.net/dome/fisheyecorrect/) - Note: it give relative distance (y) based on theta in radians (x)
  if (lens=="Omnitech-ORIFL190-3"){
    VZAbins<-seq(startVZA,endVZA,(endVZA-startVZA)/nrings)
    rset<-round(rc*(0.7418*(VZAbins*pi/180)-0.0309*(VZAbins*pi/180)^2)-0.0022*(VZAbins*pi/180)^3-0.0176*(VZAbins*pi/180)^4)
  }
  
  #From Paul Bourke (http://www.paulbourke.net/dome/fisheyecorrect/) - Note: it give relative distance (y) based on theta in radians (x)
  if (lens=="Aico-ACHIR01028B10M"){
    VZAbins<-seq(startVZA,endVZA,(endVZA-startVZA)/nrings)
    rset<-round(rc*(0.4896*(VZAbins*pi/180)-0.054*(VZAbins*pi/180)^2)+0.0651*(VZAbins*pi/180)^3-0.0208*(VZAbins*pi/180)^4)
  }
  
  #From Paul Bourke (http://www.paulbourke.net/dome/fisheyecorrect/) - Note: it give relative distance (y) based on theta in radians (x)
  if (lens=="Aico-ACHIR01420B9M"){
    VZAbins<-seq(startVZA,endVZA,(endVZA-startVZA)/nrings)
    rset<-round(rc*(0.6194*(VZAbins*pi/180)-0.0094*(VZAbins*pi/180)^2)+0.0126*(VZAbins*pi/180)^3-0.0041*(VZAbins*pi/180)^4)
  }
  
  # From Paul Bourke (http://www.paulbourke.net/dome/fisheyecorrect/) - Note: it give relative distance (y) based on theta in radians (x)
  # if (lens=="Raynox-DCR-CF187"){
  #   VZAbins<-seq(startVZA,endVZA,(endVZA-startVZA)/nrings)
  #   rset<-round(rc*(0.9442*(VZAbins*pi/180)+0.0138*(VZAbins*pi/180)^2)-0.167*(VZAbins*pi/180)^3+0.0213*(VZAbins*pi/180)^4)
  # }
  
  #From Paul Bourke (http://www.paulbourke.net/dome/fisheyecorrect/) - Note: it give relative distance (y) based on theta in radians (x)
  if (lens=="SMTEC-SL-190"){
    VZAbins<-seq(startVZA,endVZA,(endVZA-startVZA)/nrings)
    rset<-round(rc*(0.7091*(VZAbins*pi/180)-0.0123*(VZAbins*pi/180)^2)-0.0013*(VZAbins*pi/180)^3-0.0126*(VZAbins*pi/180)^4)
  }
  
  #From Paul Bourke (http://www.paulbourke.net/dome/fisheyecorrect/) - Note: it give relative distance (y) based on theta in radians (x)
  if (lens=="Evetar-E3307"){
    VZAbins<-seq(startVZA,endVZA,(endVZA-startVZA)/nrings)
    rset<-round(rc*(0.6929*(VZAbins*pi/180)+0.003*(VZAbins*pi/180)^2)+0.0541*(VZAbins*pi/180)^3-0.0548*(VZAbins*pi/180)^4)
  }
  
  #From Paul Bourke (http://www.paulbourke.net/dome/fisheyecorrect/) - Note: it give relative distance (y) based on theta in radians (x)
  if (lens=="Evetar-E3267A"){
    VZAbins<-seq(startVZA,endVZA,(endVZA-startVZA)/nrings)
    rset<-round(rc*(.5171*(VZAbins*pi/180)-0.055*(VZAbins*pi/180)^2)+0.054*(VZAbins*pi/180)^3-0.019*(VZAbins*pi/180)^4)
  }
  
  #From Paul Bourke (http://www.paulbourke.net/dome/fisheyecorrect/) - Note: it give relative distance (y) based on theta in radians (x)
  if (lens=="Evetar-E3279"){
    VZAbins<-seq(startVZA,endVZA,(endVZA-startVZA)/nrings)
    rset<-round(rc*(0.4914*(VZAbins*pi/180)-0.0616*(VZAbins*pi/180)^2)+0.1315*(VZAbins*pi/180)^3-0.0323*(VZAbins*pi/180)^4)
  }
  
  #From Paul Bourke (http://www.paulbourke.net/dome/fisheyecorrect/) - Note: it give relative distance (y) based on theta in radians (x)
  if (lens=="ArduCam-M25156H18"){
    VZAbins<-seq(startVZA,endVZA,(endVZA-startVZA)/nrings)
    rset<-round(rc*(0.6826*(VZAbins*pi/180)-0.0284*(VZAbins*pi/180)^2)+0.0292*(VZAbins*pi/180)^3-0.0185*(VZAbins*pi/180)^4)
  }
  
  #From Paul Bourke (http://www.paulbourke.net/dome/fisheyecorrect/) - Note: it give relative distance (y) based on theta in radians (x)
  if (lens=="Bosch-Flexidome-7000"){
    VZAbins<-seq(startVZA,endVZA,(endVZA-startVZA)/nrings)
    rset<-round(rc*(.7823*(VZAbins*pi/180)-0.0102*(VZAbins*pi/180)^2)-0.0325*(VZAbins*pi/180)^3-0.0127*(VZAbins*pi/180)^4)
  }
  
  if (lens=="equisolid") {
    VZAbins<-seq(startVZA,endVZA,(endVZA-startVZA)/nrings)
    maxVZAmm<-2*sin(maxVZA/2*pi/180)
    rset<-round(rc*2*sin(VZAbins/2*pi/180)/maxVZAmm)
  }
  
  if (lens=="stereographic") {
    VZAbins<-seq(startVZA,endVZA,(endVZA-startVZA)/nrings)
    maxVZAmm<-2*tan(maxVZA/2*pi/180)
    rset<-round(rc*2*tan(VZAbins/2*pi/180)/maxVZAmm)
  }
  
  if (lens=="orthographic") {
    VZAbins<-seq(startVZA,endVZA,(endVZA-startVZA)/nrings)
    maxVZAmm<-2*sin(maxVZA*pi/180)
    rset<-round(rc*2*sin(VZAbins*pi/180)/maxVZAmm)
  }
  ######################################################
  
  
  # Create a new column 'ring' by cutting the 'r' column into segments
  imgdf$ring <- cut(imgdf$r, rset, include.lowest = TRUE,
                    labels = seq(startVZA + ((endVZA - startVZA) / 2 / nrings),
                                 endVZA - ((endVZA - startVZA) / 2 / nrings),
                                 (endVZA - startVZA) / nrings))
  imgdf <- imgdf[!is.na(imgdf$ring), ]
  
  # Convert 'ring' to numeric
  imgdf$ring <- as.numeric(as.character(imgdf$ring))
  
  # Create a new column 'azimuth_group' by cutting the 'azimuth' column into segments
  imgdf$azimuth_group <- cut(imgdf$azimuth, seq(0, 2 * pi, 2 * pi / nseg) * 180 / pi,
                             include.lowest = TRUE,
                             labels = rep(seq(2 * pi / nseg, 2 * pi, 2 * pi / nseg)) * 180 / pi)
  
  # Convert 'alpha.to' to numeric
  imgdf$azimuth_group <- as.numeric(as.character(imgdf$azimuth_group))
  
  #reduce outputs
  if(reduce_outputs==T){
    idx<-which(colnames(imgdf)%in%c('Gap','ring','azimuth_group'))
    imgdf<-imgdf[,idx]
  }
  
  return(imgdf)
}
deWitt_G<-function(){
  deWit<-data.frame(
    mu=c(2.77,1.172,3.326,.433,1,1.101),
    nu=c(1.172,2.77,3.326,.433,1,1.93),
    LAD=c('Planophile','Erectophile','Plagiophile','Extremophile','Uniform','Spherical')
  )
  
  out<-list()
  for(idx in 1:nrow(deWit)){
    mu <- deWit$mu[idx]
    nu <- deWit$nu[idx]
    
    rad=pi/180
    
    realG<-NULL
    zen_angles<-seq(0.5,89.5,by=1)
    
    
    ft=matrix(NA,ncol=length(zen_angles))
    for (i in 1:length(zen_angles)){
      z<-zen_angles[i]
      t=2*z*rad/pi
      a=(1-t)^(mu-1)
      b=t^(nu-1)
      ft[i]=2/pi*(1/beta(mu,nu))*a*b# probability density function
    }
    frel=ft/sum(ft)
    
    Gi <- NULL
    for(i in 1:length(zen_angles)){
      z<-zen_angles[i]
      
      for(j in 1:length(zen_angles)){
        z2<-zen_angles[j]
        (cot_i <- 1/tan(z*rad))
        (cot_fl <- 1/tan(z2*rad))
        
        if(abs(cot_i*cot_fl)>1) {
          A = cos(z*rad)*cos(z2*rad)
        } else {
          A = cos(z*rad)*cos(z2*rad)*(1+(2/pi)*((tan(acos(cot_i*cot_fl)))-acos(cot_i*cot_fl)))
        }
        phi=A*frel[i]
        
        Gi <- c(Gi, phi)
      }
    }
    
    Gmat <- matrix(Gi, ncol=90)
    Gfun<-apply(Gmat, 1, sum)# G-function
    realG<-rbind(realG,round(Gfun,3))
    
    out[[idx]]<-data.frame(zenith=seq(0.5,89.5,by=1),pdf=round(as.numeric(frel),4),G=as.numeric(realG),type=deWit$LAD[idx])
    names(out)[idx]<-deWit$LAD[idx]
  }
  return(out)
}

getPAINils_fromImg <- function(fcs,Code,
                       sat_PAI=8, # user defined inputs for saturated cells
                       G_funct='Spherical',Z_max=70,G_list=deWitt_G()){   # grouping. can be either an image or site code..

  G_focus<-G_list[[which(names(G_list)==G_funct)]]
  # get the G function
  ring_bins<-unique(fcs$ring)
  G_focus$ring_group<-NA
  for(k in 1:length(G_focus$zenith)){
    G_focus$ring_group[k]<- ring_bins[which((ring_bins - G_focus$zenith[k])^2 == min((ring_bins - G_focus$zenith[k])^2))]
  }
  G_0<-G_focus[G_focus$zenith<Z_max,]
  G_0$zenith_rad<-(G_0$zenith* pi) / 180
  
  PAI_by_ring<-list()
  
  # calculate weights and add a gap to saturated pixels based on saturated PAI vals
  zenith<-(as.numeric(unique(ring_bins))* pi) / 180
  zenith<-zenith[order(zenith)]
  
  max_Pgap<-list()
  for(rngs in 1:length(zenith)){
    max_Pgap[[rngs]] <- exp(sum(-sat_PAI*G_0$G[G_0$ring_group==ring_bins[rngs]])/ sum(cos(G_0$zenith_rad[G_0$ring_group==ring_bins[rngs]])))
  }
  
  # get gap fractions
  zen_group<-as.numeric(unique(as.character(fcs$ring)))
  zen_group<-zen_group[order(zen_group)]
  for (rngs in 1:length(unique(fcs$ring))){
    idx_a<-which(fcs$ring==zen_group[rngs])
    gap_by_azm<-list()
    fcs_2<-fcs[idx_a,]
    
    
    
    # calculate Pgap; 
    # also calculate Pgap for LX clumping that gives numbers to saturated azum segs
    satP_azm<-list()
    for (azm in 1:length(unique(fcs_2$azimuth_group))){
      azumith_seg<-fcs_2[fcs_2$azimuth_group==unique(fcs_2$azimuth_group)[azm],]
      
      PgapT<-sum(azumith_seg$Gap)/nrow(azumith_seg)
      PgapT_LX<-PgapT
      if(PgapT==0&nrow(azumith_seg)>0){ 
        PgapT_LX<-max_Pgap[[rngs]]
        satP_azm[[azm]]<-1
      }else{satP_azm[[azm]]<-0}
      
      gap_by_azm[[azm]]<-as.data.frame(cbind(PgapT,PgapT_LX))
    }
    gap_by_azm<-do.call(rbind,gap_by_azm)
    satP_azm<-sum(do.call(rbind,satP_azm))
    
    if(satP_azm>0){gap_by_azm$PgapT[gap_by_azm$PgapT==0]<-max_Pgap[[rngs]]}
    
    
    #effective LAI, PAI and WAI calcs
    PAIe_rng<-sum(-log(mean(na.omit(gap_by_azm$PgapT)))*cos(G_0$zenith_rad[G_0$ring_group==ring_bins[rngs]]))/(0.5*nrow(G_0[G_0$ring_group==ring_bins[rngs],]))

    
    
    #if any of them are 0s the above code will shit itself so this gives it a 0 value for LAI - meaning the LAI is very low...
    if(!is.finite(PAIe_rng)&mean(na.omit(gap_by_azm$PgapT))==0){PAIe_rng<-0}

    #LX clumping adjusted LAI calcs
    PAI_rng<-sum(mean(-log(na.omit(gap_by_azm$PgapT_LX)))*cos(G_0$zenith_rad[G_0$ring_group==ring_bins[rngs]]))/(0.5*nrow(G_0[G_0$ring_group==ring_bins[rngs],]))
    
    # edit saturated segments
    if(mean(na.omit(gap_by_azm$PgapT))==1){
      Psat_ring<-1
    }else{Psat_ring<-0}
    

    PAI_LXG1_rng=NA
    PAI_LXG2_rng=NA
    
    
    PAI_by_ring[[rngs]]<-as.data.frame(cbind(PAIe_rng,
                                             PAI_rng,
                                             PAI_LXG1_rng,
                                             PAI_LXG2_rng,
                                             Psat_ring,
                                             satP_azm
    ))
  }
  PAI_by_ring<-do.call(rbind,PAI_by_ring)
  PAI_by_img_temp<-as.data.frame(t(as.matrix(colMeans(PAI_by_ring[,1:4]))))
  PAI_by_img_temp<-cbind(PAI_by_img_temp,as.data.frame(t(as.matrix(colSums(PAI_by_ring[,5:6])))))
  PAI_by_img_temp$Code<-as.vector(Code)

  colnames(PAI_by_img_temp)[1:4]<-c('PAIe',
                                    'PAI_LX',
                                    'PAI_LXG1',
                                    'PAI_LXG2')
  return(PAI_by_img_temp)
}
