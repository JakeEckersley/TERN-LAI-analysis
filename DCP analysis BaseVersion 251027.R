########### upward facing dcp photos
source('hemispheric_PAI_functions_251020update.R')
library(terra)
library(hemispheR)
library(ggplot2)
library(patchwork)
library(ggpubr)
library(dplyr)
library(gtools)
library(autothresholdr)
library(parallel)
library(future.apply)
library(progressr)


# for tern - create a script that looks in a folder and just runs through all images
# create an initial run script (i.e., no QC, just the list of images and paths then the initial run)
# add a blank 'manual threshold' column
# change the QC script to check for a manual threshold and to use that as img_thd where > 0 else PAIe == 0

QC<-read.csv("J:/TERN LAI final doc/TERN_dcp_LAI_QC_complete.csv")
QC$file.path<-gsub('D:/',"J:/",QC$file.path)
QC$file.dirs<-gsub('D:/',"J:/",QC$file.dirs)
img_output_path<-"J:/BinarisedImagesBase/DCP upward/"
dir.create(img_output_path,showWarnings = F)
colnames(QC)[which(colnames(QC)=="img_thd")]<-"AutoThd"

# this is just a temporary fix...
QC$Comments[QC$AutoThd == 501] <- 'No plant material'
QC$AutoThd[QC$AutoThd == 501] <- 0
QC$ManualThd <- NA
QC$ManualThd[QC$thresh_type=='Manual']<-QC$AutoThd[QC$thresh_type=='Manual']
QC$AutoThd[QC$thresh_type=='Manual']<-NA
QC<-QC[QC$path_id!='clpm-lai-core1ha-20170502-upward-DSC04285',]

#DCP_analysis(QC[20,])

DCP_analysis<-function(ImgDat,img_outputs=img_output_path,PrintImg=T,GLA=F){
  
  start<-Sys.time()
  
  #####################################
  # read image
  #####################################
  dcp<-ImgDat$file.path
  Code<-tools::file_path_sans_ext(basename(dcp))
  rast_im<-read_B_and_gamma_correct(dcp,rgb=T,stretch_im = T)
  
  
  #####################################
  # Check QC
  #####################################  
  qc_fail<-ifelse(ImgDat$QC==1,F,T)
  if(is.na(qc_fail)){
    ImgDat$QC = NA
    ImgDat$Comments='QC needed'
  }
  
  #####################################
  # Check img isn't weird
  #####################################  
  if(!nlyr(rast_im)%in%c(1,3)){
    message('odd value for ',Code)
    outputs <- data.frame(path_id=Code,PAIe=NA,img_thd=NA)
    ImgDat <- cbind(ImgDat[c('Site','Image.type','Plot','Date')],outputs,ImgDat[,c('Image.orientation.2','Image.code','Method','Camera','Download.code','file.path','QC','Comments','thresh_type','Duplicate')])
    ImgDat$Comments <- 'Odd number of layers in image'
    ImgDat$QC <- 0
    return(list(ImgDat,NULL))
  }
  
  
  ##############################################
  # select the blue band or create GLA img
  ##############################################
  if(GLA){
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
  }else{
    rast_im<-rast_im[[which(names(rast_im)=='Blue')]]
  }
  
  #####################################
  # Determine threshold
  #####################################
  img.mat <- matrix(terra::values(rast_im),nrow = nrow(rast_im),byrow=T)
  
  if(is.na(ImgDat$AutoThd)){ #check if there is a number in thresh_QC
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
  
  #####################################
  # Binarise the image
  #####################################
  img.bw <- terra::classify(rast_im, rbind(c(-Inf,img_thd,0),c( img_thd,Inf,1)))
  base::names(img.bw) <- base::names(rast_im)
  mk<- is.na(img.bw)
  img.bw[mk]<-NA
  
  #############################
  # Generate outputs
  #############################
  if(img_thd!=0){
    outputs <- getPAIdcp(img.bw,ImgID=Code,GapListOut=T)
    DataDetails<-data.frame(
      th_channel = "Blue",
      th_method = 'Otsu'
    )
    GapList <- outputs[[2]]
    outputs <- cbind(outputs[[1]],DataDetails)
  }else{
    outputs <- data.frame(Code=Code,
                          PAIe=0,
                          k=0.5,      
                          th_channel = "Blue",
                          th_method = 'Otsu')
    GapList<-getPAIdcp(img.bw,ImgID=Code,GapListOut=T)[[2]]
  }
  
  # update QC csv
  ImgDat <- cbind(ImgDat[,c('Site','Image.type','Plot','Date')],outputs,ImgDat[,c('Image.orientation.2','Image.code','Method','Camera','Download.code','file.path','thresh_type','AutoThd','ManualThd','QC','Duplicate','Comments')])
  GapList$Code <- paste(ImgDat[,c('Site','Plot','Date')],collapse='_')
  
  
  #############################
  # Plot QC img
  #############################
  if(PrintImg){
    png(paste0(img_outputs,Code,'.png'), width = 1600, height = 1200, res = 150)
    
    # Set 2x2 layout (use 3 panels)
    par(mfrow = c(2, 2))
    par(mar = c(4, 4, 2, 1))  # reduce margins a bit
    
    # re-read rgb if it is downward
    if(GLA){
      rast_im<-read_B_and_gamma_correct(dcp,rgb=T,stretch_im = T)
    }
    
    # --- plot #1 rgb
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
         main = paste0('thresh = ',img_thd),
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
    
    ## --- 4. Blank
    plot.new()
    text(0.5, 0.5, " ")
    
    dev.off()
  }
  
  rm(img.bw,vals,rast_im)
  
  end<-Sys.time()
  
  message(Code,' ',end-start)
  
  return(list(ImgDat,GapList))
  
  terra::tmpFiles(remove = TRUE)
  
}


###################################
# Apply the function
###################################
# Split QC into chunks (for each survey)
split_QC <- split(QC, list(QC$Site,QC$Plot,QC$Date))
split_QC<-split_QC[unlist(lapply(split_QC,function(x)nrow(x)>0))]

#################################################
###### run in parallel
#################################################
# # the following parallelises the application of the dcp function. the same can be achieved by running the lapply function in the 'sequential processing' section

handlers(global = TRUE)   # enable global progress reporting



FastApply <- function(X) {
  # define progressor of total number of rows across all chunks
  total_rows <- sum(sapply(X, nrow))
  p <- progressr::progressor(steps = total_rows)
  
  future_lapply(
    X,
    function(df_chunk) {
      results <- vector("list", nrow(df_chunk))
      
      for (i in seq_len(nrow(df_chunk))) {
        results[[i]] <- DCP_analysis(df_chunk[i, ])
        p(message = sprintf("Row %d of %d", i, nrow(df_chunk)))
      }
      return(results)
      
    },
    future.seed = TRUE,
    future.packages = c('terra','hemispheR','ggplot2','patchwork','ggpubr','gtools')
  )
}



plan(multisession, workers = 10)

  LAI_chunks <- FastApply(split_QC)

plan(sequential)

LAI_df<-do.call(smartbind,
                lapply(LAI_chunks,function(x){
                  x <- do.call(smartbind,lapply(x,function(y)y<-y[[1]]))
                  }
                )
  )

saveRDS(LAI_chunks,'J:/TERN_dcp_LAI_QC_complete_251012update.rds')

#########################################
## Run in sequence
#########################################
# LAI_chunks <- lapply(split_QC,function(df_chunk) {
#     results <- vector("list", nrow(df_chunk))
#     for (i in seq_len(nrow(df_chunk))) {
#       results[[i]] <- DCP_analysis(df_chunk[i, ])
#     }
#     return(results)
#   })
# 
# LAI_df<-do.call(smartbind,lapply(LAI_chunks,function(x){
#       x <- do.call(smartbind,lapply(x,function(y)y<-y[[1]]))
#     })
#   )

#########################################################
# Add site level uncertainty
#########################################################
# data
LAI_chunks<-readRDS('J:/TERN_dcp_LAI_QC_complete_251012update.rds')
LAI_df<-read.csv('J:/TERN_dcp_LAI_QC_complete_251021update.csv')

# update the dodgy ones belinda pointed out...
ReDo<-read.csv("J:/Tern images/re-do subset.csv")

LAI_df<-LAI_df[-which(LAI_df$file.path%in%gsub('D:/','J:/',ReDo$file.path)),]
write.csv(LAI_df,'J:/TERN_dcp_LAI_QC_complete_251021update.csv',na='',row.names = F)


s57<-list()
for(i in 1:nrow(ReDo)){
  s57[[i]]<-paste(ReDo[i,c('Site','Plot','Date')],collapse = '.')
}
LAI_chunks<-LAI_chunks[-which(names(LAI_chunks)%in%unique(unlist(s57)))]


LAI_df<-split(LAI_df,list(LAI_df$Site,LAI_df$Plot,LAI_df$Date))
LAI_df<-LAI_df[which(unlist(lapply(LAI_df,function(x)nrow(x)>0)))]


# get site level vals with site level clumping
LAI_gaps <- lapply(LAI_chunks, function(chunk) {
  
  # loop through sub-elements
  out <- lapply(seq_along(chunk), function(i) {
    df <- chunk[[i]][[2]]           # extract the gaps part
    cbind(ImgIdx = i, df)           # add img index
  })
  
  do.call(rbind, out)               # bind into one dataframe
})


names(LAI_gaps)

############################
# Apply EQ accross sites
############################
outList<-list()
for(i in 1:length(LAI_gaps)){
  message(paste0(round(i/length(LAI_gaps)*100,2),' %'))
  outList[[i]]<-getPAIdcp.df(LAI_gaps[[i]])
  cbind(data.frame("Code" = names(LAI_gaps)[i]),outList[[i]])
}

outList<-do.call(rbind,outList)

###################################################################
# add the summary statistics from the individual img analysis
##################################################################
sumCols<-c('PAIe','PAI_Mcf','ClumpMcf','ClumpEffMcf',
           'CanopyFrac','SmallGapMcf','LargeGapMcf','GapFrac',
           'CrownPorMcf','CrownCovMcf','PAI_Alv',
           'ClumpAlv','CanopyFracAlv','SmallGapAlv',
           'LargeGapAlv','CrownPorAlv',
           'CrownCovAlv','GapThAlv')
LAI_list<-lapply(LAI_df,function(chunk){
    SumRows<-list()
    for(i in 1:length(sumCols)){
      if(sumCols[i]%in%colnames(chunk)){
        SumRows[[i]]<-data.frame(mean(chunk[,sumCols[i]],na.rm=T),sd(chunk[,sumCols[i]],na.rm=T))
        colnames(SumRows[[i]])<-c(paste0(sumCols[i],"_IndMean"),paste0(sumCols[i],"_IndSD"))
      }
    }
    return(do.call(cbind,SumRows))
  })
LAI_list<-do.call(rbind,LAI_list)

results<-cbind(outList,LAI_list)
write.csv(results,'J:/TERN_dcp_SiteAggretateLAI_QC_complete_251025update.csv')

























GetImgBw <- function(PAIrow){
  start<-Sys.time()
  dcp<-PAIrow$file.path
  
  # read in img
  rast_im <- read_B_and_gamma_correct(dcp,rgb=T,stretch_im = T)
  
  # select the blue band
  rast_im <- rast_im[[which(names(rast_im)=='Blue')]]
  
  #####################################
  # Binarise the image
  #####################################
  img.mat <- matrix(terra::values(rast_im),nrow = nrow(rast_im),byrow=T)
  
  #check if there is a number in ManualThd
  if(!is.na(PAIrow$ManualThd)){
    img_thd <- PAIrow$ManualThd
  }else if(is.na(PAIrow$AutoThd)){ 
    img_thd <- unlist(autothresholdr::auto_thresh(round(img.mat), method='Otsu',ignore_na=TRUE))[1]
    PAIrow$AutoThd <- img_thd
    PAIrow$Comments<-'QC needed'
    PAIrow$QC<-NA
    PAIrow$thresh_type<-'Auto'
  }else{
    img_thd <- PAIrow$AutoThd
  }
  
  img.bw <- terra::classify(rast_im, rbind(c(-Inf,img_thd,0),c( img_thd,Inf,1)))
  base::names(img.bw) <- base::names(rast_im)
  mk<- is.na(img.bw)
  img.bw[mk]<-NA
  return(img.bw)
}

#SitePlotDateDf = dcp_gaps[[1]]

SitePlotDatePAI<-function(SitePlotDateDf) {
  init = T
  for (j in 1:nrow(SitePlotDateDf)) {
    p(message = sprintf("Row %d of %d", i, nrow(SitePlotDateDf)))
    if(init){
      img.bw <- GetImgBw(SitePlotDateDf[j, ])
      init = F
    }else{
      img.bw <- c(img.bw,GetImgBw(SitePlotDateDf[j, ]))
    }
  }
  
  outputs <- getPAIdcp(img.bw,
                       ImgID=as.vector(apply(SitePlotDateDf[, c('Site','Plot','Date')],1,function(x) paste(x, collapse = "_")))
                       )
  DataDetails<-data.frame(
    th_channel = "Blue",
    th_method = 'Otsu'
  )
  outputs <- cbind(outputs,DataDetails)
  return(outputs)
}

FastApplyPAI <- function(x) {
  
  # define progressor of total number of rows across all chunks
  total_rows <- sum(sapply(x, nrow))
  p <- progressr::progressor(steps = total_rows)
  
  SitePlotDatePAI<-function(SitePlotDateDf) {
    init = T
    for (j in 1:nrow(SitePlotDateDf)) {
      p(message = sprintf("Row %d of %d", i, nrow(SitePlotDateDf)))
      if(init){
        img.bw <- GetImgBw(SitePlotDateDf[j, ])
        init = F
      }else{
        img.bw <- c(img.bw,GetImgBw(SitePlotDateDf[j, ]))
      }
    }
    
    outputs <- getPAIdcp(img.bw,
                         ImgID=as.vector(apply(SitePlotDateDf[, c('Site','Plot','Date')],1,function(x) paste(x, collapse = "_")))
    )
    DataDetails<-data.frame(
      th_channel = "Blue",
      th_method = 'Otsu'
    )
    outputs <- cbind(outputs,DataDetails)
    return(outputs)
  }
  
  future_lapply(
    x,
    SitePlotDatePAI,
    future.seed = TRUE,
    future.packages = c('terra','ggplot2')
  )
}

plan(multisession, workers = 20)

  SiteDF <- FastApplySitePAI(dcp_gaps)

plan(sequential)

