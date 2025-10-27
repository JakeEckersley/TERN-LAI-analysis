# apply clumping to df
source('hemispheric_PAI_functions.R')
library(terra)
library(hemispheR)
library(gtools)


PAI_df<-read.csv('J:/TERN_dhp_LAI_QC_complete_251012update.csv')

# define columns to summarise
colsub <- c("PAIe_57","PAI_LX_57","Psat_ring_57","satP_azm_57","CI_57",
          "PAIe_70","PAI_LX_70","Psat_ring_70","satP_azm_70","CI_70")
coltots <- c("Psat_ring_57","satP_azm_57","Psat_ring_70","satP_azm_70")

# split into site and date
PAI_df<-split(PAI_df,list(PAI_df$Site,PAI_df$Plot,PAI_df$Date))

# filter out dates without samples
PAI_df <- PAI_df[unlist(lapply(PAI_df,function(x)nrow(x)>0))]

summarisePAI<-  function(x) {
  header <- x[1,c('Site','Plot','Date')]
  header$n_imgs <- nrow(x)
  x <- x[x[,'QC'] == 1 & x[,'Duplicate']  == 1, ]
  header$n_retained <- nrow(x)
  results<-list()
  for(i in 1:length(colsub)){
    results[[i]]<-data.frame(mean(x[,colsub[i]], na.rm=TRUE), sd(x[,colsub[i]], na.rm=TRUE), sum(x[,colsub[i]], na.rm=TRUE))
    colnames(results[[i]])<-c(paste0(colsub[i],'_mean'),paste0(colsub[i],'_sd'),paste0(colsub[i],'_sum'))
    if(!colsub[i]%in%coltots){
      results[[i]]<-results[[i]][,-3]
    }
  }
  results<-do.call(cbind,results)
  results<-cbind(header,results)
  return(results)
}

PAI_sum<-do.call(rbind,lapply(PAI_df,summarisePAI))
rownames(PAI_sum)<-NULL

out <- as.data.frame(aggregate(PAI_df[, colsub],
                 by = list(PAI_df$Site,PAI_df$Plot,PAI_df$Date),
                 FUN = function(x) c(mean = mean(x, na.rm=TRUE), sd = sd(x, na.rm=TRUE), sum = sum(x, na.rm=TRUE))))



######################################################
# Summary stats for the Gaps_70 df
######################################################
Gaps_70<-readRDS('J:/TERN_dhp_RawGapFrac_70.rds')
Gaps_70<-split(Gaps_70,list(Gaps_70$Site,Gaps_70$Plot,Gaps_70$Date))
LAI_dat<-list()
for(i in 1:length(Gaps_70)){
  if(nrow(Gaps_70[[i]])==0){
    LAI_dat[[i]]<-NULL
  }else{
    Code<-names(Gaps_70)[i]
    LAI_dat[[i]]<-getPAINils_fromImg(Gaps_70[[i]],Code)
  }
}
LAI_70<-do.call(rbind,LAI_dat[-which(unlist(lapply(LAI_dat,is.null)))])
code_split<-strsplit(LAI_70$Code,'.',fixed=T)
LAI_70$CI_70<-(LAI_70$PAIe/LAI_70$PAI_LX)
LAI_70$Site<-unlist(lapply(code_split,function(x)x<-x[1]))
LAI_70$Plot<-unlist(lapply(code_split,function(x)x<-x[2]))
LAI_70$Date<-unlist(lapply(code_split,function(x)x<-x[3]))
colnames(LAI_70)[1:2]<-paste0(colnames(LAI_70)[1:2],'_70')


Gaps_57<-readRDS('J:/TERN_dhp_RawGapFrac_57.rds')
Gaps_57<-split(Gaps_57,list(Gaps_57$Site,Gaps_57$Plot,Gaps_57$Date))
LAI_dat<-list()
for(i in 1:length(Gaps_57)){
  if(nrow(Gaps_57[[i]])==0){
    LAI_dat[[i]]<-NULL
  }else{
    Code<-names(Gaps_57)[i]
    LAI_dat[[i]]<-getPAINils_fromImg(Gaps_57[[i]],Code)
  }
}
LAI_57<-do.call(rbind,LAI_dat[-which(unlist(lapply(LAI_dat,is.null)))])
code_split<-strsplit(LAI_57$Code,'.',fixed=T)
LAI_57$CI_57<-(LAI_57$PAIe/LAI_57$PAI_LX)
LAI_57$Site<-unlist(lapply(code_split,function(x)x<-x[1]))
LAI_57$Plot<-unlist(lapply(code_split,function(x)x<-x[2]))
LAI_57$Date<-unlist(lapply(code_split,function(x)x<-x[3]))
colnames(LAI_57)[1:2]<-paste0(colnames(LAI_57)[1:2],'_57')

LAI_combined<-cbind(LAI_70[,c(1:2,8)],LAI_57[,c(1:2,8)],LAI_70[,-c(1:4,7,8)])

write.csv(LAI_combined,'J:/TERN_dhp_LAI_QC_complete_SiteLevelClumping_251006update.csv',na='',row.names = F)


####### bind those together.






#############################################
## do the same for the dcp images
#############################################
dcp_df<-read.csv('J:/TERN_dcp_LAI_QC_complete_251015update.csv')
dcp_gaps<-split(dcp_df,list(dcp_df$Site,dcp_df$Plot,dcp_df$Date))
LAI_dat<-list()
for(i in 1:length(dcp_gaps)){
  if(nrow(dcp_gaps[[i]])==0){
    LAI_dat[[i]]<-NULL
  }else{
    Code<-names(dcp_gaps)[i]
    code_split<-strsplit(Code,'.',fixed=T)
    LAI_dat[[i]]<-data.frame(Site=code_split[1],Plot=code_split[2],Date=code_split[3])
    k <- 0.5 # used for euc forest. not equal to spherical LAD! g(angle)/cos(angle)
    cc <- 1-dcp_gaps[[i]]$Large_gaps # this is canopy cover. check that large gaps is already a proportion
    Gapfrac <- 1-dcp_gaps[[i]]$Canopy # is this foliar cover?? I think so...
    CP = dcp_gaps[[i]]$Canopy
    Clump<-(1-(1-mean(ff/fc)))*log(1-mean(ff))/log(1-mean(ff/fc))*mean(ff)
    PAI_cc <- -mean(fc)*log(1-mean(ff/fc))/k # 1-ff/fc is crown porosity
    PAI_lx <- mean(log(1-fc))/k # 
    PAI_clx <- -mean(fc)*mean(log(1-ff/fc))/k
    PAIe<-log(1-mean(fc))/k
    PAI_cc_check<-Clump*log(1-mean(ff/fc))/k
    ClumpEff <- PAIe/PAI
  }
}
DCP_agg<-do.call(rbind,LAI_dat[-which(unlist(lapply(LAI_dat,is.null)))])
code_split<-strsplit(DCP_agg$Code,'.',fixed=T)
DCP_agg$Site<-unlist(lapply(code_split,function(x)x<-x[1]))
DCP_agg$Plot<-unlist(lapply(code_split,function(x)x<-x[2]))
DCP_agg$Date<-unlist(lapply(code_split,function(x)x<-x[3]))
DCP_agg$


LAI_57$CI_57<-(LAI_57$PAIe/LAI_57$PAI_LX)
