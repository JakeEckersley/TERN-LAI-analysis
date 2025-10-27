# create a file I can label
file_locations<-c("J:/Tern images")
im_folders<-list.dirs(file_locations,recursive = F)
im_folders<-list.dirs(im_folders)
im_folders<-im_folders[grepl('/images',im_folders)]
im_files<-list.files(im_folders,full.names = T)


file_details<-do.call(rbind,strsplit(tools::file_path_sans_ext(basename(im_files)),'-'))
colnames(file_details)<-c('Site','Image type','Plot','Date','Image orientation','Image code')
im_details<-lapply(strsplit(dirname(im_files),'/',fixed=T),function(x){t(as.matrix(unlist(strsplit(x[3],"_",fixed=T))))})
im_details[!which(lapply(im_details,ncol)==5)]
im_details<-do.call(rbind,im_details)[,-1]
colnames(im_details)<-c('Method','Camera','Folder suffix','Download code')

im_details<-cbind(file_details,im_details,
                  data.frame(file.dirs=dirname(im_files),
                            file.path=im_files))
im_details$AutoThd<-NA
im_details$ManualThd<-NA
im_details$QC<-NA
im_details$Duplicate<-NA
im_details$Comments<-NA

# write a csv
write.csv(im_details,paste0(file_locations,'/Tern image details.csv'),row.names = F,na='')
