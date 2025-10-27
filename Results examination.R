library(ggplot2)
library("exifr")
library(gtools)
# Define path
img_dir <- "E:/Tern images"

# List all files recursively
all_files <- list.files(path = img_dir, recursive = TRUE, full.names = TRUE)
length(all_files)
lai_data_dcp <- read.csv("E:/TERN LAI final doc/TERN_dcp_LAI_QC_complete.csv")
lai_data_57 <- read.csv("E:/TERN_dcp_57_LAI_QC_updated250926.csv")
lai_data_dhp <- read.csv("E:/TERN LAI final doc/TERN_dhp_LAI_QC_complete.csv")
lai_data_57_old <- read.csv("E:/TERN LAI final doc/TERN_dcp_57_LAI_QC_complete.csv")

(
  nrow(lai_data_dcp[lai_data_dcp$thresh_type%in%c('manual','Manual'),])+
  nrow(lai_data_dhp[lai_data_dhp$thresh_type%in%c('manual','Manual'),])+
  nrow(lai_data_57[lai_data_57$thresh_type%in%c('manual','Manual'),])+
  nrow(lai_data_dcp[lai_data_dcp$QC!=1,])+
  nrow(lai_data_dhp[lai_data_dhp$QC!=1,])+
  nrow(lai_data_57[lai_data_57$QC!=1,])
  
)/(
  nrow(lai_data_dcp)+
  nrow(lai_data_dhp)+
  nrow(lai_data_57)
)

lai_data_dcp$Date <- as.Date(as.character(lai_data_dcp$Date), format = "%Y%m%d")


# filter out if QC is bad
lai_data_dcp<-lai_data_dcp[lai_data_dcp$QC==1,]

# check LAI
unique(lai_data_dcp_dhp$Site)
lai_focus<-lai_data_dcp[lai_data_dcp$Site=="alic"&lai_data_dcp$Plot=="tower_footprint",]

# Aggregate PAIe by Date (mean per day)
daily_mean <- aggregate(PAIe ~ Date, data = lai_focus, FUN = mean, na.rm = TRUE)

# Plot
ggplot(daily_mean, aes(x = Date, y = PAIe)) +
  geom_point(color = "blue", size = 2) +
  geom_line(color = "darkblue") +
  theme_bw() +
  labs(x = "Date", y = "Mean PAIe", title = "Mean PAIe")



#######################
# check with updated 57 degrees
#######################
library(gtools)
lai_data_57<-smartbind(lai_data_57_old,lai_data_57)
lai_data_57$Date <- as.Date(as.character(lai_data_57$Date), format = "%Y%m%d")
for(i in 1:nrow(lai_data_57)){
  if(is.na(lai_data_57$QC[i])){
    lai_data_57$QC[i]<-lai_data_57$QC_1[i]
  }
}
# filter out if QC is bad
lai_data_57<-lai_data_57[lai_data_57$QC==1,]

# check LAI
unique(lai_data_57$Plot)
lai_focus<-lai_data_57

# Aggregate PAIe by Date (mean per day)
daily_mean <- aggregate(PAIe ~ Date, data = lai_focus, FUN = mean, na.rm = TRUE)

# Plot
ggplot(daily_mean, aes(x = Date, y = PAIe)) +
  geom_point(color = "blue", size = 2) +
  geom_line(color = "darkblue") +
  theme_bw() +
  labs(x = "Date", y = "Mean PAIe", title = "Alic Tower footprint")




#######################
# check whroo..
# filter out if QC is bad
lai_data_dhp$Date <- as.Date(as.character(lai_data_dhp$Date), format = "%Y%m%d")

# check LAI
unique(lai_data_dhp$Site)
lai_focus<-lai_data_dhp[lai_data_dhp$Site=="vicd_whroo"&lai_data_dhp$Plot=="core1ha",]
lai_focus<-lai_focus[lai_focus$QC==1,]

# Aggregate PAIe by Date (mean per day)
daily_mean <- aggregate(PAIe_70 ~ Date, data = lai_focus, FUN = mean, na.rm = TRUE)

# Plot
ggplot(daily_mean, aes(x = Date, y = PAIe_70)) +
  geom_point(color = "blue", size = 2) +
  geom_line(color = "darkblue") +
  theme_bw() +
  labs(x = "Date", y = "Mean PAIe", title = "Mean PAIe")


