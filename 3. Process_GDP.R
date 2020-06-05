
# Clean up R
rm(list=ls())

# Load libraries
library(plyr)
library(readxl)

# Working directory
wd<-"/Users/alistairsenior/Dropbox/Human lifetables and Nutrition"
#wd<-"/Users/alistairsenior/Dropbox (Sydney Uni)/Work-Home Share/Human lifetables and Nutrition"

# Read in all lifetables
setwd(paste0(wd, "/rawdata/Maddison"))

# Read in the data
data<-read_excel("mpd2018.xlsx", sheet=2)
head(data)

# Throw out any data with missing GDP - note we will have a strategy for dealing with these IF they have correlates - this just means we will not end upo with any data with all missing covariates
data<-data[-which(is.na(data$cgdppc) == T),]

# Throw out any data for country groupings that no longer exist 
data<-data[-which(data$countrycode == "SUN" & data$year > 1991),]
data<-data[-which(data$countrycode == "YUG" & data$year > 1991),]
data<-data[-which(data$countrycode == "CSK" & data$year > 1991),]

# Save the clean table as a .csv file
setwd(paste0(wd, "/cleandata"))
write.table(data, file="Clean_GDP_Mad.csv", sep=",", row.names=F, col.names=names(data))

# Read in all lifetables
setwd(paste0(wd, "/rawdata/WB"))

data_i<-read.table("GDP_perCapita.csv", sep=",", skip=4, header=T)
data_i<-data_i[,-c(3:4, dim(data_i)[2])]

# Get the years
years<-names(data_i)[-c(1:2)]
years<-unlist(lapply(strsplit(years, "X"), "[[", 2))

# Tranform to wide format
long_i<-reshape(data_i, direction="long", idvar=c("Country.Name", "Country.Code"), varying=list(3:dim(data_i)[2]), v.names="X")
long_i$Year<-as.numeric(years[long_i$time])
long_i$predictor<-"GDP_WB"
long<-long_i	

# Remove the 'time' column
long<-long[,-3]

# Reshape to be wide format
wide<-reshape(long, direction="wide", idvar=c("Country.Name", "Country.Code", "Year"), v.names="X", timevar="predictor")
wide<-wide[order(wide$Country.Name),]
names(wide)[-c(1:3)]<-unlist(lapply(strsplit(names(wide)[-c(1:3)], ".", fixed=T), "[[", 2))

# For some of the older LTs we have Czeck data coded as CSK Lets use the CZE data for that 
tag<-which(wide$Country.Code == "CZE" & wide$Year < 1991)
CZ_data<-wide[tag,]
CZ_data$Country.Name<-"Czechoslovakia"
CZ_data$Country.Code<-"CSK"
wide<-rbind(wide, CZ_data)


# For some of the older LTs we have Russian data coded as SUN Lets use the RUS data for that 
tag<-which(wide$Country.Code == "RUS" & wide$Year < 1991)
RU_data<-wide[tag,]
RU_data$Country.Name<-"USSR"
RU_data$Country.Code<-"SUN"
wide<-rbind(wide, RU_data)

# As above drop any missing data
wide<-wide[-which(is.na(wide$GDP_WB) == T),]

# Save the clean table as a .csv file
setwd(paste0(wd, "/cleandata"))
write.table(wide, file="Clean_GDP_WB.csv", sep=",", row.names=F, col.names=names(wide))

