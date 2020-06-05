
# Clean up R
rm(list=ls())

# Load libraries
library(plyr)
library(countrycode)

# Working directory
wd<-"/Users/alistairsenior/Dropbox/Work-Home Share/Human lifetables and Nutrition"

# Read in all lifetables
setwd(paste0(wd, "/rawdata/FAO"))
new<-read.csv("NormalizedFBS_new.csv")
old<-read.csv("NormalizedFBS_historic.csv")
head(new)
head(old)

# combine the two
combined<-rbind(old, new)

# We only need a fraction of these data - total supply of food, fat protein and alcohol, with carbs taken as the difference - note protein and fat must be converted from g to kcal later
group1<-which(combined$Item == "Grand Total" & combined$Element == "Food supply (kcal/capita/day)")
group2<-which(combined$Item == "Grand Total" & combined$Element == "Protein supply quantity (g/capita/day)")
group3<-which(combined$Item == "Grand Total" & combined$Element == "Fat supply quantity (g/capita/day)")
group4<-which(combined$Item == "Alcoholic Beverages" & combined$Element == "Food supply (kcal/capita/day)")

combined$Food_group<-NA
combined$Food_group[group1]<-"Total.kcal"
combined$Food_group[group2]<-"Protein.g"
combined$Food_group[group3]<-"Fat.g"
combined$Food_group[group4]<-"Alcohol.kcal"

filtered<-combined[c(group1, group2, group3, group4),]

# Cut out some of the other stuff that is not needed 
filtered<-filtered[,c("Area", "Year", "Food_group", "Value")]

# Reshape to be wide format
wide<-reshape(filtered, direction="wide", idvar=c("Area", "Year"), v.names="Value", timevar="Food_group")
wide$Protein.kcal<-wide$Value.Protein.g * 4
wide$Fat.kcal<-wide$Value.Fat.g * 9
wide$Alcohol.kcal<-wide$Value.Alcohol.kcal
# Note for the UAE Alcohol is listed as NA as the country is dry - lets set it to 0
wide$Alcohol.kcal[which(is.na(wide$Alcohol.kcal) == T)]<-0
wide$Carbo.kcal<-wide$Value.Total.kcal - wide$Protein.kcal - wide$Fat.kcal - wide$Alcohol.kcal

# Clean out some more
wide<-wide[,c("Area", "Year", "Protein.kcal", "Carbo.kcal", "Fat.kcal", "Alcohol.kcal")]

# We have to split and duplicate the belgium-luxembourg data (which was pooled until 2000) as we have seperatre lifetables for these
tag<-which(wide$Area == "Belgium-Luxembourg")
BL_data<-wide[tag,]
wide$Area[tag]<-"Belgium"
BL_data$Area<-"Luxembourg"
wide<-rbind(wide, BL_data)

# We have to split and duplicate the serbia-montenegro data (which was pooled in mid 2000s) as we have seperatre lifetables for these
tag<-which(wide$Area == "Serbia and Montenegro")
SM_data<-wide[tag,]
wide$Area[tag]<-"Serbia"
SM_data$Area<-"Montenegro"
wide<-rbind(wide, SM_data)

# We have to split out the data on Yugoslav SFR as we have seperate lifetables for bosnia, serbia, slovenia, croatia, macedonia and montenegro 
tag<-which(wide$Area == "Yugoslav SFR")
YG_data<-wide[tag,]
YG_data$Area<-"Bosnia & Herzegovina"
wide<-rbind(wide, YG_data)
YG_data$Area<-"Serbia"
wide<-rbind(wide, YG_data)
YG_data$Area<-"Slovenia"
wide<-rbind(wide, YG_data)
YG_data$Area<-"Macedonia"
wide<-rbind(wide, YG_data)
YG_data$Area<-"Montenegro"
wide<-rbind(wide, YG_data)
YG_data$Area<-"Croatia"
wide<-rbind(wide, YG_data)

# We have to split out the data on Czechoslovakia as we have seperate lifetables for Czechia, and slovakia 
tag<-which(wide$Area == "Czechoslovakia")
CZ_data<-wide[tag,]
CZ_data$Area<-"Czechia"
wide<-rbind(wide, CZ_data)
CZ_data$Area<-"Slovakia"
wide<-rbind(wide, CZ_data)

# We have to split out the data on USSR as we have seperate lifetables for Russia, Lithuania, Tajikistan and Estonia 
tag<-which(wide$Area == "USSR")
US_data<-wide[tag,]
US_data$Area<-"Lithuania"
wide<-rbind(wide, US_data)
US_data$Area<-"Estonia"
wide<-rbind(wide, US_data)
US_data$Area<-"Russian Federation"
wide<-rbind(wide, US_data)
US_data$Area<-"Tajikistan"
wide<-rbind(wide, US_data)

# Sort by country
wide<-wide[order(wide$Area),]

# Convert to iso3c country names as used in the LT data
wide$Country<-countrycode(wide$Area, origin="country.name", destination="iso3c")

# Add a country code for YUG - seems to be missing from the function and we have one or two LTs for Yugoslavia as a whole
tag<-which(wide$Area == "Yugoslav SFR")
wide$Country[tag]<-"YUG"

# Add a country code for CSK - seems to be missing from the function, and we have LTs for Czechoslovakia as a whole
tag<-which(wide$Area == "Czechoslovakia")
wide$Country[tag]<-"CSK"

# Add a country code for SUN - seems to be missing from the function, and we have LTs for USSR as a whole
tag<-which(wide$Area == "USSR")
wide$Country[tag]<-"SUN"

# Save the clean table as a .csv file
setwd(paste0(wd, "/cleandata"))
write.table(wide, file="Clean_FB.csv", sep=",", row.names=F, col.names=names(wide))

