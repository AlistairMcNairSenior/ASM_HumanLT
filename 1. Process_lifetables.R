
# Clean up R
rm(list=ls())

# Load libraries
library(plyr)

# Working directory
wd<-"/Users/alistairsenior/Dropbox (Sydney Uni)/Human lifetables and Nutrition" # Work iMac
setwd(wd)

# Read in the standard from Wilmoth et al 2012
standard<-read.csv("clean_data/wilmoth_standard.csv")

# Read in all lifetables
setwd(paste0(wd, "/rawdata/HMD"))
lifetables<-read.csv("hld.csv")
head(lifetables)

# For lifetables with a year-range calculate the mid point
lifetables$YearMid<-floor((lifetables$Year1 + lifetables$Year2) / 2)

# Cull out those ages not in the standard
lifetables<-lifetables[-which(is.na(match(lifetables$Age, standard$Age)) == T),]
head(lifetables)

# Cull out any lifetable data not coming from the whole country
lifetables<-lifetables[-which(lifetables$Region != 0),]
lifetables<-lifetables[-which(lifetables$Residence != 0),]
lifetables<-lifetables[-which(lifetables$Ethnicity != 0),]
lifetables<-lifetables[-which(lifetables$SocDem != 0),]

# Add a lifetable ID for the year, country, sex and reference
lifetables$LT_ID<-paste(lifetables$Country, lifetables$YearMid, lifetables$Sex, lifetables$Ref.ID, sep="_")

# For any one reference for a country x year x sex we want the minimal TypeLT - (1, 2, or 4)
IDS<-unique(lifetables$LT_ID)
for(i in 1:length(IDS)){
	sub<-lifetables[which(lifetables$LT_ID == IDS[i]),]
	sub<-sub[which(sub$TypeLT == min(sub$TypeLT)),]
	if(i == 1){
		lifetables_trim<-sub
	}else{
		lifetables_trim<-rbind(lifetables_trim, sub)
	}
}

# Do we still have several 'versions' for anyone ID
check<-ddply(lifetables_trim, .(LT_ID), summarise, length(unique(Version)))
check[which(check$..1 != 1),"LT_ID"]

# Yes there is this one from hongkong in 1973 with both sexes - the data seem identical so lets just use version1
tag<-which(is.na(match(lifetables_trim$LT_ID, check[which(check$..1 != 1),"LT_ID"])) == F)
hk_1973<-lifetables_trim[tag,]
lifetables_trim<-lifetables_trim[-tag,]
lifetables_trim<-rbind(lifetables_trim, hk_1973[which(hk_1973$Version == 1),])

# Do we still have several 'versions' for anyone reference
check<-ddply(lifetables_trim, .(LT_ID), summarise, length(unique(Version)))
check[which(check$..1 != 1),"LT_ID"]
# No

# Now within each lifetable lets normalise the l.x. to its cohort size
lifetables_trim$l.x.prop<-0
IDS<-unique(lifetables_trim$LT_ID)
for(i in 1:length(IDS)){
	tag<-which(lifetables_trim$LT_ID == IDS[i])
	sub<-lifetables_trim[tag,]
	lifetables_trim$l.x.prop[tag]<-sub$l.x. / sub$l.x.[1]
}

# Lets check that the first entry is 1 for all lifetables
check<-ddply(lifetables_trim, .(LT_ID), summarise, l.x.prop[1])
check[which(check$..1 != 1),"LT_ID"]
# OK, looking good

# Now need to check whether there are multiple references per year x country x sex
# Add a lifetable ID for the year, country, sex
lifetables_trim$LT_ID2<-paste(lifetables_trim$Country, lifetables_trim$YearMid, lifetables_trim$Sex, sep="_")
check<-ddply(lifetables_trim, .(LT_ID2), summarise, length(unique(LT_ID)))
check[which(check$..1 != 1),"LT_ID2"]
# Yes there are lots, so lets average for those

IDS<-unique(lifetables_trim$LT_ID2)
for(i in 1:length(IDS)){
	tag<-which(lifetables_trim$LT_ID2 == IDS[i])
	sub<-lifetables_trim[tag,]
	sub_mu<-ddply(sub, .(Age), summarise, l.x.=mean(l.x.prop))
	sub_mu$Country<-sub$Country[1]
	sub_mu$Year<-sub$YearMid[1]
	sub_mu$Sex<-sub$Sex[1]
	if(i == 1){
		lifetables_clean<-sub_mu
	}else{
		lifetables_clean <-rbind(lifetables_clean, sub_mu)
	}
}

# Drop any LT data prior to 1961 - when the food balance sheets start
lifetables_clean<-lifetables_clean[-which(lifetables_clean$Year < 1961),]

# Drop any LT data after 2017 - when the food balance sheets end
lifetables_clean<-lifetables_clean[-which(lifetables_clean$Year > 2017),]

# Save the clean table as a .csv file
setwd(paste0(wd, "/clean_data"))
write.table(lifetables_clean, file="Clean_LT.csv", sep=",", row.names=F, col.names=names(lifetables_clean))

