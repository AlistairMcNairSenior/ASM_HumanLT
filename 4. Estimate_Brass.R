
rm(list=ls())

# Working directory
wd<-"/Users/alistairsenior/Dropbox (Sydney Uni)/Human lifetables and Nutrition" # Work iMac
#wd<-"/Users/asenior/Dropbox (Sydney Uni)/Human lifetables and Nutrition"
setwd(wd)

# Load libraries
library(arm)
library(plyr)
library(ggplot2)
library(lmerTest)
library(countrycode)
library(rnaturalearth)
library(mgcv)
library(gridExtra)
library(MortalityLaws)
library(Cairo)
library(doSNOW)
library(readr)
library(dplyr)
library(lme4)
library(mice)
library(miceadds)
library(mitml)
library(GGally)

# The header functions written for this analysis
source("scripts/0. Header_Functions.R")

# Read in the standard from Wilmoth 2012
standard<-read.csv("clean_data/wilmoth_standard.csv")

# Read in the cleaned LT data
LT_data<-read.csv("clean_data/Clean_LT.csv")
LT_data$Sex<-as.factor(LT_data$Sex)

# Add the LT_ID - country x year
LT_data$LT_ID<-paste(LT_data$Country, LT_data$Year, sep="_")

# Ditch the year 0 and 110 -  these have lx 1 and 0 (inf and -inf) in the predictors
LT_data<-LT_data[-which(LT_data$Age == 0 ),]
LT_data<-LT_data[-which(LT_data$Age == 110 ),]

# Ditch any countries for which ages 5 and 60 are not present - we need these years to model properly
LT_data$sex_ID<-paste0(LT_data$LT_ID, "_", LT_data$Sex)
check<-ddply(LT_data, .(sex_ID), summarise, age_check=sum(is.na(match(Age, c(5, 60))) == F))
drop<-check[which(check$age_check != 2),"sex_ID"]
LT_data<-LT_data[-which(is.na(match(LT_data$sex_ID, drop)) == F),]

# Check whether any of the data have l.x. == 1 or 0
which(LT_data$l.x. == 1)
which(LT_data$l.x. == 0)
# There are some for 0
LT_data<-LT_data[-which(LT_data$l.x. == 0),]

# Logit transform the y
LT_data$y<-logit(LT_data$l.x.)

# Prepare the predictors
# Note this is done in a sex-specific way
LT_data$v1<-NA
LT_data$v2<-NA
LT_data$v3<-NA

# For males
LT_males<-LT_data[which(LT_data$Sex == 1),]
tag<-match(LT_males$Age, standard$Age)

# V1 is the logit survvial in the age class in the standard pattern
LT_males$v1<-logit(standard$lx_male[tag])

# V2 is gamma_x * (1 - logit(age5)ratio)
age5_ratio<-LT_males[which(LT_males$Age == 5), c("y", "v1", "LT_ID")] 
age5_ratio$ratio<-(1 - (age5_ratio$y / age5_ratio$v1))
LT_males$v2<-standard$gamma_male[tag] * age5_ratio$ratio[match(LT_males$LT_ID, age5_ratio$LT_ID)]

# V3 is theta_x * (1 - logit(age60)ratio)
age60_ratio<-LT_males[which(LT_males$Age == 60), c("y", "v1", "LT_ID")] 
age60_ratio$ratio<-(1 - (age60_ratio$y / age60_ratio$v1))
LT_males$v3<-standard$theta_male[tag] * age60_ratio$ratio[match(LT_males$LT_ID, age60_ratio$LT_ID)]

# For females
LT_females<-LT_data[which(LT_data$Sex == 2),]
tag<-match(LT_females$Age, standard$Age)

# V1 is the logit survival in the age class in the standard pattern
LT_females$v1<-logit(standard$lx_female[tag])

# V2 is gamma_x * (1 - logit(age5)ratio)
age5_ratio<-LT_females[which(LT_females$Age == 5), c("y", "v1", "LT_ID")] 
age5_ratio$ratio<-(1 - (age5_ratio$y / age5_ratio$v1))
LT_females$v2<-standard$gamma_female[tag] * age5_ratio$ratio[match(LT_females$LT_ID, age5_ratio$LT_ID)]

# V3 is theta_x * (1 - logit(age60)ratio)
age60_ratio<-LT_females[which(LT_females$Age == 60), c("y", "v1", "LT_ID")] 
age60_ratio$ratio<-(1 - (age60_ratio$y / age60_ratio$v1))
LT_females$v3<-standard$theta_female[tag] * age60_ratio$ratio[match(LT_females$LT_ID, age60_ratio$LT_ID)]

# Rebind together the data now we have sex-specific predictors
LT_data<-rbind(LT_males, LT_females)

# Run the random-regression on the LT_data
model<-lmer(y ~ v1 + Sex + v1:Sex + (1 + v1 + Sex + v1:Sex|LT_ID) + offset(v2) + offset(v3), data=LT_data, control=lmerControl(optimizer="bobyqa"))
summary(model)

# Save the coefficients from the model
table_name<-"tables/Table_S8_(LMM).csv"
write.table(Sys.time(), file=table_name, sep=",", row.names=F, col.names=F)
write.table("", file=table_name, sep=",", row.names=F, col.names=F, append=T)
write.table(summary(model)$methTitle, file=table_name, sep=",", row.names=F, col.names=F, append=T)
write.table(as.character(summary(model)$call)[c(2, 4)], file=table_name, sep=",", row.names=F, col.names=F, append=T)
fixed_effects<-cbind(row.names(summary(model)$coeff), as.data.frame(summary(model)$coeff))
names(fixed_effects)[1]<-"Coefficient"
write.table("", file=table_name, sep=",", row.names=F, col.names=F, append=T)
write.table("Fixed Effects", file=table_name, sep=",", row.names=F, col.names=F, append=T)
write.table(fixed_effects, file=table_name, sep=",", row.names=F, col.names=names(fixed_effects), append=T)
write.table("", file=table_name, sep=",", row.names=F, col.names=F, append=T)
write.table("Random Effects VCV Matrix", file=table_name, sep=",", row.names=F, col.names=F, append=T)
random_effects<-as.data.frame(summary(model)$varcor[[1]][c(1:4), c(1:4)])
random_effects<-cbind(row.names(random_effects), random_effects)
names(random_effects)[1]<-" "
write.table(random_effects, file=table_name, sep=",", row.names=F, col.names=names(random_effects), append=T)
write.table("", file=table_name, sep=",", row.names=F, col.names=F, append=T)
write.table(paste0("Data came from ", length(unique(LT_data$LT_ID)), " cohorts (Country:Year)"), file=table_name, sep=",", row.names=F, col.names=F, append=T)

# Add the random and fixed coefficients to get specific LT coefficients
random<-ranef(model)$LT_ID
beta<-model@beta
for(i in 1:4){
	random[,i]<-random[,i] + beta[i]
}

# Get the sex specific brass coefficients
estimates_males<-random[,c(1,2)]
estimates_females<-estimates_males + random[,c(3,4)]
estimates_males$Sex<-as.factor("Males")
estimates_females$Sex<-as.factor("Females")
estimates_males$LT<-row.names(estimates_males)
estimates_females$LT<-row.names(estimates_females)
names(estimates_males)<-c("alpha", "beta", "Sex", "LT")
names(estimates_females)<-names(estimates_males)

# OK now we have all the available brass coefs going back as far 1961, lets see what we have interms of nutritional data

# Get the Nutrition data
FB_data<-read.csv("clean_data/Clean_FB.csv")

# In the FB data there are regional supplies (incl. for the world) - these are labelled as country NA, and we will drop them
FB_data<-FB_data[-which(is.na(FB_data$Country) == T),]

# Make the ID
FB_data$Food_ID<-paste0(FB_data$Country, "_", FB_data$Year)

# Lets add in the GDP predictors from Maddison project
SES_data<-read.csv("clean_data/Clean_GDP_Mad.csv")
head(SES_data)

# Be sure to throw out any data before 1961
SES_data<-SES_data[-which(SES_data$year < 1961),]

# Make the ID
SES_data$SES_ID<-paste0(SES_data$countrycode, "_", SES_data$year)

# Lets add in the GDP predictors from the WB
SES_data2<-read.csv("clean_data/Clean_GDP_WB.csv")
head(SES_data2)

# Be sure to throw out any data before 1961
SES_data2<-SES_data2[-which(SES_data2$Year < 1961),]

# Make the ID
SES_data2$SES_ID<-paste0(SES_data2$Country.Code, "_", SES_data2$Year)

# Aggregate all the IDs and lets see what we have data for
All_ID<-unique(c(estimates_males$LT, FB_data$Food_ID, SES_data$SES_ID, SES_data2$SES_ID))

# Create a data frame
all_males<-data.frame(ID=All_ID)
all_males$Country<-unlist(lapply(strsplit(as.character(all_males$ID), "_"), "[[", 1))
all_males$Year<-unlist(lapply(strsplit(as.character(all_males$ID), "_"), "[[", 2))

# Add in the FB data
all_males<-cbind(all_males, FB_data[match(all_males$ID, FB_data$Food_ID), c("Protein.kcal", "Carbo.kcal", "Fat.kcal", "Prop_animal_prot")])

# Add in the GDP data from the Maddison PJ
all_males<-cbind(all_males, SES_data[match(all_males$ID, SES_data$SES_ID), c("cgdppc")])
names(all_males)[dim(all_males)[2]]<-"GDP_perCapita"

# Add in the GDP data from the WB
all_males<-cbind(all_males, SES_data2[match(all_males$ID, SES_data2$SES_ID), c("GDP_WB")])
names(all_males)[dim(all_males)[2]]<-"GDP_perCapita_WB"

# Create a copy for the females
all_females<-all_males

# Add in the sex to each data frame
all_males$Sex<-"Males"
all_females$Sex<-"Females"

# Now match to the sex specific alpha and beta estimates
all_males<-cbind(all_males, estimates_males[match(all_males$ID, estimates_males$LT), c("alpha", "beta")])
all_females<-cbind(all_females, estimates_females[match(all_females$ID, estimates_females$LT), c("alpha", "beta")])

# Sort by country and year
all_males<-all_males[order(all_males$Country, all_males$Year),]
all_females<-all_females[order(all_females$Country, all_females$Year),]

# Bind the sexes together
full_data<-rbind(all_males, all_females)

# Now lets calculate the life-expectancy and SD in age at death for each lifetable - while we are doing that we also get predictions for lx in the LT data based on each country's alpha and beta

# Recode sex to be as in the full data
tag<-which(LT_data$Sex == 1)
LT_data$Sex2<-"Females"
LT_data$Sex2[tag]<-"Males"
LT_data$Sex2<-as.factor(LT_data$Sex2)
LT_data$Sex2<-relevel(LT_data$Sex2, ref="Males")

# Add in some NA columns to hold data
full_data$e0<-NA
full_data$l5<-NA
full_data$l60<-NA
full_data$s0<-NA
LT_data$pred_l.x.<-NA

for(j in 1:dim(full_data)[1]){
	
	# Get the jth lifetable
	tag<-which(LT_data$LT_ID == full_data$ID[j] & LT_data$Sex2 == full_data$Sex[j])
	if(length(tag != 0)){
		subset<-LT_data[tag,]
		class_0<-subset[1,]
		class_0$Age<-0
		class_0$l.x.<-1
		subset<-rbind(class_0, subset)
		
		# Get the full lifetable
		sexes<-c("male", "female")
		LT_j<-LifeTable(x=subset$Age, lx=subset$l.x., sex=sexes[subset$Sex[1]])$lt
		LT_j$sx<-0
		n_class<-dim(LT_j)[1]
		# Get the SD in life expectancy
		for(i in 1:n_class){
			LT_j$sx[i]<-sqrt(sum(LT_j$mx[c(i:n_class)] * ((LT_j$x[c(i:n_class)] + LT_j$ax[c(i:n_class)]) - LT_j$ex[1])^2))
		}
		
		# Add in the life expectancy and survival at ages 5 and 60
		full_data$e0[j]<-LT_j$ex[1]
		full_data$l5[j]<-subset$l.x.[which(subset$Age == 5)]
		full_data$l60[j]<-subset$l.x.[which(subset$Age == 60)]
		full_data$s0[j]<-LT_j$sx[1]
		
		# Predict the lx from the alpha and beta values
		LT_data$pred_l.x.[tag]<-invlogit(full_data$alpha[j] + full_data$beta[j] * LT_data$v1[tag] + LT_data$v2[tag] + LT_data$v3[tag])
	}
}

# Number of countries
length(unique(full_data$Country))
length(unique(full_data$ID))

# This analysis is performed on the complete cases, so I need to drop those with missing data in any of the following columns
check<-c("Protein.kcal", "Carbo.kcal", "Fat.kcal", "GDP_perCapita", "alpha", "beta")
missing<-apply(full_data[,check], 2, is.na)
complete_cases<-full_data[-which(apply(missing, 1, sum) > 0),]

# Save the full dataset minus any mising data
write.table(complete_cases, file="brass_data/Brass_complete_cases.csv", sep=",", row.names=F, col.names=names(complete_cases))

# R1 suggested splitting a random subset of the data and repeating the analysis - to do this I will randomly select the 50% of the LT IDs then match in each sex (so we have the same random assemblage of years/countries for each sex)

# Select a random subset of LT IDs
set.seed(123)
IDs<-unique(complete_cases$ID)
random_1<-IDs[sample(seq(1, length(IDs), 1), round(length(IDs)/2))]

# Pull out those selected
tag<-which(is.na(match(complete_cases$ID, random_1)) == F)
sub1<-complete_cases[tag,]

# And the other half
sub2<-complete_cases[-tag,]

# Save them
write.table(sub1, file="brass_data/Brass_random1.csv", sep=",", row.names=F, col.names=names(sub1))
write.table(sub2, file="brass_data/Brass_random2.csv", sep=",", row.names=F, col.names=names(sub2))

# Subset the data for those in Lieberman et al. 2020 AJCN, which are Australia, Canada, Czech Republic, Finland, France, Germany, Ireland, Japan, Korea, Netherlands, Norway, Poland, UK and USA
countries_lieb<-c("AUS", "CAN", "CZE", "FIN", "FRA", "DEU", "IRL", "JPN", "KOR", "NLD", "NOR", "POL", "GBR", "USA")
sub_data<-complete_cases[which(is.na(match(complete_cases$Country, countries_lieb)) == F),]
sub_data<-droplevels(sub_data)

# Save the sub dataset
write.table(sub_data, file="brass_data/Brass_subset.csv", sep=",", row.names=F, col.names=names(sub_data))

# Lets check correlations between l.x. and predicted l.x.

# Cull out the missing data from the LT data to establish correltations
LT_data<-LT_data[-which(is.na(match(LT_data$LT_ID, complete_cases$ID)) == T),]

r_sqm<-round(cor(LT_data$l.x.[which(LT_data$Sex2=="Males")], LT_data$pred_l.x.[which(LT_data$Sex2=="Males")])^2, 3)
r_sqf<-round(cor(LT_data$l.x.[which(LT_data$Sex2=="Females")], LT_data$pred_l.x.[which(LT_data$Sex2=="Females")])^2, 3)

# I had plotted this but removed it for space, lets just report the r2
r_sqm
r_sqf

# Lastly lets create a dataset with imputations for the missing data
# This last chunk of code below to do the imputation of missing data was done Shinichi

# data with missingness - note I am excluding the propotion of animal protein from this as we do not analyse that directly
dat <- full_data[,-7]

# reformat to wide format
wide_dat<-data.frame(ID = unique(dat$ID))
wide_dat<-cbind(wide_dat, dat[match(wide_dat$ID, dat$ID), c(2:8)])
males<-dat[which(dat$Sex == "Males"),]
females<-dat[which(dat$Sex == "Females"),]
wide_dat<-cbind(wide_dat, males[match(wide_dat$ID, males$ID), c(10:11)])
names(wide_dat)[c(9:10)]<-c("alpha_m", "beta_m")
wide_dat<-cbind(wide_dat, females[match(wide_dat$ID, females$ID), c(10:11)])
names(wide_dat)[c(11:12)]<-c("alpha_f", "beta_f")

names(wide_dat)

md.pattern(wide_dat)
summary(wide_dat)
dim(wide_dat)

# just imputing just alpha and beta
# make character into factor 
# needs to make Country into intergers
# standarise for brining onto the same scale (good for model) but need to turn it back to the orginal (remember original SD and mean)

# Log the GDP
wide_dat[,c(7:8)]<-apply(wide_dat[,c(7:8)], 2, log)

# save the means and sds, which will be needed to backtransform
mus<-apply(wide_dat[,-c(1:3)], 2, mean, na.rm=T)
sds<-apply(wide_dat[,-c(1:3)], 2, sd, na.rm=T)

# and standardise
sdat<-wide_dat
sdat[,-c(1:3)]<-apply(sdat[,-c(1:3)], 2, scale)
sdat$Country<-as.numeric(as.factor(sdat$Country))

dim(sdat)
md.pattern(sdat)
  
#ggpairs(sdat)

# set up redictor matrix and imputation methods
pred_matrix <- make.predictorMatrix(sdat)
imp_method <- make.method(sdat)

# slope and random slope (2)
# think this is fine for all
pred_matrix[ , "Year"] <- 2

# Setting as 3 will include cluster mean 
pred_matrix[,  c("Protein.kcal", "Carbo.kcal", "Fat.kcal", "GDP_perCapita", "GDP_perCapita_WB", "alpha_m", "beta_m", "alpha_f", "beta_f")] <- 1

# # our cluster (-2)
pred_matrix[, "Country"] <- -2

# stting 0 for non-missing data, and also GDP_WB - we dont need this
no_missing <- c("ID", "Country", "Year")
pred_matrix[no_missing, ] <- 0
pred_matrix[, "ID"] <- 0

# also put 0 for diag
diag(pred_matrix) <- 0

# checking
pred_matrix

# setting methods 
imp_method[c("Protein.kcal", "Carbo.kcal", "Fat.kcal", "GDP_perCapita", "GDP_perCapita_WB", "alpha_m", "beta_m", "alpha_f", "beta_f")] <- "2l.pmm" # obs level
imp_method

# checking
# needs fast computer!!!!
m<-50
imp <- mice(sdat, 
            m = m, 
            maxit = 20,
            method = imp_method,
            predictorMatrix = pred_matrix)

# look what happens
# looks mostly good
# densityplot(imp, ~ Protein.kcal) 
# densityplot(imp, ~ Carbo.kcal) 
# densityplot(imp, ~ Fat.kcal)
# densityplot(imp, ~ GDP_perCapita)
# densityplot(imp, ~ GDP_perCapita_WB)
# densityplot(imp, ~ alpha_m)
# densityplot(imp, ~ beta_m)
# densityplot(imp, ~ alpha_f)
# densityplot(imp, ~ beta_f)


# some are good but not others
# plot(imp) 

# get imputed data as a list
# mean of these imputed will be probably match with the best estimate
GDP <- mids2mitml.list(imp)

# Average over the lists
set_i<-GDP[[1]][,-c(1:3)] * (1/m)
for(i in 2:m){
	set_i<-set_i + GDP[[i]][,-c(1:3)] * (1/m)
}

# Unscale
for(i in 1:dim(set_i)[2]){
	set_i[,i]<-set_i[,i] * sds[i] + mus[i]
}

# Unlog the GDP
set_i[,c("GDP_perCapita", "GDP_perCapita_WB")]<-exp(set_i[,c("GDP_perCapita", "GDP_perCapita_WB")])

# Add back in
wide_dat[,-c(1:3)]<-set_i

# Turn back in to long format
males<-wide_dat[,-c(11,12)]
females<-wide_dat[,-c(9,10)]
names(males)[c(9,10)]<-c("alpha", "beta")
names(females)[c(9,10)]<-c("alpha", "beta")
males$Sex<-"Males"
females$Sex<-"Females"

# Check that the original complete obs match the imputed
full_imp<-rbind(males, females)
plot(full_imp$Protein.kcal, full_data$Protein.kcal)
plot(full_imp$GDP_perCapita_WB, full_data$GDP_perCapita_WB)
# Yep

full_imp$imputed<-is.na(match(full_imp$ID, full_data$ID))

# Work out which entries have been imputed - useful to know
data_missing<-as.data.frame(is.na(full_data[,c("Protein.kcal", "Carbo.kcal", "Fat.kcal", "GDP_perCapita", "GDP_perCapita_WB", "alpha", "beta")]))
names(data_missing)<-paste0(names(data_missing), "_imp")
full_imp<-cbind(full_imp, data_missing)

# Drop any imputed alpha and beta values - performance here is V bad and we are imputing too much 
full_imp2<-full_imp[-which(full_imp$alpha_imp == T | full_imp$beta_imp == T),]

write.table(full_imp2, file="brass_data/Brass_imputed.csv", sep=",", row.names=F, col.names=names(full_imp2))
