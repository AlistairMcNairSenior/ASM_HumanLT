
# Clean up 
rm(list=ls())

# Working directory
wd<-"/Users/alistairsenior/Dropbox (Sydney Uni)/Human lifetables and Nutrition" # Work iMac
#wd<-"/Users/asenior/Dropbox (Sydney Uni)/Human lifetables and Nutrition"
#wd<-"/Users/asenior/Dropbox (Sydney Uni)/Human lifetables and Nutrition" # Home iMac
setwd(wd)

# Load libraries
library(mgcv)
library(doSNOW)

source("scripts/0. Header_Functions.R")

#################################################
############# GAMs on alpha and beta ############
#################################################

# We will fit all of the following models in each sex seperately
formulas_list<-list()
formulas_list[[1]]<-" ~ 1 + s(Country, bs=\"re\")"
formulas_list[[2]]<-" ~ s(Protein.kcal, Carbo.kcal, Fat.kcal, k=k_nut) + s(Country, bs=\"re\")"
formulas_list[[3]]<-" ~ s(Year, k=10, bs=\"cr\") + s(Country, bs=\"re\")"
formulas_list[[4]]<-" ~ s(GDP_perCapita, k=10, bs=\"cr\") + s(Country, bs=\"re\")"
formulas_list[[5]]<-" ~ s(Protein.kcal, Carbo.kcal, Fat.kcal, k=k_nut) + s(Year, k=10, bs=\"cr\") + s(Country, bs=\"re\")"
formulas_list[[6]]<-" ~ s(Protein.kcal, Carbo.kcal, Fat.kcal, k=k_nut) + s(GDP_perCapita, k=10, bs=\"cr\") + s(Country, bs=\"re\")"
formulas_list[[7]]<-" ~ s(Year, k=10, bs=\"cr\") + s(GDP_perCapita, k=10, bs=\"cr\") + s(Country, bs=\"re\")"
formulas_list[[8]]<-" ~ te(Protein.kcal, Carbo.kcal, Fat.kcal, Year, bs=c(\"tp\", \"cr\"), d=c(3,1), k=c(k_nut, 7)) + s(Country, bs=\"re\")"
formulas_list[[9]]<-" ~ te(Protein.kcal, Carbo.kcal, Fat.kcal, GDP_perCapita, bs=c(\"tp\", \"cr\"), d=c(3,1), k=c(k_nut, 7)) + s(Country, bs=\"re\")"
formulas_list[[10]]<-" ~ te(Year, GDP_perCapita, k=10) + s(Country, bs=\"re\")"
formulas_list[[11]]<-" ~ te(Protein.kcal, Carbo.kcal, Fat.kcal, Year, bs=c(\"tp\", \"cr\"), d=c(3,1), k=c(k_nut, 7)) + s(GDP_perCapita, k=10, bs=\"cr\") + s(Country, bs=\"re\")"
formulas_list[[12]]<-" ~ te(Protein.kcal, Carbo.kcal, Fat.kcal, GDP_perCapita, bs=c(\"tp\", \"cr\"), d=c(3,1), k=c(k_nut, 7)) + s(Year, k=10, bs=\"cr\") + s(Country, bs=\"re\")"
formulas_list[[13]]<-" ~ te(Year, GDP_perCapita, k=10) + s(Protein.kcal, Carbo.kcal, Fat.kcal, k=k_nut) + s(Country, bs=\"re\")"

# The sexes over which we run the two models
sex<-c("Males", "Females")

#################################################################
################## Complete Cases Analysis ######################
#################################################################

## Start with the full complete cases dataset
full_data<-read.csv("brass_data/Brass_complete_cases.csv", stringsAsFactors=T)

# Models for macronutrient supply and alpha and beta
# We will set gamma as the log(n)/2, which is BIC-like (see ?bam) - however note that n is the rows of full_data/2 because sexes are analysed seperately
n<-dim(full_data)[1]/2
gamma<-log(n)/2

# Run in parallel using doSNOW - one model per core, then wash rinse and repeat for the other sexes via the loop
cl<-makeCluster(length(formulas_list), outfile="")
registerDoSNOW(cl)

# Open the loop for the sexes
for(i in 1:length(sex)){
	
	# Get the ith dataset for sex i
	data_i<-full_data[which(full_data$Sex == sex[i]),]
	
	# Run in parallel
	models_list<-foreach(p = 1:length(formulas_list)) %dopar% {
		require(mgcv)
		k_nut<-100
		form<-formulas_list[[p]]
		GAM<-gam(list(as.formula(paste0("alpha", form)), as.formula(paste0("beta", form))), data=data_i, family=mvn(d=2), gamma=gamma)
		return(GAM)	
	}
	
	# Save the models
	save(models_list, file=paste0("models/Complete_cases_GAMS_", sex[i], ".rdata"))
		
}

# Load the males and estimate the deviance explained, AIC, delta AIC and weights
load("models/Complete_cases_GAMS_Males.rdata")
#lapply(models_list, gam.check)
summaries<-lapply(models_list, summary)
dev_males<-unlist(lapply(summaries, "[[", 14)) * 100
AIC_males<-unlist(lapply(models_list, AIC))
delta_males<-AIC_males - min(AIC_males)
weights_males<-exp(-0.5 * delta_males) / sum(exp(-0.5 * delta_males))
models_males<-models_list

# Repeat for females
load("models/Complete_cases_GAMS_Females.rdata")
#lapply(models_list, gam.check)
summaries<-lapply(models_list, summary)
dev_females<-unlist(lapply(summaries, "[[", 14)) * 100
AIC_females<-unlist(lapply(models_list, AIC))
delta_females<-AIC_females - min(AIC_females)
weights_females<-exp(-0.5 * delta_females) / sum(exp(-0.5 * delta_females))
models_females<-models_list

# Put the AIC tables together and save
AIC_table<-cbind(dev_males, AIC_males, delta_males, weights_males, dev_females, AIC_females, delta_females, weights_females)
write.table(AIC_table, file="tables/Table_2.csv", sep=",", row.names=F, col.names=colnames(AIC_table))

# In both cases model 11 is favoured
AIC_favoured_models<-list(models_males[[11]], models_females[[11]])
save(AIC_favoured_models, file="models/Complete_cases_AIC_GAMS.rdata")

# Create Table S6
write.GAM(AIC_favoured_models[[1]], csv.file="tables/Table_S6_males(complete_cases).csv")
write.GAM(AIC_favoured_models[[2]], csv.file="tables/Table_S6_females(complete_cases).csv")

#################################################################
####################### Lieberman subet #########################
#################################################################

# Subset the data for those in Lieberman et al. 2020 AJCN, which are Australia, Canada, Czech Republic, Finland, France, Germany, Ireland, Japan, Korea, Netherlands, Norway, Poland, UK and USA
full_data<-read.csv("brass_data/Brass_subset.csv", stringsAsFactors=T)

# Refit model 11
form<-formulas_list[[11]]

# We will set gamma as the log(n)/2, which is BIC-like (see ?bam)
n<-dim(full_data)[1]/2
gamma<-log(n)/2

# Run with doSNOW
cl<-makeCluster(length(sex), outfile="")
registerDoSNOW(cl)
AIC_favoured_models<-foreach(p = 1:length(sex)) %dopar% {
	
	# Load the details locally	
	require(mgcv)
	k_nut<-30
	
	# Get the ith dataset for sex i
	data_i<-full_data[which(full_data$Sex == sex[p]),]
	
	# Run the GAM
	GAM<-gam(list(as.formula(paste0("alpha", form)), as.formula(paste0("beta", form))), data=data_i, family=mvn(d=2), gamma=gamma)
	
	# Return the GAM
	return(GAM)		
}

# Save the models
save(AIC_favoured_models, file="models/Subset_AIC_GAMS.rdata")
load("models/Subset_AIC_GAMS.rdata")

# Create Table S7
write.GAM(AIC_favoured_models[[1]], csv.file="tables/Table_S7_males(subset).csv")
write.GAM(AIC_favoured_models[[2]], csv.file="tables/Table_S7_females(subset).csv")

# NOTE table S8 is the LMM coefficients in the script 4. Estimate_Brass

# tables S9 and S10 are the associations between GDP/nutrient supply and time in the script for Figure S4

#################################################################
####################### Imputed dataset #########################
#################################################################

# Repeat with the imputed dataset
full_data<-read.csv("brass_data/Brass_imputed.csv", stringsAsFactors=T)

# Refit model 11
form<-formulas_list[[11]]

# We will set gamma as the log(n)/2, which is BIC-like (see ?bam) - however note that n is the rows of full_data/2 because sexes are analysed seperately
n<-dim(full_data)[1]/2
gamma<-log(n)/2

# Run with doSNOW
cl<-makeCluster(length(sex), outfile="")
registerDoSNOW(cl)
AIC_favoured_models<-foreach(p = 1:length(sex)) %dopar% {
	
	# Load the details locally	
	require(mgcv)
	k_nut<-100
	
	# Get the ith dataset for sex i
	data_i<-full_data[which(full_data$Sex == sex[p]),]
	
	# Run the GAM
	GAM<-gam(list(as.formula(paste0("alpha", form)), as.formula(paste0("beta", form))), data=data_i, family=mvn(d=2), gamma=gamma)
	
	# Return the GAM
	return(GAM)		
}

# Save the models
save(AIC_favoured_models, file="models/Imputed_AIC_GAMS.rdata")
load("models/Imputed_AIC_GAMS.rdata")

write.GAM(AIC_favoured_models[[1]], csv.file="tables/Table_S11_males(imputed).csv")
write.GAM(AIC_favoured_models[[2]], csv.file="tables/Table_S11_females(imputed).csv")

#################################################################
####################### Random dataset1 #########################
#################################################################

# Repeat with the random dataset 1
full_data<-read.csv("brass_data/Brass_random1.csv", stringsAsFactors=T)

# Refit model 11
form<-formulas_list[[11]]

# We will set gamma as the log(n)/2, which is BIC-like (see ?bam) - however note that n is the rows of full_data/2 because sexes are analysed seperately
n<-dim(full_data)[1]/2
gamma<-log(n)/2

# Run with doSNOW
cl<-makeCluster(length(sex), outfile="")
registerDoSNOW(cl)
AIC_favoured_models<-foreach(p = 1:length(sex)) %dopar% {
	
	# Load the details locally	
	require(mgcv)
	k_nut<-50
	
	# Get the ith dataset for sex i
	data_i<-full_data[which(full_data$Sex == sex[p]),]
	
	# Run the GAM
	GAM<-gam(list(as.formula(paste0("alpha", form)), as.formula(paste0("beta", form))), data=data_i, family=mvn(d=2), gamma=gamma)
	
	# Return the GAM
	return(GAM)		
}


# Save the models
save(AIC_favoured_models, file="models/Random1_AIC_GAMS.rdata")
load("models/Random1_AIC_GAMS.rdata")

write.GAM(AIC_favoured_models[[1]], csv.file="tables/Table_S12_males(random1_subset).csv")
write.GAM(AIC_favoured_models[[2]], csv.file="tables/Table_S12_females(random1_subset).csv")

#################################################################
####################### Random dataset2 #########################
#################################################################

# Repeat with the random dataset 2
full_data<-read.csv("brass_data/Brass_random2.csv", stringsAsFactors=T)

# Refit model 11
form<-formulas_list[[11]]

# We will set gamma as the log(n)/2, which is BIC-like (see ?bam) - however note that n is the rows of full_data/2 because sexes are analysed seperately
n<-dim(full_data)[1]/2
gamma<-log(n)/2

# Run with doSNOW
cl<-makeCluster(length(sex), outfile="")
registerDoSNOW(cl)
AIC_favoured_models<-foreach(p = 1:length(sex)) %dopar% {
	
	# Load the details locally	
	require(mgcv)
	k_nut<-50
	
	# Get the ith dataset for sex i
	data_i<-full_data[which(full_data$Sex == sex[p]),]
	
	# Run the GAM
	GAM<-gam(list(as.formula(paste0("alpha", form)), as.formula(paste0("beta", form))), data=data_i, family=mvn(d=2), gamma=gamma)
	
	# Return the GAM
	return(GAM)		
}


# Save the models
save(AIC_favoured_models, file="models/Random2_AIC_GAMS.rdata")
load("models/Random2_AIC_GAMS.rdata")

write.GAM(AIC_favoured_models[[1]], csv.file="tables/Table_S13_males(random2_subset).csv")
write.GAM(AIC_favoured_models[[2]], csv.file="tables/Table_S13_females(random2_subset).csv")
