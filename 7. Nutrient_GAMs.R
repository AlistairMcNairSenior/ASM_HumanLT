
# Clean up 
rm(list=ls())

# Working directory
wd<-"/Users/alistairsenior/Dropbox/Human lifetables and Nutrition"
wd<-"/Users/alistairsenior/Dropbox (Sydney Uni)/Human lifetables and Nutrition"
#wd<-"/Users/asenior/Dropbox (Sydney Uni)/Human lifetables and Nutrition"
setwd(wd)

# Load libraries
library(mgcv)

source("0.Header_Functions.R")

#################################################
############# GAMs on alpha and beta ############
#################################################

# Full data
full_data<-read.csv("Brass_complete_cases.csv")

# Models for macronutrient supply and alpha and beta
# We will set gamma as the log(n)/2, which is BIC-like (see ?bam)
gamma<-log(dim(full_data)[1])/2

# Additive model with time
form<- " ~ s(Protein.kcal, Carbo.kcal, Fat.kcal, by=Sex, k=125) + s(Year, k=10, bs=\"cr\") + Sex + s(Country, bs=\"re\")"
GAM<-gam(list(as.formula(paste0("alpha", form)), as.formula(paste0("beta", form))), data=full_data, family=mvn(d=2), gamma=gamma)
save(GAM, file="GAM_add.rdata")
load("GAM_add.rdata")
summary(GAM)
gam.check(GAM)
write.GAM(GAM, csv.file="Table_S6_(supply+time).csv")
aic_3d<-AIC(GAM)

# Interaction with time
form<- " ~ te(Protein.kcal, Carbo.kcal, Fat.kcal, Year, by=Sex, bs=c(\"tp\", \"cr\"), d=c(3,1), k=c(125, 7)) + Sex + s(Country, bs=\"re\")"
GAM<-gam(list(as.formula(paste0("alpha", form)), as.formula(paste0("beta", form))), data=full_data, family=mvn(d=2), gamma=gamma)
save(GAM, file="GAM_int.rdata")
load("GAM_int.rdata")
summary(GAM)
gam.check(GAM)
write.GAM(GAM, csv.file="Table_S7_(supply:time).csv")
aic_4d<-AIC(GAM)

# Difference in AIC between models
aic_4d - aic_3d

# Interaction with time + GDP
form<- " ~ te(Protein.kcal, Carbo.kcal, Fat.kcal, Year, by=Sex, bs=c(\"tp\", \"cr\"), d=c(3,1), k=c(125, 7)) + Sex + s(Country, bs=\"re\") + s(GDP_perCapita, k=10)"
GAM<-gam(list(as.formula(paste0("alpha", form)), as.formula(paste0("beta", form))), data=full_data, family=mvn(d=2), gamma=gamma)
save(GAM, file="GAM_int_GDP.rdata")
load("GAM_int_GDP.rdata")
summary(GAM)
gam.check(GAM)
write.GAM(GAM, csv.file="Table_S8_(supply:time+gdp).csv")
aic_4d_GDP<-AIC(GAM)

# Difference in AIC
aic_4d_GDP - aic_4d

# We should also check the model without macros to see if we are getting any extrac information from there
form<- " ~ te(Year, GDP_perCapita, by=Sex, k=10) + Sex + s(Country, bs=\"re\")"
GAM<-gam(list(as.formula(paste0("alpha", form)), as.formula(paste0("beta", form))), data=full_data, family=mvn(d=2), gamma=gamma)
save(GAM, file="GAM_Year_GDP.rdata")
load("GAM_Year_GDP.rdata")
summary(GAM)
gam.check(GAM)
write.GAM(GAM, csv.file="Table_S9_(gdp:time).csv")
aic_GDP<-AIC(GAM)

# Difference in AIC
aic_4d_GDP - aic_GDP

# Subset the data for those in Lieberman et al. 2020 AJCN, which are Australia, Canada, Czech Republic, Finland, France, Germany, Ireland, Japan, Korea, Netherlands, Norway, Poland, UK and USA
sub_data<-read.csv("Brass_subset.csv")

# Refit the models
# We will set gamma as the log(n)/2, which is BIC-like (see ?bam)
gamma<-log(dim(sub_data)[1])/2

# Additive model with time
form<- " ~ s(Protein.kcal, Carbo.kcal, Fat.kcal, by=Sex, k=50) + s(Year, k=7, bs=\"cr\") + Sex + s(Country, bs=\"re\")"
GAM<-gam(list(as.formula(paste0("alpha", form)), as.formula(paste0("beta", form))), data=sub_data, family=mvn(d=2), gamma=gamma)
save(GAM, file="GAM_add_sub.rdata")
load("GAM_add_sub.rdata")
summary(GAM)
gam.check(GAM)
# Note I have increased k for nutrients all the way up to 200, with out edf substantially changing or the model predictions substantially changing
write.GAM(GAM, csv.file="Table_S10_(supply+time_SUBSET).csv")
aic_3d<-AIC(GAM)

# Interaction with time
form<- " ~ te(Protein.kcal, Carbo.kcal, Fat.kcal, Year, by=Sex, bs=c(\"tp\", \"cr\"), d=c(3,1), k=c(50, 5)) + Sex + s(Country, bs=\"re\")"
GAM<-gam(list(as.formula(paste0("alpha", form)), as.formula(paste0("beta", form))), data=sub_data, family=mvn(d=2), gamma=gamma)
save(GAM, file="GAM_int_sub.rdata")
load("GAM_int_sub.rdata")
summary(GAM)
gam.check(GAM)
write.GAM(GAM, csv.file="Table_S11_(supply:time_SUBSET).csv")
aic_4d<-AIC(GAM)

# Difference in AIC between models
aic_4d - aic_3d

# Interaction with time + gdp
form<- " ~ te(Protein.kcal, Carbo.kcal, Fat.kcal, Year, by=Sex, bs=c(\"tp\", \"cr\"), d=c(3,1), k=c(50, 5)) + Sex + s(Country, bs=\"re\") + s(GDP_perCapita, k=10)"
GAM<-gam(list(as.formula(paste0("alpha", form)), as.formula(paste0("beta", form))), data=sub_data, family=mvn(d=2), gamma=gamma)
save(GAM, file="GAM_int_GDP_sub.rdata")
load("GAM_int_GDP_sub.rdata")
summary(GAM)
gam.check(GAM)
write.GAM(GAM, csv.file="Table_S12_(supply:time+gdp_SUBSET).csv")
aic_4d_GDP<-AIC(GAM)

aic_4d_GDP - aic_4d

# GDP interaction with time 
form<- " ~ te(Year, GDP_perCapita, by=Sex, k=10) + Sex + s(Country, bs=\"re\")"
GAM<-gam(list(as.formula(paste0("alpha", form)), as.formula(paste0("beta", form))), data=sub_data, family=mvn(d=2), gamma=gamma)
save(GAM, file="GAM_Year_GDP_sub.rdata")
load("GAM_Year_GDP_sub.rdata")
summary(GAM)
gam.check(GAM)
write.GAM(GAM, csv.file="Table_S13_(gdp:time_SUBSET).csv")
aic_GDP<-AIC(GAM)

# Difference in AIC
aic_4d_GDP - aic_GDP

# NOTE table S14 is the LMM coefficients in the script 4. Estimate_Brass

# tables S15 and S16 are the associations between GDP/nutrient supply and time in the script 12. Figure S1 
# Repeat with the imputed dataset
full_data<-read.csv("Brass_imputed.csv")

# Models for macronutrient supply and alpha and beta
# We will set gamma as the log(n)/2, which is BIC-like (see ?bam)
gamma<-log(dim(full_data)[1])/2

# #Additive model with time
form<- " ~ s(Protein.kcal, Carbo.kcal, Fat.kcal, by=Sex, k=125) + s(Year, k=10, bs=\"cr\") + Sex + s(Country, bs=\"re\")"
GAM<-gam(list(as.formula(paste0("alpha", form)), as.formula(paste0("beta", form))), data=full_data, family=mvn(d=2), gamma=gamma)
save(GAM, file="GAM_add_imp.rdata")
load("GAM_add_imp.rdata")
summary(GAM)
gam.check(GAM)
write.GAM(GAM, csv.file="Table_S17_(supply+time_IMP).csv")
aic_3d<-AIC(GAM)

# Interaction with time
form<- " ~ te(Protein.kcal, Carbo.kcal, Fat.kcal, Year, by=Sex, bs=c(\"tp\", \"cr\"), d=c(3,1), k=c(125, 7)) + Sex + s(Country, bs=\"re\")"
GAM<-gam(list(as.formula(paste0("alpha", form)), as.formula(paste0("beta", form))), data=full_data, family=mvn(d=2), gamma=gamma)
save(GAM, file="GAM_int_imp.rdata")
load("GAM_int_imp.rdata")
summary(GAM)
gam.check(GAM)
write.GAM(GAM, csv.file="Table_S18_(supply:time_IMP).csv")
aic_4d<-AIC(GAM)

# Difference in AIC between models
aic_4d - aic_3d

# Interaction with time + GDP
form<- " ~ te(Protein.kcal, Carbo.kcal, Fat.kcal, Year, by=Sex, bs=c(\"tp\", \"cr\"), d=c(3,1), k=c(125, 7)) + Sex + s(Country, bs=\"re\") + s(GDP_perCapita, k=10)"
GAM<-gam(list(as.formula(paste0("alpha", form)), as.formula(paste0("beta", form))), data=full_data, family=mvn(d=2), gamma=gamma)
save(GAM, file="GAM_int_GDP_imp.rdata")
load("GAM_int_GDP_imp.rdata")
summary(GAM)
gam.check(GAM)
write.GAM(GAM, csv.file="Table_S19_(supply:time+gdp_IMP).csv")
aic_4d_GDP<-AIC(GAM)

# Difference in AIC
aic_4d_GDP - aic_4d

# We should also check the model without macros to see if we are getting any extrac information from there
form<- " ~ te(Year, GDP_perCapita, by=Sex, k=10) + Sex + s(Country, bs=\"re\")"
GAM<-gam(list(as.formula(paste0("alpha", form)), as.formula(paste0("beta", form))), data=full_data, family=mvn(d=2), gamma=gamma)
save(GAM, file="GAM_Year_GDP_imp.rdata")
load("GAM_Year_GDP_imp.rdata")
summary(GAM)
gam.check(GAM)
write.GAM(GAM, csv.file="Table_S20_(supply:gdp_IMP).csv")
aic_GDP<-AIC(GAM)

# Difference in AIC
aic_4d_GDP - aic_GDP

