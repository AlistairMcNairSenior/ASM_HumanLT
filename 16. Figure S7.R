
# Clean up
rm(list=ls())

# Working directory
wd<-"/Users/alistairsenior/Dropbox (Sydney Uni)/Human lifetables and Nutrition" # Work iMac
#wd<-"/Users/asenior/Dropbox (Sydney Uni)/Human lifetables and Nutrition" # Work Macbook
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
library(ggcorrplot)

source("scripts/0. Header_Functions.R")

#################################################
#################### FIGURE S7 ##################
#################################################

full_data<-read.csv("brass_data/Brass_imputed.csv", stringsAsFactors=T)
head(full_data)
range(full_data$Year)
# Lets start by making two maps, that summarise data representation for the periods 1961-1990 and 1991-2017

# World data
world<-ne_countries(scale="medium", returnclass="sf")
world$CC<-countrycode(world$name_sort, origin="country.name", destination="iso3c")

# For the sake of visualisation put the YUG data in to SRB, SUN in to RUS, and CSK in to CZE
full_data_map<-full_data
tag<-which(full_data_map$Country == "SUN")
full_data_map$Country[tag]<-"RUS"

tag<-which(full_data_map$Country == "YUG")
full_data_map$Country[tag]<-"SRB"

tag<-which(full_data_map$Country == "CSK")
full_data_map$Country[tag]<-"CZE"

# Get the number of lifetables by country
summary_61<-ddply(full_data_map[which(full_data_map$Year < 1991),], .(Country), summarise, n_LT=length(unique(Year)))
summary_91<-ddply(full_data_map[which(full_data_map$Year >= 1991),], .(Country), summarise, n_LT=length(unique(Year)))

# Match the data in to the world data
world$LT_61<-summary_61$n_LT[match(world$CC, summary_61$Country)]
world$LT_91<-summary_91$n_LT[match(world$CC, summary_91$Country)]

# Countries with missing data have 0 LTs
world$LT_61<-world$LT_61/30 * 100
world$LT_91<-world$LT_91/27 * 100

# Map for 1961-1990
map_61<-ggplot(data=world) +
			geom_sf(aes(fill=LT_61), lwd=0) + 
			theme_bw() +
			theme(panel.grid.major=element_blank(), panel.grid.minor=element_blank(), legend.title=element_blank(), title=element_text(size=15)) + 
			labs(subtitle="1961 - 1990") + 
			scale_fill_gradient(low="linen", high="red", limits=c(1, 100)) + 
			coord_sf(ylim=c(-57, 80)) + 
			theme(legend.position=c(0.075, 0.25), legend.key.size=unit(0.4, "cm"))
			
# Map for 1991-2017
map_91<-ggplot(data=world) +
			geom_sf(aes(fill=LT_91), lwd=0) + 
			theme_bw() +
			theme(panel.grid.major=element_blank(), panel.grid.minor=element_blank(), title=element_text(size=15)) + 
			labs(subtitle="1991 - 2017") + 
			scale_fill_gradient(low="linen", high="red", limits=c(1, 100)) + 
			coord_sf(ylim=c(-57, 80)) +
			theme(legend.position="none")


# Fit GAMs fpr each outcome over time and add to the plot
gamma<-log(dim(full_data)[1])/2

# National changes in nutrition supplies over time - just need data from one sex
males<-full_data[which(full_data$Sex=="Males"),]
long<-data.frame(Supply=c(males$Protein.kcal, males$Carbo.kcal, males$Fat.kcal), Year=rep(males$Year, 3), Country=rep(males$Country, 3), Nutrient=as.factor(c(rep("Protein", dim(males)[1]), rep("Carbohydrate", dim(males)[1]), rep("Fat", dim(males)[1]))))

plot_nut<-ggplot(data=long, aes(x=Year, y=Supply, color=Nutrient)) + 
			geom_point(size=0.2, alpha=0.3) +
			theme_bw() +
			xlim(1960, 2020) +
			theme(axis.text=element_text(size=15), axis.title=element_text(size=15), legend.background=element_blank()) +
			theme(legend.position=c(0.5, 0.9), legend.key.size=unit(0.4, "cm"), legend.text=element_text(size=12), legend.title=element_text(size=15), legend.direction="horizontal", title=element_text(size=15)) +
			labs(color="", y="Supply kcal/capita/day") +
			guides(color=guide_legend(override.aes=list(size=4))) +
			ylim(0, 2625)

# GAM for change in nutrient supply over time
GAM_supply<-gam(Supply ~ s(Year, by=Nutrient, k=50) + Nutrient + s(Country, bs="re"), data=long, gamma=gamma)
summary(GAM_supply)
gam.check(GAM_supply)
# Some heteorgeneity in residuals (likely unaccounted for covariates, maybe should reformulate as MVN model as heteroscedastic residuals) - k fine
write.GAM(GAM_supply, csv.file="tables/Table_S9_(supply_time_IMP).csv")

# Predictions for the nutrients
Predictors<-data.frame(Year=rep(seq(1961, 2017, 1), 3), Nutrient=c(rep(as.character("Protein"), 57), rep(as.character("Carbohydrate"), 57), rep(as.character("Fat"), 57)))

prediction<-predict(GAM_supply, newdata=Predictors, exclude="s(Country)", newdata.guaranteed=T, se.fit=T, type="response")
Predictors$y<-prediction$fit
Predictors$se<-prediction$se.fit
Predictors$l.se<-Predictors$y - Predictors$se * 1.96
Predictors$u.se<-Predictors$y + Predictors$se * 1.96

plot_nut<-plot_nut +
			geom_ribbon(data=Predictors, aes(x=Year, y=y, ymin=l.se, ymax=u.se, fill=Nutrient), alpha=0.5, linetype=0) +
			geom_line(data=Predictors, aes(x=Year, y=y, color=Nutrient), size=0.1) + 
			guides(fill=F)

# Get the long format data from 1970 and 2010
long_70<-long[which(long$Year == 1970),]
long_10<-long[which(long$Year == 2010),]

# Make some densities
plot_dens70<-ggplot(long_70, aes(Supply, color=Nutrient)) +
			geom_density(aes(fill=Nutrient), alpha=0.5, linetype=0) + 
			theme_bw() + 
			theme(legend.position="none") +
			labs(x="Supply kcal/capita/day", y="Density", subtitle="1970") +
			theme(axis.text=element_text(size=15), axis.title=element_text(size=15), title=element_text(size=15)) +
			xlim(0, 2500)

plot_dens10<-ggplot(long_10, aes(Supply, color=Nutrient)) +
			geom_density(aes(fill=Nutrient), alpha=0.5, linetype=0) + 
			theme_bw() + 
			theme(legend.position="none") +
			labs(x="Supply kcal/capita/day", y="Density", subtitle="2010") +
			theme(axis.text=element_text(size=15), axis.title=element_text(size=15), title=element_text(size=15)) +
			xlim(0, 2500)			

# GDP over time
plot_GDP<-ggplot(males, aes(x=Year, y=GDP_perCapita)) +
			geom_point(size=0.2, alpha=0.2) +
			theme_bw() +
			xlim(1960, 2020) +
			theme(axis.text=element_text(size=15), axis.title=element_text(size=15)) +
			theme(title=element_text(size=15)) + 
			labs(y="GDP per Capita (2011 USD)")

# GAM for change in nutrient supply over time
GAM_GDP<-gam(GDP_perCapita ~ s(Year, k=50) + s(Country, bs="re"), data=males, gamma=gamma, family=gaussian(link="log"))
summary(GAM_GDP)
gam.check(GAM_GDP)
# Some heteorgeneity in residuals (likely unaccounted for covariates) - k fine
write.GAM(GAM_GDP, csv.file="tables/Table_S10_(GDP_time_IMP).csv")

# Predictions for the nutrients
Predictors<-data.frame(Year=rep(seq(1961, 2017, 1), 1))

prediction<-predict(GAM_GDP, newdata=Predictors, exclude="s(Country)", newdata.guaranteed=T, se.fit=T, type="link")
ilink<-family(GAM_GDP)$linkinv
Predictors$y<-prediction$fit
Predictors$se<-prediction$se.fit
Predictors<-transform(Predictors, l.se=ilink(y - 1.96 * se), u.se=ilink(y + 1.96 * se), y=ilink(y))

plot_GDP<-plot_GDP +
			geom_ribbon(data=Predictors, aes(x=Year, y=y, ymin=l.se, ymax=u.se), alpha=0.5, linetype=0) +
			geom_line(data=Predictors, aes(x=Year, y=y), size=0.1) + 
			guides(fill=F)
			

# Lets make correlograms for the following variables in 1970 and 2010: Lets also add in total energy as the sum of PCF
full_data$Total<-full_data$Protein.kcal + full_data$Carbo.kcal + full_data$Fat.kcal
cor_var<-c("Total", "Protein.kcal", "Carbo.kcal", "Fat.kcal", "GDP_perCapita")
cor_1970<-cor(full_data[which(full_data$Year == 1970), cor_var])
rownames(cor_1970)<-c("Total", "Protein", "Carb.", "Fat", "GDP")
colnames(cor_1970)<-c("Total", "Protein", "Carb.", "Fat", "GDP")
cor_2010<-cor(full_data[which(full_data$Year == 2010), cor_var])
rownames(cor_2010)<-c("Total", "Protein", "Carb.", "Fat", "GDP")
colnames(cor_2010)<-c("Total", "Protein", "Carb.", "Fat", "GDP")

plot_cor1<-ggcorrplot(cor_1970, lab=T, show.legend=F, type="lower") + labs(subtitle="1970") +
			theme(axis.text=element_text(size=15), axis.title=element_text(size=15), title=element_text(size=15))
plot_cor2<-ggcorrplot(cor_2010, lab=T, show.legend=F, type="lower") + labs(subtitle="2010") +
			theme(axis.text=element_text(size=15), axis.title=element_text(size=15), title=element_text(size=15))

# Now lets arrange all those plots
pdf("figures/Figure_S7.pdf", height=10, width=12)

grid.arrange(map_61+labs(title="A"), map_91+labs(title="B"), plot_nut+labs(title="C"), plot_dens70+labs(title="D"), plot_dens10+labs(title="E"), plot_GDP+labs(title="F"), plot_cor1+labs(title="G"), plot_cor2+labs(title="H"),
									layout_matrix=rbind(c(1,1,1,2,2,2),
			  					  					    c(3,3,4,4,5,5),
								  					    c(6,6,7,7,8,8)))

dev.off()

# Check the sample sizes for whole dataset
dim(males)
length(unique(full_data$Country))

# Check for 1970 and 2010
tag<-which(males$Year == 1970)
dim(males[tag,])

tag<-which(males$Year == 2010)
dim(males[tag,])
