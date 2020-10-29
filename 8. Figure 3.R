
# Clean up
rm(list=ls())

# Working directory
wd<-"/Users/alistairsenior/Dropbox (Sydney Uni)/Human lifetables and Nutrition" # Work iMac
#wd<-"/Users/asenior/Dropbox (Sydney Uni)/Human lifetables and Nutrition" # Macbook
setwd(wd)

# Load libraries
library(arm)
library(plyr)
library(ggplot2)
library(mgcv)
library(gridExtra)
library(MortalityLaws)
library(Cairo)

source("scripts/0. Header_Functions.R")

# Read in the standard from Wilmoth et al 2012
standard<-read.csv("clean_data/wilmoth_standard.csv")

#################################################
#################### FIGURE 3 ###################
#################################################

# Load the full data
full_data<-read.csv("brass_data/Brass_complete_cases.csv")

# So the model with macros * time + GDP has the best fit
load("models/Complete_cases_AIC_GAMS.rdata")
model_aic<-AIC_favoured_models[[2]]

# Get year of choice surfaces for PCF in females
year_plot<-2016
dataset_plot<-full_data[which(full_data$Year == year_plot), ]
med_GDP<-round(median(dataset_plot$GDP_perCapita, na.rm=T))
predict_val<-data.frame(Year=year_plot, GDP_perCapita=med_GDP, Sex=as.factor("Females"))

# Specify the desired layout for the surfaces
XYZ_list<-list()
XYZ_list[[1]]<-c("Protein.kcal", "Carbo.kcal", "Fat.kcal")
XYZ_list[[2]]<-c("Protein.kcal", "Fat.kcal", "Carbo.kcal")
XYZ_list[[3]]<-c("Carbo.kcal", "Fat.kcal", "Protein.kcal")

# Limits for the y.axis
y_limits<-list()
y_limits[[1]]<-c(1000, 2100)
y_limits[[2]]<-c(400, 1600)
y_limits[[3]]<-c(750, 1600)

labels_list<-c("Protein kcal/capita/day", "Carbohydrate kcal/capita/day", "Fat kcal/capita/day")

# Find the min and max set of predictions
mins<-array(NA, c(3,2))
maxs<-mins
for(j in 1:3){
	
	# Set the parameters for XYZ set j
	XYZ<-XYZ_list[[j]]
	z.val<-round(quantile(dataset_plot[,XYZ[3]])[3])
		
	gg_surfaces<-ggSurface(GAM=model_aic, data=dataset_plot, XYZ=XYZ, labels=XYZ, exclude=c("s(Country)", "s.1(Country)"), predict_val=predict_val, traits=c("alpha", "beta"), z.val=z.val, y.limits=y_limits[[j]])
	
	# Save the min and max values for scaling
	mins[j,1]<-min(gg_surfaces[[1]]$data$fit)
	maxs[j,1]<-max(gg_surfaces[[1]]$data$fit)	
	mins[j,2]<-min(gg_surfaces[[2]]$data$fit)
	maxs[j,2]<-max(gg_surfaces[[2]]$data$fit)	

}


# Find the min and max
min_use<-apply(mins, 2, min)
max_use<-apply(maxs, 2, max)

# Now refit to normalised scaling across surfaces and save those for presentation
surfaces_list_F<-list()
for(j in 1:3){
	
	# Set the parameters for XYZ set j
	XYZ<-XYZ_list[[j]]
	labels<-labels_list[match(XYZ, XYZ_list[[1]])]
	z.val<-round(quantile(dataset_plot[,XYZ[3]])[3])
	
	# Remake the surfaces sacles by the corss-surface min and max
	surfaces<-ggSurface(GAM=model_aic, data=dataset_plot, XYZ=XYZ, labels=labels, exclude=c("s(Country)", "s.1(Country)"), predict_val=predict_val, traits=c("alpha", "beta"), surf_min=min_use, surf_max=max_use, subtitle=paste0(labels[3], " = ", z.val), z.val=z.val, y.limits=y_limits[[j]])
	
	# Annotate
	surfaces[[1]]<-surfaces[[1]] + annotate("text", x = floor(min(dataset_plot[,XYZ[1]])), y = max(y_limits[[j]]), label = expression(italic(alpha)~~Females), hjust = 0, vjust = 1, size = 7)
	surfaces[[2]]<-surfaces[[2]] + annotate("text", x = floor(min(dataset_plot[,XYZ[1]])), y = max(y_limits[[j]]), label = expression(italic(beta)~~Females), hjust = 0, vjust = 1, size = 7)
	
	# Save them				
	surfaces_list_F[[j]]<-surfaces
}

################# REPEAT for males

# Get year of choice surfaces for PCF
model_aic<-AIC_favoured_models[[1]]
predict_val<-data.frame(Year=year_plot, GDP_perCapita=med_GDP, Sex=as.factor("Males"))

# Find the min and max set of predictions
mins<-array(NA, c(3,2))
maxs<-mins
for(j in 1:3){
	
	# Set the parameters for XYZ set j
	XYZ<-XYZ_list[[j]]
	z.val<-round(quantile(dataset_plot[,XYZ[3]])[3])
		
	gg_surfaces<-ggSurface(GAM=model_aic, data=dataset_plot, XYZ=XYZ, labels=XYZ, exclude=c("s(Country)", "s.1(Country)"), predict_val=predict_val, traits=c("alpha", "beta"), z.val=z.val, y.limits=y_limits[[j]])
	
	# Save the min and max values for scaling
	mins[j,1]<-min(gg_surfaces[[1]]$data$fit)
	maxs[j,1]<-max(gg_surfaces[[1]]$data$fit)	
	mins[j,2]<-min(gg_surfaces[[2]]$data$fit)
	maxs[j,2]<-max(gg_surfaces[[2]]$data$fit)	

}

# Find the min and max
min_use<-apply(mins, 2, min)
max_use<-apply(maxs, 2, max)

# Now refit to normalised scaling across surfaces and save those for presentation
surfaces_list_M<-list()
for(j in 1:3){
	
	# Set the parameters for XYZ set j
	XYZ<-XYZ_list[[j]]
	labels<-labels_list[match(XYZ, XYZ_list[[1]])]
	z.val<-round(quantile(dataset_plot[,XYZ[3]])[3])
	
	# Remake the surfaces sacles by the corss-surface min and max
	surfaces<-ggSurface(GAM=model_aic, data=dataset_plot, XYZ=XYZ, labels=labels, exclude=c("s(Country)", "s.1(Country)"), predict_val=predict_val, traits=c("alpha", "beta"), surf_min=min_use, surf_max=max_use, subtitle=paste0(labels[3], " = ", z.val), z.val=z.val, y.limits=y_limits[[j]])
	
	# Annotate
	surfaces[[1]]<-surfaces[[1]] + annotate("text", x = floor(min(dataset_plot[,XYZ[1]])), y = max(y_limits[[j]]), label = expression(italic(alpha)~~Males), hjust = 0, vjust = 1, size = 7)
	surfaces[[2]]<-surfaces[[2]] + annotate("text", x = floor(min(dataset_plot[,XYZ[1]])), y = max(y_limits[[j]]), label = expression(italic(beta)~~Males), hjust = 0, vjust = 1, size = 7)
	
	# Save them				
	surfaces_list_M[[j]]<-surfaces
}

################################################

# Now lets arrange all those plots
CairoPDF("figures/Figure_3.pdf", height=20, width=15)

grid.arrange(surfaces_list_F[[1]][[1]]+labs(title="A"), 
				surfaces_list_F[[2]][[1]]+labs(title="B"), 
				surfaces_list_F[[3]][[1]]+labs(title="C"), 
				surfaces_list_F[[1]][[2]]+labs(title="D"), 
				surfaces_list_F[[2]][[2]]+labs(title="E"), 
				surfaces_list_F[[3]][[2]]+labs(title="F"), 
				surfaces_list_M[[1]][[1]]+labs(title="G"), 
				surfaces_list_M[[2]][[1]]+labs(title="H"), 
				surfaces_list_M[[3]][[1]]+labs(title="I"), 
				surfaces_list_M[[1]][[2]]+labs(title="J"), 
				surfaces_list_M[[2]][[2]]+labs(title="K"), 
				surfaces_list_M[[3]][[2]]+labs(title="L"),
									layout_matrix=rbind(c(1,2,3),
														c(4,5,6),
			  					  					    c(7,8,9),
			  					  					    c(10,11,12)))

dev.off()

