
# Clean up
rm(list=ls())

# Working directory
wd<-"/Users/alistairsenior/Dropbox (Sydney Uni)/Human lifetables and Nutrition" # Work iMac
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
#################### FIGURE S6 ##################
#################################################

# Load the full data
full_data<-read.csv("brass_data/Brass_subset.csv")

# So the model with macros * time + GDP has the best fit
load("models/Subset_AIC_GAMS.rdata")
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
y_limits[[1]]<-c(1300, 1800)
y_limits[[2]]<-c(700, 1500)
y_limits[[3]]<-c(1100, 1500)

# List of axes labels
labels_list<-c("Protein kcal/capita/day", "Carbohydrate kcal/capita/day", "Fat kcal/capita/day")

# Find the min and max set of predictions
mins<-array(NA, c(3,2))
maxs<-mins
surfaces_list<-list()
for(j in 1:3){
	
	# Set the parameters for XYZ set j
	XYZ<-XYZ_list[[j]]
	z.val<-round(quantile(dataset_plot[,XYZ[3]])[3])
	
	# Generate alpha and beta surfaces for slice i	
	surfaces_list[[j]]<-ggSurface(GAM=model_aic, data=dataset_plot, XYZ=XYZ, labels=XYZ, exclude=c("s(Country)", "s.1(Country)"), predict_val=predict_val, traits=c("alpha", "beta"), z.val=z.val, y.limits=y_limits[[j]])
	
}

# Now for each pair of surfaces (alpha and beta) we need to convert the alpha and beta values to LTs and estimate life expectancy at birth (e0)
predictions_list<-list()

# To begin with do the conversion and find the mins and maxs across all surfaces
mins<-array(NA, c(3,1))
maxs<-mins
for(k in 1:3){
	
	# Pull out the alpha and betas from the plots
	predictions<-surfaces_list[[k]][[1]]$data
	predictions$alpha<-predictions$fit
	predictions$beta<-surfaces_list[[k]][[2]]$data$fit
	
	# Add in an NA columns to hold data
	predictions$e0<-NA

	for(j in 1:dim(predictions)[1]){
		
		# Which sex are we making predictions for
		convert1<-c("Males", "Females")
		sex_used<-match(predict_val$Sex, convert1)
				
		# Get the jth lifetable
		tag<-which(standard$Age == 5)
		ls5<-standard[tag, grep("lx", names(standard))[sex_used]]
		l5<-invlogit(predictions$alpha[j] + predictions$beta[j] * logit(ls5))
		tag<-which(standard$Age == 60)
		ls60<-standard[tag, grep("lx", names(standard))[sex_used]]
		l60<-invlogit(predictions$alpha[j] + predictions$beta[j] * logit(ls60))
		l.x.<-invlogit(predictions$alpha[j] + predictions$beta[j] * logit(standard[, grep("lx", names(standard))[sex_used]]) + standard[, grep("gamma", names(standard))[sex_used]] * (1 - (logit(l5)/logit(ls5))) + standard[, grep("theta", names(standard))[sex_used]] * (1 - (logit(l60)/logit(ls60))))
		
		# Get the full lifetable
		convert2<-c("male", "female")
		LT_j<-LifeTable(x=standard$Age, lx=l.x., sex=convert2[sex_used])$lt
		predictions$e0[j]<-LT_j$ex[which(LT_j$x == 0)]
	}
	
	# Save the predictions
	predictions_list[[k]]<-predictions
	mins[k,1]<-min(predictions$e0)
	maxs[k,1]<-max(predictions$e0)
						
}

# Find the min and max
min_use<-apply(mins, 2, min)
max_use<-apply(maxs, 2, max)

# Now make the plots
expectancy_list_F<-list()

# To begin with do the conversion to find the mins and maxs across all surfaces
for(k in 1:3){
	
	# Pull out the alpha and betas from the plots
	predictions<-predictions_list[[k]]
	
	# List to hold the kth set of life expectancies
	expectancy_k<-list()
	
	# Get the right labels for slice i
	labels<-labels_list[match(XYZ_list[[k]], XYZ_list[[1]])]
	
	# This specifies the color scheme for surface	
	rgb.palette<-colorRampPalette(c("blue","cyan","yellow","red"), space="Lab", interpolate="linear")
	map<-rgb.palette(256)
	mn<-min_use[1]
	mx<-max_use[1]
	locs<-(range(predictions$e0, na.rm=TRUE) - mn) / (mx-mn) * 256
	contour_use<-signif((max(predictions$e0, na.rm=T)-min(predictions$e0, na.rm=T))/5, 1)
	
	expectancy_list_F[[k]]<-ggplot(predictions, aes(x=x, y=y)) +
				geom_raster(aes(fill=e0), show.legend=F, interpolate=F, na.rm=T) +
				scale_fill_gradientn(colors=map[locs[1]:locs[2]]) +
				geom_contour(data=predictions, aes(x=x, y=y, z=e0), na.rm=T, color="black", binwidth=contour_use) +	
				geom_label_contour(data=predictions, aes(x=x, y=y, z=e0), size=3, binwidth=contour_use, skip=1) +
				theme_bw() +
				labs(x = labels[1], y = labels[2], subtitle=paste0(labels[3], " = ", round(quantile(dataset_plot[, XYZ_list[[k]][3]])[3]))) +
				theme(axis.text=element_text(size=15), axis.title=element_text(size=15)) +
				theme(title=element_text(size=15))	+
				ylim(y_limits[[k]])	+
				annotate("text", x = floor(min(predictions[,1])), y = max(y_limits[[k]]), label = expression(italic(e[0])~~Females), hjust = 0, vjust = 1, size = 7)
									
}

#######################################

# Repeat for males
model_aic<-AIC_favoured_models[[1]]
predict_val<-data.frame(Year=year_plot, GDP_perCapita=med_GDP, Sex=as.factor("Males"))

# Find the min and max set of predictions
mins<-array(NA, c(3,2))
maxs<-mins
surfaces_list<-list()
for(j in 1:3){
	
	# Set the parameters for XYZ set j
	XYZ<-XYZ_list[[j]]
	z.val<-round(quantile(dataset_plot[,XYZ[3]])[3])
	
	# Generate alpha and beta surfaces for slice i	
	surfaces_list[[j]]<-ggSurface(GAM=model_aic, data=dataset_plot, XYZ=XYZ, labels=XYZ, exclude=c("s(Country)", "s.1(Country)"), predict_val=predict_val, traits=c("alpha", "beta"), z.val=z.val, y.limits=y_limits[[j]])
	
}

# Now for each pair of surfaces (alpha and beta) we need to convert the alpha and beta values to LTs and estimate life expectancy at birth (e0)
predictions_list<-list()

# To begin with do the conversion and find the mins and maxs across all surfaces
mins<-array(NA, c(3,1))
maxs<-mins
for(k in 1:3){
	
	# Pull out the alpha and betas from the plots
	predictions<-surfaces_list[[k]][[1]]$data
	predictions$alpha<-predictions$fit
	predictions$beta<-surfaces_list[[k]][[2]]$data$fit
	
	# Add in an NA columns to hold data
	predictions$e0<-NA

	for(j in 1:dim(predictions)[1]){
		
		# Which sex are we making predictions for
		convert1<-c("Males", "Females")
		sex_used<-match(predict_val$Sex, convert1)
				
		# Get the jth lifetable
		tag<-which(standard$Age == 5)
		ls5<-standard[tag, grep("lx", names(standard))[sex_used]]
		l5<-invlogit(predictions$alpha[j] + predictions$beta[j] * logit(ls5))
		tag<-which(standard$Age == 60)
		ls60<-standard[tag, grep("lx", names(standard))[sex_used]]
		l60<-invlogit(predictions$alpha[j] + predictions$beta[j] * logit(ls60))
		l.x.<-invlogit(predictions$alpha[j] + predictions$beta[j] * logit(standard[, grep("lx", names(standard))[sex_used]]) + standard[, grep("gamma", names(standard))[sex_used]] * (1 - (logit(l5)/logit(ls5))) + standard[, grep("theta", names(standard))[sex_used]] * (1 - (logit(l60)/logit(ls60))))
		
		# Get the full lifetable
		convert2<-c("male", "female")
		LT_j<-LifeTable(x=standard$Age, lx=l.x., sex=convert2[sex_used])$lt
		predictions$e0[j]<-LT_j$ex[which(LT_j$x == 0)]
	}
	
	# Save the predictions
	predictions_list[[k]]<-predictions
	mins[k,1]<-min(predictions$e0)
	maxs[k,1]<-max(predictions$e0)
						
}

# Find the min and max
min_use<-apply(mins, 2, min)
max_use<-apply(maxs, 2, max)

# Now make the plots
expectancy_list_M<-list()

# To begin with do the conversion to find the mins and maxs across all surfaces
for(k in 1:3){
	
	# Pull out the alpha and betas from the plots
	predictions<-predictions_list[[k]]
	
	# List to hold the kth set of life expectancies
	expectancy_k<-list()
	
	# Get the right labels for slice i
	labels<-labels_list[match(XYZ_list[[k]], XYZ_list[[1]])]
	
	# This specifies the color scheme for surface	
	rgb.palette<-colorRampPalette(c("blue","cyan","yellow","red"), space="Lab", interpolate="linear")
	map<-rgb.palette(256)
	mn<-min_use[1]
	mx<-max_use[1]
	locs<-(range(predictions$e0, na.rm=TRUE) - mn) / (mx-mn) * 256
	contour_use<-signif((max(predictions$e0, na.rm=T)-min(predictions$e0, na.rm=T))/5, 1)
	
	expectancy_list_M[[k]]<-ggplot(predictions, aes(x=x, y=y)) +
				geom_raster(aes(fill=e0), show.legend=F, interpolate=F, na.rm=T) +
				scale_fill_gradientn(colors=map[locs[1]:locs[2]]) +
				geom_contour(data=predictions, aes(x=x, y=y, z=e0), na.rm=T, color="black", binwidth=contour_use) +	
				geom_label_contour(data=predictions, aes(x=x, y=y, z=e0), size=3, binwidth=contour_use, skip=1) +
				theme_bw() +
				labs(x = labels[1], y = labels[2], subtitle=paste0(labels[3], " = ", round(quantile(dataset_plot[, XYZ_list[[k]][3]])[3]))) +
				theme(axis.text=element_text(size=15), axis.title=element_text(size=15)) +
				theme(title=element_text(size=15))	+
				ylim(y_limits[[k]])	+
				annotate("text", x = floor(min(predictions[,1])), y = max(y_limits[[k]]), label = expression(italic(e[0])~~Females), hjust = 0, vjust = 1, size = 7)
									
}

##################################

# Now lets arrange all those plots
CairoPDF("figures/Figure_S6.pdf", height=10, width=15)

grid.arrange(expectancy_list_F[[1]]+labs(title="A"), 
				expectancy_list_F[[2]]+labs(title="B"), 
				expectancy_list_F[[3]]+labs(title="C"),
				expectancy_list_M[[1]]+labs(title="D"), 
				expectancy_list_M[[2]]+labs(title="E"), 
				expectancy_list_M[[3]]+labs(title="F"),	 
								layout_matrix=rbind(c(1,2,3),
													c(4,5,6)))

dev.off()

