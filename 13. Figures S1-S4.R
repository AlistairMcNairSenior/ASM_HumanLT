
# Clean up
rm(list=ls())

# Working directory
wd<-"/Users/alistairsenior/Dropbox (Sydney Uni)/Human lifetables and Nutrition" # Work iMac
#wd<-"/Users/asenior/Dropbox (Sydney Uni)/Human lifetables and Nutrition" # Macbook
wd<-"/Users/asenior/Dropbox (Sydney Uni)/Human lifetables and Nutrition" # Home iMac
setwd(wd)

# Load libraries
library(arm)
library(plyr)
library(ggplot2)
library(mgcv)
library(gridExtra)
library(MortalityLaws)
library(Cairo)
library(doSNOW)
library(metR)

source("scripts/0. Header_Functions.R")

# Read in the standard from Wilmoth et al 2012
standard<-read.csv("clean_data/wilmoth_standard.csv")

#################################################
#################### FIGURE S2 ##################
#################################################

# Load the full data
full_data<-read.csv("brass_data/Brass_complete_cases.csv")

# So the model with macros * time + GDP has the best fit
load("models/Complete_cases_AIC_GAMS.rdata")

# Just reorganising to put females first
GAM<-list()
GAM[[1]]<-AIC_favoured_models[[2]]
GAM[[2]]<-AIC_favoured_models[[1]]

# Sepcify the sexes
sexes<-c("Females", "Males")

##################################################
### PART 1: AVERAGE SURFACES ACROSS ALL TIMES ####
##################################################

# The only way to get a year independent surface is to calculate a surface for each year that data are available, then average over them
years<-ddply(full_data, .(Year), summarise, n.LT=length(ID))
years$weight<-years$n.LT / sum(years$n.LT)

# We will predict over the entire range of full dataset for the median of the whole dataset
dataset_plot<-full_data
med_GDP<-round(median(dataset_plot$GDP_perCapita, na.rm=T))

# Specify the desired layout for the surfaces
XYZ_list<-list()
XYZ_list[[1]]<-c("Protein.kcal", "Carbo.kcal", "Fat.kcal")
XYZ_list[[2]]<-c("Protein.kcal", "Fat.kcal", "Carbo.kcal")
XYZ_list[[3]]<-c("Carbo.kcal", "Fat.kcal", "Protein.kcal")

# Nice labels
labels_list<-c("Protein kcal/capita/day", "Carbohydrate kcal/capita/day", "Fat kcal/capita/day")

# Set some common limits (found by trial and error)
y_limits<-list()
y_limits[[1]]<-c(1000, 2250)
y_limits[[2]]<-c(500, 1600)
y_limits[[3]]<-c(500, 1600)
x_limits<-list()
x_limits[[1]]<-c(200, 525)
x_limits[[2]]<-c(150, 550)
x_limits[[3]]<-c(1000, 2500)

# # Do the avergaing using doSNOW per sex tp speed things up
# cl<-makeCluster(length(sexes), outfile="")
# registerDoSNOW(cl)

# # Do SNOW for parallel sexes
# average_list<-foreach(s = 1:length(sexes)) %dopar% {

	# # Load mgcv and MASS locally
	# require(mgcv)
	# require(MASS)

	# # List to save the average alphas and betas
	# average_alpha<-list()
	# average_beta<-list()
	
	# # Go through each year and estmate a surface for alpha and beta and calculate the weighted average
	# for(y in 1:dim(years)[1]){
		
		# # Loop for the 3 slices
		# for(j in 1:3){
			
			# # Pull out the yth year
			# predict_val<-data.frame(GDP_perCapita=med_GDP, Sex=as.factor(sexes[s]), Year=years$Year[y])
			
			# # Specify the desired layout for the surfaces j
			# XYZ<-XYZ_list[[j]]
			# x.limits<-c(floor(min(dataset_plot[,XYZ[1]])), ceiling(max(dataset_plot[,XYZ[1]])))
			# y.limits<-c(floor(min(dataset_plot[,XYZ[2]])), ceiling(max(dataset_plot[,XYZ[2]])))
			# z.val<-round(quantile(dataset_plot[,XYZ[3]])[3])
		
			# gg_surfaces<-ggSurface(GAM=GAM[[s]], data=dataset_plot, XYZ=XYZ, labels=XYZ[c(1,2)], exclude=c("s(Country)", "s.1(Country)"), predict_val=predict_val, traits=c("alpha", "beta"), subtitle=paste0(XYZ[3], " = ", z.val), z.val=z.val, x.limits=x.limits, y.limits=y.limits)
			
			# # Weight by year and add in (save on the first pass)
			# if(y == 1){
				# average_alpha[[j]]<-gg_surfaces[[1]]$data
				# average_beta[[j]]<-gg_surfaces[[2]]$data
				# average_alpha[[j]]$fit<-0
				# average_beta[[j]]$fit<-0
			# }
			# average_alpha[[j]]$fit<-average_alpha[[j]]$fit + gg_surfaces[[1]]$data$fit * years$weight[y]
			# average_beta[[j]]$fit<-average_beta[[j]]$fit + gg_surfaces[[2]]$data$fit * years$weight[y]
			
		# }
	# }
	
	# # Return the averaged lists for the sth sex
	# return(list(average_alpha, average_beta))

# }

# save(average_list, file="models/Temporal_average_surface.rdata")
load("models/Temporal_average_surface.rdata")

# Now make the PCF surfaces for alpha and beta in each sex
# List for each of the sexes
average_surfaces<-list()

# Loop for the sexes
for(s in 1:length(sexes)){

	# List for sex_s
	average_s<-list(list(), list())
	
	# Find the min and max on each
	mins<-array(NA, c(3,2))
	maxs<-mins
	for(i in 1:3){
		mins[i,1]<-min(average_list[[s]][[1]][[i]]$fit, na.rm=T)
		mins[i,2]<-min(average_list[[s]][[2]][[i]]$fit, na.rm=T)
		maxs[i,1]<-max(average_list[[s]][[1]][[i]]$fit, na.rm=T)
		maxs[i,2]<-max(average_list[[s]][[2]][[i]]$fit, na.rm=T)
	}	
	min_use<-apply(mins, 2, min)
	max_use<-apply(maxs, 2, max)

	for(k in 1:3){
		
		# Get the right labels for slice i
		labels<-labels_list[match(XYZ_list[[k]], XYZ_list[[1]])]
		
		# Pull out the right alpha for the plots
		predictions<-average_list[[s]][[1]][[k]]
		
		# This specifies the color scheme for surface	
		rgb.palette<-colorRampPalette(c("blue","cyan","yellow","red"), space="Lab", interpolate="linear")
		map<-rgb.palette(256)
		mn<-min_use[1]
		mx<-max_use[1]
		locs<-(range(predictions$fit, na.rm=TRUE) - mn) / (mx-mn) * 256
		contour_use<-signif((max(predictions$fit, na.rm=T)-min(predictions$fit, na.rm=T))/5, 1)
		
		# Plot alpha
		average_s[[1]][[k]]<-ggplot(predictions, aes(x=x, y=y)) +
					geom_raster(aes(fill=fit), show.legend=F, interpolate=F, na.rm=T) +
					scale_fill_gradientn(colors=map[locs[1]:locs[2]]) +
					geom_contour(data=predictions, aes(x=x, y=y, z=fit), na.rm=T, color="black", binwidth=contour_use) +	
					geom_label_contour(data=predictions, aes(x=x, y=y, z=fit), size=3, binwidth=contour_use, skip=1) +
					theme_bw() +
					labs(x = labels[1], y = labels[2], subtitle=paste0(labels[3], " = ", round(quantile(dataset_plot[, XYZ_list[[k]][3]])[3]))) +
					theme(axis.text=element_text(size=15), axis.title=element_text(size=15)) +
					theme(title=element_text(size=15))	+
					ylim(y_limits[[k]])	+
					xlim(x_limits[[k]]) +
					annotate("text", x = min(x_limits[[k]]), y = max(y_limits[[k]]), label = paste0("\U03B1 ", sexes[s]), hjust = 0, vjust = 1, size = 7)
					
		# Repeat for beta
		predictions<-average_list[[s]][[2]][[k]]
		
		# This specifies the color scheme for surface	
		rgb.palette<-colorRampPalette(c("blue","cyan","yellow","red"), space="Lab", interpolate="linear")
		map<-rgb.palette(256)
		mn<-min_use[2]
		mx<-max_use[2]
		contour_use<-signif((max(predictions$fit, na.rm=T)-min(predictions$fit, na.rm=T))/5, 1)
		
		# Plot alpha
		average_s[[2]][[k]]<-ggplot(predictions, aes(x=x, y=y)) +
					geom_raster(aes(fill=fit), show.legend=F, interpolate=F, na.rm=T) +
					scale_fill_gradientn(colors=map[locs[1]:locs[2]]) +
					geom_contour(data=predictions, aes(x=x, y=y, z=fit), na.rm=T, color="black", binwidth=contour_use) +	
					geom_label_contour(data=predictions, aes(x=x, y=y, z=fit), size=3, binwidth=contour_use, skip=1) +
					theme_bw() +
					labs(x = labels[1], y = labels[2], subtitle=paste0(labels[3], " = ", round(quantile(dataset_plot[, XYZ_list[[k]][3]])[3]))) +
					theme(axis.text=element_text(size=15), axis.title=element_text(size=15)) +
					theme(title=element_text(size=15))	+
					ylim(y_limits[[k]])	+
					xlim(x_limits[[k]]) +
					annotate("text", x = min(x_limits[[k]]), y = max(y_limits[[k]]), label = paste0("\U03B2 ", sexes[s]), hjust = 0, vjust = 1, size = 7)
					
										
	}
	
	# Save the sth set of surfaces
	average_surfaces[[s]]<-average_s
}

##################################################
######### PART 2: YEAR-SPECIFIC SURFACES #########
##################################################

# Get year of choice surfaces for PCF in females
year_plot<-c(1970, 1990, 2010)
mins_years<-array(NA, c(length(year_plot), 2))
maxs_years<-mins_years
labels_list<-c("Prot.", "Carb.", "Fat")

# Loop for years
for(y in 1:length(year_plot)){
	
	# Pull put dataset y	
	dataset_plot<-full_data[which(full_data$Year == year_plot[y]), ]
	predict_val<-data.frame(Year=year_plot[y], GDP_perCapita=med_GDP, Sex=as.factor("Females"))
	
	# Find the min and max set of predictions to normalise across time
	mins<-array(NA, c(3,2))
	maxs<-mins
	
	# Loop for each 	slice
	for(i in 1:3){
		
		# Get the right XYZ set
		XYZ<-XYZ_list[[i]]	
		z.val<-round(quantile(dataset_plot[,XYZ[3]])[3])
		
		# Get the right labels for slice i
		labels<-labels_list[match(XYZ_list[[i]], XYZ_list[[1]])]

		gg_surfaces<-ggSurface(GAM=GAM[[1]], data=dataset_plot, XYZ=XYZ, labels=labels, exclude=c("s(Country)", "s.1(Country)"), predict_val=predict_val, traits=c("alpha", "beta"), subtitle=paste0(labels[3], " = ", z.val), z.val=z.val, x.limits=x_limits[[i]], y.limits=y_limits[[i]], skip=2)
		
		# Save the min and max values for scaling
		mins[i,1]<-min(gg_surfaces[[1]]$data$fit)
		maxs[i,1]<-max(gg_surfaces[[1]]$data$fit)	
		mins[i,2]<-min(gg_surfaces[[2]]$data$fit)
		maxs[i,2]<-max(gg_surfaces[[2]]$data$fit)	
		
	}
	
	# Get the min for each year
	mins_years[y,]<-apply(mins, 2, min)
	maxs_years[y,]<-apply(maxs, 2, max)
	
}

# Find the min and max across years
min_use<-apply(mins_years, 2, min)
max_use<-apply(maxs_years, 2, max)

# Now refit the appropriately scaled surfaces
year_surfaces_F<-list()
for(y in 1:length(year_plot)){
	
	# Pull put dataset y	
	dataset_plot<-full_data[which(full_data$Year == year_plot[y]), ]
	predict_val<-data.frame(Year=year_plot[y], GDP_perCapita=med_GDP, Sex=as.factor("Females"))
	
	# Loop for each 	slice
	plots_y<-list()
	for(i in 1:3){
		
		# Get the right XYZ set
		XYZ<-XYZ_list[[i]]	
		z.val<-round(quantile(dataset_plot[,XYZ[3]])[3])
		
		# Get the right labels for slice i
		labels<-labels_list[match(XYZ_list[[i]], XYZ_list[[1]])]

		plots_y[[i]]<-ggSurface(GAM=GAM[[1]], data=dataset_plot, XYZ=XYZ, labels=labels, exclude=c("s(Country)", "s.1(Country)"), predict_val=predict_val, traits=c("alpha", "beta"), subtitle=paste0(labels[3], " = ", z.val), z.val=z.val, x.limits=x_limits[[i]], y.limits=y_limits[[i]], surf_min=min_use, surf_max=max_use, skip=2)
		
		# Annotate with alpoha and beta
		plots_y[[i]][[1]]<-plots_y[[i]][[1]]	 +
							annotate("text", x = min(x_limits[[i]]), y = max(y_limits[[i]]), label = "\U03B1 Females", hjust = 0, vjust = 1, size = 7)
		plots_y[[i]][[2]]<-plots_y[[i]][[2]]	 +
							annotate("text", x = min(x_limits[[i]]), y = max(y_limits[[i]]), label = "\U03B2 Females", hjust = 0, vjust = 1, size = 7)									 
		
	}
	
	# Save set y
	year_surfaces_F[[y]]<-plots_y
	
}


# Repeat for males

mins_years<-array(NA, c(length(year_plot), 2))
maxs_years<-mins_years

# Loop for years
for(y in 1:length(year_plot)){
	
	# Pull put dataset y	
	dataset_plot<-full_data[which(full_data$Year == year_plot[y]), ]
	predict_val<-data.frame(Year=year_plot[y], GDP_perCapita=med_GDP, Sex=as.factor("Males"))
	
	# Find the min and max set of predictions to normalise across time
	mins<-array(NA, c(3,2))
	maxs<-mins
	
	# Loop for each 	slice
	for(i in 1:3){
		
		# Get the right XYZ set
		XYZ<-XYZ_list[[i]]	
		z.val<-round(quantile(dataset_plot[,XYZ[3]])[3])
		
		# Get the right labels for slice i
		labels<-labels_list[match(XYZ_list[[i]], XYZ_list[[1]])]

		gg_surfaces<-ggSurface(GAM=GAM[[2]], data=dataset_plot, XYZ=XYZ, labels=labels, exclude=c("s(Country)", "s.1(Country)"), predict_val=predict_val, traits=c("alpha", "beta"), subtitle=paste0(labels[3], " = ", z.val), z.val=z.val, x.limits=x_limits[[i]], y.limits=y_limits[[i]], skip=2)
		
		# Save the min and max values for scaling
		mins[i,1]<-min(gg_surfaces[[1]]$data$fit)
		maxs[i,1]<-max(gg_surfaces[[1]]$data$fit)	
		mins[i,2]<-min(gg_surfaces[[2]]$data$fit)
		maxs[i,2]<-max(gg_surfaces[[2]]$data$fit)	
		
	}
	
	# Get the min for each year
	mins_years[y,]<-apply(mins, 2, min)
	maxs_years[y,]<-apply(maxs, 2, max)
	
}

# Find the min and max across years
min_use<-apply(mins_years, 2, min)
max_use<-apply(maxs_years, 2, max)

# Now refit the appropriately scaled surfaces
year_surfaces_M<-list()
for(y in 1:length(year_plot)){
	
	# Pull put dataset y	
	dataset_plot<-full_data[which(full_data$Year == year_plot[y]), ]
	predict_val<-data.frame(Year=year_plot[y], GDP_perCapita=med_GDP, Sex=as.factor("Males"))
	
	# Loop for each 	slice
	plots_y<-list()
	for(i in 1:3){
		
		# Get the right XYZ set
		XYZ<-XYZ_list[[i]]	
		z.val<-round(quantile(dataset_plot[,XYZ[3]])[3])
		
		# Get the right labels for slice i
		labels<-labels_list[match(XYZ_list[[i]], XYZ_list[[1]])]

		plots_y[[i]]<-ggSurface(GAM=GAM[[2]], data=dataset_plot, XYZ=XYZ, labels=labels, exclude=c("s(Country)", "s.1(Country)"), predict_val=predict_val, traits=c("alpha", "beta"), subtitle=paste0(labels[3], " = ", z.val), z.val=z.val, x.limits=x_limits[[i]], y.limits=y_limits[[i]], surf_min=min_use, surf_max=max_use, skip=2)
		
		# Annotate with alpoha and beta
		plots_y[[i]][[1]]<-plots_y[[i]][[1]]	 +
							annotate("text", x = min(x_limits[[k]]), y = max(y_limits[[k]]), label = "\U03B1 Males", hjust = 0, vjust = 1, size = 7)
		plots_y[[i]][[2]]<-plots_y[[i]][[2]]	 +
							annotate("text", x = min(x_limits[[k]]), y = max(y_limits[[k]]), label = "\U03B2 Males", hjust = 0, vjust = 1, size = 7)
		
	}
	
	# Save set y
	year_surfaces_M[[y]]<-plots_y
	
}


##################################################
########### PART 3: CREATE A TIMELINE ############
##################################################

timeline<-data.frame(time=c(1970, 1990, 2010), y=0)
timeline$lower<--1
timeline$upper<-0

time_course<-ggplot(timeline, aes(x=time, y=y)) +
	geom_line(color="purple") +
	theme(axis.line=element_blank(),axis.text.x=element_blank(),
          axis.text.y=element_blank(),axis.ticks=element_blank(),
          axis.title.x=element_blank(),
          axis.title.y=element_blank(),legend.position="none",
          panel.background=element_blank(),panel.border=element_blank(),panel.grid.major=element_blank(),
          panel.grid.minor=element_blank(),plot.background=element_blank()) +
     geom_errorbar(aes(ymin=lower, ymax=upper), width=0.01, color="purple") +
     ylim(c(-3,1)) +
     xlim(c(1964, 2015)) +
     annotate("text", x=timeline$time, y=timeline$lower-0.75, label=timeline$time, size=10)
     
    


##################################################
################# PART 4: ARRANGE ################
##################################################


# Set the array for layout
mylayout<-array(NA, c(5, 20))
mylayout[1,]<-c(1,1,2,2,3,3, NA, 7,7,8,8,9,9, NA, 13,13,14,14,15,15)
mylayout[2,]<-c(1,1,2,2,3,3, NA, 7,7,8,8,9,9, NA, 13,13,14,14,15,15)
mylayout[3,]<-19
mylayout[4,]<-c(4,4,5,5,6,6, NA, 10,10,11,11,12,12, NA, 16,16,17,17,18,18)
mylayout[5,]<-c(4,4,5,5,6,6, NA, 10,10,11,11,12,12, NA, 16,16,17,17,18,18)


# Now lets arrange all those plots
CairoPDF("figures/Figure_S1.pdf", height=10, width=35)

grid.arrange(year_surfaces_F[[1]][[1]][[1]],
			 year_surfaces_F[[1]][[2]][[1]],
			 year_surfaces_F[[1]][[3]][[1]],
			 year_surfaces_F[[1]][[1]][[2]],
			 year_surfaces_F[[1]][[2]][[2]],
			 year_surfaces_F[[1]][[3]][[2]],
			 year_surfaces_F[[2]][[1]][[1]],
			 year_surfaces_F[[2]][[2]][[1]],
			 year_surfaces_F[[2]][[3]][[1]],
			 year_surfaces_F[[2]][[1]][[2]],
			 year_surfaces_F[[2]][[2]][[2]],
			 year_surfaces_F[[2]][[3]][[2]],
			 year_surfaces_F[[3]][[1]][[1]],
			 year_surfaces_F[[3]][[2]][[1]],
			 year_surfaces_F[[3]][[3]][[1]],
			 year_surfaces_F[[3]][[1]][[2]],
			 year_surfaces_F[[3]][[2]][[2]],
			 year_surfaces_F[[3]][[3]][[2]],
			 time_course,
				layout_matrix=mylayout)

dev.off()


# Now lets arrange all those plots
CairoPDF("figures/Figure_S2.pdf", height=10, width=35)

grid.arrange(year_surfaces_M[[1]][[1]][[1]],
			 year_surfaces_M[[1]][[2]][[1]],
			 year_surfaces_M[[1]][[3]][[1]],
			 year_surfaces_M[[1]][[1]][[2]],
			 year_surfaces_M[[1]][[2]][[2]],
			 year_surfaces_M[[1]][[3]][[2]],
			 year_surfaces_M[[2]][[1]][[1]],
			 year_surfaces_M[[2]][[2]][[1]],
			 year_surfaces_M[[2]][[3]][[1]],
			 year_surfaces_M[[2]][[1]][[2]],
			 year_surfaces_M[[2]][[2]][[2]],
			 year_surfaces_M[[2]][[3]][[2]],
			 year_surfaces_M[[3]][[1]][[1]],
			 year_surfaces_M[[3]][[2]][[1]],
			 year_surfaces_M[[3]][[3]][[1]],
			 year_surfaces_M[[3]][[1]][[2]],
			 year_surfaces_M[[3]][[2]][[2]],
			 year_surfaces_M[[3]][[3]][[2]],
			 time_course,
				layout_matrix=mylayout)

dev.off()


# Now lets arrange all those plots
CairoPDF("figures/Figure_S3.pdf", height=10, width=15)

grid.arrange(average_surfaces[[1]][[1]][[1]],
			 average_surfaces[[1]][[1]][[2]],
			 average_surfaces[[1]][[1]][[3]],
			 average_surfaces[[1]][[2]][[1]],
			 average_surfaces[[1]][[2]][[2]],
			 average_surfaces[[1]][[2]][[3]],
				layout_matrix=rbind(c(1,2,3),c(4,5,6)))

dev.off()


# Now lets arrange all those plots
CairoPDF("figures/Figure_S4.pdf", height=10, width=15)

grid.arrange(average_surfaces[[2]][[1]][[1]],
			 average_surfaces[[2]][[1]][[2]],
			 average_surfaces[[2]][[1]][[3]],
			 average_surfaces[[2]][[2]][[1]],
			 average_surfaces[[2]][[2]][[2]],
			 average_surfaces[[2]][[2]][[3]],
				layout_matrix=rbind(c(1,2,3),c(4,5,6)))

dev.off()