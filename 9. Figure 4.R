
# Clean up
rm(list=ls())

# Working directory
wd<-"/Users/alistairsenior/Dropbox/Human lifetables and Nutrition"
wd<-"/Users/alistairsenior/Dropbox (Sydney Uni)/Human lifetables and Nutrition"
wd<-"/Users/asenior/Dropbox (Sydney Uni)/Human lifetables and Nutrition"
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
library(MASS)
library(metR)

source("0.Header_Functions.R")

# Read in the standard from Wilmoth et al 2012
standard<-read.csv("wilmoth_standard.csv")

#################################################
#################### FIGURE 4 ###################
#################################################

# # # Load the full data
full_data<-read.csv("Brass_complete_cases.csv")

# # Load the AIC favoured model which had Macro*time + GDP
# # So the model with macros * time has the best fit
load("GAM_int_GDP.rdata")

# Year and dataset to use for plotting
year_plot<-2016
dataset_plot<-full_data[which(full_data$Year == year_plot), ]
med_GDP<-round(median(dataset_plot$GDP_perCapita, na.rm=T))

# We will run each sex over 1 core using doSNOW
sexes<-c("Females", "Males")

# For each age class we will make 6 surfaces, PC surfaces for min/max with the F axis held at low med high
surface_order<-list()
surface_order[[1]]<-c("Protein.kcal", "Carbo.kcal", "Fat.kcal")
surface_order[[2]]<-c("Protein.kcal", "Carbo.kcal", "Fat.kcal")
surface_order[[3]]<-c("Protein.kcal", "Carbo.kcal", "Fat.kcal")

labels<-c("Protein kcal/cap/day", "Carbohydrate kcal/cap/day", "Fat kcal/cap/day")

slice_3<-quantile(dataset_plot$Fat.kcal)[c(2,3,4)]

# This specifies the color scheme for surface	
rgb.palette<-colorRampPalette(c("blue","cyan","yellow","red"), space="Lab", interpolate="linear")
map<-rgb.palette(256)

# Which ages do we want to plot
ages_plot<-c(5, 50, 75)

# List to hold the plots
plots_list<-list()

for(s in 1:length(sexes)){
	
	# Ages_list to hold the plots for all the ages
	ages_list<-list()
	
	# For each age class
	for(a in 1:length(ages_plot)){
		
		# List to hold the 3 slices					
		slices_list<-list()
		
		# Objects for the min and max on each surface
		surf_min<-NA
		surf_max<-NA
			
			# For each PCF surface on the minimum
			for(i in 1:length(surface_order)){
				
				# The order for the ith surface
				surface_i<-surface_order[[i]]
				
				# Lets work out the protein, carb and fat supply values for which we will get lifetables
				energy_1<-seq(floor(min(dataset_plot[,surface_i[1]])), ceiling(max(dataset_plot[,surface_i[1]])), length=101)
				energy_2<-seq(floor(min(dataset_plot[,surface_i[2]])), ceiling(max(dataset_plot[,surface_i[2]])), length=101)
				energy_3<-slice_3[i]
				
				# Get all of the combinations
				predict_values<-expand.grid(energy_1, energy_2, energy_3)
				names(predict_values)<-surface_order[[i]]			
				
				# Cull out any unobserved energy supplies
				in.poly<-as.numeric(inhull(predict_values[,c("Protein.kcal", "Carbo.kcal", "Fat.kcal")], dataset_plot[,c("Protein.kcal", "Carbo.kcal", "Fat.kcal")]) != -1)
				predict_values<-predict_values[which(in.poly==1),]
				
				# Now convert each value to the ASM
				predict_values$qx<-NA
				for(v in 1:dim(predict_values)[1]){
					predict_values$qx[v]<-as.numeric(convert_brass_general(predict_values[v, surface_order[[1]]], sex=sexes[s], age=ages_plot[a], GDP=med_GDP, year=year_plot, standard=standard, GAM=GAM))
				}
				
				# Save the ith set of predicted values
				slices_list[[i]]<-predict_values
				surf_min<-min(surf_min, predict_values$qx, na.rm=T)
				surf_max<-max(surf_max, predict_values$qx, na.rm=T)
			}	
			
			# Re loop through to actually plot
			for(i in 1:length(surface_order)){
			
				# Make the surface
				predict_values<-slices_list[[i]]
				locs<-(range(predict_values$qx, na.rm=TRUE) - surf_min) / (surf_max-surf_min) * 256
				contour_use<-signif((max(predict_values$qx, na.rm=T)-min(predict_values$qx, na.rm=T))/5, 1)
				names(predict_values)[c(1,2)]<-c("x", "y")
				plot<-ggplot(predict_values, aes(x=x, y=y)) +
						geom_raster(aes(fill=qx), show.legend=F, interpolate=F, na.rm=T) +
						scale_fill_gradientn(colors=map[locs[1]:locs[2]]) +
						geom_contour(data=predict_values, aes(x=x, y=y, z=qx), na.rm=T, color="black", binwidth=contour_use) +	
						geom_label_contour(data=predict_values, aes(x=x, y=y, z=qx), size=3, binwidth=contour_use, skip=1) +
						theme_bw() +
						labs(x = labels[1], y = labels[2], subtitle=paste0(labels[3], " = ", round(slice_3[i]), ", Age = ", ages_plot[a], " to ", ages_plot[a]+5)) +
						theme(axis.text=element_text(size=15), axis.title=element_text(size=15)) +
						theme(title=element_text(size=15)) + 
						xlim(range(energy_1)) +
						ylim(range(energy_2)) 
											
				slices_list[[i]]<-plot	
	
			}
		
		# Add the surfaces to the figure	
		ages_list[[a]]<-slices_list	

	}

	# Save the slices for the sth sex
	plots_list[[s]]<-ages_list

}


pdf("Figure 4.pdf", height=15, width=15)


grid.arrange(plots_list[[1]][[1]][[1]]+labs(title="A"),
						plots_list[[1]][[1]][[2]]+labs(title="B"), 
						plots_list[[1]][[1]][[3]]+labs(title="C"),
						plots_list[[1]][[2]][[1]]+labs(title="D"),
						plots_list[[1]][[2]][[2]]+labs(title="E"), 
						plots_list[[1]][[2]][[3]]+labs(title="F"),
						plots_list[[1]][[3]][[1]]+labs(title="G"),
						plots_list[[1]][[3]][[2]]+labs(title="H"), 
						plots_list[[1]][[3]][[3]]+labs(title="I"),
#						plots_list[[1]][[4]][[1]]+labs(title="J"),
#						plots_list[[1]][[4]][[2]]+labs(title="K"), 
#						plots_list[[1]][[4]][[3]]+labs(title="L"),
						layout_matrix=rbind(c(1,2,3),
											c(4,5,6),
											c(7,8,9)))
#											c(10,11,12)))

dev.off()
