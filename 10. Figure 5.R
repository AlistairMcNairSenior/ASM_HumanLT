
# Clean up
rm(list=ls())

# Working directory
wd<-"/Users/alistairsenior/Dropbox (Sydney Uni)/Human lifetables and Nutrition" # work iMac
#wd<-"/Users/asenior/Dropbox (Sydney Uni)/Human lifetables and Nutrition"
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

source("scripts/0. Header_Functions.R")

# Read in the standard from Wilmoth et al 2012
standard<-read.csv("clean_data/wilmoth_standard.csv")

#################################################
#################### FIGURE 5 ###################
#################################################

# # # Load the full data
full_data<-read.csv("brass_data/Brass_complete_cases.csv")

# # Load the AIC favoured model which had Macro*time + GDP
# # So the model with macros * time has the best fit
load("models/Complete_cases_AIC_GAMS.rdata")

# Rearrange the GAMs to put females first
GAM<-list()
GAM[[1]]<-AIC_favoured_models[[2]]
GAM[[2]]<-AIC_favoured_models[[1]]

# Year and dataset to use for plotting
year_plot<-2016
dataset_plot<-full_data[which(full_data$Year == year_plot), ]
med_GDP<-round(median(dataset_plot$GDP_perCapita, na.rm=T))

# We will run with x starting values for the optimiser based on mv-normal distribution, and take the best of all observed
set.seed(123)
# Put together in a list
init_list<-list()
for(i in 1:500){
	init_list[[i]]<-mvrnorm(1, apply(dataset_plot[,c(4:6)], 2, mean), cov(dataset_plot[,c(4:6)]))
}

# We will run each sex over 1 core using doSNOW
sexes<-c("Females", "Males")

# Set the cluster
cl<-makeCluster(length(sexes), outfile="")
registerDoSNOW(cl)

# List to hold the predictions for each sex
predictions_list<-foreach(j = 1:length(sexes)) %dopar% {
	
	# Create some variables to hold the max/min starting with the ages that we will search through, note I am dropping 105 and 110 here as preliminary analyses suggested that there are many many PCFs that maximise/minimise these
	ages<-c(0, 1, seq(5, 100, 5))
	max_expect<-list()
	min_expect<-max_expect
	
	# Run the optimiser for each age class
	for(i in 1:length(ages)){
		
		# Go through the initiating list
		res_temp_min<-data.frame(val=rep(NA, length(init_list)), P=NA, C=NA, F=NA) 
		res_temp_max<-res_temp_min
		
		# Now optimise for each starting value
		for(o in 1:length(init_list)){
			
			# Pull out the oth set of starting values
			pars<-init_list[[o]]
			
			# First get the combination estimated to mimise qx at age i and then repeat to get the max within intial values o
			optim_results<-optim(par=pars, fn=convert_brass, sex=sexes[j], age=ages[i], GDP=med_GDP, year=year_plot, data=dataset_plot, standard=standard, GAM=GAM[[j]])
			res_temp_min[o,]<-c(optim_results$value, optim_results$par)
			
			optim_results<-optim(par=pars, fn=convert_brass, sex=sexes[j], age=ages[i], GDP=med_GDP, year=year_plot, data=dataset_plot, standard=standard, GAM=GAM[[j]], aim="maximise")
			res_temp_max[o,]<-c(optim_results$value, optim_results$par)
		}
			
		# Save the results in to the list
		max_expect[[i]]<-res_temp_max
		min_expect[[i]]<-res_temp_min	
		
	}
	
	# Save the list to be returned
	names(max_expect)<-ages
	names(min_expect)<-ages
	predictions<-list(min_expect, max_expect)
	return(predictions)

}

# Save and reload (takes a while to get these data)
save(predictions_list, file="models/Predictions.Rdata")
load("models/Predictions.Rdata")

# Go through each in the list and find the minimal value
list_new<-list()
for(s in 1:length(predictions_list)){
	list_s<-predictions_list[[s]]
	list_temp<-list()
	for(m in 1:length(list_s)){
		list_m<-list_s[[m]]
		data_temp<-data.frame(ages=as.numeric(names(list_m)), qx=NA, Protein.kcal=NA, Carbo.kcal=NA, Fat.kcal=NA)
		for(a in 1:length(list_m)){
			list_a<-list_m[[a]]
			hit<-list_a[which(list_a$val == min(list_a$val)),]
			if(m == 2){
				# Inverting the max values
				hit$val<-hit$val*-1
			}
			data_temp[a,-1]<-hit
		}
		list_temp[[m]]<-data_temp
	}
	names(list_temp)<-c("Min", "Max")
	list_new[[s]]<-list_temp	
}
names(list_new)<-sexes
predictions_list<-list_new

# Go through each and plot the composition and total energy for each sex minimising
min_list<-list()

for(i in 1:length(sexes)){
	
	# Pull out the ith data set
	data_i2<-predictions_list[[i]][[1]]
	data_i2$type<-"Min"
		
	# Add a categorical predictor for each row to color for age class
	data_i2$age<-standard$Age[c(1:dim(data_i2)[1])]
	
	# Limit for the % and total kj axis
	ylim_range=c(5, 75)
	zlim_max<-4000
	margin<-c(0.5,0.8,0.5,0.8)
	
	# Caluclate as a %
	data_i2$total<-data_i2$Fat.kcal + data_i2$Protein.kcal + data_i2$Carbo.kcal
	data_i2$Fat.kcal<-data_i2$Fat.kcal / data_i2$total * 100
	data_i2$Carbo.kcal<-data_i2$Carbo.kcal / data_i2$total * 100
	data_i2$Protein.kcal<-data_i2$Protein.kcal / data_i2$total * 100
	data_i2$total.plot<-data_i2$total / zlim_max * ylim_range[2]
	
				
	# Plot the min 	
	min_list[[i]]<-ggplot(data_i2, aes(x=ages, y=Fat.kcal, color="Fat")) + 
			geom_point(size=3) + 
			geom_line() +
			geom_point(aes(x=ages, y=Carbo.kcal, color="Carbohydrate"), size=3) +
			geom_line(aes(x=ages, y=Carbo.kcal, color="Carbohydrate")) +
			geom_point(aes(x=ages, y=Protein.kcal, color="Protein"), size=3) +
			geom_line(aes(x=ages, y=Protein.kcal, color="Protein")) +
			geom_point(aes(x=ages, y=total.plot, color="Total Energy"), size=3) +
			geom_line(aes(x=ages, y=total.plot, color="Total Energy")) +	
			labs(x="Age", y="% of Energy", subtitle=paste0("Lowest Mortality in ", sexes[i]), color=element_blank()) + 
			theme_bw() +
			theme(axis.line = element_line(colour = "black")) + 
			theme(axis.text=element_text(size=15), axis.title=element_text(size=15), plot.subtitle=element_text(size=15)) +
			scale_y_continuous(limits=ylim_range, sec.axis = sec_axis(~. /ylim_range[2]*zlim_max, name="Total Energy (kcal/cap/day)")) +
			theme(plot.margin = unit(margin, "cm"))
	
	if(i == 1){
		min_list[[i]]<-min_list[[i]] + 
			theme(legend.position=c(0.5, 0.94), legend.text=element_text(size=13), legend.background=element_blank(), legend.direction="horizontal")
	}else{
		min_list[[i]]<-min_list[[i]] + 
			theme(legend.position="none")
	}
				
}


# Lastly lets look at manipulating the median and seeing how that affects ASM
ratio_list<-list()

# The median supply refernece diet and check it is actually within what is observed
median_diet.kcal<-round(c(median(dataset_plot$Protein.kcal), median(dataset_plot$Carbo.kcal), median(dataset_plot$Fat.kcal)))
inhull(array(median_diet.kcal, c(1,3)), dataset_plot[,c("Protein.kcal", "Carbo.kcal", "Fat.kcal")])

# Estimate the probability of ASM at the median supply
delta_data<-data.frame(age=standard$Age[c(1:dim(data_i2)[1])], qx_f=NA, qx_m=NA)
for(x in 1:dim(delta_data)[1]){
	delta_data$qx_f[x]<-convert_brass_general(x=median_diet.kcal, sex=sexes[1], age=delta_data$age[x], GDP=med_GDP, year=year_plot, standard=standard, GAM=GAM[[1]], stat="qx")
	delta_data$qx_m[x]<-convert_brass_general(x=median_diet.kcal, sex=sexes[2], age=delta_data$age[x], GDP=med_GDP, year=year_plot, standard=standard, GAM=GAM[[2]], stat="qx")
}


# We will first +/- X kcal across the board, then +/1 X at a time to each macornutrient
X<-0.1
manipulations<-array(0, c(3,3))
diag(manipulations)<-X
manipulations<-rbind(rep(0.1, 3), manipulations)
labels<-c("Total Energy", "Protein", "Carbohydrate", "Fat")
color_pallette<-c("blue", "red", "lightskyblue", "rosybrown2")
text_add<-paste0(c("Protein = ", "Carbohydrate = ", "Fat = "), median_diet.kcal)
text_add<-c("Median Supply (kcal/cap/day)", text_add)
text_x<-80
text_y<-seq(1.15, 1.09, length=4)
text_data<-data.frame(x=text_x, y=text_y, label=text_add)

for(m in 1:dim(manipulations)[1]){
	
	# Add the energy
	delta_m<-delta_data
	for(x in 1:dim(delta_data)[1]){
		delta_m$qx_f[x]<-convert_brass_general(x=median_diet.kcal * (1+manipulations[m,]), sex=sexes[1], age=delta_data$age[x], GDP=med_GDP, year=year_plot, standard=standard, GAM=GAM[[1]], stat="qx")
		delta_m$qx_m[x]<-convert_brass_general(x=median_diet.kcal * (1+manipulations[m,]), sex=sexes[2], age=delta_data$age[x], GDP=med_GDP, year=year_plot, standard=standard, GAM=GAM[[2]], stat="qx")
	}
	delta_m$type<-paste0("+", X*100, "%")
	
	# Subtract the energy
	delta_m2<-delta_data
	for(x in 1:dim(delta_data)[1]){
		delta_m2$qx_f[x]<-convert_brass_general(x=median_diet.kcal * (1-manipulations[m,]), sex=sexes[1], age=delta_data$age[x], GDP=med_GDP, year=year_plot, standard=standard, GAM=GAM[[1]], stat="qx")
		delta_m2$qx_m[x]<-convert_brass_general(x=median_diet.kcal * (1-manipulations[m,]), sex=sexes[2], age=delta_data$age[x], GDP=med_GDP, year=year_plot, standard=standard, GAM=GAM[[2]], stat="qx")
	}
	delta_m2$type<-paste0("-", X*100, "%")
	
	# Bind together
	delta_m<-rbind(delta_m, delta_m2)
	
	# Make relative
	delta_m$qx_f<-delta_m$qx_f/delta_data$qx_f
	delta_m$qx_m<-delta_m$qx_m/delta_data$qx_m
	
	margin<-c(0.5,1.6,0.5,0.8)
	ratio<-ggplot(delta_m, aes(x=age, y=qx_f, color=paste0("Females " , delta_m$type))) +
					geom_point(size=3) + 
					geom_line() +
					geom_point(aes(x=age, y=qx_m, color=paste0("Males ", delta_m$type)), size=3) + 
					geom_line(aes(x=age, y=qx_m, color=paste0("Males ", delta_m$type))) + 
					labs(y = expression(italic(q[x])~~Relative~~to~~Median), x=expression(Age~~italic(x)), subtitle=labels[m], color=element_blank()) +
					theme_bw() +
					ylim(c(0.9, 1.15)) +
					theme(axis.text=element_text(size=15), axis.title=element_text(size=15) , plot.subtitle=element_text(size=15)) +
					geom_abline(intercept=1, slope=0) +
					scale_color_manual(values=color_pallette) +
					theme(plot.margin = unit(margin, "cm"))
						
	# Add a legend on pass 1
	if(m == 1){
		ratio<-ratio + 
				theme(legend.position=c(0.2, 0.18), legend.text=element_text(size=15), legend.background=element_blank()) + 
				annotate("text", x=text_data$x, y=text_data$y, label=text_data$label)
	}else{
		ratio<-ratio + theme(legend.position="none")
	}
	
	ratio_list[[m]]<-ratio	

}
	 	
			
# Now lets arrange all those plots
pdf("figures/Figure_5.pdf", height=15, width=15)

grid.arrange(min_list[[1]]+labs(title="A"),
				min_list[[2]]+labs(title="B"),
				ratio_list[[1]]+labs(title="C"),
				ratio_list[[2]]+labs(title="D"),
				ratio_list[[3]]+labs(title="E"),
				ratio_list[[4]]+labs(title="F"),
									layout_matrix=rbind(c(1,2),
														c(3,4),
														c(5,6)))

dev.off()

