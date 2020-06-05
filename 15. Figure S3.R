
# Clean up
rm(list=ls())

# Working directory
wd<-"/Users/alistairsenior/Dropbox/Human lifetables and Nutrition"
wd<-"/Users/alistairsenior/Dropbox (Sydney Uni)/Human lifetables and Nutrition"
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

source("0.Header_Functions.R")

# Read in the standard from Wilmoth et al 2012
standard<-read.csv("wilmoth_standard.csv")

#################################################
#################### FIGURE S3 ##################
#################################################

# # # Load the full data
full_data<-read.csv("Brass_imputed.csv")

# # Load the AIC favoured model which had Macro*time + GDP
# # So the model with macros * time + GDP has the best fit
load("GAM_int_GDP_imp.rdata")


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

# # List to hold the predictions for each sex
# predictions_list<-foreach(j = 1:length(sexes)) %dopar% {
	
	# # Create some variables to hold the max/min starting with the ages that we will search through, note I am dropping 105 and 110 here as preliminary analyses suggested that there are many many PCFs that maximise/minimise these
	# ages<-c(0, 1, seq(5, 100, 5))
	# max_expect<-list()
	# min_expect<-max_expect
	
	# # Run the optimiser for each age class
	# for(i in 1:length(ages)){
		
		# # Go through the initiating list
		# res_temp_min<-data.frame(val=rep(NA, length(init_list)), P=NA, C=NA, F=NA) 
		# res_temp_max<-res_temp_min
		
		# # Now optimise for each starting value
		# for(o in 1:length(init_list)){
			
			# # Pull out the oth set of starting values
			# pars<-init_list[[o]]
			
			# # First get the combination estimated to mimise qx at age i and then repeat to get the max within intial values o
			# optim_results<-optim(par=pars, fn=convert_brass, sex=sexes[j], age=ages[i], GDP=med_GDP, year=year_plot, data=dataset_plot, standard=standard, GAM=GAM)
			# res_temp_min[o,]<-c(optim_results$value, optim_results$par)
			
			# optim_results<-optim(par=pars, fn=convert_brass, sex=sexes[j], age=ages[i], GDP=med_GDP, year=year_plot, data=dataset_plot, standard=standard, GAM=GAM, aim="maximise")
			# res_temp_max[o,]<-c(optim_results$value, optim_results$par)
		# }
			
		# # Save the results in to the list
		# max_expect[[i]]<-res_temp_max
		# min_expect[[i]]<-res_temp_min	
		
	# }
	
	# # Save the list to be returned
	# names(max_expect)<-ages
	# names(min_expect)<-ages
	# predictions<-list(min_expect, max_expect)
	# return(predictions)

# }

# # Save and reload (takes a while to get these data)
# save(predictions_list, file="Predictions_imp.Rdata")
load("Predictions_imp.Rdata")

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

# Go through each and plot the composition and total energy
min_list<-list()
max_list<-list()
ratio_list<-list()
 
for(i in 1:length(sexes)){
	
	# Pull out the ith data set
	data_i1<-predictions_list[[i]][[2]]
	data_i1$type<-"Max"
	data_i2<-predictions_list[[i]][[1]]
	data_i2$type<-"Min"
		
	# Add a categorical predictor for each row to color for age class
	data_i1$age<-standard$Age[c(1:dim(data_i1)[1])]
	data_i2$age<-standard$Age[c(1:dim(data_i2)[1])]
	
	# Plot the max 	
	max_list[[i]]<-ggplot(data_i1, aes(x=ages, y=Fat.kcal, color="Fat")) + 
			geom_point(size=3) + 
			geom_line() +
			geom_point(aes(x=ages, y=Carbo.kcal, color="Carbohydrate"), size=3) +
			geom_line(aes(x=ages, y=Carbo.kcal, color="Carbohydrate")) +
			geom_point(aes(x=ages, y=Protein.kcal, color="Protein"), size=3) +
			geom_line(aes(x=ages, y=Protein.kcal, color="Protein")) +
			ylim(c(0, 2500)) + 
			labs(x="Age", y="kcal/cap/day", subtitle=paste0("High Mortality, ", sexes[i])) + 
			theme_bw() +
			theme(axis.line = element_line(colour = "black")) + 
			theme(axis.text=element_text(size=15), axis.title=element_text(size=15)) + 
			theme(legend.position="none")
				
	# Plot the min 	
	min_list[[i]]<-ggplot(data_i2, aes(x=ages, y=Fat.kcal, color="Fat")) + 
			geom_point(size=3) + 
			geom_line() +
			geom_point(aes(x=ages, y=Carbo.kcal, color="Carbohydrate"), size=3) +
			geom_line(aes(x=ages, y=Carbo.kcal, color="Carbohydrate")) +
			geom_point(aes(x=ages, y=Protein.kcal, color="Protein"), size=3) +
			geom_line(aes(x=ages, y=Protein.kcal, color="Protein")) +
			ylim(c(0, 2500)) + 
			labs(x="Age", y="kcal/cap/day", subtitle=paste0("Low Mortality, ", sexes[i]), color=element_blank()) + 
			theme_bw() +
			theme(axis.line = element_line(colour = "black")) + 
			theme(axis.text=element_text(size=15), axis.title=element_text(size=15))
	
	if(i == 1){
		min_list[[i]]<-min_list[[i]] + 
			theme(legend.position=c(0.23, 0.85), legend.text=element_text(size=15))
	}else{
		min_list[[i]]<-min_list[[i]] + 
			theme(legend.position="none")
	}

						
	# Estimate the log ratio of probability of ASM
	q_max<-data_i1[,"qx"]
	q_min<-data_i2[,"qx"]
	
	# Put together as a dataset to plot
	delta_data<-data.frame(age=standard$Age[c(1:dim(data_i1)[1])], max=q_max, min=q_min)
	delta_data$Ratio<-delta_data$min / delta_data$max
	
	ratio<-ggplot(delta_data, aes(x=age, y=max)) +
						geom_point(aes(color="max"), size=3) + 
						geom_line(aes(color="max")) +
						geom_point(aes(x=age, y=min, color="min"), size=3)	+
						geom_line(aes(x=age, y=min, color="min")) + 
						geom_point(aes(x=age, y=Ratio, color="ratio"), size=3)	+
						geom_line(aes(x=age, y=Ratio, color="ratio")) + 
						labs(y ="Probability / Ratio", x=expression(Age~~italic(x)), subtitle=sexes[i], color=element_blank()) +
						theme_bw() +
						theme(axis.text=element_text(size=15), axis.title=element_text(size=15)) +
						scale_colour_manual(values=c("max" = "red", "min" = "blue", "ratio" = "black"), labels=c(expression(q[x]^high), expression(q[x]^low), expression(q[x]^low / q[x]^high))) 
						
	# Add a legend on pass 1
	if(i == 1){
		ratio<-ratio + theme(legend.position=c(0.23, 0.333), legend.text=element_text(size=15))
	}else{
		ratio<-ratio + theme(legend.position="none")
	}
	
	ratio_list[[i]]<-ratio	
				
}

# Now lets arrange all those plots
pdf("Figure_S3.pdf", height=10, width=15)

grid.arrange(min_list[[1]]+labs(title="A"),
				max_list[[1]]+labs(title="B"), 
				ratio_list[[1]]+labs(title="C"),
				min_list[[2]]+labs(title="D"),
				max_list[[2]]+labs(title="E"), 
				ratio_list[[2]]+labs(title="F"),
									layout_matrix=rbind(c(1,2,3),
														c(4,5,6)))

dev.off()


