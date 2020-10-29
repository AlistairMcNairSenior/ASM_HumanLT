
# FUnction to work out whether points in nutrient spaace fall within the hull given by the raw data

inhull <- function(testpts, calpts, hull=convhulln(calpts), tol=mean(mean(abs(calpts)))*sqrt(.Machine$double.eps)) { 

	require(sp)
	require(geometry)

	# https://tolstoy.newcastle.edu.au/R/e8/help/09/12/8784.html
	calpts <- as.matrix(calpts) 
	testpts <- as.matrix(testpts) 
	p <- dim(calpts)[2] 
	cx <- dim(testpts)[1] # rows in testpts
	nt <- dim(hull)[1] # number of simplexes in hull 
	nrmls <- matrix(NA, nt, p)
	
	degenflag <- matrix(TRUE, nt, 1) 
	for (i in 1:nt){ 
		nullsp<-t(Null(t(calpts[hull[i,-1],] - matrix(calpts[hull[i,1],],p-1,p, byrow=TRUE))))  
		if (dim(nullsp)[1] == 1){
			nrmls[i,]<-nullsp
			degenflag[i]<-FALSE
		}
	}

	if(length(degenflag[degenflag]) > 0) warning(length(degenflag[degenflag])," degenerate faces in convex hull")
	nrmls <- nrmls[!degenflag,] 
	nt <- dim(nrmls)[1] 
	
	center = apply(calpts, 2, mean) 
	a<-calpts[hull[!degenflag,1],] 
	nrmls<-nrmls/matrix(apply(nrmls, 1, function(x) sqrt(sum(x^2))), nt, p)
	
	dp <- sign(apply((matrix(center, nt, p, byrow=TRUE)-a) * nrmls, 1, sum))
	nrmls <- nrmls*matrix(dp, nt, p)
	
	aN <- diag(a %*% t(nrmls)) 
	val <- apply(testpts %*% t(nrmls) - matrix(aN, cx, nt, byrow=TRUE), 1,min) 
	
	val[abs(val) < tol] <- 0 
	as.integer(sign(val)) 
}


# Function to drop all missing data for a set of variables from a predictor
drop.missing<-function(data, check, verbose=T){
	
	if(verbose == T){
		print("n.rows before:")
		print(dim(data)[1])
	}
		
	for(i in 1:length(check)){
		missing<-which(is.na(data[,check[i]]) == T)
		if(length(missing) > 0){
			data<-data[-missing,]
		}
	}
	
	if(verbose == T){
		print("n.rows after:")
		print(dim(data)[1])
	}
	
	return(data)
}


# Write GAM

write.GAM<-function(GAM, csv.file){
	
	# create the csv file to write to, if you want one
	write.table(Sys.time(), file=csv.file, sep=",", row.names=F, col.names=F)
	write.table(" ", file=csv.file, sep=",", row.names=F, col.names=F, append=T)
	
	# Add the formula
	write.table(as.character(GAM$formula), file=csv.file, sep=",", row.names=F, col.names=F, append=T)
	write.table(" ", file=csv.file, sep=",", row.names=F, col.names=F, append=T)
	
	# Get the summary
	summary_GAM<-summary(GAM)
	
	# Write the linear terms
	p.table<-round(summary_GAM$p.table, 4)
	p.table<-as.data.frame(cbind(row.names(p.table), p.table))
	names(p.table)[1]<-"Coef."
	suppressWarnings(write.table(p.table, file=csv.file, sep=",", row.names=F, col.names=names(p.table), append=T))
	write.table(" ", file=csv.file, sep=",", row.names=F, col.names=F, append=T)
	
	# Write the smooth terms
	s.table<-round(summary_GAM$s.table, 4)
	s.table<-as.data.frame(cbind(row.names(s.table), s.table))
	names(s.table)[1]<-"Coef."
	suppressWarnings(write.table(s.table, file=csv.file, sep=",", row.names=F, col.names=names(s.table), append=T))	
	write.table(" ", file=csv.file, sep=",", row.names=F, col.names=F, append=T)
	
	# Write the n and deviance explained
	dev.expl<-paste0("n = ", summary_GAM$n, ": % Dev. Explained = ", round(summary_GAM$dev.expl * 100, 2), ": AIC = ", round(AIC(GAM)))
	write.table(dev.expl, file=csv.file, sep=",", row.names=F, col.names=F, append=T)		
	write.table(" ", file=csv.file, sep=",", row.names=F, col.names=F, append=T)


}


ggSurface<-function(GAM, data, XYZ, labels, predict_val, surf_min=NA, surf_max=NA, x.limits=NA, y.limits=NA, z.val=NA, exclude, subtitle="", traits=NA, lab_size=3, nlevels=5, contour_at=NA, skip=0){
	
	require(ggplot2)
	require(sp)
	require(geometry)
	require(mgcv)
	require(metR)
	
	# This specifies the color scheme for surface
	rgb.palette<-colorRampPalette(c("blue","cyan","yellow","red"), space="Lab", interpolate="linear")
	map<-rgb.palette(256)
	
	# What are the outcomes being modelled, if not specified
	if(is.na(traits)[1]){
		traits<-unlist(lapply(strsplit(as.character(summary(GAM)$formula), " ~ "), "[[", 1))
	}
	
	# List to hold the plots
	plots_list<-list()

	# List for the order of plots
	nutrient.order<-XYZ[c(1,2,3)]
	
	# List for the labels
	labels.order<-labels[c(1,2,3)]
				
	# Values to predict over, if they are unspecified
	if(is.na(x.limits)[1] == T){
		x.limits<-c(floor(min(data[,nutrient.order[1]])), ceiling(max(data[,nutrient.order[1]])))
	}
	if(is.na(y.limits)[1] == T){
		y.limits<-c(floor(min(data[,nutrient.order[2]])), ceiling(max(data[,nutrient.order[2]])))
	}		
			
	# If we do not specify values to slice at, use the 25, 50, and 75 %ile
	if(is.na(z.val) == T){
		z.val<-round(median(data[,nutrient.order[3]]))
	}
			
	# Fitted list to hold some results for later
	x.new<-seq(min(x.limits, na.rm=T), max(x.limits, na.rm=T), len=501)
	y.new<-seq(min(y.limits, na.rm=T), max(y.limits, na.rm=T), len=501)
	z.new<-z.val
	predictors<-as.data.frame(expand.grid(x.new, y.new, z.new))
	names(predictors)<-nutrient.order
	in.poly<-as.numeric(inhull(predictors[,c(1:3)], data[,names(predictors)]) != -1)
			
	# Add the predictors for the additional 'confounders'
	predictors<-cbind(predictors, predict_val)
	predictors<-predictors[-which(in.poly == 0),]
			
	# Do the predictions
	predictions<-predict(GAM, newdata=predictors, type="response", exclude=exclude, newdata.guaranteed=T)

	# Loop for the proteins
	for(k in 1:length(traits)){
		
		# Get the proedictions for the kth trait					
		predictions_k<-predictions[,k]
								
		# Find the min and max values across all predictions
		mn<-surf_min[k]
		mx<-surf_max[k]
		if(is.na(mn)==T){
			mn<-min(predictions_k, na.rm=T)
		}
		if(is.na(mx)==T){
			mx<-max(predictions_k, na.rm=T)
		}
		locs<-(range(predictions_k, na.rm=TRUE) - mn) / (mx-mn) * 256	
		
		plot_data<-predictors
		plot_data$fit<-predictions_k
		plot_data$x<-plot_data[,nutrient.order[1]]
		plot_data$y<-plot_data[,nutrient.order[2]]
		
		# Set the contour
		if(is.na(contour_at)[1] == T){
			contour_use<-signif((max(predictions_k, na.rm=T)-min(predictions_k, na.rm=T))/nlevels, 1)
		}else{
			contour_use<-contour_at	
		}
		
		# Make the plot
		plot<-ggplot(plot_data, aes(x=x, y=y)) +
				geom_raster(aes(fill=fit), show.legend=F, interpolate=F, na.rm=T) +
				scale_fill_gradientn(colors=map[locs[1]:locs[2]]) +
				geom_contour(data=plot_data, aes(x=x, y=y, z=fit), na.rm=T, color="black", binwidth=contour_use) +	
				geom_label_contour(data=plot_data, aes(x=x, y=y, z=fit), size=lab_size, binwidth=contour_use, skip=skip) +
				theme_bw() +
				labs(x = labels.order[1], y = labels.order[2], subtitle=subtitle) +
				theme(axis.text=element_text(size=15), axis.title=element_text(size=15)) +
				theme(title=element_text(size=15)) + 
				xlim(x.limits) +
				ylim(y.limits)
				
		# Save the plot		
		plots_list[[k]]<-plot
	}
	
	return(plots_list)

}

# Function to estimate ASM for a given age, sex, year, and PCF - designed to work with an optimiser so we can get the PCF associated with maximal/minimal ASM

convert_brass<-function(x, sex, age, GDP, year, data, standard, GAM, stat="qx", aim="minimise"){
		
	# Make sure the required packages are loaded on the cores
	require(arm)
	require(mgcv)
	require(MortalityLaws)
	
	# Create a copy for analysing
	predict_i<-data.frame(Protein.kcal=x[1], Carbo.kcal=x[2], Fat.kcal=x[3], Sex=sex, GDP_perCapita=GDP, Year=year)
	
	# check it is within the observed data
	in.poly<-as.numeric(inhull(predict_i[,c("Protein.kcal", "Carbo.kcal", "Fat.kcal")], data[,c("Protein.kcal", "Carbo.kcal", "Fat.kcal")]) != -1)
	
	# If it is observed do the converison, else return 1000 + ED from the mean supply
	if(in.poly == 1){
		
		# Get the predictions
		predictions<-as.data.frame(predict(GAM, newdata=predict_i, exclude=c("s(Country)", "s.1(Country)"), newdata.guaranteed=T))
		names(predictions)<-c("alpha", "beta")
		predictions<-cbind(predict_i, predictions)
		
			
		# Which sex are we making predictions for
		convert1<-c("Males", "Females")
		sex_used<-match(predictions$Sex[1], convert1)
					
		# Get the jth lifetable
		tag<-which(standard$Age == 5)
		ls5<-standard[tag, grep("lx", names(standard))[sex_used]]
		l5<-invlogit(predictions$alpha + predictions$beta * logit(ls5))
		tag<-which(standard$Age == 60)
		ls60<-standard[tag, grep("lx", names(standard))[sex_used]]
		l60<-invlogit(predictions$alpha + predictions$beta * logit(ls60))
		l.x.<-invlogit(predictions$alpha + predictions$beta * logit(standard[, grep("lx", names(standard))[sex_used]]) + standard[, grep("gamma", names(standard))[sex_used]] * (1 - (logit(l5)/logit(ls5))) + standard[, grep("theta", names(standard))[sex_used]] * (1 - (logit(l60)/logit(ls60))))
			
		# Get the full lifetable
		convert2<-c("male", "female")
		LT<-LifeTable(x=standard$Age, lx=l.x., sex=convert2[sex_used])$lt
		
		# What are we looking for
		target<-LT[which(LT$x == age), stat]
		
		# If we want to maxise the qx multiply by -1
		if(aim == "maximise"){
			target<-target * -1
		}
		
	}else{
		# If we are outside the boundary, return 1000 + ED from the centre - helps guide the optimiser back there
		target<-1000 + sqrt((mean(data$Protein.kcal) - x[1])^2 + (mean(data$Carbo.kcal) - x[2])^2 + (mean(data$Fat.kcal) - x[3])^2)
	}
	
	return(target)
			
}


# A version of above but for general use

convert_brass_general<-function(x, sex, age, GDP, year, standard, GAM, stat="qx"){
		
	# Make sure the required packages are loaded on the cores
	require(arm)
	require(mgcv)
	require(MortalityLaws)
	
	# Create a copy for analysing
	predict_i<-data.frame(Protein.kcal=x[1], Carbo.kcal=x[2], Fat.kcal=x[3], Sex=sex, GDP_perCapita=GDP, Year=year)
		
	# Get the predictions
	predictions<-as.data.frame(predict(GAM, newdata=predict_i, exclude=c("s(Country)", "s.1(Country)"), newdata.guaranteed=T))
	names(predictions)<-c("alpha", "beta")
	predictions<-cbind(predict_i, predictions)
	
		
	# Which sex are we making predictions for
	convert1<-c("Males", "Females")
	sex_used<-match(predictions$Sex[1], convert1)
				
	# Get the jth lifetable
	tag<-which(standard$Age == 5)
	ls5<-standard[tag, grep("lx", names(standard))[sex_used]]
	l5<-invlogit(predictions$alpha + predictions$beta * logit(ls5))
	tag<-which(standard$Age == 60)
	ls60<-standard[tag, grep("lx", names(standard))[sex_used]]
	l60<-invlogit(predictions$alpha + predictions$beta * logit(ls60))
	l.x.<-invlogit(predictions$alpha + predictions$beta * logit(standard[, grep("lx", names(standard))[sex_used]]) + standard[, grep("gamma", names(standard))[sex_used]] * (1 - (logit(l5)/logit(ls5))) + standard[, grep("theta", names(standard))[sex_used]] * (1 - (logit(l60)/logit(ls60))))
		
	# Get the full lifetable
	convert2<-c("male", "female")
	LT<-LifeTable(x=standard$Age, lx=l.x., sex=convert2[sex_used])$lt
	
	# What are we looking for
	target<-LT[which(LT$x == age), stat]
			
	return(target)
			
}

