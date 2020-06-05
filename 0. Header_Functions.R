
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


# A function to plot MR gam surfaces

plot_MR.GAM<-function(GAM, data, XYZ_list=NA, predict_val=NA, exclude=NULL, csv.file=NA, pdf.file=NA, slice_at=NA, fit.resolution=101, no.cols=256, nlev=8, include_se=F, markers=NA, scale_surface=NA, cex.axis=2, labels_list=XYZ_list, cex.lab=2, direction=1){
	
	require(sp)
	require(geometry)
	require(mgcv)
	
	# If we want to plot the surfaces
	if(is.na(pdf.file) == F){
		# Set the layout
		# Open the pdf file for plotting
		if(include_se == T){
			# DO you want surfaces for SE
			if(direction == 1){
				# Set the layout for direction 1 - reading left to right: lower, middle, upper slice
				pdf(pdf.file, height=6 * 5, width=3 * 5)
				par(mfrow=c(6, 3), mar=c(6,6,5,1))
			}else{
				# Set the layout for direction 2 - reading top to bottom: lower, middle, upper slice
				pdf(pdf.file, height=3 * 5, width=6 * 5)
				par(mar=c(6,6,5,1))
				layout(as.matrix(array(seq(1, 6*3, 1), c(3, 6))))
			}	
		}else{
			pdf(pdf.file, height=3 * 5, width=3 * 5)
			if(direction == 1){
				par(mfrow=c(3, 3), mar=c(6,6,5,1))
			}else{
				par(mar=c(6,6,5,1))
				layout(as.matrix(array(seq(1, 3*3, 1), c(3, 3))))
			}		
		}
	}
	
	# Order the markers
	markers_list<-list()
	if(is.list(markers) == T){
		markers_list[[1]]<-list(markers[[1]], markers[[2]])
		markers_list[[2]]<-list(markers[[1]], markers[[3]])
		markers_list[[3]]<-list(markers[[2]], markers[[3]])
	}
	
	# This specifies the color scheme for surface
	rgb.palette<-colorRampPalette(c("blue","cyan","yellow","red"), space="Lab", interpolate="linear")
	map<-rgb.palette(no.cols)
	
	# create the csv file to write to, if you want one
	if(is.na(csv.file) == F){
		write.table(Sys.time(), file=csv.file, sep=",", row.names=F, col.names=F)
		write.table(" ", file=csv.file, sep=",", row.names=F, col.names=F, append=T)
		
		# Write the linear terms
		p.table<-round(summary(GAM)$p.table, 4)
		p.table<-as.data.frame(cbind(row.names(p.table), p.table))
		names(p.table)[1]<-"Coef."
		suppressWarnings(write.table(p.table, file=csv.file, sep=",", row.names=F, col.names=names(p.table), append=T))
		write.table(" ", file=csv.file, sep=",", row.names=F, col.names=F, append=T)
		
		# Write the smooth terms
		s.table<-round(summary(GAM)$s.table, 4)
		s.table<-as.data.frame(cbind(row.names(s.table), s.table))
		names(s.table)[1]<-"Coef."
		suppressWarnings(write.table(s.table, file=csv.file, sep=",", row.names=F, col.names=names(s.table), append=T))	
		write.table(" ", file=csv.file, sep=",", row.names=F, col.names=F, append=T)
		
		# Write the n and deviance explained
		dev.expl<-paste0("n = ", summary(GAM)$n, ": % Dev. Explained = ", round(summary(GAM)$dev.expl * 100, 2), ": AIC = ", round(AIC(GAM)))
		write.table(dev.expl, file=csv.file, sep=",", row.names=F, col.names=F, append=T)		
		write.table(" ", file=csv.file, sep=",", row.names=F, col.names=F, append=T)
	}

	# What are the outcomes being modelled
	traits<-unlist(lapply(strsplit(as.character(summary(GAM)$formula), " ~ "), "[[", 1))
		
	# Formatting the progress bar
 	pb <- txtProgressBar(min = 0, max = length(traits), style = 3)
	progress<-0
	
	# Loop for the proteins
	for(k in 1:length(traits)){
		
		# If we are plotting the surface make the surfaces
		if(is.na(pdf.file) == F){
			
			# List for the order of plots
			list.order<-list()
			list.order[[1]]<-XYZ_list[[k]][c(1,2,3)]
			list.order[[2]]<-XYZ_list[[k]][c(1,3,2)]
			list.order[[3]]<-XYZ_list[[k]][c(2,3,1)]
			
			# List for the labels
			labels.order<-list()
			labels.order[[1]]<-labels_list[[k]][c(1,2,3)]
			labels.order[[2]]<-labels_list[[k]][c(1,3,2)]
			labels.order[[3]]<-labels_list[[k]][c(2,3,1)]
			
			# Lists to hold the predictions for the nth combinations
			predictors.list<-list()
			predictions.list<-list()
			xyz.list<-list() 
			cv.list<-list()
			
			# Go through the nutrient orders and get the predicted values, note we will do the predictions for all three combinations, then find the minimla and maximal values, then go back though the nutrients orders and plot out
			for(n in 1:3){
			
				# Order to plot the nutrients
				nutrient.order<-list.order[[n]]
				
				# Values to predict over
				x.limits<-c(floor(min(data[,nutrient.order[1]])), ceiling(max(data[,nutrient.order[1]])))
				y.limits<-c(floor(min(data[,nutrient.order[2]])), ceiling(max(data[,nutrient.order[2]])))
				
				# If we do not specify values to slice at, use the 25, 50, and 75 %ile
				if(is.list(slice_at) == F){
					z.vals<-round(quantile(data[,nutrient.order[3]])[c(2:4)])
				}else{
					z.vals<-slice_at[[4-n]]	
				}
				
				# Fitted list to hold some results for later
				x.new<-seq(min(x.limits, na.rm=T), max(x.limits, na.rm=T), len=fit.resolution)
				y.new<-seq(min(y.limits, na.rm=T), max(y.limits, na.rm=T), len=fit.resolution)
				z.new<-z.vals
				predictors<-as.data.frame(expand.grid(x.new, y.new, z.new))
				names(predictors)<-nutrient.order
				in.poly<-as.numeric(inhull(predictors[,c(1:3)], data[,names(predictors)]) != -1)
				
				# Add the predictors for the additional 'confounders'
				predictors<-cbind(predictors, predict_val)
				
				# Do the predictions
				predictions<-predict(GAM, newdata=predictors, type="response", exclude=exclude, se.fit=T, newdata.guaranteed=T)
				
				# Edit out based on the marker list
				predictions$fit[which(in.poly == 0),]<-NA
				predictions$se.fit[which(in.poly == 0),]<-NA
				
				# Save the nth set of predictions
				predictions.list[[n]]<-predictions$fit[,k]
				predictors.list[[n]]<-predictors
				xyz.list[[n]]<-list(x.new, y.new, z.new)
				cv.list[[n]]<-predictions$se.fit[,k]
			}
			
			# Find the min and max values across all predictions
			mn<-min(unlist(predictions.list), na.rm=T)
			mx<-max(unlist(predictions.list), na.rm=T)
			
			# If no color scale is specified for the outcomes scale by the predicted values
			if(sum(is.na(scale_surface)) < length(scale_surface)){
				# Find the absolute max values across all predictions, and the scale
				upp_abs<-max(abs(c(scale_surface, mn, mx)))	
				mn<-(-upp_abs)
				mx<-upp_abs
			}
			
			# now do the coefficient of the error
			mn.cv<-min(unlist(cv.list), na.rm=T)
			mx.cv<-max(unlist(cv.list), na.rm=T)

			# Now go back though the predictions and plot
			for(n in 1:3){
				
				# Order to plot the nutrients
				nutrient.order<-list.order[[n]]
				labs<-labels.order[[n]]
				
				# Pull out the nth set of predictors and predictions
				predictors<-predictors.list[[n]]
				predictions<-predictions.list[[n]]
				x.new<-xyz.list[[n]][[1]]
				y.new<-xyz.list[[n]][[2]]
				z.new<-xyz.list[[n]][[3]]
					
				# Do the 3 quantiles for the predictions
				for(i in 1:length(z.new)){
								
					# Subset for the ith quantile
					ith_Quantile<-predictions[which(predictors[, nutrient.order[3]] == z.new[i])]
											
					surf<-matrix(ith_Quantile, nrow=fit.resolution)
								
					locs<-round((range(surf, na.rm=TRUE) - mn) / (mx-mn) * no.cols)
					image(x.new, y.new, surf, col=map[locs[1]:locs[2]], xlab="", ylab="", axes=FALSE)
					mtext(paste0(labs[3], " = ", z.new[i]), line=1, cex=cex.lab)
					mtext(labs[1], side=1, line=4, cex=cex.lab)
					mtext(labs[2], side=2, line=4, cex=cex.lab)
					if(i == 3 & n == 1){
						mtext(traits[k], line=2, font=1, cex=1, at = max(x.new)*0.9)
					}
					axis(1, cex.axis=cex.axis)
					axis(2, cex.axis=cex.axis)
					contour(x.new, y.new, surf, add=TRUE, levels=pretty(range(mn, mx), nlev), labcex=1, lwd=3)			
					# Add any markers
					if(is.list(markers) == T){
						abline(v=markers_list[[n]][[1]], col="grey")
						abline(h=markers_list[[n]][[2]], col="grey")
					}
					
				}
				
				# Now the 3 quantiles for the errors if we want those
				if(include_se == T){
					
					# Pull out the nth set of errors
					predictions<-cv.list[[n]]
					
					# GO through each quantile
					for(i in 1:length(z.new)){
									
						# Subset for the ith quantile
						ith_Quantile<-predictions[which(predictors[, nutrient.order[3]] == z.new[i])]
												
						surf<-matrix(ith_Quantile, nrow=fit.resolution)
									
						locs<-round((range(surf, na.rm=TRUE) - mn.cv) / (mx.cv-mn.cv) * no.cols)
						image(x.new, y.new, surf, col=map[locs[1]:locs[2]], xlab="", ylab="", axes=FALSE)
						mtext(paste0("se"), line=1, cex=cex.lab)
						mtext(labs[1], side=1, line=4, cex=cex.lab)
						mtext(labs[2], side=2, line=4, cex=cex.lab)
						axis(1, cex.axis=cex.axis)
						axis(2, cex.axis=cex.axis)
						contour(x.new, y.new, surf, add=TRUE, levels=pretty(range(mn.cv, mx.cv), nlev), labcex=1, lwd=3)
									
					}
				}
			}
		}		
		
		# Update the progress
		progress<-progress + 1
		setTxtProgressBar(pb, progress)
	}
	
	if(is.na(pdf.file) == F){
		# Close the file 
		dev.off()
	}
	
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

# Function to test residuals for a list of models via GAM
test.resid<-function(models, pdf.file, traits){
	
	# Open the file for plotting
	pdf(pdf.file)
	par(mfrow=c(1,1))
	
	# For each model in the list		
	for(i in 1:length(models)){
		
		# Get the predictions and pearson residuals	
		x<-predict(models[[i]])
		y<-resid(models[[i]], type="pearson")
		
		# Plot them
		plot(x, y, xlab="Predicted Values", ylab="Pearson Residuals", main=traits[i])
		abline(h=0, col="gray")
		
		# Fit a GAM to check
		if(sd(x) > 0){
			test<-gam(y ~ s(x))
			new.x<-seq(min(x), max(x), len=1000)
			new.y<-predict(test, newdata=data.frame(x=new.x), se=T)
			lines(new.x, new.y$fit, col=2, lwd=2)
			lines(new.x, new.y$fit + new.y$se.fit*1.96, col=2, lwd=2, lty=2)
			lines(new.x, new.y$fit - new.y$se.fit*1.96, col=2, lwd=2, lty=2)
			mtext(paste0("p = ", as.character(round(summary(test)$s.table), 3)[4]))	
		}	
	}
	
	# Close the file
	dev.off()

}

# Function that compares k different models for comparison - the most common would be a null model, and a model including some variable
# models are passed in a list of increasing complexity. E.g. in the most simple example, the first element would be the null model, and the second the more complex model
# Each element itself can be n.traits long for different outcomes/traits being modelled
# A list of n.traits with the LRTs for the pairwise comparisons is returned

model.select<-function(models, traits, alpha=0.05){
	
	# Get the pairwise comparisons amongst traits
	model_combn<-combn(length(models), 2)
	test_list<-list()
	
	# GO through each trait in turn
	for(i in 1:length(traits)){
		
		# Create a copy of combinations list
		model_combn_i<-model_combn
		model_combn_i<-rbind(model_combn_i, NA, NA)
		row.names(model_combn_i)<-c("M1", "M2", "p", "dev")
				
		# Now get all the model comparisons for this trait, and do a LRT and save the p-value
		for(j in 1:dim(model_combn_i)[2]){
			
			# Pull out those models
			model1ij<-models[[model_combn_i[1,j]]][i][[1]]
			model2ij<-models[[model_combn_i[2,j]]][i][[1]]
			
			# DO the LRT on the models
			test<-anova(model1ij, model2ij, test="Chisq")	
			
			# Add in the pvalues for the LRT
			if(is.na(test$Pr[2]) == T){
				model_combn_i[3,j]<-1
			}else{
				model_combn_i[3,j]<-test$Pr[2]
			}
			model_combn_i[4,j]<-test$Dev[2]
		}
		
		# Save those results for the ith trait
		test_list[[i]]<-model_combn_i
		names(test_list)[i]<-traits[i]
	}

	# Return the list of model comparisons
	return(test_list)
}



# Function to return and index calculated from life table parameters - it is a subset of the above code 

return_index<-function(q, x, type, l=NA){
	
	require(Matrix)
	require(signal)
	
	# Convert l to q if necessary	
	if(is.na(q) == T){
		p<-c(l[-1], 0) / l
		q<-(1 - p)
	}
	
	#------------------------------------------
	# Preliminaries
	#------------------------------------------
	# defining survival probabilities
	  p <- 1-q
	# number of transient states + 1, (in construction of U first row is zero)
	  s <- length(p)+1
	# number of transient states only
	  s2 <- length(p)
	# other things we will need later
	  I <- diag(rep(1,s))  # identity matrix
	  I <- as(I,"sparseMatrix")
	  e <- rep(1,s) # vector of ones for summations
	  e1 <- c(1,rep(0,s-1))
	  age <- x
	# C matrix is for calculating cumulative sums
	  C <- Matrix(0,nrow=s,ncol=s)
	    for (i in 1:s){
	      C[,i] <- c(rep(0,i),rep(1,s-i))
	      }
	  C <- as(C,"sparseMatrix")
	#--------------------------------------------------
	# Markov chain formulation of longevity
	#--------------------------------------------------
	# U matrix describes transient states in the Markov chain
	  U <- subdiag(p,-1)
	  U <- as(U,"sparseMatrix")
	# N matrix, where (i,j) entry is the mean time spent in each age class i,
	# conditional upon starting in age class j
	  N <- solve(I-U)
	  N <- as(N,"sparseMatrix")
	# M matrix has probability of death at each age on diagonal
	  M <- diag(c(q,1))
	  M <- as(M,"sparseMatrix")
	# The distribution of ages at death
	  B <- M %*% N  # the complete distribution of ages at death
	  f <- B %*% e1 # the distribution of ages at death from birth (1st age class)
	  B <- as(B,"sparseMatrix")
	# survivorship (alternatively ell <- e-C%*%f )
	  ell <- N %*% e1
	  ell <- as(ell,"sparseMatrix")
	# remaining life expectancy at each age
	  mean_eta <- colSums(N)-0.5
	# life expectancy at birth (or first age class)
	  eta <- mean_eta[1]
	# NB: in Markov chain formulations, the life expectancy at birth is always
	# 0.5 years higher than that found by conventional life table methods,
	# which is why we subtract 0.5 years
	#------------------------------------
	# Indices of lifespan variation
	#------------------------------------
	# variance in lifespan
	  V <- colSums(N) %*% (2*N-I) - mean_eta*mean_eta
	# standard deviation in lifespan
	  S <- sqrt(V)
	  
	  if(type=="S"){
	  	out<-S[1]
	  }
	  
	  if(type == "e0"){
	  	out<-eta
	  }
	  
	  return(out)
	  
}

#---------------------------------
#---- special functions
#---------------------------------
#--- function for making subdiagonals (from Bill Venables)
  subdiag <- function (v, k) {
    n <- length(v) + abs(k)
    x <- matrix(0, n, n)
    if (k == 0)
        diag(x) <- v
    else if (k < 0)
      { ## sub-diagonal
        j <- 1:(n+k)
        i <- (1 - k):n
        x[cbind(i, j)] <- v
} 
x } 


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

