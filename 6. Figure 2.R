
# Clean up 
rm(list=ls())

# Working directory
wd<-"/Users/alistairsenior/Dropbox/Human lifetables and Nutrition"
wd<-"/Users/alistairsenior/Dropbox (Sydney Uni)/Human lifetables and Nutrition"
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
library(doSNOW)

source("0.Header_Functions.R")

#################################################
#################### FIGURE 2 ###################
#################################################

# Full data
full_data<-read.csv("Brass_complete_cases.csv")

# Start with densities for alpha and beta by sex
dens_alpha<-ggplot(full_data, aes(alpha, color=Sex)) +
			geom_density(aes(fill=Sex), alpha=0.5, linetype=0) + 
			theme_bw() + 
			theme(legend.position=c(0.15, 0.85), legend.key.size=unit(0.4, "cm"), legend.title=element_blank(), legend.text=element_text(size=15), legend.background=element_blank()) +
			labs(x=expression(italic("\U03B1")), y="Density") +
			theme(axis.text=element_text(size=15), axis.title=element_text(size=15), title=element_text(size=15)) +
			xlim(-3, 3) +
			scale_fill_manual(values=c("#d8b365", "#5ab4ac"))

dens_beta<-ggplot(full_data, aes(beta, color=Sex)) +
			geom_density(aes(fill=Sex), alpha=0.5, linetype=0) + 
			theme_bw() + 
			theme(legend.position="none") +
			labs(x=expression(italic("\U03B2")), y="Density") +
			theme(axis.text=element_text(size=15), axis.title=element_text(size=15), title=element_text(size=15)) +
			xlim(0, 2) +
			scale_fill_manual(values=c("#d8b365", "#5ab4ac"))

# Lets check sex-specific correlations between alpha and e0
tag_f<-which(full_data$Sex == "Females")
tag_m<-which(full_data$Sex == "Males")

r_m<-round(cor(full_data$alpha[tag_m], full_data$e0[tag_m])^2, 2)
r_f<-round(cor(full_data$alpha[tag_f], full_data$e0[tag_f])^2, 2)
val<-c(substitute(expression(italic(r^2) == v), list(v=r_m)), substitute(expression(italic(r^2) == v), list(v=r_f)))
ys<-c(43, 35)
cols<-c("#5ab4ac", "#d8b365")

plot_e0_alpha<-ggplot(data=full_data, aes(x=alpha, y=e0, color=Sex)) +
				geom_point(alpha=0.1, size=0.5) +
				theme_bw() +
				theme(axis.text=element_text(size=15), axis.title=element_text(size=15)) +
				theme(title=element_text(size=15)) + 
				labs(x=expression(italic("\U03B1")), y=expression(italic(e[0]))) +
				theme(legend.position="none") +
				scale_color_manual(values=c("#d8b365", "#5ab4ac"))

for(i in 1:2){
	plot_e0_alpha<-plot_e0_alpha + annotate("text", label=lapply(val, eval)[[i]], x=1.5, y=ys[i], size=8, color=cols[i])
}						

# Lets check sex-specific correlations between beta and e0

r_m<-round(cor(full_data$beta[tag_m], full_data$e0[tag_m])^2, 2)
r_f<-round(cor(full_data$beta[tag_f], full_data$e0[tag_f])^2, 2)
val<-c(substitute(expression(italic(r^2) == v), list(v=r_m)), substitute(expression(italic(r^2) == v), list(v=r_f)))
ys<-c(43, 35)
cols<-c("#5ab4ac", "#d8b365")

plot_e0_beta<-ggplot(data=full_data, aes(x=beta, y=e0, color=Sex)) +
				geom_point(alpha=0.1, size=0.5) +
				theme_bw() +
				theme(axis.text=element_text(size=15), axis.title=element_text(size=15)) +
				theme(title=element_text(size=15)) + 
				labs(x=expression(italic("\U03B2")), y=expression(italic(e[0]))) +
				theme(legend.position="none") +
				scale_color_manual(values=c("#d8b365", "#5ab4ac"))

for(i in 1:2){
	plot_e0_beta<-plot_e0_beta + annotate("text", label=lapply(val, eval)[[i]], x=1.75, y=ys[i], size=8, color=cols[i])
}						

# Lets check sex-specific correlations between alpha and l5

r_m<-round(cor(full_data$alpha[tag_m], full_data$l5[tag_m])^2, 2)
r_f<-round(cor(full_data$alpha[tag_f], full_data$l5[tag_f])^2, 2)
val<-c(substitute(expression(italic(r^2) == v), list(v=r_m)), substitute(expression(italic(r^2) == v), list(v=r_f)))
ys<-c(0.765, 0.725)
cols<-c("#5ab4ac", "#d8b365")

plot_l5_alpha<-ggplot(data=full_data, aes(x=alpha, y=l5, color=Sex)) +
				geom_point(alpha=0.1, size=0.5) +
				theme_bw() +
				theme(axis.text=element_text(size=15), axis.title=element_text(size=15)) +
				theme(title=element_text(size=15)) + 
				labs(x=expression(italic("\U03B1")), y=expression(italic(l[5]))) +
				ylim(0.7, 1) +
				theme(legend.position="none") +
				scale_color_manual(values=c("#d8b365", "#5ab4ac"))

for(i in 1:2){
	plot_l5_alpha<-plot_l5_alpha + annotate("text", label=lapply(val, eval)[[i]], x=1.5, y=ys[i], size=8, color=cols[i])
}						

# Lets check sex-specific correlations between alpha and l60

r_m<-round(cor(full_data$alpha[tag_m], full_data$l60[tag_m])^2, 2)
r_f<-round(cor(full_data$alpha[tag_f], full_data$l60[tag_f])^2, 2)
val<-c(substitute(expression(italic(r^2) == v), list(v=r_m)), substitute(expression(italic(r^2) == v), list(v=r_f)))
ys<-c(0.2, 0.075)

plot_l60_alpha<-ggplot(data=full_data, aes(x=alpha, y=l60, color=Sex)) +
				geom_point(alpha=0.1, size=0.5) +
				theme_bw() +
				theme(axis.text=element_text(size=15), axis.title=element_text(size=15)) +
				theme(title=element_text(size=15)) + 
				labs(x=expression(italic("\U03B1")), y=expression(italic(l[60]))) +
				ylim(0, 1) +
				theme(legend.position="none") +
				scale_color_manual(values=c("#d8b365", "#5ab4ac"))

for(i in 1:2){
	plot_l60_alpha<-plot_l60_alpha + annotate("text", label=lapply(val, eval)[[i]], x=1.5, y=ys[i], size=8, color=cols[i])
}						

# Lets check sex-specific correlations between beta and l5
			
r_m<-round(cor(full_data$beta[tag_m], full_data$l5[tag_m])^2, 2)
r_f<-round(cor(full_data$beta[tag_f], full_data$l5[tag_f])^2, 2)
val<-c(substitute(expression(italic(r^2) == v), list(v=r_m)), substitute(expression(italic(r^2) == v), list(v=r_f)))
ys<-c(0.765, 0.725)

plot_l5_beta<-ggplot(data=full_data, aes(x=beta, y=l5, color=Sex)) +
				geom_point(alpha=0.1, size=0.5) +
				theme_bw() +
				theme(axis.text=element_text(size=15), axis.title=element_text(size=15)) +
				theme(title=element_text(size=15)) + 
				labs(x=expression(italic("\U03B2")), y=expression(italic(l[5]))) +
				ylim(0.7, 1) +
				theme(legend.position="none") +
				scale_color_manual(values=c("#d8b365", "#5ab4ac"))

for(i in 1:2){
	plot_l5_beta<-plot_l5_beta + annotate("text", label=lapply(val, eval)[[i]], x=1.75, y=ys[i], size=8, color=cols[i])
}						

# Lets check sex-specific correlations between beta and l60

r_m<-round(cor(full_data$beta[tag_m], full_data$l60[tag_m])^2, 2)
r_f<-round(cor(full_data$beta[tag_f], full_data$l60[tag_f])^2, 2)
val<-c(substitute(expression(italic(r^2) == v), list(v=r_m)), substitute(expression(italic(r^2) == v), list(v=r_f)))
ys<-c(0.2, 0.075)

plot_l60_beta<-ggplot(data=full_data, aes(x=beta, y=l60, color=Sex)) +
				geom_point(alpha=0.1, size=0.5) +
				theme_bw() +
				theme(axis.text=element_text(size=15), axis.title=element_text(size=15)) +
				theme(title=element_text(size=15)) + 
				labs(x=expression(italic("\U03B2")), y=expression(italic(l[60]))) +
				ylim(0, 1) +
				theme(legend.position="none") +
				scale_color_manual(values=c("#d8b365", "#5ab4ac"))

for(i in 1:2){
	plot_l60_beta<-plot_l60_beta + annotate("text", label=lapply(val, eval)[[i]], x=1.75, y=ys[i], size=8, color=cols[i])
}	

# Now lets arrange all those plots
cairo_pdf("Figure_2.pdf", height=15, width=10)

grid.arrange(dens_alpha+labs(title="A"), dens_beta+labs(title="B"), plot_e0_alpha+labs(title="C"), plot_e0_beta+labs(title="D"), plot_l5_alpha+labs(title="E"), plot_l60_alpha+labs(title="F"), plot_l5_beta+labs(title="G"), plot_l60_beta+labs(title="H"), 
									layout_matrix=rbind(c(1,2),
														c(3,4),
			  					  					    c(5,6),
			  					  					    c(7,8)))

dev.off()

