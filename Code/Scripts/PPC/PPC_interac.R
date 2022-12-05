#PPC for interaction models

library(MCMCvis)
library(rstan)
library(dplyr)
library(scales)
library(Rlab)

#Season of interest
seas <- "spring"

#Read in model and associated data
model_in <- readRDS("Fit_models_and_data/Models/Spring_interaction.rds")
spec_df <- readRDS("Fit_models_and_data/Data/Spring_interaction_df.rds")

#Plot parameter fits
MCMCplot(model_in,params = c("mu_omega","mu_beta","mu_theta"))

#Save species list separately
spec_list <- unique(data_in$spec_name)

x1.rep <- spec_df$Ap_norm
x2.rep <- spec_df$SS_norm

#Number of samples of posterior
m <- 50

#Number of iterations
ni <- 2000

#Select chains to sample for each draw
chains <- As.mcmc.list(model_in)
samp_chain <- sample(1:4,m,replace = TRUE)

#Select which iteration to sample for each
samp_iter <- sample(1:ni,m,replace = TRUE)

#Set up plot
par(mfrow=c(1,1))

#Plot the vagrancy index values
par(las=1,mar=c(6,6,6,6))
par(pty="s")
break_algo <- seq(from=0, to=100000, by=10)
hg <-hist(spec_df$vagrancy_index,xlim=c(0,200),breaks=break_algo,prob=TRUE)
plot(hg, col = alpha("dodgerblue4",.2),xlim=c(0,200),ylab="",xlab="Vagrancy Index",main="") # Plot 1st histogram using a transparent color
title(ylab = "Frequency", line = 4)     
par(new=TRUE)
plot(density(spec_df$vagrancy_index,bw=2,from = 0,to=200),lwd=1.5,axes = FALSE,main="",xlab="",ylab="")

#Get the species, and years of each record
spec.rep <- spec_df$spec_name_fct
year.rep <- spec_df$YEAR_ind

#Get alpha, beta index
yr_sp.rep <- spec_df$sp_yr

mean_rec <- c()
median_rec <- c()
sd_rec <- c()
#For each of the samples,
for (i in 1:m){
  
  #Get parameter values, indexing for each specific year, species, for each record (vector)
  alpha.rep <- chains[[samp_chain[i]]][samp_iter[i],paste("alpha[",yr_sp.rep,"]",sep = "")]
  beta.rep <- chains[[samp_chain[i]]][samp_iter[i],paste("beta[",spec.rep,"]",sep = "")]
  theta.rep <- chains[[samp_chain[i]]][samp_iter[i],paste("theta[",spec.rep,"]",sep = "")]
  omega.rep <- chains[[samp_chain[i]]][samp_iter[i],paste("theta[",spec.rep,"]",sep = "")]
  shp.rep <- chains[[samp_chain[i]]][samp_iter[i],paste("shp[",spec.rep,"]",sep = "")]
  
  #Linear predictor (vector) 
  #Reference:
  #lp = alpha[sp_yr] + beta[ii] .* x_ap + theta[ii] .* x_ss + omega[ii] .* x_ap .* x_ss; 
  
  lp.rep <- as.numeric(alpha.rep) + as.numeric(beta.rep) * as.numeric(x1.rep) +  as.numeric(theta.rep) * as.numeric(x2.rep) + 
    as.numeric(omega.rep) * as.numeric(x1.rep) * as.numeric(x2.rep)
  
  
  #Pull from gamma distribution with shape and rate
  y_rep <- rgamma(length(lp.rep),shape = shp.rep,rate = (shp.rep/exp(lp.rep)))
  
  #Add line, calculate median, mean and sd
  lines(density(y_rep,bw=2,from = 0,to=200),lwd=1.5,main="",xlab="",ylab="",col=alpha("seagreen3",.2))
  mean_rec <- c(mean_rec,mean(y_rep))
  median_rec <- c(median_rec,median(y_rep))
  sd_rec <- c(sd_rec,sd(y_rep))
}

mean(spec_df$vagrancy_index)
median(spec_df$vagrancy_index)
sd(spec_df$vagrancy_index)

range(mean_rec)
range(median_rec)
range(sd_rec)



