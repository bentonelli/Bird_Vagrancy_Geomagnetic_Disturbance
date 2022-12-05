#PPC 2-9-21

library(MCMCvis)
library(rstan)
library(dplyr)
library(scales)
library(Rlab)

#Season of interest
seas <- "fall"

#Read in model and associated data
model_in <- readRDS("Fit_models_and_data/Models/Fall_Ap.rds")
spec_df <- readRDS("Fit_models_and_data/Data/Fall_Ap_df.rds")

#Number of samples of posterior
m <- 50

#Plot parameter fits
MCMCplot(model_in,params = c("omega","mu_nu","psi","eta"))

#Save species list separately
spec_list <- unique(data_in$spec_name)

#Get Ap /SS values
if (seas=="fall"){
  x.rep <- spec_df$SS_21
} else {
  x.rep <- spec_df$SS_norm
}

#Number of iterations
if (seas=="fall"){
  ni <- 3000 
}else{
  ni <- 2000
}


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

#Get age values
age.rep <- spec_df$AGE_CODE

mean_rec <- c()
median_rec <- c()
sd_rec <- c()
#For each of the samples,
for (i in 1:m){
  
  if(seas == "fall"){
    kappa.rep <- chains[[samp_chain[i]]][samp_iter[i],"kappa_age"]
    
    #Assign ages randomly to data
    age.rep[which(is.na(age.rep))] <- rbern(length(which(is.na(age.rep))),kappa.rep)
  }
  
  #Get parameter values, indexing for each specific year, species, for each record (vector)
  alpha.rep <- chains[[samp_chain[i]]][samp_iter[i],paste("alpha[",yr_sp.rep,"]",sep = "")]
  beta.rep <- chains[[samp_chain[i]]][samp_iter[i],paste("beta[",yr_sp.rep,"]",sep = "")]
  lambda.rep <- chains[[samp_chain[i]]][samp_iter[i],paste("lambda[",spec.rep,"]",sep = "")]
  nu.rep <- chains[[samp_chain[i]]][samp_iter[i],paste("nu[",spec.rep,"]",sep = "")]
  
  #get shape value for chain sample (scalar)
  shp.rep <- chains[[samp_chain[i]]][samp_iter[i],paste("shp[",spec.rep,"]",sep = "")]

  #Linear predictor (vector) 
  #Reference:
  #lp[j] = alpha[sp_yr[j]] + beta[sp_yr[j]] * x[j] + lambda[ii[j]] * age_obs[j] + nu[ii[j]] * age_obs[j] * x[j]; 
  
  lp.rep <- as.numeric(alpha.rep) + as.numeric(beta.rep) * as.numeric(x.rep) + 
    as.numeric(lambda.rep) * as.numeric(age.rep) + as.numeric(nu.rep) * as.numeric(age.rep) * as.numeric(x.rep)
  
  #Pull from gamma distribution with shape and rate
  y_rep <- rgamma(length(lp.rep),shape = shp.rep,rate = (shp.rep/exp(lp.rep)))
  
  #Plot y_rep density, save mean, median and sd
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



