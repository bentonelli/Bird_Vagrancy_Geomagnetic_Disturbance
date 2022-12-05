#Post-hoc phylogenetic analysis

#load packages
library(picante)
library(ape)
library(dplyr)
library(ggplot2)
library(ggtree)
library(readr)
library(MCMCvis)
library(phytools)

#Define whether model is an interaction model or not
interac <- FALSE

#Define which season is being tested
season_of_interest <- "fall"

#Read in appropriate model and associated data
model_in <- readRDS("Fit_models_and_data/Models/Fall_SS.rds")
data_in <- readRDS("Fit_models_and_data/Data/Fall_SS_df.rds")

#Read in appropriate parameter values based on model type
if (interac){
  model_summary <- MCMCsummary(model_in,params=c("beta","theta","omega","mu_theta","mu_beta","mu_omega","psi","eta"))
} else {
  model_summary <- MCMCsummary(model_in,params=c("mu_beta","omega","psi","eta"))
}

#read in species data
if (season_of_interest == "fall"){
  spec_traits_f <- read_csv("Data/Spec_data/Fall_spec_traits.csv", col_names = TRUE)
  colnames(spec_traits_f)[1] <- "spec_name"
  spec_traits_f <- filter(spec_traits_f, spec_name %in% unique(data_in$spec_name))
  #Species traits, in order, from model run
  spec_traits_tax <- spec_traits_f
} else {
  spec_traits_s <- read_csv("Data/Spec_data/Spring_spec_traits.csv", col_names = TRUE)
  colnames(spec_traits_s)[1] <- "spec_name"
  spec_traits_s <- filter(spec_traits_s, spec_name %in% unique(data_in$spec_name))
  #Species traits, in order, from model run
  spec_traits_tax <- spec_traits_s
}

if (interac){
  spec_traits_tax$mu_diff_t <- NA
  spec_traits_tax$mu_diff_b <- NA
  spec_traits_tax$mu_diff_o <- NA
  #For each species, get difference in mu.theta - theta
  mu_theta_in <- model_summary["mu_theta","50%"]
  mu_beta_in <- model_summary["mu_beta","50%"]
  mu_omega_in <- model_summary["mu_omega","50%"]
  theta_in <- MCMCsummary(model_in,params="theta")
  beta_in <- MCMCsummary(model_in,params="beta")
  omega_in <- MCMCsummary(model_in,params="omega")
  spec_traits_tax$theta <- theta_in$`50%`
  spec_traits_tax$beta <- beta_in$`50%`
  spec_traits_tax$omega <- omega_in$`50%`
  
} else {
  omega_in <- model_summary["omega","50%"]
  nu_in <-MCMCsummary(model_in,params="nu")
  mu_beta_in <- MCMCsummary(model_in,params="mu_beta")
  spec_traits_tax$mu_beta <- mu_beta_in$`50%`
  spec_traits_tax$nu <- nu_in$`50%`
}


if (interac){
  for (i in 1:nrow(spec_traits_tax)){
    #Get taxon
    theta_in <- model_summary[paste("theta","[",i, "]",sep=""),"50%"]
    spec_traits_tax$mu_diff_t[i] <- theta_in - mu_theta_in
    
    #Get taxon
    beta_in <- model_summary[paste("beta","[",i, "]",sep=""),"50%"]
    spec_traits_tax$mu_diff_b[i] <- beta_in - mu_beta_in
    
    omega_in <- model_summary[paste("omega","[",i, "]",sep=""),"50%"]
    spec_traits_tax$mu_diff_o[i] <- omega_in - mu_omega_in
    
  }
}else{
  for (i in 1:nrow(spec_traits_tax)){
    #Get taxon
    beta_in <- model_summary[paste("mu_beta","[",i, "]",sep=""),"50%"]
    spec_traits_tax$mu_diff_b[i] <- beta_in - omega_in
  }
}



#Filter to species available in the tree - 
# Note here that some species are excluded due to lack of tree data
tree <- ape::read.nexus("Data/Spec_data/Phylo/full_tree_2_22/output.nex")
spec_tree_list <- gsub("_"," ",tree$tree_8214$tip.label)
spec_traits_tax_f <- filter(spec_traits_tax,`Scientific Name` %in% spec_tree_list)
inc_ind <- which(spec_traits_tax$`Scientific Name` %in% spec_tree_list)
drop_tips <- spec_tree_list[!(spec_tree_list %in% spec_traits_tax_f$`Scientific Name`)]
drop_tips <- gsub(" ","_",drop_tips)
trees<-lapply(tree,drop.tip,tip=drop_tips)

#Look at a random tree
plt_tree <- trees[[3]]
p <- ggtree(plt_tree)

#Process species name data for analysis
bn <- spec_traits_tax_f$`Scientific Name`
bn2 <- gsub(' ', '_', bn)
bn3 <- data.frame(name = bn2, num = 1:length(bn2))
j_idx <- dplyr::left_join(data.frame(name = plt_tree$tip.label), bn3)

#Real variables
#For Test 3 - Do a MLR to take into account the effect of migration length and breeding latitude
if (interac){
  param_name <- "omega" #change below too
  lmfit1 <- (lm(spec_traits_tax$omega~
                  spec_traits_tax$`Breeding Lat`+spec_traits_tax$`Migration Length`))
  pred_sens <- lmfit1$coefficients[1] + lmfit1$coefficients[2] * spec_traits_tax$`Breeding Lat` +
    lmfit1$coefficients[3] * spec_traits_tax$`Migration Length`
  pred_sens <- pred_sens[inc_ind]
} else {
  param_name <- "nu" #change below too
  lmfit1 <- (lm(spec_traits_tax$mu_beta~
                  spec_traits_tax$`Breeding Lat`+spec_traits_tax$`Migration Length`))
  pred_sens <- lmfit1$coefficients[1] + lmfit1$coefficients[2] * spec_traits_tax$`Breeding Lat` +
    lmfit1$coefficients[3] * spec_traits_tax$`Migration Length`
  pred_sens <- pred_sens[inc_ind]
}

#For Test 2,3 - use bootstrapping to extract chain values randomly, 
#one for each of the trials done below.
# The idea here is to incorporate the uncertainty of the model median estimates

if (interac){
  in_chains <- MCMCchains(model_in,params=c(param_name))
  in_chains <- in_chains[,inc_ind]
  sens_trait_v <- c()
  for (j in 1:length(trees)){
    sens_trait_add <- in_chains[sample(1:nrow(in_chains),1),] - pred_sens
    sens_trait_v <- rbind(sens_trait_v,sens_trait_add)
  }
} else {
  beta_chains <- MCMCchains(model_in,params=c("nu"))
  beta_chains <- beta_chains[,inc_ind]
  sens_trait_v <- c()
  for (j in 1:length(trees)){
    sens_trait_add <- beta_chains[sample(1:nrow(beta_chains),1),] - pred_sens
    sens_trait_v <- rbind(sens_trait_v,sens_trait_add)
  }
}

#Plot bootstrapped measure to make sure everything looks good
sens_trait_test <- as.numeric(sens_trait_v[1,])

p <- ggtree(plt_tree)
d1 <- data.frame(id=plt_tree$tip.label, val=sens_trait_test[j_idx$num])
p2 <- facet_plot(p, panel="dot", data=d1, geom=geom_point, aes(x=val), color=alpha('darkorchid4',.5))
p2

#calculate Bloomberg's K (and optionally, Pagel's lambda) for each tree
out.df <- data.frame(K = rep(NA, length(trees), PIC.var.P = NA))
lambda_rec <- c()
for (i in 1:length(trees)){
  #i <- 1
  #Optional, bootstrapping approach
  sens_trait <- as.numeric(sens_trait_v[i,])
  names(sens_trait) <- bn3$name
  tree_n <- trees[[i]]
  phy_res <- picante::phylosignal(sens_trait, tree_n,reps=1000)
  lambda_pagel <- phytools::phylosig(tree_n,sens_trait,method = "lambda")
  lambda_rec <- c(lambda_rec,lambda_pagel$lambda)
  out.df$K[i] <- phy_res$K
  out.df$PIC.var.P[i] <- phy_res$PIC.variance.P
}

#summarize output
hist(out.df$PIC.var.P)
sum(out.df$PIC.var.P < 0.05)/length(trees)
mean(out.df$K)
sd(out.df$K)

hist(lambda_rec)
median(lambda_rec)

