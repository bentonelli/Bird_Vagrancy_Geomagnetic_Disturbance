#Caterpillar plots for interaction, and base models

library(MCMCvis)

### Base models ####
model_in <- readRDS("Fit_models_and_data/Models/Spring_Ap.rds")

mb_in <- MCMCsummary(model_in,params = "mu_beta",probs = c(.25,.75))
sig_mb <- which(mb_in$`25%` > 0 | mb_in$`75%` < 0)
col_in <- rep("grey60",nrow(mb_in))
col_in[sig_mb] <- "black"
MCMCplot(model_in,params = "mu_beta",rank = TRUE,guide_axis = TRUE, col = col_in,ylab=NULL,labels=NULL,
         sz_thick = 2,sz_thin = 1,sz_med = 1.1,sz_labels = .5,xlab=NULL)

### For interaction plots ####
model_in <- readRDS("Fit_models_and_data/Models/Fall_interaction.rds")

#Get betas - AP
mb <- MCMCsummary(model_in,params=c("beta"),probs = c(.25,.75))
#Thetas - SS
mt <- MCMCsummary(model_in,params=c("theta"),probs = c(.25,.75))
#Omegas - AP x SS
mo <- MCMCsummary(model_in,params=c("omega"),probs = c(.25,.75))

#Get median estimates of mu_beta, mu_theta, mu_omega
median_top_level <- MCMCsummary(model_in,params = c("mu_beta","mu_theta","mu_omega"))

sig_mb <- which(mb$`25%` > 0 | mb$`75%` < 0)
col_in <- rep("grey60",nrow(mb))
col_in[sig_mb] <- "black"
#Plot each separately
MCMCplot(model_in,params = "beta",rank = TRUE,
         guide_axis = TRUE, col = col_in,ylab=NULL,labels=NULL,xlab=NULL,
         sz_thick = 2,sz_thin = 1,sz_med = 1.1,sz_labels = .5)

sig_mt <- which(mt$`25%` > 0 | mt$`75%` < 0)
col_in <- rep("grey60",nrow(mt))
col_in[sig_mt] <- "black"
MCMCplot(model_in,params = "theta",rank = TRUE,
         guide_axis = TRUE, col = col_in,ylab=NULL,labels=NULL,xlab=NULL,
         sz_thick = 2,sz_thin = 1,sz_med = 1.1,sz_labels = .5)

sig_mo <- which(mo$`25%` > 0 | mo$`75%` < 0)
col_in <- rep("grey60",nrow(mo))
col_in[sig_mo] <- "black"
MCMCplot(model_in,params = "omega",rank = TRUE,
         guide_axis = TRUE, col = col_in,ylab=NULL,labels=NULL,xlab=NULL,
         sz_thick = 2,sz_thin = 1,sz_med = 1.1,sz_labels = .5)

