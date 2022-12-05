#Figure 4

#Plot individual predicted species vagrancy under high and low predictions

library(readr)
library(ggplot2)
library(MCMCvis)
library(dplyr)

#Read in model and data
model_in <- readRDS("Fit_models_and_data/Models/Fall_interaction.rds")
geo_data <- readRDS("Fit_models_and_data/Data/Fall_interaction_df.rds")
spec_list <- unique(geo_data$spec_name)

model_chains <- data.frame(MCMCchains(model_in, "mu_beta"),
                           MCMCchains(model_in, "mu_theta"),
                           MCMCchains(model_in, "mu_omega"),
                           MCMCchains(model_in, "mu_shp"),
                           MCMCchains(model_in,"mu_alpha"),
                           MCMCchains(model_in,"beta"),
                           MCMCchains(model_in,"theta"),
                           MCMCchains(model_in,"omega"),
                           MCMCchains(model_in, "gamma"),
                           MCMCchains(model_in, "shp"))

#Get species geo-indices ranges (95% of records range)
geo_range_spec <- geo_data %>%
  group_by(spec_name) %>%
  summarise(lower=quantile(Ap_21,probs=c(.025)),upper=quantile(Ap_21,probs=c(.975)))

sa_range_spec <- geo_data %>%
  group_by(spec_name) %>%
  summarise(lower=quantile(SS_norm,probs=c(.025)),upper=quantile(SS_norm,probs=c(.975)))

#Set sunspot value
ss_val <- 1.5

#If species should be removed, add here
#FOR HIGH SS ESTIMATES ONLY: Species that have been only captured in recent years have very few records at high
# sunspot counts, so remove those here to avoid showing very extrapolated relationships
to_remove <- as.numeric(which(sa_range_spec$upper<ss_val))
#to_remove <- c()
#Limit to 3 species
not_three_species <- which(!spec_list %in% c("Lark Bunting","Cerulean Warbler","Broad-tailed Hummingbird"))
to_remove <- c(to_remove,not_three_species) 

remove_list <- c(to_remove)

#Set range of normalized Ap values
intra_x <- seq(from = -1.5, to = 5,by=.05)
intra_x_trans <- intra_x * sd(geo_data$Ap_21) + mean(geo_data$Ap_21)



ss_val * sd(geo_data$SS_21) + mean(geo_data$SS_21)


ix_rec <- list()
pred_vagr_rec <- list()

#Set up plot
(par(pty = "s"))
par(las=1)
plot(x=NULL, y=NULL,type="l",ylim=c(2,5),xlim=c(4,30),
     lwd=4,ylab="",xlab="",cex.lab=1.5,cex.axis=1.5)
sig_pos_count <- 0
sig_neg_count <- 0
not_sig_count <- 0

#Limit species list to example species

for (nn in 1:length(spec_list)){
  if (!nn %in% remove_list){
    
    #Get 95% range for Ap
    #Don't plot past the range of the real data
    spec_upper <- geo_range_spec$upper[nn]
    spec_lower <- geo_range_spec$lower[nn]
    xx <- which(intra_x_trans < spec_upper & intra_x_trans > spec_lower)
    intra_x_spec <- intra_x[which(intra_x_trans <spec_upper & intra_x_trans > spec_lower)]
    
    pred.seq.intra_x <- array(dim = c(length(model_chains$gamma), length(intra_x_spec)))
    alpha_ind <- which(colnames(model_chains)==paste("mu_alpha.",nn,".",sep=""))
    beta_ind <- which(colnames(model_chains)==paste("beta.",nn,".",sep=""))
    theta_ind <- which(colnames(model_chains)==paste("theta.",nn,".",sep=""))
    omega_ind <- which(colnames(model_chains)==paste("omega.",nn,".",sep=""))
    shp_ind <- which(colnames(model_chains)==paste("shp.",nn,".",sep=""))
    
    for (i in 1:nrow(model_chains)){
      sens <-  model_chains[[alpha_ind]][i] + model_chains[[beta_ind]][i] * intra_x_spec + 
        model_chains[[theta_ind]][i] * ss_val + model_chains[[omega_ind]][i] * intra_x_spec * ss_val
      spec_shp <- model_chains[[shp_ind]][i]
      exp_vagr <- rep(NA,length(sens))
      
      for (ll in 1:length(sens)){
        exp_vagr[ll] <- qgamma(.5,shape = spec_shp,rate = spec_shp/exp(sens[ll]))
      }
      pred.seq.intra_x[i,] <- exp_vagr
    }
    
    #Get predicted vagrancy
    pred_vagr <- apply(pred.seq.intra_x, 2, mean)
    #Get 95% CrI for slope
    pred_vagr_upper <- apply(pred.seq.intra_x, 2, quantile,probs=c(.975))
    pred_vagr_lower <- apply(pred.seq.intra_x, 2, quantile,probs=c(.025))
    
    sig_pos <- log(pred_vagr_lower[length(pred_vagr_lower)]) - log(pred_vagr_upper[1])
    sig_neg <- log(pred_vagr_lower[1]) - log(pred_vagr_upper[length(pred_vagr_upper)])
    
    xx_in <- intra_x_trans[xx]
    yy_in <- pred_vagr[xx]
    yy_up <- pred_vagr_upper[xx]
    yy_low <- pred_vagr_lower[xx]
    #Highlight species
    if (nn == 29){
      points(xx_in, log(yy_in),type="l",col=alpha("darkorchid3",.6),ylim=c(0,400),xlim=c(0,35),
             lwd=6,ylab="Vagrancy % Increase",xlab="Geomagnetic Disturbance",cex.lab=1.25,cex.axis=1.25,main=spec_list[nn])
      pg_ind <- which(!is.na(yy_in))
      polygon(x = c(xx_in[pg_ind], rev(xx_in[pg_ind])),
              y = c(log(yy_low[pg_ind]), rev(log(yy_up[pg_ind]))),col = alpha("darkorchid3",.3), border = NA)
    } else if (nn == 40){
      points(xx_in, log(yy_in),type="l",col=alpha("deepskyblue3",.6),ylim=c(0,400),xlim=c(0,35),
             lwd=6,ylab="Vagrancy % Increase",xlab="Geomagnetic Disturbance",cex.lab=1.25,cex.axis=1.25,main=spec_list[nn])
      pg_ind <- which(!is.na(yy_in))
      polygon(x = c(xx_in[pg_ind], rev(xx_in[pg_ind])),
              y = c(log(yy_low[pg_ind]), rev(log(yy_up[pg_ind]))),col = alpha("deepskyblue3",.3), border = NA)
    } else if (nn == 73){
      points(xx_in, log(yy_in),type="l",col=alpha("black",.6),ylim=c(0,400),xlim=c(0,35),
             lwd=6,ylab="Vagrancy % Increase",xlab="Geomagnetic Disturbance",cex.lab=1.25,cex.axis=1.25,main=spec_list[nn])
      pg_ind <- which(!is.na(yy_in))
      polygon(x = c(xx_in[pg_ind], rev(xx_in[pg_ind])),
              y = c(log(yy_low[pg_ind]), rev(log(yy_up[pg_ind]))),col = alpha("black",.3), border = NA)
    } else {
      if (sig_pos > 0){
        points(intra_x_trans[xx], log(pred_vagr[xx]),type="l",col=alpha("darkorange",.1),ylim=c(0,400),xlim=c(0,35),
               lwd=6,ylab="Vagrancy % Increase",xlab="Geomagnetic Disturbance",cex.lab=1.25,cex.axis=1.25,main=spec_list[nn])
        sig_pos_count <- sig_pos_count + 1
      } else if (sig_neg > 0){
        points(intra_x_trans[xx], log(pred_vagr[xx]),type="l",col=alpha("springgreen3",.1),ylim=c(0,400),xlim=c(0,35),
               lwd=6,ylab="Vagrancy % Increase",xlab="Geomagnetic Disturbance",cex.lab=1.25,cex.axis=1.25,main=spec_list[nn])
        sig_neg_count <- sig_neg_count + 1
      } else{
        points(intra_x_trans[xx], log(pred_vagr[xx]),type="l",col=alpha("gray67",.1),ylim=c(0,400),xlim=c(0,35),
               lwd=6,ylab="Vagrancy % Increase",xlab="Geomagnetic Disturbance",cex.lab=1.25,cex.axis=1.25,main=spec_list[nn])
        not_sig_count <- not_sig_count + 1
      }
    }
    ix_rec <- c(ix_rec,intra_x_trans[xx])
    pred_vagr_rec <- c(ix_rec,intra_x_trans[xx])
  }
}

pred.seq.intra_x <- array(dim = c(length(model_chains$gamma), length(intra_x)))
for (i in 1:nrow(model_chains)){
  sens <-  model_chains$gamma[i] + model_chains$mu_beta[i] * intra_x + 
    model_chains$mu_theta[i] * ss_val + model_chains$mu_omega[i] * intra_x * ss_val
  mu_shp <- model_chains$mu_shp[i]
  #Draw from Gamma, for reference:
  #lp = alpha[sp_yr] + beta[ii] .* x_ap + theta[ii] .* x_ss + omega[ii] .* x_ap .* x_ss; 
  #shp_over_lp_exp = shp[ii] ./ exp(lp);
  #y ~ gamma(shp[ii], shp_over_lp_exp);
  exp_vagr <- rep(NA,length(sens))
  for (ll in 1:length(sens)){
    exp_vagr[ll] <- qgamma(.5,shape = mu_shp,rate = mu_shp/exp(sens[ll]))
  }
  pred.seq.intra_x[i,] <- exp_vagr
}

pred_vagr <- apply(pred.seq.intra_x, 2, mean)

points(intra_x_trans[1:90], log(pred_vagr[1:90]),type="l",col=alpha("seagreen4",.9),ylim=c(20,60),
       lwd=7,ylab="Vagrancy",xlab="Geomagnetic Disturbance",cex.lab=1.75,cex.axis=1.75,lty=5)
lower_25 <- apply(pred.seq.intra_x, 2, quantile, probs = 0.025)
upper_975 <- apply(pred.seq.intra_x, 2, quantile, probs = 0.975)
#lower_25 <- (lower_25 - vagr_0)/vagr_0
#upper_975 <- (upper_975 - vagr_0)/vagr_0

polygon(x = c(intra_x_trans[1:90], rev(intra_x_trans[1:90])),
        y = c(log(lower_25[1:90]), rev(log(upper_975[1:90]))),col = alpha("seagreen4",.3), border = NA)

#Compute slope of overall line with CI
(log(pred_vagr[90])-log(pred_vagr[1]))/(intra_x_trans[90]-intra_x_trans[1])

#Upper CI
(log(upper_975[90])-log(lower_25[1]))/(intra_x_trans[90]-intra_x_trans[1])

#Lower CI
(log(lower_25[90])-log(upper_975[1]))/(intra_x_trans[90]-intra_x_trans[1])
