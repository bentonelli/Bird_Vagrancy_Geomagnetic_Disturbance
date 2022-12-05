#Figure 3

library(MCMCvis)
library(shinystan)
library(scales)
library(readr)
library(lubridate)
library(dplyr)
library(ggplot2)
library(ebirdst)

#Read in model and data
model_in <- readRDS("Fit_models_and_data/Models/Fall_Ap.rds")
model_df <- readRDS("Fit_models_and_data/Data/Fall_Ap_df.rds")
spec_list <- unique(model_df$spec_name)

#Read in species data
spec_traits_f <- read_csv("Data/Spec_data/Fall_spec_traits.csv", col_names = TRUE)
#Filter to included species
spec_traits_f <- filter(spec_traits_f,spec_traits_f$`Species Name` %in% spec_list)

m_sum_mu_betas <- MCMCsummary(model_in,params="mu_beta")

model_chains <- data.frame(MCMCchains(model_in, "mu_gamma"),
                           MCMCchains(model_in, "psi"),
                           MCMCchains(model_in, "eta"),
                           MCMCchains(model_in, "mu_nu"),
                           MCMCchains(model_in, "mu_shp"),
                           MCMCchains(model_in, "omega"),
                           MCMCchains(model_in,"mu_lambda"))

#Get standardized migration length (as it is formatted in the models)
ml_stand <- (spec_traits_f$`Migration Length` - mean(spec_traits_f$`Migration Length`))/sd(spec_traits_f$`Migration Length`)
xx <- seq(min(ml_stand),max(ml_stand),by=.1)
yy <- mean(model_chains$omega) + mean(model_chains$psi) * xx

#Get 95% CrI
cri_mat <- c()
for (each_x in xx){
  temp <- model_chains$omega + each_x * model_chains$psi
  cri_mat <- cbind(cri_mat,temp)
}

uppers <- c()
lowers <- c()
for (each_col in 1:ncol(cri_mat)){
  uppers <- c(uppers,quantile(cri_mat[,each_col],.975))
  lowers <- c(lowers,quantile(cri_mat[,each_col],.025))
}

par(las=1,pty="s")
breed_lats <- round(spec_traits_f$`Breeding Lat`)
breed_lats <- breed_lats - min(breed_lats) + 1
cols_ramp <- colorRampPalette(c("sienna2","grey","turquoise3"),bias=1)
cols_ramp <- cols_ramp(max(breed_lats))
col_breed_lats<-cols_ramp[breed_lats]

plot(spec_traits_f$`Migration Length`,m_sum_mu_betas$`50%`,pch=19,col=alpha(col_breed_lats,.8),
     xlab="Migration Length (km)",ylab="Sensitivity to Geom. Dist.",cex=2,cex.axis=1.25)
polygon(x = c(xx_trans, rev(xx_trans)),
        y = c((lowers), rev(uppers)),col = alpha("tan4",.25), border = NA)
xx_trans <- xx * sd(spec_traits_f$`Migration Length`) + mean(spec_traits_f$`Migration Length`)
points(xx_trans,yy,type="l",col="tan4",lty=1,lwd="5")
abline(h=0,lty=2,lwd=2,col="darkgrey")

legend("bottomright",legend=rev(c("67°","60°","50°","40°","30°","20°")),cex=1.2,pch=19,col=rev(cols_ramp[c(47,37,27,17,7,1)]),
       title = "Breeding Latitude",horiz = TRUE)


