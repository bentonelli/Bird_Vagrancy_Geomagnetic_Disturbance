# Read in model 1, get predicted increase in vagrancy of +1SD
library(MCMCvis)


model_ap_fall <- readRDS("Fit_models_and_data/Models/Fall_Ap.rds")
model_chains <- MCMCchains(model_ap_fall,params=c("mu_gamma",
                                                  "omega",
                                                  "mu_shp"
                                                  ))
MCMCsummary(model_ap_fall,params=c("mu_gamma",
                                    "omega",
                                    "mu_shp"))
model_chains <- as.data.frame(model_chains)
xx <- seq(-3,3,by=1)

lp <-c()
lp[1] <- mean(model_chains$mu_gamma + 0 * model_chains$omega)
lp[2] <- mean(model_chains$mu_gamma + 2 * model_chains$omega)
mean_mu_shp <- mean(model_chains$mu_shp)

pred_1 <- rgamma(1000000,mean_mu_shp,(mean_mu_shp/exp(lp[1])))
pred_2 <- rgamma(1000000,mean_mu_shp,(mean_mu_shp/exp(lp[2])))

100 * sum(pred_1 > 100)/1000000
100 * sum(pred_2 > 100)/1000000

median(pred_1)
median(pred_2)

par(las=1)
par(pty="s")
par(mar=c(6,6,3,3))
plot(density(pred_1,bw=5),lwd=4,col="black",main="",xlim=c(0,150),
     xlab="Predicted Vagrancy Index",cex.lab=1.2,cex.axis=1)
points(density(pred_2,bw=5),type="l",lwd=4,col="dodgerblue3")
abline(v=100,lty=2,col="grey",lwd=4)
legend("topright",legend = c("Average Geomagnetic Disturbance","+2 Std. Dev. Geomagnetic Disturbance"),lty = 1,lwd=4,col=c("black","dodgerblue3"))

