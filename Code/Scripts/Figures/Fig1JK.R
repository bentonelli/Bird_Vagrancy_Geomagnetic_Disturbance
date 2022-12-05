# Figure 1:
# Long term vagrancy plots

library(dplyr)
library(readr)
library(ggplot2)

#Read in data
spec_data <- readRDS("Fit_models_and_data/Data/Fall_Ap_df.rds")
spec_traits <- read_csv("Data/Spec_data/Fall_spec_traits.csv")
geo_data <- read_csv("Data/Geo_data/geo_data.csv")

par(mfrow=c(2,1),mai = c(.5, .5, 0.2, 1),mgp = c(4, 1, 0))
head(spec_data)


#Get fall data
par(las=1)
plot(geo_data$DATE,geo_data$Ap_21,type="l",col=alpha("skyblue3",.75),xlab="",cex.axis=1.5,cex.lab=1.5,yaxt="n",ylab="")
#par(las=3)
axis(4,cex.axis=1.5)
par(las=3)
par(las=1)
plot(geo_data$DATE,geo_data$SS_21,type="l",col=alpha("tomato3",.75),xlab="Year",cex.axis=1.5,cex.lab=1.5,yaxt="n",ylab="")
axis(4,cex.axis=1.5)
par(las=3)



