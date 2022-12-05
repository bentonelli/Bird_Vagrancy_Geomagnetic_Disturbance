#Look at non-breeding latitude sensitivity to geomag. disturbance 

library(MCMCvis)
library(readr)
library(ggplot2)
library(qdapRegex)

spring_list <- readRDS("Fit_models_and_data/Data/Spring_Ap_df.rds")
spring_list <- unique(spring_list$spec_name)

spring_data <- read_csv("Data/Spec_data/Spring_spec_traits.csv")
colnames(spring_data)[1] <- "spec" 

s_model <- readRDS("Fit_models_and_data/Models/Spring_Ap.rds")

#Get betas
s_betas <- MCMCsummary(s_model,params="mu_beta")

#Merge medians
s_medians <- as.data.frame(s_betas$`50%`)

s_medians$spec <- spring_list

#Get the NB lat estimate 
filenames <- list.files("Data/Spec_NB_lats/", pattern="*.rds", full.names=TRUE)
spec_list <- unlist(qdapRegex::ex_between(filenames, "Data/Spec_NB_lats//", ".rds"))
spec_list <- gsub("_", " ", spec_list)

spec_lat_rec <- c()
for (each_file in filenames){
  spec_centroids <- readRDS(each_file)
  spec_centroids <- as.data.frame(spec_centroids)
  spec_nb <- spec_centroids[which(spec_centroids$V2 == "nonbreeding"),]
  spec_mean_lat <- mean(as.numeric(spec_nb$V4))
  spec_mean_lon <- mean(as.numeric(spec_nb$V3))  
  spec_lat_rec <- c(spec_lat_rec,spec_mean_lat)
}
spec_lat_rec <- as.data.frame(spec_lat_rec)
spec_lat_rec$spec <- spec_list

spr <- merge(spring_data,s_medians)
colnames(spr)[6] <- "mu_beta"

spec_w_nb <- merge(spr,spec_lat_rec,by="spec",all=TRUE)

#A Few species need to be manually estimated, because eBird doesn't estimate them,
#or the code to get nb centroids otherwise failed.

#Black-billed Cuckoo
spec_w_nb$spec_lat_rec[11] <- -2.5
#Cassin's Vireo
spec_w_nb$spec_lat_rec[36] <- 23
#Cerulean Warbler
spec_w_nb$spec_lat_rec[37] <- 4
#Connecticut Warbler
spec_w_nb$spec_lat_rec[42] <- -2.5
#Purple Martin
spec_w_nb$spec_lat_rec[91] <- -7

#Remove CAHU which wasn't in this spring dataset
spec_w_nb <- spec_w_nb[-c(32),]

#Fit segmented regression
library(SiZer)
library(ggplot2)

model <- piecewise.linear(spec_w_nb$spec_lat_rec,spec_w_nb$mu_beta, CI=FALSE)
#group according to the changepoint
spec_w_nb$grp = factor(ifelse(spec_w_nb$spec_lat_rec > model$change.point,1,0))

ggplot(spec_w_nb,aes(x=spec_lat_rec,y=mu_beta,group=grp),size=4) + 
  theme_bw() + theme(panel.grid.major = element_blank(),
                     panel.grid.minor = element_blank(),
                     axis.text=element_text(size=18),text = element_text(size=60),
                     axis.title=element_text(size=18,face="bold")) +
  xlab("Non-breeding Latitude") + ylab("Sensitivity to Geom. Disturbance") +
  theme(plot.margin = unit(c(1,1,1,1), "cm")) +
  theme(aspect.ratio=1) +
  geom_vline(xintercept = 24,color="firebrick3", 
             linetype="dashed", size=1.5) +
  geom_point(size = 4,col=alpha("dodgerblue4",.6)) + 
  geom_smooth(method="lm",formula=y~x,col="black")

