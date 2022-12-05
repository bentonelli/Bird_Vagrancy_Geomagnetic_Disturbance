#Plot to show the distribution of records from fall and spring seasons by date

library(readr)
library(dplyr)
library(lubridate)
library(ebirdst)
library(ggplot2)

#Read in fall data
fall_data <- read_csv("Data/Spec_data/Fall_spec_traits.csv")
spec_list <- fall_data$`Species Name`


days_since_rec_f <- c()
for (targ_spec in spec_list){
  #Read in species file
  spec_file <- read_csv(paste("Data/Processed_Species_Data/",
                              gsub(" ", "_", targ_spec),"_","processed.csv",sep=""))
  #Get the ebird dates for fall migration
  migr_dates <- yday(as.matrix(ebirdst_runs[which(ebirdst_runs$common_name==targ_spec),16:17]))
  spec_file <- filter(spec_file,yday(DATE) >= migr_dates[1] & yday(DATE) <= migr_dates[2])
  days_since <- yday(spec_file$DATE) - migr_dates[1]
  days_since_rec_f <- c(days_since_rec_f,days_since)
}


#For spring
#Read in spring data
spring_data <- read_csv("Data/Spec_data/Spring_spec_traits.csv")
spec_list <- spring_data$`Spec name`

days_since_rec_s <- c()
for (targ_spec in spec_list){
  #Read in species file
  spec_file <- read_csv(paste("Data/Processed_Species_Data/",
                              gsub(" ", "_", targ_spec),"_","processed.csv",sep=""))
  #Get the ebird dates for fall migration
  migr_dates <- yday(as.matrix(ebirdst_runs[which(ebirdst_runs$common_name==targ_spec),20:21]))
  spec_file <- filter(spec_file,yday(DATE) >= migr_dates[1] & yday(DATE) <= migr_dates[2])
  days_since <- yday(spec_file$DATE) - migr_dates[1]
  days_since_rec_s <- c(days_since_rec_s,days_since)
}

# Now compare the two, fall, and spring
par(mfrow=c(2,1))
par(mar = c(5, 5, 1, 5))
par(las=1)
options(scipen=5)
bp <- seq(0,175,by=7)
hist(days_since_rec_s,breaks=bp,xlim=c(0,170),col=alpha("springgreen4",.8),main="",xlab="",ylab="",cex.axis=1.25,right=F)
abline(v=median(days_since_rec_s),lty=3,lwd=6)
hist(days_since_rec_f,breaks=bp,xlim=c(0,170),col=alpha("sienna3",.8),main="",xlab="",ylab="",cex.axis=1.25,right=F)
abline(v=median(days_since_rec_f),lty=3,lwd=6)
median(days_since_rec_s)
median(days_since_rec_f)
