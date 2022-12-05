#Combine banding, encounter data

library(ebirdst)
library(dplyr)
library(readr)
library(lubridate)

#### Combine all banding, encounter data ####
#Import banding/encounter data - note that the fields are different between the two sets

#Add all banding records together
Y_band_records <- c()

#Add data csvs as variables
#ADD HERE

#Combine all banding files, save only certain columns
band_files_list <- list()
for (band_files in band_files_list){
  band_files <- subset(band_files,select=c(LAT_DECIMAL_DEGREES,LON_DECIMAL_DEGREES,COORD_PRECISION,BANDING_DATE,BIRD_STATUS,AGE_CODE,SEX_CODE,SPECIES_NAME))
  Y_band_records <- rbind(Y_band_records,band_files)
}

#Add column for banding or encounter
Y_band_records$B_or_E <- "B"
#Redo column names
colnames(Y_band_records)[1:8] <- c("LAT_DECIMAL_DEGREES","LON_DECIMAL_DEGREES",
                                   "COORD_PRECISION","DATE","BIRD_STATUS","AGE_CODE",
                                   "SEX_CODE","SPECIES_NAME")
#Add all encounters together 
Y_enc_records <- c()
#Add data csvs as variables
#ADD HERE

enc_files_list <- list()
for (enc_files in enc_files_list){
  #Account for unknown ages in encounter data (not in the dataset, so we just set to unknown here)
  enc_files$AGE_CODE <- 0
  #subset data to the columns needed
  enc_files <- subset(enc_files,select=c(E_LAT_DECIMAL_DEGREES,E_LON_DECIMAL_DEGREES,E_COORD_PRECISION,ENCOUNTER_DATE,B_BIRD_STATUS,AGE_CODE,B_SEX_CODE,B_SPECIES_NAME))
  #Save to master list
  Y_enc_records <- rbind(Y_enc_records,enc_files)
}
Y_enc_records$B_or_E <- "E"
#Rename columns to match banding data
colnames(Y_enc_records)[1:8] <- c("LAT_DECIMAL_DEGREES","LON_DECIMAL_DEGREES",
                                  "COORD_PRECISION","DATE","BIRD_STATUS","AGE_CODE",
                                  "SEX_CODE","SPECIES_NAME")

#Combine banding/encounter records to one MASSIVE df
All_records_unfiltered <- rbind(Y_band_records,Y_enc_records)

#Reformat date column
All_records_unfiltered$DATE <- as.Date(All_records_unfiltered$DATE,"%m/%d/%Y")

#### Filter banding/encounter data ####
Target_spec_list <- read_csv("") #Need to add csv with list of species here!

#Filter based on target species list
t_spec_list <- as.list(Target_spec_list$`Spec name`)
filt_B_E_data <- filter(All_records_unfiltered, SPECIES_NAME %in% t_spec_list)

#Filter data to period 1960-2019
filt_B_E_data <- filter(filt_B_E_data,DATE > as.Date("12/31/1959","%m/%d/%Y") & DATE < as.Date("01/01/2020","%m/%d/%Y"))

#Location needs to be complete (non-zero). There are a few locations that are zeros in the data
filt_B_E_data <- filter(filt_B_E_data, LAT_DECIMAL_DEGREES != 0 & LON_DECIMAL_DEGREES != 0)

#Location accuracy needs to be 0 (exact), 1 (1 minute block), 10 (10 minute block). All other are not precise enough for our purposes.
filt_B_E_data <- filter(filt_B_E_data,COORD_PRECISION == 0 | COORD_PRECISION == 1 | COORD_PRECISION == 10)

#BBS code needs to be 3 - for wild birds, not transported or experimentally altered.
filt_B_E_data <- filter(filt_B_E_data,BIRD_STATUS == 3)

#Locations are specified to either be exact, 1 minute blocks, or 10 minute blocks. For 1 and 10 precision, the lat long is the center point of the block.
#For locations where the precision is not exact, choose a random location within the block and assign that as new location.

#Create random numbers for each row, multiply by the coordinate precision, convert from minutes to decimals, add to existing
filt_B_E_data$LAT_DECIMAL_DEGREES <- runif(nrow(filt_B_E_data),-.5,.5)*filt_B_E_data$COORD_PRECISION/60 + filt_B_E_data$LAT_DECIMAL_DEGREES
filt_B_E_data$LON_DECIMAL_DEGREES <- runif(nrow(filt_B_E_data),-.5,.5)*filt_B_E_data$COORD_PRECISION/60 + filt_B_E_data$LON_DECIMAL_DEGREES

#Write!
#write_csv(filt_B_E_data,"")


