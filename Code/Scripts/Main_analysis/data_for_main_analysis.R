#Process species-specific dataframes to create new 

library(MCMCvis)
library(coda)
library(rstan)
library(readr)
library(rlist)
library(dplyr)
library(lubridate)
library(ebirdst)
library(parallel)

#Get species list with traits
spec_traits <- read_csv("Data/Spec_data/Fall_spec_traits.csv", col_names = TRUE)

#Optional to trim dataset to a certain number of species - for testing
#spec_traits <- spec_traits[1:40,]

#Get Splits/Lumps exclusion list
excl_list <- read_csv("Data/Spec_data/S_L_Exclusions.csv", col_names = TRUE)

# Read in geo data
geo_data <- read_csv("Data/Geo_data/geo_data.csv")

#Add info to column names
colnames(spec_traits)[1:4] <- c("spec_name","breed_lat","migr_len","taxon") 

#Sort alphabetically
spec_traits <- spec_traits %>% arrange(spec_name)

#Save spec list separately
spec_list <- spec_traits$spec_name

#Normalize species trait variables to percentage of mean
spec_traits$migr_len <- (spec_traits$migr_len - mean(spec_traits$migr_len))/sd(spec_traits$migr_len)
spec_traits$breed_lat <- (spec_traits$breed_lat - mean(spec_traits$breed_lat))/sd(spec_traits$breed_lat)

#Define season
season_of_interest <- "fall"

#Define metric - Ap or SS
metric_of_interest <- "SS"

### Preprocess data, combine records from all species into single dataframe ####

#DF for banding records data
all_spec_df <- c()
#DF for average SS/Ap by year
geo_seas_yr <- c()

counter <- 0
#For each species, add to working df
for (targ_spec in spec_list){
  counter <- counter + 1
  print(targ_spec)
  
  #Read in species file
  spec_file <- read_csv(paste("Data/Processed_Species_Data/",gsub(" ", "_", targ_spec),"_","processed.csv",sep=""))
  
  #Delete records based on dates on exclusion list
  if (targ_spec %in% excl_list$`Spec name`){
    print(targ_spec)
    #Get exlcusion start, end
    excl_rw <- which(excl_list$`Spec name` == targ_spec)
    exc_start <- as.numeric(excl_list[excl_rw,"Exclude_start"])
    exc_end <- as.numeric(excl_list[excl_rw,"Exclude_end"])
    
    #remove those years from the spec file
    spec_file <- filter(spec_file,year(spec_file$DATE) < exc_start | year(spec_file$DATE) > exc_end)
  }
  
  # Pick season of interest - note the season needs to be defined below as well
  if (season_of_interest == "fall"){
    # Get migration timing (fall = 16,17;spring=20,21
    migr_dates <- yday(as.matrix(ebirdst_runs[which(ebirdst_runs$common_name==targ_spec),16:17]))
  } else {
    migr_dates <- yday(as.matrix(ebirdst_runs[which(ebirdst_runs$common_name==targ_spec),20:21]))
  }
  
  #If season skips through new year, things need to be done differently 
  #First, filter data to migration season
  if (migr_dates[1] < migr_dates[2]){
    #Filter data to season
    spec_data <- filter(spec_file,yday(DATE) >= migr_dates[1] & yday(DATE) <= migr_dates[2])
    geo_data_loaded <- filter(geo_data,yday(DATE) >= migr_dates[1] & yday(DATE) <= migr_dates[2])
    spec_data$YEAR <- year(spec_data$DATE)
    geo_data_loaded$YEAR <- year(geo_data_loaded$DATE)
  } else {
    print("Mult-year season")
    spec_data <- filter(spec_file,yday(DATE) >= migr_dates[1] | yday(DATE) <= migr_dates[2])
    geo_data_loaded <- filter(geo_data,yday(DATE) >= migr_dates[1] | yday(DATE) <= migr_dates[2])
    
    #Because this season crosses over the new year, an adjustment needs to be made to the year,
    #such that it reflects the season. In this case, for each date that is in the new(er) year,
    #subtract one year. For an example, if the season is between Dec 2018 and Jan 2019, year for 
    #all points should be 2018
    rep_ind <- yday(spec_data$DATE) < migr_dates[1]
    spec_data$YEAR <- year(spec_data$DATE)
    spec_data$YEAR[rep_ind] <- spec_data$YEAR[rep_ind] - 1
    
    #Same with geo information
    rep_ind <- yday(geo_data_loaded$DATE) < migr_dates[1]
    geo_data_loaded$YEAR <- year(geo_data_loaded$DATE)
    geo_data_loaded$YEAR[rep_ind] <- geo_data_loaded$YEAR[rep_ind] - 1
    
    #Lastly, change the age code - birds in the new year that are marked as second year birds are
    #effectively hatch year birds (because of the way the banding codes code by calendar year)
    rep_ind <- yday(spec_data$DATE) < migr_dates[1] & spec_data$AGE_CODE == 5
    spec_data$AGE_CODE[rep_ind] <- 2
  }
  
  #Next up is to get SS/Ap as the SS/Ap anomaly, now SS/Ap is defined by its value relative to the year
  #Make two new dfs, where SS_21/Ap_21 is scaled based on the year. 
  
  if (metric_of_interest == "Ap"){
    print("Metric: Ap")
    
    #Select the important info from the species data
    vagrancy_anomaly <- dplyr::select(spec_data,vagrancy_index,Ap_21,YEAR,AGE_CODE)
    
    #Get year index
    vagrancy_anomaly$YEAR_ind <- vagrancy_anomaly$YEAR - 1959
    
  } else {
    print("Metric: SS")
    
    #Select the important info from the species data
    vagrancy_anomaly <- dplyr::select(spec_data,vagrancy_index,SS_21,YEAR,AGE_CODE)
    
    #Get year index
    vagrancy_anomaly$YEAR_ind <- vagrancy_anomaly$YEAR - 1959
  }
  
  #Add species
  vagrancy_anomaly$spec_name <- targ_spec
  
  
  if (season_of_interest == "fall"){
    #Filter by age - unknowns (0) get assigned NA, after hatching year (1) and second year plus (5,6,7,8)
    #become 0, hatching year (2) becomes 1
    
    #Filter out local and juvenile codes ages - obsolete codes/not migratory birds
    vagrancy_anomaly <- filter(vagrancy_anomaly, vagrancy_anomaly$AGE_CODE != 3)
    vagrancy_anomaly <- filter(vagrancy_anomaly, vagrancy_anomaly$AGE_CODE != 4)
    
    #Get locations of hatch year birds, non-hatch year birds
    unk_birds <- which(vagrancy_anomaly$AGE_CODE %in% c(0))
    nhy_ind <- which(vagrancy_anomaly$AGE_CODE %in% c(1,5,6,7,8))
    hy_ind <- which(vagrancy_anomaly$AGE_CODE %in% c(2))
    
    #Change age codes
    vagrancy_anomaly[nhy_ind,'AGE_CODE'] <- 0
    vagrancy_anomaly[hy_ind,'AGE_CODE'] <- 1
    vagrancy_anomaly[unk_birds,'AGE_CODE'] <- NA
    
  } else {
    #Filter by age - unknowns (0) and after hatching year (1) get assigned NA,and third year plus (6,7,8)
    #become 0, hatching year (2) becomes 1
    
    #Filter out local and juvenile codes ages - obsolete codes/not migratory birds
    #Also filter out hatching year birds, if these records exist (shouldn't be possible)
    vagrancy_anomaly <- filter(vagrancy_anomaly, vagrancy_anomaly$AGE_CODE != 2)
    vagrancy_anomaly <- filter(vagrancy_anomaly, vagrancy_anomaly$AGE_CODE != 3)
    vagrancy_anomaly <- filter(vagrancy_anomaly, vagrancy_anomaly$AGE_CODE != 4)
    
    #Get locations of hatch year birds, non-hatch year birds - note hy here stands in for SY
    unk_birds <- which(vagrancy_anomaly$AGE_CODE %in% c(0,1))
    nhy_ind <- which(vagrancy_anomaly$AGE_CODE %in% c(6,7,8))
    hy_ind <- which(vagrancy_anomaly$AGE_CODE %in% c(5))
    
    #Change age codes
    vagrancy_anomaly[nhy_ind,'AGE_CODE'] <- 0
    vagrancy_anomaly[hy_ind,'AGE_CODE'] <- 1
    vagrancy_anomaly[unk_birds,'AGE_CODE'] <- NA
  }
  
  #Take a random sample based on size of dataset
  if (nrow(vagrancy_anomaly) > 20000){
    
    num_samples <- 20000
    va_ind <- 1:nrow(vagrancy_anomaly)
    
    #Take sample probabilistically so that the balance between years is as equal as possible.
    yr_ind <- as.data.frame(vagrancy_anomaly$YEAR_ind)
    colnames(yr_ind)[1] <- "YEAR_ind"
    yr_count <- vagrancy_anomaly %>%
      group_by(YEAR_ind) %>%
      summarise(n_yr = n())
    yr_probs <- left_join(yr_ind,yr_count,by="YEAR_ind")
    yr_probs <- 1/yr_probs$n_yr
    yr_probs <- yr_probs/sum(yr_probs)
    
    vagrancy_anomaly_trim <- vagrancy_anomaly[sample(va_ind,num_samples,prob=yr_probs,replace = FALSE),]
    vagrancy_anomaly <- vagrancy_anomaly_trim
  }
  
  #Add to all species df
  all_spec_df <- rbind(all_spec_df,vagrancy_anomaly)
  
}

#Add species names as factors
all_spec_df$spec_name_fct <- as.numeric(as.factor(all_spec_df$spec_name))

