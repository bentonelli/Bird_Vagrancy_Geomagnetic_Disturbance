#### Functions ####

### Function for creating dataframes for species, adding vagrancy index,seasons and geo data
get_spec_VI_geo_season <- function(targ_spec,all_banding_points,geo_data){
  library(ebirdst)
  library(dplyr)
  library(readr)
  library(lubridate)
  library(sp)
  library(raster)
  library(sf)
  library(circular)
  library(geosphere)
  library(nngeo)
  library(pryr)
  #Import the banding data for target species
  spec_targ_b <- filter(all_banding_points,SPECIES_NAME == targ_spec)
  
  #Add week of the year for comparison to ebird data
  spec_targ_b$WEEK <- week(spec_targ_b$DATE)
  
  #Download the abundance data for a target species
  sp_path <- ebirdst_download(species = targ_spec,force=TRUE)
  
  abd <- get_species_path(targ_spec) %>% 
    load_raster("abundance", .)
  #Set raster options up front
  rasterOptions(todisk = FALSE,maxmemory = 8e+09,chunksize = 5e+08,memfrac = .8) 
  
  ### Construct weekly relative abundance maps ####
  
  spec_targ_b_week_final <- data.frame(NULL)
  #For each week,
  for (each_week in 1:52){
    
    print(each_week)
    #Split banding data to the week
    spec_targ_b_week <- filter(spec_targ_b,WEEK==each_week)
    
    #Get locations of banding data, check if they are outside the eBird model range
    locations <- as.data.frame(dplyr::select(spec_targ_b_week,LON_DECIMAL_DEGREES,LAT_DECIMAL_DEGREES))
    
    print("Splitting Data")
    #Split weekly eBird abundance data
    w_abd_r <- abd[[each_week]]
    
    #Rarely, eBird may have no modeled occupancy during a week, in this case the max value will be zero, and we can skip the rest
    if (maxValue(w_abd_r) > 0){
      print("Trimming")
      #print(mem_used())
      w_abd_r <- trim(w_abd_r)
      canProcessInMemory(w_abd_r,verbose=TRUE)
      print("Projecting")
      #Project to workable lat/lon format
      w_abd_r_proj <- projectRaster(w_abd_r, crs = "+init=epsg:4326",method = "ngb")
      
      print("writing to points")
      #Create relative abundance csv in lat/long/abund format
      rel_abd_spatial <- rasterToPoints(w_abd_r_proj)
      colnames(rel_abd_spatial) <- c("longitude", "latitude", "abundance_umean")
      rel_abd_spatial <- as.data.frame(rel_abd_spatial)
      print("Doing the rest...")

      #Remove 0s
      rel_abd_spatial <- filter(rel_abd_spatial,abundance_umean > 0)
      
      #Create column for relative abundance
      rel_abd_spatial$perc_occupancy <- rel_abd_spatial$abundance_umean/sum(rel_abd_spatial$abundance_umean)
      
      #Pull 10k points to represent weekly abundance
      nn_point_centers <- rel_abd_spatial[sample(1:nrow(rel_abd_spatial),10000,replace=TRUE,prob=rel_abd_spatial$perc_occupancy),1:2]
      
      #These are just center points of 2.96x2.96km squares, so we want to expand the points by a random factor
      #Use a fast estimate of distance here, note that this only causes very minor distortion
      # 111,111 meters = 1 degree of latitude, and 111,111 meters * cos(latitude) = 1 degree longitude, so correct for that here
      rand_lon_change <- runif(nrow(nn_point_centers),-1480,1480)/(111111*cos(rad(nn_point_centers$latitude)))
      rand_lat_change <- runif(nrow(nn_point_centers),-1480,1480)/111111
      
      #Also note how minuscule these changes are...
      
      nn_point_centers$longitude <- nn_point_centers$longitude + rand_lon_change
      nn_point_centers$latitude <- nn_point_centers$latitude + rand_lat_change
      
      #Filter banding/encounter data outside model range
      #Get locations of banding data, check if they are outside the eBird model range
      
      brt <- as.data.frame(extract(w_abd_r_proj,locations))
      outside_range_index <- which(is.na(brt$V1))
      
      #Eliminate those records from the dataset
      if (length(outside_range_index) > 0){
        spec_targ_b_week <- spec_targ_b_week[-outside_range_index,]
        paste(length(outside_range_index)," points eliminated for being outside range",sep="")
      }
      
      #Create vagrancy index values
      sp_band_points <- dplyr::select(spec_targ_b_week,LON_DECIMAL_DEGREES,LAT_DECIMAL_DEGREES)
      spec_targ_b_week$vagrancy_index <- ref_to_NN_group(nn_point_centers,sp_band_points,10)
      
      #Add to all week dataframe
      spec_targ_b_week_final <- rbind(spec_targ_b_week_final,spec_targ_b_week)
    }
  }
  
  ### ADD geo to species frame ###
  
  row_matches <- match(spec_targ_b_week_final$DATE,geo_data$DATE)
  spec_w_geo <- cbind(spec_targ_b_week_final,geo_data[row_matches,4:13])
  
  ### ADD seasons ###
  
  #Import season dates as weeks
  season_weeks <- week(as.matrix(ebirdst_runs[which(ebirdst_runs$common_name==targ_spec),c(8,9,12,13,16,17,20,21)]))
  #1=Breeding,2=nonbreeding,3=postbreeding,4=prebreeding
  
  #Correct for missing seasons. Important to note that these seasons are not modeled because of data quality issues
  #In this dataset, only one season is missing from datasets, and only b/nb are missing
  count <- 0
  for (date in season_weeks){
    count <- count + 1
    if(is.na(date)){
      if (count == 1){
        #breeding start = prebreeding + 1
        season_weeks[1] <- season_weeks[8]+1
      } else if (count == 2){
        #breeding end = post breeding start - 1
        season_weeks[2] <- season_weeks[5]-1
      } else if (count == 3){
        #nb start = post breeding end + 1
        season_weeks[3] <- season_weeks[6]+1
      } else if (count == 4){
        #nb end = pre breeding start - 1
        season_weeks[4] <- season_weeks[7]-1
      }
    }
  }
  
  #Find the season that overlaps the new year, mostly should be non-breeding
  for (n in 1:4){
    season_overlap <- season_weeks[(n*2-1)] > season_weeks[(n*2)]
    if (season_overlap == TRUE){
      if (n == 1){
        season_targ <- "breeding"
      } else if (n == 2){
        season_targ <- "nonbreeding"
      } else if (n == 3){
        season_targ <- "postbreeding_migration"
      } else if (n == 4){
        season_targ <- "prebreeding_migration"
      }
    }
  }
  #Create table of weeks, seasons
  wk_num <- 1:52
  seas_wks <- c()
  for (wks in wk_num){
    if (wks >= season_weeks[1] & wks <= season_weeks[2]){
      #breeding
      seas_wks <- c(seas_wks,"breeding")
    } else if (wks >= season_weeks[3] & wks <= season_weeks[4]){
      #non-breeding
      seas_wks <- c(seas_wks,"nonbreeding")
    } else if (wks >= season_weeks[5] & wks <= season_weeks[6]){
      #post-breeding
      seas_wks <- c(seas_wks,"postbreeding_migration")
    } else if (wks >= season_weeks[7] & wks <= season_weeks[8]){
      #pre-breeding
      seas_wks <- c(seas_wks,"prebreeding_migration")
    } else {
      #Overlap season
      seas_wks <- c(seas_wks,season_targ)
    }
  }
  
  wk_w_season <- as.data.frame(cbind(wk_num,seas_wks))
  
  #Match weeks to banding data
  row_matches <- match(spec_w_geo$WEEK,wk_w_season$wk_num)
  spec_w_geo$season <- wk_w_season[row_matches,2]
  
  
  ### Scrap week, coord precision
  spec_w_geo$COORD_PRECISION <- NULL
  spec_w_geo$WEEK <- NULL
  
  return(spec_w_geo)
}

### NN matching function
ref_to_NN_group <- function(reference_mat,banding_points,x){
  #Set up flagging system for what happens when there are 0 banding points
  if (nrow(banding_points)==0){
    nn_record <- c()
    return(nn_record)
  }
  reference_mat <- as.matrix(reference_mat)
  banding_points <- as.matrix(banding_points)
  #For each banding point
  nn_record <- c()
  for (n in 1:nrow(banding_points)){
    #pull in banding points
    test_point <- banding_points[n,]
    #Create list of distance to all points
    comp_nn <- spDistsN1(reference_mat,test_point,longlat=TRUE)
    #Get average of all x number of NN
    nn_avg <- mean(sort(comp_nn,decreasing = FALSE)[1:x])
    #Save to vector
    nn_record <- c(nn_record,nn_avg) 
  } 
  return(nn_record)
}


#Here is the code to process a list of species, attach geo data, calculate VI
spec_list <- read_csv("Data/Spec_data/Fall_spec_traits.csv", 
                                   col_names = FALSE)
geo_data <- read_csv("Data/Geo_data/geo_data.csv")
All_Spec_B_E_Combined <- read_csv("Data/Compiled_Banding_Records/All_Spec_B_E_Combined.csv")

spec_list <- spec_list$X1

for (targ_spec in spec_list){
  print(targ_spec)
  spec_df <- get_spec_VI_geo_season(targ_spec,All_Spec_B_E_Combined,geo_data)
  #write_csv(spec_df,paste(gsub(" ", "_", targ_spec),"_processed.csv",sep=""))
}

