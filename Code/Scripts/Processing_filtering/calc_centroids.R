# Code to get lat-long points for each week of a species
# These can then be used to get the centroids of breeding, non-breeding ranges

library(ebirdst)
library(raster)
library(dplyr)
library(lubridate)

#Using ebirdst, get the centroid of range for each week.
#Download the abundance data for a target species

#Ebird access key is needed!
#set_ebirdst_access_key("89ba0jbnqg6l",overwrite = FALSE)

#Get species list used in analysis - here as fall
spec_list <- readRDS("Fit_models_and_data/Data/Fall_interaction_df.rds")
spec_list <- unique(spec_list$spec_name)

#For each species, 
for (targ_spec in spec_list[1:length(spec_list)]){
  print(targ_spec)
  
  
  sp_path <- ebirdst_download(species = targ_spec,force=TRUE)
  
  abd <- load_raster("abundance", path = sp_path)
  #Set raster options up front
  rasterOptions(todisk = FALSE,maxmemory = 8e+09,chunksize = 5e+08,memfrac = .8)
  
  #Get dates of seasons
  breed_weeks <- week(as.matrix(ebirdst_runs[which(ebirdst_runs$common_name==targ_spec),8:9]))
  nbreed_weeks <- week(as.matrix(ebirdst_runs[which(ebirdst_runs$common_name==targ_spec),12:13]))
  postbreed_weeks <- week(as.matrix(ebirdst_runs[which(ebirdst_runs$common_name==targ_spec),16:17]))
  prebreed_weeks <- week(as.matrix(ebirdst_runs[which(ebirdst_runs$common_name==targ_spec),20:21]))
  
  if (!is.na(nbreed_weeks[1])){
    #Assign to all weeks. Note, this willl get messed up if some season other than nb crosses
    # over into the new year.
    if (nbreed_weeks[1]-nbreed_weeks[2]>0){
      week_season <- rep("nonbreeding",52)
      week_season[prebreed_weeks[1]:prebreed_weeks[2]] <- "prebreeding_migration"
      week_season[breed_weeks[1]:breed_weeks[2]] <- "breeding"
      week_season[postbreed_weeks[1]:postbreed_weeks[2]] <- "postbreeding_migration"
      week_include <- c(nbreed_weeks[1]:52,1:nbreed_weeks[2])
    } else {
      print("Alert")
      week_include <- c(nbreed_weeks[1]:nbreed_weeks[2])
    }
    
    
    week_centroid <- c()
    for (each_week in week_include){
      
      print(each_week)
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
        
        print("Writing to points")
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
        week_centroid <- rbind(week_centroid,c(each_week,week_season[each_week],mean(nn_point_centers$longitude),mean(nn_point_centers$latitude)))
        
      }
    }
    #Save as RDS
    saveRDS(week_centroid,paste(gsub(" ", "_", targ_spec),"_centroids.rds",sep=""))
  } 
}





