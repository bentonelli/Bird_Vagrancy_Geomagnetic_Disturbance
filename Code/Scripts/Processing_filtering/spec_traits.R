#Code for getting migratory distance and breeding latitude for species

library(dplyr)
library(ebirdst)
library(lubridate)
library(raster)
library(circular)
library(geosphere)
library(rgdal)


spec_list <- c("Yellow-breasted Chat")

migr_dist_all <- c()
breed_lat_all <- c()
for (targ_spec in spec_list){
  print(targ_spec)
  
  #Download the abundance data for a target species
  sp_path <- ebirdst_download(species = targ_spec,force=TRUE)
  
  abd <- get_species_path(targ_spec) %>% 
    load_raster("abundance")
  
  #Set raster options up front
  #rasterOptions(todisk = FALSE,maxmemory = 8e+09,chunksize = 5e+08,memfrac = .8)
  
  #Using the eBird abundance data, calculate the weighted average breeding latitude of the species
  b_weeks <- week(as.matrix(ebirdst_runs[which(ebirdst_runs$common_name==targ_spec),8:9]))
  
  #For seasons that aren't explicitly defined by eBird, we need to find the unmodeled weeks
  #that are in between the other seasons
  if (is.na(b_weeks[1])){
    print("breeding season undefined")
    b_weeks <- week(as.matrix(ebirdst_runs[which(ebirdst_runs$common_name==targ_spec),c(21,16)])) + c(1,-1)
  }
  
  
  if (b_weeks[1] < b_weeks[2]){
    b_weeks <- b_weeks[1]:b_weeks[2]
    
  } else {
    print("Mult-year season")
    b_weeks <- c(b_weeks[1]:52,1:b_weeks[2])
  }
  
  
  
  
  #For each week, draw 10k points, add to list
  all_points_breeding <- c()
  for (each_week in b_weeks){
    
    w_abd_r <- abd[[each_week]]
    
    #Only take data if the week is modeled by eBird
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
      
      #Add to all week dataframe
      all_points_breeding <- rbind(all_points_breeding,nn_point_centers)
    }
  }
  
  #Do same for non-breeding
  all_points_n_breeding <- c()
  nb_weeks <- week(as.matrix(ebirdst_runs[which(ebirdst_runs$common_name==targ_spec),12:13]))
  
  #For seasons that aren't explicitly defined by eBird, we need to find the unmodeled weeks
  #that are in between the other seasons
  if (is.na(nb_weeks[1])){
    print("non-breeding season undefined")
    nb_weeks <- week(as.matrix(ebirdst_runs[which(ebirdst_runs$common_name==targ_spec),c(17,20)])) + c(1,-1)
  }
  
  
  if (nb_weeks[1] < nb_weeks[2]){
    nb_weeks <- nb_weeks[1]:nb_weeks[2]
    
  } else {
    print("Mult-year season")
    nb_weeks <- c(nb_weeks[1]:52,1:nb_weeks[2])
  }
  
  #For each week, draw 10k points, add to list
  all_points_n_breeding <- c()
  for (each_week in nb_weeks){
    
    w_abd_r <- abd[[each_week]]
    
    #Only take data if the week is modeled by eBird
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
      
      #Add to all week dataframe
      all_points_n_breeding <- rbind(all_points_n_breeding,nn_point_centers)
    }
  }
  
  #Get mean long, latitude - breeding
  mean_breeding_lat <- mean(all_points_breeding$latitude)
  mean_breding_lon <- mean(all_points_breeding$longitude)
  migr_distance <- NA
  if (length(all_points_n_breeding)>0){
    
    #Get mean long, latitude - nb
    mean_n_breeding_lat <- mean(all_points_n_breeding$latitude)
    mean_n_breding_lon <- mean(all_points_n_breeding$longitude)
    
    #Get migration distance, kilometers
    migr_distance <- distHaversine(c(mean_breding_lon,mean_breeding_lat),c(mean_n_breding_lon,mean_n_breeding_lat))/1000
    
    
  }
  migr_dist_all <- c(migr_dist_all,migr_distance)
  breed_lat_all <- c(breed_lat_all,mean_breeding_lat)
  
}

#Output all species list
out_spec_traits <- as.data.frame(spec_list)
out_spec_traits$migr_dist <- migr_dist_all

out_spec_traits
out_spec_traits$breed_lat <- breed_lat_all