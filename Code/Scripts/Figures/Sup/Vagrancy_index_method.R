#Figure 5 (+ a few helpful extra figures!)

# Plot for vagrancy metric workflow
library(readr)
library(ggmap)
library(ebirdst)
library(raster)
library(dplyr)
library(circular)
library(lubridate)
library(rnaturalearth)
library(rnaturalearthdata)
library(ggplot2)

#Get map plots of single bird

targ_spec <- "Lark Bunting"
spec_file <- read_csv(paste("Data/Processed_Species_Data/",gsub(" ", "_", targ_spec),"_","processed.csv",sep=""))

spec_file <- filter(spec_file,season=="postbreeding_migration")

#Sort randomly to avoid sorting, and plotting, by date
rows <- sample(nrow(spec_file))
spec_file <- spec_file[rows,]

world <- ne_countries(scale = "medium", returnclass = "sf")

mappoints <- ggplot(data = world) +
  geom_sf() +
  coord_sf(crs = "+init=epsg:3035")

### Full species data - not just one week like in Fig 5 ####

mappoints <- ggplot(data = world) +
  geom_sf(fill="grey95") +
  theme(panel.background = element_rect(fill = "white"),panel.grid = element_blank()) +
  coord_sf(xlim = c(-130, -62), ylim = c(22, 63), expand = FALSE) +
  geom_point(data = spec_file,aes(x=LON_DECIMAL_DEGREES,y=LAT_DECIMAL_DEGREES,size=6),color="darkblue") +
  xlab("") + ylab("") + theme(axis.text.x=element_blank(),axis.text.y=element_blank(),axis.ticks=element_blank(),legend.position="none") + 
  scale_color_viridis_c(option = "magma") + theme(plot.margin = unit(c(0,0,0,0), "cm"))
mappoints

mappoints <- ggplot(data = world) +
  geom_sf(fill= "grey95") +
  theme(panel.background = element_rect(fill = "white"),panel.grid = element_blank()) +
  coord_sf(xlim = c(-130, -62), ylim = c(22, 63), expand = FALSE) + 
  geom_point(data = spec_file,aes(x=LON_DECIMAL_DEGREES,y=LAT_DECIMAL_DEGREES,color=log(vagrancy_index)),size=6,alpha=.9) +
  xlab("") + ylab("") + scale_color_viridis_c(name="Vagrancy Index (log)",option="viridis",direction=1) + 
  theme(axis.text.x=element_blank(),axis.text.y=element_blank(),axis.ticks=element_blank(),legend.position="bottom") +
  theme(legend.key.size = unit(1.5, 'cm')) +
  theme(legend.title = element_text(size=20)) + 
  theme(legend.text = element_text(size=16)) +
  theme(plot.margin = unit(c(1,1,1,1), "cm"))
mappoints

#Now plot Ap index
mappoints <- ggplot(data = world) +
  geom_sf(fill= "grey95") +
  theme(panel.background = element_rect(fill = "white"),panel.grid = element_blank()) +
  coord_sf(xlim = c(-125, -65), ylim = c(22, 60), expand = FALSE) + 
  geom_point(data = spec_file,aes(x=LON_DECIMAL_DEGREES,y=LAT_DECIMAL_DEGREES,color=Ap_21),size=6,alpha=.8) +
  xlab("") + ylab("") + scale_color_viridis_c(name="Geomagnetic Disturbance",option="inferno",direction=1) + 
  theme(axis.text.x=element_blank(),axis.text.y=element_blank(),axis.ticks=element_blank(),legend.position="bottom") +
  theme(legend.key.size = unit(1.5, 'cm')) +
  theme(legend.title = element_text(size=20)) + 
  theme(legend.text = element_text(size=16)) +
  theme(plot.margin = unit(c(1,1,1,1), "cm"))
mappoints

#Now plot SS index
mappoints <- ggplot(data = world) +
  geom_sf(fill= "grey95") +
  theme(panel.background = element_rect(fill = "white"),panel.grid = element_blank()) +
  coord_sf(xlim = c(-125, -65), ylim = c(22, 60), expand = FALSE) +  
  geom_point(data = spec_file,aes(x=LON_DECIMAL_DEGREES,y=LAT_DECIMAL_DEGREES,color=SS_21),size=6,alpha=.8) +
  xlab("") + ylab("") + scale_color_viridis_c(name="21-Day Sunspot Mean",option="magma",direction=-1) + 
  theme(axis.text.x=element_blank(),axis.text.y=element_blank(),axis.ticks=element_blank(),legend.position="bottom") +
  theme(legend.key.size = unit(1.5, 'cm')) +
  theme(legend.title = element_text(size=22)) + 
  theme(legend.text = element_text(size=18)) +
  theme(plot.margin = unit(c(0,0,0,0), "cm"))
mappoints


### One week plots ####

print(targ_spec)
targ_week <- 33

#Download the abundance data for a target species
#Need ebird access key!
#set_ebirdst_access_key("pef7q1rc52av",overwrite = TRUE)
sp_path <- ebirdst_download(species = targ_spec)

abd <- load_raster("abundance",path=sp_path)
w_abd_r <- abd[[35]]

#Set raster options up front
rasterOptions(todisk = FALSE,maxmemory = 8e+09,chunksize = 5e+08,memfrac = .8)

#Get points
print("Trimming")
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


### Now Plot 10K points ####
mappoints <- ggplot(data = world) +
  geom_sf(fill= "grey95") +
  theme(panel.background = element_rect(fill = "white"),panel.grid = element_blank()) +
  coord_sf(xlim = c(-125, -65), ylim = c(22, 60), expand = FALSE) +  
  geom_point(data = nn_point_centers,aes(x=longitude,y=latitude),size=1,alpha=.25,pch=19,color="gold1") +
  xlab("") + ylab("") + 
  theme(axis.text.x=element_blank(),axis.text.y=element_blank(),axis.ticks=element_blank(),legend.position="bottom") +
  theme(legend.key.size = unit(1.5, 'cm')) +
  theme(legend.title = element_text(size=22)) + 
  theme(legend.text = element_text(size=18)) +
  theme(plot.margin = unit(c(0,0,0,0), "cm"))
mappoints


### Lastly, plot banding points from same week ####
spec_file <- read_csv(paste("Data/Processed_Species_Data/",gsub(" ", "_", targ_spec),"_","processed.csv",sep=""))
spec_file <- filter(spec_file, week(DATE)==targ_week)

mappoints <- ggplot(data = world) +
  geom_sf(fill= "grey95") +
  theme(panel.background = element_rect(fill = "white"),panel.grid = element_blank()) +
  coord_sf(xlim = c(-125, -65), ylim = c(22, 60), expand = FALSE) +  
  geom_point(data = spec_file,aes(x=LON_DECIMAL_DEGREES,y=LAT_DECIMAL_DEGREES),size=6,alpha=.8) +
  theme(axis.text.x=element_blank(),axis.text.y=element_blank(),axis.ticks=element_blank()) +
  xlab("") + ylab("") + 
  theme(plot.margin = unit(c(0,0,0,0), "cm"))
mappoints

mappoints <- ggplot(data = world) +
  geom_sf(fill= "grey95") +
  theme(panel.background = element_rect(fill = "white"),panel.grid = element_blank()) +
  coord_sf(xlim = c(-125, -65), ylim = c(22, 60), expand = FALSE) +  
  geom_point(data = spec_file,aes(x=LON_DECIMAL_DEGREES,y=LAT_DECIMAL_DEGREES,color=log(vagrancy_index)),size=6,alpha=.8) +
  xlab("") + ylab("") + scale_color_viridis_c(name="Vagrancy Index (log)",direction=-1) + 
  theme(axis.text.x=element_blank(),axis.text.y=element_blank(),axis.ticks=element_blank()) +
  theme(legend.position = "none") +
  theme(plot.margin = unit(c(0,0,0,0), "cm"))
mappoints

