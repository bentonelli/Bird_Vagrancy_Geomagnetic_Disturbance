#Figure 2

# Plot for species-specific vagrancy, ap, and ss indices
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
library(dichromat)
library(cowplot)

#Get map plots of single bird
targ_spec  <- "Cerulean Warbler"

#get the global average vagrancy first
all_data <- readRDS("Fit_models_and_data/Data/Fall_interaction_df.rds")
avg_ap <- mean(all_data$Ap_21)
avg_ss <- mean(all_data$SS_21)

spec_file <- read_csv(paste("Data/Processed_Species_Data/",gsub(" ", "_", targ_spec),"_","processed.csv",sep=""))

spec_file <- filter(spec_file,season=="postbreeding_migration")
#Sort randomly to avoid sorting, and plotting, by date
rows <- sample(nrow(spec_file))
spec_file <- spec_file[rows,]
nrow(spec_file)

world <- ne_countries(scale = "medium", returnclass = "sf")

mappoints <- ggplot(data = world) +
  geom_sf() +
  coord_sf(crs = "+init=epsg:3035")

par(mfrow=c(3,1))
#Vagrancy plot
p1 <- ggplot(data = world) +
  geom_sf(fill= "white") +
  theme(panel.background = element_rect(fill = "grey95"),panel.grid = element_blank()) +
  coord_sf(xlim = c(-125, -65), ylim = c(22, 60), expand = TRUE) +
  geom_point(data = spec_file,aes(x=LON_DECIMAL_DEGREES,y=LAT_DECIMAL_DEGREES,color=log(vagrancy_index)),size=6,alpha=.7) +
  xlab("") + ylab("") + scale_color_viridis_c(name="Vagrancy Index (log)",option="viridis",direction=-1,limits=c(0,8.6)) + 
  theme(axis.text.x=element_blank(),axis.text.y=element_blank(),axis.ticks=element_blank(),legend.position="none") +
  theme(legend.key.size = unit(1.75, 'cm')) +
  theme(legend.title = element_blank()) + 
  theme(legend.text = element_text(size=18)) +
  theme(plot.margin = unit(c(0,0,0,0), "cm"))

#Now plot Ap index
p2 <- ggplot(data = world) +
  geom_sf(fill= "white") +
  theme(panel.background = element_rect(fill = "grey95"),panel.grid = element_blank()) +
  coord_sf(xlim = c(-125, -65), ylim = c(22, 60), expand = TRUE) +
  geom_point(data = spec_file,aes(x=LON_DECIMAL_DEGREES,y=LAT_DECIMAL_DEGREES,color=Ap_21),size=6,alpha=.7) +
  xlab("") + ylab("") + 
  scale_colour_gradient2(low = "burlywood", mid = "lightgrey", high = "deepskyblue4", midpoint = avg_ap,limits=c(0,42)) +
  #scale_color_viridis_c(name="Geomagnetic Disturbance",option="inferno",direction=1) + 
  theme(axis.text.x=element_blank(),axis.text.y=element_blank(),axis.ticks=element_blank(),legend.position="none") +
  theme(legend.key.size = unit(1.75, 'cm')) +
  theme(legend.title = element_blank()) + 
  theme(legend.text = element_text(size=18)) +
  theme(plot.margin = unit(c(0,0,0,0), "cm"))

#Now plot SS index
p3 <- ggplot(data = world) +
  geom_sf(fill= "white") +
  theme(panel.background = element_rect(fill = "grey95"),panel.grid = element_blank()) +
  coord_sf(xlim = c(-125, -65), ylim = c(22, 60), expand = TRUE) +  
  geom_point(data = spec_file,aes(x=LON_DECIMAL_DEGREES,y=LAT_DECIMAL_DEGREES,color=SS_21),size=6,alpha=.7) +
  xlab("") + ylab("") + 
  scale_colour_gradient2(low = "burlywood4", mid = "lightgrey", high = "tomato2", midpoint = avg_ss,limits=c(0,222)) +
  theme(axis.text.x=element_blank(),axis.text.y=element_blank(),axis.ticks=element_blank(),legend.position="none") +
  theme(legend.key.size = unit(1.75, 'cm')) +
  theme(legend.title = element_blank()) + 
  theme(legend.text = element_text(size=18)) +
  theme(legend.text.align = 0) +
  theme(plot.margin = unit(c(0,0,0,0), "cm"))

plot_grid(p1, p2,p3,nrow = 3)

