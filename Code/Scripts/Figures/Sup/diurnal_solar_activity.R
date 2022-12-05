#Diurnal, nocturnal differences in sens. to solar activity

library(readr)
library(lubridate)
library(ebirdst)
library(dplyr)
library(MCMCvis)
library(shinystan)
library(ggpubr)
library(MASS)

spec_list <- read_csv("Data/Spec_data/spec_fall_d_n_class.csv")
model_in <- readRDS("Fit_models_and_data/Models/Fall_interaction.rds")

ddd <- filter(spec_list,`Diurnal Migrant` == "D")
nnn <- filter(spec_list,`Diurnal Migrant` == "N")

t.test(ddd$theta_50,nnn$theta_50)

#Reorder
spec_list$diurn <- factor(spec_list$`Diurnal Migrant` , levels=c("D", "N", "B", "U"))

my_comparisons <- list( c("D","N"))
p <- ggplot(spec_list, aes(x=diurn, y=theta_50)) + 
  theme_bw(base_size = 22) +
  theme(legend.position = "none") +
  theme(axis.text.x = element_text(angle = 45,vjust=.7,hjust=.7)) + 
  scale_x_discrete(labels= c("Diurnal","Nocturnal","Both","Unknown")) + 
  labs(x = "Migration Timing", y = "Effect of Solar Activity on Vagrancy") + 
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank()) +
  geom_boxplot(outlier.shape = NA) +
  geom_jitter(aes(colour = factor(diurn),size=4),position=position_jitter(0.2)) +
  stat_compare_means(comparisons = my_comparisons,label.y=.4,method = "t.test")
p + scale_color_manual(values=alpha(c("goldenrod", "darkblue","seagreen4","honeydew4"),.6))


