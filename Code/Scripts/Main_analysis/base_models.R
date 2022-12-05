#Code to run base models

library(rstan)
library(readr)

#Set metric
metric_of_interest = "Ap"

#Read in data
all_spec_df <- readRDS("Fit_models_and_data/Data/Fall_Ap_df.rds")

#Read in species traits
spec_traits <- read_csv("Data/Spec_data/Fall_spec_traits.csv")
spec_traits <- filter(spec_traits,`Species Name` %in% unique(all_spec_df$spec_name))

#Normalize specis traits
spec_traits$migr_len <- (spec_traits$`Migration Length`-mean(spec_traits$`Migration Length`))/sd(spec_traits$`Migration Length`)
spec_traits$breed_lat <- (spec_traits$`Breeding Lat`-mean(spec_traits$`Breeding Lat`))/sd(spec_traits$`Breeding Lat`)

#Add species names as factors
all_spec_df$spec_name_fct <- as.numeric(as.factor(all_spec_df$spec_name))

#Get missing age indices, use to create an array age_miss that has 0 for non-missing, 1 for missing
age_miss_ind <- which(is.na(all_spec_df$AGE_CODE))
age_miss <- rep(0,length(all_spec_df$AGE_CODE))
age_miss[age_miss_ind]  <- 1

#Reformat age for integration in the model. Missing variables become -1, from NA
age_obs <- all_spec_df$AGE_CODE
age_obs[age_miss_ind]  <- -1

yr_sp_in <- c()
for (spec_ind in 1:max(all_spec_df$spec_name_fct)){
  yr_sp_in <- c(yr_sp_in,rep(spec_ind,60))
}

### Create and run model ####
if (metric_of_interest == "Ap") {
  vsw_data <- list(
    N = nrow(all_spec_df),
    x = all_spec_df$Ap_21, # Note here that Ap_21 is normalized
    y = all_spec_df$vagrancy_index,
    Nyrs = max(all_spec_df$YEAR_ind),
    tt = as.numeric(as.factor(all_spec_df$YEAR)),
    Nsp = max(all_spec_df$spec_name_fct),
    ii = all_spec_df$spec_name_fct,
    migr_len = spec_traits$migr_len,
    breed_lat = spec_traits$breed_lat,
    age_obs = age_obs,
    age_miss = age_miss,
    N_sp_yr = max(all_spec_df$spec_name_fct)*max(all_spec_df$YEAR_ind),
    sp_yr = all_spec_df$sp_yr,
    yr_sp_in = yr_sp_in
  )
} else {
  vsw_data <- list(
    N = nrow(all_spec_df),
    x = all_spec_df$SS_21,
    y = all_spec_df$vagrancy_index,
    Nyrs = max(all_spec_df$YEAR_ind),
    tt = as.numeric(as.factor(all_spec_df$YEAR)),
    Nsp = max(all_spec_df$spec_name_fct),
    ii = all_spec_df$spec_name_fct,
    migr_len = spec_traits$migr_len,
    breed_lat = spec_traits$breed_lat,
    age_obs = age_obs,
    age_miss = age_miss,
    N_sp_yr = max(all_spec_df$spec_name_fct)*max(all_spec_df$YEAR_ind),
    sp_yr = all_spec_df$sp_yr,
    yr_sp_in = yr_sp_in
  )
}


fit1 <- stan(
  file = "Code/Stan_models/fall_base.stan",  # Stan program
  data = vsw_data,    # named list of data
  chains = 4,             # number of Markov chains
  warmup = 1000,          # number of warmup iterations per chain
  iter = 4000,            # total number of iterations per chain
  cores = 4,              # naumber of cores (could use one per chain)
  control = list(adapt_delta = .8,
                 max_treedepth = 10,
                 stepsize=.06),
  init_r = 1
)

#saveRDS(fit1,"")
#saveRDS(all_spec_df,"")
