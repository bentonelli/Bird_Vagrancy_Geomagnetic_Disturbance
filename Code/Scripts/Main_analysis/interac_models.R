# Code to run interaction models
# Data here is already "prepackaged"
library(rstan)

#Read in data
all_spec_df <- readRDS("Fit_models_and_data/Data/Fall_interaction_df.rds")

#Create yr,sp index
yr_sp_in <- c()
for (spec_ind in 1:max(all_spec_df$spec_name_fct)){
  yr_sp_in <- c(yr_sp_in,rep(spec_ind,60))
}

vsw_data <- list(
  N = nrow(all_spec_df),
  x_ap = all_spec_df$Ap_norm,
  x_ss = all_spec_df$SS_norm,
  y = all_spec_df$vagrancy_index,
  Nyrs = max(all_spec_df$YEAR_ind),
  tt = as.numeric(as.factor(all_spec_df$YEAR)),
  Nsp = max(all_spec_df$spec_name_fct),
  ii = all_spec_df$spec_name_fct,
  N_sp_yr = max(all_spec_df$spec_name_fct)*max(all_spec_df$YEAR_ind),
  sp_yr = all_spec_df$sp_yr,
  yr_sp_in = yr_sp_in
)
fit1 <- stan(
  file = "Code/Stan_models/fall_interac.stan",  # Stan program
  data = vsw_data,    # named list of data
  chains = 4,             # number of Markov chains
  warmup = 1000,          # number of warmup iterations per chain
  iter = 3000,            # total number of iterations per chain
  cores = 4,              # number of cores (could use one per chain)
  control = list(adapt_delta = .8,
                 max_treedepth = 10,
                 stepsize = .06)
)

#saveRDS(fit1,"")
#saveRDS(all_spec_df,"")