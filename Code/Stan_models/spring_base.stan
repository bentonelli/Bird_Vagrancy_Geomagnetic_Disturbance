// Spring base
data {
  int<lower=0> N;
  vector[N] x;
  vector[N] y;
  int<lower=1> Nyrs;
  int<lower = 1, upper = Nyrs> tt[N];
  int<lower=1> Nsp;
  int<lower = 1, upper = Nsp> ii[N];
  vector[Nsp] migr_len;
  vector[Nsp] breed_lat;
  vector[N] age_obs;  
  int<lower=1> N_sp_yr;
  int<lower = 1, upper = N_sp_yr> sp_yr[N];
  int<lower = 1, upper = Nsp> yr_sp_in[N_sp_yr];
}

parameters {
  
  real<lower=0> mu_shp;
  real<lower=0> rate_shp;
  vector<lower=0>[Nsp] shp;
  
  vector[N_sp_yr] alpha;
  vector[N_sp_yr] beta;
  
  vector[Nsp] gamma;
  real mu_gamma;
  
  real omega;
  real psi;
  real eta;
  
    real mu_nu;
  real mu_lambda;
  
  real<lower=0> sigma_alpha;
  real<lower=0> sigma_beta;
  real<lower=0> sigma_mu_beta;
  real<lower=0> sigma_nu;
  real<lower=0> sigma_lambda;
  real<lower=0> sigma_gamma;
  
  //For NC params
  vector[Nsp] nu_raw;
  vector[Nsp] lambda_raw;

  vector[Nsp] mu_beta_raw;

}

transformed parameters {
  vector[Nsp] lambda;
  vector[Nsp] nu;

  vector[Nsp] mu_beta;
  vector[Nsp] delta;
  
  //For non-centered params.
  nu = mu_nu + nu_raw * sigma_nu;
  lambda = mu_lambda + lambda_raw * sigma_lambda;
  
  delta = omega + psi * migr_len + eta * breed_lat;
  
  mu_beta = delta + sigma_mu_beta * mu_beta_raw;
  
}

model {
  
  vector[N] lp;
  vector[N] shp_over_lp_exp;
  
  mu_shp ~ normal(1,1);
  rate_shp ~ normal(1,1);
  
  mu_gamma ~ normal(3,1);
  
  sigma_beta ~ std_normal();
  sigma_alpha ~ std_normal();
  sigma_mu_beta ~ std_normal();
  sigma_gamma ~ std_normal();
  sigma_nu ~ std_normal();
  sigma_lambda ~ std_normal();
  
  nu_raw ~ std_normal();
  lambda_raw ~ std_normal();

  mu_beta_raw ~ std_normal();

  gamma ~ normal(mu_gamma,sigma_gamma);
  
  omega ~ std_normal();
  psi ~ std_normal();
  eta ~ std_normal();
  
  mu_lambda ~ std_normal();
  mu_nu ~ std_normal();
  
  shp ~ gamma(mu_shp,rate_shp);
  
  alpha ~ normal(gamma[yr_sp_in],sigma_alpha);
  beta ~ normal(mu_beta[yr_sp_in],sigma_beta);
  
  lp = alpha[sp_yr] + beta[sp_yr] .* x + lambda[ii] .* age_obs + nu[ii] .* age_obs .* x;
  
  shp_over_lp_exp = shp[ii] ./ exp(lp);
  y ~ gamma(shp[ii], shp_over_lp_exp);
  
}
