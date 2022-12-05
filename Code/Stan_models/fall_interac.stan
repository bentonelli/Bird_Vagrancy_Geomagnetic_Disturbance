//Fall interaction model
data {
  int<lower=0> N;
  vector[N] x_ap;
  vector[N] x_ss;
  vector[N] y;
  int<lower=1> Nyrs;
  int<lower = 1, upper = Nyrs> tt[N];
  int<lower=1> Nsp;
  int<lower = 1, upper = Nsp> ii[N];
  int<lower=1> N_sp_yr;
  int<lower = 1, upper = N_sp_yr> sp_yr[N];
  int<lower = 1, upper = Nsp> yr_sp_in[N_sp_yr];
}

parameters {
  
  real<lower=0> mu_shp;
  real<lower=0> rate_shp;
  vector<lower=0>[Nsp] shp;
  
  real mu_beta;
  real mu_theta;
  real mu_omega;
  real gamma;
  
  vector[N_sp_yr] alpha;
  vector[Nsp] beta;
  vector[Nsp] theta;
  vector[Nsp] omega;
  vector[Nsp] mu_alpha;

  
  real<lower=0> sigma_alpha;
  real<lower=0> sigma_beta;
  real<lower=0> sigma_theta;
  real<lower=0> sigma_omega;
  real<lower=0> sigma_mu_alpha;


}

model {
  
  vector[N] lp;
  vector[N] shp_over_lp_exp;
  
  mu_shp ~ normal(1,1);
  rate_shp ~ normal(1,1);
  gamma ~ normal(3,1);
  
  sigma_alpha ~ std_normal();
  sigma_beta ~ std_normal();
  sigma_theta ~ std_normal();
  sigma_omega ~ std_normal();

  mu_alpha ~ normal(gamma,sigma_mu_alpha);
  mu_beta ~ std_normal();
  mu_theta ~ std_normal();
  mu_omega ~ std_normal();
  
  shp ~ gamma(mu_shp,rate_shp);
  alpha ~ normal(mu_alpha[yr_sp_in],sigma_alpha);
  beta ~ normal(mu_beta,sigma_beta);
  theta ~ normal(mu_theta,sigma_theta);
  omega ~ normal(mu_omega,sigma_omega);
  
  lp = alpha[sp_yr] + beta[ii] .* x_ap + theta[ii] .* x_ss + omega[ii] .* x_ap .* x_ss; 
  shp_over_lp_exp = shp[ii] ./ exp(lp);
  y ~ gamma(shp[ii], shp_over_lp_exp);
  
}
