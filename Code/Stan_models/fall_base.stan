// Fall Ap/SS
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
  int<lower = -1, upper = 1> age_obs[N];                // for missing cov
  int<lower = 0, upper = 1> age_miss[N];       // for missing cov
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

  real<lower=0> mu_gamma;
  vector<lower=0>[Nsp] gamma;

  real omega;
  real psi;
  real eta;

  real mu_nu;
  real mu_lambda;

  real<lower=0> sigma_alpha;
  real<lower=0> sigma_beta;
  real<lower=0> sigma_mu_beta;
  real<lower=0> sigma_gamma;
  real<lower=0> sigma_nu;
  real<lower=0> sigma_lambda;

  real<lower = 0, upper = 1> kappa_age;               // for missing cov

  // For non-centered params.
  vector[Nsp] lambda_raw;
  vector[Nsp] nu_raw;

  vector[Nsp] mu_beta_raw;
}

transformed parameters{

  vector[Nsp] lambda;
  vector[Nsp] nu;

  vector[Nsp] mu_beta;

  vector[Nsp] delta;

  lambda = mu_lambda + sigma_lambda * lambda_raw;
  nu = mu_nu + sigma_nu * nu_raw;

  delta = omega + psi * migr_len + eta * breed_lat;

  mu_beta = delta + sigma_mu_beta * mu_beta_raw;
}

model {

  vector[N] lp;
  vector[N] lp0;
  vector[N] lp1;
  vector[N] shp_over_lp_exp;
  vector[N] shp_over_lp_exp0;
  vector[N] shp_over_lp_exp1;

  sigma_beta ~ std_normal();
  sigma_alpha ~ std_normal();
  sigma_mu_beta ~ std_normal();
  sigma_gamma ~ std_normal();
  sigma_nu ~ std_normal();
  sigma_lambda ~ std_normal();

  mu_shp ~ gamma(1,1);
  rate_shp ~ gamma(1,1);

  mu_gamma ~ std_normal();

  gamma ~ normal(mu_gamma,sigma_gamma);

  omega ~ std_normal();
  psi ~ std_normal();
  eta ~ std_normal();

  kappa_age ~ uniform(0,1);

  mu_lambda ~ std_normal();
  mu_nu ~ std_normal();

  shp ~ gamma(mu_shp,rate_shp);

  //For non-centered parameters
  lambda_raw ~ std_normal();
  nu_raw ~ std_normal();
  mu_beta_raw ~ std_normal();

  alpha ~ normal(gamma[yr_sp_in],sigma_alpha);
  beta ~ normal(mu_beta[yr_sp_in],sigma_beta);

  //Below code works to impute the missing age parameters
  for (j in 1:N){
      if (age_miss[j] == 1) {
            // age missing
            // as if draw was a 1
            lp1[j] = alpha[sp_yr[j]] + beta[sp_yr[j]] * x[j] + lambda[ii[j]] + nu[ii[j]] * x[j];
            shp_over_lp_exp1[j] = shp[ii[j]] ./ exp(lp1[j]);

            // as if draw was a 0
            lp0[j] = alpha[sp_yr[j]] + beta[sp_yr[j]] * x[j];
            shp_over_lp_exp0[j] = shp[ii[j]] ./ exp(lp0[j]);

            target += log_mix(kappa_age,
                    gamma_lpdf(y[j] | shp[ii[j]], shp_over_lp_exp1[j]),
                    gamma_lpdf(y[j] | shp[ii[j]], shp_over_lp_exp0[j]));
        } else {
            // age not missing
            age_obs[j] ~ bernoulli(kappa_age);

            lp[j] = alpha[sp_yr[j]] + beta[sp_yr[j]] * x[j] + lambda[ii[j]] * age_obs[j] + nu[ii[j]] * age_obs[j] * x[j];
            shp_over_lp_exp[j] = shp[ii[j]] ./ exp(lp[j]);
            y[j] ~ gamma(shp[ii[j]], shp_over_lp_exp[j]);

        }
  }
}
