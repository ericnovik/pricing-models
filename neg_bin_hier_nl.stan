data {
  int<lower=0> N;       // total number of observations  
  int<lower=0> qty[N];  // response variable
  int<lower=1> J;       // number of titles
  int<lower=1,upper=J> prod_key[N]; // specific title indicator: product key
  vector[N] time;
  vector[N] price;
}
transformed data {
  vector[N] price_sqr;
  
  for (i in 1:N) {
    price_sqr[i] <- price[i]^2;
  }
}
parameters {
  real alpha;
  vector[J] alpha_prod_raw;
  real<lower=0> sigma_alpha_prod;
  real<lower=0> sigma_beta_price;
  real<lower=0> sigma_beta_time;
  real mu_beta_price;
  real mu_beta_time;
//  real beta_price;
//  real beta_time;
  vector[J] beta_time_raw; 
  vector[J] beta_price_raw;
  vector<lower=0>[J] phi_raw;  
  real<lower=0> mu_phi;
  real<lower=0> sigma_phi;
}
transformed parameters {
  vector[J] alpha_prod;
  vector[J] beta_time; 
  vector[J] beta_price;
  vector[J] phi;

  alpha_prod <- sigma_alpha_prod * alpha_prod_raw;
  beta_time <- mu_beta_time + sigma_beta_time * beta_time_raw;
  beta_price <- mu_beta_price + sigma_beta_price * beta_price_raw;
  for (j in 1:J)
    phi[j] <- exp(phi_raw[j]);
}
model {
  vector[N] eta;     
  vector[N] phis;
  // linear predictor
  for (n in 1:N) {
    eta[n] <- alpha + alpha_prod[prod_key[n]] + exp(2 * time[n] * beta_time[prod_key[n]]) + 
              price[n] / 10 * beta_price[prod_key[n]];
    phis[n] <- phi[prod_key[n]];
  }
  
  // sigma priors
  sigma_alpha_prod ~ normal(0, 1);
  sigma_beta_time ~ normal(0, 1);
  sigma_beta_price ~ normal(0, 1);

  mu_beta_price ~ normal(0, 1);
  mu_beta_time ~ normal(0, 1);
//  beta_price ~ normal(0, 1);
//  beta_time ~ normal(0, 1);

  alpha ~ normal(7, 1);
  
  alpha_prod_raw ~ normal(0, 1);
  beta_time_raw ~ normal(0, 1);
  beta_price_raw ~ normal(0, 1);

  // overdispersion prior spec
  mu_phi ~ cauchy(0, 2.5);
  phi_raw ~ normal(mu_phi, sigma_phi);
  sigma_phi ~ cauchy(0, 1);
  
  // likelihood
  qty ~ neg_binomial_2_log(eta, phis);
}
generated quantities {
  real etas[N];
  real phis[N];
  
  for (n in 1:N) {
    etas[n] <- alpha + alpha_prod[prod_key[n]] + exp(2 * time[n] * beta_time[prod_key[n]]) + 
             price[n] / 10  * beta_price[prod_key[n]];
    phis[n] <- phi[prod_key[n]];
  }
}
