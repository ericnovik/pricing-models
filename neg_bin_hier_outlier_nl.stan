data {
  int<lower=0> N;       // total number of observations  
  int<lower=0> qty[N];  // response variable
  int<lower=1> J;       // number of titles
  int<lower=1,upper=J> prod_key[N]; // specific title indicator: product key
  vector[N] time;
  vector[N] price;
  vector[N] outlier;
  real prior_global_mean;
  real<lower=0> prior_global_sd;
}
transformed data {
  vector[N] price_sqr;
  vector[N] time_sqr;
  
  for (i in 1:N) {
    price_sqr[i] <- price[i]^2;
    time_sqr[i] <- time[i]^2 / 10;
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
  real mu_beta_price_sqr;
  vector[J] beta_time_raw; 
  vector[J] beta_price_raw;
  real<lower=0> phi;  
  real<lower=0> outlier_effect;
}
transformed parameters {
  vector[J] alpha_prod;
  vector[J] beta_time; 
  vector[J] beta_price;

  alpha_prod <- sigma_alpha_prod * alpha_prod_raw;
  beta_time <- mu_beta_time + sigma_beta_time * beta_time_raw;
  beta_price <- mu_beta_price + sigma_beta_price * beta_price_raw;
}
model {
  vector[N] eta;     
  // linear predictor
  for (n in 1:N) {
    eta[n] <- alpha + alpha_prod[prod_key[n]] 
      + exp(2 * time[n] * beta_time[prod_key[n]])
      + price[n] / 10 * beta_price[prod_key[n]]
      + price_sqr[n] * mu_beta_price_sqr
      + outlier[n] * outlier_effect;
  }
  
  // sigma priors
  sigma_alpha_prod ~ normal(0, 1);
  sigma_beta_time ~ normal(0, 1);
  sigma_beta_price ~ normal(0, 1);

  mu_beta_price ~ normal(-1, 1);
  mu_beta_time ~ normal(-1, 1);
  mu_beta_price_sqr ~ normal(0, 1);

  alpha ~ normal(prior_global_mean, 
                 prior_global_sd);
  
  alpha_prod_raw ~ normal(0, 1);
  beta_time_raw ~ normal(0, 1);
  beta_price_raw ~ normal(0, 1);

  // overdispersion prior spec
  phi ~ cauchy(0, 2.5);

  // outlier effect
  outlier_effect ~ cauchy(0, 1);
  
  // likelihood
  qty ~ neg_binomial_2_log(eta, phi);
}
generated quantities {
  real etas[N];
  real phis[N];
  
  for (n in 1:N) {
    etas[n] <- alpha + alpha_prod[prod_key[n]]
      + exp(2 * time[n] * beta_time[prod_key[n]])
      + price[n] / 10  * beta_price[prod_key[n]]
      + price_sqr[n] * mu_beta_price_sqr
      + outlier[n] * outlier_effect;
    phis[n] <- phi;
  }
}
