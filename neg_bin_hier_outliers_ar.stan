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
  int year_month_ind[N];
  int num_outlier;
  vector[N] promote;
  vector[N] L1_qty;
}
transformed data {
  vector[N] price_sqr;
  vector[N] time_sqr;
  vector[N] L1_log_qty;
  
  for (i in 1:N) {
    price_sqr[i] <- price[i]^2;
    time_sqr[i] <- time[i]^2 / 10;
    L1_log_qty[i] <- if_else(L1_qty[i] == -999, 0, if_else(L1_qty[i] == 0, log(0.5),
                             log(L1_qty[i])));
  }
}
parameters {
  real alpha;
  vector[J] alpha_prod_raw;
  real<lower=0> sigma_alpha_prod;
  real<lower=0> sigma_beta_price;
  real<lower=0> sigma_beta_time;
  real<lower=0> sigma_beta_ar;
  real mu_beta_price;
  real mu_beta_time;
  real mu_beta_ar;
  real promote_price;
  real promote_level;
  vector[J] beta_time_raw; 
  vector[J] beta_price_raw;
  vector[J] beta_ar_raw;
  real<lower=0> phi;  
  vector<lower=0>[num_outlier] outlier_effects;
  real<lower=0> month_sd;
  vector[14] month_raw;
}
transformed parameters {
  vector[J] alpha_prod;
  vector[J] beta_time; 
  vector[J] beta_price;
  vector[J] beta_ar;
  vector[14] month;
  vector[N] etas;     

  alpha_prod <- sigma_alpha_prod * alpha_prod_raw;
  month <- month_sd * month_raw;
  beta_time <- mu_beta_time + sigma_beta_time * beta_time_raw;
  beta_price <- mu_beta_price + sigma_beta_price * beta_price_raw;
  beta_ar <- mu_beta_ar + sigma_beta_ar * beta_ar_raw;

  {
    int outlier_counter;
    real price_effect;
    // linear predictor
    outlier_counter <- 1;
    for (n in 1:N) {
      if (promote[n] == 0) {
        price_effect <- price[n] / 10 * beta_price[prod_key[n]];
      } else {
        price_effect <- price[n] / 10 * promote_price + promote_level;
      }
      if (outlier[n] > 0) {
        etas[n] <- alpha + alpha_prod[prod_key[n]] 
          + price_effect
          + 2 * time[n] * beta_time[prod_key[n]]
          + month[year_month_ind[n]]
          + outlier_effects[outlier_counter]
          + L1_log_qty[n] * beta_ar[prod_key[n]];
        outlier_counter <- outlier_counter + 1;
      } else {
        etas[n] <- alpha + alpha_prod[prod_key[n]]
          + price_effect
          + 2 * time[n] * beta_time[prod_key[n]]
          + month[year_month_ind[n]]
          + L1_log_qty[n] * beta_ar[prod_key[n]];
      }
    }
  }
}
model {
  // sigma priors
  sigma_alpha_prod ~ normal(0, 1);
  sigma_beta_time ~ normal(0, 1);
  sigma_beta_price ~ normal(0, 1);
  sigma_beta_ar ~ normal(0, 1);
  month_sd ~ normal(0, 1);

  mu_beta_price ~ normal(-1, 1);
  mu_beta_time ~ normal(-1, 1);
  mu_beta_ar ~ normal(0, 1);
  month_raw ~ normal(0, 1);

  alpha ~ normal(prior_global_mean, 
                 prior_global_sd);
  
  alpha_prod_raw ~ normal(0, 1);
  beta_time_raw ~ normal(0, 1);
  beta_price_raw ~ normal(0, 1);
  beta_ar_raw ~ normal(0, 1);

  // promotional priors
  promote_price ~ normal(-1, 1);
  promote_level ~ normal(1, 1);

  // overdispersion prior spec
  phi ~ cauchy(0, 2.5);

  // outlier effect
  outlier_effects ~ normal(0, 2.5);
  
  // likelihood
  qty ~ neg_binomial_2_log(etas, phi);
}
generated quantities {
  real phis[N];
  
  for (n in 1:N) 
    phis[n] <- phi;
}
