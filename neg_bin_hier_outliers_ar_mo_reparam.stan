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
  real prior_global_alpha_sd;
  real<lower=0> prior_global_alpha_sd_sd;
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
  real<lower=0,upper=1> month_ar_raw;
  vector[14] month_raw;
  real<lower=0> ar_sd;
}
transformed parameters {
  vector[J] alpha_prod;
  vector[J] beta_time; 
  vector[J] beta_price;
  vector[J] beta_ar;
  real month_ar;
  vector[14] month;
  vector[N] etas;     
  real alpha_sq;

  month_ar <- 2 * month_ar_raw - 1;
  alpha_prod <- sigma_alpha_prod * alpha_prod_raw;
  alpha_sq <- ar_sd / (1.0 - pow(month_ar,2.0));
  
  beta_time <- mu_beta_time + sigma_beta_time * beta_time_raw;
  beta_price <- mu_beta_price + sigma_beta_price * beta_price_raw;
  beta_ar <- mu_beta_ar + sigma_beta_ar * beta_ar_raw;
  {
    matrix[14,14] month_cov;
    matrix[14,14] L_month_cov;
    for (n in 1:13) {
      for (m in (n+1):14) {
        month_cov[n, m] <- pow(month_ar,m-n) * alpha_sq;
        month_cov[m, n] <- month_cov[n, m];
      }
      month_cov[n, n] <- alpha_sq;
    }
    month_cov[14, 14] <- alpha_sq;
    L_month_cov <- cholesky_decompose(month_cov);
    month <- L_month_cov * month_raw;
  }

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
  sigma_alpha_prod ~ normal(prior_global_alpha_sd, prior_global_alpha_sd_sd);
  sigma_beta_time ~ normal(0, 1);
  sigma_beta_price ~ normal(0, 1);
  sigma_beta_ar ~ normal(0, 1);
  month_ar_raw ~ beta(4, 2);
  ar_sd ~ normal(0, 1);

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
