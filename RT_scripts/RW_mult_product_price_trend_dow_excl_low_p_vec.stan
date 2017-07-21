functions {
  int neg_bi_safe_rng(real eta, real phi) {
    real phi_div_exp_eta = phi / exp(eta);
    real gam_draw;
    if (phi == 0 || phi_div_exp_eta == 0 || is_inf(phi) || is_inf(phi_div_exp_eta))
      return -9;
    gam_draw = gamma_rng(phi, phi_div_exp_eta);
    if (gam_draw > 1e9 || is_nan(gam_draw)) {
      return -9;
    } else {
      return poisson_rng(gam_draw);
    }
  }
  vector standardize(vector x) {
    return (x - mean(x)) ./ sd(x);
  }
  matrix standardize_mat(matrix x) {
    return (x - mean(to_vector(x))) ./ sd(to_vector(x));
  }
  int[] which(int[] x, int ret_size, int inp_size) {
    int which_vec[ret_size];
    int inds_cntr = 0;
    for (n in 1:inp_size)
      if (x[n] == 1) {
        inds_cntr = inds_cntr + 1;
        which_vec[inds_cntr] = n;
      }
    return which_vec;
  }
}
data {
  int<lower=1> N;
  int<lower=1> J;
  int y[J, N];
  matrix[J,N] price;
  int dow_ind[N];
  int N_hypo_prices;
  vector[N_hypo_prices] hypo_prices;
}
transformed data {
  matrix[N, J] price_std;
  vector[N] inds_std;
  int n_oos = 30;
  vector[n_oos] inds_for;
  vector[N_hypo_prices] hypo_prices_std;
  int wday_oos[n_oos];
  int last_wday = dow_ind[N] - 1;
  int include_indicator[J, N];
  int n_obs[J];

  for (j in 1:J) {
    for (n in 1:N)
      if ((price[j,n] >= 4.99 && n == 1) || (price[j,n] >= 4.99 && price[j,n - 1] >= 4.99)) {
        include_indicator[j,n] = 1;
      } else {
        include_indicator[j,n] = 0;
      }
    n_obs[j] = sum(include_indicator[j]);
  }

  price_std = standardize_mat(price)';
  {
    vector[N] inds;
    for (n in 1:N)
      inds[n] = n;
    {
      real incr_t = 1 / sd(inds); 
      inds_std = standardize(inds);
      for (i in 1:rows(inds_for))
        inds_for[i] = i * incr_t + inds_std[N];
    }
  }
  for (n in 1:N_hypo_prices) {
    hypo_prices_std[n] = (hypo_prices[n] - mean(price)) / sd(price);
  }
  for (n in 1:n_oos) 
    wday_oos[n] = (last_wday + n) % 7 + 1;
}
parameters {
  //
  vector[J] eta_phi;
  real<lower=0> sigma_phi;
  real mu_phi;
 
  real<lower=0> sigma_sigma;
  real mu_sigma;
  row_vector[J] eta_sigma;
  
  vector[J] beta_t;
  real alpha_t;
  real<lower=0> sigma_beta_t;
  matrix[N,J] eta;
  
  real alpha;
  real<lower=0> sigma_mu;
  vector[J] mu;
  
  vector[J] beta_p;
  real alpha_p_eta;
  real<lower=0> sigma_beta_p;

  vector[N] day_eta;
  real<lower=0> sigma_day;
  vector[7] dow;
}
transformed parameters {
  matrix[N,J] theta;
  row_vector[J] sigma;
  vector[J] phi;
  vector[N] day;
  real alpha_p;
  
  sigma = exp(mu_sigma + sigma_sigma * eta_sigma);
  phi = exp(mu_phi + sigma_phi * eta_phi);
  theta[1] = sigma .* eta[1];
  for (n in 2:N)
    theta[n] = theta[n-1] + sigma .* eta[n];
  day = alpha + sigma_day * day_eta + dow[dow_ind];
  alpha_p = -0.05 + 0.1 * alpha_p_eta;
}
model {
  to_vector(eta) ~ normal(0, 1);
  
  // random intercepts
  alpha ~ normal(2, 2);
  sigma_mu ~ normal(0, 1);
  mu ~ normal(0, 1);
  
  // random elasticities
  alpha_p_eta ~ normal(0, 1);
  sigma_beta_p ~ normal(0, 1);
  beta_p ~ normal(0, 1);
  
  // random time trends
  alpha_t ~ normal(-0.25, 1);
  sigma_beta_t ~ normal(0, 1);
  beta_t ~ normal(0, 1);
  
  // random observation noise
  eta_phi ~ normal(0, 1);
  sigma_phi ~ normal(0, 1);
  mu_phi ~ normal(0, 1);
  
  // random walk noise
  mu_sigma ~ normal(0, 1);
  sigma_sigma ~ normal(0, 1);
  eta_sigma ~ normal(0, 1);

  day_eta ~ normal(0, 1);
  sigma_day ~ normal(0, 1);
  dow ~ normal(0, 1);
  
  for (j in 1:J) {
    int which_vec[n_obs[j]] = which(include_indicator[j], n_obs[j], N);
    y[j,which_vec] ~ neg_binomial_2_log((day + sigma_mu*mu[j]
                                + (alpha_p + sigma_beta_p * beta_p[j]) * col(price_std,j)
                                + (alpha_t + sigma_beta_t * beta_t[j]) * inds_std
                                + col(theta,j))[which_vec], phi[j]);
  }
}
generated quantities {
  int y_rep[J, N];
  int y_forward[J, n_oos];
  matrix[J, n_oos] rev_forward;
  matrix[J, n_oos] theta_for;
  matrix[J, N_hypo_prices] y_hypo;
  matrix[J, N_hypo_prices] rev_hypo;
  vector[n_oos] new_day;

  for (n in 1:n_oos)
    new_day[n] = normal_rng(dow[wday_oos[n]], sigma_day);
  for (j in 1:J) {
    theta_for[j,1] = theta[N,j] + normal_rng(0, sigma[j]);
    y_forward[j,1] = neg_bi_safe_rng(alpha + new_day[1] + sigma_mu*mu[j] + (alpha_p + sigma_beta_p * beta_p[j]) * price_std[N,j] + (alpha_t + sigma_beta_t * beta_t[j]) * inds_for[1] + theta_for[j,1], phi[j]);
    rev_forward[j,1] = y_forward[j, 1] * price[j, N];
  }
  for (j in 1:J)
    for (i in 2:rows(inds_for)) {
      theta_for[j,i] = theta_for[j,i - 1] + normal_rng(0, sigma[j]);
      y_forward[j,i] = neg_bi_safe_rng(alpha + new_day[i] + sigma_mu*mu[j] + (alpha_p + sigma_beta_p * beta_p[j]) * price_std[N,j] + (alpha_t + sigma_beta_t * beta_t[j]) * inds_for[i] + theta_for[j,i], phi[j]);
      rev_forward[j,i] = y_forward[j, i] * price[j, N];
    }

  for (j in 1:J) {
    for (p in 1:N_hypo_prices) {
      int y_hypo_prod[n_oos];
      for (i in 1:n_oos) {
          y_hypo_prod[i] = neg_bi_safe_rng(alpha + new_day[i] + sigma_mu*mu[j] + (alpha_p + sigma_beta_p * beta_p[j]) * hypo_prices_std[p] + (alpha_t + sigma_beta_t * beta_t[j]) * inds_for[i] + theta_for[j,i], phi[j]);
      }
      y_hypo[j,p] = sum(y_hypo_prod);
      rev_hypo[j,p] = y_hypo[j,p] * hypo_prices[p];
    }
  }
}
