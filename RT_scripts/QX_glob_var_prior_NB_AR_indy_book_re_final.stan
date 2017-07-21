functions {
  int[] count_obs(int[] inds) {
    int N = size(inds);
    int J = max(inds);
    int counter = 1;
    int cur_ind = inds[1];
    int n_obs[J];
    for (n in 2:N) {
      if (inds[n] == inds[n - 1]) {
        counter = counter + 1;
      } else {
        n_obs[cur_ind] = counter;
        counter = 1;
        cur_ind = inds[n];
      }
    }
    n_obs[inds[N]] = counter;
    return n_obs;
  }
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
  real pred(vector beta, real ar_1,
            real wk_ind, real price, 
            real spike) {
              return beta[1] * ar_1
                     + beta[2] * wk_ind
                     + beta[3] * price
                     + beta[4] * spike
                     + beta[5];
            }
}
data {
  int N; // number of data points
  int D; // number of predictors
  int J_g_1; // Number of groups
  int J_g_2; // Number of groups
  int N_wk; // number of weeks observed
  int N_mo; // number of months
  int N_hypo; // number of months
  int wk_ind[N];
  int mo_cnt[N];
  int new_mo_fore;
  vector[N] prices_raw;
  matrix[N, D + 1] Q; // scaled Q matrix, with last column of 1's
  matrix[N, D + 1] des_mat; // raw X matrix, with last column of 1's
  int inds_1[N]; // vector of integer group indices for each obs
  int inds_2[N]; // vector of integer group indices for each obs
  int inds_j_2[J_g_1]; // vector of integer group 2 indices for each group 1
  int y[N]; // vector of observations
  matrix[D, D] R_inv; // Inverse of scaled R from QR decomp
  vector[D] means_fore;
  vector[D] means;
  vector[D] scales;
  vector[N_hypo] hypo_prices;
}
transformed data {
  matrix[D + 1, N] t_X;
  matrix[D + 1, N] t_Q;
  vector[D] prior_mean;
  vector[2 * D + 5] prior_counts;
  vector[J_g_1] last_obs;
  real last_wk_ind;
  real last_sales_spike;
  vector[N_hypo] hypo_prices_std;
  int n_obs[J_g_1] = count_obs(inds_1);
  vector[J_g_1] last_price;
  vector[J_g_1] last_price_raw;
  int pcounts_rows;
  int num_sds;
  real mean_log1p_qty_sld;
  t_Q = Q';
  t_X = des_mat';
  {
    real csum;
    real qty;
    csum = 0;
    for (n in 1:N) {
      qty = y[n];
      csum = csum + log1p(qty);
    }
    mean_log1p_qty_sld = csum / N;
  }
  for (j in 1:J_g_1) {
    last_obs[j] = (log1p(y[sum(n_obs[1:j])]) - means_fore[1]) / scales[1];
    last_price[j] = des_mat[sum(n_obs[1:j]), 3];
    last_price_raw[j] = prices_raw[sum(n_obs[1:j])];
  }
  last_wk_ind = (max(wk_ind) + 1 - means_fore[2]) / scales[2];
  last_sales_spike = -means_fore[4] / scales[4];
  for (i in 1:N_hypo) {
    hypo_prices_std[i] = (hypo_prices[i] - means_fore[3]) / scales[3];
  }
  prior_mean = rep_vector(0, D);
  num_sds = 2 * D + 2;
  for (d in 1:(num_sds/2 - 1)) 
    prior_counts[d] = 2;
  prior_counts[num_sds/2] = 20;
  for (d in (num_sds/2 + 1):(num_sds - 1))
    prior_counts[d] = 2;
  prior_counts[num_sds] = 4;
  prior_counts[num_sds + 1] = 4;
  prior_counts[num_sds + 2] = 10;
  prior_counts[num_sds + 3] = 2;
  pcounts_rows = rows(prior_counts);
}
parameters {
  row_vector[D + 1] fixed;
  matrix[D + 1, J_g_1] eta_1;
  matrix[D + 1, J_g_2] eta_2;
  real<lower=0> scale;
  simplex[pcounts_rows] prop_var;
  real<lower=0> sd_obs;
  vector[N_mo] mo_eta;
  vector[N_wk] wk_eta;
  cholesky_factor_corr[D + 1] L_omega_2;
  cholesky_factor_corr[D + 1] L_omega_1;
}
transformed parameters {
  matrix[D + 1, J_g_1] group_beta_1;
  matrix[D + 1, J_g_2] group_beta_2;
  vector[D + 1] sds_1;
  vector[D + 1] sds_2;
  vector[pcounts_rows] vars;
  vector[pcounts_rows] sds;
  real mo_sd;
	real wk_sd;
  row_vector[N] prod_mu;
  vector[N_mo] mo;
  vector[N_wk] wk;
  vars = rows(prior_counts) * square(scale) * prop_var;
  for (d in 1:(pcounts_rows)) 
    sds[d] = sqrt(vars[d]);
  sds_1[1:(D + 1)] = sds[1:(D + 1)];
  sds_2 = sds[(D + 2):num_sds];
  mo_sd = sds[num_sds + 1];
  wk_sd = sds[num_sds + 2];

  mo = mo_sd * mo_eta;
  wk = wk_sd * wk_eta;

  group_beta_1 = diag_pre_multiply(sds_1, L_omega_1) * eta_1;
  group_beta_2 = diag_pre_multiply(sds_2, L_omega_2) * eta_2;

  for (n in 1:N) {
    prod_mu[n] = col(group_beta_1,inds_1[n])' * col(t_X,n)
                 + col(group_beta_2,inds_2[n])' * col(t_X,n)
                 + mo[mo_cnt[n]]
                 + wk[wk_ind[n]];
    
  }
  prod_mu = prod_mu + fixed * t_Q;
}
model {
  fixed[D + 1] ~ normal(mean_log1p_qty_sld, 1);
  fixed[1:D] ~ normal(prior_mean, 1);
  scale ~ gamma(2, 8);
  prop_var ~ dirichlet(prior_counts);
  L_omega_1 ~ lkj_corr_cholesky(3.0);
  L_omega_2 ~ lkj_corr_cholesky(3.0);
  to_vector(eta_1) ~ normal(0, 1);
  to_vector(eta_2) ~ normal(0, 1);
  sd_obs ~ gamma(2, 2);
  mo_eta ~ normal(0, 1);
  wk_eta ~ normal(0, 1);

  y ~ neg_binomial_2_log(prod_mu, 1/sd_obs);
}
generated quantities {
  vector[D+1] fixed_r;
  int y_forward[J_g_1];
  vector[J_g_1] rev_forward;
  vector[J_g_1] mu_forward;
  int y_hypo[J_g_1, N_hypo];
  matrix[J_g_1, N_hypo] rev_hypo;
  matrix[J_g_1, N_hypo] mu_hypo;
  matrix[D + 1, J_g_1] prod_beta;
  matrix[D + 1, D + 1] Omega_1;
  matrix[D + 1, D + 1] Omega_2;
  real port_rev_forward;
  int port_forward;
  real new_wk; 
  real new_mo;
  fixed_r[1:D] = R_inv * fixed[1:D]';
  fixed_r[D+1] = fixed[D+1] - fixed_r[1:D]' * means;

  new_wk = normal_rng(0, wk_sd);
  new_mo = normal_rng(0, mo_sd);
  for (j in 1:J_g_1) {
    prod_beta[1:D,j] = fixed_r[1:D] + group_beta_1[1:D,j] + group_beta_2[1:D,inds_j_2[j]];
    prod_beta[D+1, j] = fixed_r[D + 1] + group_beta_1[D+1,j] + group_beta_2[D+1,inds_j_2[j]];
    mu_forward[j] = pred(col(prod_beta,j),
                         last_obs[j],
                         last_wk_ind,
                         last_price[j],
                         last_sales_spike)
                    + new_wk
                    + (new_mo_fore == 1 ? new_mo : mo[N_mo]);
    y_forward[j] = neg_bi_safe_rng(mu_forward[j], 1 / sd_obs);
    rev_forward[j] = y_forward[j] * last_price_raw[j]; 
    for (p in 1:N_hypo) {
      mu_hypo[j,p] = pred(col(prod_beta,j),
                          last_obs[j],
                          last_wk_ind,
                          hypo_prices_std[p],
                          last_sales_spike)
                      + new_wk
                      + (new_mo_fore == 1 ? new_mo : mo[N_mo]);
      y_hypo[j, p] = neg_bi_safe_rng(mu_hypo[j,p], 1 / sd_obs);
      rev_hypo[j, p] = y_hypo[j, p] * hypo_prices[p];
    }
  }
  port_rev_forward = sum(rev_forward);
  port_forward = sum(y_forward);
  Omega_1 = L_omega_1 * L_omega_1';
  Omega_2 = L_omega_2 * L_omega_2';
}
