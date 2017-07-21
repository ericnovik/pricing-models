library(ggplot2)
library(rstan)
library(dplyr)

# This file calculates the out-of-sample forecasts for each book
# under the observed prices and the daily, random walk model's parameters
# Allows for the model validation and criticism of the dail model.

source('daily_model_data_prep.R')
round_to_99 <- function(x) {
  (round(x*100,digits = -2) - 1) / 100
}

# Read in data used to fit the daily model
book_data <- readRDS('class_123_data_update_20170521_w_dow_out_mo_ind.RDS')

# Read in the file used to push price recommendations 
# only needed for the list of PROD_KEYs used to build a portfolio of books 
recos <- read.csv('hi_vol_recommended_prices_20170524.csv', sep = '|')

# Read in samples from daily model, likely created in the "summarise_stan_fit.R" file
s <- readRDS('daily_run_dow_model_fit.RDS')

# Read in the newest raw data (assign step in pipeline) for the books 
# of interest, likely to be all books in recos file
new_data <- readRDS('assigned_converted_netezza.us.randomhouse.com_tableau_apps_dev_fetch_data.sql.jinja_pull_hi_vol_full_alt_20170630.json_24975eade1786a3735c926f1165da639.csv.RDS')

# Read in the raw data (assign step in pipeline) used
# to fit the model
old_data <- readRDS('assigned_converted_netezza.us.randomhouse.com_tableau_apps_dev_fetch_data.sql.jinja_pull_7k_full_alt_20170521.json_3903d10718558ed0e9afab8a67fd4dd2.csv.RDS')
s_date <- min(old_data$date)

# Number of days forward the forecast should be run for
# Limited only by the lenght of your new data and the number
# of out-of-sample days calculated by your Stan model
days <- as.integer(max(new_data$date) - max(old_data$date))
# Change stan model
days <- ifelse(days > book_data$n_oos, book_data$n_oos, days)

# Create new data from book data
# Filter data before model data starts
new_data <- new_data %>% filter(date >= s_date)
book_dat_update <- data_updater(book_data$keys, new_data)

# Calculate observed daily revenue for each book
book_dat_update$rev <- book_dat_update$y * book_dat_update$price

# Get list of dates
dates <- unique(new_data$date)

# Get list of dates for out-of-sample period
oos_dates <- dates[(book_data$N + 1):(book_data$N + days)]

# Expose neg_bi_safe_rng from stan model for use in forecast
mod <- stan_model('~/Dropbox/PricingEngine/RW_mult_product_price_trend_dow_excl_low_p_vec.stan')
funcs <- expose_stan_functions(mod)

# Get out-of-sample price array for each book and standardize the 
# prices
price_arr <- book_dat_update$price[,(book_data$N + 1):(book_data$N + days)]
prices_std <- (price_arr - mean(book_data$price)) / sd(book_data$price)

# Get number of posterior samples
N_iter <- nrow(s$eta_phi)

# Prep array for revenue and quantity sold out-of-sample forecasts
obs_oos <- array(NA, c(book_data$J, days, N_iter, 2))

# Quantities needed for linear-time-trend extrapolation
incr_t <- 1 / sd(1:book_data$N)
mu_t <- mean(1:book_data$N)

# Quantities needed for day-of-week effect extrapolation
dow_oos <- book_dat_update$dow_ind[(book_data$N + 1):(book_data$N + days)]
dow_e_oos <- s$dow[,dow_oos] + sweep(matrix(rnorm(N_iter * days),N_iter, days),1,s$sigma_day,FUN = '*')

# Loop for generating quantity sold and revenue forecasts
# NB: This loop is tied to the structure of the DAILY model
# so if that structure changes, this loop needs to change
# Specifically, lines 81 and 86 and 87 need to change
for (book in 1:book_data$J) {
  rw <- matrix(NA, N_iter, days)
  for (ind in 1:N_iter) {
    rw[ind, ] <- s$theta[ind, book_data$N, book] + cumsum(rnorm(days, 0, s$sigma[ind, book])) +
      dow_e_oos[ind,]
  }

  for (ind in 1:N_iter) {
    mu_units <- s$alpha[ind] + s$sigma_mu[ind] * s$mu[ind, book] + (s$alpha_p[ind] + s$sigma_beta_p[ind] * s$beta_p[ind, book]) * prices_std[book,] +
      (s$alpha_t[ind] + s$sigma_beta_t[ind] * s$beta_t[ind, book]) * (incr_t * (book_data$N - mu_t) + incr_t * (1:days)) + rw[ind,]
    obs_oos[book, , ind, 1] <- sapply(mu_units, function(mu) neg_bi_safe_rng(mu, s$phi[ind,book]))
    obs_oos[book, , ind, 2] <- obs_oos[book, , ind, 1] * price_arr[book,]
  }
  print(paste0("book #", book, " with prod key = ", book_data$keys[book], " finished"))
}

# data.frame for out-of-sample quantity sold forecasts from cmdstan
# representing no change in price, matched to observed quantity sold
oos_fores_vol <- data.frame(med = as.vector(t(apply(s$y_forward,c(2,3),median))),
                            hi = as.vector(t(apply(s$y_forward,c(2,3),quantile,0.9))),
                            lo = as.vector(t(apply(s$y_forward,c(2,3),quantile,0.1))),
                            key = as.vector(sapply(book_data$keys, rep, days)),
                            date = rep(oos_dates,49),
                            price_change = rep(1:days,49),
                            obs = as.vector(t(book_dat_update$y[,(book_data$N + 1):(book_data$N + days)])))

# Plot product by product forecast vs. observed quantity sold
oos_fores_vol %>% filter(key %in% c(recos$PROD_KEY)) %>% ggplot() +
  geom_line(aes(x = date, y = med)) + geom_ribbon(aes(x = date, ymin = lo, ymax = hi), alpha = 0.3) + 
  geom_point(aes(x = date, y = obs)) + geom_vline(xintercept = as.numeric(as.Date('2017-05-29'))) + 
  facet_wrap( ~ key, scales = 'free')

# data.frame for out-of-sample book by day revenue forecasts from cmdstan
# representing no change in price, matched to observed revenue
oos_fores_rev <- data.frame(med = as.vector(t(apply(s$rev_forward,c(2,3),median))),
                            hi = as.vector(t(apply(s$rev_forward,c(2,3),quantile,0.9))),
                            lo = as.vector(t(apply(s$rev_forward,c(2,3),quantile,0.1))),
                            key = as.vector(sapply(book_data$keys, rep, days)),
                            date = rep(oos_dates,49),
                            price_change = rep(1:days,49),
                            obs = as.vector(t(book_dat_update$rev[,(book_data$N + 1):(book_data$N + days)])))

# data.frame for out-of-sample portfolio revenue forecast
# representing no change in price, matched to observed revenue
oos_fores_rev_tot <- data.frame(med = median(apply(s$rev_forward[,,1:days],1,sum)),
                                hi = quantile(apply(s$rev_forward[,,1:days],1,sum),0.95),
                                lo = quantile(apply(s$rev_forward[,,1:days],1,sum),0.1),
                                obs = sum(book_dat_update$rev[,(book_data$N + 1):(book_data$N + days)]))

# Plot product by product forecast vs. observed revenue
oos_fores_rev %>% filter(key %in% c(recos$PROD_KEY)) %>% ggplot() +
  geom_line(aes(x = date, y = med)) + geom_ribbon(aes(x = date, ymin = lo, ymax = hi), alpha = 0.3) + 
  geom_point(aes(x = date, y = obs)) + geom_vline(xintercept = as.numeric(as.Date('2017-05-29'))) + 
  facet_wrap( ~ key, scales = 'free')

# Pick out revenue forecasts from obs_oos array
rev_forward <- aperm(obs_oos[,,,2],c(3,1,2))

# data.frame for out-of-sample book by book total revenue forecast for
# out-of-sample period with observed price, matched to total observed revenue by book
oos_sum_rev <- data.frame(med = apply(apply(rev_forward[,,15:days],c(1,2),sum),c(2),median),
                          hi = apply(apply(rev_forward[,,15:days],c(1,2),sum),c(2),quantile,0.975),
                          lo = apply(apply(rev_forward[,,15:days],c(1,2),sum),c(2),quantile,0.025),
                          key = as.character(book_data$keys),
                          obs = apply(book_dat_update$rev[,(book_data$N + 15):(book_data$N + days)],1,sum))

# observed revenue
obs_rev <- apply(book_dat_update$rev[,(book_data$N + 1):(book_data$N + days)],1,sum)

# Daily observed revenue by book
rev_tot <- apply(rev_forward[,,1:days],c(1,2),sum)
plt_df <- data.frame(rev = as.vector(rev_tot), key = as.vector(sapply(book_data$keys, rep, N_iter)),
                     obs_rev = as.vector(sapply(obs_rev, rep, N_iter)))

# Book by book plot, excluding outliers 110144361 and 30780604
plt_df %>% filter(key %in% c(as.character(recos$PROD_KEY)),
                        key != as.character(110144361),
                        key != as.character(30780604)) %>% ggplot() +
  geom_histogram(aes(x = rev)) +
  geom_vline(aes(xintercept = obs_rev)) +
  facet_wrap(~ key, scales = 'free')


# Book by book plot, excluding outliers 110144361 and 30780604
oos_sum_rev %>% filter(key %in% c(as.character(recos$PROD_KEY)),
                        key != as.character(110144361),
                        key != as.character(30780604)) %>% arrange(lo) %>% ggplot() +
  geom_errorbar(aes(x = key, ymin = lo, ymax = hi)) + geom_point(aes(x = key, y = obs)) +
  theme(axis.text.x = element_text(angle = 90, hjust = 1))


# Generate a list of PROD_KEYS which did have prices
# changed according to optimal price recommendation
# Excluding outliers
oos_sum_rev %>% filter(key %in% c(as.character(recos$PROD_KEY)),
                       !(key %in% as.character(c(110144361,30780604,30726475,110154444)))) %>% .$key %>% as.character %>% as.integer -> keys_sel

# Match stan index with PROD_KEYs in keys_sel
key_sel_vec <- match(keys_sel, book_data$keys)

# Match stan index with PROD_KEYS not in keys_sel
# For use as control titles
titles_comp <- match(setdiff(oos_sum_rev$key, union(recos$PROD_KEY,c(110144361,30780604,30726475,110154444))), book_data$keys)

# Total portfolio sum of revenue for out-of-sample period using
# observed prices
oos_sum_rev_tot <- data.frame(med = median(apply(rev_forward[,key_sel_vec,5:days],c(1),sum)),
                              hi = quantile(apply(rev_forward[,key_sel_vec,5:days],c(1),sum), 0.95),
                              lo = quantile(apply(rev_forward[,key_sel_vec,5:days],c(1),sum), 0.05),
                              obs = sum(book_dat_update$rev[key_sel_vec,(book_data$N + 5):(book_data$N + days)]))

# Total portfolio sum of revenue for out-of-sample period
# under no pricing change
oos_fores_rev_tot <- data.frame(med = median(apply(s$rev_forward[,key_sel_vec,5:days],1,sum)),
                                hi = quantile(apply(s$rev_forward[,key_sel_vec,5:days],1,sum),0.95),
                                lo = quantile(apply(s$rev_forward[,key_sel_vec,5:days],1,sum),0.05),
                                obs = sum(book_dat_update$rev[key_sel_vec,(book_data$N + 5):(book_data$N + days)]))

# Calculate distribution of change in revenue due to price change
change_in_rev_dist <- apply(rev_forward[,key_sel_vec,5:days],c(1),sum) / apply(s$rev_forward[,key_sel_vec,5:days],1,sum)

# Calculate book by book change in revenue at observed prices vs. 
# scenario where no prices were changed
oos_change_rev_tot <- data.frame(med = median(),
                                hi = quantile(apply(rev_forward[,key_sel_vec,5:days],c(1),sum) / apply(s$rev_forward[,key_sel_vec,5:days],1,sum),0.95),
                                lo = quantile(apply(rev_forward[,key_sel_vec,5:days],c(1),sum) / apply(s$rev_forward[,key_sel_vec,5:days],1,sum),0.05))

# Obesrved lift for titles_comp and books set to optimal prices
obs_lift <- sum(book_dat_update$rev[key_sel_vec,(book_data$N + 5):(book_data$N + days)]) / sum(book_dat_update$rev[key_sel_vec,(book_data$N - 27):(book_data$N + 4)]) - 1
obs_lift <- sum(book_dat_update$rev[titles_comp,(book_data$N + 5):(book_data$N + days)]) / sum(book_dat_update$rev[titles_comp,(book_data$N - 27):(book_data$N + 4)]) - 1

(mean(book_dat_update$rev[key_sel_vec,(book_data$N + 5):(book_data$N + 34)]) / mean(book_dat_update$rev[key_sel_vec,(book_data$N - 25):(book_data$N + 4)]) - 1) - (mean(book_dat_update$rev[titles_comp,(book_data$N + 5):(book_data$N + 34)]) / mean(book_dat_update$rev[titles_comp,(book_data$N - 25):(book_data$N + 4)]) - 1)
