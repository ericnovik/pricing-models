setwd('~/RH')
library(ggplot2)
library(rstan)
library(dplyr)
library(reshape2)
library(TTR)
source('hir_data_sim.R')

dat_myst <- readRDS('data_2016-04-08 12:24:30_low_vol_Mystery.rds')
dat_myst <- arrange(dat_myst,prod_key,year,week)
dat_myst <- dat_myst %>% mutate(year_month = interaction(factor(year),factor(month),drop=T),
                                promote = price <= 3.99)
dat_myst$year_month_ind <- dat_myst %>% group_indices(year_month)
dat_myst$dsd <- dat_myst$ysd * 365
dat_myst$dsd_mod <- round(dat_myst$dsd,0)
dat_myst[dat_myst$year == 2016,]$dsd_mod <- dat_myst[dat_myst$year == 2016,]$dsd_mod - 1L
dat_myst$wsd <- floor(dat_myst$dsd_mod / 7)

gpd_dat <- group_by(dat_myst,prod_key)
small_timers <- filter(gpd_dat, mean(w_qty_sld) <= 20 & !any(wsd < 100))
big_timers <- filter(gpd_dat,  mean(w_qty_sld) > 20 & !any(wsd < 100))

ggplot(data = small_timers,aes(x=wsd,y=w_qty_sld,colour = prod_key_factor )) + geom_point() + geom_line() 
ggplot(data = big_timers,aes(x=wsd,y=w_qty_sld,colour = prod_key_factor )) + geom_point() + geom_line() 

data_clean <- function(df) {
  gpd <- group_by(df, prod_key)
  gpd$prod_key_stan <- as.integer(factor(gpd$prod_key))
  df_clean <- gpd %>% arrange(prod_key,year,week) %>% 
    mutate(run_median = TTR::runMedian(w_qty_sld, cumulative = TRUE,n=1),
           run_mad = TTR::runMAD(w_qty_sld,cumulative=TRUE,n=1),
           run_z = (w_qty_sld - run_median) / (1.4826 * run_mad),
           outlier = ifelse(ifelse(is.na(run_z),0,run_z) > 3, 1, 0),
           L1_qty = c(-999,head(w_qty_sld,-1)))
  df_clean <- arrange(df_clean,prod_key_stan,year,week)
  df_clean$y_obs_num <- as.vector(unlist(sapply(1:length(unique(df_clean$prod_key_stan)),
                                     function(x) 1:sum(df_clean$prod_key_stan == x))))
  stan_data <- with(df_clean,
                    list(N = length(prod_key),
                         qty = w_qty_sld,
                         time = ysd,
                         price = price,
                         J = length(unique(prod_key)),
                         prod_key = prod_key_stan,
                         outlier = outlier,
                         year_month_ind = year_month_ind,
                         num_outlier = sum(outlier),
                         promote = promote,
                         L1_qty = L1_qty))
  return(list(stan_data = stan_data,
              df_clean = df_clean))
}

small_clean <- data_clean(small_timers)
small_stan <- small_clean$stan_data
small_stan$prior_global_mean = 2
small_stan$prior_global_sd = 1
small_stan$prior_global_alpha_sd = 0
small_stan$prior_global_alpha_sd_sd = 1

big_clean <- data_clean(big_timers)
big_stan <- big_clean$stan_data
big_stan$prior_global_mean = 5
big_stan$prior_global_sd = 1
big_stan$prior_global_alpha_sd = 2
big_stan$prior_global_alpha_sd_sd = 1
#test_fit <- stan('neg_bin_rd_gp_var_nohier.stan', data = data_real, chains = 2, cores = 2, iter = 400, warmup = 300)
#fit_hier <- stan('neg_bin_hier.stan', data = data_real, chains = 2, cores = 2, iter = 400, warmup = 300)
#fit_hier_nl <- stan('neg_bin_hier_nl.stan', data = data_real, chains = 2, cores = 2, iter = 400, warmup = 300)
#fit_hier_nl_no_first <- stan('neg_bin_hier_outlier_nl.stan', data = data_real, chains = 2, cores = 2, iter = 2000, warmup = 1500,control=list(adapt_delta=0.9))
fit_hier_nl_no_first <- stan('neg_bin_hier_outlier_nl.stan', data = big_stan, 
                             chains = 2, cores = 2, iter = 500, warmup = 300)

mnth_mod <- stan('neg_bin_hier_outlier_nl_gp_var.stan', data = big_stan, 
                             chains = 2, cores = 2, iter = 500, warmup = 300)

mnth_mod <- stan('neg_bin_hier_outlier_nl_gp_var.stan', data = big_stan, 
                             chains = 2, cores = 2, iter = 500, warmup = 300)

mnth_mod_small <- stan('neg_bin_hier_outliers_ar.stan', data = small_stan, 
                             chains = 4, cores = 2, iter = 1000, warmup = 700)

mnth_mod_ar_big <- stan('neg_bin_hier_outliers_ar.stan', data = big_stan, 
                             chains = 4, cores = 2, iter = 1000, warmup = 700)

mnth_mod_ar_big_ar <- stan('neg_bin_hier_outliers_ar_mo.stan', data = big_stan, 
                             chains = 4, cores = 2, iter = 1000, warmup = 700)

fit_hier_nl_no_first <- stan('neg_bin_hier_outlier_nl.stan', data = big_stan, 
                             chains = 2, cores = 2, iter = 2000, warmup = 1500)
pp_gen <- function(loc, scale=NULL) {
  stopifnot(all(dim(loc) == dim(scale)))
  n_gen <- dim(loc)[1]
  n_y <- dim(loc)[2]
  
  if (is.null(scale)) {
    y_rep <- sapply(1:n_gen,function(n) rpois(n_y, lambda=exp(loc[n,])))
  } else {
    y_rep <- sapply(1:n_gen,function(n) rnegbin(n_y, mu=exp(loc[n,]),
                                                theta = scale[n,]))
  }
  return(y_rep)
}

post_loc_prec <- function(stan_fit) {
  loc_draws <- rstan::extract(stan_fit,pars=c('etas'))[[1]]
  prec_draws <- rstan::extract(stan_fit,pars=c('phis'))[[1]]
  return(list(loc=loc_draws,
              prec=prec_draws))
}

pp_fit <- function(stan_fit) {
  post_loc <- post_loc_prec(stan_fit)
  reps <- pp_gen(post_loc$loc,
                 post_loc$prec)
  return(reps)
}

pp_fit_pois <- function(stan_fit) {
  post_loc <- rstan::extract(stan_fit,pars=c('etas'))[[1]]
  reps <- pp_gen(post_loc)
  return(reps)
}

y_rep_gp <- function(stan_fit, data) {
  y_rep <- pp_fit(stan_fit)
  jn_yrep <- data.frame(cbind(y_rep,data$prod_key))
  names(jn_yrep) <- c(paste('y_rep_',1:200,sep=''),'prod_key')
  return(jn_yrep)
}

plt_fun <- function(reps, gpd, key, ...) {
  print(plot(apply(reps[reps$prod_key == key,1:200],1,mean), type='l', ...))
  print(points(1:length(gpd[gpd$prod_key == key,]$y), gpd[gpd$prod_key == key,]$y))
}

log_fun <- function(x) {
  x_filt <- ifelse(x == 0, log(0.5), log(x))
  return(x_filt)
}

#
# @param y,yrep,group Validated y, yrep, and group objects from the user.
# @param stat Either NULL or a string naming a function.
# @value If \code{stat} is NULL, a molten data frame grouped by group and
#   variable. If \code{stat} specifies a function then a summary table created
#   by dplyr::summarise.
#

pp_means_by_gp_by_x <- function(y, y_rep, gp, x) {
  lower_25 <- apply(y_rep,2,quantile,0.05)
  median <- apply(y_rep,2,quantile,0.5)
  upper_75 <- apply(y_rep,2,quantile,0.95)
  df <- data.frame(
    lower_25 = lower_25,
    y = median,
    upper_75 = upper_75,
    gp = gp,
    x = x)
  df_data <- data.frame(y = y,
                       x = x,
                       gp = gp)
  p <- ggplot(data = df, aes(x = x, y = y, ymin = lower_25,
                             ymax = upper_75)) + 
    geom_smooth(stat = "identity",
                fill = "#DCBCBC",
                color = "#C79999") + 
    
    geom_point(data = df_data, color = 'black', size = 0.25) +
    facet_wrap(facets = "gp", scales = "free_y")
  print(p)
}

                    

pp_means_by_group <- function(y_rep, data) {
  gpd_y <- group_by(data.frame(y = data$qty,
                               prod_key = data$prod_key,
                               y_obs_num = as.vector(unlist(sapply(1:max(data$prod_key),
                                              function(x) 1:sum(data$prod_key == x))))),prod_key)
  gpd_y_rep <- data.frame(cbind(y_rep,data$prod_key,gpd_y$y_obs_num))
  names(gpd_y_rep) <- c(paste('y_rep_',1:dim(y_rep)[2],sep=''),'prod_key','y_obs_num')
  gpd_y_rep <- filter(gpd_y_rep, y_obs_num >= 40)
  means <- summarise(gpd_y_rep, mean(y))
  means <- dplyr::rename(means, obs_mean=`mean(y)`)
  melted <- melt(gpd_y_rep, id.vars = c('prod_key','y_obs_num'))
  gpd_key_rep <- group_by(melted, prod_key, variable)
  gpd_key_obs <- group_by(melted, prod_key, y_obs_num)
  mean_dist_by_group <- summarise(gpd_key_rep,mean(value))
  jnd <- left_join(means, mean_dist_by_group, by = 'prod_key')
  jnd <- dplyr::rename(jnd, value = `mean(value)`)
  obs_mean_by_group <- summarise(gpd_key_obs,mean(value))
  jnd_indiv <- left_join(gpd_y, obs_mean_by_group, by = c('prod_key','y_obs_num'))
  jnd_indiv <- dplyr::rename(jnd_indiv, value = `mean(value)`)
  return(list(jnd=jnd, 
              jnd_indiv= jnd_indiv))
}

plt_dens_means <- function(df_plt) {
  plt <- ggplot(data = df_plt, aes(x = value)) + 
       geom_density(fill = "orange", alpha = 0.5) + # Make it pretty
       facet_wrap(~ prod_key, scales = "free") +
       geom_vline(aes(xintercept = obs_mean), colour = "red") +
       ggtitle("Actual observed means vs. modeled means\n") + theme_bw()
  return(plt)
}

plt_mean_v_data <- function(df_plt) {
  plt <- ggplot(data = df_plt, aes(x = y_obs_num, y = value)) +
       geom_line(aes(colour = 'red')) + # Make it pretty
       facet_wrap(~ prod_key, scales = "free") +
       geom_line(aes(x = y_obs_num, y = y),size = 0.1) + geom_point(aes(x = y_obs_num, y = y),size = 0.1) + 
       ggtitle("Modeled mean by week vs. observed data by week\n") + theme_bw()
  return(plt)
}

pp <- pp_fit(fit_hier_nl_no_first)
pp <- t(pp_fit(mnth_mod_small))
pp <- t(pp_fit(mnth_mod_ar_big))
pp_small <- t(pp_fit(mnth_mod_small))
pp_means_by_gp_by_x(y = small_stan$qty, gp = small_stan$prod_key, y_rep = pp_small, x = rep(1:60,118))
pp_means_by_gp_by_x(y =big_stan$qty, gp = big_stan$prod_key, y_rep = pp, x = rep(1:60,30))
pp <- pp_fit_pois(fit_hier_pois)
id_y <- pp_means_by_group(pp, big_stan)
p

plt_mean_v_data(id_y$jnd_indiv)
plt_dens_means(id_y$jnd)
