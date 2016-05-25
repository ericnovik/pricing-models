setwd('~/RH')
library(ggplot2)
library(rstan)
library(dplyr)
library(reshape2)
source('hir_data_sim.R')

dat_myst <- readRDS('data_2016-04-08 12:24:30_low_vol_Mystery.rds')
dat_myst <- arrange(dat_myst,prod_key,year,week)


plot_ysd_price <- function(prod_key, dat) {
  return(plot(dat[dat$prod_key == prod_key,'ysd'],
              dat[dat$prod_key == prod_key,'price']))
}

plot_price_qty <- function(prod_key, dat) {
  return(plot(dat[dat$prod_key == prod_key,'price'],
              log(dat[dat$prod_key == prod_key,'w_qty_sld'])))
}

plot_ysd_qty <- function(prod_key=NULL, dat, plog=TRUE,...) {
  if (!is.null(prod_key)) {
    dat_plt <- filter(dat,prod_key == prod_key)
  } else {
    dat_plt <- dat
  }
  qty <- dat_plt$w_qty_sld
  ysd <- dat_plt$ysd
  if (plog) {
    qty <- ifelse(qty == 0, log(0.5), log(qty))
  }
  return(plot(as.data.frame(ysd)[,1],
              as.data.frame(qty)[,1], ...))
}

plot_qty_acf <- function(prod_key=NULL, dat, plog=TRUE,...) {
  if (!is.null(prod_key)) {
    dat_plt <- filter(dat,prod_key == prod_key)
  } else {
    dat_plt <- dat
  }
  qty <- dat_plt$w_qty_sld
  ysd <- dat_plt$ysd
  if (plog) {
    qty <- ifelse(qty == 0, log(0.5), log(qty))
  }
  return(pacf(as.data.frame(qty)[,1], ...))
}

dat_sub <- filter(group_by(dat_myst,prod_key), w_qty_sld[1] > 20)
dat_sub$prod_key_stan <- as.integer(factor(dat_sub$prod_key))
dat_sub <- arrange(dat_sub,prod_key_stan,year,week)
dat_sub$y_obs_num <- unlist(sapply(1:length(unique(dat_sub$prod_key_stan)),
                                   function(x) 1:sum(dat_sub$prod_key_stan == x)))
dat_sub <- dat_sub[dat_sub$y_obs_num != 1,]
data_real <- with(dat_sub,
                  list(N = length(prod_key),
                       qty = w_qty_sld,
                       time = ysd,
                       price = price,
                       J = length(unique(prod_key)),
                       prod_key = prod_key_stan))

#test_fit <- stan('neg_bin_rd_gp_var_nohier.stan', data = data_real, chains = 2, cores = 2, iter = 400, warmup = 300)
#fit_hier <- stan('neg_bin_hier.stan', data = data_real, chains = 2, cores = 2, iter = 400, warmup = 300)
#fit_hier_nl <- stan('neg_bin_hier_nl.stan', data = data_real, chains = 2, cores = 2, iter = 400, warmup = 300)
fit_hier_nl_no_first <- stan('neg_bin_hier_nl.stan', data = data_real, chains = 2, cores = 2, iter = 2000, warmup = 1500,control=list(adapt_delta=0.9))
fit_hier_ <- stan('poisson_hier_nl.stan', data = data_real, chains = 2, cores = 1, iter = 2000, warmup = 1500,control=list(adapt_delta=0.9, max_treedepth=13))

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
  post_loc <- 
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

pp_means_by_group <- function(y_rep, data) {
  gpd_y <- group_by(data.frame(y = data$qty,
                               prod_key = data$prod_key,
                               y_obs_num = unlist(sapply(1:max(data$prod_key),
                                              function(x) 1:sum(data$prod_key == x)))),prod_key)
  means <- summarise(gpd_y, mean(y))
  means <- dplyr::rename(means, obs_mean=`mean(y)`)
  gpd_y_rep <- data.frame(cbind(y_rep,data$prod_key,gpd_y$y_obs_num))
  names(gpd_y_rep) <- c(paste('y_rep_',1:dim(y_rep)[2],sep=''),'prod_key','y_obs_num')
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

pp <- pp_fit(fit_hier_gp_no_first)
id_y <- pp_means_by_group(pp, data_real)

plt_mean_v_data(id_y$jnd_indiv)
plt_dens_means(id_y$jnd)
