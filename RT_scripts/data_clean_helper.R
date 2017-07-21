library(dplyr)
library(rstan)
library(lubridate)
library(loo)
library(reshape2)
library(parallel)

# Transform daily data from transform step
# in pipeline to weekly data
# for use in the weekly model
#' @param df data from transform step
#' @return tibble  of weekly aggregates
prep_raw_daily_data <- function(df) {
  ## temporarily skip this block
  # go the weekly
  df <- df %>% 
    mutate(
      qty_sld1 = ifelse(is.na(qty_sld),0,qty_sld),
      qty_sld2 = ifelse(qty_sld1 >= 30000, prior_units, qty_sld1) # extreme outliers
    ) %>% distinct(prod_key, date, .keep_all=TRUE)
  n_days_prods <- df %>%
    group_by(prod_key) %>%
    summarise(n_days = length(date))
  max_days <- max(n_days_prods$n_days)
  excl_prods <- n_days_prods %>% filter(n_days < max_days) %>% .$prod_key
  incl_prods <- setdiff(n_days_prods$prod_key,excl_prods)
 
  first_day <- min(df$date)
  first_weekday <- lubridate::wday(first_day)
  first_day <- ifelse(first_weekday == 2, first_day, first_day + (7 - (first_weekday - 2)))
  last_day <- max(df$date)
  last_weekday <- lubridate::wday(last_day)
  last_day <- ifelse(last_weekday == 1, last_day, last_day - (last_weekday - 1))
  df <- df %>% filter(date >= first_day & date <= last_day, prod_key %in% incl_prods)
  
  dw <- df %>% 
    group_by(prod_key, calendar_week_year, calendar_week) %>%
    summarise(
      map = first(map),
      qty_sld = sum(qty_sld2,na.rm=T),
      years_since_osd = median(years_since_osd),
      years_since_lppe = median(years_since_lppe),
      author = first(author_short),
      low_p_wk = as.integer(min(price) < 4.99),
      price = ifelse(sum(qty_sld2) > 0,weighted.mean(price, qty_sld2),mean(price)),
      calendar_month = round(median(calendar_month)),
      repeat_bestseller = first(repeat_bestseller),
      week_max_date = week(max(date)),
      year_max_date = year(max(date)),
      n_days = n(),
      wday_str = paste('7',week_max_date,year_max_date,sep='-'),
      #week_of_sunday = lubridate::as_date(wday_str,format='%u-%W-%Y')
      week_of_sunday = as.Date(wday_str,format='%u-%W-%Y')
    )
  return(dw)
}

# Used to filter data prior to creating Stan data
# according to set of critera applied at the book level
#' @param data weekly aggregated data 
#' @param qty_filt condition on which to exclude books below a certain weekly
#'                 median quantity threshold
#' @param price_filt filters for max and min price observed
#' @param time_filt filter books that have gone on sale no earlier than one year ago
#'                  at any point during the data time frame
filt_data <- function(data, qty_filt = 'med_qty > 20', price_filt =  'max_p <= 19.99 & min_p >= 4.99', time_filt = 'min_ysd > 1') {
  data %>% group_by(prod_key) %>% summarise(med_qty = median(qty_sld),
                                          max_p = max(price),
                                          min_p = min(price),
                                          min_ysd = min(years_since_osd)) %>% filter_(paste(c(qty_filt, price_filt, time_filt), collapse = ' & ')) -> prods_keep
  
  data %>% filter(prod_key %in% prods_keep$prod_key) -> sub_fit
  return(sub_fit)
}

# Helper function for padding zeros
#' @param x integer
padder <- function(x) {
  nchar_x <- nchar(x)
  if (nchar_x < 2)
    return(paste0('0',x))
  else
    return(paste0(x))
}

# Take results of prep_raw_daily_data
# and append stan related indexing variables
# for weekly and monthly random effects
#' @param df tibble from prep_raw_daily_data
#' @return tibble ready to generate stan input
prep_stan_dat <- function(df) {
  df <- df %>% ungroup()
  year_mo_df <- data.frame(week_of_sunday = sort(unique(df$week_of_sunday))) %>%
    mutate(mo_nm = sapply(month(week_of_sunday),padder),
           year = year(week_of_sunday),
           year_mo = paste(year,mo_nm, sep='_'),
           fac_year_mo = factor(year_mo),
           mo_cnt = as.integer(fac_year_mo),
           year_ind = as.integer(factor(isoyear(week_of_sunday)))) %>%
    select(week_of_sunday,mo_cnt,year_ind)
  df %>% left_join(year_mo_df,by='week_of_sunday') -> df
  
  dw_T <- df %>% ungroup() %>%
    mutate(stan_prod = as.integer(factor(prod_key)),
           stan_subj = as.integer(factor(map))) %>%
    arrange(prod_key, calendar_week_year, calendar_week) %>% group_by(prod_key) %>%
    mutate(
      L1_qty = c(-999, head(qty_sld, -1)),
      z2 = (qty_sld - median(qty_sld, na.rm=T)) / (mad(qty_sld,na.rm=T)),
      sales_spike = ifelse(ifelse(is.na(z2),0,z2) > 3, 1, 0),
      L1_spike = c(-999, head(sales_spike, -1)),
      L1_low_p = c(-999, head(low_p_wk, -1))
    ) %>% ungroup() %>%
    filter(L1_qty != -999, low_p_wk == 0, L1_low_p == 0) %>% 
     mutate(
      price_sqr = price ^ 2,
      wk_ind = as.integer(factor(mapply(function(x, y) paste0(x,'_',y), calendar_week_year, calendar_week)))
     ) %>% 
    group_by(prod_key) %>%
    mutate(
      log1p_L1_qty = log1p(L1_qty),
      log1p_qty = log1p(qty_sld)
    )
  return(dw_T)
}

# Function taking tibble from prep_stan_dat and generating a list to
# pass to rstan::stan_rdump
#' @param df tibble from prep_stan_dat
#' @param scale_X boolean controlling whether columns of X are scaled to sd = 1
#' @param y_trans function to transform quantity sold (outcome) from integer to some other space
#' @param coefs character vector indicating which variables to include as fixed effects
#' @param eta_prior prior for lkj_corr priors over group-level correlation matrices
#' @return list with stan_data element to pass to stan_rdump, also contains
#'          map from stan book index to prod_key and stan subject index to subject
mk_stan_data_two_groups <- function(df, scale_X = T, y_trans = function(x) log1p(x), 
                                    coefs = c('price','years_since_osd','price_squared'), eta_prior = 5) {
  X = df[,coefs]
  subj_map <- unique(df[,c('stan_prod','stan_subj')]) %>% arrange(stan_prod) %>% select(stan_subj)
  means = colMeans(X)
  X_cent = sweep(X,2,means)
  scales = apply(X_cent,2,sd)
  if (scale_X) {
    X_cent = sweep(X_cent,2,scales,FUN = '/')
  }
  last_obs = df %>% group_by(prod_key) %>% summarise(l_p = last(price),
                                                     obs = last(qty_sld))
  prod_inds <- df %>% group_by(prod_key) %>% group_indices()
  
  QR = qr(X_cent)
  Q = qr.Q(QR)
  R = qr.R(QR)
  sqrtN = sqrt(nrow(X) - 1)
  Q_scale = Q * sqrtN
  R_scale = R * 1 / sqrtN
  R_inv = solve(R_scale)

  price_hypo = seq(4.99,14.99,by=1)
  
  max_date <- lubridate::ymd(max(df$week_of_sunday))
  next_wk <- max_date + 7
  new_mo_ind <- as.integer(lubridate::month(max_date) != lubridate::month(next_wk))
  stan_data <- with(df, list(N = length(qty_sld),
                             D = dim(X_cent)[2],
                             scales = scales,
                             new_mo_fore = new_mo_ind,
                             J_g_1 = length(unique(stan_prod)),
                             J_g_2 = length(unique(stan_subj)),
                             wk_ind = wk_ind,
                             N_mo = length(unique(mo_cnt)),
                             N_wk = max(wk_ind),
                             prices_raw = df$price,
                             mo_cnt = mo_cnt,
                             y=sapply(qty_sld,y_trans),
                             inds_1=stan_prod,
                             inds_2=stan_subj,
                             inds_j_2=subj_map$stan_subj,
                             mo_ind = calendar_month,
                             Q = cbind(Q_scale,1),
                             des_mat = cbind(X_cent,1),
                             hypo_prices = price_hypo,
                             N_hypo = length(price_hypo),
                             means = rep(0,length(means)),
                             means_fore = means,
                             R_inv = R_inv))
  prod_map <- unique(df[,c('prod_key','stan_prod')])
  subj_map <- unique(df[,c('map','stan_subj')])
  return(list(stan_data=stan_data,
              prod_map = prod_map,
              subj_map = subj_map,
              means = means,
              R_inv = R_inv))
}

# Wrapper for filter and stan data creation
#' @param df result of prep_raw_daily_data
#' @param qty_filter book-level quantity filter
#' @param y_trans transformation for outcome (quantity sold by week)
#' @param coefs coefs to include as fixed effects
#' @param price_filter book-level price filter
#' @return list to save to disk containing all relevant info prior to running cmdstan
full_data_pipe <- function(df, qty_filter, ytrans, coefs, price_filter) {
  stan_data_prepped <- filt_data(df, qty_filt = qty_filter, price_filt = price_filter) %>%
     prep_stan_dat()
  stan_data <- stan_data_prepped %>%
    mk_stan_data_two_groups(coefs = coefs,scale_X = T,y_trans = ytrans,eta_prior = 1)
  stan_data$raw_data <- stan_data_prepped
  return(stan_data)
}

# Wrapper for saving model data for input into cmdstan
#' @param mdl_data data returned by full_data_pipe
#' @param mdl_sv_pth full path to directory to save the model data
#' @param mdl_prefix name to prepend to .RDS and .dat files for use in stan model
sv_mdl_data <- function(mdl_data, mdl_sv_pth, mdl_prefix) {
  saveRDS(mdl_data,paste0(mdl_sv_pth,mdl_prefix,'.RDS'))
  rstan::stan_rdump(names(mdl_data$stan_data),file = paste0(mdl_sv_pth,mdl_prefix,'.dat'),envir = list2env(mdl_data$stan_data))
}

# Wrapper for running data prep pipeline
#' @param mdl_data tibble returned by prep_raw_daily_data
#' @param mdl_sv_pth full path to directory to save the model data
#' @param mdl_prefix name to prepend to .RDS and .dat files for use in stan model
#' @param qty_filter book-level quantity filter
#' @param y_trans transformation for outcome (quantity sold by week)
#' @param coefs coefs to include as fixed effects
#' @param price_filter book-level price filter
data_prep_wrapper <- function(mdl_data, mdl_sv_pth, mdl_prefix, qty_filt, coefs = c("log1p_L1_qty","wk_ind","price","sales_spike" ), ytrans = I, price_filter = 'max_p <= 19.99 & min_p >= 0.01') {
  mdl_data <- mdl_data %>%
    mutate(
     years_since_lppe = days_since_lppe_osd / 365.25
  ) %>% prep_raw_daily_data()
  mdl_data %>% group_by(prod_key) %>% summarise(n_nas = sum(is.na(years_since_lppe))) %>% filter(n_nas < 10) %>% .$prod_key -> nas_incl ## revisit recent release logic
  jnd_cln <- mdl_data %>% ungroup() %>% filter(prod_key %in% nas_incl)
  jnd_cln %>% ungroup() %>% group_by(prod_key) %>% summarise(sd_lppe = min(diff(years_since_lppe),na.rm = T)) %>% filter(sd_lppe < 0) -> lppe_release_products

  jnd_cln <- filter(jnd_cln, !is.na(map),!(prod_key %in% lppe_release_products))
  stan_dat <- full_data_pipe(df = jnd_cln, qty_filter = qty_filt, ytrans = ytrans,coefs = coefs, price_filter = price_filter)
  sv_mdl_data(mdl_data = stan_dat, mdl_sv_pth = mdl_sv_pth, mdl_pre = mdl_prefix)
}
