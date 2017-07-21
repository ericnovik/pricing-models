library(rstan)
library(dplyr)
library(lubridate)

filter_missing_data <- function(df) {
  n_rows_pass = length(unique(df$date))
  passing_prods = df %>% group_by(prod_key) %>% 
    tally() %>%
    filter(n == n_rows_pass)
  df_clean = df %>% semi_join(passing_prods)
}

#' @param keys vector of PROD_KEYS to model
#' @param full_data data from pipeline "assign"step 
#'                  containing a superset of the 
#'                  PROD_KEYS from the first argument
#' @param n_oos number of days to forecast revenue and quantity
#'              sold under different prices and no price change
#' @param hypo_prices vector of prices at which to forecast revenue
#'                    and quantity sold
#' @return list to be fed to rstan::stan_rdump for cmdstan modeling
data_updater <- function(keys, full_data, n_oos = 30,
                         hypo_prices = seq(from = 4.99, 
                                           to = 14.99, by = 1)) {
  date_set <- unique(full_data$date)
  N_new <- length(date_set)
  J <- length(keys)
  qtys <- matrix(NA_real_,J, N_new)
  prices <- matrix(NA_real_,J, N_new)
  full_data %>% group_by(prod_key) %>% arrange(date) -> full_data
  dow_ind <- wday(date_set)
  
  for (key_i in seq_along(keys)) {
    key <- keys[key_i]
    df <- full_data %>% filter(prod_key == key)
    if(length(df$units) > N_new) {
      df %>% group_by(date) %>% summarise(units = first(units), 
                                          price = first(price)) -> df
      print(key)
    }
    else if (length(df$units) < N_new) {
      print(key)
    }
    qty <- df$units
    price <- df$price
    qty[is.na(qty)] <- 0
    qty[qty < 0] <- 0
    qtys[key_i,] <- qty
    prices[key_i,] <- price
  }
  mo_fac <- apply(cbind(lubridate::year(date_set),lubridate::month(date_set)),1,function(x) paste(x,collapse='_'))
  wk_fac <- apply(cbind(lubridate::isoyear(date_set),lubridate::isoweek(date_set)),1,function(x) paste(x,collapse='_'))
  mo_lvls <- unique(mo_fac)
  wk_lvls <- unique(wk_fac)
  mo_ind <- as.integer(factor(mo_fac, levels = mo_lvls))
  wk_ind <- as.integer(factor(wk_fac, levels = wk_lvls))
  return(list(N = N_new, J = J, y = qtys, 
              keys = keys,
              price = prices, dow_ind = dow_ind,
              N_hypo_prices = length(hypo_prices),
              hypo_prices = hypo_prices,
              N_wk = length(wk_lvls),
              n_oos = n_oos))
}

plotter <- function(dat) {
  bad_keys <- c()
  for (j in 1:dat$J) {
    plot(dat$y[j,])
    bad_key <- readline(prompt = 'Is this a bad key: ')
    if (bad_key == '1')
      bad_keys <- c(bad_keys, dat$keys[j])
  }
  return(setdiff(dat$keys,bad_keys))
}
