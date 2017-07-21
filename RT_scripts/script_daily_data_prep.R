
# Daily model run
source('daily_model_data_prep.R')

daily_raw_data <- readRDS('/analytics/data/tests/for_eric/assigned_converted_netezza.us.randomhouse.com_tableau_apps_dev_fetch_data.sql.jinja_amz_midvol.json_625715384dc3d8de9acd53574a656812.csv.RDS')
p_keys <- unique(daily_raw_data$primary_prod_key)
daily_stan_data <- data_updater(p_keys, daily_raw_data)

saveRDS(daily_stan_data, file = '/analytics/data/tests/for_eric/daily_model_test_20170721.RDS')

rstan::stan_rdump(names(daily_stan_data), 
                  file = '/analytics/data/tests/for_eric/daily_model_test_20170721.dat',
                  envir = list2env(daily_stan_data))

# For model compilation, run the following command from inside the toplevel directory for cmdstan
# make rhome/rtrangucci/PricingEngine/RW_mult_product_price_trend_dow_excl_low_p_vec
