
# weekly model prep

source('data_clean_helper.R')

weekly_model_raw_data <- readRDS('/analytics/data/tests/for_eric/transformed_sample4.RDS')
weekly_model_raw_data$map <- weekly_model_raw_data$rhsubject_desc

data_prep_wrapper(mdl_data = weekly_model_raw_data, 
                  mdl_sv_pth = '/analytics/data/tests/for_eric/weekly_model_run_20170721/',
                  mdl_prefix = 'wkly_mod_20170721', qty_filt = 'med_qty > 0')


# For model compilation, run the following command from inside the toplevel directory for cmdstan
# make rhome/rtrangucci/PricingEngine/QX_glob_var_prior_NB_AR_indy_book_re_final
