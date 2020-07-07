## Estimate Rt using epinow2
## Test if the assumed delay to case or death observation is three days too long and much more variable than the true delay

## Load dependencies and set parameters ------------------------------------------------------
rm(list = ls())
source('00-load_packages.R')
source('00-util.R')
source('00-run_test.R')
ggplot2::theme_set(theme_bw())

## Check if synthetic data already exists.
## If so, load from cache
## If not, make the data
if(!file.exists('true_pars.rds')){source('01-make_SEIR_data.R')}
parlist <- readRDS('true_pars.rds')
# Synthetic data is loaded using:
# get_sim_df()

## Set parameters for EpiNow2 test
testpars <- list(
  last_obs_time = 150,
  output_folder = 'misspec-delay-pars',
  ## True delays
  true_mean_case_delay = 5,
  true_sd_case_delay = 1.7,
  true_mean_death_delay = 15,
  true_sd_death_delay = 1.5,
  true_mean_inc = exp(EpiNow2::covid_incubation_period[1, ]$mean),
  true_sd_inc = exp(EpiNow2::covid_incubation_period[1, ]$sd))
## Delays specified in model
testpars$input_mean_case_delay = testpars$true_mean_case_delay + 3
testpars$input_sd_case_delay = testpars$true_sd_case_delay + 2
testpars$input_mean_death_delay = testpars$true_mean_death_delay
testpars$input_sd_death_delay = testpars$true_sd_death_delay
testpars$input_mean_inc = testpars$true_mean_inc
testpars$input_sd_inc = testpars$true_sd_inc
testpars$input_mean_gi = parlist$true_mean_GI
testpars$input_sd_gi = sqrt(parlist$true_var_GI)
dir_check(testpars$output_folder)

run_test(parlist, 
         testpars)

