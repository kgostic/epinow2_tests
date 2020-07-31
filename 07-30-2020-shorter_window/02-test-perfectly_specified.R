## Estimate Rt using epinow2
## Test if the assumed delay to case or death observation is three days too long and much more variable than the true delay

## Load dependencies and set parameters ------------------------------------------------------
rm(list = ls())
source('../00-load_packages.R')
source('../00-util.R')
source('../00-run_test.R')
ggplot2::theme_set(theme_bw())
cat(getwd())
this_dir = gsub('.+/epinow2_tests/(.+)', '\\1', getwd())

## Check if synthetic data already exists.
## If so, load from cache
## If not, make the data
parlist <- load_parlist()
# Synthetic data is loaded using:
# get_sim_df()

## Set parameters for EpiNow2 test
testpars <- list(
  last_obs_time = 150,
  output_folder = 'perfectly_specified',
    # Cases
  true_log_mean_case_delay = log(5),
  true_log_sd_case_delay = log(1.7),
    # Deaths
  true_log_mean_death_delay = log(15),
  true_log_sd_death_delay = log(1.5),
    # Incubation
  true_log_mean_inc = EpiNow2::covid_incubation_period[1, ]$mean,
  true_log_sd_inc = EpiNow2::covid_incubation_period[1, ]$sd,
    # Generation interval
  true_mean_gi = parlist$true_mean_GI,
  true_sd_gi = sqrt(parlist$true_var_GI))

## Delays specified as priors, may or may not match truth
    # Cases
testpars$input_lmean_case_delay = testpars$true_log_mean_case_delay
testpars$input_lsd_case_delay = testpars$true_log_sd_case_delay
    # Deaths
testpars$input_lmean_death_delay = testpars$true_log_mean_death_delay
testpars$input_lsd_death_delay = testpars$true_log_sd_death_delay
    # Incubation
testpars$input_lmean_inc = testpars$true_log_mean_inc
testpars$input_lsd_inc = testpars$true_log_sd_inc
    # Generation int
testpars$input_mean_gi = testpars$true_mean_gi
testpars$input_sd_gi = testpars$true_sd_gi

## If output directory does not exist, creat it.
dir_check(testpars$output_folder)

## Run epinow() ------------------------------ 
run_test(parlist, testpars, prior_smoothing_window = 3)
## Code called by run_test
print(run_test)




## Genearte a report using outputs of test run
rmarkdown::render(input = "../Make_test_report.Rmd", output_file = sprintf('%s-results.html', testpars$output_folder), output_dir = '../reports/',
                  params = list(path = "perfectly_specified", replot = TRUE, dir_date = this_dir))
