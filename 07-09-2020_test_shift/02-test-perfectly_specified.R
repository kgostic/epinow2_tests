## Estimate Rt using epinow2
## Test if the assumed delay to case or death observation is three days too long and much more variable than the true delay

## Load dependencies and set parameters ------------------------------------------------------
rm(list = ls())
source('../00-load_packages.R')
source('../00-util.R')
source('../00-run_test.R')
ggplot2::theme_set(theme_bw())

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
  ## True delays are used to generate data
  true_log_mean_case_delay = log(5),
  true_log_sd_case_delay = log(1.7),
  true_log_mean_death_delay = log(15),
  true_log_sd_death_delay = log(1.5),
  true_log_mean_inc = EpiNow2::covid_incubation_period[1, ]$mean,
  true_log_sd_inc = EpiNow2::covid_incubation_period[1, ]$sd,
  true_mean_gi = parlist$true_mean_GI,
  true_sd_gi = sqrt(parlist$true_var_GI))

## Delays specified as priors, may or may not match truth
testpars$input_lmean_case_delay = testpars$true_log_mean_case_delay
testpars$input_lsd_case_delay = testpars$true_log_sd_case_delay
testpars$input_lmean_death_delay = testpars$true_log_mean_death_delay
testpars$input_lsd_death_delay = testpars$true_log_sd_death_delay
testpars$input_lmean_inc = testpars$true_log_mean_inc
testpars$input_lsd_inc = testpars$true_log_sd_inc
testpars$input_mean_gi = testpars$true_mean_gi
testpars$input_sd_gi = testpars$true_sd_gi
## If output directory does not exist, creat it.
dir_check(testpars$output_folder)



## Set up inputs for estimate_infections ---------------------------------------

## Set delay distributions for input into EpiNow2 
generation_time <- list(mean = testpars$input_mean_gi, # gamma
                        mean_sd = 1,
                        sd = testpars$input_sd_gi,
                        sd_sd = 1,
                        max = 30)

incubation_period <- list(mean = testpars$input_lmean_inc, # lognormal
                          mean_sd = EpiNow2::covid_incubation_period[1, ]$mean_sd,
                          sd = testpars$input_lsd_inc,
                          sd_sd = EpiNow2::covid_incubation_period[1, ]$sd_sd,
                          max = 30)

obs_rep_delay <- list(mean = testpars$input_lmean_case_delay, ## lognormal, assume mean of 5 days for outpatient testing 
                      mean_sd = 1,   ## (48h from symptom onset + 48h for reporting + report the following day)
                      sd = testpars$input_lsd_case_delay,
                      sd_sd = 1,
                      max = 30)

death_rep_delay <- list(mean = testpars$input_lmean_death_delay, ## lognorma 
                        mean_sd = 1, #For death, assume a mean of 15 days from symptom onset and an sd of about 7d
                        sd = testpars$input_lsd_death_delay,
                        sd_sd = 1,
                        max = 30)


# Plot the specified distributions
png(sprintf('figs/specified_distributions-%s.png', testpars$output_folder), width = 7, height = 7, units = 'in', res = 300)
par(mfrow = c(2,2))
xx = seq(0, 30, by = 0.01)

true_case_dist <- function(x){dlnorm(x, testpars$true_log_mean_case_delay, testpars$true_log_sd_case_delay)}
true_death_dist <- function(x){dlnorm(x, testpars$true_log_mean_death_delay, testpars$true_log_sd_death_delay)} # Additional delay from symptoms -> observation

## Delays to case observation -----
# Input
plot(xx, dlnorm(xx, obs_rep_delay$mean, obs_rep_delay$sd), 
     type = 'l', main = 'assumed delay from onset to case detection', xlab = 'days', ylab = 'dens')
# Truth
lines(xx, dlnorm(xx, true_case_dist(xx)), col = 'red', lty = 2)
legend('topright', c('input', 'true'), col = c(1, 'red'), lty = 1)

## Delays to death ------
plot(xx, dlnorm(xx, death_rep_delay$mean, death_rep_delay$sd),
     type = 'l', main = 'assumed delay from onset to death', xlab = 'days', ylab = 'dens')
lines(xx, true_death_dist(xx), col = 'red', lty = 2)
legend('topright', c('input', 'true'), col = c(1, 'red'), lty = 1)

## Generation interval ---
plot(xx, 
     dgamma(xx, shape = with(generation_time, get_shape(mean, sd^2)), 
            rate = with(generation_time, get_rate(mean, sd^2))), 
     type = 'l', main = 'assumed generation time', xlab = 'days', ylab = 'dens')
lines(xx, dgamma(xx, get_shape(testpars$true_mean_gi, testpars$true_sd_gi^2), 
                 get_rate(testpars$true_mean_gi, testpars$true_sd_gi^2)), 
      col = 'red', lty = 2)
legend('topright', c('input', 'true'), col = c(1, 'red'), lty = 1)
## Incubation period
plot(xx, 
     dlnorm(xx, incubation_period$mean, incubation_period$sd), 
     type = 'l', main = 'assumed incubation time', xlab = 'days', ylab = 'dens')
lines(xx, 
      dlnorm(xx, testpars$true_log_mean_inc, testpars$true_log_sd_inc), 
      col = 'red', lty = 2)
legend('topright', c('input', 'true'), col = c(1, 'red'), lty = 1)
dev.off()





## Generate synthetic observations  -------------------------------------------

# Define functions to draw n samples from each delay distribution
r_inc_dist <- function(n){rlnorm(n, testpars$true_log_mean_inc, testpars$true_log_sd_inc)}
# Define delays to observation. Assume lognormal form if an alternate form is not specified on input
r_true_case_delay <- function(n){rlnorm(n, testpars$true_log_mean_case_delay, testpars$true_log_sd_case_delay)}
r_true_death_delay <- function(n){rlnorm(n, testpars$true_log_mean_death_delay, testpars$true_log_sd_death_delay)} # Additional delay from symptoms -> observation

outpatient_delay_dist <- function(nn){r_inc_dist(nn) + r_true_case_delay(nn)}
death_delay_dist <- function(nn){r_inc_dist(nn) + r_true_death_delay(nn)}

source('../00-inf_to_obs.R')
## Wrapper function to get synthetic data in a data frame that includes times of observation
get_obs_ts <- function(mt = NULL, cdelay = outpatient_delay_dist, ddelay = death_delay_dist){ # "Max time" - time at which to end the time series of observations
  if(length(mt) == 0){mt = max(get_sim_df()$time)}
  ## Sample from delays to get times of observation
  get_sim_df() %>%
    merge(
      get_tObs_from_tInf(get_sim_df()$incidence, 
                         get_sim_df()$time, 
                         cdelay, 
                         return_times = T),
      by = 'time', all = TRUE) %>% rename(obs_outpatient = n) %>%
    merge(
      get_tObs_from_tInf(get_sim_df()$incidence, 
                         get_sim_df()$time, 
                         ddelay, 
                         return_times = T),
      by = 'time', all = TRUE) %>% rename(obs_deaths = n) %>%
    as.tbl() -> obs_df
  
  ## Truncate at time mt
  obs_df <- obs_df %>% filter(time <= mt)
  return(obs_df)
}

## Generate synthetic times of observation up to t = 150
obs_df <- get_obs_ts(max_time) %>%
  mutate(date = Sys.Date()-rev(time)) ## Add calendar dates for input into epinow

## Plot
obs_df %>%
  select(date, incidence, contains('obs')) %>%
  mutate(shifted_observed = lead(obs_outpatient, n = round(exp(incubation_period$mean+incubation_period$sd^2/2)+exp(obs_rep_delay$mean+obs_rep_delay$sd^2/2))),
         shifted_deaths = lead(obs_deaths, n = round(exp(incubation_period$mean+incubation_period$sd^2/2)+exp(death_rep_delay$mean+death_rep_delay$sd^2/2)))) %>%
  pivot_longer(-date, names_to = 'data_type', values_to = 'count') %>%
  mutate(data_type = factor(data_type, 
                            levels = c('incidence', 'obs_outpatient', 'obs_deaths', 'shifted_observed', 'shifted_deaths'),
                            labels = c('latent (true) infections\nlater observed as...', 'cases', 'deaths', 'shifted cases', 'shifted deaths'))) %>%
  ggplot() +
  geom_line(aes(x = date, y = count, color = data_type)) +
  ggtitle('Synthetic Data') +
  scale_color_viridis_d(direction = -1) -> obs
obs + guides(color = 'legend') + theme(legend.position = 'right')
ggsave2(sprintf('figs/%s-synthetic_obs.png', testpars$output_folder), width = 6, height = 4, units = 'in', dpi = 300)



## Input into inference model  -------------------------------------------
## Write a wrapper to reformat the desired synthetic data for input into epiEstim
get_in_dat <- function(obs_colname, odf = obs_df){
  odf[,c('date', obs_colname)] %>%
    setnames(c('date', 'confirm')) %>%
    as.data.table() %>%
    return()
}







## Set intputs into estimate_infections
reported_cases = get_in_dat('obs_outpatient')
family = 'negbin'
generation_time = generation_time  ## Generation time priors
incubation_period = incubation_period ## Incubation priors
reporting_delay = obs_rep_delay ## Reporting dealy priors
gp = list(basis_prop = 0.5, boundary_scale = 2)
rt_prior = list(mean = 1, sd = 1) 
adapt_delta = 0.99
max_treedepth = 15
cores = 4
chains = 4 
samples = 2000 
warmup = 500
estimate_rt = TRUE
horizon = 7
verbose = TRUE
return_fit = TRUE
debug = FALSE