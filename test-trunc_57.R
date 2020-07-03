## Estimate Rt using epinow2

## Load dependencies and set parameters ------------------------------------------------------
rm(list = ls())
source('load_packages.R')
source('util.R')
if(!file.exists('true_pars.rds')){source('01-make_SEIR_data.R')}
parlist <- readRDS('true_pars.rds')
ggplot2::theme_set(theme_bw())

# Load synthetic data using:
# get_sim_df()


## Set parameters
testpars <- list(
  last_obs_time = 57,
  output_folder = 'truncated_57',
  ## True delays
  true_mean_case_delay = 5,
  true_sd_case_delay = 1.7,
  true_mean_death_delay = 15,
  true_sd_death_delay = 1.5,
  true_mean_inc = exp(EpiNow2::covid_incubation_period[1, ]$mean),
  true_sd_inc = exp(EpiNow2::covid_incubation_period[1, ]$sd))
## Delays specified in model
testpars$input_mean_case_delay = testpars$true_mean_case_delay
testpars$input_sd_case_delay = testpars$true_sd_case_delay
testpars$input_mean_death_delay = testpars$true_mean_death_delay
testpars$input_sd_death_delay = testpars$true_sd_death_delay
testpars$input_mean_inc = testpars$true_mean_inc
testpars$input_sd_inc = testpars$true_sd_inc
testpars$input_mean_gi = parlist$true_mean_GI
testpars$input_sd_gi = sqrt(parlist$true_var_GI)

dir_check(testpars$output_folder)



## If aiming to input a perfectly specified generation interval, calculate the appropriate lognormal parameters -------------------------------------------
if(testpars$input_mean_gi == parlist$true_mean_GI & testpars$input_sd_gi == sqrt(parlist$true_var_GI)){
  # The stan model assumes a lognormal GI. The true GI is gamma-distributed. Fit the true GI to a lognormal for input into stan.
  # Draw samples from the true GI
  gi_samples = rgamma(1000, get_shape(parlist$true_mean_GI, parlist$true_var_GI), get_rate(parlist$true_mean_GI, parlist$true_var_GI)) 
  ## Fit the true generation time to a lognormal
  lnormal_gi_pars <- optim(par = c(mm = log(8), sd = 1), fn = function(pars){
    -sum(dlnorm(gi_samples, pars['mm'], pars['sd'], log = TRUE))
  })
testpars$input_mean_gi = exp(lnormal_gi_pars$par[1])
testpars$input_sd_gi = exp(lnormal_gi_pars$par[2])
}



## Set delay distributions for input into EpiNow2 -------------------------------------------
generation_time <- list(mean = log(testpars$input_mean_gi),
                        mean_sd = 1,
                        sd = log(testpars$input_sd_gi),
                        sd_sd = 1,
                        max = 30)

incubation_period <- list(mean = log(testpars$input_mean_inc),
                          mean_sd = EpiNow2::covid_incubation_period[1, ]$mean_sd,
                          sd = log(testpars$input_sd_inc),
                          sd_sd = EpiNow2::covid_incubation_period[1, ]$sd_sd,
                          max = 30)

obs_rep_delay <- list(mean = log(testpars$input_mean_case_delay), ## assume mean of 5 days for outpatient testing 
                      mean_sd = 1,   ## (48h from symptom onset + 48h for reporting + report the following day)
                      sd = log(testpars$input_sd_case_delay),
                      sd_sd = 1,
                      max = 30)

death_rep_delay <- list(mean = log(testpars$input_mean_death_delay), ## For death, assume a mean of 15 days from symptom onset and an sd of about 7d
                        mean_sd = 1,
                        sd = log(testpars$input_sd_death_delay),
                        sd_sd = 1,
                        max = 30)


write_rds(testpars, sprintf('%s/testpars.rds', testpars$output_folder))
write_rds(generation_time, sprintf('%s/gen_interval.rds', testpars$output_folder))
write_rds(incubation_period, sprintf('%s/incubation_pd.rds', testpars$output_folder))
write_rds(obs_rep_delay, sprintf('%s/case_delay.rds', testpars$output_folder))
write_rds(death_rep_delay, sprintf('%s/death_delay.rds', testpars$output_folder))


# Plot the specified distributions
pdf(sprintf('%s/specified_distributions.pdf', testpars$output_folder))
xx = seq(0, 30, by = 0.01)

plot(xx, dlnorm(xx, obs_rep_delay$mean, obs_rep_delay$sd), type = 'l', main = 'assumed delay from onset to case detection', xlab = 'days', ylab = 'dens')
lines(xx, dlnorm(xx, log(testpars$true_mean_case_delay), log(testpars$true_sd_case_delay)), col = 'red', lty = 2)
legend('topright', c('input', 'true'), col = c(1, 'red'), lty = 1)

plot(xx, dlnorm(xx, death_rep_delay$mean, death_rep_delay$sd), type = 'l', main = 'assumed delay from onset to death', xlab = 'days', ylab = 'dens')
lines(xx, dlnorm(xx, log(testpars$true_mean_death_delay), log(testpars$true_sd_death_delay)), col = 'red', lty = 2)
legend('topright', c('input', 'true'), col = c(1, 'red'), lty = 1)

plot(xx, dlnorm(xx, generation_time$mean, generation_time$sd), type = 'l', main = 'assumed generation time', xlab = 'days', ylab = 'dens')
  lines(xx, 
        dgamma(xx, get_shape(parlist$true_mean_GI, parlist$true_var_GI), 
               get_rate(parlist$true_mean_GI, parlist$true_var_GI)), 
        col = 'red', lty = 2)
  legend('topright', c('input', 'true'), col = c(1, 'red'), lty = 1)

plot(xx, dlnorm(xx, incubation_period$mean, incubation_period$sd), type = 'l', main = 'assumed incubation time', xlab = 'days', ylab = 'dens')
lines(xx, dlnorm(xx, log(testpars$true_mean_inc), log(testpars$true_sd_inc)), col = 'red', lty = 2)
legend('topright', c('input', 'true'), col = c(1, 'red'), lty = 1)
dev.off()





## Generate synthetic observations  -------------------------------------------

# Define functions to draw n samples from each delay distribution
r_inc_dist <- function(n){rlnorm(n, meanlog = incubation_period$mean, sdlog = incubation_period$sd)} 

outpatient_delay_dist <- function(nn, inc_dist = r_inc_dist){
  r_sym_to_obs_dist <- function(n){rlnorm(n, obs_rep_delay$mean, obs_rep_delay$sd)}
  inc_dist(nn) + r_sym_to_obs_dist(nn)
}

death_delay_dist <- function(nn, inc_dist = r_inc_dist){
  r_sym_to_obs_dist <- function(n){rlnorm(n, death_rep_delay$mean, death_rep_delay$sd)} # Additional delay from symptoms -> observation
  inc_dist(nn) + r_sym_to_obs_dist(nn)
}

source('inf_to_obs.R')
## Wrapper function to get synthetic data in a data frame that includes times of observation
get_obs_ts <- function(mt = NULL){ # "Max time" - time at which to end the time series of observations
  if(length(mt) == 0){mt = max(get_sim_df()$time)}
  ## Sample from delays to get times of observation
  get_sim_df() %>%
    merge(
      get_tObs_from_tInf(get_sim_df()$incidence, 
                         get_sim_df()$time, 
                         outpatient_delay_dist, 
                         return_times = T),
      by = 'time', all = TRUE) %>% rename(obs_outpatient = n) %>%
    merge(
      get_tObs_from_tInf(get_sim_df()$incidence, 
                         get_sim_df()$time, 
                         death_delay_dist, 
                         return_times = T),
      by = 'time', all = TRUE) %>% rename(obs_deaths = n) %>%
    as.tbl() -> obs_df
  
  ## Truncate at time mt
  obs_df <- obs_df %>% filter(time <= mt)
  return(obs_df)
}

## Generate synthetic times of observation up to t = 150
obs_df <- get_obs_ts(testpars$last_obs_time) %>%
  mutate(date = Sys.Date()-rev(time)) ## Add calendar dates for input into epinow

## Plot
obs_df %>%
  select(time, incidence, contains('obs')) %>%
  pivot_longer(-time, names_to = 'data_type', values_to = 'count') %>%
  mutate(is_imputed = ifelse(!grepl(pattern = 'obs', x = data_type), 'observed', 'underlying'),
         data_type = factor(data_type, 
                            levels = c('incidence', 'obs_outpatient', 'imputed_hospital', 'obs_deaths'),
                            labels = c('latent (true) infections\nlater observed as...', 'cases', 'hospital admissions', 'deaths'))) %>%
  ggplot() +
  geom_line(aes(x = time, y = count, color = data_type, lty = is_imputed)) +
  ggtitle('Synthetic Data') +
  scale_color_viridis_d(direction = -1)+
  guides(color = 'none')+
  theme(legend.position = 'bottom')+
  labs(color = '', linetype = 'Type') -> obs
obs + guides(color = 'legend') + theme(legend.position = 'right')
ggsave2(sprintf('%s/synthetic_obs.png', testpars$output_folder), width = 7, height = 5, units = 'in', dpi = 300)







## Input into inference model  -------------------------------------------
## Write a wrapper to reformat the desired synthetic data for input into epiEstim
get_in_dat <- function(obs_colname, odf = obs_df){
  odf[,c('date', obs_colname)] %>%
    setnames(c('date', 'confirm')) %>%
    as.data.table() %>%
    return()
}



## Fit to synthetic case observations
est_from_cases <- EpiNow2::epinow(reported_cases = get_in_dat('obs_outpatient'), 
                                  generation_time = generation_time,  ## Generation time priors
                                  incubation_period = incubation_period, ## Incubation priors
                                  reporting_delay = obs_rep_delay, ## Reporting dealy priors
                                  rt_prior = list(mean = 1, sd = 1), horizon = 7,
                                  samples = 2000, warmup = 500, cores = 4,
                                  chains = 4, verbose = TRUE,
                                  target_folder = paste0(testpars$output_folder, '/cases'))

est_from_deaths <- EpiNow2::epinow(reported_cases = get_in_dat('obs_deaths'), 
                                   generation_time = generation_time,  ## Generation time priors
                                   incubation_period = incubation_period, ## Incubation priors
                                   reporting_delay = death_rep_delay, ## Reporting dealy priors
                                   rt_prior = list(mean = 1, sd = 1), horizon = 7,
                                   samples = 2000, warmup = 500, cores = 4,
                                   chains = 4, verbose = TRUE,
                                   target_folder = paste0(testpars$output_folder, '/deaths'))