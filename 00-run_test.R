## Function to input parameters, plot true inputs vs. inputs specified for epinow2, and then run inference in stan.


run_test <- function(parlist,  # List of parameters used to generate synthetic data
                     testpars, # List of parameters used as inputs to epinow2 (may be misspecified relative to parlist)
                     max_time = 150,
                     prior_smoothing_window,
                     debug = FALSE)
  {
  
  
  ## Set delay distributions for input into EpiNow2 -------------------------------------------
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
  
  
  
  write_rds(testpars, sprintf('%s/testpars.rds', testpars$output_folder))
  write_rds(generation_time, sprintf('%s/gen_interval.rds', testpars$output_folder))
  write_rds(incubation_period, sprintf('%s/incubation_pd.rds', testpars$output_folder))
  write_rds(obs_rep_delay, sprintf('%s/case_delay.rds', testpars$output_folder))
  write_rds(death_rep_delay, sprintf('%s/death_delay.rds', testpars$output_folder))
  
  
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
  lines(xx, true_case_dist(xx), col = 'red', lty = 2)
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
  
  
  if(!debug){
  ## Fit to synthetic case observations
  est_from_cases <- EpiNow2::epinow(reported_cases = get_in_dat('obs_outpatient'), 
                                    generation_time = generation_time,
                                    delays = list(reporting_delay = obs_rep_delay,
                                                  incubation_period = incubation_period), 
                                    prior_smoothing_window = prior_smoothing_window,
                                    rt_prior = list(mean = 2, sd = 1), horizon = 7,
                                    samples = 2000, warmup = 500, cores = 4,
                                    chains = 4, verbose = TRUE,
                                    target_folder = paste0(testpars$output_folder, '/cases'))
  
  est_from_deaths <- EpiNow2::epinow(reported_cases = get_in_dat('obs_deaths'), 
                                     generation_time = generation_time,  ## Generation time priors
                                     delays = list(reporting_delay = obs_rep_delay,
                                                   incubation_period = incubation_period), 
                                     rt_prior = list(mean = 2, sd = 1), horizon = 7,
                                     samples = 2000, warmup = 500, cores = 4,
                                     chains = 4, verbose = TRUE,
                                     prior_smoothing_window = prior_smoothing_window,
                                     target_folder = paste0(testpars$output_folder, '/deaths'))
  }
  
}


## For debugging
# r_case_dist = function(nn){runif(n = nn, min = 0, max = 10)}
# r_death_dist = function(nn){runif(n = nn, min = 5, max = 25)}
# d_case_dist = function(xx) ifelse(xx >=0 & xx <= 10, 1/10, 0)
# d_death_dist = function(xx) ifelse(xx >=5 & xx <= 25, 1/20, 0)
# r_case_dist = NULL 
# r_death_dist = NULL
# d_case_dist = NULL
# d_death_dist = NULL



#run_test(parlist = parlist, testpars = testpars)
