## Function to input parameters, plot true inputs vs. inputs specified for epinow2, and then run inference in stan.


run_test <- function(parlist,  # List of parameters used to generate synthetic data
                     testpars, # List of parameters used as inputs to epinow2 (may be misspecified relative to parlist)
                     max_time = 150, # Late timepoint in observations
                     r_case_dist = NULL, # Function that draws samples from the delay from symptom onset to case  or death observation. Default is lognormal, with params given by testpars. If another form is desired (e.g. uniform), specify an alternte function here.
                     r_death_dist = NULL,
                     d_case_dist = NULL,
                     d_death_dist = NULL
){
  
  
  ## Set delay distributions for input into EpiNow2 -------------------------------------------
  generation_time <- list(mean = testpars$input_mean_gi, # gamma
                          mean_sd = 1,
                          sd = testpars$input_sd_gi,
                          sd_sd = 1,
                          max = 30)
  
  incubation_period <- list(mean = log(testpars$input_mean_inc), # lognormal
                            mean_sd = EpiNow2::covid_incubation_period[1, ]$mean_sd,
                            sd = log(testpars$input_sd_inc),
                            sd_sd = EpiNow2::covid_incubation_period[1, ]$sd_sd,
                            max = 30)
  
  obs_rep_delay <- list(mean = log(testpars$input_mean_case_delay), ## lognormal, assume mean of 5 days for outpatient testing 
                        mean_sd = .5,   ## (48h from symptom onset + 48h for reporting + report the following day)
                        sd = log(testpars$input_sd_case_delay),
                        sd_sd = .5,
                        max = 30)
  
  death_rep_delay <- list(mean = log(testpars$input_mean_death_delay), ## lognorma 
                          mean_sd = .5, #For death, assume a mean of 15 days from symptom onset and an sd of about 7d
                          sd = log(testpars$input_sd_death_delay),
                          sd_sd = .5,
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
  
  if(length(d_case_dist) == 0){d_case_dist <- function(x){dlnorm(x, log(testpars$true_mean_case_delay), log(testpars$true_sd_case_delay))}}
  if(length(d_death_dist) == 0){d_death_dist <- function(x){dlnorm(x, log(testpars$true_mean_death_delay), log(testpars$true_sd_death_delay))}} # Additional delay from symptoms -> observation
  
  ## Delays to case observation
  plot(xx, d_case_dist(xx), 
       type = 'l', main = 'assumed delay from onset to case detection', xlab = 'days', ylab = 'dens')
  lines(xx, 
        dlnorm(xx, log(testpars$true_mean_case_delay), log(testpars$true_sd_case_delay)), 
        col = 'red', lty = 2)
  legend('topright', c('input', 'true'), col = c(1, 'red'), lty = 1)
  ## Delays to death
  plot(xx, d_death_dist(xx),
       type = 'l', main = 'assumed delay from onset to death', xlab = 'days', ylab = 'dens')
  lines(xx, 
        dlnorm(xx, log(testpars$true_mean_death_delay), log(testpars$true_sd_death_delay)), 
        col = 'red', lty = 2)
  legend('topright', c('input', 'true'), col = c(1, 'red'), lty = 1)
  ## Generation interval
  plot(xx, 
       dgamma(xx, shape = with(testpars, get_shape(input_mean_gi, input_sd_gi^2)), 
              rate = with(testpars, get_rate(input_mean_gi, input_sd_gi^2))), 
       type = 'l', main = 'assumed generation time', xlab = 'days', ylab = 'dens')
  lines(xx, 
        dgamma(xx, get_shape(parlist$true_mean_GI, parlist$true_var_GI), 
               get_rate(parlist$true_mean_GI, parlist$true_var_GI)), 
        col = 'red', lty = 2)
  legend('topright', c('input', 'true'), col = c(1, 'red'), lty = 1)
  ## Incubation period
  plot(xx, 
       dlnorm(xx, incubation_period$mean, incubation_period$sd), 
       type = 'l', main = 'assumed incubation time', xlab = 'days', ylab = 'dens')
  lines(xx, 
        dlnorm(xx, log(testpars$true_mean_inc), log(testpars$true_sd_inc)), 
        col = 'red', lty = 2)
  legend('topright', c('input', 'true'), col = c(1, 'red'), lty = 1)
  dev.off()
  
  
  
  
  
  ## Generate synthetic observations  -------------------------------------------
  
  # Define functions to draw n samples from each delay distribution
  r_inc_dist <- function(n){rlnorm(n, meanlog = log(testpars$true_mean_inc), sdlog = log(testpars$true_sd_inc))}
  # Define delays to observation. Assume lognormal form if an alternate form is not specified on input
  if(length(r_case_dist) == 0){r_case_dist <- function(n){rlnorm(n, log(testpars$true_mean_case_delay), log(testpars$true_sd_case_delay))}}
  if(length(r_death_dist) == 0){r_death_dist <- function(n){rlnorm(n, log(testpars$true_mean_death_delay), log(testpars$true_sd_death_delay))}} # Additional delay from symptoms -> observation
  
  outpatient_delay_dist <- function(nn){r_inc_dist(nn) + r_case_dist(nn)}
  death_delay_dist <- function(nn){r_inc_dist(nn) + r_death_dist(nn)}
  
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
    pivot_longer(-date, names_to = 'data_type', values_to = 'count') %>%
    mutate(is_imputed = ifelse(!grepl(pattern = 'obs', x = data_type), 'observed', 'underlying'),
           data_type = factor(data_type, 
                              levels = c('incidence', 'obs_outpatient', 'imputed_hospital', 'obs_deaths'),
                              labels = c('latent (true) infections\nlater observed as...', 'cases', 'hospital admissions', 'deaths'))) %>%
    ggplot() +
    geom_line(aes(x = date, y = count, color = data_type, lty = is_imputed)) +
    ggtitle('Synthetic Data') +
    scale_color_viridis_d(direction = -1)+
    guides(color = 'none')+
    theme(legend.position = 'bottom')+
    labs(color = '', linetype = 'Type') -> obs
  obs + guides(color = 'legend') + theme(legend.position = 'right')
  ggsave2(sprintf('figs/%s-synthetic_obs.png', testpars$output_folder), width = 6, height = 4, units = 'in', dpi = 300)
  
  
  # obs_df %>%
  #   mutate(shifted_observed = lead(obs_outpatient, n = round(exp(incubation_period$mean)+exp(obs_rep_delay$mean)))) %>%
  #   pivot_longer(c(incidence, contains('obs'))) %>%
  #   ggplot()+
  #   geom_line(aes(x = date, y = value, color = name))
  # 
  # obs_df %>%
  #   mutate(shifted_observed = lead(obs_deaths, n = round(exp(incubation_period$mean)+exp(death_rep_delay$mean)))) %>%
  #   pivot_longer(c(incidence, contains('obs'))) %>%
  #   ggplot()+
  #   geom_line(aes(x = date, y = value, color = name))
  
  
  
  obs_df %>%
    mutate(shifted_observed = lag(obs_outpatient, n = round(exp(incubation_period$mean+incubation_period$sd^2/2)+exp(obs_rep_delay$mean+obs_rep_delay$sd^2/2)))) %>%
    pivot_longer(c(incidence, contains('obs'))) %>%
    ggplot()+
    geom_line(aes(x = date, y = value, color = name))
  
  xx %>%
    filter(grepl('infections', variable) | grepl('report', variable)) %>%
    ggplot(aes(x = date)) +
    geom_line(aes(y = mean, color = variable))+
    geom_ribbon(aes(ymin = lower, ymax = upper, fill = variable), alpha = .5)
  
  
  
  
  
  
  
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
                                    rt_prior = list(mean = 2, sd = 1), horizon = 7,
                                    samples = 2000, warmup = 500, cores = 4,
                                    chains = 4, verbose = TRUE,
                                    target_folder = paste0(testpars$output_folder, '/cases'))
  
  est_from_deaths <- EpiNow2::epinow(reported_cases = get_in_dat('obs_deaths'), 
                                     generation_time = generation_time,  ## Generation time priors
                                     incubation_period = incubation_period, ## Incubation priors
                                     reporting_delay = death_rep_delay, ## Reporting dealy priors
                                     rt_prior = list(mean = 2, sd = 1), horizon = 7,
                                     samples = 2000, warmup = 500, cores = 4,
                                     chains = 4, verbose = TRUE,
                                     target_folder = paste0(testpars$output_folder, '/deaths'))
  
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
