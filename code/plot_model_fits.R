#!/usr/bin/env Rscript
library(reshape2)
library(tidyr)
library(dplyr)
library(cowplot)
library(ggplot2)
library(assertthat)
library(stringr)

# Load model functions
source('model_functions.R')
source('merge_data.R')
# Load fixed basic parameters
source('basic_parameters.R')

# Function for calculating infection history probabilities
source('imprinting_probabilites.R')

# Function for generating stochastic multinomial realizations of model predictions
source('synth_data_functions.R')

# Plotting functions
source('base_plotting_functions.R')

# Alpha and Number of replicates for bootstrapping C.Is
n_CI_replicates <- 1000
CI_alpha <- 0.05

args = commandArgs(trailingOnly = T)

input_file = args[1] # Input file specifying case data, covariate data, etc.
# (should be in the same directory as likelihood profiles)

# Example of input file:
#case_data_path='../results/processed_data/case_data_nz_all_surveillance_untyped_assigned.csv'
#demographic_data_path='../results/processed_data/demographic_data.csv'
#intensity_scores_path='../results/processed_data/intensity_scores.csv'
#lineage_frequencies_path='../results/processed_data/lineage_frequencies_gisaid-genbank_noVicin1990s.csv'
#season_incidence_curves_path='../results/processed_data/season_incidence_curves.csv'
#start_birth_year=1959
#subset_region='New_Zealand'
#reporting_age_cutoff=1
#precomputed_history_probs_path='NA'
#initial_pars='initial_parameters/three_class_atk_rate_model_noVicin1990s.csv'

constrained_pars = args[2] # Comma-separated string of parameters to constrain, if any
constrained_par_values = args[3] # Comma-separated values of constrained parameters, if any
constraint_type = args[4] # Soft: look for best values of remaining parameters
                          # Hard: keep remaining parameters at global MLE

constrained_pars = str_split(constrained_pars, ',')[[1]]
constrained_par_values = as.numeric(str_split(constrained_par_values, ',')[[1]])

if(all(is.na(constrained_pars)) & !all(is.na(constrained_par_values)) |
   (!all(is.na(constrained_pars)) & all(is.na(constrained_par_values)))){
  stop('Parameter constraints incorrectly specified.')
}

n_pars_constrained <- length(constrained_pars)

source(input_file)
if(precomputed_history_probs_path == 'NA'){precomputed_history_probs_path <- NA}

profiles_directory <- paste0(dirname(input_file),'/')
plot_directory <- paste0(dirname(profiles_directory), '/model_fit_plots/')

if(!all(is.na(constrained_pars))){
  constraint_label <- paste0(constraint_type, '_constrained_',
                             paste0(paste(constrained_pars, constrained_par_values, sep = '-'),
                                    collapse = '_')
                             )
  plot_directory <- paste0(plot_directory, constraint_label, '/')
}

dir.create(plot_directory, showWarnings = FALSE)

main <- function(){
  case_data <- as_tibble(read.csv(case_data_path))
  case_data <- set_min_birth_year(case_data, start_birth_year)
  
  demographic_data <- as_tibble(read.csv(demographic_data_path))
  demographic_data <- normalize_demographic_data(demographic_data, start_birth_year)
  
  intensity_scores <- as_tibble(read.csv(intensity_scores_path))
  lineage_frequencies <- as_tibble(read.csv(lineage_frequencies_path))
  season_incidence_curves <- as_tibble(read.csv(season_incidence_curves_path))
  
  # Combine demographic and case data
  dem_plus_case_data <- merge_data(demographic_data, case_data)
  
  if(exists('min_age')){
    dem_plus_case_data <- impose_min_age(dem_plus_case_data, min_age)
  }
  
  # Read combined likelihood profile
  combined_profile <- as_tibble(read.csv(paste0(profiles_directory,
                                                'combined_likelihood_profiles.csv'), header = T))
  # Read model name from profile tibble
  model <- as.character(unique(combined_profile$model))
  stopifnot(length(model) == 1)
  
  # Model function
  model_function <- get(paste(model, '_model', sep = ''))
  
  # Parameter names
  model_par_names <- get(paste(model, '_model_par_names', sep = ''))  
  
  mle_pars <- get_MLE_parameters(combined_profile)
  
  stopifnot(sum(is.na(mle_pars)) == 0)
  
  write.csv(tibble(par = model_par_names, value = mle_pars),
            paste0(plot_directory,'pars.csv'), row.names = F)
  
  # If no precomputed history probabilities given, compute them
  if(is.na(precomputed_history_probs_path)){
    # Compute and export history probabilities under the MLE
    R_V <- mle_pars[1]
    R_Y <- mle_pars[2]
    chi_V <- mle_pars[3]
    chi_Y <- mle_pars[4]
    gamma_VY <- mle_pars[5]
    gamma_YV <-  mle_pars[6]
    gamma_AV <- mle_pars[7]
    gamma_AY <- mle_pars[8]
    beta1 <- mle_pars[9]
    beta2 <-  mle_pars[10]
    beta3 <-  mle_pars[11]
    
    if(model == 'main'){
      reporting_factor_vic = mle_pars[model_par_names == 'reporting_factor']
      reporting_factor_yam = mle_pars[model_par_names == 'reporting_factor']
    }else{
      stopifnot(model == 'two_rhos')
      reporting_factor_vic = mle_pars[model_par_names == 'reporting_factor_vic']
      reporting_factor_yam = mle_pars[model_par_names == 'reporting_factor_yam']
    }
    
    chi_VY = chi_Y * gamma_VY
    chi_YV = chi_V * gamma_YV
    chi_AV = chi_V * gamma_AV
    chi_AY = chi_Y * gamma_AY
    
    base_iprobs_tibble <- as_tibble(expand.grid(
      observation_year = seq(min(dem_plus_case_data$observation_year), max(dem_plus_case_data$observation_year)),
      cohort_value = seq(min(dem_plus_case_data$cohort_value), max(dem_plus_case_data$cohort_value)))) %>%
      mutate(cohort_type = 'age', country = unique(dem_plus_case_data$country),
             region = NA, lineage = NA, n_cases = NA, CLY_total_cases = NA, rel_pop_size = NA) %>%
      arrange(observation_year)
    
    
     history_probs <- calculate_iprobs(dem_plus_case_data = base_iprobs_tibble,
                                      lineage_frequencies = lineage_frequencies,
                                      intensity_scores = intensity_scores,
                                      chi_VY = chi_VY, chi_YV = chi_YV,
                                      chi_AV = chi_AV, chi_AY = chi_AY,
                                      beta1 = beta1, beta2 = beta2, beta3 = beta3,
                                      cutoff_age = cutoff_age,
                                      maternal_ab_duration = maternal_ab_duration,
                                      season_incidence_curves = season_incidence_curves,
                                      school_start_age = school_start_age,
                                      oldest_atk_rate_age = oldest_atk_rate_age,
                                      birth_year_cutoff = birth_year_cutoff)
    
    
    history_probs <- history_probs %>% select(country, observation_year, cohort_type, cohort_value,
                                              matches('P_')) %>%
      select(-rel_pop_size) %>% unique()
    
    write.csv(history_probs, paste0(dirname(profiles_directory),'/MLE_history_probs.csv'), row.names = F)
    history_probs_constrained <- NULL # Don't precompute for constrained parameters
    
  }else{
    # If fit was done with precomputed history probabilities (i.e., non-sentinel or sentinel specific fits)
    history_probs <- as_tibble(read.csv(precomputed_history_probs_path, header = T))
    history_probs_constrained <- history_probs 
  }
  
  
  # Get model predictions
  predictions <- model_function(parameters = mle_pars, dem_plus_case_data, lineage_frequencies, intensity_scores,
                                cutoff_age, maternal_ab_duration, season_incidence_curves, school_start_age,
                                oldest_atk_rate_age, reporting_age_cutoff, precomputed_history_probs = history_probs)
 
  predictions = list(predictions)
  
  # Predictions with constrained parameters, if any
  if(!all(is.na(constrained_pars))){
    if(constraint_type == 'soft'){
      # Soft constrained parameters: a parameter is fixed, others are conditional MLEs given constraint
      constrained_pars = get_MLE_parameters(combined_profile,
                                            constrained_pars, constrained_par_values)
    }else{
      stopifnot(constraint_type == 'hard')
      constrained_pars = get_hard_constrained_pars(mle_pars, constrained_pars, constrained_par_values, model_par_names)
    }
    
    write.csv(tibble(par = model_par_names, global_mle = get_MLE_parameters(combined_profile),
                     constrained_values = constrained_pars),
              paste0(plot_directory,'pars.csv'), row.names = F)
    
    constrained_pars[is.na(constrained_pars)] <- 0
    
    constrained_predictions <- model_function(parameters = constrained_pars, dem_plus_case_data, lineage_frequencies,
                                                   intensity_scores,cutoff_age, maternal_ab_duration, season_incidence_curves,
                                                   school_start_age, oldest_atk_rate_age, reporting_age_cutoff,
                                                   precomputed_history_probs = history_probs_constrained)
    
    # Merge global MLE and constrained predictions into single predictions list object
    predictions <- list(predictions[[1]], constrained_predictions)
    names(predictions) <- c('MLE',constraint_label)
    
    # Log-likelihood by observation year, country and lineage, for the global MLE and for the constrained mle
    loglik_by_obs_year_global_MLE <- neg_loglik_function_by_obs_year(parameters = mle_pars, model_function, dem_plus_case_data,
                                                          lineage_frequencies, intensity_scores,cutoff_age, maternal_ab_duration,
                                                          season_incidence_curves, school_start_age, oldest_atk_rate_age, reporting_age_cutoff,
                                                          precomputed_history_probs = history_probs)
    loglik_by_obs_year_constrained <- neg_loglik_function_by_obs_year(parameters = constrained_pars, model_function, dem_plus_case_data,
                                                                     lineage_frequencies, intensity_scores,cutoff_age, maternal_ab_duration,
                                                                     season_incidence_curves, school_start_age, oldest_atk_rate_age,
                                                                     reporting_age_cutoff, precomputed_history_probs = history_probs_constrained)
    # Comparison of log-likelihood between MLE and constrained parameters, year by year.
    loglik_by_obs_year <- left_join(loglik_by_obs_year_global_MLE %>% rename(loglik_global_MLE = loglik),
                                    loglik_by_obs_year_constrained %>% rename(loglik_constrained = loglik),
                                    by = c('country', 'lineage', 'observation_year')) %>%
      mutate(LRS = 2*(loglik_global_MLE - loglik_constrained)) %>%
      mutate(p = 1 - pchisq(LRS, n_pars_constrained)) %>%
      mutate(global_LRS_without_this_line = sum(LRS) - LRS)
    
    write.csv(loglik_by_obs_year,
              paste0(plot_directory,'loglik_by_obs_year.csv'), row.names = F)
    
    write.csv(loglik_by_obs_year %>% mutate(loglik_dif = LRS/2) %>% 
                group_by(lineage) %>%
                summarise(loglik_dif = sum(loglik_dif),
                          loglik_global_MLE = sum(loglik_global_MLE),
                          loglik_constrained = sum(loglik_constrained)),
              paste0(plot_directory,'fit_comparison.csv'),
              row.names = F)
    
  }
  
  # Plot with fraction of cases aggregated across obs. years and countries
  # (We're no longer fitting to multiple countries simultaneously)
  fraction_cases_pooled <- plot_pooling_countries(predictions, n_CI_replicates, CI_alpha,plot_predictions = T,
                                                  plot_fraction = T, demographic_normalization = F) 
  save_plot(paste0(plot_directory,'fraction_cases_pooled.pdf'),
            fraction_cases_pooled + theme(legend.position = 'None') +
              scale_x_continuous(breaks = c(1960,1980,2000,2020),
                                 limits = c(1958,2020)),
            base_height = 6, base_width =11)
  
  fraction_cases_pooled_normalized <- plot_pooling_countries(predictions, n_CI_replicates, CI_alpha,plot_predictions = T,
                                                             plot_fraction = T, demographic_normalization = T) +
    geom_hline(yintercept = 1, linetype = 2)
  
  
  save_plot(paste0(plot_directory,'fraction_cases_pooled_normalized.pdf'),
            fraction_cases_pooled_normalized +
              ylab('Fraction of cases divided by fraction\nof the population in cohort'),
            base_height = 4, base_width = 8)
  
  # Cases by observation year
  cases_by_obs_year <- plot_by_obs_year(predictions, n_CI_replicates, CI_alpha, plot_predictions = T,
                                           plot_excess_cases = F, plot_fraction = F, demographic_normalization = F)
  save_plot(paste0(plot_directory,'cases_by_obs_year.pdf'),
            cases_by_obs_year[[1]],
            base_height = 10, base_width = 5)
  
  # Cases by observation year by age  
  cases_by_obs_year_by_age <- plot_by_obs_year(predictions, n_CI_replicates, CI_alpha, plot_predictions = T,
                                           plot_excess_cases = F, plot_fraction = F, demographic_normalization = F,
                                           plot_cohort_type = 'age')
  save_plot(paste0(plot_directory,'cases_by_obs_year_by_age.pdf'),
            cases_by_obs_year_by_age[[1]] + xlab('Age (years)'),
            base_height = 10, base_width = 5)
  
  # Excess cases by observation year
  excess_cases_by_obs_year <- plot_by_obs_year(predictions, n_CI_replicates, CI_alpha,plot_predictions = F,
                                               plot_excess_cases = T, plot_fraction = F, demographic_normalization = F)
  
  save_plot(paste0(plot_directory,'excess_cases_by_obs_year.pdf'),
            plot_grid(plotlist = excess_cases_by_obs_year, ncol = length(excess_cases_by_obs_year)),
            base_height = 13, base_width =6.5 * length(excess_cases_by_obs_year))
  
  save_plot(paste0(plot_directory,'fraction_cases_pooled_panel.pdf'),
            plot_grid(fraction_cases_pooled + ylab(''),
                      fraction_cases_pooled_normalized + ylab('') +
                        scale_y_continuous(breaks = c(0,2,4,6),
                                           labels = c('0.00','2.00','4.00','6.00'),
                                           limits = c(0,6)),
                        #scale_y_continuous(breaks = seq(-0.03,0.06,0.03)) +
                        #theme(axis.text.y = element_text(size = 8.7)),
                      nrow = 2),
            base_height = 10, base_width =9)
  
  # Plots with raw (absolute number of cases) distributions
  raw_dist_by_obs_year <- plot_by_obs_year(predictions, n_CI_replicates = 1,CI_alpha,plot_predictions = F,
                                           plot_fraction = F, demographic_normalization = F)
  
  # Plots with fraction of cases normalized by demography by obs. year
  #fraction_by_obs_year <- plot_by_obs_year(predictions, n_CI_replicates, CI_alpha,plot_predictions = T,
  #                                         plot_fraction = T, demographic_normalization = T)
  # Plots with number of cases normalized by demography by obs. year
  cases_by_obs_year <- plot_by_obs_year(predictions, n_CI_replicates, CI_alpha,plot_predictions = T,
                                           plot_fraction = T, demographic_normalization = T)
  
  # Plots with number of cases normalized by demography by interval
  #cases_by_interval <- plot_by_interval(predictions, n_CI_replicates, CI_alpha,plot_predictions = T,
  #                                      plot_fraction = F, demographic_normalization = T)
  
  # Save all lists with separate plots by country
  for(i in 1:length(raw_dist_by_obs_year)){
    country <- names(raw_dist_by_obs_year)[[i]]
    save_plot(paste0(plot_directory,'raw_cases_by_obs_year_',
                     country,'.pdf'),
              raw_dist_by_obs_year[[i]],
              base_height = 10, base_width =7)
    #save_plot(paste0(plot_directory,'cases_by_interval_',
    #                 country,'.pdf'),
    #          cases_by_interval[[i]],
    #          base_height = 7, base_width =7)
  }
  
  
}

main()

  
