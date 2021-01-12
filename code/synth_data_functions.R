source('model_functions.R')
library(dplyr)

# Function for generating stochastic realization based on predicted case probabilities
generate_stochastic_realization <- function(predictions){
  grouping_vars <- c('country','observation_year','lineage')
  # If SurveillanceType present as a variable in predictions tibble, add it as a grouping var
  if("SurveillanceType" %in% names(predictions)){
    grouping_vars <- c(grouping_vars, "SurveillanceType")
  }
  
  realization <- predictions %>%
    group_by_at(.vars = vars(grouping_vars)) %>%
    # Simulate from model assuming same number of total cases per country/obs.year/lineage as observed
    mutate(imprinting_predicted_cases = rmultinom(n = 1, size = CLY_total_cases[1], prob = pred_case_prob)) %>%
    ungroup() 
  
  # Check that total number of simulated cases in each country/obs.year/lineage (CYL) is the same as the observed one
  check <- realization %>%  group_by_at(.vars = vars(grouping_vars))  %>%
    summarise(total_predicted_cases = sum(imprinting_predicted_cases), total_obs_cases = CLY_total_cases[1]) %>%
    ungroup() %>%
    summarise(check = sum(total_predicted_cases == total_obs_cases), n_combinations = n()) %>%
    summarise(sums_match = check/n_combinations) %>%
    summarise(check = (1 - sums_match) < 1e-10) %>% pull(check)
  stopifnot(check)
  return(realization)
}


generate_synthetic_data <- function(model_name, model_parameters, dem_plus_case_data, lineage_frequencies,
                                    intensity_scores, cutoff_age,maternal_ab_duration, season_incidence_curves,
                                    school_start_age, oldest_atk_rate_age, reporting_age_cutoff,
                                    precomputed_history_probs, stochastic){
  # Get model function from model name passed as an argument
  model_function <- get(paste(model_name, 'model', sep ='_'))
  
  # Get case probabilities and expected numbers of cases under the model
  synthetic_data <- model_function(parameters = model_parameters, dem_plus_case_data, lineage_frequencies,
                                   intensity_scores, cutoff_age, maternal_ab_duration, season_incidence_curves,
                                   school_start_age, oldest_atk_rate_age, reporting_age_cutoff,
                                   precomputed_history_probs)
  
  # If stochastic == T, replace model expectation by a stochastic realization of the model 
  if(stochastic == T){
    synthetic_data <- generate_stochastic_realization(synthetic_data)
  }else{
    # If stochastic == F, synthetic data are the rounded expected number of cases under the model 
    synthetic_data <- mutate(synthetic_data, imprinting_predicted_cases = round(imprinting_predicted_cases)) 
    # Rounding causes the total n. of predicted cases to be slightly different than the actual total number.
  }
  synthetic_data <- synthetic_data %>%
    rename(n_cases = obs_cases) %>%
    # Retain only variables present in original data, plus imprinting_predicted_cases
    select(one_of(c(colnames(dem_plus_case_data),
                    'imprinting_predicted_cases'))) %>%
    # Remove true number of observed cases
    select(-n_cases) %>% 
    # Rename "imprinting_predicted_cases" as "obs_cases" (i.e. treat simulation as synthetic data)
    rename(n_cases = imprinting_predicted_cases) %>%
    # Add minimum birth year for each cohort in each observation year
    mutate(minimum_birth_year = ifelse(cohort_type == 'birth_year',
                                       cohort_value, observation_year - cohort_value - 1))
  return(synthetic_data)
}


# Export model parameters and base data to a pars file:
# pars_file <- file(paste(output_directory,replicate, '_syn_data_',model_name, '_pars.txt', sep = ''))
# writeLines(c(paste("base_data_file = ", base_data_file_path, sep =''),
#              paste("model_name = '", model_name, "'", sep = ''),
#              paste(c("model_parameters = ", paste(model_parameters, collapse =',')), collapse = ''),
#              paste("stochastic = ", stochastic, sep = '')
#              ),
#            pars_file)
# close(pars_file)

