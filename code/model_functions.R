library(dplyr)
library(stringr)

parameter_bounds <- tibble(
  parameter = c("R_V","R_Y","chi_V","chi_Y","gamma_VY","gamma_YV","gamma_AV", "gamma_AY","beta1","beta2",
                "beta3", "reporting_factor_us","reporting_factor_aus","reporting_factor_nz"),
  lower_bound = c(rep(0,11), rep(1e-4,3)),
  upper_bound = c(rep(1,8), rep(0.75,3), rep(50,3))
)

# Function for constraining data to minimum birth year (e.g. 1988 or 2000)
set_min_birth_year <- function(case_data, min_birth_year){
  case_data <- case_data %>%
    filter(minimum_birth_year >= min_birth_year) %>%
    select(-minimum_birth_year)
  return(case_data)
} 

# Function for constraining demographic data by minimum birth year and normalizing rel. fractions of each b.year
normalize_demographic_data <- function(demographic_data, min_birth_year){
  # Constrain by minimum birth year
  demographic_data <- demographic_data %>%
    filter(minimum_birth_year >= min_birth_year) %>%
    select(-minimum_birth_year)
  
  # Calculate fraction of population in each b.year/age/age group for each year in each country
  demographic_data <- demographic_data %>% group_by(observation_year, country) %>% mutate(year_total = sum(n_persons)) %>% 
    mutate(fraction = n_persons / year_total) %>%
    ungroup()
  # Check that the fractions in each class sum to one in each year in each country
  fraction_sums <- demographic_data %>% group_by(observation_year, country) %>% summarise(total = sum(fraction)) %>% pull(total)
  stopifnot(sum(abs(fraction_sums - 1) < 1e-7) == length(fraction_sums))
  return(demographic_data)
}

# Function for calculating min. and max. possible age for cohorts identified by age, b. year or age group
calculate_age_bounds <- function(data){
  modified_data <- data %>%
    # Calculate minimum and maximum possible age for each cohort 
    mutate(max_age = NA, min_age = NA) %>% rowwise() %>%
    mutate(max_age = ifelse(cohort_type == 'age_group',  # If cohort identified by age group (Canadian data)...
                            as.numeric(str_replace(str_extract(cohort_value, '-[0-9]+'), '-','')),
                            max_age),
           min_age = ifelse(cohort_type == 'age_group',  # If cohort identified by age group (Canadian data)...
                            as.numeric(str_replace(str_extract(cohort_value, '[0-9]+-'), '-','')),
                            min_age),
           # If cohort identified by age, use "cohort_value":
           max_age = ifelse(cohort_type == 'age',as.numeric(cohort_value), max_age),
           min_age = ifelse(cohort_type == 'age',as.numeric(cohort_value), min_age),
           # If id'ed. by b. year, estimate as max. age
           max_age = ifelse(cohort_type == 'birth_year', observation_year - as.numeric(cohort_value), max_age),
           min_age = ifelse(cohort_type == 'birth_year', observation_year - as.numeric(cohort_value) - 1, min_age)
    ) %>% ungroup()
  return(modified_data)
}

# Function for normalizing case probs. so they sum to 1 for each country/lineage/obs. year combination
normalize_case_probs <- function(predictions){
  grouping_vars <- c('country','observation_year','lineage')
  # If SurveillanceType present as a variable in predictions tibble, add it as a grouping var
  if("SurveillanceType" %in% names(predictions)){
    grouping_vars <- c(grouping_vars, "SurveillanceType")
  }
  # predictions: tibble passed by predicted_case_probability or predicted_case_probability_with_reassortment
  normalized_tibble <- predictions %>% 
    group_by_at(.vars = vars(grouping_vars)) %>%
    mutate(pred_case_prob = pred_case_prob / sum(pred_case_prob)) %>% ungroup() %>%
    # Calculate predicted case numbers by country/lineage/observation year totals
    mutate(demography_predicted_cases = rel_pop_size * CLY_total_cases,
           imprinting_predicted_cases = pred_case_prob * CLY_total_cases) %>%
    ungroup()
  return(normalized_tibble)
}

# Function for checking that total model probabilties by country/lineage/observation year sum to 1
check_total_probs <- function(predictions){
  grouping_vars <- c('country','observation_year','lineage')
  # If SurveillanceType present as a variable in predictions tibble, add it as a grouping var
  if("SurveillanceType" %in% names(predictions)){
    grouping_vars <- c(grouping_vars, "SurveillanceType")
  }
  # predictions: tibble passed by predicted_case_probability or predicted_case_probability_with_reassortment
  total_probs <- unique(predictions %>% group_by_at(.vars = vars(grouping_vars)) %>% 
                          summarise(total_prob = sum(pred_case_prob)) %>% pull(total_prob))
  stopifnot(all(abs(total_probs - 1) < 1e-10 ))
}

# Function for calculating relative susceptibility given probabilities of different exposure histories
calculate_relative_susceptibility <- function(data_plus_protection){
  # data_with_protection: tibble with cases, impr. probs and protection effects distributed
  output_tibble <- data_plus_protection %>%
    rowwise() %>%
    mutate(
      rel_susceptibility = ifelse(lineage == 'B/Victoria',
                                  P_0 +
                                    P_A0 * (1 - chi_AV) +
                                    P_Y0 * (1 - chi_YV) + 
                                    P_AY0 * (1 - max(chi_AV, chi_YV)) + 
                                    P_AV0 * (1 - chi_V) +
                                    P_AVY * (1 - chi_V) +
                                    P_YV * (1 - chi_V) +
                                    P_V0 * (1 - chi_V) * (1 - R_V) +
                                    P_VY * (1 - chi_V) * (1 - R_V),
                                  P_0 +
                                    P_A0 * (1 - chi_AY) +
                                    P_V0 * (1 - chi_VY) + 
                                    P_AV0 * (1 - max(chi_AY, chi_VY)) + 
                                    P_AY0 * (1 - chi_Y) +
                                    P_AVY * (1 - chi_Y) +
                                    P_VY * (1 - chi_Y) +
                                    P_Y0 * (1 - chi_Y) * (1 - R_Y) +
                                    P_YV * (1 - chi_Y) * (1 - R_Y)
                                  )
      ) %>%
    ungroup()
  return(output_tibble)
}

# Function for adding country-specific age effects to model predictions (renormalizes resulting probs)
add_country_reporting_effect <- function(predictions, reporting_factor_us, reporting_factor_aus,
                                         reporting_factor_nz, reporting_age_cutoff){
  
  country_effects <- c(reporting_factor_us, reporting_factor_aus, reporting_factor_nz)
  names(country_effects) <- c('United States', 'Australia', 'New Zealand')
  
  country_reporting_effect_predictions <- calculate_age_bounds(predictions) %>%
    # Distribute country-specific reporting effects
    mutate(country_reporting_effect = as.numeric(country_effects[as.character(country)])) %>%
    # Add risk modifiers for the relevant ages
    mutate(risk_modifier = ifelse(max_age <= reporting_age_cutoff, country_reporting_effect, 1)) %>%
    # Multiply predicted probabilities by age-based risk increase
    mutate(pred_case_prob = pred_case_prob*risk_modifier)  %>%
    select(-country_reporting_effect)
  
  # Re-normalize predicted probabilities
  country_reporting_effect_predictions <- normalize_case_probs(country_reporting_effect_predictions)
  check_total_probs(country_reporting_effect_predictions)
  return(country_reporting_effect_predictions)
}

### ==== Multiply age specific attack rate into multinomial probabilities
# beta1 if age < school_start_year, beta2 if school_start_year <= age < oldest_atk_rate_age, beta3 otherwise
add_age_specific_attack_rates <- function(predictions, beta1, beta2, beta3, school_start_age, oldest_atk_rate_age){
  
  # Distribute school start age
  predictions <- left_join(predictions, school_start_age, by = 'country')
  
  predictions <- calculate_age_bounds(predictions) 
  
 
  
  predictions <- predictions %>%
    rowwise() %>%
    mutate(age_specific_attack_rate = ifelse(max_age < school_start_age,
                                             beta1,
                                             ifelse(max_age < oldest_atk_rate_age, beta2, beta3))) %>%
    ungroup() %>%
    mutate(pred_case_prob = pred_case_prob * age_specific_attack_rate) %>%
    ungroup()
  
  predictions <- normalize_case_probs(predictions)
  check_total_probs(predictions)
  return(predictions)
}

# Function for calculating case multinomial probabilities
compute_case_probs <- function(data_plus_protection, beta1, beta2, beta3,reporting_factor_us, reporting_factor_aus,
                               reporting_factor_nz, school_start_age, oldest_atk_rate_age, reporting_age_cutoff){
  # data_plus_protection: tib. w/ demog. and case data, impr. probs. and protec. pars.; passed by model functions
  predictions <- data_plus_protection %>%
    mutate(pred_case_prob = rel_pop_size*rel_susceptibility)

  # Add age-specific attack rates
  predictions <- add_age_specific_attack_rates(predictions, beta1 = beta1, beta2 = beta2, beta3 = beta3,
                                              school_start_age = school_start_age, oldest_atk_rate_age = oldest_atk_rate_age)
  # Add reporting effect
  predictions <- add_country_reporting_effect(predictions, reporting_factor_us = reporting_factor_us,
                                              reporting_factor_aus = reporting_factor_aus,
                                              reporting_factor_nz = reporting_factor_nz,
                                              reporting_age_cutoff = reporting_age_cutoff)
  predictions <- normalize_case_probs(predictions)
  check_total_probs(predictions)
  return(predictions)
}

# Function for generating initial parameter values given model parameter names
init_pars <- function(par_names, initial_par_bounds){
  # Check that all parameter names are listed in initial_par_bounds
  stopifnot(all(par_names %in% initial_par_bounds$parameter))
  sample_init_value <- function(par, initial_par_bounds){
    bounds <- initial_par_bounds %>% filter(parameter == par) %>% select(lower_bound, upper_bound)
    return(runif(1, bounds$lower_bound, bounds$upper_bound))
  }
  return(sapply(par_names, sample_init_value, initial_par_bounds = initial_par_bounds))
}

# Functions for generating lower and upper bounds for parameter values given vector of parameter names
lower_bounds <- function(par_names, parameter_bounds){
  # Check that all parameter names are listed in parameter_bounds
  stopifnot(all(par_names %in% parameter_bounds$parameter))
  lbound <- function(par, parameter_bounds){
    bounds <- parameter_bounds %>% filter(parameter == par) %>% select(lower_bound, upper_bound)
    return(bounds$lower_bound)
  }
  return(sapply(par_names, lbound, parameter_bounds = parameter_bounds))
}

upper_bounds <- function(par_names, parameter_bounds){
  # Check that all parameter names are listed in parameter_bounds
  stopifnot(all(par_names %in% parameter_bounds$parameter))
  ubound <- function(par, parameter_bounds){
    bounds <- parameter_bounds %>% filter(parameter == par) %>% select(lower_bound, upper_bound)
    return(bounds$upper_bound)
  }
  return(sapply(par_names, ubound, parameter_bounds = parameter_bounds))
}

# Function implementing main model. Passed to loglik functions.
main_model <- function(parameters, dem_plus_case_data, lineage_frequencies,
                       intensity_scores, cutoff_age, maternal_ab_duration, season_incidence_curves,
                       school_start_age, oldest_atk_rate_age, reporting_age_cutoff){
  
  R_V <- parameters[1]
  R_Y <- parameters[2]
  chi_V <- parameters[3]
  chi_Y <- parameters[4]
  gamma_VY <- parameters[5]
  gamma_YV <- parameters[6]
  gamma_AV <- parameters[7]
  gamma_AY <- parameters[8]
  beta1 <- parameters[9]
  beta2 <- parameters[10]
  beta3 <- parameters[11]
  reporting_factor_us <- parameters[12]
  reporting_factor_aus <- parameters[13]
  reporting_factor_nz <- parameters[14]
  
  chi_VY = chi_Y * gamma_VY
  chi_YV = chi_V * gamma_YV
  chi_AV = chi_V * gamma_AV
  chi_AY = chi_Y * gamma_AY
  
  # Calculate exposure history probabilities
  data_plus_improbs <- calculate_iprobs(dem_plus_case_data = dem_plus_case_data,
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
  
  # Distribute chi's and R's as columns
  data_plus_protection <- data_plus_improbs %>%
    mutate(R_V, R_Y, chi_V, chi_Y, chi_VY, chi_YV, chi_AV, chi_AY)

  # Calculate relative susceptibilities
  data_plus_protection <- calculate_relative_susceptibility(data_plus_protection)
  
  # Compute case multinomial probabilities
  predictions <- compute_case_probs(data_plus_protection, beta1 = beta1, beta2 = beta2, beta3 = beta3,
                                    reporting_factor_us = reporting_factor_us, 
                                    reporting_factor_aus = reporting_factor_aus,
                                    reporting_factor_nz = reporting_factor_nz,
                                    school_start_age = school_start_age,
                                    oldest_atk_rate_age = oldest_atk_rate_age,
                                    reporting_age_cutoff = reporting_age_cutoff)
  
  return(predictions)
}

main_model_par_names <-   c("R_V","R_Y","chi_V","chi_Y","gamma_VY","gamma_YV","gamma_AV","gamma_AY","beta1",
                            "beta2", "beta3", "reporting_factor_us","reporting_factor_aus","reporting_factor_nz")


# Log-likelihood function by country / observation year / lineage combination (and surveillance type, if applicable)
neg_loglik_function_by_obs_year <- function(parameters, model_function, dem_plus_case_data,
                                            lineage_frequencies, intensity_scores,cutoff_age, maternal_ab_duration,
                                            season_incidence_curves, school_start_age, oldest_atk_rate_age, reporting_age_cutoff){
  
  multinomial_draw_unit <- c('country', 'lineage', 'observation_year')
  if("SurveillanceType" %in% names(dem_plus_case_data)){
    multinomial_draw_unit <- c(multinomial_draw_unit, 'SurveillanceType')
  }
  
  neg_loglik_by_obs_year <- model_function(parameters = parameters, dem_plus_case_data = dem_plus_case_data,
                                           lineage_frequencies = lineage_frequencies, intensity_scores = intensity_scores,
                                           cutoff_age = cutoff_age, maternal_ab_duration = maternal_ab_duration,
                                           season_incidence_curves = season_incidence_curves, school_start_age = school_start_age,
                                           oldest_atk_rate_age = oldest_atk_rate_age, reporting_age_cutoff = reporting_age_cutoff) %>% 
    # Calculate log-likelihood separately for each country/lineage/observation year
    group_by_at(vars(multinomial_draw_unit)) %>%
    # ...under a multinomial distribution
    summarise(loglik = dmultinom(x = obs_cases, prob = pred_case_prob, log = T)) %>%
    ungroup()
  return(neg_loglik_by_obs_year)
  
}

# Negative log-likelihood function (optim minimizes function by default) 
neg_loglik_function <- function(parameters, model_function, dem_plus_case_data,
                                lineage_frequencies, intensity_scores,cutoff_age, maternal_ab_duration,
                                season_incidence_curves, school_start_age, oldest_atk_rate_age, reporting_age_cutoff){
  
  # Compute predicted probabilities
  loglik <- neg_loglik_function_by_obs_year(parameters, model_function, dem_plus_case_data,
                                            lineage_frequencies, intensity_scores,cutoff_age, maternal_ab_duration,
                                            season_incidence_curves, school_start_age, oldest_atk_rate_age,
                                            reporting_age_cutoff) %>%
    # Sum log-likelihoods across countries/lineages/observation years
    summarise(loglik = sum(loglik)) %>% 
    # Pull total log-likelihood
    pull(loglik)
  gc()
  return(-loglik)
}



# Constrained negative log-lik function with specified parameters constrained to specified values
constrained_loglik_function <- function(free_parameters, model_function, model_par_names, selected_par_names,
                                        selected_par_values,dem_plus_case_data, lineage_frequencies,
                                        intensity_scores,cutoff_age, maternal_ab_duration, season_incidence_curves,
                                        school_start_age, oldest_atk_rate_age, reporting_age_cutoff, parameter_bounds){
  
  free_parameter_names = model_par_names[model_par_names %in% selected_par_names == F]
  # If any free parameter value is outside of bounds, set loglik to NA
  if(any(free_parameters - lower_bounds(free_parameter_names, parameter_bounds) < -1e-7) |
     any(free_parameters - upper_bounds(free_parameter_names, parameter_bounds) > 1e-7)){
    loglik <- NA # Optim's default is to minimize, not maximize, thus returning +Inf
  }else{
    full_parameters <- rep(NA, length(model_par_names))
    
    for(i in 1:length(selected_par_names)){
      full_parameters[model_par_names == selected_par_names[i]] <- selected_par_values[i]
    }
    full_parameters[model_par_names %in% selected_par_names == F] <- free_parameters
    
    loglik <- do.call(neg_loglik_function, args = list(parameters = full_parameters,
                                                       model_function = model_function,
                                                       dem_plus_case_data = dem_plus_case_data,
                                                       lineage_frequencies = lineage_frequencies,
                                                       intensity_scores = intensity_scores,
                                                       cutoff_age = cutoff_age,
                                                       maternal_ab_duration = maternal_ab_duration,
                                                       season_incidence_curves = season_incidence_curves,
                                                       school_start_age = school_start_age,
                                                       oldest_atk_rate_age = oldest_atk_rate_age,
                                                       reporting_age_cutoff = reporting_age_cutoff))
  }
  gc()
  return(loglik)
}
