# Functions for plotting model fits. Called by plot_model_fits.R and exploratory_prediction_plots.R
# Function for generating stochastic multinomial realizations of model predictions
source('generate_synthetic_data.R')

# Function for adjusting predictions so all cohorts are plotted by specified cohort type (age or b. year)
set_data_to_cohort_type <- function(predictions, plot_cohort_type){
  # If plot cohort type is birth year, convert age to minimum possible birth year
  if(plot_cohort_type == 'birth_year'){
    predictions <- predictions %>% rowwise() %>%
      mutate(
        cohort_value = ifelse(cohort_type == 'age',
                              observation_year - cohort_value - 1,
                              cohort_value)
      ) %>% ungroup()
  }
  # If plot cohort type is age, convert birth year to maximum possible age
  if(plot_cohort_type == 'age'){
    predictions <- predictions %>% rowwise() %>%
      mutate(
        cohort_value = ifelse(cohort_type == 'birth_year',
                              observation_year - cohort_value,
                              cohort_value)
      ) %>% ungroup()
  }
  return(predictions)
}

# Function to add binomial quantiles to the predicted number of cases
add_binomial_quantiles <- function(predictions, alpha = 0.05){
  predictions <- predictions %>% 
    mutate(imprinting_predicted_cases_lower = qbinom(p = alpha/2, size = CLY_total_cases, prob = pred_case_prob),
           imprinting_predicted_cases_upper = qbinom(p = 1 - alpha/2, size = CLY_total_cases, prob = pred_case_prob))
  return(predictions)
}

# Function for computing the fraction of cases of a lineage by birth year
compute_case_fractions <- function(cases){
  # This automatically determines the level of aggregation 
  grouping_variables <- names(cases)[names(cases) %in% c('lineage','country','observation_year','interval',
                                                         'SurveillanceType')]
  stopifnot('lineage' %in% grouping_variables)
  
  cases <- cases %>% 
    group_by_at(vars(one_of(grouping_variables)))
  
  cases <- cases %>%
    mutate(total_cases = sum(obs_cases),
           fraction_obs_cases = obs_cases / total_cases,
           demography_predicted_fraction = demography_predicted_cases / total_cases,
           imprinting_predicted_fraction = imprinting_predicted_cases / total_cases) %>%
    ungroup()
  
  return(cases)
}

get_MLE_parameters <- function(combined_profile, constrained_par_name = NA, constrained_par_value = NA){
  model <- as.character(unique(combined_profile$model))
  stopifnot(length(model) == 1)
  model_par_names <- get(paste(model, '_model_par_names', sep = '')) 
  # Parameters that were effectively estimated (i.e., no ancestral imprinting if post 1988, or US if focusing on AUS/NZ)
  inferred_par_names <- model_par_names[model_par_names %in% colnames(combined_profile)]
  
  if(!all(is.na(constrained_par_name))){
    stopifnot(!is.na(constrained_par_value))
    for(i in 1:length(constrained_par_name)){
      combined_profile <- combined_profile %>%
        filter_(paste0(constrained_par_name[i],'==',constrained_par_value[i]))
    }
    if(nrow(combined_profile) == 0){
      stop('Constrained parameter values not available in profile')}
  }
  
  mle_pars <- combined_profile %>% filter(loglik == max(loglik))
  if(nrow(mle_pars) > 1){
    mle_pars <- mle_pars[1,]
  }
  mle_pars <- mle_pars %>% select(inferred_par_names) %>% 
    unlist(., use.names=T)
  
  # Fill in irrelavant parameters (e.g. ancestral imprinting if post 1988, NZ if focusing on AUS)
  mle_pars <- left_join(tibble(par = model_par_names), 
                        melt(mle_pars) %>% mutate(par = rownames(melt(mle_pars))), by = 'par')
  
  return(mle_pars$value)
}

get_hard_constrained_pars <- function(mle_pars, constrained_par_name, constrained_par_value, model_par_names){
  constrained_mle_pars <- mle_pars
  for(i in 1:length(constrained_par_name)){
    constrained_mle_pars[model_par_names == constrained_par_name[i]] <- constrained_par_value[i]
  }
  return(constrained_mle_pars)
}


# Function for aggregating predictions across observation years (by interval or across all years)
# Will aggregate across countries if any are present.
# Function for by country plots will pass country-specific predictions
aggregate_predictions <- function(predictions, grouping_variables = NULL, interval_breaks = NULL){
  # Interval breaks: tibble with country and year of interval break (1st y of 2nd interval)
  if('interval' %in% grouping_variables & is.null(interval_breaks)){
    stop('Breaks must be specified if aggregating within time intervals')
  }
  if('interval' %in% grouping_variables & ('country' %in% grouping_variables == F)){
    stop('Aggregating within time intervals while pooling across countries not implemented')
  }
  if('observation_year' %in% grouping_variables & ('country' %in% grouping_variables == F)){
    stop('Aggregating across countries by observation year not implemented.')
  }
  stopifnot(('interval' %in% grouping_variables & 'observation_year' %in% grouping_variables) == F)
  
  grouping_variables <- c('lineage','cohort_value', grouping_variables)
  
  if('interval' %in% grouping_variables){
    # For aggregating within time intervals
    range_obs_years <- predictions %>% group_by(country) %>% 
      summarise(min_obs_year = min(observation_year),
                max_obs_year = max(observation_year))
    
    agg_predictions <- left_join(predictions, interval_breaks, by = 'country')
    agg_predictions <- left_join(agg_predictions, range_obs_years, by = 'country')
    
    agg_predictions <- agg_predictions %>%
      # Divide data into intervals
      mutate(interval = ifelse(observation_year < interval_break,
                               paste0(min_obs_year, '-', interval_break -1),
                               paste0(interval_break, '-', max_obs_year))) 
  }else{
    agg_predictions <- predictions 
  }
  agg_predictions <- agg_predictions %>%
    group_by_at(vars(one_of(grouping_variables))) %>%
    summarise_at(c('obs_cases','demography_predicted_cases','imprinting_predicted_cases'),
                 sum) %>%
    ungroup()
  
  return(agg_predictions)
}

add_boostrap_CI <- function(predictions, n_replicates = 1000, alpha = 0.05, grouping_variables = NULL,
                            interval_breaks = NULL){
  
  # Generate first realization and aggregate it at chosen level
  realizations <- generate_stochastic_realization(predictions)
  
  # If grouping_variables inc. 'observation_year', no aggreg. is done (and must include 'country')
  realizations <- aggregate_predictions(realizations, grouping_variables, interval_breaks)  
  realizations <- compute_case_fractions(realizations) %>% 
    mutate(replicate = 1)
  
  # Generate replicate realizations
  for(i in 2:n_replicates){
    new_realization <- generate_stochastic_realization(predictions)
    new_realization <- aggregate_predictions(new_realization,grouping_variables, interval_breaks) 
    new_realization <- compute_case_fractions(new_realization) %>%
      mutate(replicate = i)
    realizations <- bind_rows(realizations,new_realization)
  }
  
  bootstrap_CIs <- realizations %>% group_by_at(vars(one_of(c(grouping_variables,'cohort_value','lineage')))) %>% 
    summarise(
      imprinting_predicted_cases_lower = quantile(imprinting_predicted_cases, alpha/2),
      imprinting_predicted_cases_upper = quantile(imprinting_predicted_cases, 1 - alpha/2),
      imprinting_predicted_fraction_lower = quantile(imprinting_predicted_fraction, alpha/2),
      imprinting_predicted_fraction_upper = quantile(imprinting_predicted_fraction, 1 - alpha/2)
    ) %>% ungroup()
  
  # Return aggregated dataalong with bootstrap C.Is
  predictions <- aggregate_predictions(predictions, grouping_variables, interval_breaks)
  
  # Calculate case fractions for observed data
  predictions <- compute_case_fractions(predictions)
  
  return(left_join(predictions,bootstrap_CIs))
}

base_plot_function <- function(data, plot_predictions, plot_excess_cases = F, country = NULL,
                               line_size = 1, point_size = 2){
  #data: list of tibbles with observed cases and associated predictions
  if(!is.null(country)){
    for(i in 1:length(data)){
      data[[i]] <-  data[[i]] %>% filter(country == !!country)  
    }
  }
  
  if(plot_excess_cases){
    stopifnot(plot_predictions == F)
    if(length(data) == 1){
      pl <- data[[1]] %>% ggplot(aes(x = cohort_value, y = predicted_value - observed_value)) +
        geom_col() +
        geom_linerange(aes(ymin = predicted_value_lower - observed_value,
                           ymax = predicted_value_upper - observed_value), size = 1, col = 'red',
                       alpha = 0.4)
    }else{
      stopifnot(length(data) ==2)
      pl <-  data[[1]] %>% ggplot(aes(x = cohort_value, y = predicted_value - observed_value)) +
        geom_point(colour = 'red', shape = 21, size = 0.8) +
        geom_line(color = 'red') +
        geom_linerange(aes(ymin = predicted_value_lower - observed_value,
                           ymax = predicted_value_upper - observed_value), size = 0.5, col = 'red',
                       alpha = 0.4) +
        geom_point(data = data[[2]], color = 'blue', size = 0.5, alpha = 0.8)
    }
    
  }else{
    pl <- data[[1]] %>% 
      ggplot(aes(x = cohort_value, y = observed_value)) + 
      geom_col()
    
    if(plot_predictions){
      stopifnot(plot_excess_cases == F)
      pl <- pl + geom_linerange(aes(ymin = predicted_value_lower,
                                    ymax = predicted_value_upper), size = 1, col = 'red',
                                alpha = 0.4) +
        geom_line(aes(y=predicted_value), size = line_size, col = 'red') +
        geom_point(aes(y=predicted_value), size = point_size, col = 'red')
      
      if(length(data) > 1){
        stopifnot(length(data) ==2)
        pl <- pl + geom_linerange(data = data[[2]],
                                  aes(ymin = predicted_value_lower,
                                      ymax = predicted_value_upper), size = 1, col = 'blue',
                                  alpha = 0.4) +
          geom_line(data = data[[2]], aes(y=predicted_value), size = line_size, col = 'blue', alpha = 0.5) +
          geom_point(data = data[[2]], aes(y=predicted_value), size = point_size, col = 'blue', alpha = 0.5)
        
        if(!is.null(names(data))){
          pl <- pl +
            # Add empty geometry just to add color guide
            geom_line(data = tibble(cohort_value = c(1980,1981), predicted_value = c(0.025,0.025),
                                    prediction_type = factor(names(data), levels = 
                                                               names(data)),
                                    lineage = c('B/Victoria','B/Yamagata')),
                      aes(y = predicted_value, x = cohort_value, color = prediction_type)) +
            scale_color_manual(name = '', values = c('red','blue'),
                               labels = names(data)) +
            theme(legend.position = 'top')
        }
      }
    }
  }
  return(pl)
}

normalize_predictions <- function(predictions,plot_fraction, demographic_normalization){
  predictions <- predictions %>% rowwise() %>%
    mutate(observed_value = ifelse(plot_fraction, fraction_obs_cases, obs_cases),
           predicted_value = ifelse(plot_fraction, imprinting_predicted_fraction,
                                    imprinting_predicted_cases),
           predicted_value_lower = ifelse(plot_fraction,imprinting_predicted_fraction_lower,
                                          imprinting_predicted_cases_lower),
           predicted_value_upper = ifelse(plot_fraction,imprinting_predicted_fraction_upper,
                                          imprinting_predicted_cases_upper)) %>%
    ungroup()
  
  if(demographic_normalization){
    predictions <- predictions %>%
      rowwise() %>%
      # Calculate observed and predicted deviations from demographic expectation
      mutate(
        observed_value = observed_value / ifelse(plot_fraction,
                                                 demography_predicted_fraction,
                                                 demography_predicted_cases),
        predicted_value = predicted_value / ifelse(plot_fraction,
                                                   demography_predicted_fraction,
                                                   demography_predicted_cases),
        predicted_value_lower = predicted_value_lower / ifelse(plot_fraction,
                                                               demography_predicted_fraction,
                                                               demography_predicted_cases),
        predicted_value_upper = predicted_value_upper / ifelse(plot_fraction,
                                                               demography_predicted_fraction,
                                                               demography_predicted_cases)
        
      ) %>% ungroup()
  }
  return(predictions)
}

# Function for making plots for a specified country
plot_by_obs_year <- function(predictions, n_CI_replicates, CI_alpha, plot_predictions = T,
                             plot_excess_cases = F, plot_fraction = F, demographic_normalization = F){
  
  grouping_variables <- c('observation_year','country')
  
  if(plot_fraction){
    if(plot_excess_cases){
      ylabel <- 'Predicted minus observed cases'
    }else{
      ylabel <- 'Fraction of cases'
    }
  }else{
    if(plot_excess_cases){
      ylabel <- 'Predicted minus observed cases'
    }else{
      ylabel <- 'Number of cases'
    }
  }
  if(demographic_normalization){
    ylabel <- paste0(ylabel, ' relative to\ndemographic expectation')
  }
  
  if(plot_predictions == T | plot_excess_cases == T){
    # Remove country/lineage/year combinations with 5 cases or less from plots with fits
    for(i in 1:length(predictions)){
      predictions[[i]] <- predictions[[i]] %>% filter(CLY_total_cases >= 5)
    }
  }
  
  for(i in 1:length(predictions)){
    # Adjust predictions so all cohorts plotted by birth year
    predictions[[i]] <- set_data_to_cohort_type(predictions[[i]], plot_cohort_type = 'birth_year')
    
    # Add bootstrap CIs.
    predictions[[i]] <- add_boostrap_CI( predictions[[i]], n_CI_replicates, CI_alpha, grouping_variables,
                                         interval_breaks = NULL)
    # Normalize
    predictions[[i]] <- normalize_predictions(predictions[[i]], plot_fraction, demographic_normalization)
  }  
  
  countries_in_data <- unique(predictions[[1]]$country)
  
  plot_list <- lapply(countries_in_data,FUN = base_plot_function,
                      data = predictions, plot_predictions = plot_predictions,
                      plot_excess_cases = plot_excess_cases,
                      line_size = 0.5, point_size = 1)
  for(i in 1:length(plot_list)){
    plot_list[[i]] <- plot_list[[i]] + 
      facet_grid(observation_year~lineage, scales = 'free_y') +
      ggtitle(label = countries_in_data[i]) +
      xlab('Year of birth') +
      ylab(ylabel) +
      background_grid(major = "xy", minor = "none") +
      theme(axis.text.x = element_text(size = 12))
  }
  names(plot_list) <- countries_in_data
  return(plot_list)
}

plot_by_country <- function(predictions, n_CI_replicates, CI_alpha, plot_predictions = T,
                            plot_fraction = F, demographic_normalization = F){
  
  grouping_variables <- c('country')
  for(i in 1:length(predictions)){
    # Adjust predictions so all cohorts plotted by birth year
    predictions[[i]] <- set_data_to_cohort_type(predictions[[i]], plot_cohort_type = 'birth_year')
    
    # Add bootstrap CIs.
    predictions[[i]] <- add_boostrap_CI(predictions[[i]], n_CI_replicates, CI_alpha, grouping_variables,
                                        interval_breaks = NULL)
    
    # Normalize
    predictions[[i]] <- normalize_predictions(predictions[[i]], plot_fraction, demographic_normalization)
  }
  
  
  pl <- base_plot_function(predictions, plot_predictions, country = NULL)
  
  pl <- pl + facet_grid(country~lineage, scales = 'free_y') + 
    xlab('Year of birth') +
    ylab(paste0(ifelse(plot_fraction,'Fraction','Number'),' of cases',
                ifelse(demographic_normalization, ' relative\nto demographic expectation',''))) +
    background_grid(major = "xy", minor = "none")
  return(pl)
}

plot_by_interval <- function(predictions, n_CI_replicates, CI_alpha, plot_predictions = T,
                             plot_fraction = F, demographic_normalization = F){
  # Hard-coded interval breaks (1st year of 2nd interval)
  interval_breaks <- tibble(country = c('Australia','New Zealand','United States'),
                            interval_break = c(2008,2008,2013))
  
  pl <- predictions 
  # Adjust predictions so all cohorts plotted by birth year
  pl <- set_data_to_cohort_type(pl, plot_cohort_type = 'birth_year')
  
  grouping_variables <- c('country','interval')
  
  # Add bootstrap CIs.
  pl <- add_boostrap_CI(pl, n_CI_replicates, CI_alpha, grouping_variables, interval_breaks)
  
  # Normalize
  pl <- normalize_predictions(pl, plot_fraction, demographic_normalization)
  countries_in_data <- unique(pl$country)
  
  plot_list <- lapply(countries_in_data,FUN = base_plot_function,
                      data = pl, plot_predictions = plot_predictions)
  for(i in 1:length(plot_list)){
    plot_list[[i]] <- plot_list[[i]] + 
      facet_grid(lineage ~ interval, scales = 'free_y') +
      ggtitle(label = countries_in_data[i]) +
      xlab('Year of birth') +
      ylab(paste0(ifelse(plot_fraction,'Fraction','Number'),' of cases',
                  ifelse(demographic_normalization, ' relative\nto demographic expectation',''))) +
      background_grid(major = "xy", minor = "none")
  }
  names(plot_list) <- countries_in_data
  return(plot_list)
}

# Function for combining Aus, NZ and US results (but not Canada, bc. of age groups)
plot_pooling_countries <- function(predictions, n_CI_replicates, CI_alpha, plot_predictions = T,
                                   plot_fraction = F, demographic_normalization = F,
                                   pool_surveillance_types = T){
  
  grouping_variables <- NULL

  if(pool_surveillance_types == F){
    stopifnot('SurveillanceType' %in% names(predictions[[1]]))
    grouping_variables <- 'SurveillanceType'
  }

  
  for(i in 1:length(predictions)){
    # Adjust predictions so all cohorts plotted by birth year
    predictions[[i]] <- set_data_to_cohort_type(predictions[[i]], plot_cohort_type = 'birth_year')
    
    # Add bootstrap CIs.
    predictions[[i]] <- add_boostrap_CI(predictions[[i]], n_CI_replicates, CI_alpha, grouping_variables,
                                        interval_breaks = NULL)
    # Normalize
    predictions[[i]] <- normalize_predictions(predictions[[i]], plot_fraction, demographic_normalization)
    
  }
  
  pl <- base_plot_function(predictions, plot_predictions, country = NULL)
  
  if(pool_surveillance_types == F){
    pl <- pl + facet_grid(SurveillanceType~lineage, scales = 'free_y')
  }else{
    pl <- pl + facet_grid(.~lineage, scales = 'free_y')
  }
  
  pl <- pl + 
    xlab('Year of birth') +
    ylab(paste0(ifelse(plot_fraction,'Fraction','Number'),' of cases',
                ifelse(demographic_normalization, ' relative\nto demographic expectation',''))) +
    background_grid(major = "xy", minor = "none")
  return(pl)
}

# Function for plotting parameter estimates for an individual model
plot_model_pars <- function(model, parameters_summary, combined_model_selection_results){
  # Check if model has true values listed in parameters_summary (synthetic data)
  has_true_values <- filter(parameters_summary, model_name == model) %>%
    summarise(has_true_values = 'true_value' %in% statistic) %>% 
    pull(has_true_values)
  
  model_par_names <- get(paste(model, '_model_par_names', sep = ''))  
  
  # par_estimates_pl <- combined_model_selection_results %>% filter(model_name == model) %>%
  #   select(model_par_names) %>%
  #   melt() %>%
  #   rename(parameter = variable) %>%
  #   ggplot(aes(x = value, y = parameter)) +
  #   geom_density_ridges2(stat = 'binline')
  
  par_estimates_pl <- ggplot(data = filter(parameters_summary,
                                           statistic %in% c('replicate_min','replicate_max'),
                                           model_name == model) %>%
                               spread(key = statistic, value = value)) + 
    geom_linerange(aes(x = par_name, ymin = replicate_min, ymax = replicate_max),
                   color = "#7570b3") +
    ylim(-0.1,2) +
    coord_flip() + 
    xlab('Parameter') + ylab('Value') +
    theme(axis.text.y = element_text(size = 9), legend.position = c(0.65,0.95))
  
  if(has_true_values){
    par_estimates_pl <- par_estimates_pl +
      geom_point(data = filter(parameters_summary, statistic %in% c('true_value', 'replicate_mean'),
                               !is.na(value), model_name == model),
                 aes(x = par_name, y = value,
                     color = factor(statistic)),
                 size = 2.5, alpha = 0.6) +
      scale_color_manual(values = c("#7570b3","#fc8d59"), name = '',
                         labels = c('Mean MLE value', 'True value'))
    
  }else{
    par_estimates_pl <- par_estimates_pl +
      geom_point(data = filter(parameters_summary, statistic %in% c('replicate_mean'), !is.na(value),
                               model_name == model),
                 aes(x = par_name, y = value),
                 color = "#7570b3", size = 2.5, alpha = 0.6) +
      scale_color_manual(values = c("#7570b3","#7570b3"), name = '')
  }
  
  save_plot(paste(plot_directory, model, '_parameters.pdf', sep = ''), par_estimates_pl,
            base_height = 5, base_width = 8)
  
}

plot_age_obs_vs_time_predictions <- function(predictions){
  predictions %>% group_by(lineage, observation_year) %>%
    summarise(mean_age = sum(cohort_value*imprinting_predicted_cases/sum(imprinting_predicted_cases))) %>%
    ggplot(aes(x = observation_year, y = mean_age)) +
    geom_point() +
    facet_grid(lineage~.)
  
  predictions %>% group_by(lineage, observation_year) %>%
    mutate(min_birth_year = observation_year - cohort_value - 1) %>%
    filter(min_birth_year <=1995) %>%
    mutate(pred_case_prob = pred_case_prob/sum(pred_case_prob)) %>%
    summarise(mean_byear = sum(min_birth_year*pred_case_prob)) %>%
    ggplot(aes(x = observation_year, y = mean_byear)) +
    geom_point() +
    facet_grid(lineage~.) + 
    geom_line()
  
  predictions %>% 
    mutate(birth_year = observation_year - cohort_value - 1) %>%
    filter(CLY_total_cases >= 50) %>%
    mutate(observation_year = factor(as.character(observation_year),
                                     levels = as.character(sort(unique(observation_year))))) %>%
    ggplot(aes(x = birth_year, y = obs_cases/CLY_total_cases, color = observation_year,
               group = observation_year)) +
    facet_grid(lineage~country) +
    geom_line() +
    theme_dark() +
    scale_color_brewer(palette = 3, name = 'Observation year') +
    theme_dark() +
    xlab('Birth year') +
    ylab('Observed fraction of cases')
  
  predictions %>% 
    mutate(birth_year = observation_year - cohort_value - 1) %>%
    filter(birth_year <= 2002) %>%
    filter(observation_year <=2010) %>%
    group_by(country, lineage, observation_year) %>%
    mutate(pred_case_prob = pred_case_prob/sum(pred_case_prob)) %>%
    ungroup() %>%
    mutate(observation_year = factor(as.character(observation_year),
                                     levels = as.character(sort(unique(observation_year))))) %>%
    ggplot(aes(x = birth_year, y = pred_case_prob, color = observation_year,
               group = observation_year)) +
    facet_grid(lineage~country) +
    geom_line() +
    scale_color_brewer(palette = 3, name = 'Observation year') +
    theme_dark() +
    xlab('Birth year') +
    ylab('Predicted fraction of cases')
  
  predictions %>% 
    mutate(birth_year = observation_year - cohort_value - 1) %>%
    filter(birth_year <= 2002) %>%
    filter(observation_year >=2002) %>%
    group_by(country, lineage, observation_year) %>%
    mutate(pred_case_prob = pred_case_prob/sum(pred_case_prob)) %>%
    summarise(fraction_cases_in_90s = sum(pred_case_prob[birth_year <= 1990]))
}
