# Fits specified model fixing specified parameters (for likelihood profiling)
library(dplyr)
library(optimParallel)
library(lubridate)

args = commandArgs(trailingOnly = T)
parameter_file = args[1] 

# Source file specifying data, model, selected parameters and subset region, minimum birth year to consider
source(parameter_file)

# Write start time of fitting
start_time = now()

# Import model functions
source('model_functions.R')
source('imprinting_probabilites.R')
source('merge_data.R')

# Load basic parameters (maternal ab duration, school start age)
source('basic_parameters.R')

# Read processed data files
case_data <- as_tibble(read.csv(case_data_path))


case_data <- set_min_birth_year(case_data, start_birth_year)

demographic_data <- as_tibble(read.csv(demographic_data_path))



demographic_data <- normalize_demographic_data(demographic_data, start_birth_year)

intensity_scores <- as_tibble(read.csv(intensity_scores_path))
lineage_frequencies <- as_tibble(read.csv(lineage_frequencies_path))
season_incidence_curves <- as_tibble(read.csv(season_incidence_curves_path))

# Combine demographic and case data
dem_plus_case_data <- merge_data(demographic_data, case_data)

if(!is.na(min_age)){
  min_age = as.integer(min_age)
  dem_plus_case_data <- impose_min_age(dem_plus_case_data, min_age)
}

# Separate string of model parameter names into character vector
selected_par_names <- strsplit(selected_par_names, split = ',')[[1]]

# Same with parameter values
selected_par_values <- as.numeric(strsplit(selected_par_values, split = ',')[[1]])
print(paste('BOUNDED PARS', bounded_par_names))

# Read tibble with bounds to sample initial parameters from
initial_par_bounds <- as_tibble(read.csv(initial_par_bounds_path))

# If pre-computed history probabilities are being used, read them
if(precomputed_history_probs_path == 'NA'){precomputed_history_probs_path = NA}
if(!is.na(precomputed_history_probs_path)){
  print(precomputed_history_probs_path)
  precomputed_history_probs = as_tibble(read.csv(precomputed_history_probs_path))
}else{
  precomputed_history_probs = NULL
}

# Initialize constrained parameter names and constrained parameter values based on profiled parameters
constrained_par_names <- selected_par_names
constrained_par_values <- selected_par_values

# Change parameter bounds (abs. bounds, not initial par bounds) if non-default bounds were provided
# (For a parameter with equal lower and upper bounds, simply add parameter to vector of constrained parameters)
bounded_par_names <- ifelse(bounded_par_names == 'NA',NA, bounded_par_names)
bounded_par_lbounds <- ifelse(bounded_par_lbounds == 'NA',NA, bounded_par_lbounds)
bounded_par_ubounds <- ifelse(bounded_par_ubounds == 'NA',NA, bounded_par_ubounds)

if(!is.na(bounded_par_names)){
  if(any(is.na(c(bounded_par_lbounds, bounded_par_ubounds)))){
    stop('Non-default parameter bounds not specified')
  }
  bounded_par_names <- strsplit(bounded_par_names, split = ',')[[1]]
  bounded_par_lbounds <- strsplit(bounded_par_lbounds, split = ',')[[1]]
  bounded_par_ubounds <- strsplit(bounded_par_ubounds, split = ',')[[1]]
  
  for(i in 1:length(bounded_par_names)){
    if(('default' %in% c(bounded_par_lbounds[i], bounded_par_ubounds[i]) == F) &
       bounded_par_lbounds[i] == bounded_par_ubounds[i]){
      constrained_par_names <- c(constrained_par_names, bounded_par_names[i])
      constrained_par_values <- c(constrained_par_values, as.numeric(bounded_par_lbounds[i]))
    }else{
      if(bounded_par_lbounds[i] != 'default'){
        parameter_bounds$lower_bound[parameter_bounds$parameter == bounded_par_names[i]] = 
          as.numeric(bounded_par_lbounds[i])
      }
      if(bounded_par_ubounds[i] != 'default'){
        parameter_bounds$upper_bound[parameter_bounds$parameter == bounded_par_names[i]] =
          as.numeric(bounded_par_ubounds[i])
      }
    }
  }
}

# Get model function and parameter names from model_functions script
model_function <- get(paste(selected_model_name, 'model', sep ='_'))
model_par_names <- get(paste(selected_model_name, '_model_par_names', sep = ''))

# Check that pars. of interested are listed in model parameter names
stopifnot(all(selected_par_names %in% model_par_names))
# Check specified # of parameter values = specified # of par. names
stopifnot(length(selected_par_values) == length(selected_par_names))

# Initialize cluster for optimParallel
cluster = makeForkCluster(nnodes = n_cores)
#cluster = makeCluster(nnodes = 9, type = 'FORK')

# If start birth year >= 1988, ancestral protection parameters are constrained to 0
if(start_birth_year >= 1988){
  constrained_par_names <- c(constrained_par_names, 'gamma_AV','gamma_AY')
  constrained_par_values <- c(constrained_par_values, 0, 0)
}

# Names of free parameters
free_parameter_names = model_par_names[model_par_names %in% constrained_par_names == F]

# Stop if constrained parameter values are outside of allowed bounds
value_outside_bound <- tibble(parameter = selected_par_names, value = selected_par_values,
                              lbound = lower_bounds(parameter, parameter_bounds),
                              ubound = upper_bounds(parameter, parameter_bounds)) %>%
  summarise(value_outside_bound = any(value < lbound) | any(value > ubound)) %>%
  pull(value_outside_bound)
stopifnot(value_outside_bound == F)

# Retrieve bounds for free parameters
lbounds <- lower_bounds(par_names = free_parameter_names, parameter_bounds)
ubounds <- upper_bounds(par_names = free_parameter_names, parameter_bounds)

# Fit model with parameters of interest fixed (and unused country-specific pars. fixed to NA if any)
convergence <- FALSE
while(convergence == FALSE){
  # Generate initial values for free parameters
  ipars <- init_pars(par_names = free_parameter_names, initial_par_bounds)
  
  model_fit <- optimParallel(par = ipars, fn = constrained_loglik_function, model_function = model_function,
                             model_par_names = model_par_names, selected_par_names = constrained_par_names,
                             selected_par_values = constrained_par_values, dem_plus_case_data = dem_plus_case_data,
                             lineage_frequencies = lineage_frequencies, intensity_scores = intensity_scores,
                             cutoff_age = cutoff_age, maternal_ab_duration = maternal_ab_duration, 
                             season_incidence_curves = season_incidence_curves, school_start_age = school_start_age,
                             oldest_atk_rate_age = oldest_atk_rate_age,
                             reporting_age_cutoff = reporting_age_cutoff,
                             precomputed_history_probs = precomputed_history_probs,
                             parameter_bounds = parameter_bounds,
                             #method = 'CG',
                             method = 'L-BFGS-B', 
                             lower = lbounds,
                             upper = ubounds,
                             control = list(maxit = 10000),
                             parallel = list(cl = cluster, forward = F))
  convergence <- model_fit$convergence
  convergence <- ifelse(convergence == 0, T, F)
}
loglik <- -model_fit$value
par_estimates <- model_fit$par
names(par_estimates) <- free_parameter_names

# Write results
results <- data.frame(rbind(c(selected_par_values, loglik, selected_model_name,
                              par_estimates,
                              constrained_par_values[(constrained_par_names %in% selected_par_names) == F]
                              )))
colnames(results) <- c(selected_par_names, 'loglik','model',
                       free_parameter_names,
                       constrained_par_names[(constrained_par_names %in% selected_par_names) == F])

write.csv(results,
          file = paste(c(output_directory, paste(selected_par_values, collapse = '_'), '.csv'),
                       collapse = ''),
          row.names = F)
# Write start and end time
end_time = now()
#time_file_conn <- file(paste(c(output_directory, paste(selected_par_values, collapse = '_'),
#                               '_fitting_time.csv'),
#                             collapse = ''))

#writeLines(c(paste('Time:', as.duration(start_time %--% end_time))),
#           time_file_conn)