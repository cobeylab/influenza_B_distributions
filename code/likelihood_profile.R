# Constructs likelihood profiles for the specified model and parameters
# Based on Ed Baskerville's code for runmany
library(stringr)

args <- commandArgs(trailingOnly = T)
args[args == 'NA'] <- NA

case_data_path = args[1] # case_data_path = '../results/processed_data/case_data_nz_all_surveillance_untyped_assigned.csv'
demographic_data_path = args[2] # demographic_data_path = '../results/processed_data/demographic_data.csv'
intensity_scores_path = args[3] # intensity_scores_path = '../results/processed_data/intensity_scores.csv'
lineage_frequencies_path = args[4] # lineage_frequencies_path = '../results/processed_data/lineage_frequencies_gisaid-genbank_noVicin1990s.csv'
season_incidence_curves_path = args[5] # season_incidence_curves_path = '../results/processed_data/season_incidence_curves.csv'
start_birth_year = as.numeric(args[6]) # Earliest birth year in case data to analyze
subset_region = args[7] # New_Zealand or Australia
reporting_age_cutoff = as.numeric(args[8]) # Differential reporting if age <= reporting_age_cutoff
precomputed_history_probs_path = args[9] # set to NA for all analyses in the paper.

selected_model_name = args[10] # set to 'main' for all analyses in the paper.
selected_par_names = args[11] # Comma-separated list of selected parameter names, e.g. "R_Y,gamma_AY" or selected_par_names"R_Y"

# Parameter names
# R_V
# R_Y
# chi_V
# chi_Y
# gamma_VY
# gamma_YV
# gamma_AY
# gamma_AV
# beta1
# beta2
# beta3
# reporting_factor

lower_limits = args[12] # Comma-separated list of lower bounds for parameter values 
upper_limits = args[13] # Comma-separated list of upper bounds for parameter values 
increments = args[14] # Comma-separated values for the increment between par1 and par2 values
output_directory <- args[15] # Output for the likelihood profile over the chosen parameter(s)
n_cores <- as.integer(args[16]) # Number of cores to use. 
initial_par_bounds_path = as.character(args[17]) # csv with intervals to sample initial values. Default in the paper: code/intial_parameters/three_class_atk_rate_model_noVicin1990s.csv
bounded_par_names <- as.character(args[18]) # Comma-separated list of pars with non-default bounds. Can be left as NA
bounded_par_lbounds <- as.character(args[19]) # Comma-separated list of lower non-default bounds. Can be left as NA
bounded_par_ubounds <- as.character(args[20]) # Comma-separated list of upper non-default bounds. Can be left as NA
min_age <- as.integer(args[21]) # Excludes all ages < min_age

if(!is.na(precomputed_history_probs_path) & is.na(bounded_par_names)){
  stop('Must constrain betas if using pre-computed infection history probabilities')
}

# Separate strings of parameter bounds into numeric vector
lower_limits <- as.numeric(strsplit(lower_limits, split = ',')[[1]])
upper_limits <- as.numeric(strsplit(upper_limits, split = ',')[[1]])
increments <- as.numeric(strsplit(increments, split = ',')[[1]])

# Function to generate dataframe with all parameter value combinations
generate_par_combinations <- function(lower_limits, upper_limits, increments, selected_par_names){
  stopifnot(length(lower_limits) == length(upper_limits))
  # Range of parameter values
  values_list <- list()
  for(i in 1:length(lower_limits)){
    values_list[[i]] <- seq(lower_limits[i], upper_limits[i], increments[i])
  }
  values_dataframe <- expand.grid(values_list)
  colnames(values_dataframe) <- strsplit(selected_par_names, ',')[[1]]
  return(values_dataframe)
}

n_parameter_sets <- nrow(generate_par_combinations(lower_limits, upper_limits, increments,
                                                   selected_par_names))

# Gen. list of files with individual combinations of par. values (to pass to to fixed_pars_likelihood.R)
generate_par_files <- function(lower_limits, upper_limits, increments, selected_par_names,
                               selected_model_name){

  values_dataframe <- generate_par_combinations(lower_limits, upper_limits, increments,
                                                selected_par_names)
  file_counter <- 0
  for(i in 1:nrow(values_dataframe)){
    row <- values_dataframe[i,] 
    file_counter <- file_counter + 1
    file_path <- paste(output_directory, 'parameter_set_', file_counter, '.R', sep = '')
    file_conn <- file(file_path)
    
    writeLines(c(paste0("case_data_path = '", case_data_path, "'"),
                 paste0("demographic_data_path = '", demographic_data_path, "'"),
                 paste0("intensity_scores_path = '", intensity_scores_path, "'"),
                 paste0("lineage_frequencies_path = '", lineage_frequencies_path, "'"),
                 paste0("season_incidence_curves_path = '", season_incidence_curves_path, "'"),
                 paste0("subset_region = '", subset_region, "'"),
                 paste0("reporting_age_cutoff = ", reporting_age_cutoff),
                 paste0("precomputed_history_probs_path = '", precomputed_history_probs_path,"'"),
                 paste0("start_birth_year = ", start_birth_year),
                 paste0("selected_model_name = '", selected_model_name,"'"),
                 paste0("selected_par_names = '", selected_par_names, "'"),
                 paste(c("selected_par_values = '", paste(values_dataframe[i,], collapse = ','), "'"),
                       collapse = ''),
                 paste0("output_directory = '", output_directory, "'"),
                 paste0("n_cores = ", n_cores, sep = ''),
                 paste0("initial_par_bounds_path = '", initial_par_bounds_path,"'"),
                 paste0("bounded_par_names = '", bounded_par_names,"'"),
                 paste0("bounded_par_lbounds = '", bounded_par_lbounds, "'"),
                 paste0("bounded_par_ubounds = '", bounded_par_ubounds, "'"),
                 paste0("min_age = ", min_age)
                  ),
                file_conn)
    close(file_conn)
  }
  
}

generate_par_files(lower_limits, upper_limits, increments, selected_par_names, selected_model_name)

# Generate sbatch files that call fixed_pars_likelihood for individual parameter_set files
sbatch_file_name <- paste0('likprof_', selected_model_name, '_model_',
                           str_replace(selected_par_names,',','_VS_'),'.sbatch')

generate_sbatch_file <- function(sbatch_file_name, output_directory, parameter_combinations){
  n_parameter_sets <- nrow(parameter_combinations)
  
  # csv files in output directory
  csv_files = list.files(output_directory, pattern = '.csv')
  
  # Generate comma-separated list of numbers corresponding to parameter sets
  array <- ''
  for(i in 1:n_parameter_sets){
    # Check if a csv file identified by parameter values already exists in output directory
    results_file_name <- paste0(paste(parameter_combinations[i,], collapse = '_'), '.csv')
    if(results_file_name %in% csv_files){
      # If so, do not submit job for that parameter combination
      warning(paste(c("Results already present for" , parameter_combinations[i,], '- Skipping'), collapse = ' '))
      file.remove(paste0(output_directory, 'parameter_set_',i,'.R'))
    }else{
      array = paste(array, i, sep = ',')  
    }
  }
  # Remove extra comma at the beginning
  array <- substring(array, 2)
  
  file_conn <- file(sbatch_file_name)
  writeLines(c("#!/bin/bash",
               paste0("#SBATCH --job-name=LikPro_", paste(selected_par_names, collapse = '_vs_')),
               paste0("#SBATCH --output=sbatch_files/out_err_files/LikPro_%A_%a_", paste(selected_par_names, collapse = '_vs_'), '.out'),
               paste0("#SBATCH --error=sbatch_files/out_err_files/LikPro_%A_%a_", paste(selected_par_names, collapse = '_vs_'), '.err'),
               paste("#SBATCH --array=", array, sep = ''),
               "#SBATCH --time=60:00:00",
               "#SBATCH --partition=cobey",
               "#SBATCH --nodes=1",
               "#SBATCH --exclude=midway2-bigmem05",
               paste("#SBATCH --ntasks=", n_cores, sep = ''),
               "#SBATCH --mem-per-cpu=1000",
               "module load R/3.4.3",
               paste("parameter_set_file=", output_directory,
                     "parameter_set_${SLURM_ARRAY_TASK_ID}.R", sep = ''),
               "Rscript fixed_pars_likelihood.R $parameter_set_file",
               "rm $parameter_set_file"),
             file_conn)
  close(file_conn)
}


parameter_combinations <- generate_par_combinations(lower_limits, upper_limits, increments,
                                                    selected_par_names)
generate_sbatch_file(sbatch_file_name, output_directory, parameter_combinations)

# Run sbatch file
system(paste0('sbatch ', sbatch_file_name))

# Remove sbatch file
file.remove(sbatch_file_name)


