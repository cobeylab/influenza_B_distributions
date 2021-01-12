source('synth_data_functions.R')
source('imprinting_probabilites.R')
source('merge_data.R')
source('basic_parameters.R')

args <- commandArgs(trailingOnly = T)

# path to file with true parameter values (output case data will be exported to the same directory)
true_pars_file_path <- args[1]
output_dir <- dirname(true_pars_file_path)

# By default, synth data will be based on these configurations
case_data_path='../results/processed_data/case_data_nz_all_surveillance_untyped_assigned.csv'
demographic_data_path='../results/processed_data/demographic_data.csv'
intensity_scores_path='../results/processed_data/intensity_scores.csv'
lineage_frequencies_path='../results/processed_data/lineage_frequencies_gisaid-genbank_noVicin1990s.csv'
season_incidence_curves_path='../results/processed_data/season_incidence_curves.csv'
start_birth_year=1959
subset_region='New_Zealand'
reporting_age_cutoff=1
model_name = 'main'
model_par_names <- get(paste(model_name, '_model_par_names', sep = ''))

# By default, will use the expectation, not a stochastic realization
stochastic <- F

# Read base data and covariates
case_data <- as_tibble(read.csv(case_data_path))
case_data <- set_min_birth_year(case_data, start_birth_year)

demographic_data <- as_tibble(read.csv(demographic_data_path))
demographic_data <- normalize_demographic_data(demographic_data, start_birth_year)

intensity_scores <- as_tibble(read.csv(intensity_scores_path))
lineage_frequencies <- as_tibble(read.csv(lineage_frequencies_path))
season_incidence_curves <- as_tibble(read.csv(season_incidence_curves_path))

# Combine demographic and case data
dem_plus_case_data <- merge_data(demographic_data, case_data)
generate_synthetic_data(model = 'main')

# Read true parameter values
true_pars <- as_tibble(read.csv(true_pars_file_path))
stopifnot(all(true_pars$par == model_par_names))
true_pars <- true_pars$value

synth_data <- generate_synthetic_data(model_name, model_parameters = true_pars, dem_plus_case_data, lineage_frequencies,
         intensity_scores, cutoff_age, maternal_ab_duration, season_incidence_curves,
         school_start_age, oldest_atk_rate_age, reporting_age_cutoff,
         precomputed_history_probs = NULL, stochastic)

synth_data <- synth_data %>% select(country, region, observation_year, cohort_type, 
                                    cohort_value, minimum_birth_year, lineage, n_cases)

write.csv(synth_data, file = paste0(output_dir,'/synth_case_data.csv'), row.names = F)


