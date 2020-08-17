#!/usr/bin/env Rscript
library(dplyr)
library(tidyr)
library(stringr)

args <- commandArgs(trailingOnly = TRUE)
main_directory <- args[1] # results_directory <- '../results/fitting_replicates/simulated_data/2018-06-11/'

update_reporting_factor_cols <- function(profile_tibble){
  # Renames reporting parameters for old runs of the single-rho (main) model, which used to be labelled according to country
  reporting_factor_cols <- colnames(profile_tibble)[grepl('reporting_factor',colnames(profile_tibble))]
  
  if(any(c("reporting_factor_aus","reporting_factor_us","reporting_factor_nz") %in% reporting_factor_cols)){
    # Retain column that's non-NA, corresponding to country the model was fitted to
    na_columns <- reporting_factor_cols[colSums(is.na(profile_tibble[reporting_factor_cols])) > 0]
    stopifnot(length(na_columns) == 2)
    retained_cols <- colnames(profile_tibble)[(colnames(profile_tibble) %in% na_columns) == F]
    profile_tibble <- profile_tibble[, retained_cols]
    # Rename reporting factor column as simply 'reporting_factor'
    names(profile_tibble)[grepl('reporting_factor', names(profile_tibble))] <- 'reporting_factor'
  }
  return(profile_tibble)
  
}



# Funciton for reading and combining results files within a single directory
combine_files <- function(results_directory){
  files <- list.files(results_directory)
  if(length(files) > 0){
    model_selection_files <- files[grepl('.csv',files) & 
                                     !grepl('fitting_time', files) &
                                     !grepl('likelihood_profile', files)]
    if(length(model_selection_files) > 0){
      file_paths <- paste(results_directory,model_selection_files, sep = '')
      
      result_tibbles_list <- lapply(file_paths, FUN = read.csv, header = T)
      result_tibbles_list <- lapply(result_tibbles_list, FUN = update_reporting_factor_cols)
    
      combined_tibble <- result_tibbles_list[[1]]
      if(length(result_tibbles_list) > 1){
        for(i in 2:length(result_tibbles_list)){
          combined_tibble <- bind_rows(combined_tibble, result_tibbles_list[[i]])
        }
      }
      
      write.csv(combined_tibble,
                paste(results_directory, 'likelihood_profile.csv', sep = ''),
                row.names = F)
    }
  }
}

combine_profiles <- function(profiles_directory){
  profile_paths <- list.dirs(profiles_directory)
  profile_paths <- str_replace(profile_paths[-1],'//','/')
  profile_paths <- paste0(profile_paths, '/likelihood_profile.csv')
  profile_paths <- profile_paths[file.exists(profile_paths)]
  
  profiles_list <- lapply(profile_paths, FUN = read.csv, header = T)
  
  combined_profiles <- profiles_list[[1]]
  if(length(profiles_list) > 1){
    for(i in 2:length(profiles_list)){
      combined_profiles <- bind_rows(combined_profiles, profiles_list[[i]])
    }
  }
  
  combined_profiles <- combined_profiles %>% select(model, loglik, everything())

  
  write.csv(combined_profiles,
            paste(profiles_directory, 'combined_likelihood_profiles.csv', sep = ''),
            row.names = F)
}

dir_list <- list.dirs(main_directory)
dir_list <- dir_list[-1]
dir_list <- str_replace(dir_list, '//','/')
dir_list <- paste0(dir_list, '/')
dir_list <- dir_list[grepl('likelihood_profiles', dir_list)]


sapply(dir_list, FUN = combine_files)
combine_profiles(main_directory)

# Remove csv files
#file.remove(file_paths)