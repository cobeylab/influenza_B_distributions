library(dplyr)
# Basic parameters
cutoff_age <- 100
maternal_ab_duration <- 26
school_start_age <- tibble(country = c('United States','Australia','New Zealand','Europe','China'),
                           school_start_age = c(5,6,6,6,6))
oldest_atk_rate_age <- 18 # Age at which the oldest attack rate bin starts.

birth_year_cutoff <- 1970 # For birth years prior to this year, repeat this year's infection history probabilities
