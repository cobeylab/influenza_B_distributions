#!/usr/bin/env Rscript

library(dplyr)

# Function for merging demographic and case data
merge_data <- function(demographic_data, case_data){
  output_vars <- c("country", "region", "observation_year", "lineage","cohort_type",
                   "cohort_value", "rel_pop_size", "n_cases", "CLY_total_cases")
  multinomial_draw_unit <- c('country', 'lineage', 'observation_year')
  full_join_vars <- c("country","region","observation_year","cohort_type","cohort_value","lineage")
  
  # Duplicate demography tibble to match Victoria and Yamagata rows in merged tibble
  double_dem_tibble <- bind_rows(mutate(demographic_data, lineage = 'B/Victoria'),
                                 mutate(demographic_data, lineage = 'B/Yamagata'))
  
  # Filter demography tibble to retain only obs. years when there were cases
  double_dem_tibble <- filter(double_dem_tibble, observation_year %in% unique(case_data$observation_year)) %>%
    filter()
  
  # Merge demography and case tibbles
  # Doing a full join with the adjusted demog. tibble ensures b. years with no cases are represented as 0s
  dem_plus_case_data <- full_join(case_data, double_dem_tibble,
                                  by = full_join_vars)
  # Replace NAs with 0s for birth years without any cases %>%
  dem_plus_case_data$n_cases[is.na(dem_plus_case_data$n_cases)] <- 0
  
  dem_plus_case_data <- dem_plus_case_data %>%
    # Add total number of cases per country/lineage/observation year (CLY_total_cases)
    group_by_at(vars(multinomial_draw_unit)) %>%
    mutate(CLY_total_cases = sum(n_cases)) %>% 
    ungroup() %>%
    # Remove lineage/country/obs. year combinations with no cases
    filter(CLY_total_cases > 0) %>%
    rename(rel_pop_size = fraction) %>%
    select(output_vars) %>%
    arrange(observation_year, lineage, cohort_value)

  return(dem_plus_case_data)
}









