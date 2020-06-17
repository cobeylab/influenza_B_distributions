library(dplyr)
library(tidyr)

# Data from the 2015 eLife paper
as_tibble(read.csv('../data/case_data/Vijaykrishna2015_Australian_cases.tsv', sep = '\t')) %>%
  rename(lineage = typesubtype) %>%
  mutate(age = floor(age)) %>%
  group_by(year, lineage, country, age) %>%
  count() %>% 
  ungroup() %>% filter(country == 'Australia') %>%
  mutate(region = 'AUSNZ', cohort_type = 'age') %>%
  mutate(minimum_birth_year = year - age - 1) %>%
  rename(cohort_value = age, observation_year = year, n_cases = n) %>%
  select(country, region, observation_year, cohort_type, cohort_value, minimum_birth_year, lineage, n_cases) %>%
  write.csv('../results/processed_data/Australian_case_data.csv')
