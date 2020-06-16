library(dplyr)
library(ggplot2)
library(cowplot)
theme_set(theme_cowplot())
library(stringr)
library(tidyr)
library(reshape2)

# New Zealand Data from ESR in 2012-2019, with type of surveillance information.
esr_nz_data <- as_tibble(read.csv('../data/case_data/nz_data_2012-2019_batch.csv', header = T)) %>%
  rename(age = age_in_years, observation_year = year, lineage = influenza_virus) %>%
  mutate(lineage = as.character(lineage))

esr_nz_data <- esr_nz_data %>% filter(!is.na(age))

unique(esr_nz_data$lineage)
esr_nz_data$lineage[esr_nz_data$lineage == 'B Yam'] <- 'B/Yamagata'
esr_nz_data$lineage[esr_nz_data$lineage == 'B Vic'] <- 'B/Victoria'
esr_nz_data$lineage[esr_nz_data$lineage == 'B'] <- 'unknown'

nz_data <- esr_nz_data %>%
  group_by(observation_year, age, lineage, severity) %>%
  summarise(n_cases = n()) %>%
  ungroup() %>%
  mutate(region = 'AUSNZ', country = 'New Zealand', cohort_type = 'age') %>%
  mutate(minimum_birth_year = observation_year - age - 1) %>%
  rename(cohort_value = age) %>%
  select(severity, country, region, observation_year, cohort_type, cohort_value, minimum_birth_year, lineage, n_cases)

# Check if cases with unknown lineage can be assigned when a lineage is clearly dominant
nz_data %>% 
  group_by(observation_year) %>%
  summarise(n_vic = sum(n_cases[lineage == 'B/Victoria']),
            n_yam = sum(n_cases[lineage == 'B/Yamagata']),
            n_unknown = sum(n_cases[lineage == 'unknown'])) %>%
  ungroup() %>%
  mutate(fraction_vic_in_assigned = n_vic / (n_vic + n_yam),
         fraction_yam_in_assigned = 1 - fraction_vic_in_assigned)
#nz_data %>% ggplot(aes(x = age, y = n_cases)) +
#  geom_col(aes(fill = lineage)) + 
#  facet_grid(observation_year~.)

# Assignining unidentified cases in seasons dominated by a single lineage
nz_data_untyped_assigned <- esr_nz_data
nz_data_untyped_assigned$lineage[nz_data_untyped_assigned$lineage == 'unknown' & nz_data_untyped_assigned$observation_year == 2013] <- 'B/Yamagata'
nz_data_untyped_assigned$lineage[nz_data_untyped_assigned$lineage == 'unknown' & nz_data_untyped_assigned$observation_year == 2014] <- 'B/Yamagata'
nz_data_untyped_assigned$lineage[nz_data_untyped_assigned$lineage == 'unknown' & nz_data_untyped_assigned$observation_year == 2017] <- 'B/Yamagata'
nz_data_untyped_assigned$lineage[nz_data_untyped_assigned$lineage == 'unknown' & nz_data_untyped_assigned$observation_year == 2019] <- 'B/Victoria'

nz_data_untyped_assigned <- nz_data_untyped_assigned %>%
  filter(lineage != 'unknown') %>%
  group_by(observation_year, age, lineage, severity) %>%
  summarise(n_cases = n()) %>%
  ungroup() %>%
  mutate(region = 'AUSNZ', country = 'New Zealand', cohort_type = 'age') %>%
  mutate(minimum_birth_year = observation_year - age - 1) %>%
  rename(cohort_value = age) %>%
  select(severity, country, region, observation_year, cohort_type, cohort_value, minimum_birth_year, lineage, n_cases)
  
# --------- Exporting raw data ----------
# all surveillance combined
nz_data %>% 
  filter(lineage != 'unknown') %>%
  group_by(country, region, observation_year, cohort_type,
           cohort_value, minimum_birth_year, lineage) %>%
  summarise(n_cases = sum(n_cases)) %>%
  write.csv('../results/processed_data/case_data_batches/case_data_nz_2012-2019_all_surveillance.csv', row.names = F)

# Sentinel (GP only)
nz_data %>% 
  filter(lineage != 'unknown', severity == 'GP consultation') %>%
  write.csv('../results/processed_data/case_data_batches/case_data_nz_2012-2019_sentinel_only.csv', row.names = F)

# Non-sentinel (hospital admission, referred to as non-sentinel for consistency with pre-2012 data but technically part of sentinel surveillance)
nz_data %>% 
  filter(lineage != 'unknown', severity == 'Hospital admission') %>%
  write.csv('../results/processed_data/case_data_batches/case_data_nz_2012-2019_nonsentinel_only.csv', row.names = F)

# --------- Exporting data with unidentified cases assigned where possible ----------
# All surveillance combined
nz_data_untyped_assigned %>%
  filter(lineage != 'unknown') %>%
  group_by(country, region, observation_year, cohort_type,
           cohort_value, minimum_birth_year, lineage) %>%
  summarise(n_cases = sum(n_cases)) %>%
  write.csv('../results/processed_data/case_data_batches/case_data_nz_2012-2019_all_surveillance_untyped_assigned.csv', row.names = F)

# Sentinel (GP only)
nz_data_untyped_assigned %>% 
  filter(lineage != 'unknown', severity == 'GP consultation') %>%
  write.csv('../results/processed_data/case_data_batches/case_data_nz_2012-2019_sentinel_only_untyped_assigned.csv', row.names = F)

# Non-sentinel (hospital admission)
nz_data_untyped_assigned %>% 
  filter(lineage != 'unknown', severity == 'Hospital admission') %>%
  write.csv('../results/processed_data/case_data_batches/case_data_nz_2012-2019_nonsentinel_only_untyped_assigned.csv', row.names = F)

