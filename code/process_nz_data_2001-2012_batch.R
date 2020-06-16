library(dplyr)
library(ggplot2)
library(cowplot)
theme_set(theme_cowplot())
library(stringr)
library(tidyr)
library(reshape2)

# New Zealand Data from ESR in 2001-2002, with type of surveillance information.
esr_nz_data <- as_tibble(read.csv('../data/case_data/nz_data_2001-2012_batch.csv', header = T)) %>%
  mutate(sAge = as.character(sAge))

# Standardize age column:
esr_nz_data$sAge[grepl('#NAME?', esr_nz_data$sAge)] <- NA
esr_nz_data$sAge[grepl('Y', esr_nz_data$sAge)] <- str_replace(esr_nz_data$sAge[grepl('Y', esr_nz_data$sAge)],'Y','')
esr_nz_data$sAge[grepl('y', esr_nz_data$sAge)] <- str_replace(esr_nz_data$sAge[grepl('y', esr_nz_data$sAge)],'y','')

# For ages listed in months, days, or weeks, set age to greatest complete year of life 
esr_nz_data$sAge[grepl('M', esr_nz_data$sAge)] <- as.character(floor(as.numeric(str_replace(esr_nz_data$sAge[grepl('M', esr_nz_data$sAge)],'M',''))/12))
esr_nz_data$sAge[grepl('m', esr_nz_data$sAge)] <- as.character(floor(as.numeric(str_replace(esr_nz_data$sAge[grepl('m', esr_nz_data$sAge)],'m',''))/12))
esr_nz_data$sAge[grepl('D', esr_nz_data$sAge)]  <- as.character(floor(as.numeric(str_replace(esr_nz_data$sAge[grepl('D', esr_nz_data$sAge)],'D',''))/365))
esr_nz_data$sAge[grepl('W', esr_nz_data$sAge)] <- as.character(floor(as.numeric(str_replace(esr_nz_data$sAge[grepl('W', esr_nz_data$sAge)],'W',''))/53))
esr_nz_data$sAge[grepl('wk', esr_nz_data$sAge)] <- as.character(floor(as.numeric(str_replace(esr_nz_data$sAge[grepl('wk', esr_nz_data$sAge)],'wk',''))/53))

esr_nz_data$sAge[esr_nz_data$sAge == '-'] <- NA
esr_nz_data$sAge[esr_nz_data$sAge == ''] <- NA

# Check no ages with characters in them remain.
stopifnot(length(esr_nz_data$sAge[grepl('[A-Z|a-z]', esr_nz_data$sAge) & !is.na(esr_nz_data$sAge)]) == 0)

# Convert all ages to complete years
esr_nz_data <- esr_nz_data %>%
  mutate(sAge = floor(as.numeric(sAge))) %>%
  rename(age = sAge,  year = iYear)
  
# Order age groups
esr_nz_data$sAgeGrp <- factor(esr_nz_data$sAgeGrp, levels = c("0<1 year", "1-4 yrs", "5-19 yrs", "20-34 yrs", "35-49 yrs","50-64 yrs", "65 + yrs","Unknown"))

# Named vector to assign lineage to test strains in tibble
test_strains <- c('B/Johannesburg/5/99 - like' = 'B/Yamagata', 
                  'B/Yamanashi/166/98 - like' = 'B/Yamagata',
                  'B/Malaysia/2506/2004-like' = 'B/Victoria',
                  'B/Wisconsin/1/2010 - like virus' = 'B/Yamagata',
                  'B/Shanghai/361/2002-like' = 'B/Yamagata',
                  'B/Florida/4/2006 - like' = 'B/Yamagata',
                  'B/Sichuan/379/99 - like' = 'B/Yamagata',
                  'Untypable' = 'unknown',
                  'Influenza B virus identified' = 'unknown',
                  'B/Hong Kong/330/2001 - like' = 'B/Victoria',
                  'Typing to follow' = 'unknown',
                  'B/Brisbane/60/2008 - like virus isolated' = 'B/Victoria'
                  )

stopifnot(all(test_strains[c('B/Johannesburg/5/99 - like','B/Johannesburg/5/99 - like', 'B/Brisbane/60/2008 - like virus isolated',
               'B/Sichuan/379/99 - like', 'B/Malaysia/2506/2004-like', 'Untypable')] ==
            c('B/Yamagata', 'B/Yamagata', 'B/Victoria','B/Yamagata','B/Victoria','unknown')))

# Assign lineages
esr_nz_data <- esr_nz_data %>% mutate(lineage = test_strains[as.character(sAntigenicStrain)]) %>%
  # Replace lineage = NA with lineage = unknown
  rowwise() %>%
  mutate(lineage = ifelse(is.na(lineage),'unknown',lineage)) %>%
  ungroup()

# A visual check that lineage assignment was correct and unambiguous.
esr_nz_data %>%
  group_by(sAntigenicStrain) %>%
  summarise(lineage = paste(unique(lineage), sep = ';'))

# ====== EXPORT DATA AGGREGATED BY AGE =======
nz_data <- esr_nz_data  %>% 
  filter(lineage != 'unknown', is.na(age) == F) %>%
  group_by(year, age, lineage, SurveillanceType) %>%
  summarise(n_cases = n()) %>% 
  ungroup() %>%
  mutate(region = 'AUSNZ', country = 'New Zealand', cohort_type = 'age') %>%
  mutate(minimum_birth_year = year - age - 1) %>%
  rename(cohort_value = age, observation_year = year) %>%
  select(SurveillanceType, country, region, observation_year, cohort_type, cohort_value, minimum_birth_year, lineage, n_cases)

# Assigns a lineage to unidentified cases in seasons with a dominant lineage
nz_data_untyped_assigned <- esr_nz_data 
nz_data_untyped_assigned$lineage[nz_data_untyped_assigned$lineage == 'unknown' & nz_data_untyped_assigned$year == 2002] <- 'B/Victoria'
nz_data_untyped_assigned$lineage[nz_data_untyped_assigned$lineage == 'unknown' & nz_data_untyped_assigned$year == 2011] <- 'B/Victoria'

nz_data_untyped_assigned <- nz_data_untyped_assigned %>% 
  filter(lineage != 'unknown', is.na(age) == F) %>%
  group_by(year, age, lineage, SurveillanceType) %>%
  summarise(n_cases = n()) %>% 
  ungroup() %>%
  mutate(region = 'AUSNZ', country = 'New Zealand', cohort_type = 'age') %>%
  mutate(minimum_birth_year = year - age - 1) %>%
  rename(cohort_value = age, observation_year = year) %>%
  select(SurveillanceType, country, region, observation_year, cohort_type, cohort_value, minimum_birth_year, lineage, n_cases)
  
# New Zealand only (all surveillance)
# Remove non-sentinel data from 2012 to avoid overlap with data batch for 2012-2019
nz_data %>%
  filter(!(SurveillanceType == 'Other Influenza Surveillance' & observation_year == 2012)) %>%
  group_by(country, region, observation_year, cohort_type,
           cohort_value, minimum_birth_year, lineage) %>%
  summarise(n_cases = sum(n_cases)) %>%
  write.csv('../results/processed_data/case_data_batches/case_data_nz_2001-2012_all_surveillance.csv', row.names = F)

# New Zealand only (sentinel only)
nz_data %>%
  filter(SurveillanceType == 'Sentinel Surveillance') %>%
  write.csv('../results/processed_data/case_data_batches/case_data_nz_2001-2012_sentinel_only.csv', row.names = F)

# New Zealand only (nonsentinel only)
# Remove non-sentinel data from 2012 to avoid overlap with data batch for 2012-2019
nz_data %>%
  filter(SurveillanceType == 'Other Influenza Surveillance', observation_year != 2012) %>%
  write.csv('../results/processed_data/case_data_batches/case_data_nz_2001-2012_nonsentinel_only.csv', row.names = F)

# ------------- New Zealand data with untyped cases assigned where possible --------------
# All surveillance, removing non-sentinel data from 2012 to avoid overlap with data batch for 2012-2019
nz_data_untyped_assigned %>%
  filter(!(SurveillanceType == 'Other Influenza Surveillance' & observation_year == 2012)) %>%
  group_by(country, region, observation_year, cohort_type,
           cohort_value, minimum_birth_year, lineage) %>%
  summarise(n_cases = sum(n_cases)) %>%
  write.csv('../results/processed_data/case_data_batches/case_data_nz_2001-2012_all_surveillance_untyped_assigned.csv', row.names = F)

# Sentinel only with untyped cases assigned where possible: (sentinel only)
nz_data_untyped_assigned %>%
  filter(SurveillanceType == 'Sentinel Surveillance') %>%
  group_by(country, region, observation_year, cohort_type,
           cohort_value, minimum_birth_year, lineage) %>%
  summarise(n_cases = sum(n_cases)) %>%
  write.csv('../results/processed_data/case_data_batches/case_data_nz_2001-2012_sentinel_only_untyped_assigned.csv', row.names = F)

# Non-sentinel only with untyped cases assigned where possible: (nonsentinel only)
nz_data_untyped_assigned %>%
  # Remove non-sentinel data from 2012 to avoid overlap with data batch for 2012-2019
  filter(SurveillanceType == 'Other Influenza Surveillance', observation_year != 2012) %>%
  group_by(country, region, observation_year, cohort_type,
           cohort_value, minimum_birth_year, lineage) %>%
  summarise(n_cases = sum(n_cases)) %>%
  write.csv('../results/processed_data/case_data_batches/case_data_nz_2001-2012_nonsentinel_only_untyped_assigned.csv', row.names = F)
