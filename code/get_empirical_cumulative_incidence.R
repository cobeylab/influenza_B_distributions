library(lubridate)
library(dplyr)
library(tidyr)

# Function for getting number of ISO weeks in a year
n_iso_weeks <- function(year){
  last_days <- paste(year, '12', seq(20,31), sep = '-')
  n_iso_weeks <- max(lubridate::isoweek(last_days))
  return(n_iso_weeks)
}

get_empirical_cumulative_incidence <- function(case_data){
  cumulative_incidence <- case_data  %>%
    filter(is.na(observation_year) == F) %>%
    mutate(observation_date = paste(observation_year, observation_month, observation_day, sep = '-'), 
           observation_date = ymd(observation_date),
           # In S. Hemisphere, season week number is simply calendar week number
           season_week = ifelse(is.na(observation_date), NA, lubridate::isoweek(observation_date))) %>%
    filter(!is.na(season_week)) %>%
    rename(season_start_year = observation_year) %>% # Identify season as "season_start_year"
    group_by(season_start_year, season_week) %>%
    # Count number of cases by week in the season (1,2,...52/53) for each season and country
    summarise(weekly_cases = n()) %>%
    ungroup() %>% 
    group_by(season_start_year) %>%
    # Cumulative number of cases in each week in each season/ country
    mutate(cumulative_n_cases = cumsum(weekly_cases),
           season_total_cases = sum(weekly_cases)) %>%
    # Cumulative fraction of cases in each week in each season/ country
    mutate(cumulative_B_fraction = cumulative_n_cases /season_total_cases) %>%
    # Add number of weeks in each year
    rowwise() %>%
    mutate(n_weeks_in_year = n_iso_weeks(season_start_year)) %>%
    # Cumulative fraction of the season duration in each week
    mutate(cumulative_season_fraction = season_week / n_weeks_in_year) %>%
    ungroup()
  return(cumulative_incidence)
}

case_data_nz_batch1 <- as_tibble(read.csv('../data/case_data/nz_data_2001-2012_batch.csv')) %>%
  select(iYear, dDateSpecimenTaken) %>%
  rename(date = dDateSpecimenTaken) %>%
  mutate(date = as.character(date)) %>%
  separate(date, into = c('observation_month','observation_day','two_digit_year'), sep = '/') %>%
  select(-two_digit_year) %>%
  rename(observation_year = iYear) %>%
  filter(observation_year != '2012') # To avoid overlap with second batch
  
case_data_nz_batch2 <- as_tibble(read.csv('../data/case_data/nz_data_2012-2019_batch.csv')) %>%
  select(year,date_specimen) %>%
  mutate(date_specimen = as.character(date_specimen)) %>%
  separate(date_specimen, into = c('observation_month','observation_day','two_digit_year')) %>%
  select(-two_digit_year) %>%
  rename(observation_year = year) 

case_data_australia <- read.table('../data/case_data/Vijaykrishna2015_Australian_cases.tsv', sep = '\t', header = T)
case_data_australia <- as_tibble(case_data_australia) %>%
  separate(date, c('observation_year','observation_month','observation_day'),sep = '-') %>%
  mutate(observation_year = as.integer(observation_year), observation_month = as.integer(observation_month)) %>%
  select(observation_year, observation_month, observation_day)

bind_rows(get_empirical_cumulative_incidence(case_data_nz_batch1),
          get_empirical_cumulative_incidence(case_data_nz_batch2)) %>%
  write.csv('../results/processed_data/raw_cumulative_incidence_curves/raw_cumulative_incidence_curves_nz.csv', row.names = F)

get_empirical_cumulative_incidence(case_data_australia) %>%
  write.csv('../results/processed_data/raw_cumulative_incidence_curves/raw_cumulative_incidence_curves_aus.csv', row.names = F)
