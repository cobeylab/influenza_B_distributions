library(dplyr)
library(reshape2)
library(stringr)
source('model_functions.R')
library(cowplot)
library(ggplot2)
theme_set(theme_cowplot())

output_directory <- '../results/processed_data/'

# Get start_birth_year from basic_parameters.R
source('basic_parameters.R')

main <- function(){
  # Import Australian and NZ demographic data
  aus_demographics <- read.csv('../data/demographic_data/australian_age_distribution.csv', header = T)
  nz_demographics <- read.csv('../data/demographic_data/nz_age_distribution.csv', header = T)
  
  # Melt into dataframe with columns year, age, n_persons
  aus_demographics <- melt(aus_demographics, 
                           id.vars = 'year', variable.name = 'age', value.name = 'n_persons')
  nz_demographics <- melt(nz_demographics,
                          id.vars = 'year', variable.name = 'age', value.name = 'n_persons')
  
  # Convert into tibbles and parse out age value from 'age' column
  aus_demographics <- as_tibble(aus_demographics) %>%
    mutate(age = as.integer(str_extract(age,'[0-9]+')),
           year = as.integer(gsub('Jun-', '',year))
           )
  nz_demographics <- as_tibble(nz_demographics) %>%
    mutate(age = as.integer(str_extract(age,'[0-9]+')))
  
  # Merge Australian and NZ tibbles
  ausnz_demographic_data <- bind_rows(
    mutate(aus_demographics, country = 'Australia'),
    mutate(nz_demographics, country = 'New Zealand')
  ) %>% 
    # Add variable "region"
    mutate(region = 'AUSNZ') %>%
    # Rename age as "cohort_value" and create new variable "cohort_type" equal to "age"
    rename(cohort_value = age) %>%
    mutate(cohort_type = 'age') %>%
    select(country, region, year, cohort_type, cohort_value, n_persons)
    
  # Merge with AUS/NZ and US data:
  demographic_data <- ausnz_demographic_data %>%
    mutate(cohort_value = as.character(cohort_value))
  
  # Add minimum possible birth year for each cohort (cohorts identified by age in a given obs. year have 2 possible b. years)
  demographic_data <- demographic_data %>%
    rename(observation_year = year) %>% # Rename year as observation_year to match case data
    mutate(minimum_birth_year = ifelse(cohort_type == 'birth_year', as.integer(cohort_value),
                                       observation_year-as.integer(cohort_value) - 1)) %>%
    select(country,region,observation_year,cohort_type,cohort_value,minimum_birth_year,n_persons)
  
  return(demographic_data)
}
demographic_data <- main()
write.csv(demographic_data, paste0(output_directory, 'demographic_data.csv'), row.names = F)

normalize_demographic_data(demographic_data %>% filter(cohort_value < 90),1500) %>%
  filter(country %in% c('Australia','New Zealand')) %>%
  mutate(cohort_value = as.numeric(cohort_value)) %>%
  group_by(country, observation_year) %>%
  summarise(mean_age = sum(cohort_value*fraction, na.rm = T)) %>%
  filter(observation_year >= 2001, observation_year <= 2013) %>%
  ggplot(aes(x = observation_year, y = mean_age, color = country)) +
  geom_line()+
  geom_point()


general_population_age <- normalize_demographic_data(demographic_data, 1900) %>%
  filter(country %in% c('Australia','New Zealand')) %>%
  mutate(cohort_value = as.numeric(cohort_value)) %>%
  filter(observation_year %in% c(2002,2005,2008,2012)) %>%
  filter(!(cohort_value >= 90 & country == 'New Zealand')) %>%
  ggplot(aes(x = cohort_value, y = fraction, color = country)) +
  geom_line() +
  facet_grid(observation_year~.) +
  theme(legend.position = 'top',
        axis.text.y = element_text(size = 10)) +
  xlab('Age (years)') +
  ylab('Fraction of general population')
save_plot('../figures/general_population_age_distribution.pdf',general_population_age,
          base_height = 7, base_width = 5)

normalize_demographic_data(demographic_data, 1952) %>%
  filter(country %in% c('Australia','New Zealand')) %>%
  mutate(cohort_value = as.numeric(cohort_value)) %>%
  mutate(min_birth_year = observation_year - cohort_value - 1) %>%
  filter(observation_year %in% c(2002,2005,2008,2011)) %>%
  ggplot(aes(x = min_birth_year, y = fraction, color = country)) +
  geom_line() +
  facet_grid(observation_year~.)
  
