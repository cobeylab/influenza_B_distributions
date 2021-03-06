library(tidyr)
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

# Function to expand aggregated age distribution into a tibble at the one-year resolution
# Distributes the number of persons equally by each year in an interval
# (Used for the aggregated Chinese demography data)
expand_age_distribution <- function(dem_data){
  expanded_tibble <- c()
  for(y in unique(dem_data$year)){
    for(group in unique(dem_data$age_group)){
      bounds <- as.numeric(str_split(group, '-')[[1]])
      expanded_tibble <- bind_rows(expanded_tibble,
                                   tibble(year = y, age_group = group, age_group_length = bounds[2] - bounds[1] + 1,
                                          age = bounds[1]:bounds[2]))
      
    }
  }
  return(left_join(expanded_tibble, dem_data %>% mutate(age_group = as.character(age_group))) %>%
           rename(total_in_age_group = n_persons) %>%
           mutate(n_persons = total_in_age_group / age_group_length) %>%
           select(year, age, n_persons))
}

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
  
  # Process EU demography
  format_eu_demography <- function(eu_demography_file_path){
    as_tibble(read.table(eu_demography_file_path, sep = '\t', header = T)) %>%
      rename(age = freq.unit.age.sex.geo.TIME_PERIOD) %>%
      mutate(age = as.character(age)) %>%
      mutate(age = ifelse(grepl('Y_LT1',age), 'Y0', age)) %>%
      mutate(age = str_extract(age,'Y[0-9]*')) %>%
      mutate(across(matches('X20'), as.character)) %>%
      pivot_longer(cols = matches('X20'),
                   names_to = 'year', values_to = 'n_persons') %>%
      mutate(n_persons = ifelse(n_persons == ': ',NA, n_persons)) %>%
      mutate(age = as.integer(str_replace(age, 'Y','')),
             year = as.integer(str_replace(year,'X','')),
             n_persons = as.integer(str_remove_all(n_persons,'[a-z]*'))) %>%
      arrange(age) %>%
      select(year, age, n_persons)
  }
  
  eu_demography_2007_2012 <- format_eu_demography('../data/demographic_data/EU_demographics_27_countries_2007-2012.tsv')
  eu_demography_2013_2020 <- format_eu_demography('../data/demographic_data/EU_demographics_28_countries_2013-2020.tsv')
  
  # Process China demography
  china_demography <- as_tibble(read.csv('../data/demographic_data/China_age_distribution.csv')) %>%
    filter(age_group != '95+')
  
  # Fill missing years using the average of neighboring years
  china_demography <- bind_rows(china_demography,
                                china_demography %>% filter(year %in% c(2009, 2011)) %>%
                                  group_by(age_group) %>% summarise(n_persons = mean(n_persons)) %>%
                                  mutate(year = 2010) %>% ungroup() %>% select(year, age_group, n_persons))
  # Use 2018 for 2019, since 2018 is the last year in the National Statistics Bureau data.
  china_demography <- bind_rows(china_demography,
                                china_demography %>% filter(year == 2018) %>% mutate(year = 2019))
  
  china_demography <- expand_age_distribution(china_demography)
  
  
  # Merge tibbles
  demographic_data <- bind_rows(
    mutate(aus_demographics, country = 'Australia', region = 'AUSNZ'),
    mutate(nz_demographics, country = 'New Zealand', region = 'AUSNZ'),
    mutate(eu_demography_2007_2012, country = 'Europe', region = 'Europe'),
    mutate(eu_demography_2013_2020, country = 'Europe', region = 'Europe'),
    china_demography %>% mutate(country = 'China', region = 'China')
    ) %>% 
    # Rename age as "cohort_value" and create new variable "cohort_type" equal to "age"
    rename(cohort_value = age) %>%
    mutate(cohort_type = 'age') %>%
    select(country, region, year, cohort_type, cohort_value, n_persons)
    
  demographic_data <- demographic_data %>%
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

normalize_demographic_data(demographic_data %>% filter(country == 'New Zealand', cohort_value < 90),1500) %>%
  mutate(cohort_value = as.numeric(cohort_value)) %>%
  mutate(min_birth_year = observation_year - cohort_value - 1) %>%
  group_by(country, observation_year) %>%
  summarise(mean_age = sum(cohort_value*fraction, na.rm = T),
            mean_birth_year = sum(min_birth_year*fraction, na.rm = T)) %>%
  filter(observation_year >= 2001, observation_year <= 2019) %>%
  ggplot(aes(x = observation_year, y = mean_age, color = country)) +
  geom_line()+
  geom_point()

general_population_age <- normalize_demographic_data(demographic_data %>% 
                                                       filter(country %in% c('Australia', 'New Zealand')),
                                                     1900) %>%
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

normalize_demographic_data(demographic_data %>% filter(country %in% c('Australia', 'New Zealand')), 1952) %>%
  mutate(cohort_value = as.numeric(cohort_value)) %>%
  mutate(min_birth_year = observation_year - cohort_value - 1) %>%
  filter(observation_year %in% c(2002,2005,2008,2011)) %>%
  ggplot(aes(x = min_birth_year, y = fraction, color = country)) +
  geom_line() +
  facet_grid(observation_year~.)


# Calculate vaccine coverage in 2016 using ESR report
age_group_levels <- c('<1','1-4','5-19','20-34','35-49','50-64','65plus')

nz_demographics_2016_by_age_group <- demographic_data %>% filter(country == 'New Zealand', observation_year == 2016) %>%
  mutate(cohort_value = as.numeric(cohort_value)) %>%
  mutate(age_group = 
           case_when(cohort_value <1 ~ '<1',
                     (cohort_value >=1 & cohort_value <= 4) ~ '1-4',
                     (cohort_value >=5 & cohort_value <= 19) ~ '5-19',
                     (cohort_value >=20 & cohort_value <= 34) ~ '20-34',
                     (cohort_value >=35 & cohort_value <= 49) ~ '35-49',
                     (cohort_value >=50 & cohort_value <= 64) ~ '50-64',
                     (cohort_value >= 65) ~ '65plus')
         ) %>%
  mutate(age_group = factor(age_group, levels = age_group_levels)) %>%
  group_by(age_group) %>%
  summarise(n_persons = sum(n_persons)) %>%
  ungroup() 

vaccine_coverage_nz_2016 <- as_tibble(read.csv('../data/vaccine_coverage_nz_2016.csv'))

vaccine_coverage_nz_2016 <- left_join(nz_demographics_2016_by_age_group, vaccine_coverage_nz_2016, by = 'age_group') %>%
  mutate(age_group = factor(age_group, levels = age_group_levels)) %>%
  mutate(coverage = n_vaccines_distributed/n_persons * 100) %>% arrange(age_group)

vaccine_coverage_nz_2016 %>%
  filter((age_group %in% c('<1','1-4','5-19'))) %>%
  summarise(sum(n_persons*coverage)/sum(n_persons))

vaccine_coverage_nz_2016 %>%
  filter((age_group %in% c('<1','1-4','5-19','65plus')) == F) %>%
  summarise(sum(n_persons*coverage)/sum(n_persons))



