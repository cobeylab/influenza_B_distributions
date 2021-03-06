library(dplyr)
library(tidyr)
library(ggplot2)
library(cowplot)
theme_set(theme_cowplot())

# Merges gisaid metadata files into a single file
years_downloaded <- as.character(seq(1980,2019))
years_downloaded <- c('pre_1980', years_downloaded)

gisaid_B_data <- tibble()

# Import and merge files
for(year in years_downloaded){
  file_path <- paste('../data/gisaid_metadata/metadata_files/', year, '_B_meta.csv', sep = '')
  year_tibble <- as_tibble(read.csv(file_path, header = T)) %>%
    select(Isolate_Id, Isolate_Name, Lineage,
           Collection_Date, Location, Passage_History, Host_Age, Host_Age_Unit) %>%
    mutate(Collection_Date = as.character(Collection_Date))
  gisaid_B_data <- bind_rows(gisaid_B_data, year_tibble)
}

# Parse collection date and location to get collection year and country.
gisaid_B_data <- gisaid_B_data %>% separate(Location, c('continent','country'), sep = ' / ',
                                            extra = 'drop', fill = 'right') %>%
  # Remove non-standard collection dates (separated by '/')
  filter(grepl('/', Collection_Date) == F) %>%
  # Parse collection dates
  separate(Collection_Date, c('year','month','day'), sep = '-', fill = 'right') %>%
  mutate(year = as.integer(year), month = as.integer(month)) 

gisaid_B_data %>% filter(continent %in% c('Europe','North America') | country == 'China', !is.na(month)) %>% 
  group_by(year) %>% 
  summarise(previous_season = sum(month < 10), next_season = sum(month >=10)) %>%
  mutate(fraction_previous_season = previous_season / (previous_season + next_season)) %>%
  ungroup() %>% 
  summarise(median(fraction_previous_season))

gisaid_B_data <- gisaid_B_data %>%
  # Convert year to season (retain name of the variable as 'year')
  rowwise() %>%
  mutate(year = ifelse(continent %in% c('Europe','North America') | country == 'China' & 
                         (month < 10 | is.na(month)),
                       year - 1, year)) %>%
  ungroup() %>%
  # Rename variables to lower case
  rename(isolate_id = Isolate_Id, isolate_name = Isolate_Name, lineage = Lineage,
         host_age = Host_Age, host_age_unit = Host_Age_Unit) %>%
  mutate(lineage = as.character(lineage)) %>%
  mutate(lineage = ifelse(lineage == '', NA, lineage))

# EXPORT data used to calculate lineage frequencies

write.csv(gisaid_B_data, '../data/gisaid_metadata/gisaid_metadata.csv', row.names = F)

# EXPORT data on isolates with lineage and age info to use as case data:

gisaid_data_by_lineage_and_age <- gisaid_B_data %>%
  filter(!is.na(lineage), host_age_unit %in% c('Y','M')) %>%
  filter(!is.na(host_age)) %>%
  mutate(lineage = paste0('B/',lineage)) %>%
  mutate(host_age = ifelse(host_age_unit == 'M', host_age/12, host_age)) %>%
  mutate(host_age = floor(host_age)) %>%
  select(isolate_name, lineage, year, continent, country, host_age) %>% 
  unique()

gisaid_data_by_lineage_and_age %>%
  group_by(continent) %>%
  count()

gisaid_data_by_lineage_and_age %>% 
  filter(continent == 'Europe') %>%
  group_by(country) %>%
  count() %>%
  arrange(desc(n)) %>%
  ungroup() %>%
  mutate(cum_sum = cumsum(n))

# EU cases:
countries_outside_EU <- c('Russian Federation','Norway', 'Ukraine','Estonia', 'Moldova, Republic of',
                      'Albania', 'Macedonia, the former Yogoslav Republic of', 'Bosnia and Herzegovina',
                      'Belarus', 'Montenegro')

gisaid_EU_data_by_lineage_and_age <- gisaid_data_by_lineage_and_age %>%
  filter(continent == 'Europe', !(country %in% countries_outside_EU)) %>%
  group_by(year, host_age, lineage) %>%
  summarise(n_cases = n()) %>%
  ungroup() %>%
  mutate(minimum_birth_year = year - host_age - 1) %>%
  rename(cohort_value = host_age, observation_year = year) %>%
  mutate(country = 'Europe', region = 'Europe', cohort_type = 'age') %>%
  select(country, region, observation_year, cohort_type, cohort_value,
         minimum_birth_year, lineage, n_cases) %>%
  arrange(observation_year, cohort_value)
  

write.csv(gisaid_EU_data_by_lineage_and_age %>% filter(observation_year >= 2008),
          '../results/processed_data/case_data_europe_gisaid.csv', row.names = F)

# China cases:
gisaid_data_by_lineage_and_age %>%
  filter(country %in% c('China','Hong Kong (SAR)')) %>%
  group_by(country) %>% count()

gisaid_data_by_lineage_and_age %>%
  filter(country %in% c('China','Hong Kong (SAR)')) %>%
  group_by(year) %>% count()

gisaid_china_data_by_lineage_and_age <- gisaid_data_by_lineage_and_age %>%
  filter(country %in% c('China','Hong Kong (SAR)')) %>%
  group_by(year, host_age, lineage) %>%
  summarise(n_cases = n()) %>%
  ungroup() %>%
  mutate(minimum_birth_year = year - host_age - 1) %>%
  rename(cohort_value = host_age, observation_year = year) %>%
  mutate(country = 'China', region = 'China', cohort_type = 'age') %>%
  select(country, region, observation_year, cohort_type, cohort_value,
         minimum_birth_year, lineage, n_cases) %>%
  arrange(observation_year, cohort_value)

write.csv(gisaid_china_data_by_lineage_and_age %>% filter(observation_year >= 2009),
          '../results/processed_data/case_data_china_gisaid.csv', row.names = F)


# PLOTS
pooled_age_dist_europe <- gisaid_EU_data_by_lineage_and_age %>%
  filter(minimum_birth_year >= 1959) %>%
  group_by(minimum_birth_year, lineage) %>%
  summarise(n_cases = sum(n_cases)) %>%
  group_by(lineage) %>%
  mutate(fraction_cases = n_cases / sum(n_cases)) %>%
  ggplot(aes(x = minimum_birth_year, y = fraction_cases)) + geom_col() + facet_grid(.~lineage) +
  xlab('Year of birth') +
  ylab('Fraction of isolates') +
  background_grid() +
  xlim(1959,2020.5)
  
save_plot('../figures/gisaid_age_dist_europe_pooled.pdf', pooled_age_dist_europe,
          base_height = 6, base_width = 12)

age_dist_by_season_europe <- gisaid_EU_data_by_lineage_and_age %>%
  mutate(fraction_cases = n_cases / sum(n_cases)) %>%
  ggplot(aes(x = minimum_birth_year, y = n_cases)) + geom_col() + 
  facet_grid(observation_year~lineage, scales = 'free') +
  xlab('Year of birth') +
  ylab('Number of isolates') +
  background_grid() 
save_plot('../figures/gisaid_age_dist_europe_by_season.pdf', age_dist_by_season_europe,
          base_height = 10, base_width = 9)

pooled_age_dist_china <- gisaid_china_data_by_lineage_and_age %>%
  filter(minimum_birth_year >= 1959) %>%
  group_by(minimum_birth_year, lineage) %>%
  summarise(n_cases = sum(n_cases)) %>%
  group_by(lineage) %>%
  mutate(fraction_cases = n_cases / sum(n_cases)) %>%
  ggplot(aes(x = minimum_birth_year, y = fraction_cases)) + geom_col() + facet_grid(.~lineage) +
  xlab('Year of birth') +
  ylab('Fraction of isolates') +
  background_grid() +
  xlim(1959,2020.5)
save_plot('../figures/gisaid_age_dist_china_pooled.pdf', pooled_age_dist_china,
          base_height = 5, base_width = 9)

age_dist_by_season_china <- gisaid_china_data_by_lineage_and_age %>%
  mutate(fraction_cases = n_cases / sum(n_cases)) %>%
  ggplot(aes(x = minimum_birth_year, y = n_cases)) + geom_col() + 
  facet_grid(observation_year~lineage, scales = 'free') +
  xlab('Year of birth') +
  ylab('Number of isolates') +
  background_grid() 
save_plot('../figures/gisaid_age_dist_china_by_season.pdf', age_dist_by_season_china,
          base_height = 10, base_width = 9)

isolate_data_by_country_EU <- gisaid_data_by_lineage_and_age %>%
  filter(continent == 'Europe', !(country %in% countries_outside_EU)) %>%
  group_by(year,country) %>%
  summarise(n_isolates = n()) %>%
  group_by(year) %>%
  mutate(total_isolates = sum(n_isolates),
         fraction = n_isolates / total_isolates) %>%
  ungroup()

isolate_data_by_country_EU  %>%
  group_by(country) %>%
  summarise(n_isolates = sum(n_isolates)) %>%
  ungroup() %>%
  arrange(desc(n_isolates)) %>%
  mutate(fraction = n_isolates / sum(n_isolates),
         cum_fraction = cumsum(fraction))


# Export distribution by country of isolates used to fit age distributions in the EU
# (for plotting together with the distribution by country of isolates used to calculate lin. freqs.)
isolate_data_by_country_EU <- isolate_data_by_country_EU %>%
  mutate(data_type = 'Isolate data used to fit age distributions (European Union only)')

write.csv(isolate_data_by_country_EU, '../results/processed_data/europe_age_dist_isolate_distribution.csv',
          row.names = F)
