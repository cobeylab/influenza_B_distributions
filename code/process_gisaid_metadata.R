library(dplyr)
library(tidyr)
# Merges gisaid metadata files into a single file


years_downloaded <- as.character(seq(1980,2019))
years_downloaded <- c('pre_1980', years_downloaded)

gisaid_B_data <- tibble()

# Import and merge files
for(year in years_downloaded){
  file_path <- paste('../data/gisaid_metadata/metadata_files/', year, '_B_meta.csv', sep = '')
  year_tibble <- as_tibble(read.csv(file_path, header = T)) %>%
    select(Isolate_Id, Isolate_Name, Lineage,
           Collection_Date, Location, Passage_History) %>%
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
  rename(isolate_id = Isolate_Id, isolate_name = Isolate_Name, lineage = Lineage) %>%
  mutate(lineage = ifelse(lineage == '', NA, lineage))

write.csv(gisaid_B_data, '../data/gisaid_metadata/gisaid_metadata.csv', row.names = F)

