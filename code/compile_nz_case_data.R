library(dplyr)
library(tidyr)
library(stringr)
library(cowplot)
theme_set(theme_cowplot())
library(ggplot2)

# Read nz data
nz_data <- lapply(as.list(list.files('../results/processed_data/case_data_batches/', full.names = T)),
                  FUN = read.csv)
names(nz_data) <- str_replace(list.files('../results/processed_data/case_data_batches/', full.names = F),'\\.csv','')

# Compile data from both surveillance types combined, separately with and without assignment of unidentified cases
bind_rows(nz_data$`case_data_nz_2001-2012_all_surveillance`,
          nz_data$`case_data_nz_2012-2019_all_surveillance`) %>%
  write.csv('../results/processed_data/case_data_nz_all_surveillance.csv', row.names = F)

bind_rows(nz_data$`case_data_nz_2001-2012_all_surveillance_untyped_assigned`,
          nz_data$`case_data_nz_2012-2019_all_surveillance_untyped_assigned`) %>%
  write.csv('../results/processed_data/case_data_nz_all_surveillance_untyped_assigned.csv', row.names = F)

# Data from sentinel cases only (meaning GP, since hospitalizations post 2012 are labelled sentinel by ESR)
# No overlap since 2012 data from second batch are all hospitalizations...
bind_rows(nz_data$`case_data_nz_2001-2012_sentinel_only_untyped_assigned`,
          nz_data$`case_data_nz_2012-2019_sentinel_only_untyped_assigned`) %>%
  write.csv('../results/processed_data/case_data_nz_sentinel_only_untyped_assigned.csv', row.names = F)

# Data from non-sentinel cases only (meaning hospitalizations, since hospitalizations post 2012 are labelled sentinel by ESR)
# 2012 non-sentinel cases from first batch removed to avoid double counting with 2012 from second batch.
bind_rows(nz_data$`case_data_nz_2001-2012_nonsentinel_only_untyped_assigned`,
          nz_data$`case_data_nz_2012-2019_nonsentinel_only_untyped_assigned`) %>%
  write.csv('../results/processed_data/case_data_nz_nonsentinel_only_untyped_assigned.csv', row.names = F)

raw_sentinel_cases <- as_tibble(bind_rows(nz_data$`case_data_nz_2001-2012_sentinel_only`,
                                nz_data$`case_data_nz_2012-2019_sentinel_only`))

raw_nonsentinel_cases <- as_tibble(bind_rows(nz_data$`case_data_nz_2001-2012_nonsentinel_only`,
                                   nz_data$`case_data_nz_2012-2019_nonsentinel_only`))


full_join(raw_sentinel_cases %>% rename(n_sentinel = n_cases),
          raw_nonsentinel_cases %>% rename(n_nonsentinel = n_cases),
          by = c('country','region','observation_year', 'cohort_type','cohort_value','minimum_birth_year','lineage')) %>%
  group_by(country,region, minimum_birth_year, lineage) %>%
  summarise(n_sentinel = sum(n_sentinel, na.rm = T),
         n_nonsentinel = sum(n_nonsentinel, na.rm = T)) %>%
  ungroup() %>%
  group_by(lineage) %>%
  mutate(total_sentinel_cases = sum(n_sentinel, na.rm = T),
         total_nonsentinel_cases = sum(n_nonsentinel, na.rm = T)) %>%
  ungroup() %>%
  filter(minimum_birth_year >= 1952) %>%
  mutate(fraction_sentinel = n_sentinel / total_sentinel_cases,
         fraction_nonsentinel = n_nonsentinel / total_nonsentinel_cases) %>%
  ggplot(aes(x = fraction_sentinel, y = fraction_nonsentinel)) +
  geom_text(aes(label = minimum_birth_year)) +
  facet_grid(.~lineage)

full_join(raw_sentinel_cases %>% rename(n_sentinel = n_cases),
          raw_nonsentinel_cases %>% rename(n_nonsentinel = n_cases),
          by = c('country','region','observation_year', 'cohort_type','cohort_value','minimum_birth_year','lineage')) %>%
  group_by(observation_year, lineage) %>%
  mutate(total_sentinel_cases = sum(n_sentinel, na.rm = T),
         total_nonsentinel_cases = sum(n_nonsentinel, na.rm = T)) %>%
  ungroup() %>%
  mutate(fraction_sentinel = n_sentinel / total_sentinel_cases,
         fraction_nonsentinel = n_nonsentinel / total_nonsentinel_cases) %>%
  filter(observation_year %in% c(2001,2002,2005,2007,2008, 2011, 2013, 2014, 2015, 2017, 2019)) %>%
  #filter(!(lineage == 'B/Victoria' & observation_year == 2007),
  #       !(lineage == 'B/Yamagata' & observation_year == 2002),
  #       !(lineage == 'B/Yamagata' & observation_year == 2003)) %>%
  ggplot(aes(x = minimum_birth_year)) +
  geom_col(aes(y = n_sentinel), color = 'red', fill = 'red', alpha = 0.7) +
  geom_col(aes(y = n_nonsentinel), fill = 'blue', alpha = 0.7) +
  facet_grid(observation_year~lineage, scales = 'free_y') +
  ggtitle('Year/lineage combinations with few cases omitted\nto improve visualization')



bind_rows(raw_sentinel_cases %>% mutate(surveillance_type = 'sentinel'),
          raw_nonsentinel_cases %>% mutate(surveillance_type = 'non-sentinel')) %>%
  group_by(country,region, minimum_birth_year, surveillance_type, lineage) %>%
  summarise(n_cases = sum(n_cases, na.rm = T)) %>%
  group_by(surveillance_type, lineage) %>%
  mutate(total_cases = sum(n_cases, na.rm = T))  %>%
  mutate(fraction_cases = n_cases / total_cases) %>%
  ungroup() %>%
  ggplot(aes(x = minimum_birth_year, y = fraction_cases,
             color = surveillance_type, group = surveillance_type)) +
  geom_line() +
  facet_grid(.~lineage)


bind_rows(raw_sentinel_cases %>% mutate(surveillance_type = 'sentinel'),
          raw_nonsentinel_cases %>% mutate(surveillance_type = 'non-sentinel')) %>%
  group_by(surveillance_type, lineage, observation_year) %>%
  mutate(total_cases = sum(n_cases, na.rm = T))  %>%
  mutate(fraction_cases = n_cases / total_cases) %>%
  ungroup() %>%
  filter(observation_year %in% c(2001,2002,2003,2005,2007, 2008, 2011)) %>%
  filter(!(lineage == 'B/Victoria' & observation_year == 2007),
         !(lineage == 'B/Yamagata' & observation_year == 2002)) %>%
  ggplot(aes(x = minimum_birth_year, y = fraction_cases,
             color = surveillance_type, fill = surveillance_type, group = surveillance_type)) +
  geom_col(position = 'dodge', alpha = 0.5) +
  facet_grid(observation_year~lineage, scales = 'free_y')
  

