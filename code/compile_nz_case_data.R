library(dplyr)
library(tidyr)
library(stringr)
library(ggplot2)
library(cowplot)
library(viridis)
theme_set(theme_cowplot())


# Read nz data
nz_data <- lapply(as.list(list.files('../results/processed_data/case_data_batches/', full.names = T)),
                  FUN = read.csv)
names(nz_data) <- str_replace(list.files('../results/processed_data/case_data_batches/', full.names = F),'\\.csv','')

# Compile data from both surveillance types combined, separately with and without assignment of unidentified cases
case_data_nz_all_surveillance <- bind_rows(nz_data$`case_data_nz_2001-2012_all_surveillance`,
                                           nz_data$`case_data_nz_2012-2019_all_surveillance`)
case_data_nz_all_surveillance %>%
  write.csv('../results/processed_data/case_data_nz_all_surveillance.csv', row.names = F)

bind_rows(nz_data$`case_data_nz_2001-2012_all_surveillance_untyped_assigned`,
          nz_data$`case_data_nz_2012-2019_all_surveillance_untyped_assigned`) %>%
  write.csv('../results/processed_data/case_data_nz_all_surveillance_untyped_assigned.csv', row.names = F)

# Case data with assignment of unidentified cases but removing children < 10 years old
bind_rows(nz_data$`case_data_nz_2001-2012_all_surveillance_untyped_assigned`,
          nz_data$`case_data_nz_2012-2019_all_surveillance_untyped_assigned`) %>%
  filter(cohort_value >= 10) %>%
  write.csv('../results/processed_data/case_data_nz_all_surveillance_untyped_assigned_no_children.csv', row.names = F)

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

# Plot with correlation between sentinel, non-sentinel

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

sentinel_vs_nonsentinel_by_obs_year <- full_join(raw_sentinel_cases %>% rename(n_sentinel = n_cases),
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
  geom_col(aes(y = n_sentinel), color = 'goldenrod1', fill = 'goldenrod1', alpha = 0.7) +
  geom_col(aes(y = n_nonsentinel), fill = 'aquamarine4', alpha = 0.7) +
  facet_grid(observation_year~lineage, scales = 'free_y') +
  xlab('Year of birth') +
  ylab('Number of cases')  +
  background_grid() +
  theme(axis.text.y = element_text(size = 8))
save_plot('../figures/sentinel_vs_nonsentinel_by_obs_year.pdf',
          sentinel_vs_nonsentinel_by_obs_year,
          base_width = 7, base_height = 8)
#+ ggtitle('Year/lineage combinations with few cases omitted\nto improve visualization')

#Plot showing complete age distributions
# Including Australia
aus_cases <- as_tibble(read.csv('../results/processed_data/Australian_case_data.csv')) %>%
  mutate(surveillance_type = 'Australia')

byear_distribution_full <- bind_rows(list(raw_sentinel_cases %>% mutate(surveillance_type = 'New Zealand (sentinel)'),
               raw_nonsentinel_cases %>% mutate(surveillance_type = 'New Zealand (non-sentinel)'),
               aus_cases)) %>%
  group_by(country,region, minimum_birth_year, surveillance_type, lineage) %>%
  summarise(n_cases = sum(n_cases, na.rm = T)) %>%
  group_by(surveillance_type, lineage) %>%
  mutate(total_cases = sum(n_cases, na.rm = T))  %>%
  mutate(fraction_cases = n_cases / total_cases) %>%
  ungroup() %>%
  ggplot(aes(x = minimum_birth_year, y = fraction_cases,
             color = surveillance_type, group = surveillance_type)) +
  geom_line() +
  facet_grid(lineage~.) +
  scale_x_continuous(breaks = seq(1920,2020,20)) +
  theme(legend.position = c(0.02,0.9)) +
  scale_color_viridis(discrete = T, name = '', option = 'D') +
  xlab('Birth year') +
  ylab('Fraction of cases')


save_plot('../figures/byear_distribution_full.pdf', 
          byear_distribution_full,
          base_width = 8, base_height = 6)

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
  

aggregated_all_surveillance_distribution <- case_data_nz_all_surveillance %>%
  group_by(lineage, minimum_birth_year) %>%
  summarise(n_cases = sum(n_cases, na.rm = T)) %>%
  ungroup() %>% group_by(lineage) %>%
  mutate(total_cases = sum(n_cases)) %>%
  mutate(fraction_cases = n_cases/total_cases) %>%
  ungroup()

aggregated_all_surveillance_distribution_by_age <- case_data_nz_all_surveillance %>%
  group_by(lineage, cohort_value) %>%
  summarise(n_cases = sum(n_cases, na.rm = T)) %>%
  ungroup() %>% group_by(lineage) %>%
  mutate(total_cases = sum(n_cases)) %>%
  mutate(fraction_cases = n_cases/total_cases) %>%
  ungroup()

  
# Alternative plot per Katia's suggestion
aggregated_all_surveillance_distribution_histograms <- aggregated_all_surveillance_distribution %>%
  ggplot(aes(x = minimum_birth_year, y = fraction_cases)) +
  geom_col(data = aggregated_all_surveillance_distribution %>% filter(lineage == 'B/Victoria'),
             fill = 'darkorange1', alpha = 0.8) +
  geom_col(data = aggregated_all_surveillance_distribution %>% filter(lineage == 'B/Yamagata'),
           fill = 'mediumpurple3', alpha = 0.7) +
  scale_x_continuous(breaks = seq(1960,2020,5), limits = c(1959,2020),
                     expand = c(0.01,0)) +
  #scale_y_continuous(expand = c(0.005,0)) +
  theme(axis.text.x = element_text(size = 12),
        legend.position = 'None') +
  ylab('Fraction of cases') + xlab ('Year of birth') +
  scale_color_manual(name = '',
                     labels = c('B/Victoria','B/Yamagata'),
                     values = c('darkorange1','mediumpurple3'))

save_plot('../figures/byear_distribution_nz_allsurveillance.pdf',aggregated_all_surveillance_distribution_histograms,
          base_height = 4, base_width = 10)

aggregated_all_surveillance_distribution_histograms_by_age <- aggregated_all_surveillance_distribution_by_age %>%
  ggplot(aes(x = cohort_value, y = fraction_cases)) +
  geom_col(data = aggregated_all_surveillance_distribution_by_age %>% filter(lineage == 'B/Victoria'),
           fill = 'darkorange1', alpha = 0.8) +
  geom_col(data = aggregated_all_surveillance_distribution_by_age %>% filter(lineage == 'B/Yamagata'),
           fill = 'mediumpurple3', alpha = 0.7) +
  scale_x_continuous(expand = c(0.01,0), limits = c(-0.5,60)) +
  #scale_y_continuous(expand = c(0.005,0)) +
  theme(axis.text.x = element_text(size = 12),
        legend.position = 'None') +
  ylab('Fraction of cases') + xlab ('Age') +
  scale_color_manual(name = '',
                     labels = c('B/Victoria','B/Yamagata'),
                     values = c('darkorange1','mediumpurple3'))

save_plot('../figures/age_distribution_nz_allsurveillance.pdf',aggregated_all_surveillance_distribution_histograms_by_age,
          base_height = 4, base_width = 10)




# ======================= CORRELATIONS BETWEEN AGE/B. YEAR AND OBSERVATION YEAR =====================
uncounted_data <- as_tibble(case_data_nz_all_surveillance) %>% uncount(n_cases) %>%
  rename(year = observation_year, age = cohort_value) %>%
  select(country, year, lineage, age) %>%
  mutate(min_birth_year = year - age - 1)

joint_tibble %>% 
  rowwise() %>%
  mutate(period = ifelse(year <= 2010, '2001-2010','2011-2019')) %>%
  ungroup() %>%
  group_by(period, lineage, age) %>%
  summarise(n_cases = n()) %>%
  group_by(lineage, period) %>%
  mutate(fraction_cases = n_cases / sum(n_cases)) %>%
  ggplot(aes(x = age, y = fraction_cases)) +
  geom_col() +
  facet_grid(period~lineage, scales = 'free') +
  geom_smooth(se = F)

joint_tibble %>% 
  rowwise() %>%
  mutate(period = ifelse(year <= 2010, '2001-2010','2011-2019')) %>%
  ungroup() %>%
  group_by(period, lineage, min_birth_year) %>%
  summarise(n_cases = n()) %>%
  group_by(lineage, period) %>%
  mutate(fraction_cases = n_cases / sum(n_cases)) %>%
  ggplot(aes(x = min_birth_year, y = fraction_cases)) +
  geom_col() +
  facet_grid(period~lineage, scales = 'free') +
  geom_smooth(se = F)
  
age_vs_obs_year_CI_NZ <- lm(age~year:lineage, data = uncounted_data)
confint(age_vs_obs_year_CI_NZ, 'year:lineageB/Victoria', 0.95)
confint(age_vs_obs_year_CI_NZ, 'year:lineageB/Yamagata', 0.95)

birth_year_vs_obs_year_CI_NZ <-  lm(min_birth_year~year:lineage, data = uncounted_data)
confint(birth_year_vs_obs_year_CI_NZ, 'year:lineageB/Victoria', 0.95)
confint(birth_year_vs_obs_year_CI_NZ, 'year:lineageB/Yamagata', 0.95)

age_vs_obs_year_pl <- uncounted_data %>%
  ggplot(aes(x = year, y = age)) +
  geom_point(position = position_jitter(width = 0.1), alpha = 0.6) +
  facet_grid(.~lineage) + geom_smooth(method = 'lm') +
  xlab('Observation year') + ylab('Age') +
  scale_x_continuous(breaks = seq(2002,2020,2))

byear_vs_obs_year_pl <- uncounted_data %>%
  filter(country == 'New Zealand') %>%
  ggplot(aes(x = year, y = min_birth_year)) +
  geom_point(position = position_jitter(width = 0.1), alpha = 0.6) +
  facet_grid(.~lineage) + geom_smooth(method = 'lm') +
  xlab('Observation year') + ylab('Birth year')+
  scale_x_continuous(breaks = seq(2002,2018,2))

save_plot('../figures/age_and_byear_vs_obs_year.pdf',
          plot_grid(age_vs_obs_year_pl, byear_vs_obs_year_pl, nrow = 2),
          base_height = 9, base_width = 10)