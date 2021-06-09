library(dplyr)
library(tidyr)
library(readr)
library(ggplot2)
library(cowplot)
library(patchwork)
theme_set(theme_cowplot())


isolate_distribution_eu_linfreqs <- read_csv('../results/processed_data/europe_lineage_freq_isolate_distribution.csv')
isolate_distribution_eu_agedist <- read_csv('../results/processed_data/europe_age_dist_isolate_distribution.csv')

isolate_distribution_east_asia_linfreqs <- read_csv('../results/processed_data/east_asia_lineage_freq_isolate_distribution.csv')

isolate_distribution_eu <- bind_rows(isolate_distribution_eu_linfreqs,
          isolate_distribution_eu_agedist) %>%
  mutate(country_group = case_when(
    country %in% c('Spain','Portugal') ~ 'Spain and Portugal',
    country %in% c('United Kingdom', 'Ireland') ~ 'UK and Ireland',
    country == 'Russian Federation' ~ 'Russia',
    country == 'France' ~ 'France',
    country == 'Italy' ~ 'Italy',
    country == 'Germany' ~ 'Germany',
    country == 'Austria' ~ 'Austria',
    country %in% c('Sweden','Norway','Dernmark', 'Finland') ~ 'Scandinavia',
    country == 'Netherlands' ~ 'Netherlands',
    country == 'Greece' ~ 'Greece',
    TRUE ~ 'Other')) 


country_levels <- as.character(isolate_distribution_eu %>%
  group_by(country_group) %>%
  summarise(total_country_group_isolates = sum(n_isolates)) %>%
  arrange(desc(total_country_group_isolates)) %>%
  filter(country_group != 'Other') %>%
  pull(country_group))

country_levels <- c(country_levels,'Other')

isolate_distribution_eu <- isolate_distribution_eu  %>%
  mutate(country_group = factor(country_group, levels = country_levels))

year_totals_eu <- isolate_distribution_eu %>%
  select(year, total_isolates, data_type) %>%
  unique()

isolate_distribution_eu_pl <- isolate_distribution_eu %>%
  ggplot(aes(x = year, y = fraction)) +
  geom_col(aes(fill = country_group)) +
  facet_wrap('data_type') +
  geom_text(data = year_totals_eu, aes(x = year, y = 1.03, label = total_isolates),
            angle = 90, size = 3,show.legend = F) +
  background_grid() +
  scale_fill_brewer(name = 'Country', type = 'div',
                    palette = 4) +
  xlab('Year') +
  ylab('Fraction')



year_totals_east_asia <- isolate_distribution_east_asia_linfreqs %>%
  select(year, total_isolates) %>% unique()


isolate_distribution_east_asia_pl <- isolate_distribution_east_asia_linfreqs %>%
  ggplot(aes(x = year, y = fraction)) +
  geom_col(aes(fill = country)) +
  facet_wrap('data_type') +
  geom_text(data = year_totals_east_asia, aes(x = year, y = 1.03, label = total_isolates),
            angle = 90, size = 3,show.legend = F) +
  background_grid() +
  scale_fill_brewer(name = 'Country', type = 'div',
                    palette = 4) +
  xlab('Year') +
  ylab('Fraction') +
  theme(legend.position = c(1,0.5)) 


combined_plot <- (isolate_distribution_eu_pl) / 
  (isolate_distribution_east_asia_pl  | plot_spacer())


save_plot('../figures/isolate_data_by_country_EU_east_asia.pdf',
          combined_plot, 
          base_height = 10.5, base_width = 16)
          
