library(cowplot)
library(ggplot2)
theme_set(theme_cowplot())
library(reshape2)
library(dplyr)
library(viridis)

# Raw Gisaid / Genbank frequencies for Aus/NZ and the U.S.
gisaid_genbank_frequencies_local <- as_tibble(read.csv('../results/processed_data/lineage_frequencies_gisaid-genbank_raw.csv',
                                                 header = T)) %>%
  mutate(type = 'genetic',
         source = 'Gisaid/Genbank (local)') %>%
  select(year, year_total, country, fraction_yamagata, fraction_victoria, type, source)


# Gisaid / Genbank frequencies estimated from all countries when data was scarce
gisaid_genbank_frequencies_global <- as_tibble(read.csv('../results/processed_data/lineage_frequencies_gisaid-genbank_global.csv',
                                                 header = T)) %>%
  rename(fraction_yamagata = fraction_yamagata_global,
         fraction_victoria = fraction_victoria_global,
         year_total = year_total_global) %>%
  mutate(type = 'genetic',
         country = 'All countries',
         source = 'Gisaid/Genbank (All countries)') %>%
  select(year,year_total, country, fraction_yamagata, fraction_victoria, type, source)

# Gisaid / Genbank frequencies estimated from all countries when data was scarce, but with Yam freq = 100% in the 90s.
gisaid_genbank_frequencies_global_noVicin1990s <- as_tibble(read.csv('../results/processed_data/lineage_frequencies_gisaid-genbank_noVicin1990s.csv',
                                                        header = T)) %>%
  mutate(type = 'genetic',
         country = 'All countries',
         source = 'Gisaid/Genbank (All countries + no Vic in 1990s)') %>%
  select(year,year_total, country, fraction_yamagata, fraction_victoria, type, source)


gisaid_genbank_frequencies <- bind_rows(gisaid_genbank_frequencies_local, gisaid_genbank_frequencies_global)
gisaid_genbank_frequencies <- bind_rows(gisaid_genbank_frequencies, gisaid_genbank_frequencies_global_noVicin1990s)


surveillance_frequencies <- bind_rows(
  as_tibble(read.csv('../data/surveillance_report_frequencies/Aus_lineage_freqs.csv')),
  as_tibble(read.csv('../data/surveillance_report_frequencies/NZ_lineage_freqs.csv')),
  as_tibble(read.csv('../data/surveillance_report_frequencies/US_lineage_freqs.csv'))
) %>% 
  rename(year_total = n_specimens_tested) %>%
  select(year, year_total, country, fraction_yamagata, fraction_victoria, method) %>%
  mutate(source = 'surveillance reports') %>%
  rename(type = method) %>%
  mutate(type = as.character(type))

lineage_frequencies <- bind_rows(gisaid_genbank_frequencies, surveillance_frequencies)

# ---------------------------------------------------- Plots ---------------------------------------------------
surveillance_vs_gisaid_genbank_freqs_pl <- ggplot(lineage_frequencies %>%
                                                    filter(is.na(fraction_yamagata) == F,
                                                           country != 'All countries',
                                                           year > 1983),
       aes(x = year, y = fraction_yamagata)) +
  geom_point(aes(color = factor(source)), alpha = 0.8) +
  geom_text(data = lineage_frequencies %>%
              filter(is.na(fraction_yamagata) == F, source == 'Gisaid/Genbank (local)', year > 1983),
                     aes(label = year_total), color = 'blue',
            position = position_nudge(x = 0, y = 0.1), size = 3) +
  geom_text(data = lineage_frequencies %>%
              filter(is.na(fraction_yamagata) == F, source == 'surveillance reports', year > 1983),
            aes(label = year_total), color = 'orange',
            position = position_nudge(x = , y = -0.1), size = 3) +
  geom_line(aes(color = factor(source)), alpha = 0.8) +
  facet_grid(country~.) +
  #ggtitle('Note: Gisaid/Genbank frequencies were estimated jointly for Australia and New Zealand') +
  xlab('Year (season start)') +
  ylab('Frequency of B/Yamagata') +
  scale_x_continuous(breaks = seq(1982,2020,2), limits = c(1982, 2020)) +
  scale_color_manual(values = c('blue','orange'), name = 'Source',
                     labels = c('Gisaid/Genbank', 'surveillance reports')) +
  theme(legend.position = 'top', axis.text.x = element_text(size = 9))

save_plot('../figures/lineage_frequencies/surveillance_vs_gisaid_genbank_freqs.pdf',
          surveillance_vs_gisaid_genbank_freqs_pl,
          base_height = 5, base_width = 10)

lineage_frequencies %>% filter(source == 'Gisaid/Genbank (All countries + no Vic in 1990s)')

# Plot with lineage frequencies used in the model, assuming all Yam in the 90s
lineage_frequencies %>% filter(source == 'Gisaid/Genbank (local)', country == 'New Zealand', year >= 1988)

local_points <- lineage_frequencies %>% filter(source == 'Gisaid/Genbank (local)', country == 'New Zealand', year >= 1990) %>%
  filter(!is.na(year_total)) %>%
  filter((year <= 2000) | (year > 2000 & year_total >= 10)) %>%
  mutate(source = 'New Zealand and Australia isolates')

global_points <- lineage_frequencies %>% filter(year >= 2001, source == "Gisaid/Genbank (All countries + no Vic in 1990s)") %>%
  filter((year %in% unique(local_points$year)) == F) %>% unique() %>%
  mutate(source = 'Isolates from all represented countries')

points <- bind_rows(local_points, global_points) %>% arrange(year)


lineage_frequencies_gisaid_genbank_noVicin1990s <- ggplot(lineage_frequencies %>% filter(source == "Gisaid/Genbank (All countries + no Vic in 1990s)",
                                      year >= 1988), aes(x = year, y = fraction_yamagata)) +
  geom_line(alpha = 0.5) + 
  geom_point(data = points, shape = 21, size = 6, (aes(fill = source))) +
  geom_text(data = points, size = 2, aes(label = year_total)) +
  xlab('Season') + ylab('Frequency of B/Yamagata') +
  scale_x_continuous(breaks = seq(1960,2020,5), limits = c(1959,2020),
                     expand = c(0.01,0)) +
  #facet_grid(region~.) +
  scale_fill_brewer(type = 'qual', palette = 3, name = 'Frequency estimates from:',
                    labels = c('All represented countries', 'Australia/New Zealand')) +
  theme(legend.position = c(0.02,0.9),
        legend.background = element_rect(fill = 'white'),
        legend.text = element_text(size = 10),
        legend.title = element_text(size = 10)) +
  geom_vline(xintercept = c(1970,1983), linetype = 2, color = 'gray30')  +
  geom_text(data = tibble(year = 1976.5,
                          #region = c('Australia / New Zealand', 'China', 'United States'),
                          #region = c('Australia / New Zealand','United States'),
                          fraction_yamagata = 0.55,
                          #text = c('','Approximate interval of\nB/Vic - B/Yam split','')),
                          text = c('Estimates of\nB/Vic - B/Yam split','')),
            aes(label = text))

save_plot('../figures/lineage_frequencies/lineage_frequencies_gisaid-genbank_ausnz_noVicin1990s.pdf',
          lineage_frequencies_gisaid_genbank_noVicin1990s,base_width = 10, base_height = 4)


# Plot with antigenic frequencies overlaid across countries
# lineage_frequencies_surveillance_raw_pl <- ggplot(lineage_frequencies %>%
#                                          filter(is.na(fraction_yamagata) == F,
#                                                 source == 'surveillance reports'),
#                                        aes(x = year, y = fraction_yamagata)) +
#   geom_point(aes(color = factor(country)), alpha = 0.9) +
#   geom_line(aes(color = factor(country)), alpha = 0.9)  +
#   xlab('Year (season start)') +
#   ylab('Frequency of B/Yamagata') +
#   scale_x_continuous(breaks = seq(1980,2020,2)) +
#   scale_color_brewer(type = 'qual', palette = 2, name = 'Country')
# 
# save_plot('../figures/lineage_frequencies/lineage_frequencies_surveillance_raw.pdf',
#           lineage_frequencies_surveillance_raw_pl,
#           base_height = 4, base_width = 14)


