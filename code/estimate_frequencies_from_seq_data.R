library(cowplot)
library(ggplot2)
theme_set(theme_cowplot())
library(reshape2)
library(dplyr)
library(lubridate)

output_directory <- '../results/processed_data/'

assign_region <- function(combined_isolate_data){
  regions <- tibble(country = c('United States','Australia','New Zealand'),
         region = c('United States', 'AUSNZ','AUSNZ'))
  return(left_join(combined_isolate_data,regions, by = 'country'))
}


main <- function(){
  genbank_B_data <- read.csv('../data/genbank_data/genbank_data.csv', header = T)
  gisaid_B_data <- read.csv('../data/gisaid_metadata/gisaid_metadata.csv')
  
  # Filter Genbank and Gisaid data to retain only isolates with known collection year and year >= 1970
  # Lower bound for split between Yam and Vic lineages
  genbank_B_data <- genbank_B_data %>% filter(is.na(year) == F, year >= 1970)
  gisaid_B_data <- gisaid_B_data %>% filter(is.na(year) == F, year >= 1970)
  
  # ------- Consolidade multiple occurrences of isolate/year(season)/country combinations -------
  # Note that years have been adjusted in the process_* scripts to represent seasons
  # And discard isolate/year/country combinations with inconsistent lineage assignments
  genbank_B_data <- genbank_B_data %>% group_by(isolate_name, country, year) %>% 
    # Find how many different lineages were assigned for each combination
    mutate(n_occurrences = n(), 
           n_lineage_assignments = length(unique(lineage)[is.na(unique(lineage)) == F])) %>%
    # Consolidate multiple occurrences
    # And discard isolate/year/country combinations with inconsistent lineage assignments
    filter(n_occurrences == 1, n_lineage_assignments < 2 ) %>%
    select(isolate_name, country, lineage, year) %>% ungroup()
  
  gisaid_B_data <- gisaid_B_data %>% group_by(isolate_name, country, year) %>% 
    # Find how many different lineages were assigned for each combination
    mutate(n_occurrences = n(), 
           n_lineage_assignments = length(unique(lineage)[is.na(unique(lineage)) == F])) %>%
    # Consolidade multiple occurrences
    # And discard isolate/year/country combinations with inconsistent lineage assignments
    filter(n_occurrences == 1, n_lineage_assignments < 2 ) %>%
    select(isolate_name, continent, country, lineage, year) %>% ungroup()
  
  # ----- Compare blast assignments with assignments from Genbank
  # Create tibble by merging...
  full_join(
    # ...Gisaid sub-tibble containing only isolates present in Genbank
    gisaid_B_data %>% filter(isolate_name %in% (genbank_B_data %>% pull(isolate_name))) %>% 
      select(isolate_name, lineage) %>% rename(gisaid_lineage = lineage),
    # ...Genbank sub-tiblle containing only isolates present in Gisaid
    genbank_B_data %>% filter(isolate_name %in% (gisaid_B_data %>% pull(isolate_name))) %>% 
      select(isolate_name, lineage) %>%
      rename(blast_lineage = lineage)
  ) %>%
    # Consider only gisaid isolates with lineage assignment
    filter(is.na(gisaid_lineage) == F) %>%
    group_by(gisaid_lineage,blast_lineage) %>% count() %>% ungroup %>% mutate(fraction = n/sum(n)) %>%
    filter(gisaid_lineage == blast_lineage) %>% summarise(agreement = sum(fraction))
  # -----
  
  # Merge tibbles by isolate name, country and year
  combined_isolate_data <- full_join(genbank_B_data %>% select(isolate_name, country, lineage, year),
                                     gisaid_B_data %>% select(isolate_name, country, lineage, year),
                                     by = c('isolate_name', 'country', 'year')) %>%
    # Remove isolates where lineage is not available in either dataset
    filter(is.na(lineage.x) == F | is.na(lineage.y) == F) %>%
    # Remove isolates where year and country information is missing
    filter(is.na(year) == F,is.na(country) == F)
  
  combined_isolate_data <- combined_isolate_data %>% 
    mutate(lineage.x = as.character(lineage.x),
           lineage.y = as.character(lineage.y)) %>%
    # Retain isolates where one OR the other database provides a lineage assignment, OR, if both do, retain only isolates for which they agree.
    filter((is.na(lineage.x) == T | is.na(lineage.y) == T) | (is.na(lineage.x) == F & is.na(lineage.y) == F & lineage.x == lineage.y)
    ) %>%
    # For isolates for which only one database provides lineage assignment, use that lineage assignment
    mutate(lineage = ifelse(is.na(lineage.x), lineage.y, lineage.x)) %>%
    select(isolate_name, country, lineage, year)
  
  # Add variable 'continent'
  combined_isolate_data <- left_join(combined_isolate_data,
                                     gisaid_B_data %>% select(country, continent) %>% 
                                       filter(!is.na(country), continent != '') %>% unique(),
                                     by = 'country')
  
  # Add variable 'region' ('United States','AUSNZ')
  combined_isolate_data <- assign_region(combined_isolate_data)
  
  # Calculate relative lineage frequencies in each year (season) pooled across all countries 
  lineage_frequency_data_global <- combined_isolate_data %>% 
    group_by(year) %>% 
    summarise(n_yamagata_global = sum(lineage == 'Yamagata'),
              n_victoria_global = sum(lineage == 'Victoria')) %>%
    mutate(year_total_global = n_victoria_global + n_yamagata_global,
           fraction_yamagata_global = n_yamagata_global / year_total_global, 
           fraction_victoria_global = 1 - fraction_yamagata_global) %>% ungroup() 
  
  # Aggregate frequencies excluding Australia and New Zealand
  lineage_frequency_data_global_minus_AUSNZ <- combined_isolate_data %>% 
    filter(country != 'Australia', country != 'New Zealand') %>%
    group_by(year) %>% 
    summarise(n_yamagata_global = sum(lineage == 'Yamagata'),
              n_victoria_global = sum(lineage == 'Victoria')) %>%
    mutate(year_total_global = n_victoria_global + n_yamagata_global,
           fraction_yamagata_global = n_yamagata_global / year_total_global, 
           fraction_victoria_global = 1 - fraction_yamagata_global) %>% ungroup() 
  
  # Calculate relative lineage frequencies in each year (season) in the United States
  lineage_frequency_data_US <- combined_isolate_data %>% 
    filter(country == "United States") %>%
    group_by(year, region) %>% 
    summarise(n_yamagata = sum(lineage == 'Yamagata'),
              n_victoria = sum(lineage == 'Victoria')) %>%
    mutate(year_total = n_victoria + n_yamagata,
           fraction_yamagata = n_yamagata / year_total, 
           fraction_victoria = 1 - fraction_yamagata) %>% ungroup()
  
  # Calculate relative lineage frequencies in Australia/NZ
  lineage_frequency_data_AUSNZ <- combined_isolate_data %>% 
    filter(region == 'AUSNZ') %>%
    group_by(year, region) %>% 
    summarise(n_yamagata = sum(lineage == 'Yamagata'),
              n_victoria = sum(lineage == 'Victoria')) %>%
    mutate(year_total = n_victoria + n_yamagata,
           fraction_yamagata = n_yamagata / year_total, 
           fraction_victoria = 1 - fraction_yamagata) %>% ungroup()
  
  # Add years not represented in each region (to be filled with U.S. or global data)
  all_years <-  seq(min(lineage_frequency_data_global$year), max(lineage_frequency_data_global$year))
  years_missing_US <- all_years[all_years %in% lineage_frequency_data_US$year == F]
  years_missing_AUSNZ <- all_years[all_years %in% lineage_frequency_data_AUSNZ$year == F]

  lineage_frequency_data_AUSNZ <- full_join(tibble(region = 'AUSNZ', year = years_missing_AUSNZ),
                                            lineage_frequency_data_AUSNZ) %>%
    arrange(year)
  
  lineage_frequency_data_US <- full_join(tibble(region = 'United States', year = years_missing_US),
                                            lineage_frequency_data_US) %>%
    arrange(year)
  
  # For years with fewer than 10 isolates, fill with data from all countries combined
  lineage_frequency_data_AUSNZ_adjusted <- left_join(lineage_frequency_data_AUSNZ,
                                                     lineage_frequency_data_global %>% 
                                                     select(year, fraction_yamagata_global, fraction_victoria_global, year_total_global),
                                                     by = 'year') %>%
    rowwise() %>%
    mutate(data_source = ifelse(year_total >=10 & is.na(year_total) == F, 'AUSNZ','all countries'),
           fraction_yamagata = ifelse(data_source == 'AUSNZ', fraction_yamagata, fraction_yamagata_global),
           fraction_victoria = ifelse(data_source == 'AUSNZ', fraction_victoria, fraction_victoria_global),
           year_total = ifelse(data_source == 'AUSNZ', year_total, year_total_global)) %>%
    ungroup()%>%
    select(-matches('_global'))
  

  # For computing purposes, duplicate AUSNZ rows into a set of Australia and a set of NZ rows
  lineage_frequency_data_AUSNZ_adjusted <- lineage_frequency_data_AUSNZ_adjusted %>% 
    rename(country = region) %>%
    mutate(country = ifelse(country == 'AUSNZ', 'Australia', country))
  lineage_frequency_data_AUSNZ_adjusted <- bind_rows(lineage_frequency_data_AUSNZ_adjusted,
                                 lineage_frequency_data_AUSNZ_adjusted %>% filter(country == 'Australia') %>%
                                   mutate(country = 'New Zealand')
  ) %>% arrange(year)

  # Assume that ancestral lineage had frequency 1 prior to 1988 (B intensity will be set to 0 prior to some origin year separately)
  split_year = 1988
  lineage_frequency_data_AUSNZ_adjusted <- full_join(as_tibble(expand.grid(country = c('Australia','United States','China','New Zealand'),
                                  year = 1900:(split_year-1))) %>%
                                    mutate(fraction_ancestor = 1),
                                  lineage_frequency_data_AUSNZ_adjusted %>%
                                    filter(year >= split_year)
                                  ) %>%
    rowwise() %>%
    mutate(fraction_yamagata = ifelse(year < split_year,0,fraction_yamagata),
           fraction_victoria = ifelse(year < split_year,0,fraction_victoria),
           fraction_ancestor = ifelse(year >= split_year,0, fraction_ancestor),
           data_source = ifelse(year < split_year, 'presumed_ancestor', data_source)) %>%
    ungroup()

  # Also return raw frequencies (without filling gaps with global data)
  lineage_frequency_data_raw <- bind_rows(lineage_frequency_data_AUSNZ, lineage_frequency_data_US) %>%
    rename(country = region) %>%
    mutate(country = ifelse(country == 'AUSNZ', 'Australia', country))
  lineage_frequency_data_raw <- bind_rows(lineage_frequency_data_raw,
                                      lineage_frequency_data_raw %>% filter(country == 'Australia') %>%
                                        mutate(country = 'New Zealand')
  ) %>% arrange(year)
  
  
  return(list(adjusted = lineage_frequency_data_AUSNZ_adjusted, raw = lineage_frequency_data_raw,
              global = lineage_frequency_data_global,
              global_minus_AUSNZ = lineage_frequency_data_global_minus_AUSNZ))
}
lineage_frequency_data <- main()
write.csv(lineage_frequency_data$adjusted, paste0(output_directory, 'lineage_frequencies_gisaid-genbank.csv'), row.names = F)
write.csv(lineage_frequency_data$raw, paste0(output_directory, 'lineage_frequencies_gisaid-genbank_raw.csv'), row.names = F)
write.csv(lineage_frequency_data$global, paste0(output_directory, 'lineage_frequencies_gisaid-genbank_global.csv'), row.names = F)


# ----------------------------- Plots with fraction Yamagata over time ----------------------
raw_freq_comparison <- bind_rows(lineage_frequency_data$raw %>%
                                   filter(country != 'New Zealand') %>%
                                   select(country, year, fraction_yamagata),
                                 lineage_frequency_data$global_minus_AUSNZ %>%
                                   mutate(country = 'All countries except Australia/New Zealand') %>%
                                   rename(fraction_yamagata = fraction_yamagata_global,
                                          fraction_victoria = fraction_victoria_global) %>%
                                   select(country, year, fraction_yamagata)) %>% 
  arrange(year) %>%
  mutate()

# Plot with adjusted fraction
yam_freq_ausnz <- ggplot(assign_region(lineage_frequency_data$adjusted) %>%
                                       mutate(region = ifelse(region == 'AUSNZ','Australia / New Zealand',region)) %>%
                                       filter(year >=1984, country == 'Australia',
                                              data_source!= 'presumed_ancestor'),
                                     aes(x = year, y = fraction_yamagata, label = year_total)) +
  geom_line(alpha = 0.5) + geom_point(shape = 21, size = 6, (aes(fill = data_source))) +
  geom_text(size = 2) +
  xlab('Season') + ylab('Frequency of B/Yamagata') +
  scale_x_continuous(breaks = seq(1970,2020,5), limits = c(1970,2020),
                     expand = c(0.01,0)) +
  #facet_grid(region~.) +
  scale_fill_brewer(type = 'qual', palette = 3, name = 'Data source:',
                    labels = c('All represented countries', 'Australia/New Zealand')) +
  theme(legend.position = c(0.02,0.9)) +
  geom_vline(xintercept = c(1970,1983), linetype = 2, color = 'gray30')  +
  geom_text(data = tibble(year = 1976.5,
                          #region = c('Australia / New Zealand', 'China', 'United States'),
                          #region = c('Australia / New Zealand','United States'),
                          fraction_yamagata = 0.55,
                          #text = c('','Approximate interval of\nB/Vic - B/Yam split','')),
                          text = c('Estimates of\nB/Vic - B/Yam split','')),
            aes(label = text))

save_plot('../figures/lineage_frequencies/lineage_frequencies_gisaid-genbank_ausnz.pdf',
          yam_freq_ausnz,base_width = 12, base_height = 4)


AUSNZ_VS_US <- left_join(lineage_frequency_data$raw %>% filter(country == 'Australia'),
                         lineage_frequency_data$raw %>% filter(country == 'United States')%>%
                           rename(fraction_yamagata_US = fraction_yamagata) %>%
                           select(year, fraction_yamagata_US), by = 'year') %>%
  filter(year >= 0) 
cor.test(AUSNZ_VS_US$fraction_yamagata,AUSNZ_VS_US$fraction_yamagata_US)

AUSNZ_VS_US_pl <- AUSNZ_VS_US %>%
  ggplot(aes(x = fraction_yamagata_US, y = fraction_yamagata)) +
  geom_text(aes(label = year), size = 3) +
  xlab('Frequency of B/Yamagata\nin the United States') +
  ylab('Frequency of B/Yamagata in Australia\nand New Zealand combined') +
  theme(axis.text = element_text(size = 14),
        axis.title = element_text(size = 16)) +
  scale_x_continuous(breaks = seq(0,1,0.2)) +
  scale_y_continuous(breaks = seq(0,1,0.2),  limits = c(0,1)) +
  geom_smooth(method = 'lm')

AUSNZ_VS_global <- left_join(lineage_frequency_data$raw %>% 
                               filter(country == 'Australia'),
          lineage_frequency_data$global_minus_AUSNZ %>%
            select(year, fraction_yamagata_global, fraction_victoria_global, year_total_global),
          by = 'year') %>%
  filter(year >= 0) 
cor.test(AUSNZ_VS_global$fraction_yamagata,AUSNZ_VS_global$fraction_yamagata_global)

AUSNZ_VS_global_pl <- AUSNZ_VS_global %>%
  ggplot(aes(x = fraction_yamagata_global, y = fraction_yamagata)) +
  geom_text(aes(label = year), size = 3)  +
  ylab('Frequency of B/Yamagata in Australia\nand New Zealand combined') +
  xlab('Frequency of B/Yamagata in all\ncountries in dataset except Australia and New Zealand') +
  theme(axis.text = element_text(size = 14),
        axis.title = element_text(size = 16)) +
  scale_x_continuous(breaks = seq(0,1,0.2)) +
  scale_y_continuous(breaks = seq(0,1,0.2), limits = c(0,1)) +
  geom_smooth(method = 'lm')

save_plot('../figures/lineage_frequencies/global_vs_local_freqs.pdf',
          plot_grid(AUSNZ_VS_global_pl, AUSNZ_VS_US_pl, nrow = 1),
          base_width = 15, base_height = 6)




 
