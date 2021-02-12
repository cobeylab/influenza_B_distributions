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
    
    modified_tibble <- left_join(combined_isolate_data,regions, by = 'country')
    
    modified_tibble <- modified_tibble %>%
      mutate(region = ifelse(continent == 'Europe', 'Europe', region))
    
    east_asia_countries <- c('China', 'Hong Kong (SAR)', 'Korea, Republic of','Japan',
                             'Mongolia','Taiwan')
    
    modified_tibble <- modified_tibble %>%
      mutate(region = ifelse(country %in% east_asia_countries, 'East Asia', region))
    
  return(modified_tibble)
}


# Base dplyr pipeline (input object will already be subset and grouped depending on the region)
estimate_lineage_freqs <- function(combined_isolate_data){
  combined_isolate_data %>% 
    summarise(n_yamagata = sum(lineage == 'Yamagata'),
            n_victoria = sum(lineage == 'Victoria')) %>%
    mutate(year_total = n_victoria + n_yamagata,
           fraction_yamagata = n_yamagata / year_total, 
           fraction_victoria = 1 - fraction_yamagata) %>% 
    ungroup() 
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
  
  # Add variable 'region' 
  combined_isolate_data <- assign_region(combined_isolate_data)
  
  # Calculate relative lineage frequencies in each year (season) pooled across all countries 
  lineage_frequency_data_global <- estimate_lineage_freqs(combined_isolate_data %>% group_by(year)) %>%
    rename(n_yamagata_global = n_yamagata, n_victoria_global = n_victoria,
           year_total_global = year_total, fraction_yamagata_global = fraction_yamagata,
           fraction_victoria_global = fraction_victoria)
  
  # Aggregate frequencies across all countries excluding Australia and New Zealand
  lineage_frequency_data_global_minus_AUSNZ <- 
    estimate_lineage_freqs(combined_isolate_data %>% filter(country != 'Australia', country != 'New Zealand') %>%
                             group_by(year)) %>%
    rename(n_yamagata_global = n_yamagata, n_victoria_global = n_victoria,
           year_total_global = year_total, fraction_yamagata_global = fraction_yamagata,
           fraction_victoria_global = fraction_victoria)
  
  # Calculate relative lineage frequencies in each year (season) in the United States
  lineage_frequency_data_US <- estimate_lineage_freqs(combined_isolate_data %>% 
                                                        filter(country == "United States") %>%
                                                        group_by(year, region))
  
  # Calculate relative lineage frequencies in Australia/NZ
  # Adding Europe and East Asia to address reviewers comments
  lineage_frequency_data_AUSNZ <- estimate_lineage_freqs(combined_isolate_data %>%
                                                           filter(region %in% c('AUSNZ','Europe', 'East Asia')) %>%
                                                           group_by(year, region))
  

  # Add years not represented in each region (to be filled with U.S. or global data)
  lineage_frequency_data_AUSNZ <- left_join(tibble(expand.grid(region = unique(lineage_frequency_data_AUSNZ$region),
                                                               year = seq(min(lineage_frequency_data_AUSNZ$year),
                                                                          max(lineage_frequency_data_AUSNZ$year)))),
            lineage_frequency_data_AUSNZ)
  
  lineage_frequency_data_US <- left_join(tibble(year = seq(min(lineage_frequency_data_US$year),
                                                                   max(lineage_frequency_data_US$year))),
                                     lineage_frequency_data_US)
  
  # For years with fewer than 10 isolates, fill with data from all countries combined
  lineage_frequency_data_AUSNZ_adjusted <- left_join(lineage_frequency_data_AUSNZ,
                                                     lineage_frequency_data_global %>% 
                                                     select(year, fraction_yamagata_global, fraction_victoria_global, year_total_global),
                                                     by = 'year') %>%
    rowwise() %>%
    mutate(data_source = ifelse(year_total >=10 & is.na(year_total) == F, region,'all countries'),
           fraction_yamagata = ifelse(data_source != 'all countries', fraction_yamagata, fraction_yamagata_global),
           fraction_victoria = ifelse(data_source != 'all countries', fraction_victoria, fraction_victoria_global),
           year_total = ifelse(data_source != 'all countries', year_total, year_total_global)) %>%
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
  lineage_frequency_data_AUSNZ_adjusted <- full_join(as_tibble(expand.grid(country = unique(lineage_frequency_data_AUSNZ_adjusted$country),
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


# Lineage frequencies assuming no Vic from 1988 to 2000 in Europe, New Zealand and Australia
lineage_frequencies_gisaid_genbank_noVicin1990s <- lineage_frequency_data$adjusted %>%
  mutate(fraction_yamagata = ifelse(year>= 1988 & year <= 2000 & country != 'East Asia' , 1, fraction_yamagata),
         fraction_victoria = ifelse(year>= 1988 & year <= 2000 & country != 'East Asia', 0, fraction_victoria),
         n_yamagata = ifelse(year>= 1988 & year <= 2000 & country != 'East Asia', NA, n_yamagata),
         n_victoria = ifelse(year>= 1988 & year <= 2000 & country != 'East Asia', NA, n_victoria),
         year_total = ifelse(year>= 1988 & year <= 2000 & country != 'East Asia', NA, year_total),
         data_source = ifelse(year>= 1988 & year <= 2000 & country != 'East Asia',
                              'assuming_no_vic_in_1990s', data_source))

# 

write.csv(lineage_frequencies_gisaid_genbank_noVicin1990s,
          paste0(output_directory, 'lineage_frequencies_gisaid-genbank_noVicin1990s.csv'), row.names = F)

# Lineage frequencies assuming Vic was dominant from 1983 to 1990, then Yam was dominant
lineage_frequencies_gisaid_genbank_Vicin80s_Yamin90s <- lineage_frequencies_gisaid_genbank_noVicin1990s %>%
  mutate(fraction_yamagata = ifelse(year>= 1983 & year < 1990 & country != 'East Asia' , 0, fraction_yamagata),
         fraction_victoria = ifelse(year>= 1983 & year < 1990 & country != 'East Asia', 1, fraction_victoria),
         fraction_ancestor = ifelse(year>= 1983 & year < 1990 & country != 'East Asia', 0, fraction_ancestor),
         n_yamagata = ifelse(year>= 1983 & year < 1990 & country != 'East Asia', NA, n_yamagata),
         n_victoria = ifelse(year>= 1983 & year < 1990 & country != 'East Asia', NA, n_victoria),
         year_total = ifelse(year>= 1983 & year < 1990 & country != 'East Asia', NA, year_total),
        data_source = ifelse(year>= 1983 & year <= 2000 & country != 'East Asia',
                     'assuming_vic_in_80s_yam_in_90s', data_source))

write.csv(lineage_frequencies_gisaid_genbank_Vicin80s_Yamin90s,
          paste0(output_directory, 'lineage_frequencies_gisaid-genbank_Vicin80s_Yamin90s.csv'), row.names = F)

# ----------------------------- Plots with fraction Yamagata over time ----------------------
# Plot with global frequencies in the 1990s
yam_freq_ausnz <- ggplot(lineage_frequency_data$adjusted %>% filter(year >=1984, country == 'New Zealand',
                                              data_source!= 'presumed_ancestor'),
                                     aes(x = year, y = fraction_yamagata, label = year_total)) +
  geom_line(alpha = 0.5) + geom_point(shape = 21, size = 6, aes(fill = data_source)) +
  geom_text(size = 2, aes(color = data_source), show.legend = F) +
  xlab('Season') + ylab('Frequency of B/Yamagata') +
  scale_x_continuous(breaks = seq(1970,2020,5), limits = c(1970,2020),
                     expand = c(0.01,0)) +
  #facet_grid(region~.) +
  scale_fill_brewer(type = 'qual', palette = 3, name = 'Data source:',
                    labels = c('All represented countries', 'New Zealand/Australia')) +
  theme(legend.position = c(0.02,0.9)) +
  geom_vline(xintercept = c(1983), linetype = 2, color = 'gray30')  +
  geom_text(data = tibble(year = 1983,
                          #region = c('Australia / New Zealand', 'China', 'United States'),
                          #region = c('Australia / New Zealand','United States'),
                          fraction_yamagata = 0.55,
                          #text = c('','Approximate interval of\nB/Vic - B/Yam split','')),
                          text = c('Estimate of\nB/Vic - B/Yam split','')),
            aes(label = text)) +
  scale_color_manual(values = c('black','white'))

save_plot('../figures/lineage_frequencies/lineage_frequencies_global_1990s.pdf',
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




 
