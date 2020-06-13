library(tidyverse)
library(cowplot)
library(reshape2)
library(stringr)


# Read B HA sequence data from Genbank:
B_seq_data <- read.csv('../data/genbank_data/genbank_B_HA_seqs.csv', header =T)
B_seq_data <- as.tibble(B_seq_data)

# Use Gisaid B data to map countries to continents and resolve country name inconsistencies
source('process_gisaid_metadata.R')
country_and_continents <- gisaid_B_data %>% select(country, continent) %>% unique()

# Read output of blastn against Vic/87 and Yam/88 strains:
blast_output <- read.csv('../results/B_seq_blast/B_seq_blast.csv', header =T)
blast_output <- as.tibble(blast_output)

# Process output into tibble with the accn # of each strain and the most lik. lineage (Yam or Vic) based on bit score
blast_output <- blast_output %>% group_by(query_id) %>% 
  filter(bit_score == max(bit_score)) %>% 
  select(query_id, subject_id) %>% 
  mutate(lineage = ifelse(grepl('Victoria',subject_id), 'Victoria', 'Yamagata')) %>% 
  select(query_id, lineage) %>% mutate(accession_number = str_extract(query_id,'[^_]+')) %>%
  mutate(accession_number = paste('>', accession_number, sep = '')) %>% 
  ungroup() %>%
  select(accession_number, lineage)

# Merge tibbles
blast_output$accession_number = factor(blast_output$accession_number, 
                                       levels = levels(B_seq_data$accession_number))
genbank_B_data <- full_join(B_seq_data, blast_output, by = 'accession_number') %>%
  # Rename "strain" variable as "isolate_name" for consistency with Gisaid tibble
  rename(isolate_name = strain)

genbank_B_data <- mutate(genbank_B_data, country = as.character(country))  
# Resolve country name inconsistencies
for(country in unique(genbank_B_data$country)){
  country <- ifelse(country == '',NA,country)
  standard_name <- switch(country,
                          'South Korea' = 'Korea, Republic of',
                          'Viet Nam' = 'Vietnam',
                          'Bolivia' = "Bolivia, Plurinationial State of",
                          'Hong Kong' = "Hong Kong (SAR)",
                          'Russia' = "Russian Federation",
                          'Czechoslovakia' =  "Czech Republic",
                          'Macau' = 'Macao',
                          'Venezuela' = 'Venezuela, Bolivarian Republic of',
                          'Iran' = 'Iran, Islamic Republic of',
                          'USA' = 'United States'
  )
  if(!is.null(standard_name)){
    genbank_B_data$country[genbank_B_data$country == country] <- standard_name
  }
}

# Add continent from Gisaid data
genbank_B_data <- left_join(genbank_B_data, country_and_continents, by = 'country')
genbank_B_data$continent[genbank_B_data$country == 'Zimbabwe'] <- 'Africa'

genbank_B_data %>% filter(continent %in% c('Europe','North America') | country == 'China', !is.na(month)) %>% 
  group_by(year) %>% 
  summarise(previous_season = sum(month < 10), next_season = sum(month >=10)) %>%
  mutate(fraction_previous_season = previous_season / (previous_season + next_season)) %>%
  ungroup() %>% 
  summarise(median(fraction_previous_season))

# Convert year to season based on continent
genbank_B_data <- genbank_B_data %>%
  rowwise() %>%
  mutate(year = ifelse(continent %in% c('Europe','North America') | country == 'China' & 
                         (month < 10 | is.na(month)),
                       year - 1, year)) %>%
  ungroup()

# Plot with # of influenza B sequences on Genbank (1980-2017)
# genbank_B_seqs_vs_time_pl <- genbank_B_data %>%
#   filter(year >= 1980, country %in% c('United States', 'Australia', 'New Zealand')) %>%
#   group_by(year, lineage, country) %>%
#   count() %>%
#   ggplot(aes(x = year, y = n, label = n, color = lineage)) +
#   geom_line(aes(color = lineage)) +
#   geom_text(size = 3) +
#   facet_grid(country~.) +
#   ylab('Number of influenza B HA sequences \non Genbank (1980-2017)') +
#   xlab('Year')
# save_plot('../figures/genbank_B_seqs_over time.pdf',
#           genbank_B_seqs_vs_time_pl, base_width = 8, base_height = 6)

write.csv(genbank_B_data, '../data/genbank_data/genbank_data.csv', row.names = F)

