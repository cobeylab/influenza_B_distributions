# Plot imprinting probabilities
library(dplyr)
library(reshape2)
library(stringr)
library(ggplot2)
library(cowplot)
theme_set(theme_cowplot())

library(viridis)
library(rlang)

args = commandArgs(trailingOnly = T)

intensity_scores_path = args[1] # intensity_scores_path = '../results/processed_data/intensity_scores.csv'
lineage_frequencies_path = args[2] # lineage_frequencies_path = '../results/processed_data/lineage_frequencies_gisaid-genbank_noVicin1990s.csv'
season_incidence_curves_path = args[3] # season_incidence_curves_path = '../results/processed_data/season_incidence_curves.csv'
demographic_data_path = args[4] # demographic_data_path = '../results/processed_data/demographic_data.csv'
birth_year_range = args[5] #String with comma-separated values for the first and last birth year to consider
birth_year_range = as.numeric(str_split(birth_year_range, ',')[[1]])
observation_year = as.numeric(args[6])
pars_file <- as.character(args[7]) 
model <- as.character(args[8]) 

parent_directory = dirname(pars_file)

# Load model functions
source('model_functions.R')

# Load fixed parameters (baseline attack rate, maternal ab protection, school start age)
source('basic_parameters.R')
# Load functions for calculating imprinting probabilities
source('imprinting_probabilites.R')

# Read intensity scores, lineage frequencies and season incidence curves
intensity_scores <- as_tibble(read.csv(intensity_scores_path))
lineage_frequencies <- as_tibble(read.csv(lineage_frequencies_path))
season_incidence_curves <- as_tibble(read.csv(season_incidence_curves_path))
demographic_data <- as_tibble(read.csv(demographic_data_path)) %>%
  filter(region == 'AUSNZ') %>%
  mutate(birth_year = observation_year - cohort_value) %>%
  filter(birth_year >= birth_year_range[1], birth_year <= birth_year_range[2],
         observation_year == !!observation_year) %>%
  group_by(country) %>%
  mutate(total_country_pop = sum(n_persons)) %>%
  ungroup() %>%
  mutate(pop_fraction = n_persons / total_country_pop)

stopifnot(all(abs(demographic_data %>% group_by(country) %>% summarise(S = sum(pop_fraction)) %>% pull(S) - 1) < 1e-6))


# Read parameters
pars <- as_tibble(read.csv(pars_file)) %>% mutate(par = as.character(par))
if("constrained_values" %in% names(pars)){
  # If parameter file is for a constrained fit, adjust:
  pars <- pars %>% select(-global_mle) %>% rename(value = constrained_values)
}
if(model == 'no_ancestor'){
  gamma_AV <- 0
  gamma_AY <- 0
}


# Replace NA values for irrelevant parameters with 0
# (e.g. reporting factor in U.S. when only looking at AUS/NZ)
pars[is.na(pars)] <- 0

# Define objects for each parameter
for(i in 1:nrow(pars)){
  eval(parse(text = paste0(pars[i,1],'=',pars[i,2])))
}
chi_VY <- chi_Y * gamma_VY
chi_YV <- chi_V * gamma_YV
chi_AV <- chi_V * gamma_AV
chi_AY <- chi_Y * gamma_AY


# Tibble with combinations of country/birth year for the specified observation year
years_tibble <- as_tibble(expand.grid(c('Australia','New Zealand'),
                                      observation_year, max(birth_year_range[1], birth_year_cutoff):birth_year_range[2])) %>%
  rename(country = Var1, observation_year = Var2, birth_year = Var3)


iprobs_by_byear <- mapply(FUN = calculate_iprobs_byear,
                          birth_year = years_tibble$birth_year,
                          country_name = as.character(years_tibble$country),
                          MoreArgs = list(max_obs_year = observation_year,
                                          min_obs_year = observation_year,
                                          beta1 = beta1,
                                          beta2 = beta2,
                                          beta3 = beta3,
                                          chi_VY = chi_VY,
                                          chi_YV = chi_YV,
                                          chi_AV = chi_AV,
                                          chi_AY = chi_AY,
                                          cutoff_age = cutoff_age,
                                          maternal_ab_duration = maternal_ab_duration,
                                          lineage_frequencies = lineage_frequencies,
                                          intensity_scores = intensity_scores,
                                          season_incidence_curves = season_incidence_curves,
                                          school_start_age = school_start_age,
                                          oldest_atk_rate_age = oldest_atk_rate_age),
                          SIMPLIFY = F
)

iprobs_by_byear <- bind_rows(iprobs_by_byear) %>%
  filter(observation_year == !!observation_year)

# If range of birth years extends before birth_year_cutoff, copy history probabilities at birth_year_cutoff
if(birth_year_range[1] < birth_year_cutoff){
  for(y in birth_year_range[1]:(birth_year_cutoff-1)){
    iprobs_by_byear <- bind_rows(iprobs_by_byear,
                                 iprobs_by_byear %>% filter(birth_year == birth_year_cutoff, observation_year != birth_year) %>%
                                   mutate(birth_year = y)) 
  }
}

suscep_yam <- calculate_relative_susceptibility(
  iprobs_by_byear %>% 
    mutate(lineage = 'B/Yamagata') %>%
    mutate(R_V, R_Y, chi_V, chi_Y, chi_VY, chi_YV,chi_AV, chi_AY)) %>%
  rename(rel_suscep_yam = rel_susceptibility) %>%
  select(country, observation_year, birth_year, rel_suscep_yam)

suscep_vic <- calculate_relative_susceptibility(
  iprobs_by_byear %>% 
    mutate(lineage = 'B/Victoria') %>%
    mutate(R_V, R_Y, chi_V, chi_Y, chi_VY, chi_YV,chi_AV, chi_AY)) %>%
  rename(rel_suscep_vic = rel_susceptibility) %>%
  select(country, observation_year, birth_year, rel_suscep_vic)

iprobs_by_byear <- left_join(iprobs_by_byear,
                             suscep_vic, by = c('country','birth_year','observation_year'))
iprobs_by_byear <- left_join(iprobs_by_byear,
                             suscep_yam, by = c('country','birth_year','observation_year'))


# Plot probabilities of infection histories
history_order = c('P_A0','P_AY0','P_AV0','P_AVY','P_Y0','P_YV','P_VY','P_V0','P_0')
history_labels = c()
for(h in history_order){
  subscript = strsplit(h, '_')[[1]][2]
  if(subscript == 'AVY'){
    subscript = 'A,{VY}'
  }else{
    if(subscript == 'A0'){
      subscript = 'A,0,0'
    }else{
      subscript = paste(strsplit(subscript,'')[[1]],collapse = ',')
    }
  }
  
  history_labels = c(history_labels,
                     bquote(.(subscript)))
}

inf_history_probs_pl <- iprobs_by_byear %>% 
  select(country, observation_year, birth_year, matches('P_', ignore.case = F)) %>%
   melt(id = c('observation_year','birth_year','country'),
        variable.name = 'history', value.name = 'probability') %>%
   mutate(history = factor(history, levels = history_order))

if(min(birth_year_range) >= 1988){
  history_order <- history_order[grepl('A', history_order) == F]
  inf_history_probs_pl <- inf_history_probs_pl %>%
  filter(grepl('A',as.character(history)) == F) %>%
    mutate(history = factor(history, levels = history_order))
  colors <- c('mediumpurple3','mediumpurple1','darkorange1','darkorange3','white')
}else{
  colors <- c('grey20','grey40','grey60','grey80','mediumpurple3','mediumpurple1','darkorange1',
              'darkorange2','white')
}

inf_history_probs_pl <- inf_history_probs_pl %>%
  filter(country == 'New Zealand') %>%
   ggplot(aes(x = birth_year, y = probability, fill = history)) +
   geom_col(color = 'black',size = 0.1) +
   #facet_grid(country ~ .) +
   #ylim(0,1+1e-1) +
   xlab('Year of birth') +
   ylab('Probability') +
   scale_fill_manual(name = 'Infection\nhistory\n',
                     labels = history_labels,
                     values = colors) +
   theme(legend.position = 'left',
         legend.text = element_text(size = 10),
         legend.spacing = unit(2,'cm'))  +
  geom_hline(yintercept = 0.5, linetype = 2) +
  #geom_vline(xintercept = c(1988.5,2000.5)) +
  scale_x_continuous(breaks = seq(birth_year_range[1],birth_year_range[2],5)) +
  guides(title.vjust = 1)

 save_plot(paste0(parent_directory, '/imprinting_probabilities_',observation_year,'.pdf'),
           inf_history_probs_pl,
            base_height = 5, base_width = 12)
 
 # Processed demographic data for numbers cited in main text
examples <- iprobs_by_byear %>% mutate(P_Y = P_AVY + P_AY0 + P_VY + P_YV + P_Y0,
                                      P_V = P_AVY + P_AV0 + P_VY + P_YV + P_V0,
                                      P_V_or_A = P_A0 + P_AY0 + P_AV0 + P_AVY +
                                        P_VY + P_V0,
                                      P_Y_or_A = P_A0 + P_AY0 + P_AV0 + P_AVY +
                                        P_VY + P_Y0,
                                      P_Y_first = P_Y0 + P_YV,
                                      P_V_first = P_V0 + P_VY,
                                      P_A_first = P_A0 + P_AVY + P_AV0 + P_AY0,
                                      P_V_not_Y = P_AV0 + P_V0,
                                      P_Y_not_V = P_AY0 + P_Y0,
                                      P_naive = P_0
                                      )  %>%
  select(country, birth_year, P_Y, P_V, P_V_or_A, P_Y_or_A, P_Y_first, 
         P_V_first, P_A_first, P_V_not_Y, P_Y_not_V, P_naive, rel_suscep_vic, rel_suscep_yam)


early_90s_cohorts <- left_join(examples %>%
            filter(birth_year >=1987, birth_year <=1993),
          demographic_data %>%
            filter(birth_year >= 1987, birth_year <= 1993) %>%
            group_by(country) %>%
            mutate(pop_fraction = pop_fraction / sum(pop_fraction)) %>%
            select(country, birth_year, pop_fraction) %>%
            ungroup(),
          by = c('country','birth_year')) %>%
  group_by(country) %>%
  summarise(P_Y = sum(P_Y*pop_fraction),
            P_V = sum(P_V*pop_fraction),
            P_Y_first = sum(P_Y_first*pop_fraction),
            P_V_first = sum(P_V_first*pop_fraction),
            P_V_not_Y = sum(P_V_not_Y*pop_fraction),
            P_Y_not_V = sum(P_Y_not_V*pop_fraction),
            P_naive = sum(P_naive *pop_fraction),
            rel_suscep_vic = sum(rel_suscep_vic*pop_fraction),
            rel_suscep_yam = sum(rel_suscep_yam*pop_fraction))

write.csv(early_90s_cohorts,
          paste0(parent_directory, '/early_90s_cohorts_susceptibility_in_',observation_year,'.csv'),
          row.names = F)

cohorts_until_mid80s <- left_join(examples %>%
            filter(birth_year <=1986),
          demographic_data %>%
            filter(birth_year <=1986) %>%
            group_by(country) %>%
            mutate(pop_fraction = pop_fraction / sum(pop_fraction)) %>%
            select(country, birth_year, pop_fraction) %>%
            ungroup(),
          by = c('country','birth_year')) %>%
  group_by(country) %>%
  summarise(P_Y = sum(P_Y*pop_fraction),
            P_V = sum(P_V*pop_fraction),
            P_A_first =  sum(P_A_first*pop_fraction),
            P_Y_first = sum(P_Y_first*pop_fraction),
            P_V_first = sum(P_V_first*pop_fraction),
            P_V_not_Y = sum(P_V_not_Y*pop_fraction),
            P_Y_not_V = sum(P_Y_not_V*pop_fraction),
            P_naive = sum(P_naive *pop_fraction),
            rel_suscep_vic = sum(rel_suscep_vic*pop_fraction),
            rel_suscep_yam = sum(rel_suscep_yam*pop_fraction))

write.csv(cohorts_until_mid80s,
          paste0(parent_directory, '/cohorts_until_mid80s_susceptibility_in_',observation_year,'.csv'),
          row.names = F)

 
