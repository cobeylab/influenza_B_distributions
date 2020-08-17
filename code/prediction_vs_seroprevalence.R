# Given specified model parameters, compares fraction of individuals 0-7 years old exposed to infuenza B
# in model predictions and seroprevalence data from Bodewes et al. (2011)

library(dplyr)
library(reshape2)
library(stringr)
library(ggplot2)
library(cowplot)
theme_set(theme_cowplot())
library(viridis)

args = commandArgs(trailingOnly = T)

intensity_scores_path = args[1] # intensity_scores_path = '../results/processed_data/intensity_scores.csv'
lineage_frequencies_path = args[2] # lineage_frequencies_path = '../results/processed_data/lineage_frequencies_gisaid-genbank_noVicin1990s.csv'
season_incidence_curves_path = args[3] # season_incidence_curves_path = '../results/processed_data/season_incidence_curves.csv'
pars_file <- as.character(args[4])

birth_year_range = c(2000,2007) #String with comma-separated values for the first and last birth year to consider
observation_year = 2007

parent_directory = dirname(pars_file)

seroprevalence_data <- read.csv('../results/processed_data/netherlands_seroprevalence.csv')


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

# Read parameters
pars <- as_tibble(read.csv(pars_file)) %>% mutate(par = as.character(par))
stopifnot(sum(is.na(pars$value)) == 0)


# Define objects for each parameter
for(i in 1:nrow(pars)){
  eval(parse(text = paste0(pars[i,1],'=',pars[i,2])))
}
chi_VY <- chi_Y * gamma_VY
chi_YV <- chi_V * gamma_YV
chi_AV <- chi_V * gamma_AV
chi_AY <- chi_Y * gamma_AY

# Tibble with combinations of country/birth year for the specified observation year
years_tibble <- as_tibble(expand.grid(c('New Zealand'),
                                      observation_year, birth_year_range[1]:birth_year_range[2])) %>%
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

probs_by_age <- iprobs_by_byear %>% 
  mutate(predicted_Vic_exp = P_V0 + P_VY + P_YV,
         predicted_Yam_exp = P_Y0 + P_VY + P_YV) %>% 
  mutate(predicted_B_exp = 1 - P_0, age = observation_year - birth_year) %>%
  select(country,age,predicted_B_exp, predicted_Vic_exp, predicted_Yam_exp) %>%
  arrange(country) %>%
  filter(age > 0)

# Add total number of children observed in Dutch data to model predictions (for binomial CIs)
probs_by_age <- left_join(probs_by_age,
          seroprevalence_data %>% select(age, n_total),
          by = 'age') %>%
  mutate(std_err_B = sqrt(predicted_B_exp*(1 - predicted_B_exp)/n_total),
         std_err_Vic = sqrt(predicted_Vic_exp*(1 - predicted_Vic_exp)/n_total),
         std_err_Yam = sqrt(predicted_Yam_exp*(1 - predicted_Yam_exp)/n_total)) %>%
  mutate(llim_B = predicted_B_exp - qnorm(0.975)*std_err_B,
         ulim_B =  predicted_B_exp + qnorm(0.975)*std_err_B,
         llim_Vic = predicted_Vic_exp - qnorm(0.975)*std_err_Vic,
         ulim_Vic =  predicted_Vic_exp + qnorm(0.975)*std_err_Vic,
         llim_Yam = predicted_Yam_exp - qnorm(0.975)*std_err_Yam,
         ulim_Yam =  predicted_Yam_exp + qnorm(0.975)*std_err_Yam)

# Set lower bounds to zero if they are negative
probs_by_age$llim_B[probs_by_age$llim_B < 0] <- 0
probs_by_age$llim_Vic[probs_by_age$llim_Vic < 0] <- 0
probs_by_age$llim_Yam[probs_by_age$llim_Yam < 0] <- 0

# Correlations
correlations <- left_join(probs_by_age, seroprevalence_data, by = 'age') %>%
  group_by(country) %>%
  summarise(cor = cor.test(predicted_B_exp, seroprevalence)[[4]],
            llim = cor.test(predicted_B_exp, seroprevalence)[9]$conf.int[1],
            ulim = cor.test(predicted_B_exp, seroprevalence)[9]$conf.int[2])


# Merge predictions with seroprevalence data for plotting
probs_by_age <- full_join(probs_by_age,
                          seroprevalence_data %>% mutate(country = 'Netherlands') %>%
                            rename(predicted_B_exp = seroprevalence,
                                   predicted_Vic_exp = Vic_positive_fraction,
                                   predicted_Yam_exp = Yam_positive_fraction)) %>%
  mutate(country = factor(country, levels = c('New Zealand','Netherlands')))





nudge_neth = position_nudge(x = 0.1)

prevalence_by_age <- ggplot(probs_by_age, aes(x = age, y = predicted_B_exp,
                                              ymin = llim_B, ymax = ulim_B)) +
  geom_linerange(data = probs_by_age %>% filter(country == 'New Zealand'),
                 color = 'orange2') +
  geom_point(data = probs_by_age %>% filter(country == 'New Zealand'),
             shape = 21, fill = 'orange2', color = 'black', size = 3) +
  geom_linerange(data = probs_by_age %>% filter(country == 'Netherlands'), position = nudge_neth,
                 color = 'purple') +
  geom_point(data = probs_by_age %>% filter(country == 'Netherlands'), position = nudge_neth,
             shape = 21, fill = 'purple3', color = 'black', size = 3) +
  # adding empty geometry just to add color guide
  geom_point(data = tibble(age = 1, predicted_B_exp = -2,llim_B = -2,ulim_B = -2,
                           country = factor(c('New Zealand','Netherlands'),
                                            levels = c('New Zealand','Netherlands'))),
             aes(fill = country), shape = 21, size = 3) +
  scale_y_continuous(limits = c(-0.02,0.85), breaks = seq(0,0.8,0.1),
                     name = 'Fraction of people previously infected') +
  scale_fill_manual(values = c('orange2','purple3'),
                     labels = c('Prediction for New Zealand',
                                'Seroprevalance in the Netherlands'),
                     name = '') +
  theme(legend.position = c(0.02,0.95)) +
  ggtitle("Any influenza B") +
  scale_x_continuous(breaks = 1:7, name  = '')

# Plots by lineage
prevalence_by_age_Vic <- ggplot(probs_by_age, aes(x = age, y = predicted_Vic_exp,
                         ymin = llim_Vic, ymax = ulim_Vic)) +
  geom_linerange(data = probs_by_age %>% filter(country == 'New Zealand'),
                 color = 'orange2') +
  geom_point(data = probs_by_age %>% filter(country == 'New Zealand'),
             shape = 21, fill = 'orange2', color = 'black', size = 3) +
  geom_linerange(data = probs_by_age %>% filter(country == 'Netherlands'), position = nudge_neth,
                 color = 'purple') +
  geom_point(data = probs_by_age %>% filter(country == 'Netherlands'), position = nudge_neth,
             shape = 21, fill = 'purple3', color = 'black', size = 3) +
  scale_y_continuous(limits = c(-0.02,0.85), breaks = seq(0,0.8,0.1), name = '')+
  theme(legend.position = 'None', axis.text.y = element_blank()) +
  ggtitle("B/Victoria") +
  scale_x_continuous(breaks = 1:7, name  = 'Age (years) in 2007')

prevalence_by_age_Yam <- ggplot(probs_by_age, aes(x = age, y = predicted_Yam_exp,
                                                  ymin = llim_Yam, ymax = ulim_Yam)) +
  geom_linerange(data = probs_by_age %>% filter(country == 'New Zealand'),
                 color = 'orange2') +
  geom_point(data = probs_by_age %>% filter(country == 'New Zealand'),
             shape = 21, fill = 'orange2', color = 'black', size = 3) +
  geom_linerange(data = probs_by_age %>% filter(country == 'Netherlands'), position = nudge_neth,
                 color = 'purple') +
  geom_point(data = probs_by_age %>% filter(country == 'Netherlands'), position = nudge_neth,
             shape = 21, fill = 'purple3', color = 'black', size = 3) +
  theme(axis.text.y = element_blank()) +
  ggtitle("B/Yamagata") +
  scale_x_continuous(breaks = 1:7, name = '') +
  scale_y_continuous(limits = c(-0.02,0.85), breaks = seq(0,0.8,0.1), name = '')

save_plot(paste0(parent_directory, '/prediction_vs_seroprevalence.pdf'),
          plot_grid(prevalence_by_age, prevalence_by_age_Vic, prevalence_by_age_Yam,
                    nrow = 1),
          base_height = 6, base_width = 18)

