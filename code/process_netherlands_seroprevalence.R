library(dplyr)


# Seroprevalence data 
seroprevalence_data <- as_tibble(read.csv('../data/seroprevalence_data/Bodewes2011_seroprevalence.csv', header = T))
seroprevalence_data <- seroprevalence_data %>%
  mutate(std_err = sqrt(seroprevalence*(1 - seroprevalence)/n_total)) %>%
  mutate(llim_B = seroprevalence - qnorm(0.975)*std_err,
         ulim_B = seroprevalence + qnorm(0.975)*std_err)

# Data by lineage
seroprevalence_by_lineage <- as_tibble(read.csv('../data/seroprevalence_data/Bodewes_data_by_lineage.csv', header = T))

seroprevalence_by_lineage <- seroprevalence_by_lineage %>%
  group_by(age) %>%
  summarise(Vic_positive = sum(Vic_positive),
            Yam_positive = sum(Yam_positive), n_total = n()) %>%
  ungroup() %>%
  mutate(Vic_positive_fraction = Vic_positive / n_total,
         Yam_positive_fraction = Yam_positive / n_total) %>%
  mutate(std_err_Vic = sqrt(Vic_positive_fraction*(1 - Vic_positive_fraction)/n_total),
         std_err_Yam = sqrt(Yam_positive_fraction*(1 - Yam_positive_fraction)/n_total)) %>%
  mutate(llim_Vic = Vic_positive_fraction - qnorm(0.975)*std_err_Vic,
         ulim_Vic =  Vic_positive_fraction + qnorm(0.975)*std_err_Vic,
         llim_Yam = Yam_positive_fraction - qnorm(0.975)*std_err_Yam,
         ulim_Yam =  Yam_positive_fraction + qnorm(0.975)*std_err_Yam)

seroprevalence_data <- left_join(seroprevalence_data, seroprevalence_by_lineage %>% select(age, matches('fraction'), matches('lim')),
                                 by = 'age')

write.csv(seroprevalence_data, '../results/processed_data/netherlands_seroprevalence.csv', row.names = F)

