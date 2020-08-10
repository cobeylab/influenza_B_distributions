#!/usr/bin/env Rscript
# Plots likelihood profile
library(dplyr)
library(cowplot)
library(viridis)
library(stringr)
library(ggridges)

pairwise_divergence <- as_tibble(read.csv('../results/sequence_divergence/pairwise_divergence.csv',
                                header = T))
pairwise_divergence$segment <- as.character(pairwise_divergence$segment)
pairwise_divergence$segment[is.na(pairwise_divergence$segment)] = 'NA'
pairwise_divergence$segment <- factor(pairwise_divergence$segment, 
                                      levels = c('head','stalk','NA', '120_loop','150_loop',
                                                 '160_loop','190_helix','head_plus_stalk'))

# Exclude isolates that had both H1N1-like and H3N2-like HA sequences on GISAID.
pairwise_divergence <- filter(pairwise_divergence,
                              grepl('A/Canterbury/58/2000', pair) == F,
                              grepl('A/Canterbury/87/2000', pair) == F,
                              grepl('A/Canterbury/55/2000', pair) == F
                              )


divergence_from_ref <- as_tibble(read.csv('../results/sequence_divergence/divergence_from_reference.csv',
                                          header = T)) 
divergence_from_ref$segment <- as.character(divergence_from_ref$segment)
divergence_from_ref$segment[is.na(divergence_from_ref$segment)] = 'NA'
divergence_from_ref$segment <- factor(divergence_from_ref$segment, 
                                      levels = c('head','stalk','NA', '120_loop','150_loop',
                                                 '160_loop','190_helix','head_plus_stalk'))


# Glycosylation
glycosylation <- as_tibble(read.csv('../results/sequence_divergence/glycosylation_over_time.csv', header = T)) %>%
  mutate(segment = as.character(segment))
glycosylation$segment[is.na(glycosylation$segment)] <- 'NA'


# Add total number of seqs. per year to look for PNGS present in > x% of seqs.
n_seqs_per_year <- divergence_from_ref %>% filter(segment %in% c('head','NA')) %>%
  mutate(segment = as.character(segment))
n_seqs_per_year$segment[n_seqs_per_year$segment == 'head'] = 'HA'
n_seqs_per_year$segment[n_seqs_per_year$segment == 'head'] = 'HA'
n_seqs_per_year <- n_seqs_per_year %>%
  ungroup() %>%
  group_by(segment,lineage,year) %>%
  summarise(total_n_seqs = length(unique(strain)))

glycosylation <- left_join(glycosylation %>% group_by(segment, lineage, year, glycosylation_site) %>%
    summarise(seqs_with_PNGS = n()) %>% ungroup(),
    n_seqs_per_year,
  by = c('year','segment','lineage')) %>%
  mutate(fraction_seqs_with_PNGS = seqs_with_PNGS/total_n_seqs)

glycosylation$segment <- factor(glycosylation$segment, 
                                levels = c('HA','NA'))
glycosylation <- glycosylation %>%
  mutate(glycosylation_site = factor(as.character(glycosylation_site),
                                     levels = as.character(sort(unique(glycosylation$glycosylation_site))))
         )


segment_labels = c('HA head', 'HA stalk', 'NA')
names(segment_labels) = c('head','stalk','NA')

# Function for downsampling pairwise divergence points for plotting
downsample_pwdiv <- function(pairwise_divergence, subsample_size = 500){
  downs_tibble <- c()
  for(year in unique(pairwise_divergence$year)){
    for(segment in unique(pairwise_divergence$segment)){
      for(pair_type in unique(pairwise_divergence$pair_type)){
        sub_tibble <- pairwise_divergence %>% filter(
          segment == !!segment,
          year == !! year,
          pair_type == !!pair_type
        )
        if(nrow(sub_tibble) > subsample_size){
          sub_tibble <- sub_tibble[sample(1:nrow(sub_tibble),subsample_size,
                                          replace = F),]
        }
        downs_tibble <- bind_rows(downs_tibble, sub_tibble)
      } 
    }
  }
  return(downs_tibble)
}
pairwise_plotting_tibble <- downsample_pwdiv(pairwise_divergence, subsample_size = 500)


pairwise_annual_averages_tibble <- pairwise_divergence %>% 
  group_by(year, segment, pair_type) %>%
  summarise(pairwise_divergence = mean(pairwise_divergence)) %>%
  ungroup()

div_from_ref_averages_tibble <- divergence_from_ref %>%
  group_by(year, segment, lineage) %>%
  summarise(divergence_from_reference = mean(divergence_from_reference,na.rm=T)) %>%
  ungroup() %>%
  filter((lineage == 'Victoria' & year >= 1987) | (lineage == 'Yamagata' & year >=1988))


# ======================= Divergence in B: within and between lineages (head, stem and NA) ================
seq_pairwise_divergence_pl <- ggplot(pairwise_plotting_tibble %>%
                                       filter(segment %in% c('head','stalk','NA')) %>%
                                       filter(grepl('H',pair_type) == F) %>%
                                       mutate(pair_type = factor(pair_type,
                                                                 levels = c('Yamagata_vs_Victoria',
                                                                            'Victoria','Yamagata'))),
  aes(x = year, y = pairwise_divergence, fill = pair_type)) +
  geom_point(alpha = 0.1, size = 2, shape = 21) +
  facet_grid(segment~.) +
  geom_line(data = pairwise_annual_averages_tibble %>%
              filter(segment %in% c('head','stalk','NA')) %>%
              filter(grepl('H',pair_type) == F) %>%
              mutate(pair_type = factor(pair_type,
                                        levels = c('Yamagata_vs_Victoria',
                                                   'Victoria','Yamagata'))),
            size = 1.5, aes(color = pair_type)) +
  theme(legend.position = 'top') +
  ylab('Pairwise divergence (fraction of amino acid sites)') +
  xlab('Year') +
  scale_color_manual(name = '', labels = c('Between lineages','B/Victoria','B/Yamagata'),
                     values = c('grey20','darkorange1', 'mediumpurple1')) +
  scale_fill_manual(name = '', labels = c('Between lineages','B/Victoria','B/Yamagata'),
                    values = c('grey20','darkorange1', 'mediumpurple1')) +
  scale_y_continuous(limits = c(0,0.20), breaks = seq(0,0.20,0.05)) +
  scale_x_continuous(limits = c(1985,2020), breaks = seq(1985,2020,5)) +
  theme(plot.margin = margin(1,0.2,1,1,'cm'),
        legend.position = c(0.05,0.98),
        legend.direction = 'vertical',
        legend.text = element_text(size = 10),
        strip.background = element_blank(),
        strip.text = element_blank()) +
  background_grid(minor = 'none')

divergence_from_ref_pl <- ggplot(divergence_from_ref %>%
                                   filter(segment %in% c('head','stalk','NA')),
                                 aes(x = year, y = divergence_from_reference, fill = lineage)) +
  geom_point(alpha = 0.5, size = 2, shape = 21) +
  facet_grid(segment~.,
             labeller = labeller(segment = segment_labels)) +
  geom_line(data = div_from_ref_averages_tibble %>%
              filter(segment %in% c('head','stalk','NA')),
            size = 1.5, aes(color = lineage)) +
  theme(legend.position = 'top') +
  #ylim(0,0.125)+
  ylab('Divergence from reference strain (fraction of amino acid sites)') +
  xlab('Year') +
  scale_color_manual(name = '', labels = c('B/Victoria','B/Yamagata'),
                     values = c('darkorange1', 'mediumpurple1')) +
  scale_fill_manual(name = '', labels = c('Between lineages','B/Victoria','B/Yamagata'),
                    values = c('darkorange1', 'mediumpurple1')) +
  scale_y_continuous(limits = c(0,0.20), breaks = seq(0,0.20,0.05)) +
  scale_x_continuous(limits = c(1985,2020), breaks = seq(1985,2020,5)) +
  theme(plot.margin = margin(1,1,1,0.1,'cm'),
        legend.position = 'none') +
  background_grid(minor = 'none')
  

  
save_plot('../figures/sequence_divergence_B.pdf',
          plot_grid(seq_pairwise_divergence_pl,
                    divergence_from_ref_pl, nrow = 1),
          base_width = 10, base_height = 10)

# ====================== Divergence from reference in canonical sites ========================

canonical_pairwise_divergence_pl <- pairwise_plotting_tibble %>%
  filter(segment %in% c('120_loop','150_loop','160_loop','190_helix')) %>%
    mutate(pair_type = factor(pair_type,
                            levels = c('Yamagata_vs_Victoria',
                                       'Victoria','Yamagata'))) %>%
  ggplot(aes(x = year, y = pairwise_divergence, fill = pair_type)) +
  geom_point(alpha = 0.5, size = 2, shape = 21) +
  facet_grid(segment~.) +
  geom_line(data = pairwise_annual_averages_tibble %>%
              filter(segment %in% c('120_loop','150_loop','160_loop','190_helix')) %>%
              filter(grepl('H',pair_type) == F) %>%
              mutate(pair_type = factor(pair_type,
                                        levels = c('Yamagata_vs_Victoria',
                                                   'Victoria','Yamagata'))),
            size = 1.5, aes(color = pair_type)) +
  theme(legend.position = 'top') +
  ylab('Pairwise divergence (number of amino acid differences)') +
  xlab('Year') +
  scale_color_manual(name = '', labels = c('Between lineages','B/Victoria','B/Yamagata'),
                     values = c('grey20','darkorange1', 'mediumpurple1')) +
  scale_fill_manual(name = '', labels = c('Between lineages','B/Victoria','B/Yamagata'),
                    values = c('grey20','darkorange1', 'mediumpurple1')) +
  scale_x_continuous(limits = c(1985,2020), breaks = seq(1985,2020,5)) +
  scale_y_continuous(limits = c(0,6), breaks = 0:6) +
  theme(plot.margin = margin(1,0.2,1,1,'cm'),
        legend.position = c(0.05,0.98),
        legend.direction = 'vertical',
        legend.text = element_text(size = 10),
        strip.background = element_blank(),
        strip.text = element_blank()) +
  background_grid(minor = 'none')
   
divergence_from_ref_canonical_pl <- ggplot(divergence_from_ref %>%
                                             filter(segment %in% c('120_loop','150_loop',
                                                                   '160_loop','190_helix')),
                                 aes(x = year, y = divergence_from_reference, fill = lineage)) +
  geom_point(alpha = 0.5, size = 2, shape = 21) +
  facet_grid(segment~.) +
  geom_line(data = div_from_ref_averages_tibble %>%
              filter(segment %in% c('120_loop','150_loop',
                                    '160_loop','190_helix')),
            size = 1.5, aes(color = lineage)) +
  theme(legend.position = 'top') +
  #ylim(0,0.125)+
  ylab('Divergence from reference strain (number of amino acid differences)') +
  xlab('Year') +
  scale_color_manual(name = '', labels = c('B/Victoria','B/Yamagata'),
                     values = c('darkorange1', 'mediumpurple1')) +
  scale_fill_manual(name = '', labels = c('Between lineages','B/Victoria','B/Yamagata'),
                    values = c('darkorange1', 'mediumpurple1')) +
  scale_x_continuous(limits = c(1985,2020), breaks = seq(1985,2020,5)) +
  scale_y_continuous(limits = c(0,6), breaks = 0:6) +
  theme(plot.margin = margin(1,1,1,0.1,'cm'),
        legend.position = 'none') +
  background_grid(minor = 'none')

save_plot('../figures/sequence_divergence_B_canonical_sites.pdf',
          plot_grid(canonical_pairwise_divergence_pl,
                    divergence_from_ref_canonical_pl, nrow = 1),
          base_width = 10, base_height = 10)


# ======================= Divergence in B stalk compared with A ================

between_divergence_A_and_B <- ggplot(pairwise_plotting_tibble %>%
         filter(segment %in% c('head','stalk','NA')) %>%
         filter(grepl('_vs_', pair_type)) %>%
         mutate(pair_type = factor(pair_type,
                                   levels = c('Yamagata_vs_Victoria',
                                              'H1N1_vs_H3N2'))),
      aes(x = year, y = pairwise_divergence, fill = pair_type)) +
  geom_point(alpha = 0.1, size = 2, shape = 21) +
  facet_grid(.~segment,
             labeller = labeller(segment = segment_labels)) +
  # geom_line(data = pairwise_plotting_tibble %>%
  #             filter(grepl('_vs_', pair_type), segment != 'NA') %>%
  #             mutate(pair_type = factor(pair_type,
  #                                       levels = c('Yamagata_vs_Victoria',
  #                                                  'H1N1_vs_H3N2'))),
  #           size = 1.5, aes(color = pair_type)) +
  # theme(legend.position = 'top') +
  ylab('Pairwise divergence\n(fraction of amino acid sites)') +
  xlab('Year') +
  scale_color_manual(name = '', labels = c('Between B/Victoria and B/Yamagata',
                                           'Between H1N1 and H3N2'),
                     values = c('grey20','dodgerblue')) +
  scale_fill_manual(name = '', labels = c('Between B/Victoria and B/Yamagata',
                                          'Between H1N1 and H3N2'),
                    values = c('grey20','dodgerblue')) +
  theme(plot.margin = margin(1,1,1,0.1,'cm'),
        legend.position = 'top') +
  background_grid(minor = 'none') +
  guides(fill = guide_legend(override.aes = list(alpha=1)))
  
save_plot('../figures/sequence_divergence_A_and_B.pdf',
          between_divergence_A_and_B,
          base_width = 12, base_height = 5)

# ggplot(pairwise_plotting_tibble %>%
#          filter(segment != 'NA', grepl('_vs_', pair_type) == F),
#        aes(x = year, y = pairwise_divergence, fill = pair_type)) +
#   #geom_point(alpha = 0.1, size = 2, shape = 21) +
#   facet_grid(segment~.) +
#   geom_line(data = pairwise_annual_averages_tibble %>%
#               filter(segment != 'NA', grepl('_vs_', pair_type) == F),
#             size = 1.5, aes(color = pair_type)) +
#   guides(fill = guide_legend(override.aes = list(alpha=1))) +
#   xlim(c(1990,2020)) +
#   ylim(c(0,0.2))

# ======================= Glycosylation over time ================
minfreq <- 0.5 # Minimum frequency of PNGS site for plotting

HA_glycosylation <- glycosylation %>% filter(segment == 'HA',
                                             fraction_seqs_with_PNGS >= minfreq, year >=1983) %>%
  ggplot(aes(x = year, y = glycosylation_site, size = lineage)) +
  geom_point(aes(fill = lineage), shape = 21,
             alpha = 1) +
  scale_size_manual(values = c(5,2.5)) +
  background_grid(major = "xy", minor = "none", colour.major = 'grey80') +
  scale_x_continuous(breaks = seq(1980,2020,5)) +
  xlab('Year') +
  ylab(paste0('Potential N-linked glygosylation site\nwith ', minfreq*100, '% frequency or higher')) +
  ggtitle('HA') +
  theme(legend.position = 'top') +
  scale_fill_manual(values = c('darkorange1', 'mediumpurple1'))

NA_glycosylation <- glycosylation %>% filter(segment == 'NA',
                                             fraction_seqs_with_PNGS >= minfreq, year >=1983) %>%
  ggplot(aes(x = year, y = glycosylation_site, size = lineage)) +
  geom_point(aes(fill = lineage), shape = 21,
             alpha = 1) +
  scale_size_manual(values = c(5,2.5)) +
  background_grid(major = "xy", minor = "none", colour.major = 'grey80') +
  scale_x_continuous(breaks = seq(1980,2020,5)) +
  xlab('Year') +
  ylab('') +
  ggtitle('NA') +
  theme(legend.position = 'top') +
  scale_fill_manual(values = c('darkorange1', 'mediumpurple1'))

save_plot('../figures/glycosylation_over_time.pdf',
          plot_grid(HA_glycosylation, NA_glycosylation),
          nrow = 1, base_width = 16, base_height = 8)


specific_sites_pl <- glycosylation %>%
  filter(segment == 'HA', year >=1990, total_n_seqs >=0,
         glycosylation_site %in% c('214', '185')) %>%
  ggplot(aes(x = year, y = fraction_seqs_with_PNGS,
             color = lineage)) +
  geom_vline(xintercept = c(2005, 2008, 2011), linetype = 2) +
  geom_point() +
  geom_text(aes(label = paste(seqs_with_PNGS,total_n_seqs, sep = '/')),
            position = position_nudge(y = 0.02), size = 2.5) +
  geom_line() + 
  scale_x_continuous(breaks = 1990:2019) +
  theme(legend.position = c(0.05,0.9),
        axis.text.x = element_text(size = 9.5)) +
  facet_grid(glycosylation_site~.) +
  ylab('Fraction of sequences with \npotential N-linked glygosylation site') +
  scale_color_manual(values = c('darkorange1', 'mediumpurple1'),
                    name = '')
  

save_plot('../figures/glycosylation_sites_214_185.pdf',
          plot_grid(specific_sites_pl),
          nrow = 1, base_width = 14, base_height = 8)

