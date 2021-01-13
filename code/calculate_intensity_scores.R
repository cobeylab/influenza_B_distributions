# Processes virological surveillance data to calculate an influenza B intensity score for Australia and NZ
library(cowplot)
library(ggplot2)
theme_set(theme_cowplot())
library(dplyr)
library(tidyr)
library(lubridate)
library(FlexParamCurve)

surveillance_data_directory = "../data/surveillance_intensity/"
output_directory = '../results/processed_data/'

# Year before which B intensity is set to 0
B_origin_year = 1940

# Function for getting number of ISO weeks in a year
n_iso_weeks <- function(year){
  last_days <- paste(year, '12', seq(20,31), sep = '-')
  n_iso_weeks <- max(lubridate::isoweek(last_days))
  return(n_iso_weeks)
}

# Initial parameters for the beta CDF fit of cumulative incidence
initial_beta_parameters = list(shape1 = 2, shape2 = 2)

# Mininum number of cases in a season (below which cumulative incidence is calculated from beta fit)
min_total_cases = 50

# Function for filling NAs in the specified column of a tibble by either averaging or sampling from non NA values
fill_nas <- function(df, variable, method){
  variable <- enquo(variable)
  variable_name <- paste(quo_name(variable))
  
  nonNA_values <- df %>% filter(is.na(!!variable) == F) %>% pull(!! variable)
  
  if(method == 'average'){
    mean_value <- mean(nonNA_values)
    df <- df %>%
      mutate(!!variable_name := ifelse(is.na(!!variable),
                                       mean_value,
                                       !!variable))
  }else{
    stopifnot(method == 'sample')
    df <- df %>% rowwise() %>%
      mutate(!!variable_name := ifelse(is.na(!!variable),
                                       sample(nonNA_values,size = 1),
                                       !!variable)) %>%
      ungroup()
  }
  return(df)
}

interpolate_cumulative_incidence <- function(season_start_year, season_week, cumulative_case_data,
                                             min_total_cases){
  # Check if season_start_year for country is represented in cumulative_case_data
  season_represented = as.logical(nrow(filter(cumulative_case_data, season_start_year == !!season_start_year)))
  if(!season_represented){
    return(NA) # If season_start_year not represented in data, return NA
  }else{
    # Check if there is a value for exactly season_week of season_start_year in country
    exact_value <- filter(cumulative_case_data, season_start_year == !!season_start_year,
                          season_week == !!season_week) %>%
      pull(cumulative_B_fraction)
    
    if(length(exact_value) == 1){
      return(exact_value)
    }else{
      stopifnot(length(exact_value) == 0) # Stop if there's > 1 entry for this week/season/country
      
      # Get cumulative incidence for latest week before and earliest week after season_week
      nearest_estimates <- filter(cumulative_case_data, season_start_year == !!season_start_year,
                                  season_week != !!season_week) %>% 
        mutate(period = ifelse(season_week < !!season_week, 'before','after')) %>%
        group_by(period) %>%
        filter(season_week == max(season_week) | season_week == min(season_week)) %>%
        ungroup()
      
      lower <- filter(nearest_estimates, period == 'before') %>% 
        filter(season_week == max(season_week)) %>%
        select(season_week, cumulative_B_fraction)
      # If no lower week is represented in the data, set lower week to 0 and lower value to 0 
      lower_week <- ifelse(nrow(lower) == 1, lower$season_week, 0)
      lower_value <- ifelse(nrow(lower) == 1, lower$cumulative_B_fraction, 0)
      
      upper <- filter(nearest_estimates, period == 'after') %>%
        filter(season_week == min(season_week)) %>%
        select(season_week, cumulative_B_fraction)
      # If no upper week is represented in the data, set week to max. week in the season and upper value to 1
      upper_week <- ifelse(nrow(upper) == 1, upper$season_week, n_iso_weeks(season_start_year))
      upper_value <- ifelse(nrow(upper) == 1, upper$cumulative_B_fraction, 1)
      
      derivative <- (upper_value - lower_value)/(upper_week - lower_week)
      estimate <- lower_value +  derivative * (season_week - lower_week)
      return(estimate)
    }
  }
}

# Function for calculating season-level intensity scores
get_intensity_scores <- function(surveillance_data_directory, case_data_directory, na_fill_method){
  # =============================================================================================
  # =========== Read and process surveillance data to calculate between season intensity ========
  # =============================================================================================
  
  # -- New Zealand ------------------------------------------------------------------------------
  
  # Intensity across seasons
  between_season_intensity_nz <- read.csv(paste0(surveillance_data_directory,'ESR_surveillance_nz.csv'), header = T)
  between_season_intensity_nz <- as_tibble(between_season_intensity_nz) %>%
    select(year, A_fraction_in_flu, B_fraction_in_flu, A_fraction_in_ILI, B_fraction_in_ILI, flu_fraction_in_ILI)
  
  # Add entries for early years
  between_season_intensity_nz <- full_join(between_season_intensity_nz,
                                           tibble(year = seq(1900,1989)), by = 'year') %>%
    arrange(year)
    
  between_season_intensity_nz <- between_season_intensity_nz %>% 
    mutate(B_fraction_in_ILI_source = ifelse(is.na(B_fraction_in_flu) | is.na(flu_fraction_in_ILI),
                                             paste0(na_fill_method,'d from empirical distribution'),
                                             'observed')) %>%
    # Fill years with missing B % in flu by either sampling or taking the average (na_fill_method 'sample' or 'average')
    fill_nas(variable = A_fraction_in_flu, method = na_fill_method) %>%
    fill_nas(variable = B_fraction_in_flu, method = na_fill_method) %>%
    # Fill years with missing flu % in ILI by either sampling or taking the average (na_fill_method 'sample' or 'average')
    fill_nas(variable = flu_fraction_in_ILI, method = na_fill_method) %>%
    # Calculate B fraction in ILI as A % in Flu * Flu % in ILI
    mutate(A_fraction_in_ILI = A_fraction_in_flu * flu_fraction_in_ILI,
           B_fraction_in_ILI = B_fraction_in_flu * flu_fraction_in_ILI) %>%
    mutate(country = 'New Zealand') %>%
    rename(season_start_year = year) %>%
    select(country, season_start_year, A_fraction_in_ILI, B_fraction_in_ILI, B_fraction_in_ILI_source)
  
  # Total ILI (estimated fraction of NZ population with ILI)
  total_ILI_nz <- read.csv(paste0(surveillance_data_directory,'total_ILI_nz.csv'), header = T)
  total_ILI_nz <- as_tibble(total_ILI_nz) %>% mutate(country = 'New Zealand') %>%
    rename(season_start_year = season) %>% select(country, season_start_year, ILI_per_100K)
  
  # Merge total ILI and fraction B in ILI
  between_season_intensity_nz <- full_join(between_season_intensity_nz, total_ILI_nz)

  between_season_intensity_nz <- between_season_intensity_nz %>%
    # Fill years with ILI per 100K by either sampling or taking the average (na_fill_method 'sample' or 'average')
    mutate(ILI_per_100K_source = ifelse(is.na(ILI_per_100K),
                                        paste0(na_fill_method,'d from empirical distribution'), 
                                        'observed')) %>%
    fill_nas(variable = ILI_per_100K, method = na_fill_method) %>%
    # Define B intensity as % B in ILI times ILI cases per 100 thousand people
    mutate(B_intensity = B_fraction_in_ILI * ILI_per_100K) %>%
    # Calculate normalized intensity scores as B intensity / mean B intensity
    mutate(intensity_score = B_intensity / mean(B_intensity, na.rm=T))
  
  # -- Australia ---------------------------------------------------------------------------------
  
  # Read Australian surveillance data from the WHO
  between_season_intensity_aus_WHO <- read.csv(paste0(surveillance_data_directory,'WHO_surveillance_australia.csv'))
  between_season_intensity_aus_WHO <- as_tibble(between_season_intensity_aus_WHO) %>%
    # Remove data from before 2003
    filter(Year >= 2003,
           # Remove weeks where the number of specimens tested cannot be determined
           # (i.e. for which neither the # of specimens tested nor the # of flu- specimens are given)
           is.na(SPEC_PROCESSED_NB) == F | is.na(ALL_INF2) == F ) %>%
    # Where # of specimens tested is not reported, calculate it from + and - specimens:
    mutate(SPEC_PROCESSED_NB = ifelse(is.na(SPEC_PROCESSED_NB), ALL_INF + ALL_INF2, SPEC_PROCESSED_NB)
    ) %>%
    # Rename variables
    rename(total_isolates_tested = SPEC_PROCESSED_NB, A_positive = INF_A, B_positive = INF_B,
           year = Year, week = Week, start_date = SDATE,
           end_date = EDATE) %>%
    mutate(start_date = as.Date(start_date), end_date = as.Date(end_date), country = 'Australia') %>%
    select(country,year,week,start_date,end_date,total_isolates_tested, A_positive,B_positive) %>%
    # Aggregate by season
    group_by(year) %>%
    summarise(total_isolates_tested = sum(total_isolates_tested),
              A_positive = sum(A_positive), B_positive = sum(B_positive)) %>%
    # Calculate A and B fractions among flu isolates and among ILI specimens
    mutate(flu_positive = A_positive + B_positive,
           A_fraction_in_flu = A_positive / flu_positive,
           B_fraction_in_flu = B_positive / flu_positive,
           A_fraction_in_ILI = A_positive / total_isolates_tested,
           B_fraction_in_ILI = B_positive / total_isolates_tested) %>%
    ungroup()
  
  # Read Australian surveillance data from Australian Department of Health
  between_season_intensity_aus_DH <- read.csv(paste0(surveillance_data_directory,'DepHealth_surveillance_aus.csv'))
  between_season_intensity_aus_DH <- as_tibble(between_season_intensity_aus_DH) %>%
    # Select A and B fractions among flu isolates
    select(year, A_fraction_in_flu, B_fraction_in_flu) %>%
    # Select only data from before 2003 (when WHO data are not available or inconsistent)
    filter(year <= 2002)
  
  # Merge surveillance data for Australia from the WHO and the Department of Health
  between_season_intensity_aus <- full_join(between_season_intensity_aus_WHO, between_season_intensity_aus_DH) %>%
    arrange(year)
  
  # Add years 1982-1993
  between_season_intensity_aus <- full_join(between_season_intensity_aus,
                                            tibble(year = seq(1900,1993,1)), 
                                                   by = 'year') %>%
    arrange(year)
  
  # Get flu fraction in ILI from A and B fractions (for years that have them)
  between_season_intensity_aus <- between_season_intensity_aus %>%
     mutate(flu_fraction_in_ILI = A_fraction_in_ILI + B_fraction_in_ILI)
  
  between_season_intensity_aus <- between_season_intensity_aus %>%
      # Fill years with missing % flu in ILI by either sampling or taking the average (na_fill_method 'sample' or 'average')
      mutate(B_fraction_in_ILI_source = ifelse(is.na(flu_fraction_in_ILI),
                                               paste0(na_fill_method,'d from empirical distribution'),
                                               'observed')) %>%
      fill_nas(variable = flu_fraction_in_ILI, method = na_fill_method) %>%
      # For years w/o B fraction in ILI but with B fraction in Flu, get B % in ILI by multiplying by the averaged/sampled % flu in ILI
      mutate(A_fraction_in_ILI = A_fraction_in_flu * flu_fraction_in_ILI,
             B_fraction_in_ILI = B_fraction_in_flu * flu_fraction_in_ILI) %>%
      # For years w/o B fraction in flu (which also don't have B or flu fraction in ILI), sample/average from available values
      fill_nas(variable = B_fraction_in_ILI, method = na_fill_method) %>%
      fill_nas(variable = A_fraction_in_ILI, method = na_fill_method) %>%
      mutate(country = 'Australia') %>%
      rename(season_start_year = year) %>%
      select(country, season_start_year, A_fraction_in_ILI, B_fraction_in_ILI, B_fraction_in_ILI_source)  %>%
      # Remove 2017 season (no case data from that period)
      filter(season_start_year != 2017)
  
  # Read ILI incidence
  ili_incidence_aus <- tibble(season = 1900:1993)
  ili_incidence_aus <- bind_rows(ili_incidence_aus,
                                 as_tibble(read.csv(paste0(surveillance_data_directory,'total_ILI_aus.csv'), header = T)) %>%
                                  # Normalize across different value types (different ILI definitions, data aggregation across surveillance systems)
                                 group_by(value_type) %>%
                                 mutate(group_average = mean(peak_ILI_per_1000_visits),
                                        peak_ILI_per_1000_visits = peak_ILI_per_1000_visits/ group_average) %>%
                                 ungroup()
                      ) %>%
    # Fill years with missing peak ILI by either sampling or taking the average (na_fill_method 'sample' or 'average')
    fill_nas(variable = peak_ILI_per_1000_visits, method = na_fill_method) %>%
    rename(season_start_year = season)
  
  
  between_season_intensity_aus <- left_join(between_season_intensity_aus, ili_incidence_aus %>% 
                                              select(season_start_year, peak_ILI_per_1000_visits),
            by = 'season_start_year') %>%
    # Define B intensity as % B in ILI times ILI cases per 100 thousand people
    mutate(B_intensity = B_fraction_in_ILI * peak_ILI_per_1000_visits) %>%
    # Calculate normalized intensity scores as B intensity / mean B intensity
    mutate(intensity_score = B_intensity / mean(B_intensity, na.rm=T))
  
  between_season_intensity_ausnz <- bind_rows(between_season_intensity_nz %>% 
                                                select(country, season_start_year, B_intensity, intensity_score),
                                              between_season_intensity_aus %>%
                                                select(country, season_start_year, B_intensity, intensity_score))
  
  return(between_season_intensity_ausnz)
}

# Function for getting cumulative fraction of B cases experienced as a function of time within a season.
get_season_incidence_curves <- function(cumulative_incidence_data){
  # Fit regularized incomplete beta function to get expected cum. fraction of cases for each week
  # Remove seasons with less than min_total_cases
  beta_fit_dataset <- cumulative_incidence_data %>% filter(season_total_cases >= min_total_cases)
  beta_model <- nls(cumulative_B_fraction ~
                      pbeta(q = cumulative_season_fraction, shape1, shape2),
                    data = beta_fit_dataset, start = initial_beta_parameters)
  
  shape1 <- coef(beta_model)[1]
  shape2 <- coef(beta_model)[2]
  
  # Init. tibble of within season intensity for 1982-2019 with the logistic fit to cumulative fraction of cases
  season_incidence_curves <- tibble(season_start_year = seq(1900,2019)) %>%
    rowwise() %>%
    mutate(n_weeks = n_iso_weeks(season_start_year)
    ) %>% ungroup()
  
  season_start_year <- c()
  season_week <- c()
  for(y in season_incidence_curves$season_start_year){
    week_range <- season_incidence_curves %>% 
      filter(season_start_year == y) %>% pull(n_weeks)
    week_range <- seq(1, week_range)
    season_week <- c(season_week, week_range)
    season_start_year <- c(season_start_year,
                        rep(y, length(week_range)))
  }
  
  season_incidence_curves <- tibble(season_start_year, season_week) %>%
    # Add calendar year and week number in calendar year (for consistency with US tibble)
    mutate(year = season_start_year, week = season_week) %>%
    # Add cumulative fraction of the season duration in each week 
    group_by(season_start_year) %>%
    mutate(cumulative_season_fraction = season_week / max(season_week))
  rm(list = c("season_start_year","season_week"))
  
  # For seasons with at least min_total_cases, fill and interpolate cumulative incidence values
  season_incidence_curves <- left_join(season_incidence_curves,
                                       cumulative_incidence_data %>%
                                                   filter(season_total_cases >= min_total_cases) %>%
                                                   mutate(cumulative_B_fraction_source = 'observed') %>%
                                                   select(season_start_year, season_week,
                                                          cumulative_B_fraction_source),
                                                 by = c("season_start_year", "season_week")) %>%
    mutate(cumulative_B_fraction_source = ifelse(is.na(cumulative_B_fraction_source),
                                                       'interpolated',
                                                       cumulative_B_fraction_source)) %>%
    rowwise() %>%
    mutate(cumulative_B_fraction =
             interpolate_cumulative_incidence(season_start_year, season_week,
                                              cumulative_incidence_data %>%
                                                # Fill or interpolate only for seasons with >= min_total_cases
                                                filter(season_total_cases >= min_total_cases))) %>%
    ungroup()
  
  # For seasons without sufficient data, use beta cumulative function fit
  season_incidence_curves <- season_incidence_curves %>%
    mutate(cumulative_B_fraction_source = ifelse(is.na(cumulative_B_fraction),
                                                 'beta function', cumulative_B_fraction_source)) %>%
    rowwise() %>%
    mutate(cumulative_B_fraction = ifelse(is.na(cumulative_B_fraction),
                                                pbeta(q = cumulative_season_fraction, shape1, shape2),
                                                cumulative_B_fraction)) %>%
    ungroup() %>%
    select(year, week, season_start_year, season_week, cumulative_season_fraction, cumulative_B_fraction,
           cumulative_B_fraction_source)
    return(season_incidence_curves)
}  

# Exporting intensity scores
ausnz_intensity_scores <- get_intensity_scores(surveillance_data_directory, case_data_directory, na_fill_method = 'average')

ausnz_intensity_scores <- ausnz_intensity_scores %>%
  # Set intensity to 0 before B_origin_year
  mutate(B_intensity = ifelse(season_start_year < B_origin_year,0, B_intensity),
         intensity_score = ifelse(season_start_year < B_origin_year,0, intensity_score))

# For additional fit to European gisaid data, assume intensity scores = 1
ausnz_intensity_scores <- bind_rows(ausnz_intensity_scores,
                                    ausnz_intensity_scores %>% filter(country == 'New Zealand') %>%
                                      mutate(country = 'Europe', B_intensity = NA, intensity_score = 1))

write.csv(ausnz_intensity_scores,'../results/processed_data/intensity_scores.csv', row.names = F)

# Exporting interpolated/estimated curves for fraction of cases experienced as time within seasons
cumulative_incidence_data_nz <- as_tibble(read.csv('../results/processed_data/raw_cumulative_incidence_curves/raw_cumulative_incidence_curves_nz.csv'))
cumulative_incidence_data_aus <- as_tibble(read.csv('../results/processed_data/raw_cumulative_incidence_curves/raw_cumulative_incidence_curves_aus.csv'))

season_cumulative_incidence_ausnz <- bind_rows(get_season_incidence_curves(cumulative_incidence_data_nz) %>%
                                                 mutate(country = 'New Zealand') %>% select(year, week, country, everything()),
                                               get_season_incidence_curves(cumulative_incidence_data_aus)%>%
                                                 mutate(country = 'Australia') %>% select(year, week, country, everything()))

# For additional fit to European gisaid data, assume season intensity curves equal to those of New Zealand
season_cumulative_incidence_ausnz <- bind_rows(season_cumulative_incidence_ausnz,
                                               season_cumulative_incidence_ausnz %>% filter(country == 'New Zealand') %>%
                                                 mutate(country = 'Europe'))

write.csv(season_cumulative_incidence_ausnz, '../results/processed_data/season_incidence_curves.csv', row.names = F)



# =============================================================================================
# =========================================== Plots ===========================================
# =============================================================================================
# Merged data for plots
# between_season_intensity_ausnz <- bind_rows(between_season_intensity_aus %>% 
#                                               select(country, season_start_year, B_fraction_in_ILI, peak_ILI_per_1000_visits, intensity_score) %>%
#                                               rename(ILI_incidence = peak_ILI_per_1000_visits), 
#                                             between_season_intensity_nz %>%
#                                               select(country, season_start_year, B_fraction_in_ILI, ILI_per_100K, intensity_score) %>%
#                                               rename(ILI_incidence = ILI_per_100K))

# Plot with B fraction within ILI over time
# B_fraction_in_ILI_ausnz <- ggplot(between_season_intensity_ausnz, aes(x = season_start_year, y = B_fraction_in_ILI)) +
#   geom_line() + geom_point() +
#   xlab('') + ylab('Fraction B+ among ILI specimens') +
#   facet_grid(country~.) +
#   scale_x_continuous(breaks = seq(1988,2016,4), limits = c(1988,2016)) +
#   theme(legend.position = 'top')

# Plot with ILI incidence over time
# total_ILI_ausnz <- ggplot(between_season_intensity_ausnz, aes(x = season_start_year, y = ILI_incidence)) +
#   geom_point() +
#   geom_line() +
#   facet_grid(country~., scales = 'free') +
#   ylab("ILI incidence metric") +
#   xlab("") +
#   theme(legend.position = 'top')+
#     scale_x_continuous(breaks = seq(1988,2016,4), limits = c(1988,2016))

# B_intensity_ausnz <- ggplot(between_season_intensity_ausnz,aes(y = intensity_score, x = season_start_year)) +
#   geom_point(size = 2) +
#   geom_line() +
#   facet_grid(country~.) +
#   ylab("Influenza B intensity score") +
#   xlab("Season (start year)") +
#   theme(legend.position = 'top', axis.text.y = element_text(size = 9)) +
#   scale_x_continuous(breaks = seq(1988,2016,4), limits = c(1988,2016))

#)
#
# plot(between_season_intensity_aus$intensity_score~ between_season_intensity_nz$intensity_score)

# Plot with cumulative fraction of cases over time

season_cumulative_incidence_ausnz_pl <- season_cumulative_incidence_ausnz %>%
  filter(cumulative_B_fraction_source != 'beta function', country != 'Europe') %>%
  ggplot(aes(x = cumulative_season_fraction, y = cumulative_B_fraction)) +
  geom_line(aes(color = factor(season_start_year))) +
  geom_point(aes(color = factor(season_start_year),
                 shape = factor(cumulative_B_fraction_source,
                                levels = c('observed','interpolated')))) +
  # Beta function predictions are the same for all years w/o data; choosing 2000
  geom_line(data = filter(season_cumulative_incidence_ausnz, year == 2000, country != 'Europe'), size = 1.5) +
  xlab('Fraction of the season') +
  ylab('Cumulative fraction of B cases') +
  scale_color_discrete(name = 'Season') +
  scale_linetype_discrete(name = 'Country') +
  scale_x_continuous(breaks = seq(0,1,0.1)) +
  scale_shape_discrete(name = 'Source') +
  facet_grid(country~.)
#
# save_plot('../figures/surveillance_data/season_cumulative_incidence_ausnz.pdf',
#           season_cumulative_incidence_ausnz_pl,
#            base_height = 6,base_width =9)



