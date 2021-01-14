library(tidyr)
library(dplyr)
library(cowplot)

# Function for combining shape and overall intensity of seasons experienced in birth year
# (In the U.S., individuals born in year X may experience part of the X-1/X season)
get_first_seasons_data <- function(birth_year, country_name,season_incidence_curves,intensity_scores){
  # season_incidence_curves: object produced by calculate_intensity_scores.R with shape of each season
  # intensity scores: object procuced by calculate_intensity_scores.R with season-level intensity
  
  first_seasons_data <- left_join(season_incidence_curves %>% 
                                    filter(year == birth_year,country == country_name),
                                  intensity_scores  %>%
                                    filter(country == country_name) %>% 
                                    select(season_start_year, intensity_score),
                                  by = 'season_start_year') 

  # If birth year spans a season for which intensity score is not available, assume intensity score is 1
  first_seasons_data <- mutate(first_seasons_data,
                               intensity_score = ifelse(is.na(intensity_score), 1, intensity_score)) %>%
    # Add number of weeks in each year
    group_by(year) %>%
    mutate(n_weeks_in_year = max(week)) %>% 
    # Consider each week as a potential birth_week
    rename(birth_week = week, birth_year = year) %>%
    ungroup()
  return(first_seasons_data)
}

# Expected atck. rates in seasons experienced during YOB, accnt. for maternal Abs and the shape and intensity of seasons
birth_year_attack_rates <- function(base_attack_rate_preschool, birth_year, country_name,
                                    maternal_ab_duration, season_incidence_curves, intensity_scores, max_allowed_rate){
  # maternal_ab_duration: number of weeks post-birth when maternal Abs protect from infection
  # intensity_scores: tibble with season-level B intensity scores
  # season_incidence_curves: tibble with cumulative fraction of B cases vs. cumulative fraction of season
  if(maternal_ab_duration >= 40){warning('Chosen duration of maternal Abs >= 40 weeks')}
  
  # Intensity and shape of seasons experienced in birth year
  first_seasons_data <- get_first_seasons_data(birth_year, country_name, season_incidence_curves, intensity_scores)
  
  # Initalize tibble with effective birth week (calendar week + duration of maternal Abs) and year
  n_weeks_in_birth_year <- max(filter(first_seasons_data, birth_year == !!birth_year) %>% pull(birth_week))
  
  results <- tibble(country = country_name, birth_year = birth_year, birth_week = seq(1,n_weeks_in_birth_year)) %>%
    mutate(effective_birth_week = birth_week + maternal_ab_duration,
           # Adjust birth year
           effective_birth_year = ifelse(effective_birth_week <= n_weeks_in_birth_year, birth_year, birth_year + 1),
           # Adjust effective birth weeks that were outside original birth year
           effective_birth_week = ifelse(effective_birth_week > n_weeks_in_birth_year,
                                         effective_birth_week - n_weeks_in_birth_year, effective_birth_week)) %>%
    select(country, birth_year, birth_week, effective_birth_week, effective_birth_year)
  
  # Ignore effective birth weeks that fall into the next year (and thus season)
  
  results <- filter(results, effective_birth_year == birth_year)
  
  
  # Map data on shape and intensity from first_seasons_data onto tibble with effective b. years and weeks
  results <- left_join(results, 
                       first_seasons_data %>% rename(effective_birth_week = birth_week, effective_birth_year = birth_year) %>% 
                         select(country, effective_birth_year, effective_birth_week, season_start_year,
                                season_week, cumulative_B_fraction, intensity_score),
                       by = c('effective_birth_year','effective_birth_week', 'country')
  ) %>%
    # calculate instantaneous attack rate before adjustment by fraction experienced
    mutate(lambda = -log(1 - base_attack_rate_preschool) * intensity_score) %>%
    # Add column with attack rate conditional on each effective birth_week
    mutate(attack_rate = 1 - exp(-lambda*(1 - cumulative_B_fraction))) %>%
    # Truncate attack rates greater than max_allowed_rate
    mutate(attack_rate = ifelse(attack_rate > max_allowed_rate, max_allowed_rate, attack_rate)) %>%
    # Integrate over effective birth dates (assuming unif. distribution) to get expected attack rate in b. year
    group_by(season_start_year) %>%
    summarise(exp_attack_rate = mean(attack_rate))
  return(results)
}

# Function for getting attack rates for all seasons experienced from birth until specified obs. year
get_attack_rates <- function(beta1, beta2, beta3, birth_year, max_obs_year, country_name, maternal_ab_duration,
                             season_incidence_curves, intensity_scores, cutoff_age, school_start_age, oldest_atk_rate_age,
                             max_allowed_rate = 0.75){
  
  # Maximum possible age at max_obs_year
  max_age_at_obs_year <- max_obs_year - birth_year
  
  # School start age in country
  country_school_age <- school_start_age %>% filter(country == country_name) %>%
    pull(school_start_age)
  
  # Get attack rate for seasons experienced during year of birth
  byear_attack_rates <- birth_year_attack_rates(base_attack_rate_preschool = beta1, birth_year, country_name, maternal_ab_duration,
                                                season_incidence_curves, intensity_scores, max_allowed_rate) %>%
    #Add maximum age = 0
    mutate(max_age = 0, mean_age = 0) %>%
    rename(attack_rate = exp_attack_rate) 
  
  # Extract remaining seasons from intensity_scores tibble:
  remaining_seasons <- intensity_scores %>% filter(season_start_year >= birth_year, country == country_name) %>% 
    mutate(max_age = season_start_year - birth_year,
           mean_age = max_age - 0.5) %>%
    filter((season_start_year %in% byear_attack_rates$season_start_year) == F) %>%
    # Get values until max_age_at_obs_year or cutoff_age, whichever is smaller
    filter(max_age <= min(c(cutoff_age, max_age_at_obs_year))) %>% arrange(max_age) 
  
  # Baseline attack rates as a function of age
  baseline_attack_rates <- ifelse(remaining_seasons$max_age < country_school_age, beta1,
                                  ifelse(remaining_seasons$max_age < oldest_atk_rate_age, beta2, beta3))
  
  # Calculate yearly attack rates for remaining seasons
  remaining_attack_rates <- remaining_seasons %>%
    mutate(baseline_attack_rate = baseline_attack_rates) %>%
    mutate(lambda = -log(1 - baseline_attack_rate) * intensity_score) %>%
    mutate(attack_rate = 1 - exp(-lambda))
  
  # Combine b. year attack rate(s) with attack rates of subsequent seasons
  attack_rates <- bind_rows(byear_attack_rates, remaining_attack_rates) %>% 
    select(-mean_age) %>%
    # Truncate attack rates at max_allowed_rate, if intensity scores made them higher than that.
    mutate(attack_rate = ifelse(attack_rate > max_allowed_rate, max_allowed_rate, attack_rate),
           birth_year = birth_year, country = country_name) %>%
    select(birth_year, country, season_start_year, max_age, attack_rate)
  
  return(attack_rates)
}

# Functions to get cum. sum and prod. of a vector up to the previous element 
# (e.g. for 3rd element, get product of 1st and 2nd, arg "first" specifies value for 1st elem.)
previous_cumprod <- function(x, first){
  cprod <- first
  if(length(x) > 1){
    cprod <- c(cprod, cumprod(x)[-length(x)])
  }
  return(cprod)
}
previous_cumsum <- function(x, first){
  csum <- first
  if(length(x) > 1){
    csum <- c(csum, cumsum(x)[-length(x)])
  }
  return(csum)
}
stopifnot(all(previous_cumprod(1:5, first = 1) == c(1, 1, 2, 6, 24)))
stopifnot(all(previous_cumsum(1:5, first = 1) == c(1,1,3,6,10)))
# Functions to get cum. sum and prod. of all elements after each element
# (e.g. for 3rd element, get product of 4th, 5th, ..., arg "last" specifies value for last elem.)
following_cumprod <- function(x, last){
  if(length(x) > 1){
    cprod <- rev(cumprod(rev(x)))[-1]
    cprod <- c(cprod, last)
  }else{
    cprod <- last
  }
  return(cprod)
}
following_cumsum <- function(x, last){
  if(length(x) > 1){
    fsum <- rev(cumsum(rev(x)))[-1]
    fsum <- c(fsum, last)
  }else{
    fsum <- last
  }
  return(fsum)
}

stopifnot(all(following_cumprod(1:5, last = 0) == c(120,60,20,5,0)))
stopifnot(all(following_cumsum(1:5, last = 0) == c(14,12,9,5,0)))


# Function for merging attack rates tibble with lin. freqs. (country implicit in attack_rates object)
merge_atk_rates_and_linfreqs <- function(attack_rates, lineage_frequencies){
  country <- attack_rates$country[1]
  data <- left_join(attack_rates, lineage_frequencies %>% rename(season_start_year = year) %>%
                      mutate(country = as.character(country)) %>%
                      # Lineage frequencies calculated from all East Asia are used for China
                      mutate(country = ifelse(country == 'East Asia', 'China', country)) %>%
                      filter(country == !!country) %>%  select(season_start_year, fraction_yamagata,
                                                               fraction_victoria, fraction_ancestor),
                    by = 'season_start_year')
  return(data)
}

# Function for calculating probability of no infection with a given lineage (or with neither V nor Y) in a single season
# If calculating prob of neither V nor Y, protection[1] is for Vic and protection[2] is for Yam
# Used in calculation of phi values
P_no_infection <- function(lineage, season, protection, data){
  if(lineage == 'all_B'){
    stopifnot(is.null(protection))
    P_no_inf_in_season <- data %>% filter(season_start_year == season) %>% 
      mutate(P_no_inf_in_season = 1 - attack_rate)
  }else{
    if(lineage == 'VY'){
      stopifnot(length(protection) == 2)
      V_protection <- protection[1]
      Y_protection <- protection[2]
      P_no_inf_in_season <- data %>% filter(season_start_year == season) %>% 
        mutate(P_no_inf_in_season = 1 - attack_rate * (fraction_victoria * (1 - V_protection) + 
                                                         fraction_yamagata * (1 - Y_protection)))
    }else{
      if(lineage == 'V'){
        stopifnot(length(protection) == 1)
        P_no_inf_in_season <- data %>% filter(season_start_year == season) %>% 
          mutate(P_no_inf_in_season = 1 - attack_rate * fraction_victoria * (1 - protection))
      }else{
        if(lineage == 'Y'){
          stopifnot(length(protection) == 1)
          P_no_inf_in_season <- data %>% filter(season_start_year == season) %>% 
            mutate(P_no_inf_in_season = 1 - attack_rate * fraction_yamagata * (1 - protection))
        }
      }
    }
  }
  
  return(P_no_inf_in_season %>% pull(P_no_inf_in_season))
}

# Function for returning Phi(i,j) knowing Phi(i,j-1)
phi_recursion_end <- function(j, phi_j_minus_1, lineage, protection, data){
  # Uses the fact that Phi(i,j) = P(no inf. in season j) *  Phi(i,j-1)
  P_no_inf_in_j = P_no_infection(lineage, season = j, protection, data)
  return(exp(log(P_no_inf_in_j) + log(phi_j_minus_1)))
}

# Function for returning Phi(i,j) knowing Phi(i-,j)
# For the start recursion, phi_i_minus_1 is a tibble with values for different end seasons
phi_recursion_start <- function(i, phi_i_minus_1, lineage, protection, data){
  # Use the fact that Phi(i,j) = Phi_(i-1,j) / P(no inf in  season i-1)
  P_no_inf_in_i_minus_1 = P_no_infection(lineage, season = i-1, protection, data)
  
  phi_i <- phi_i_minus_1 %>% 
    filter(end_season >= i) %>%
    mutate(phi = exp(log(phi) - log(P_no_inf_in_i_minus_1)),
           start_season = start_season + 1)
  return(phi_i)
}

# Calculates phi for a given lineage for a range of start and end seasons
# All end seasons are considered from min_start_season to max_start_season
phi = function(lineage, min_start_season, max_start_season, max_end_season, protection, data, split_year){
  end_seasons_vector <- min_start_season:max_end_season
  
  pre_split_ends = end_seasons_vector[end_seasons_vector < split_year]
  if(length(pre_split_ends) > 0){
    pre_split_starts = min_start_season:max(pre_split_ends)
  }else{
    pre_split_starts = c()
  }
  
  post_split_ends = end_seasons_vector[end_seasons_vector >= split_year]
  if(max_start_season >= split_year){
    post_split_starts = max(split_year, min_start_season):max_start_season
  }else{
    post_split_starts = c()
  }
  
  # For any i,j pairs that are both prior to split, phi is simply 1
  if(length(pre_split_ends) > 0){
    P_pre_split <- as_tibble(expand.grid(start_season = pre_split_starts,
                                         end_season = pre_split_ends)) %>%
      filter(end_season >=start_season) %>% arrange(start_season, end_season) %>%
      mutate(lineage = lineage, phi = 1)
  }else{
    P_pre_split = c()
  }
  
  # If there are pairs contained entirely after the split, calculate phi for them...
  if(length(post_split_starts) > 0){
    earliest_split_start = post_split_starts[1]
    
    # ...starting with the first start season in the split:
    P_post_split <- P_no_infection(lineage, earliest_split_start,
                                   protection, data)
    # ...expanding to subsequent seasons:
    for(end in post_split_ends[-1]){
      P_post_split <- c(P_post_split, phi_recursion_end(end, P_post_split[length(P_post_split)],
                                                        lineage,protection, data ))
    }
    P_post_split = tibble(start_season = earliest_split_start,
                          end_season = post_split_ends,
                          lineage, phi = P_post_split)
    
    # ... then expanding to subsequent start seasons, if any, up to max_start_season
    if(max_start_season > earliest_split_start){
      for(start in (earliest_split_start):max_start_season){
        P_post_split <- bind_rows(P_post_split,
                                  phi_recursion_start(start, P_post_split %>%
                                                        filter(start_season == start-1),
                                                      lineage, protection, data))
      }
    }
    
    if(length(pre_split_starts) > 0){
      # Finally, initialize all intervals that start before the split and end after
      P_pre_post <- as_tibble(expand.grid(start_season = pre_split_starts,
                                          end_season = post_split_ends)) %>%
        filter(end_season >= start_season) %>% arrange(start_season, end_season) 
      
      # If those intervals exist, split year will be earliest post split start,
      # and ... phi for those intervals is simply phi from 1988 to j:
      P_pre_post <- left_join(P_pre_post,
                              P_post_split %>% filter(start_season == split_year) %>%
                                select(-start_season), by = 'end_season')
    }else{
      P_pre_post = c()
    }
  }else{
    P_post_split = c()
  }
  
  # Combine all tibbles
  P <- bind_rows(P_pre_split, P_post_split, P_pre_post) %>%
    arrange(start_season, end_season)
  return(P)
}

# Probability of no infection with lineage (incl. 'VY') after imprinting_lineage exp. in each possible season of 1st infection
P_no_inf_after_impr <- function(lineage, imprinting_lineage, birth_year, max_obs_year, protection, data, split_year){
  
  stopifnot(imprinting_lineage %in% c('V','Y','A'))
  if(lineage == 'VY'){
    stopifnot(length(protection) == 2)
  }else{
    stopifnot(lineage %in% c('V','Y'))
    stopifnot(length(protection) == 1)
  }
  
  prob_name = paste('P_no',lineage,'since',imprinting_lineage, sep = '_')
  
  # Initialize tibble with values when first infection is the last possible given the obs year (i.e., obs year - 1)
  # In all of those cases, prob. of **NO** infection after first infection is 1
  P <- tibble(season_1st_infection_i = birth_year:(max_obs_year-1)) %>%
    mutate(observation_year = season_1st_infection_i + 1) %>%
    mutate(!!prob_name := 1) %>%
    select(observation_year, everything()) %>%
    arrange(observation_year)
  
  # Calculate probabilities for remaining combinations of season 1st. infection and obs. years.
  if(max_obs_year > birth_year + 1){
    remaining_probs <- phi(lineage, min_start_season = birth_year + 1, max_start_season = max_obs_year,
                           max_end_season = max_obs_year - 1,protection, data, split_year) %>%
      mutate(season_1st_infection_i = start_season - 1) %>%
      mutate(observation_year = end_season + 1) %>%
      select(-start_season, -end_season) %>%
      rename(!!prob_name := phi) %>%
      select(observation_year, season_1st_infection_i, !!prob_name) %>%
      arrange(observation_year)
    P <- bind_rows(P,remaining_probs)
  }
  
  
  P <- P %>% arrange(observation_year, season_1st_infection_i)
  return(P)
}

# Prob of no VY infection BEFORE each season j after 1st. infection with A in season i, for all obs. years
P_no_VY_before_j_given_A <- function(birth_year, max_obs_year, protection, data, split_year){
  stopifnot(length(protection) == 2)
  
  # If there's only one observation year post first season experienced
  if(max_obs_year == birth_year + 1){
    final_P <- tibble(observation_year = max_obs_year, 
                      season_1st_infection_i = birth_year, subsequent_season_j = birth_year,
                      P_no_VY_before_j_given_A = 1)
  }else{
    # Initialize tibble with probs. when j = i+1 (in which case prob no inf. after i and before j = 1)
    # Cases with j = i will be ignored when computing final probabilities. Give then prob no infection before j = 1)
    P <- tibble(season_1st_infection_i = c(birth_year:(max_obs_year-1), birth_year:(max_obs_year-1)),
                subsequent_season_j = c(birth_year:(max_obs_year-1), birth_year:(max_obs_year-1) +1)) %>%
      mutate(P_no_VY_before_j_given_A = 1) %>%
      filter(subsequent_season_j != max_obs_year)
    
    if(max_obs_year > birth_year + 2){
      # Calculate probabilities for max_obs_year
      remaining_probs <- phi('VY', min_start_season = birth_year + 1, max_start_season = max_obs_year - 1,
                             max_end_season = max_obs_year -2, protection, data, split_year) %>%
        mutate(season_1st_infection_i = start_season -1, subsequent_season_j = end_season + 1) %>%
        rename(P_no_VY_before_j_given_A = phi) %>%
        select(season_1st_infection_i, subsequent_season_j, P_no_VY_before_j_given_A)
      
      P <- bind_rows(P,remaining_probs) %>% arrange(season_1st_infection_i, subsequent_season_j)
      
      final_P <- P %>% mutate(observation_year = max_obs_year) %>%
        select(observation_year, everything())
      
      # For the remaining obs. years
      for(obs_year in (birth_year+1):(max_obs_year-1)){
        final_P <- bind_rows(final_P,
                             P %>% filter(subsequent_season_j < obs_year) %>%
                               mutate(observation_year = obs_year) %>% select(observation_year, everything()))
      }
    }else{
      final_P <- bind_rows(P %>% mutate(observation_year = max_obs_year -1),
                           P %>% mutate(observation_year = max_obs_year)) %>%
        filter(subsequent_season_j < observation_year) %>%
        select(observation_year, everything())
    }
  }
  return(final_P %>% arrange(observation_year, season_1st_infection_i))
}

# Calculates prob of no V or no Y inf. AFTER each season j since 1st. inf. with A in season i, for all obs. years
P_no_inf_after_j_given_exp <- function(lineage, exposure, birth_year, min_obs_year, max_obs_year, protection, data, split_year){
  
  # This function's used to calculate probability of e.g. no Y since A in i and V in j
  # So two protection parameters are required
  stopifnot(length(protection) == 2)
  # an effective protection parameter is calculated
  # Under a multiplicative model: (1 - chi_AY)(1 - chi_VY) = 1 - (chi_AY + chi_VY - chi_AY*chi_VY)
  #effective_protection = protection[1] + protection[2] - protection[1]*protection[2]
  # Under a strongest-protection only model
  effective_protection = max(protection[1], protection[2])
  
  prob_name = paste('P_no', lineage, 'since_j_given',exposure, sep = '_')
  
  # If there's only one observation year post first season experienced
  if(max_obs_year == birth_year + 1){
    final_P <- tibble(observation_year = max_obs_year, 
                      season_1st_infection_i = birth_year, subsequent_season_j = birth_year,
                      !!prob_name := 1)
  }else{
    # Initialize tibble with probs. when subsequent season j = obs_year - 1 (thus prob of no inf. after j = 1)
    P <- tibble(observation_year = (birth_year+1):max_obs_year, 
                subsequent_season_j = (birth_year+1):max_obs_year - 1) %>%
      mutate(!!prob_name := 1)
    
    # Probs. for all possible subsequent seasons and observation years
    remaining_probs <- phi(lineage, min_start_season = birth_year + 1, max_start_season = max_obs_year - 1,
                           max_end_season = max_obs_year -1, protection = effective_protection, data, split_year) %>%
      mutate(subsequent_season_j = start_season -1, observation_year = end_season + 1) %>%
      rename(!!prob_name := phi) %>%
      select(observation_year, subsequent_season_j, !!prob_name)
    
    P <- bind_rows(P,remaining_probs) %>% arrange(observation_year, subsequent_season_j)
    
    final_P <- c()
    
    for(obs_year in min_obs_year:max_obs_year){
      obs_year_possible_is <- birth_year:(obs_year -1)
      final_P <- bind_rows(final_P, tibble(observation_year = obs_year, season_1st_infection_i = birth_year:(obs_year -1)))   
    }
    
    final_P <- left_join(final_P, P, by = c('observation_year')) %>%
      filter(subsequent_season_j >= season_1st_infection_i)
    
    final_P <- final_P %>% arrange(observation_year, season_1st_infection_i, subsequent_season_j)
  }
  return(final_P)
}


# Function for calculating probs. of exposure histories given b. year for multiple observation years
calculate_iprobs_byear <- function(birth_year, min_obs_year, max_obs_year, country_name,chi_VY, chi_YV, chi_AV, chi_AY,
                                   beta1, beta2, beta3, cutoff_age, maternal_ab_duration,
                                   lineage_frequencies, intensity_scores, season_incidence_curves, school_start_age, oldest_atk_rate_age){
  
  #print(paste(country_name, birth_year))
  
  # Read split year off of lineage frequency data
  split_year = lineage_frequencies %>% filter(fraction_yamagata == 0, fraction_victoria == 0) %>%
    summarise(M = max(year)) %>% pull(M) +1
  
  # Get attack rates for seasons following birth
  attack_rates <- get_attack_rates(beta1 = beta1, beta2 = beta2, beta3 = beta3,
                                   birth_year = birth_year, max_obs_year = max_obs_year, country_name = country_name,
                                   maternal_ab_duration = maternal_ab_duration,
                                   season_incidence_curves = season_incidence_curves, intensity_scores = intensity_scores, 
                                   cutoff_age = cutoff_age, school_start_age = school_start_age,
                                   oldest_atk_rate_age = oldest_atk_rate_age)
  
  # Merge attack rates and lineage frequencies:
  data <- merge_atk_rates_and_linfreqs(attack_rates, lineage_frequencies)
  
  
  # If seasons were experienced other than the observation year  
  if(min(data$season_start_year) < max_obs_year){
    
    # Calculate fully naive probabilities
    naive_probs <- data %>% mutate(
      P_0 = cumprod(1 - attack_rate)
    ) %>% mutate(observation_year = season_start_year + 1) %>%
      select(observation_year, birth_year, P_0) %>%
      filter(observation_year <= max_obs_year, observation_year >= min_obs_year)
    
    # Calculate probabilities of imprinting by season of first exposure
    first_infection_probs <- data %>%
      mutate(
        P_no_exp_before_i = exp(cumsum(log(1 - attack_rate)) - log(1 - attack_rate)),
        P_V_impr_in_season_i = exp(log(attack_rate) + log(fraction_victoria) + log(P_no_exp_before_i)),
        P_Y_impr_in_season_i = exp(log(attack_rate) + log(fraction_yamagata) + log(P_no_exp_before_i)),
        P_A_impr_in_season_i = exp(log(attack_rate) + log(fraction_ancestor) + log(P_no_exp_before_i))
      ) %>%
      rename(season_1st_infection_i = season_start_year) %>%
      filter(season_1st_infection_i < max_obs_year) %>%
      select(country, birth_year, season_1st_infection_i, matches('P_'))
    
    # Calculate probability of no subsequent V or Y exposure after imprinting with A in each possible season
    if(birth_year < split_year){
      P_no_VY_given_A <- P_no_inf_after_impr('VY','A', birth_year, max_obs_year, protection = c(chi_AV, chi_AY),
                                             data, split_year) %>% filter(observation_year >= min_obs_year)
    }else{
      # Skip calculation if birth year after split (no A imprinting possible)
      P_no_VY_given_A <- as_tibble(expand.grid(observation_year = (birth_year+1):max_obs_year,
                                               season_1st_infection_i = birth_year:(max_obs_year -1))) %>%
        arrange(observation_year) %>% filter(observation_year > season_1st_infection_i, observation_year >= min_obs_year) %>% 
        mutate(P_no_VY_since_A = 1)
    }
    
    # Calculate probability of no subsequent Y exposure after imprinting with V in each possible season 
    P_no_Y_given_V <- P_no_inf_after_impr('Y','V', birth_year, max_obs_year, protection = chi_VY, data, split_year) %>%
      filter(observation_year >= min_obs_year)
    
    # Calculate probability of no subsequent V exposure after imprinting with Y in each possible season 
    P_no_V_given_Y <- P_no_inf_after_impr('V','Y', birth_year, max_obs_year, protection = chi_YV, data, split_year) %>%
      filter(observation_year >= min_obs_year)
    
    # For all seasons of 1st infection i and subsequent seasons j for all obs. years, calculate:
    # ...probability no V and no Y exposures before j given A exposure in i
    subseq_season_probs <- P_no_VY_before_j_given_A(birth_year, max_obs_year,protection = c(chi_AV, chi_AY),
                                                    data, split_year) %>%
      filter(observation_year >= min_obs_year)
    # ...probability of no Y exposure after j given A exposure in i and V exposure in j
    subseq_season_probs <- left_join(subseq_season_probs,
                                     P_no_inf_after_j_given_exp('Y', 'AV', birth_year, min_obs_year, max_obs_year,
                                                                protection = c(chi_AY, chi_VY), data, split_year),
                                     by = c("observation_year", "season_1st_infection_i", "subsequent_season_j"))
    # ... and probability of no V exposure after j given A exposure in i and Y exposure in j
    subseq_season_probs <- left_join(subseq_season_probs,
                                     P_no_inf_after_j_given_exp('V', 'AY', birth_year, min_obs_year, max_obs_year,
                                                                protection = c(chi_AV, chi_YV), data, split_year),
                                     by = c("observation_year", "season_1st_infection_i", "subsequent_season_j"))
    # Add attack rates and lineage frequencies during each season j
    subseq_season_probs <- left_join(subseq_season_probs,
                                     data %>% rename(subsequent_season_j = season_start_year,
                                                     attack_rate_j = attack_rate) %>%
                                       select(subsequent_season_j, attack_rate_j, matches('fraction')),
                                     by = 'subsequent_season_j')
    
    # Integrate probabilities conditional on A exposure in i across all j
    probs_given_A_impr <- subseq_season_probs %>% 
      # Set probability of second infection in j to 0 if j = season of 1st infection
      mutate(
        attack_rate_j = ifelse(season_1st_infection_i == subsequent_season_j,0,attack_rate_j),
        prob_V_in_j_given_A = exp(log(attack_rate_j) + log(fraction_victoria) + log(1 - chi_AV)),
        prob_Y_in_j_given_A = exp(log(attack_rate_j) + log(fraction_yamagata) + log(1 - chi_AY))
      ) %>%
      group_by(observation_year, season_1st_infection_i) %>%
      summarise(
        P_V0_after_A_in_i = sum(exp(log(P_no_VY_before_j_given_A) + log(prob_V_in_j_given_A) + log(P_no_Y_since_j_given_AV))),
        P_Y0_after_A_in_i = sum(exp(log(P_no_VY_before_j_given_A) + log(prob_Y_in_j_given_A) + log(P_no_V_since_j_given_AY))),
        P_V_then_Y_after_A_in_i = sum(P_no_VY_before_j_given_A * prob_V_in_j_given_A * (1 - P_no_Y_since_j_given_AV)),
        P_Y_then_V_after_A_in_i = sum(P_no_VY_before_j_given_A * prob_Y_in_j_given_A * (1 - P_no_V_since_j_given_AY))
      )
    
    # Merge probabilities of the different 1st exposures with history probabilities conditional on them
    joint_probs <- left_join(first_infection_probs, P_no_VY_given_A, by = "season_1st_infection_i") %>%
      select(country, birth_year, observation_year, season_1st_infection_i, everything()) 
    joint_probs <- left_join(joint_probs, P_no_Y_given_V, by = c("season_1st_infection_i","observation_year"))
    joint_probs <- left_join(joint_probs, P_no_V_given_Y, by = c("season_1st_infection_i","observation_year"))
    joint_probs <- left_join(joint_probs, probs_given_A_impr, by = c("season_1st_infection_i","observation_year"))%>%
      arrange(observation_year, season_1st_infection_i)
    
    # Compute final probabilities
    joint_probs <- joint_probs %>%
      group_by(observation_year) %>%
      summarise(P_A0 = sum(P_A_impr_in_season_i * P_no_VY_since_A),
                P_AV0 = sum(P_A_impr_in_season_i * P_V0_after_A_in_i),
                P_AY0 = sum(P_A_impr_in_season_i * P_Y0_after_A_in_i),
                P_A_V_then_Y = sum(P_A_impr_in_season_i * P_V_then_Y_after_A_in_i),
                P_A_Y_then_V = sum(P_A_impr_in_season_i * P_Y_then_V_after_A_in_i),
                P_V0 = sum(P_V_impr_in_season_i * P_no_Y_since_V),
                P_Y0 = sum(P_Y_impr_in_season_i * P_no_V_since_Y),
                P_VY = sum(P_V_impr_in_season_i * (1 - P_no_Y_since_V)),
                P_YV = sum(P_Y_impr_in_season_i * (1 - P_no_V_since_Y))
      ) %>%
      # Calculate marginal P_AVY (i.e., regardless of order)
      mutate(P_AVY = P_A_V_then_Y + P_A_Y_then_V) %>%
      select(-P_A_V_then_Y, - P_A_Y_then_V)
    # Add probability of being fully naive
    joint_probs <- left_join(naive_probs, joint_probs, by = 'observation_year')
    
    # Individuals observed in the season starting in their birth year are assumed to be fully naive
    birth_year_row = c(birth_year, birth_year, 1, rep(0,length(names(joint_probs)[-1]) -2))
    names(birth_year_row) = names(joint_probs)
    joint_probs <- bind_rows(as_tibble(rbind(birth_year_row)), joint_probs)
    
    joint_probs <- joint_probs %>% mutate(country = country_name) %>% select(country, everything())
    
  }else{
    # If max. obs year is the first and only season experienced in life, individual is fully naive
    joint_probs <- tibble(country = country_name, observation_year = max_obs_year, birth_year = birth_year,
                          P_0=1, P_A0=0, P_AV0=0, P_AY0=0, P_V0=0, P_Y0=0, P_VY=0, P_YV=0, P_AVY=0) 
  }
  
  # Check probabilities sum to 1
  prob_sums <- joint_probs %>% group_by(observation_year) %>%
    mutate(S =  P_0 + P_A0 + P_AV0 + P_AY0 + P_V0 + P_Y0 + P_VY + P_YV + P_AVY) %>% pull(S)
  stopifnot(all(abs(prob_sums - 1) < 1e-7))
  return(joint_probs)
}


# Imprinting probs. when only the age, not year of birth, is known
# (Integrates over possible years of birth)
calculate_iprobs_age <- function(age, observation_year, country_name, iprobs_by_byear){
  age <- as.numeric(age)
  # Given age and obs. year, 2 birth years are possible. Assume they're equally likely
  possible_byears <- c(observation_year - age - 1 , observation_year - age)
  
  probs <- iprobs_by_byear %>% filter(country == country_name, birth_year %in% possible_byears,
                                      observation_year == !!observation_year) %>%
    select(matches('P_')) %>%
    summarise_all(mean) %>%
    mutate(country = country_name, age = age, observation_year = observation_year) %>%
    select(age, country, observation_year, matches('P_'))
  # Check if probs sum to 1
  stopifnot(abs(rowSums(probs %>% select(matches('P_'))) - 1) < 1e-7)
  return(probs)
}

# Master function for calculating & distributing probs given a tibble with combined demog. & case data
calculate_iprobs  <- function(dem_plus_case_data, lineage_frequencies, intensity_scores, chi_VY, chi_YV,
                              chi_AV, chi_AY, beta1, beta2, beta3, cutoff_age, maternal_ab_duration,
                              season_incidence_curves, school_start_age, oldest_atk_rate_age, birth_year_cutoff){
  # Vector of observation years in the data
  observation_years <- unique(dem_plus_case_data$observation_year)
  
  # Min. and Max. observation year
  max_obs_year = max(observation_years)
  min_obs_year = min(observation_years)
  
  # Observation years for case data reported by age
  obs_years_age_data <- unique(filter(dem_plus_case_data, cohort_type == 'age') %>% pull(observation_year))
  
  # Range of birth years in the data by country
  birth_year_range <- dem_plus_case_data %>%
    mutate(minimum_possible_byear = ifelse(cohort_type == 'birth_year', cohort_value,
                                           observation_year - cohort_value -1)) %>%
    group_by(country) %>%
    summarise(
      byear_llim = min(minimum_possible_byear),
      byear_ulim = max(minimum_possible_byear)
    ) %>% 
    mutate(country = as.character(country)) %>% ungroup()
  
  # Vectors to pass to mapply for calculating probs. by country/birth year combination
  country_vector <- c()
  byear_vector <- c()
  for(i in 1:nrow(birth_year_range)){
    llim = max(birth_year_range$byear_llim[i], birth_year_cutoff)
    
    country_vector <- c(country_vector, rep(birth_year_range$country[i],
                                            birth_year_range$byear_ulim[i] - llim + 1))
    byear_vector <- c(byear_vector, seq(llim, birth_year_range$byear_ulim[i]))
  }
  
  # Combinations of country/birth-year
  country_byear_combinations <- tibble(
    country = country_vector,
    birth_year = byear_vector
  )
  rm("country_vector", "byear_vector")
  
  
  # Imprinting probabilities by country/birth-year for the range of birth years observed
  iprobs_by_byear <- mapply(FUN = calculate_iprobs_byear,
                            birth_year = country_byear_combinations$birth_year,
                            country_name = country_byear_combinations$country,
                            MoreArgs = list(min_obs_year = min_obs_year,
                                            max_obs_year = max_obs_year, chi_VY = chi_VY, chi_YV = chi_YV,
                                            chi_AV = chi_AV, chi_AY = chi_AY,
                                            beta1 = beta1, beta2 = beta2, beta3 = beta3,
                                            cutoff_age = cutoff_age,maternal_ab_duration = maternal_ab_duration,
                                            lineage_frequencies = lineage_frequencies, intensity_scores = intensity_scores,
                                            season_incidence_curves = season_incidence_curves,
                                            school_start_age = school_start_age, oldest_atk_rate_age = oldest_atk_rate_age),
                            SIMPLIFY = F)
  
  
  iprobs_by_byear <- bind_rows(iprobs_by_byear)

  # If range of birth years extends before birth_year_cutoff, copy history probabilities at birth_year_cutoff
  if(birth_year_range$byear_llim[1] < birth_year_cutoff){
    for(y in birth_year_range$byear_llim[1]:(birth_year_cutoff-1)){
      iprobs_by_byear <- bind_rows(iprobs_by_byear,
                iprobs_by_byear %>% filter(birth_year == birth_year_cutoff, observation_year != birth_year) %>%
                  mutate(birth_year = y)) 
    }
  }

  # Calculate imprinting probabilies when only age is given
  iprobs_by_age <- dem_plus_case_data %>% filter(cohort_type == 'age') %>%
    select(country, observation_year, cohort_type, cohort_value) %>%
    distinct() %>%
    rowwise() %>%
    do(
      calculate_iprobs_age(age = .$cohort_value, observation_year = .$observation_year,
                           country_name = .$country, iprobs_by_byear = iprobs_by_byear)
    ) %>%
    ungroup()
  rm()
  # Distribute imprinting probabilities for rows with cases by birth year
  combined_data_byear <- left_join(dem_plus_case_data %>% filter(cohort_type == 'birth_year'),
                                   iprobs_by_byear %>% rename(cohort_value = birth_year),
                                   by = c('country', 'observation_year', 'cohort_value'))
  
  if(length(obs_years_age_data) > 0){
    # Distribute imprinting probabilities for rows with cases by age
    combined_data_age <- left_join(dem_plus_case_data %>% filter(cohort_type == 'age'),
                                   iprobs_by_age %>% rename(cohort_value = age),
                                   by = c('country', 'observation_year', 'cohort_value')
    )
    # Combine all data
    combined_data <- bind_rows(combined_data_age, combined_data_byear)
  }else{
    combined_data <- combined_data_byear
  }
  
  combined_data <- combined_data %>%
    rename(obs_cases = n_cases) %>%
    select(matches('SurveillanceType'), country, region, observation_year, cohort_type, cohort_value, lineage, obs_cases, CLY_total_cases, rel_pop_size, matches('P_'))
  return(combined_data)
}