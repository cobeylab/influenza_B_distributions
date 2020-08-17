library(dplyr)
library(ggplot2)
library(cowplot)
theme_set(theme_cowplot())
library(viridis)

main <- function(){
  # Note that maternal ab dur. here is in units of years, not weeks as in basic_parameters.R
  maternal_ab_duration <- 0.5
  
  # Age when school starts in the Netherlands
  school_start_age <- 4
  
  # Import B seroprevalence data:
  seroprevalence_data <- as_tibble(read.csv('../results/processed_data/netherlands_seroprevalence.csv', header = T) %>%
    select(age, n_positive, n_total, seroprevalence))
  
  prob_seronegative <- function(age, attack_rate_preschool, attack_rate,
                                maternal_ab_duration, school_start_age){
    # Calculate instantaneous attack rate from annual attack rate
    alpha_0 <- -log(1 - attack_rate_preschool)
    alpha_1 <- -log(1 - attack_rate)
    # To facilitate reading the code with the equation in mind:
    m <- maternal_ab_duration
    
    if(age <= school_start_age){
      prob <- exp(-alpha_0*(age - m)) * (1 - exp(-alpha_0))/alpha_0
    }else{
      prob <- exp(-alpha_0*(school_start_age - m)) * exp(-alpha_1*(age - school_start_age)) *
        (1 - exp(-alpha_1))/alpha_1
    }
    stopifnot(prob < 1 + 1e-7)
    return(prob)
  }
  # Test: ff 2 rates are the same and mAbs = 0, then P(A) = (1 - atk_rate)^a * (1 -exp(-alpha))/alpha
  alpha <- -log(1 - 0.15)
  stopifnot(abs(prob_seronegative(10,0.15,0.15,0,4) - (1 - 0.15)^10 * (1 - exp(-alpha))/alpha) < 1e-7)
  
  
  loglik_function <- function(age, n_seropositive,n_total,attack_rate_preschool, attack_rate,
                              maternal_ab_duration, school_start_age){
    
    naive_prob <- prob_seronegative(age = age, attack_rate_preschool = attack_rate_preschool,
                                    attack_rate = attack_rate, maternal_ab_duration = maternal_ab_duration,
                                    school_start_age = school_start_age)
    loglik <- dbinom(x = n_seropositive, size = n_total, prob = 1 - naive_prob, log = T)
    return(loglik)
  }

  # Master log-likelihood function for computing likelihood for different ages
  master_loglik_function <- function(serology_data, attack_rate_preschool, attack_rate,
                                     maternal_ab_duration, school_start_age){
    # Calculate log-likelihood for each age group, them sum
    total_loglik <- serology_data %>%
      rowwise() %>%
      mutate(loglik = loglik_function(age=age,
                                      n_seropositive = n_positive,
                                      n_total = n_total,
                                      attack_rate_preschool = attack_rate_preschool,
                                      attack_rate = attack_rate,
                                      maternal_ab_duration = maternal_ab_duration,
                                      school_start_age = school_start_age)) %>%
      ungroup() %>%
      summarise(total_loglik = sum(loglik, na.rm =T)) %>%
      pull(total_loglik)
    return(total_loglik)
  }
  
  # Construct likelihood curve over the interval 0-1 for attack rate in one-class model
  attack_rate <- seq(0.001,1,0.001)
  loglik_surface <- tibble(attack_rate = attack_rate) %>%
    rowwise() %>%
    mutate(
      loglik = master_loglik_function(serology_data = seroprevalence_data, attack_rate_preschool = attack_rate,
                                      attack_rate = attack_rate, maternal_ab_duration = maternal_ab_duration,
                                      school_start_age = school_start_age)
    ) %>%
    ungroup()
  mle_attack_rate_oneclass <- loglik_surface %>%  summarise(mle = attack_rate[loglik == max(loglik)]) %>% pull(mle)
  mle_loglik_oneclass <- loglik_surface %>%  summarise(mle_loglik_oneclass = max(loglik)) %>% pull(mle_loglik_oneclass)
  
  loglik_surface <- loglik_surface %>% 
    mutate(
      # Likelihood-ratio statistic
      LRS = 2 * (mle_loglik_oneclass - loglik),
      # Chi-squared test
      p_value = 1 - pchisq(LRS, df = 1),
      significant = ifelse(p_value < 0.05,T,F)
    ) 
  
  CI_limits_oneclass <- loglik_surface %>% filter(significant == F) %>%
    summarise(lower_limit = min(attack_rate),
              upper_limit = max(attack_rate))
  
  # -------------- two-class model
  twoclass_loglik_surface <- as_tibble(expand.grid(seq(0.01,0.5,0.01), seq(0.01,0.5,0.01))) %>%
    rename(attack_rate_preschool = Var1, attack_rate = Var2) %>%
    rowwise() %>%
    mutate(
      loglik = master_loglik_function(serology_data = seroprevalence_data,
                                               attack_rate_preschool = attack_rate_preschool,
                                               attack_rate = attack_rate,
                                               maternal_ab_duration = maternal_ab_duration,
                                               school_start_age = school_start_age)
    ) %>%
    ungroup() %>% 
    mutate(
      # Likelihood-ratio statistic
      LRS = 2 * (max(loglik) - loglik),
      # Chi-squared test
      p_value = 1 - pchisq(LRS, df = 2),
      significant = ifelse(p_value < 0.05,T,F)
    ) 
  
  twoclass_profile_atk_rate_preschool <- twoclass_loglik_surface %>%
    group_by(attack_rate_preschool) %>%
    summarise(loglik = max(loglik)) %>% 
    ungroup() %>%
    mutate(
      # Likelihood-ratio statistic
      LRS = 2 * (max(loglik) - loglik),
      # Chi-squared test
      p_value = 1 - pchisq(LRS, df = 1),
      significant = ifelse(p_value < 0.05,T,F)
    ) 
  CI_limits_twoclass_atk_rate_preschool <- twoclass_profile_atk_rate_preschool %>%
    filter(significant == F) %>%
    summarise(lower_limit = min(attack_rate_preschool),
              upper_limit = max(attack_rate_preschool))
  
  twoclass_profile_atk_rate <- twoclass_loglik_surface %>%
    group_by(attack_rate) %>%
    summarise(loglik = max(loglik)) %>% 
    ungroup() %>%
    mutate(
      # Likelihood-ratio statistic
      LRS = 2 * (max(loglik) - loglik),
      # Chi-squared test
      p_value = 1 - pchisq(LRS, df = 1),
      significant = ifelse(p_value < 0.05,T,F)
    ) 
  CI_limits_twoclass_atk_rate <- twoclass_profile_atk_rate %>%
    filter(significant == F) %>%
    summarise(lower_limit = min(attack_rate),
              upper_limit = max(attack_rate))
  
  twoclass_mle_attack_rate <- twoclass_loglik_surface %>% filter(loglik == max(loglik))

  
  # Comparison between one-class and two-class models
  AIC_oneclass <- 2 - 2*mle_loglik_oneclass
  AIC_twoclass <- 4 - 2*twoclass_mle_attack_rate$loglik
  
  AIC_twoclass - min(c(AIC_oneclass, AIC_twoclass))
  AIC_oneclass - min(c(AIC_oneclass, AIC_twoclass))
  
  
  # -------------- Plots -------------
  
  # Plot with likelihood surface
  B_attack_rate_lik <- ggplot(filter(loglik_surface, loglik != -Inf,
                                     attack_rate <0.3, attack_rate > 0.05),
         aes(x = attack_rate, y = loglik)) + geom_point(aes(colour = factor(significant))) +
    ylab('Log-likelihood') + xlab('Baseline attack rate') +
    geom_vline(xintercept = mle_attack_rate_oneclass, linetype = 2, color = 'red') +
    annotate('text',x = 0.25, y = -20, 
             label = paste('MLE = ', mle_attack_rate_oneclass, ' [',
                  CI_limits_oneclass$lower_limit, ' - ',
                  CI_limits_oneclass$upper_limit, ']',
                  sep = '')) +
    ggtitle('Constant attack rate model') +
    scale_x_continuous(limits = c(0.05,0.3), breaks = seq(0.05,0.3,0.05)) +
    scale_color_brewer(type = 'qual', palette = 6) +
    theme(legend.position = 'none', plot.title = element_text(size = 12))+
    geom_hline(yintercept = mle_loglik_oneclass - qchisq(0.95,1)/2) +
    geom_vline(xintercept = c(CI_limits_oneclass$lower_limit, CI_limits_oneclass$upper_limit), linetype = 2)
  
  # Plot with predicted vs. observed seroprevalence
  pred_vs_obs_seroprevalence_oneclass <- seroprevalence_data %>% rowwise() %>%
    mutate(pred_seroprevalence = prob_seronegative(age=age,
                                                 attack_rate_preschool = mle_attack_rate_oneclass,
                                                 attack_rate = mle_attack_rate_oneclass,
                                                 maternal_ab_duration = maternal_ab_duration,
                                                 school_start_age = school_start_age)) %>% 
    ungroup() %>%
    mutate(pred_seroprevalence = 1 - pred_seroprevalence) %>%
    ggplot(aes(x = age, y = seroprevalence)) + geom_point() +
    ggtitle('Constant attack rate model') +
    geom_line(aes(x = age, y = pred_seroprevalence), colour = 'red') +
    ylab('Fraction positive for B') + xlab('Age (years)') +
    scale_x_continuous(breaks = seq(1,7)) +
    scale_y_continuous(breaks = seq(0,0.7,0.1))
  
  # Plot with deviations from predicted seroprevalence
  deviations_pl_oneclass <- seroprevalence_data %>% rowwise() %>%
    mutate(pred_seroprevalence = prob_seronegative(age=age,
                                                   attack_rate_preschool = mle_attack_rate_oneclass,
                                                   attack_rate = mle_attack_rate_oneclass,
                                                   maternal_ab_duration = maternal_ab_duration,
                                                   school_start_age = school_start_age)) %>% 
    ungroup() %>%
    mutate(pred_seroprevalence = 1 - pred_seroprevalence) %>%
    mutate(deviation = seroprevalence - pred_seroprevalence) %>%
    ggplot(aes(x = age, y = deviation)) + geom_point() +
    geom_hline(yintercept = 0, linetype = 2) +
    ylab('Deviation from predicted seroprevalence')
    
  # Plot with predicted vs. observed seroprevalence (two-class model)
  pred_vs_obs_seroprevalence_twoclass <- seroprevalence_data %>% rowwise() %>%
    mutate(pred_seroprevalence = prob_seronegative(age=age,
                                                 attack_rate_preschool = 
                                                   twoclass_loglik_surface %>% 
                                                   filter(loglik == max(loglik)) %>% 
                                                   pull(attack_rate_preschool),
                                                 attack_rate = twoclass_loglik_surface %>% 
                                                   filter(loglik == max(loglik)) %>% 
                                                   pull(attack_rate),
                                                 maternal_ab_duration = maternal_ab_duration,
                                                 school_start_age = school_start_age)) %>%
    ungroup() %>%
    mutate(pred_seroprevalence = 1 - pred_seroprevalence) %>%
    ggplot(aes(x = age, y = seroprevalence)) + 
    geom_point() +
    geom_line(aes(x = age, y = pred_seroprevalence), colour = 'red') +
    ylab('Fraction positive for B') + xlab('Age (years)') +
    scale_x_continuous(breaks = seq(1,7)) +
    scale_y_continuous(breaks = seq(0,0.7,0.1))  +
    ggtitle('Two-class model')
  
  deviation_pl_twoclass <- seroprevalence_data %>% rowwise() %>%
    mutate(pred_seroprevalence = prob_seronegative(age=age,
                                                   attack_rate_preschool = 
                                                     twoclass_loglik_surface %>% 
                                                     filter(loglik == max(loglik)) %>% 
                                                     pull(attack_rate_preschool),
                                                   attack_rate = twoclass_loglik_surface %>% 
                                                     filter(loglik == max(loglik)) %>% 
                                                     pull(attack_rate),
                                                   maternal_ab_duration = maternal_ab_duration,
                                                   school_start_age = school_start_age)) %>%
    ungroup() %>%
    mutate(pred_seroprevalence = 1 - pred_seroprevalence) %>%
    mutate(deviation = seroprevalence - pred_seroprevalence)  %>%
    ggplot(aes(x = age, y = deviation)) + geom_point() +
    geom_hline(yintercept = 0, linetype = 2) +
    ylab('Deviation from predicted seroprevalence')
  
  # Combined plot
  combined_pl <- plot_grid(pred_vs_obs_seroprevalence_oneclass, pred_vs_obs_seroprevalence_twoclass,
                           deviations_pl_oneclass, deviation_pl_twoclass,
                            nrow = 2)
   
   pdf('../figures/baseline_attack_rate_seroprevalence.pdf', height = 12, width = 12)
   plot(combined_pl)
   dev.off()
}
mle_attack_rate_oneclass <- main()


