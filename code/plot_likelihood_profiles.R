#!/usr/bin/env Rscript
# Plots likelihood profile
library(dplyr)
library(ggplot2)
library(cowplot)
theme_set(theme_cowplot())
library(viridis)
library(stringr)
library(reshape2)
library(fields)
#library(directlabels)

source('model_functions.R')

args <- commandArgs(trailingOnly = T)
main_directory <- args[1]
is_synthetic <- as.logical(args[2])

get_true_par_values <- function(true_parameter_values_path, profiled_parameters){
  par_file_lines <- readLines(true_parameter_values_path)
  # Read object model_name
  eval(parse(text = par_file_lines[1]))
  
  # Get model parameter names (from model functions script)
  model_par_names = get(paste0(model_name,'_model_par_names'))
  
  true_par_values <- str_split(par_file_lines[2], '=')[[1]][2]
  true_par_values <- str_replace(true_par_values, ' ','')
  true_par_values <- as.numeric(str_split(true_par_values, ',')[[1]])
  
  # Retain values of profiled parameters
  true_par_values <- true_par_values[model_par_names %in% profiled_parameters]
  return(true_par_values)
}

add_LRS_test <- function(lik_profile){
  # All column names before 'loglik' are the names of fixed parameters.
  parameters <- colnames(lik_profile)[1:(which(colnames(lik_profile) == 'loglik') -1)]
  # Check there's at most 2 parameters in the profile
  stopifnot(length(parameters) == 1 | length(parameters) == 2)
  
  # Calculate confidence interval based on LRS
  dfreedom <- length(parameters)
  mle_loglik <- max(lik_profile$loglik)
  
  # Add 'significant' variable to lik_profile object
  lik_profile <- lik_profile %>% 
    mutate(
      # Likelihood-ratio statistic
      LRS = 2 * (mle_loglik - loglik),
      # Chi-squared test
      p_value = 1 - pchisq(LRS, df = dfreedom),
      significant = ifelse(p_value < 0.05,T,F)
    ) 
  return(lik_profile)
}

# Interpolation function using loess (for bivariate profiles)
interpolate_profile_loess <- function(lik_profile, parameters, step_size = 0.005){
  xgrid <- seq(min(pull(lik_profile[parameters[1]])), max(pull(lik_profile[parameters[1]])), step_size)
  if(length(parameters) == 1){
    profile_loess = loess(loglik ~ x,
                          data = tibble(loglik = lik_profile$loglik,
                                        x = pull(lik_profile[parameters[1]])))
    grid = data.frame(x = xgrid)
  }else{
    stopifnot(length(parameters) == 2)
    profile_loess = loess(loglik ~ x * y,
                          data = tibble(loglik = lik_profile$loglik,
                                        x = pull(lik_profile[parameters[1]]),
                                        y = pull(lik_profile[parameters[2]])))
    ygrid <-  seq(min(pull(lik_profile[parameters[2]])), max(pull(lik_profile[parameters[2]])), step_size)
    grid = expand.grid(x = xgrid, y = ygrid)
  }
  
  interp_profile = predict(profile_loess, newdata = grid)
  
  if(length(parameters) == 1){
    interp_profile <- tibble(x = xgrid, loglik = interp_profile)
    names(interp_profile) <- c(parameters[1],'loglik')
  }else{
    # melt to long format
    interp_profile <- melt(interp_profile, measure.vars = 'loglik')
    # Return data to numeric form
    interp_profile[,1] <- as.numeric(str_sub(interp_profile[,1],
                                             str_locate(interp_profile[,1], "=")[1,1] + 1))
    interp_profile[,2] <- as.numeric(str_sub(interp_profile[,2],
                                             str_locate(interp_profile[,2], "=")[1,1] + 1))
    interp_profile = as_tibble(interp_profile)
    colnames(interp_profile) <- c(parameters, 'loglik')
  }
  return(interp_profile)
}

# Interpolation function using linear/bilinear interpolation
interpolate_profile_linear <- function(lik_profile, parameters, step_size = 0.001){
  # For bivariate profiles, use bilinear interpolation
  if(length(parameters) == 2){
    # Construct x, y and z objects for interp.surface (package fields) from observed surface
    x <- pull(unique(lik_profile[parameters[1]]))
    y <- pull(unique(lik_profile[parameters[2]]))
    z <- matrix(NA, nrow = length(x), ncol = length(y))
    for(i in 1:length(x)){
      for(j in 1:length(y)){
        llik <- lik_profile %>% filter((!!sym(parameters[1])) == x[i],
                                         (!!sym(parameters[2])) == y[j]) %>%
          pull(loglik)
        if(length(llik) > 0){
          z[i,j] <- llik
        }
      }
    }
    obj <- list(x =x, y = y, z = z)
    
    xgrid <- seq(min(pull(lik_profile[parameters[1]])), max(pull(lik_profile[parameters[1]])), step_size)
    ygrid <-  seq(min(pull(lik_profile[parameters[2]])), max(pull(lik_profile[parameters[2]])), step_size)
    
    loc <- make.surface.grid(list(xgrid,ygrid))
  
    interp_surface <- interp.surface( obj, loc)
    interp_surface <- as.surface( loc, interp_surface)
    # Convert interpolated surface into tibble
    interp_profile <- interp_surface$z
    rownames(interp_profile) <- interp_surface$x
    colnames(interp_profile) <- interp_surface$y
    interp_profile <- as_tibble(melt(interp_profile))
    names(interp_profile) <- c(parameters[1], parameters[2], 'loglik')
  }else{
    # For univariate profiles, use loess
    x <- pull(lik_profile[parameters[1]])
    y <- lik_profile$loglik
    interp_loglik <- approx(x,y, n = 1000)
    
    interp_profile <- tibble(interp_loglik$x, interp_loglik$y)
    names(interp_profile) <- c(parameters[1], 'loglik')
  }
  return(interp_profile)
}

plot_likprof <- function(profile_csv_path, true_parameter_values_path, convert_gammas = F){
  # Read csv with profile
  lik_profile <- as_tibble(read.csv(profile_csv_path), header = T)
  max_loglik = lik_profile %>% filter(loglik == max(loglik)) %>% pull(loglik)
  
  # All column names before 'loglik' are the names of fixed parameters.
  parameters <- colnames(lik_profile)[1:(which(colnames(lik_profile) == 'loglik') -1)]
  # Check there's at most 2 parameters in the profile
  stopifnot(length(parameters) == 1 | length(parameters) == 2)
  
  # If this is a bivariate profile with gammas's, convert to chi's
  if(convert_gammas){
    # Calculate chis from gammas
    lik_profile <- lik_profile %>% mutate(chi_VY = gamma_VY * chi_Y, 
                                          chi_YV = gamma_YV * chi_V,
                                          chi_AV = gamma_AV * chi_V,
                                          chi_AY = gamma_AY * chi_Y)
    parameters <- str_replace(parameters,'gamma','chi')
  }
    
  # If path to true parameter values provided (synth. data), retrieve them
  if(is.na(true_parameter_values_path) == F){
    true_par_values <- get_true_par_values(true_parameter_values_path, profiled_parameters = parameters)
  }
  
  # Interpolate
  # Use LOESS if bivariate profile, linear interpolation if univariate
  interpolate_profile_function <- ifelse(length(parameters) == 1, interpolate_profile_linear,
                                interpolate_profile_loess)
  interp_profile = interpolate_profile_function(lik_profile, parameters, step_size = 0.001)
  
  # Plots
  # If profile is a bi-variate profile
  if(length(parameters) == 2){
    par1_mle <- lik_profile %>% filter(loglik == max(loglik)) %>% pull(parameters[1])
    par2_mle <- lik_profile %>% filter(loglik == max(loglik)) %>% pull(parameters[2])
    
    lik_profile_pl <- ggplot(add_LRS_test(lik_profile),
           aes_string(x = parameters[1], y = parameters[2], z = "loglik")) +
      #geom_point(aes(fill = loglik, shape = significant), size = 2.5, alpha =0.9) +
      #scale_shape_manual(values = c(21,22)) +
      geom_point(aes(fill = loglik), size = 2.5, alpha =0.9, shape = 21) +
      scale_fill_viridis(name = 'Log-likelihood',option = 'viridis') +
      stat_contour(data = interp_profile, aes_string(x = parameters[1], y = parameters[2],
                                                     z = 'loglik'),
                   breaks = max_loglik - qchisq(0.95,2)/2,
                   color = 'black') +
      geom_point(data = filter(lik_profile, loglik == max(loglik)), color = 'red',
                 shape = 21, size = 5) +
      theme(legend.position = c(0.04,1.01),
            legend.title = element_text(size = 10),
            legend.text = element_text(size = 9),
            plot.margin = margin(20,5,5,20,'pt')) +
      guides(fill = guide_colourbar(barwidth = 9, barheight = 0.9,
                                    direction = 'horizontal', title.position = 'top')) 
    
    if(!is.na(true_parameter_values_path)){
      lik_profile_pl <- lik_profile_pl + 
        geom_point(aes(x=true_par_values[1], y = true_par_values[2]), color = 'blue')
    }
    
  }else{ # If profile is univariate
    CI <- interp_profile %>%  
      filter(loglik >= max_loglik - qchisq(0.95,1)/2) %>%
      pull(!!parameters[1])
    CI_limits = c(min(CI),max(CI))
    
    par1_mle <- lik_profile %>% filter(loglik == max(loglik)) %>% pull(parameters[1])
    
    lik_profile_pl <- ggplot(lik_profile,
           aes_string(x = parameters, y = 'loglik')) + 
      geom_line(data = interp_profile) +
      geom_point(data = lik_profile) +
      geom_hline(yintercept = max(lik_profile$loglik) - qchisq(0.95, 1)/2, linetype = 2) +
      geom_vline(xintercept = CI_limits, linetype = 2) +
      ggtitle(paste('MLE = ', round(par1_mle,2), '; 95% CI ',
                    round(CI_limits[1],2), '-',
                    round(CI_limits[2],2), '',
                    sep = '')
      ) + ylab('Log-likelihood')
    
    # If synthetic data, plot true parameter value as vertical line
    if(is.na(true_parameter_values_path) == F){
      lik_profile_pl <- lik_profile_pl +  
        geom_vline(xintercept = true_par_values, color = 'blue') +
        annotate("text", x = 1.04*true_par_values,
                 y = min(lik_profile$loglik),
                 label = paste(true_par_values)) 
    }
  }
  if('chi_VY' %in% parameters & 'chi_YV' %in% parameters){
    lik_profile_pl <- lik_profile_pl + geom_abline(slope = 1, intercept = 0, linetype = 2)
  }
  return(lik_profile_pl + theme(title = element_text(size = 12)))
} 

main <- function(main_directory, is_synthetic){
  profile_directories <- list.dirs(paste0(main_directory,'likelihood_profiles/'),
                                   recursive = F)
  profile_directories <- str_replace(profile_directories,'//','/')
  if(is_synthetic){
    true_parameter_values_path <- list.files(main_directory, pattern = 'pars.txt')
    true_parameter_values_path <- paste0(main_directory, true_parameter_values_path)
  }else{
    true_parameter_values_path <- NA
  }
  profile_csv_paths <- paste0(profile_directories, '/likelihood_profile.csv')
  
  # Parameter names  
  profiled_parameter_names <- str_split(profile_directories, '/')
  profiled_parameters <- c()
  for(split_list in profiled_parameter_names){
    profiled_parameters <- c(profiled_parameters, split_list[length(split_list)])
  }
  rm(profiled_parameter_names)
  profiled_parameters <- profiled_parameters[file.exists(profile_csv_paths)]
  
  # Plots
  profile_csv_paths <- profile_csv_paths[file.exists(profile_csv_paths)]
  profile_plots <- lapply(profile_csv_paths, FUN = plot_likprof,
              true_parameter_values_path = true_parameter_values_path)
  names(profile_plots) <- profiled_parameters

  main_panel <- plot_grid(
    profile_plots$chi_V + 
      xlab(expression("Homologous protection from prior B/Vic exposure, "~italic(chi[VV]))),
    profile_plots$chi_Y + 
      xlab(expression("Homologous protection from prior B/Yam exposure, "~italic(chi[YY]))),
    profile_plots$gamma_YV + 
      xlab(expression("Protection against B/Vic from B/Yam (relative to"~italic(chi[VV])~"),"~italic(gamma[YV]))),
    profile_plots$gamma_VY + 
      xlab(expression("Protection against B/Yam from B/Vic (relative to"~italic(chi[YY])~"),"~italic(gamma[VY]))),
    profile_plots$R_V + 
      xlab(expression("Imprinting protection against B/Vic, "~ italic(R[V]))),
    profile_plots$R_Y + 
      xlab(expression("Imprinting protection against B/Yam, "~ italic(R[Y]))),
    ncol = 2
  )
  
  save_plot(paste0(main_directory,
                   'likelihood_profiles/main_panel.pdf'),
            main_panel,
            base_width = 12, base_height = 13)
  
  rates_factors_plot <- plot_grid(profile_plots$beta1 + 
                                    xlab('Attack rate for 0-5-year-olds') +
                                    ylab('Log-likelihood'),
                                  profile_plots$beta2 + 
                                    xlab('Attack rate for 6-17-year-olds') +
                                    ylab(''),
                                  profile_plots$beta3 + 
                                    xlab('Attack rate for 18+ year-olds') +
                                    ylab(''),
                                  #profile_plots$reporting_factor_aus + 
                                  #  xlab('Reporting factor for children 0-4\n(Australia)'),
                                  profile_plots$reporting_factor_nz+ 
                                    xlab('Reporting factor for children 0-4 (New Zealand)'),
                                  ncol = 2)
  save_plot(paste0(main_directory,
                   'likelihood_profiles/attack_rates_and_reporting_factors.pdf'),
            rates_factors_plot ,
            base_width = 12, base_height = 8)
  
  
  # Panel with protection only from pre 1988 exposure
  pre1988_protection_list = list(profile_plots$gamma_AV +
                                   xlab(expression("Protection against B/Vic from pre-1988 exposure (relative to"~italic(chi[VV])~"),"~italic(gamma[AV]))),
                                 profile_plots$gamma_AY +
                                   xlab(expression("Protection against B/Vic from pre-1988 exposure (relative to"~italic(chi[YY])~"),"~italic(gamma[AY]))) +
                                   ylab('')
  )
  
  pre1988_protection <- plot_grid(plotlist = pre1988_protection_list,
                                  nrow = 1)
  
  save_plot(paste0(main_directory,'likelihood_profiles/pre1988_protection.pdf'),
            pre1988_protection,
            base_height = 5, base_width = 13)
  
  # chi_Y_vs_R_Y <-  profile_plots$R_Y_VS_chi_Y +
  #   ylim(0.1,1) +
  #        ylab(expression("      Homologous protection\nfrom prior B/Yam exposure, "~italic(chi[YY]))) +
  #        xlab(expression("Imprinting protection against B/Yam, "~ italic(R[Y]))) +
  #        theme(legend.position = c(0.03,0.095),
  #              legend.background = element_rect(fill='white', linetype = 1))
  # 
  # save_plot(paste0(main_directory,'likelihood_profiles/R_Y_VS_chi_Y_broad.pdf'),
  #           chi_Y_vs_R_Y + guides(shape = F),
  #           base_height = 5, base_width = 6)
  # 
  # chi_Y_vs_R_Y_narrow <-  profile_plots$R_Y_VS_chi_Y_narrow +
  #   ylab(expression("      Homologous protection\nfrom prior B/Yam exposure, "~italic(chi[YY]))) +
  #   xlab(expression("Imprinting protection against B/Yam, "~ italic(R[Y]))) +
  #   theme(legend.position = c(0.03,0.097),
  #         legend.background = element_rect(fill='white', linetype = 1))
  # 
  # save_plot(paste0(main_directory,'likelihood_profiles/imprinting_vs_homologous_protection_Yam.pdf'),
  #           chi_Y_vs_R_Y_narrow + guides(shape = F),
  #           base_height = 5, base_width = 6)


  
  
  
  # save_plot(paste0(main_directory,
  #                  'likelihood_profiles/gamma_VY_VS_gamma_AY.pdf'),
  #           profile_plots$gamma_VY_VS_gamma_AY +
  #             xlab(expression("Protection against B/Yam from prior B/Vic exposure, "~italic(gamma[VY]))) +
  #             ylab(expression("Protection against B/Yam from pre-1988 strains, "~italic(gamma[AY]))) +
  #             theme(legend.position = 'top',
  #                   legend.background = element_rect(fill='white', linetype = 1),
  #                   legend.text = element_text(size = 6)),
  #           base_width = 8, base_height = 7
  #           )
  # 
  
  # save_plot(paste0(main_directory,
  #                  'likelihood_profiles/chi_V_VS_gamma_VY.pdf'),
  #           profile_plots$chi_Y_VS_gamma_VY +
  #             xlab(expression("Homologous protection from prior B/Yam exposure, "~italic(chi[YY]))) +
  #             ylab(expression("Protection against B/Yam from prior B/Vic exposure, "~italic(gamma[VY]))) +
  #             theme(legend.position = 'top',
  #                   legend.background = element_rect(fill='white', linetype = 1),
  #                   legend.text = element_text(size = 6)),
  #           base_width = 8, base_height = 7
  # )
  
 # Bi-variate profiles with Vic imprinting protection vs. previous exposure, same for Yam
 # Plus bi-variate profile for heterologous cross-protection
  # main_panel <- plot_grid(
  #   profile_plots$R_V_VS_chi_V_narrow +
  #     #ylab('Homologous protection\nfrom previous B/Vic exposure') +
  #     ylab(expression("      Homologous protection\nfrom prior B/Vic exposure, "~italic(chi[VV]))) +
  #     xlab(expression("Imprinting protection against B/Vic, "~ italic(R[V]))) +
  #     theme(legend.position = c(0.55,0.095),
  #           legend.background = element_rect(fill='white', linetype = 1)),
  #   profile_plots$R_Y_VS_chi_Y_narrow +
  #     ylab(expression("      Homologous protection\nfrom prior B/Yam exposure, "~italic(chi[YY]))) +
  #     xlab(expression("Imprinting protection against B/Yam, "~ italic(R[Y]))) +
  #     theme(legend.position = c(0.55,0.095),
  #           legend.background = element_rect(fill='white', linetype = 1)),
  #   profile_plots$gamma_VY_VS_gamma_YV +
  #     ylab(expression("      Protection against B/Vic\nfrom prior B/Yam exposure, "~italic(chi[YV]))) +
  #     xlab(expression("Protection against B/Yam from prior B/Vic exposure, "~italic(chi[VY]))) +
  #     theme(legend.position = c(0.55,0.095),
  #           legend.background = element_rect(fill='white', linetype = 1)),
  #   nrow = 2,
  #   labels = c('A)','B)','C)'),
  #   label_size = 14
  # )
  # 
  # save_plot(paste0(main_directory,
  #                  'likelihood_profiles/main_panel.pdf'),
  #           main_panel,
  #           base_width = 11, base_height = 11)
  
 
  # save_plot(paste0(main_directory,
  #                  'likelihood_profiles/Yam_imprinting_VS_previous_Yam_exp_effect.pdf'),
  #           profile_plots$R_Y_VS_chi_Y +
  #             ylab('Homologous protection\nfrom previous B/Yam exposure') +
  #             xlab('Imprinting protection against B/Yam'),
  #           base_width = 8, base_height = 6)
  # 
  # # Profile with absolute cross-lineage protection
  # save_plot(paste0(main_directory,'likelihood_profiles/absolute_crosslineage_protection.pdf'),
  #           profile_plots$gamma_VY_VS_gamma_YV +
  #             ylab(expression("Protection against B/Vic from prior B/Yam exposure, "~italic(chi[YV]))) +
  #             xlab(expression("Protection against B/Yam from prior B/Vic exposure, "~italic(chi[VY]))),
  #           base_height = 7, base_width = 8)
  # 
  # 
  # 
  # 
  # 
  # 
  # save_plot(paste0(main_directory,'likelihood_profiles/relative_cross_lineage_protection.pdf'),
  #           plot_grid(profile_plots$gamma_VY + ylab('Log-likelihood') +
  #                       xlab(expression("Protection against B/Yam from prior B/Vic exposure, "~italic(gamma[VY]))),
  #                     profile_plots$gamma_YV + ylab('Log-likelihood') +
  #                       xlab(expression("Protection against B/Vic from prior B/Yam exposure, "~italic(gamma[YV])))),
  #           nrow = 1,
  #           base_height = 6, base_width = 12)
  # 
  # save_plot(paste0(main_directory,'likelihood_profiles/gamma_VY.pdf'),
  #           profile_plots$gamma_VY + ylab('Log-likelihood') +
  #                       xlab(expression("Protection against B/Yam from prior B/Vic exposure, "~italic(gamma[VY]))),
  #           nrow = 1,
  #           base_height = 5, base_width = 6)
  # save_plot(paste0(main_directory,'likelihood_profiles/gamma_YV.pdf'),
  #           profile_plots$gamma_YV + ylab('Log-likelihood') +
  #             xlab(expression("Protection against B/Vic from prior B/Yam exposure, "~italic(gamma[YV]))),
  #           nrow = 1,
  #           base_height = 5, base_width = 6)
  # 

}
main(main_directory, is_synthetic)



