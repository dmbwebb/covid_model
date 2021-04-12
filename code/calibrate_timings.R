
# Set seed ----------------------------------------------------------------

    set.seed(12345)

# CALIBRATE TIMING VARIABLES ----------------------------------------------
    
    load("data/processed/data_save.RData")    


    

# Symptoms ----------------------------------------------------------------
    
    data_save$params_symptom_timing <- list(meanlog = 1.63, sdlog = 0.5) # from BMJ McLoon et al
  
    
# Serial interval parameters ---------------------------------------------------------
    

    
    
    # From 
    # He, X., Lau, E. H. Y., Wu, P., Deng, X., Wang, J., Hao, X., Lau, Y. C., Wong, J. Y., Guan, Y., Tan, X., Mo, X., 
    # Chen, Y., Liao, B., Chen, W., Hu, F., Zhang, Q., Zhong, M., Wu, Y., Zhao, L., … Leung, G. M. (2020). 
    # Temporal dynamics in viral shedding and transmissibility of COVID-19. Nature Medicine, 26(5), 672–675. 
    # https://doi.org/10.1038/s41591-020-0869-5
    
    # https://github.com/ehylau/COVID-19
    
    # See their script file: Fig1c_Rscript.R, and data is Fig1c_data.xlsx
    
    # Import He et al 2020 data
    reference_date <- as.Date("2020-01-01")
    
    serial_intervals_raw <- read_excel("data/temporal dynamics in viral shedding/Fig1c_data.xlsx",
                                       col_types = c("numeric", "date", "date", "date")) %>% 
      mutate(across(c(x.lb, x.ub, y), 
                    ~ as.numeric(as.Date(.x) - reference_date)))
    
    
    rerun_secondary_timing_data <- TRUE
    
    
    
    if (rerun_secondary_timing_data) {
      
      print("Calculating parameters for infectiousness profile")
      
      # Function for inferring the distribution of the infectiousness profile
      p.Z  = function(lb, ub, gpar, lnpar) {
        
        length(lb)
        
        #--- infectiousness, gamma distribution ---
        # gpar[1:2]: hyper-parameters (gamma)
        # x        : infection time of infectee w.r.t onset time of infector
        f.Xc = function(x, gpar) { dgamma(x, gpar[1], gpar[2]) }
        
        #--- incubation, from Li et al NEJM 2020 ---
        # lnpar[1:2]: hyper-parameter (logNormal)
        # y         : length of incubation period of infectee
        f.Y  = function(y, lnpar) { dlnorm(y, lnpar[1], lnpar[2]) }
        
        #--- convolution between incubation and infectiousness profile ---
        # gpar[3]: shift c days before symptom onset of infector
        # z      : length of serial interval
        f.Z = function(z, gpar, lnpar) {
          integrate(
            f = function(x, z, gpar, lnpar) { f.Y(z+gpar[3]-x, lnpar)*f.Xc(x, gpar) },
            lower = -Inf, 
            upper = Inf,
            z     = z,
            gpar  = gpar,
            lnpar = lnpar
          )$value
        } 
        f.Z2 = Vectorize(f.Z, vectorize.args = "z")
        
        # print(lb)
        
        #--- p.Z ---
        integrate(
          f = function(x, gpar, lnpar) { f.Z2(x, gpar, lnpar) },
          lower = lb,
          upper = ub,
          gpar  = gpar,
          lnpar = lnpar
        )$value
      }
      p.Z2 = Vectorize(p.Z, vectorize.args = c("lb", "ub"))
      
      
      lli.fx = function(gpar, x.lb, x.ub, y, lnpar) {
        # print(x.ub)
        
        lli = log(p.Z2(lb = y-(x.ub+0.5), ub = y-(x.lb-0.5), gpar, lnpar))
        out <- -sum(lli)
        print(out)
        return(out)
        
      }
      
      #--- incubation period ---
      # ln.par1 = 1.434065
      # ln.par2 = 0.6612
      
      # Use BMJ values
      ln.par1 = data_save$params_symptom_timing$meanlog # values from BMJ McAloon et al 
      ln.par2 = data_save$params_symptom_timing$sdlog
      
      
      fit = optim(
        par = c(10, 0.5, 20), 
        fn = lli.fx, 
        method = "L-BFGS-B",
        lower  = 1e-3 + c(1,0,4),
        x.lb = serial_intervals_raw$x.lb, 
        x.ub = serial_intervals_raw$x.ub,
        y    = serial_intervals_raw$y,
        lnpar = c(ln.par1, ln.par2)
      )
      
      inf.par <- fit$par
      
      inf_profile_params <- list(
        shape = inf.par[[1]],
        rate = inf.par[[2]],
        shift = inf.par[[3]]
      )
      
    }
    
    # HARD RECORD VALUES (so don't have to rerun)
    # inf_profile_params <- list(
    #   shape = 20.51651,
    #   rate = 1.592124,
    #   shift = 12.27248
    # )
    
    
    # Plot to check
    # inf_profile_data <- rgamma(1000000,
    #                            shape = serial_parameters$shape,
    #                            rate = serial_parameters$rate) - inf_profile_params$shift
    
    
    # ggplot(data = tibble(x = inf_profile_data)) + 
    # geom_density(aes(x = x))
    
    
    
    
    # allow for negative serial intervals for (shifted) gamma distribution
    min.serial <- min((serial_intervals_raw$y-serial_intervals_raw$x.ub))
    
    # Log likelihood function for gamma distribution
    lli.g = function(par, x.lb, x.ub, y) {
      lli = log(pgamma(y - (x.lb - 0.5) - min.serial, par[1], par[2]) -
                  pgamma(pmax(y - (x.ub + 0.5) - min.serial, 0.1), par[1], par[2]))
      return(-sum(lli))
    }
    
    # Calculate gamma distribution parameters for serial distribution
    ser.fit <- optim(
      c(5, 1),
      lli.g,
      x.lb = serial_intervals_raw$x.lb,
      x.ub = serial_intervals_raw$x.ub,
      y = serial_intervals_raw$y
    )
    
    
    serial_parameters <- list(
      shape = ser.fit$par[1],
      rate = ser.fit$par[2],
      shift = min.serial - 0.5
    )
    
    # 8.123758      0.6361684       -7.5
    # serial_parameters <- ser.fit$par
    
    # serial_intervals_data <- tibble(
    #   serial = rgamma(1000000,
    #                   shape = serial_parameters$shape,
    #                   rate = serial_parameters$rate) + serial_parameters$shift
    # ) %>%
    #   hist_basic(x = serial)

    # Mean of gamma = shape / rate
    serial_mean <- (serial_parameters$shape / serial_parameters$rate) + serial_parameters$shift
    
    # serial_intervals_data$serial %>% mean()
    data_save$serial_mean <- serial_mean
    
    
    # Add to data_save
    data_save$params_serial <- serial_parameters
    
    
    
    

# Other timing parameters -------------------------------------------------


    
    
    # Allow multivariate skewed normal to be drawn with an error safely (e.g. if Omega is not positive definite)
    rmsn_safe <- safely(sn::rmsn)
    

    
    # FUNCTION for creating df of timings by drawing multivariate skewed normal 
    draw_timings <- function(params, n_draw, serial_params, symptom_params) {
      
      # FOR TESTING:
      # params <- c(3, 3, 1, 2, 0, 0, 0, 5); n_draw <- 1000
      # serial_params <- data_save$serial_parameters
      # symptom_params <- list(meanlog = 1.63, sdlog = 0.5)
      
      # Extract paramter names
      var_symp <- params[[1]]
      var_serial <- params[[2]]
      cov_serial_primary <- params[[3]]
      cov_serial_secondary <- params[[4]]
      alpha_primary <- params[[5]]
      alpha_secondary <- params[[6]]
      alpha_serial <- params[[7]]
      
      Omega <- matrix(
        c(var_symp,           0,                    cov_serial_primary,
          0       ,           var_symp,             cov_serial_secondary,
          cov_serial_primary, cov_serial_secondary, var_serial),
        nrow = 3
      )
      
      # Generate basic data by drawing from multivariate skewed normal
      gen_data <- rmsn_safe(
        n = n_draw,
        xi = c(0, 0, 0),
        Omega = Omega,
        alpha = c(alpha_primary, alpha_secondary, alpha_serial)
      )
      
      # If it throws an error, return infinity
      if (!is.null(gen_data$error)) {
        print(gen_data$error)
        return(Inf)
      } else {
        gen_data <- gen_data$result
      }
      
      
      # Calculate distributions (used in secondary timings function, in draw secondary cases)
      joint_dist <- sn::makeSECdistr(
        dp = list(xi = c(0, 0, 0), 
                  Omega = Omega, 
                  alpha = c(alpha_primary, alpha_secondary, alpha_serial)),
        family = "SN",
        compNames = c("primary_symp", "secondary_symp", "serial")
      )
      
      marg_dist_primary <- sn::marginalSECdistr(joint_dist, comp = 1)
      marg_dist_secondary <- sn::marginalSECdistr(joint_dist, comp = 2)
      marg_dist_serial <- sn::marginalSECdistr(joint_dist, comp = 3)
      
      
      # Convert the distributions
      gen_data_clean <- tibble(
        primary = gen_data[, 1],
        secondary = gen_data[, 2],
        serial = gen_data[, 3]
      ) %>% 
        # CONVERT TO Qs
        mutate(
          primary_symptoms_timing_q =  sn::psn(x = - primary, dp = marg_dist_primary@dp), # note the minus sign  - to generate negative correlation between primary symptoms and serial interval!!!
          secondary_symptoms_delay_q = sn::psn(x = secondary, dp = marg_dist_secondary@dp),
          serial_interval_q =          sn::psn(x = serial, dp = marg_dist_serial@dp)
        ) %>% 
        
        # CONVERT TO actual distributions
        mutate(
          primary_symptoms = qlnorm(primary_symptoms_timing_q, meanlog = symptom_params$meanlog, sdlog = symptom_params$sdlog),
          secondary_symptoms_delay = qlnorm(secondary_symptoms_delay_q, meanlog = symptom_params$meanlog, sdlog = symptom_params$sdlog),
          serial_interval = qgamma(p = serial_interval_q, shape = serial_params$shape, rate = serial_params$rate) + serial_params$shift
        ) %>% 
        mutate(
          secondary_timing = primary_symptoms + serial_interval - secondary_symptoms_delay
        )
      
      return(gen_data_clean)
      
    }
    
    
    
    infer_secondary_timing <- function(params, min_infection = 1, n_draw, serial_params, symptom_params, secondary_mean) {
      
      # params <- c(1, 2, 0.5, 0.5, 0, 0, 0)
      
      # FIRST DRAW
      gen_data_clean <- draw_timings(params, n_draw, serial_params, symptom_params)
      
      if (is_atomic(gen_data_clean)) return(Inf) # returns Infinity if draw_timings returns Infinity
      else gen_data_clean <- gen_data_clean %>% mutate(row_id = row_number())
      
      # Set up while loop - need to replace cases where secondary_timing is below min_infection
      n_replace <- sum(gen_data_clean$secondary_timing < min_infection | is.na(gen_data_clean$secondary_timing))
      gen_data_replaced <- gen_data_clean
      
      while (n_replace > 0) {
        print(n_replace / nrow(gen_data_clean))
        new_row_ids <- gen_data_replaced %>% filter(secondary_timing < min_infection) %>% .$row_id
        
        
        gen_data_new <- draw_timings(params, n_draw, serial_params, symptom_params) %>% 
          mutate(row_id = row_number()) %>% 
          filter(row_id %in% new_row_ids)
        
        
        gen_data_replaced <- gen_data_replaced %>% 
          filter(secondary_timing >= min_infection) %>% 
          bind_rows(gen_data_new) %>% 
          arrange(row_id)
        
        n_replace <- sum(gen_data_replaced$secondary_timing < min_infection | is.na(gen_data_replaced$secondary_timing))
        
        if (is.na(n_replace)) {
          browser() # bug if n_replace is NA
        }
        
      }
      
      # Set up gamma di
      gamma_shape <- params[[8]]
      
      secondary_timing_data <- rgamma(
        n = n_draw,
        shape = gamma_shape,
        scale = secondary_mean/gamma_shape    # set mean generation time to be equal to serial interval mean  (used to be 6)
      )
      
      # Kolomogorov-Smirnov test to match distributions
      ks_secondary <- suppressWarnings(ks.test(x = gen_data_replaced$secondary_timing,
                                               y = secondary_timing_data))
      
      test_stat <- abs(ks_secondary$statistic)
      
      # print(test_stat)
      return(test_stat)
      
    }
    
    
    # OPTIMISE USING THE FUNCTION
    library("nloptr")
    set.seed(12345)
    fit <- nloptr(
      x0 = c(3, 3, 1, 2, 0, 0, 0, 5),
      eval_f = infer_secondary_timing,
      lb = 1e-5 + c(0, 0, 0, 0,       # Omega matrix
                    -Inf, -Inf, -Inf, # alphas
                    0),               # gamma parameter
      opts = list(
        "algorithm" = "NLOPT_LN_SBPLX",
        "xtol_abs"=1.0e-10,
        "maxeval" = 200,
        "print_level" = 3
      ),
      n_draw = 10000,
      serial_params = serial_parameters,
      symptom_params = list(meanlog = 1.63, sdlog = 0.5),
      min_infection = 1,
      secondary_mean = data_save$serial_mean
    )
    
    
    
    

# Export to data_save -----------------------------------------------------

    # WHEN RERUNNING
    data_save$params_timing <- fit$solution
    
    # HARD CODED (avoid rerunning)
    # data_save$params_timing <- c(2.84771514,  4.38193344,  0.96935893,  2.12849543, -0.13965082,  0.02000124,  1.01396150,  2.77771473)

    save(data_save, file = "data/processed/data_save.RData")
    
    


    

# Generate timings distribution data --------------------------------------

  
    params_test <- data_save$params_timing
    serial_params <-  data_save$params_serial
    symptom_params <-  list(meanlog = data_save$params_symptom_timing$meanlog, sdlog = data_save$params_symptom_timing$sdlog)
    # params_test <- c(0.5, 0, 4, 5, 6)
    
    var_symp <- params_test[[1]]
    var_serial <- params_test[[2]]
    cov_serial_primary <- params_test[[3]]
    cov_serial_secondary <- params_test[[4]]
    alpha_primary <- params_test[[5]]
    alpha_secondary <- params_test[[6]]
    alpha_serial <- params_test[[7]]
    
    N <- 100000
    
    set.seed(12345)
    gen_data_clean <- draw_timings(params_test, n_draw = N, serial_params = serial_params, symptom_params = symptom_params) %>% 
      mutate(row_id = row_number())
    
    n_replace <- sum(gen_data_clean$secondary_timing < 1)
    n_replace / N
    gen_data_replaced <- gen_data_clean
    
    while (n_replace > 0) {
      new_row_ids <- gen_data_replaced %>% filter(secondary_timing < 1) %>% .$row_id
      
      gen_data_new <- draw_timings(params = params_test, n_draw = N, serial_params = serial_params, symptom_params = symptom_params) %>% 
        mutate(row_id = row_number()) %>% 
        filter(row_id %in% new_row_ids)
      
      gen_data_replaced <- gen_data_replaced %>% 
        filter(secondary_timing >= 1) %>% 
        bind_rows(gen_data_new) %>% 
        arrange(row_id)
      
      n_replace <- sum(gen_data_replaced$secondary_timing < 1)
      
      if (is.na(n_replace)) {
        browser()
      }
      
    }
    
    gen_data_clean <- gen_data_replaced
    
    save(gen_data_clean, file = "data/processed/gen_data_clean_timings.RData")
    
    
    
    
    # PLOT 2D DENSITY PLOTS
    # plot_density <- function(.data, x, y) {
    #   ggplot(.data, aes(x = {{x}}, y = {{y}})) + 
    #     geom_point(alpha = 0.05) + 
    #     stat_density_2d(aes(fill = ..level..), geom = "polygon")
    # }
    # 
    
    # plot_density(gen_data_clean, x = primary_symptoms, y = serial_interval)
    # 
    # plot_density(gen_data_clean, x = secondary_symptoms_delay, y = serial_interval)
    # 
    # plot_density(gen_data_clean, x = primary_symptoms, y = secondary_symptoms_delay)
    # 
    # plot_density(gen_data_clean, x = primary_symptoms, y = secondary_timing)
    # 
    # plot_density(gen_data_clean, x = secondary_symptoms_delay, y = secondary_timing)
    # 
    # plot_density(gen_data_clean, x = primary_symptoms, y = secondary_timing - primary_symptoms)
    
    
    # PLOT COMPARISONS OF DISTRIBUTIONS
    # gen_data_clean %>% hist_basic(x = secondary_timing - primary_symptoms)
    # gen_data_clean %>% hist_basic(x = secondary_timing, binwidth = 1)
    
    
    
    # 
    # compare_dist <- function(data, generated) {
    #   ggplot() + 
    #     geom_density(data = tibble(x = data), aes(x = x), colour = "indianred", linetype = "dashed") + 
    #     geom_density(data = tibble(x = generated), aes(x = x), colour = "skyblue")
    # }
    # 
    # 
    # 
    # 
    # 
    # secondary_timing_data <- rgamma(n = N,  
    #                                 shape = params_test[[8]],
    #                                 scale = 6/params_test[[8]])
    # 
    # 
    # compare_dist(data = secondary_timing_data,
    #              generated = gen_data_clean$secondary_timing)
    # 
    # 
    # 
    # 
    # compare_dist(data = rlnorm(10000, meanlog = 1.63, sdlog = 0.5),
    #              generated = gen_data_clean$primary_symptoms)   
    # 
    # compare_dist(data = rgamma(10000, shape = serial_parameters$shape, rate = serial_parameters$rate) + min.serial - 0.5,
    #              generated = gen_data_clean$serial_interval)   
    # 
    # 
    # 
    # # PLOT DISTRIBUTIONS OF SECONDARY TIMING AND SYMPTOMS
    # gen_data_clean %>% pivot_longer(c(secondary_timing, primary_symptoms)) %>% 
    #   ggplot(aes(x = value, fill = name)) + 
    #   geom_density(alpha = 0.6)  + 
    #   facet_wrap(~ name, nrow = 2) +
    #   scale_fill_discrete(labels = c("Symptom onset", "Secondary infection")) + 
    #   labs(
    #     x = "Days relative to start of infection",
    #     y = "Density"
    #     # fill = element_blank()
    #   ) + 
    #   theme_custom() + 
    #   coord_cartesian(xlim = c(0, 20)) + 
    #   theme(
    #     legend.position = c(.97, .85),
    #     legend.justification = "right",
    #     legend.background = element_rect(fill="white",
    #                                      size=0.5, linetype="solid", 
    #                                      colour ="darkgrey"),
    #     legend.title = element_blank(),
    #     strip.background = element_blank(),
    #     strip.text.y = element_blank(),
    #     strip.text.x = element_blank()
    #   )
    # 
    # 
    # ggsave("timings_dist.png", width = 6, height = 5.5, scale = 0.6)
    # 
    # 
    # 