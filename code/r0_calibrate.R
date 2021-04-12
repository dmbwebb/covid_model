    # source("code/individual_model_functions.R")
    # load("data/processed/data_save.RData")
    load("data/processed/gen_data_clean_timings.RData")

# Test  -------------------------------------------------------------------
    
    # 
    # print(str_glue("Assuming k_scale_factor of {k_scale_factor}"))
    # 
    # data_save$k_scale_factor <- k_scale_factor
    # 
    # save(data_save, file = "data/processed/data_save.RData")

    
    RNGkind("L'Ecuyer-CMRG"); set.seed(12345)
    
    
    n_pop_total <- 100000
    k_matrix_basic <- data_save$k_matrix * (n_pop_total / data_save$k_matrix_pop)
    n_initial_cases <- round(data_save$group_props * n_pop_total / 5000) # 5000th of the pop is infected at the model start
    group_props <- c(
      round(data_save$group_props[1:3] * n_pop_total) / n_pop_total,
      1 - sum(round(data_save$group_props[1:3] * n_pop_total) / n_pop_total)
    )
    
    k_mobility_0 <- 1 %>% 
      map(~ list(k_matrix = k_matrix_basic * .x)) %>% 
      set_names(as.character(1:length(.)))
    
    

    
    sim_for_r0 <- outbreak_sims_policy_change_parallel(
      n_sims = n_sims,
      n_cores = n_cores,
      n_iterations = n_iterations,
      print_detail = TRUE,
      constant_params = list(
        n_pop_total = n_pop_total, 
        n_initial_cases = ceiling(data_save$group_props * n_pop_total * 0.001), 
        group_props = data_save$group_props, 
        dt_approx = 1,
        recov_val = 10, # people recover 10 days after symptoms show (cuts off only ~ 1% of infections at the upper end of timing distribution)
        params_timing = data_save$params_timing,
        params_symptom_timing = data_save$params_symptom_timing,
        params_serial = data_save$params_serial,
        hh_size_data = data_save$hh_data_bogota,
        test_delay_data = data_save$test_delay_data,
        # test_choice_delay_data = data_save$test_choice_delay_data %>% filter(i_group %in% 1:3), 
        # test_results_delay_data = data_save$test_results_delay_data %>% filter(i_group %in% 1:3), 
        test_sensitivity_data = data_save$test_sensitivity_data, 
        ct_delay_data = data_save$ct_delay_data,
        probs_self_test_df = data_save$probs_self_test,
        probs_isolate_symptoms_df = data_save$probs_isolate_symptoms,
        probs_isolate_test_df = data_save$probs_isolate_test,
        probs_isolate_ct_df =  data_save$probs_isolate_ct,
        p_contact_if_isolated_home = c(1, 1, 1, 1), 
        p_contact_traced = data_save$p_contact_traced,
        p_hh_quarantine = data_save$p_hh_quarantine,
        # k_matrix = data_save$k_matrix * k_scale_factor * (n_pop_total / data_save$k_matrix_pop),
        contact_dispersion = data_save$contact_dispersion,
        sar_home = data_save$sar_home_input, 
        sar_out = data_save$sar_out_input, 
        infectiousness_by_symptoms = data_save$infectiouness_by_symptoms,
        # alpha = c(0.1, 0.1, 0.1)
        alpha = c(0, 0, 0, 0)
      ),
      policy_start_times = 0,
      time_varying_params = k_mobility_0
    )
    
    
    
    
    # read_params(k_vec = k_matrix_to_vec(data_save$k_matrix),
    #             mu = rowSums(data_save$k_matrix) / (data_save$group_props * data_save$k_matrix_pop),
    #             group_props = data_save$group_props,
    #             n_pop = data_save$k_matrix_pop)
    
    

# Calculate implied exponential ----------------------------------------------------

    
    # Plot to check time window
    # plot(sim_for_r0)
    # # plot_by_i_group(sim_for_r0, prop = TRUE)
    # plot_detected_by_i_group(sim_for_r0)
    # plot_detected(sim_for_r0)
    
    
    
    new_cases <- sim_for_r0 %>% 
      group_by(sim_id, t) %>% 
      # sum up across i_groups
      summarise(across(c(n_new_cases, n_confirmed_cases), sum)) %>% 
      group_by(sim_id) %>% 
      mutate(
        new_cases_actual = n_new_cases,
        new_cases_detected =  n_confirmed_cases  - lag(n_confirmed_cases)
      ) %>% 
      mutate(
        across(c(new_cases_actual, new_cases_detected), list(log = log))
      ) %>% 
      ungroup
    
    
    min_t <- 20
    max_t <- 60
    
    
    # REGRESSION
    # (First try with actual cases, then use detected cases)
    reg <- lm(new_cases_actual_log ~ t + as.factor(sim_id), data = new_cases, 
              subset = t <= max_t & t >= min_t)
    
    
    # Coefficients
    r <- reg$coefficients[[2]]
    r_conf <- confint(reg)[2, ]
    
    (r_lab <- str_glue("r = {round(r, 3)}\n[{round(r_conf[[1]], 3)}, {round(r_conf[[2]], 3)}]"))
    
    # PLOT
    plot <- new_cases %>% 
      
      # filter(date_results >= ymd("2020-03-14")) %>% 
      # filter(new_cases != 0) %>% 
      
      ggplot(aes(x = t, y = new_cases_actual_log))  + 
      geom_point(aes(colour = sim_id), alpha = 0.2) + 
      geom_smooth(data = new_cases %>% filter(t <= max_t & t >= min_t), method = "lm") +
      theme_custom(panel.grid = element_blank()) + 
      labs(x = "Days since 1st March 2020", y = "ln(Daily Confirmed New Cases)") + 
      annotate(geom = "text", x = 0, y = 7.5, label = r_lab, vjust = "inward", hjust = "inward")
    # geom_text(label = r_lab)
    
    # Check number of susceptibles
    min_t <- 20
    max_t <- 60
    
    susceptible <- sim_for_r0 %>% filter(t >= min_t, t <= max_t) %>% 
      group_by(sim_id, t) %>% 
      summarise(prop_susceptible = sum(n_susceptible) / n_pop_total) %>% 
      group_by(t) %>% 
      summarise(prop_susceptible = mean(prop_susceptible))
    
    susceptible %$% mean(prop_susceptible)
    
    ggplot(susceptible, aes(x = t, y = prop_susceptible)) + 
      geom_line()
      
    
    
  
    print(plot)
    
    
    
# Import generation interval dist -----------------------------------------
    
    
    
    
    # ggplot(gen_data_clean, aes(x = secondary_timing)) + 
    #   geom_density(colour = "seagreen4", fill = "seagreen4", alpha = 0.5) + 
    #   scale_y_continuous(labels = scales::percent_format(accuracy = 1), minor_breaks = NULL) + 
    #   labs(y = "Density", x = "Generation Interval (days)") + 
    #   coord_cartesian(xlim = c(0, 20)) +
    #   theme_classic()
    
    
    # gen_data_clean$secondary_timing
    
    dens <- density(gen_data_clean$secondary_timing)
    
    
    
    
    gen_int <- tibble(x = dens$x, y = dens$y) %>% 
      filter(x > 0) 
    
    
    # ggplot(gen_int, aes(x = x, y = y)) + geom_line()
    
    
    
    
    
    
# Combine with exponential -----------------------------------------------
    
    
    r <- reg$coefficients[[2]]
    # r <- confint(reg)[2, ][[1]]
    
    integrand_df <- gen_int %>% 
      rename(t = x, g_x = y) %>% 
      mutate(exp_term = exp(-r * t),
             integrand = exp_term * g_x)
    
    
    
    
    
    
# Numerically integrate ---------------------------------------------------
    
    # Create a trapezoid function that approximates the integrand
    integrand_trapezoid <- approxfun(integrand_df$t,
                                     integrand_df$integrand)
    
    
    r0_inv <- integrate(integrand_trapezoid, min(integrand_df$t), max(integrand_df$t))$value
    
    r0 <- 1 / r0_inv
    
    r0_store[[as.character(outside_work_factor)]] <- list(r0 = r0, r = r, plot = plot)
    
    print(str_glue("R0 is estimated to be {round(r0, 2)}"))
    
    
    
    
    
    
    
    
    

# CALIBRATE for probability of testing ------------------------------------

    # sim_for_r0$time_series %>% 
    #   group_by(sim_id, t) %>% 
    #   # sum up across i_groups
    #   summarise(across(-c(i_group, prop_susceptible), sum)) %>% 
    #   view() %>% 
    #   mutate(n_new_cases_confirmed = n_confirmed_cases - lag(n_confirmed_cases)) %>% 
    #   relocate(t, n_new_cases, n_new_cases_confirmed, n_false_neg_self) %>% 
    #   summarise(across(c(n_new_cases, n_new_cases_confirmed, n_false_neg_self), sum, na.rm = TRUE)) %>% 
    #   # mutate(prop_confirmed = n_confirmed_cases / n_cases_cum) %>% 
    #   print_all
    # 
    
    # Check prob of false negs
    
    check_false_negs <- sim_for_r0 %>%
      group_by(sim_id, t) %>%
      # sum up across i_groups
      summarise(across(-c(i_group, policy, prop_susceptible), sum)) %>%
      mutate(prop_confirmed = eventually_confirmed_self / (eventually_false_neg_self + eventually_confirmed_self)) %>%
      select(-c(n_isolators:n_confirmed_cases), n_cases_cum, n_confirmed_cases) %>%
      mutate(prop_self_test_yn = all_self_test_yn / n_cases_cum) %>%
      tail()
    
    
    # about 22% of cases end up with false negatives... (given the timing of testing and symptoms)
    
    
    # sim_for_r0$time_series %>% relocate(n_false_neg_self) %>% 
    #   print_all
    
    