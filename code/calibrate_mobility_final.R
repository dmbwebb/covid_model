
# Parameters --------------------------------------------------------------

    RNGkind("L'Ecuyer-CMRG"); set.seed(12345)
    
    n_pop_total <- 100000
    
    # Calculate k, initial cases, group_props based on n_pop_total
    k_matrix_basic <- data_save$k_matrix * (n_pop_total / data_save$k_matrix_pop)
    n_initial_cases <- round(data_save$group_props * n_pop_total / 5000) # 5000th of the pop is infected at the model start
    group_props <- c(
      round(data_save$group_props[1:3] * n_pop_total) / n_pop_total,
      1 - sum(round(data_save$group_props[1:3] * n_pop_total) / n_pop_total)
    )
    
    policy_start_times_0 <- 0
    mobility_factors_0 <- 1
    
    k_mobility_0 <- mobility_factors_0 %>% 
      map(~ list(k_matrix = k_matrix_basic * .x)) %>% 
      set_names(as.character(1:length(.)))
    
    lower_bound_date <- ymd("2020-06-01")
    # n_loops = 15                 # 10
    # n_sims_per_round <- 35       # 35
    detection_delay = 21
    roll_smooth = 14
    # t_shift = 5
    # inc_shift = 0.1
    # inc_length = 21
    # conf_level <- 0.
    shift_scale = 1/3    # conversion from log deviation to change in mobility (e.g. 0.5 log deviation becomes 0.1 change in mobility factor)
    alteration_min = 0.02 # lower bound when stops altering
    k_alteration_smooth = 7
    # t_max <- 300



# 0. Tweaking function (wrt MEAN) ---------------------------------------

    # REAL DATA
    # confirmed_cases_by_group
    # confirmed_cases_total
    
    # And change to parallel
    
    time_series_to_case_count <- function(time_series) {
      time_series %>% 
        group_by(sim_id, t, policy) %>% 
        # sum up across i_groups
        summarise(across(c(n_new_cases, n_confirmed_cases), sum)) %>% 
        group_by(sim_id) %>% 
        rename(
          n_new_model = n_new_cases,
          n_confirmed_cum_model = n_confirmed_cases
        ) %>% 
        mutate(
          n_confirmed_new_model = n_confirmed_cum_model - lag(n_confirmed_cum_model),
          n_cum_model = cumsum(n_new_model)
        ) %>% 
        ungroup %>% 
        relocate(sim_id, t, policy, n_new_model, n_cum_model, n_confirmed_cum_model, n_confirmed_new_model) %>% 
        mutate(
          across(starts_with("n_"), list(pc = ~ . / n_pop_total)),
        )
    }
    
    
      
      
    tweak_mean <- function(
      sim_time_series, 
      lower_bound_date,
      detection_delay,
      roll_smooth,
      policy_start_times,
      mobility_factors,
      real_cases_total,
      # inc_shift,
      # inc_length
      shift_scale = 1/5,    # conversion from log deviation to change in mobility (e.g. 0.5 log deviation becomes 0.1 change in mobility factor)
      alteration_min = 0.02, # lower bound when stops altering
      k_alteration_smooth = 7
    ) {
      
      # PARAMETERS for testing
      # sim_time_series <- first_sim
      # # t_shift <- -10
      # lower_bound_date <- ymd("2020-06-01")
      # delay <- 21
      # roll_smooth <- 14
      # policy_start_times <- c(0, 96, 135, 176, 218)
      # mobility_factors <- c(1, 1.23, 0.9, 1.3, 1.6)
      # inc_shift <- 0.05
      # inc_length <- 21
      # shift_scale = 1/10
      # real_cases_total <- confirmed_cases_total
      
      # mobility_factors <- mobility_factors_new
      # policy_start_times <- policy_start_times_new
      
      # (1) COMBINE REAL AND MODEL DATA
      
      # REAL DATA
      # Get the rolling average pc in real data
      real_cases_slim <- real_cases_total %>%
        filter(n_confirmed_new != 0) %>% 
        mutate(across(c(n_confirmed_new_pc), list(roll = ~ zoo::rollmean(.x, k = roll_smooth, na.pad = TRUE, align = "center")))) %>% 
        mutate(real_pc_roll = n_confirmed_new_pc_roll)
      
      # What's the value of pc_roll on lower bound date
      lower_bound_date_value <- real_cases_slim %>% filter(date_results == lower_bound_date) %>% .$real_pc_roll
      
      
      # MODEL DATA
      model_cases_total <- time_series_to_case_count(sim_time_series)
    
      model_cases_slim <- model_cases_total %>% 
        mutate(across(c(n_confirmed_new_model_pc), list(roll = ~ zoo::rollmean(.x, k = roll_smooth, na.pad = TRUE, align = "center")))) %>% 
        mutate(model_pc_roll = n_confirmed_new_model_pc_roll) %>% 
        select(sim_id, t, model_pc_roll) %>% 
        group_by(t) %>% 
        quantile_summarise(model_pc_roll, conf_level = 0.9)
      
      # lower bound date is equivalent to which t?
      lower_bound_t <- model_cases_slim %>% filter(model_pc_roll_mean >= lower_bound_date_value) %>% .$t %>% first()
      
      model_data_with_date <- model_cases_slim %>% 
        mutate(date_model = days(t) + lower_bound_date - days(lower_bound_t))
    
      
      # MERGE REAL AND MODEL
      match_model <- real_cases_slim %>% 
        left_join(model_data_with_date, by = c("date_results" = "date_model"))
      
      max_date <- max(match_model$date_results, na.rm = TRUE) - days(detection_delay)
      max_t <- match_model$t[match_model$date_results == max_date & !is.na(match_model$t)]
    
      
      # (2) PLOT RESULTS
      log_plot <- ggplot(match_model, aes(x = date_results)) + 
        geom_line(aes(y = log(real_pc_roll)), size = 1.5, colour = "skyblue") + 
        geom_ribbon(aes(ymin = log(model_pc_roll_lower), ymax = log(model_pc_roll_upper)), fill = "indianred", alpha = 0.2) +
        geom_line(aes(y = log(model_pc_roll_median)), size = 1, colour = "indianred") +
        geom_point(aes(y = log(model_pc_roll_mean))) +
        # geom_vline(xintercept = model_deviations$date_results[model_deviations$t == t_deviate]) + 
        theme_custom() + 
        scale_x_date(breaks = "months", date_labels = "%b %d")
      
      print(log_plot)
      
      # Deviation plot
      # ggplot(match_model, aes(x = date_results)) + 
        # geom_line(aes(y = log(real_pc_roll) - log(model_pc_roll_median)))
      
      # Non log plot
      # plot <- ggplot(match_model, aes(x = date_results)) + 
      #   geom_line(aes(y = real_pc_roll), size = 1.5, colour = "skyblue") + 
      #   geom_ribbon(aes(ymin = model_pc_roll_lower, ymax = model_pc_roll_upper), fill = "indianred", alpha = 0.2) +
      #   geom_line(aes(y = model_pc_roll_median), size = 1.5, colour = "grey") +
      #   geom_point(aes(y = model_pc_roll_median)) +
      #   # geom_vline(xintercept = model_deviations$date_results[model_deviations$t == t_deviate]) + 
      #   theme_custom() + 
      #   scale_x_date(breaks = "months", date_labels = "%b %d")
      # 
      # print(plot)
      
      # (3) CALCULATE WHERE DEVIATIONS OCCUR
      # FIND THE LOG DEVIATIONS FROM THE MEAN
      model_deviations <- match_model %>% 
        filter(date_results >= lower_bound_date) %>% 
        filter(!is.na(real_pc_roll)) %>% 
        ungroup %>% 
        mutate(median_dev = log(real_pc_roll) - log(model_pc_roll_median),
               mean_dev = log(real_pc_roll) - log(model_pc_roll_mean)) %>% 
        mutate(t_alter = t - detection_delay) %>% 
        mutate(alteration = mean_dev * shift_scale,
               alteration = if_else(abs(alteration) < alteration_min, 0, alteration)) %>% 
        mutate(alteration_smooth = zoo::rollmean(alteration, k = k_alteration_smooth, na.pad = TRUE)) %>%
        # Possibly smooth here?
        select(t_alter, alteration, alteration_smooth)
      
      # Count how many correct values we have:
      correct <- model_deviations %>% filter(t_alter <= max_t) %$% sum(alteration == 0)
      total <- model_deviations %>% filter(t_alter <= max_t) %>% nrow()
      prop_correct <- correct/total
      print(str_glue("{correct} / {total} periods matched successfully (= {round(prop_correct, 2)})"))
      
      

  
      # (4) CREATE NEW POLICY VECTOR
      # Convert policy to vector
      diff_times <- lead(policy_start_times) - policy_start_times
      diff_times <- if_else(is.na(diff_times), 500, diff_times)
      # if (length(mobility_factors) > 1) {
      #   mobility_factors[length(mobility_factors)] <- NULL # adjust last value (=1) to be deleted
      # }
      
      policy_vec <- rep(mobility_factors, times = diff_times)
      
    
      

      # Alter the policy vec
      policy_tibble <- tibble(t = 0:(length(policy_vec)-1), 
                              policy_vec = policy_vec) 
      
      
      
      
      # end_policy <- mean(policy_tibble$policy_vec[policy_tibble$t >= max_t - 14], na.rm = TRUE)
      
      policy_tibble_with_alter <- policy_tibble %>% 
        # mutate(policy_vec = if_else(t > max_t, NA_real_, policy_vec)) %>%
        left_join(model_deviations, by = c("t" = "t_alter")) %>% 
        # mutate(alteration_smooth)
        mutate(alteration_smooth = if_else(is.na(alteration_smooth) & t < max(t[!is.na(alteration_smooth)], na.rm = TRUE), 0, alteration_smooth)) %>%
        mutate(policy_vec = policy_vec + alteration_smooth) %>% 
        filter_track(!is.na(policy_vec))

      policy_vec_plot <- tibble(t = 0:(length(policy_tibble_with_alter$policy_vec)-1),
             policy_vec = policy_tibble_with_alter$policy_vec) %>%
        ggplot(aes(x =  t, y = policy_vec)) + geom_line()
      
      # Convert back to start times
      mobility_factors_new <- c(first(policy_tibble_with_alter$policy_vec), policy_tibble_with_alter$policy_vec[policy_tibble_with_alter$policy_vec != lag(policy_tibble_with_alter$policy_vec)][-1])
      policy_start_times_new <- c(0, which(policy_tibble_with_alter$policy_vec != lag(policy_tibble_with_alter$policy_vec)) - 1)
      k_manual_new <- mobility_factors_new %>% 
        map(~ list(k_matrix = k_matrix_basic * .x)) %>% 
        set_names(as.character(1:length(.)))
      
      # Output results
      out <- list(
        mobility_factors = mobility_factors_new,
        k_mobility = k_manual_new,
        policy_start_times = policy_start_times_new,
        # lower_bound_date = lower_bound_date_new,
        lower_bound_t = lower_bound_t,
        prop_correct = prop_correct,
        plot = log_plot,
        policy_vec_plot = policy_vec_plot
      )
      
      return(out)
      
    }
    
    
    

    
    # tweak_median(
    #   first_sim,
    #   lower_bound_date = ymd("2020-06-01"),
    #   detection_delay = 21,
    #   roll_smooth = 14,
    #   policy_start_times = 0,
    #   mobility_factors = 1,
      # shift_scale = 1/5,    # conversion from log deviation to change in mobility (e.g. 0.5 log deviation becomes 0.1 change in mobility factor)
      # alteration_min = 0.02, # lower bound when stops altering
      # k_alteration_smooth = 7
    # )
    
    
    

# 3. First sim  ----------------------------------------------------------------

# SETUP model -------------------------------------------------------------
    

    
    first_sim <- outbreak_sims_policy_change_parallel(
      n_sims = n_sims_per_round,
      n_cores = n_cores,
      n_iterations = t_max,
      keep_all_data = FALSE,
      print_detail = TRUE,
      constant_params = list(
        n_pop = n_pop_total, 
        hh_size_data = data_save$hh_data_bogota,
        n_initial_cases = n_initial_cases, # calculated within function if remains NULL
        group_props = group_props,     # ditto
        dt_approx = 1,
        recov_val = 10, # people recover 10 days after symptoms show (cuts off only ~ 1% of infections at the upper end of timing distribution)
        params_timing = data_save$params_timing,
        params_symptom_timing = data_save$params_symptom_timing,
        params_serial = data_save$params_serial,
        test_delay_data = data_save$test_delay_data,
        test_sensitivity_data = data_save$test_sensitivity_data, 
        ct_delay_data = data_save$ct_delay_data,
        probs_self_test_df = data_save$probs_self_test,
        probs_isolate_symptoms_df = data_save$probs_isolate_symptoms,
        probs_isolate_test_df = data_save$probs_isolate_test,
        probs_isolate_ct_df =   data_save$probs_isolate_ct,
        p_contact_if_isolated_home = rep(1, 4), 
        p_contact_traced = data_save$p_contact_traced, # TAKEN FROM ABOVE
        p_hh_quarantine = data_save$p_hh_quarantine,
        
        contact_dispersion = data_save$contact_dispersion,
        sar_out = data_save$sar_out_input,
        sar_home = data_save$sar_home_input, 
        infectiousness_by_symptoms = data_save$infectiouness_by_symptoms,
        alpha = c(0, 0, 0, 0)
      ),
      policy_start_times = policy_start_times_0,
      time_varying_params = k_mobility_0
    )
    
    updated_params <-  tweak_mean(
      first_sim,
      real_cases_total = confirmed_cases_total,
      policy_start_times = policy_start_times_0,
      mobility_factors = mobility_factors_0,
      detection_delay = detection_delay,
      roll_smooth = roll_smooth,
      lower_bound_date = lower_bound_date,
      shift_scale = shift_scale,    # conversion from log deviation to change in mobility (e.g. 0.5 log deviation becomes 0.1 change in mobility factor)
      alteration_min = alteration_min, # lower bound when stops altering
      k_alteration_smooth = k_alteration_smooth
    )
      
    print(updated_params$plot)
    ggsave(paste0("figures/tweaking_process/tweak_", 0, ".png"), width = 5, height = 5)
    
    # print(updated_params$policy_vec_plot)
      
    # lower_bound_date <- updated_params$lower_bound_date
    policy_start_times <-  updated_params$policy_start_times
    mobility_factors <-  updated_params$mobility_factors
    k_mobility <- updated_params$k_mobility
    
    print("Mobility factors:")
    print(mobility_factors)
    print("Policy start times:")
    print(policy_start_times)
    
    # RECORD
    tweaks <- list(
      "0" = list(
        # lower_bound_date = lower_bound_date,
        policy_start_times = policy_start_times,
        mobility_factors = mobility_factors
        # plot = updated_params$plot
      )
    )
    
    

# 4. Loop (1)  -----------------------------------------------------------------

    
    
    RNGkind("L'Ecuyer-CMRG"); set.seed(12345)
    
    for (i in 1:n_loops) {
      
      print(str_glue("i = {i}"))
      
      later_sim <- outbreak_sims_policy_change_parallel(
        n_sims = n_sims_per_round,
        n_cores = n_cores,
        n_iterations = t_max,
        keep_all_data = FALSE,
        print_detail = FALSE,
        constant_params = list(
          n_pop = n_pop_total, 
          hh_size_data = data_save$hh_data_bogota,
          n_initial_cases = n_initial_cases, # calculated within function if remains NULL
          group_props = group_props,     # ditto
          dt_approx = 1,
          recov_val = 10, # people recover 10 days after symptoms show (cuts off only ~ 1% of infections at the upper end of timing distribution)
          params_timing = data_save$params_timing,
          params_symptom_timing = data_save$params_symptom_timing,
          params_serial = data_save$params_serial,
          test_delay_data = data_save$test_delay_data,
          test_sensitivity_data = data_save$test_sensitivity_data, 
          ct_delay_data = data_save$ct_delay_data,
          probs_self_test_df = data_save$probs_self_test,
          probs_isolate_symptoms_df = data_save$probs_isolate_symptoms,
          probs_isolate_test_df = data_save$probs_isolate_test,
          probs_isolate_ct_df =   data_save$probs_isolate_ct,
          p_contact_if_isolated_home = rep(1, 4), 
          p_contact_traced = data_save$p_contact_traced, # TAKEN FROM ABOVE
          p_hh_quarantine = data_save$p_hh_quarantine,
          
          contact_dispersion = data_save$contact_dispersion,
          sar_out = data_save$sar_out_input,
          sar_home = data_save$sar_home_input, 
          infectiousness_by_symptoms = data_save$infectiouness_by_symptoms,
          alpha = c(0, 0, 0, 0)
        ),
        policy_start_times = policy_start_times,
        time_varying_params = k_mobility
      )
      
      
      updated_params <-  tweak_mean(
        later_sim,
        real_cases_total = confirmed_cases_total,
        policy_start_times = policy_start_times,
        mobility_factors = mobility_factors,
        lower_bound_date = lower_bound_date,
        detection_delay = detection_delay,
        roll_smooth = roll_smooth,
        shift_scale = shift_scale,    # conversion from log deviation to change in mobility (e.g. 0.5 log deviation becomes 0.1 change in mobility factor)
        alteration_min = alteration_min, # lower bound when stops altering
        k_alteration_smooth = k_alteration_smooth
      )
      
      print(updated_params$plot)
      ggsave(paste0("figures/tweaking_process/tweak_", i, ".png"), width = 5, height = 5)
      
      # print(updated_params$policy_vec_plot)
      
      # UPDATE PARAMETERS for next round
      policy_start_times <-  updated_params$policy_start_times
      mobility_factors <-  updated_params$mobility_factors
      k_mobility <- updated_params$k_mobility  
      
      # RECORD
      tweaks[[as.character(i)]] <- list(
        # lower_bound_date = lower_bound_date,
        policy_start_times = policy_start_times,
        mobility_factors = mobility_factors
        # plot = updated_params$plot
      )
      
      tweaks_df <- tweaks %>% map(enframe) %>% map(pivot_wider) %>% 
        bind_rows(.id = "iteration")
      
      
      tweaks_df_slim <- tweaks_df
      save(tweaks_df_slim, file = "data/processed/tweaks_mean.RData")
      
    }
    
    
    
    

# # 5. Loop (2) -------------------------------------------------------------
# 
#     # LOAD PARAMETERS FROM SAVE
#     load("data/processed/tweaks_mean.RData")
#     policy_start_times <- tweaks_df_slim %>% .$policy_start_times %>% last()
#     mobility_factors <- tweaks_df_slim$mobility_factors %>% last()
#     k_mobility <- mobility_factors %>%
#       map(~ list(k_matrix = k_matrix_basic * .x)) %>%
#       set_names(as.character(1:length(.)))
#     
#     n_loops_2 <- 5
#     
#     RNGkind("L'Ecuyer-CMRG"); set.seed(12345)
# 
#     for (i in (1:n_loops_2 + n_loops)) {
# 
#       print(str_glue("i = {i}"))
# 
#       later_sim <- outbreak_sims_policy_change_parallel(
#         n_sims = n_sims_per_round,
#         n_cores = n_cores,
#         n_iterations = t_max,
#         keep_all_data = FALSE,
#         print_detail = FALSE,
#         constant_params = list(
#           n_pop = n_pop_total,
#           hh_size_data = data_save$hh_data_bogota,
#           n_initial_cases = n_initial_cases, # calculated within function if remains NULL
#           group_props = group_props,     # ditto
#           dt_approx = 1,
#           recov_val = 10, # people recover 10 days after symptoms show (cuts off only ~ 1% of infections at the upper end of timing distribution)
#           params_timing = data_save$params_timing,
#           params_symptom_timing = data_save$params_symptom_timing,
#           params_serial = data_save$params_serial,
#           test_delay_data = data_save$test_delay_data,
#           test_sensitivity_data = data_save$test_sensitivity_data,
#           ct_delay_data = data_save$ct_delay_data,
#           probs_self_test_df = data_save$probs_self_test,
#           probs_isolate_symptoms_df = data_save$probs_isolate_symptoms,
#           probs_isolate_test_df = data_save$probs_isolate_test,
#           probs_isolate_ct_df =   data_save$probs_isolate_ct,
#           p_contact_if_isolated_home = rep(1, 4),
#           p_contact_traced = data_save$p_contact_traced, # TAKEN FROM ABOVE
#           p_hh_quarantine = data_save$params_data$p_hh_quarantine,
# 
#           contact_dispersion = data_save$contact_dispersion,
#           sar_out = data_save$params_data$sar_out_input,
#           sar_home = data_save$params_data$sar_home_input,
#           infectiousness_by_symptoms = data_save$infectiouness_by_symptoms,
#           alpha = c(0, 0, 0, 0)
#         ),
#         policy_start_times = policy_start_times,
#         time_varying_params = k_mobility
#       )
# 
# 
#       updated_params <-  tweak_mean(
#         later_sim,
#         real_cases_total = confirmed_cases_total,
#         policy_start_times = policy_start_times,
#         mobility_factors = mobility_factors,
#         lower_bound_date = lower_bound_date,
#         detection_delay = detection_delay,
#         roll_smooth = roll_smooth,
#         shift_scale = shift_scale,    # conversion from log deviation to change in mobility (e.g. 0.5 log deviation becomes 0.1 change in mobility factor)
#         alteration_min = alteration_min, # lower bound when stops altering
#         k_alteration_smooth = k_alteration_smooth
#       )
# 
#       print(updated_params$plot)
#       ggsave(paste0("figures/tweaking_process/tweak_", i, ".png"), width = 5, height = 5)
# 
#       # UPDATE PARAMETERS for next round
#       policy_start_times <-  updated_params$policy_start_times
#       mobility_factors <-  updated_params$mobility_factors
#       k_mobility <- updated_params$k_mobility
# 
#       # RECORD
#       tweaks[[as.character(i)]] <- list(
#         # lower_bound_date = lower_bound_date,
#         policy_start_times = policy_start_times,
#         mobility_factors = mobility_factors
#         # plot = updated_params$plot
#       )
# 
#       tweaks_df <- tweaks %>% map(enframe) %>% map(pivot_wider) %>%
#         bind_rows(.id = "iteration")
# 
# 
#       tweaks_df_slim <- tweaks_df
#       save(tweaks_df_slim, file = "data/processed/tweaks_mean.RData")
# 
#     }
    
    
    
    save(later_sim, file = "data/processed/tweak_final_sim.RData")
    
    
    

    
