    load("data/processed/tweaks_mean.RData", verbose = TRUE)
    
    
    
    final_policy_start_times <- tweaks_df_slim$policy_start_times[tweaks_df_slim$iteration == tweak_to_choose][[1]]
    final_mobility_factors <- tweaks_df_slim$mobility_factors[tweaks_df_slim$iteration == tweak_to_choose][[1]]
    
    # STRIP OFF THE END that drops due to no data
    # end_adjustment <- which(final_mobility_factors == max(final_mobility_factors))
    # 
    # final_policy_start_times_trimmed <- final_policy_start_times[1:end_adjustment]
    # final_mobility_factors_trimmed <- final_mobility_factors[1:end_adjustment]
    
    data_save$policy_start_times <- final_policy_start_times
    data_save$mobility_factors <- final_mobility_factors
    
    k_matrix_basic <- data_save$k_matrix * (n_pop_total / data_save$k_matrix_pop)

    data_save$k_mobility <- data_save$mobility_factors %>% 
      map(~ list(k_matrix = k_matrix_basic * .x)) %>% 
      set_names(as.character(1:length(.)))
    
    
    
    

# RUN A FINAL SIM TO CHECK ------------------------------------------------

    RNGkind("L'Ecuyer-CMRG"); set.seed(12345)
    
    
    # n_sims_per_round <- 1
    
    final_sim <- outbreak_sims_policy_change_parallel(
      n_sims = n_sims_per_round,
      n_cores = n_cores,
      n_iterations = n_iterations,
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
      policy_start_times = data_save$policy_start_times,
      time_varying_params = data_save$k_mobility
    )
    
    save(final_sim, file = "data/processed/tweak_final_sim_v2.RData")
    save(data_save, file = "data/processed/data_save.RData")
    