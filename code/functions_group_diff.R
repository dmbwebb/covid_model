


# .................... ---------------------------------------------------------------------
# MAIN RESULTS - GROUP DIFFS - 5+ groups --------------------------------------------------



# Group diff default function ---------------------------------------------



group_diff <- function(n_sims = 1,
                       n_iterations = 1000,
                       n_cores,
                       keep_all_data = FALSE,
                       print_detail = TRUE,
                       n_pop = 100000, 
                       hh_size_data = data_save$hh_data_bogota,
                       n_initial_cases = NULL, # calculated within function if remains NULL
                       group_props = data_save$group_props,     # ditto
                       dt_approx = 1,
                       recov_val = 10, # people recover 10 days after symptoms show (cuts off only ~ 1% of infections at the upper end of timing distribution)
                       params_timing = data_save$params_timing,
                       params_symptom_timing = data_save$params_symptom_timing,
                       params_serial = data_save$params_serial,
                       test_delay_data = data_save$test_delay_data,
                       # test_choice_delay_data = data_save$test_choice_delay_data,
                       # test_results_delay_data = data_save$test_results_delay_data, 
                       test_sensitivity_data = data_save$test_sensitivity_data, 
                       ct_delay_data = data_save$ct_delay_data,
                       probs_self_test_df = data_save$probs_self_test,
                       probs_isolate_symptoms_df = data_save$probs_isolate_symptoms,
                       probs_isolate_test_df = data_save$probs_isolate_test,
                       probs_isolate_ct_df =   data_save$probs_isolate_ct,
                       p_contact_if_isolated_home = rep(1, 4), 
                       p_contact_traced = data_save$p_contact_traced,
                       p_hh_quarantine = data_save$p_hh_quarantine,
                       k_matrix = data_save$k_matrix * data_save$k_scale_factor, # multiply by new n_pop here,
                       k_matrix_pop = data_save$k_matrix_pop,
                       contact_dispersion = data_save$contact_dispersion,
                       # beta_matrix = beta_matrix_6,
                       # r0_group_1 = 5, 
                       # contacts_mean = contacts_mean_6,
                       p_immune = c(0, 0, 0, 0),
                       sar_out = data_save$sar_out_input,
                       sar_home = data_save$sar_home_input, 
                       infectiousness_by_symptoms =  data_save$infectiouness_by_symptoms,
                       alpha = c(0, 0, 0, 0)) {
  
  # Calculate number of initial cases and group proportions
  if (is.null(n_initial_cases)) {
    n_initial_cases <- round(group_props * n_pop / 5000) # one 5000th of the pop is infected at the model start
    
    if (sum(n_initial_cases == 0) > 0) stop("n_initial_cases contains 0s, this causes bugs")
  }
  
  
  
  group_props <- c(
    round(group_props[1:3] * n_pop) / n_pop,
    1 - sum(round(group_props[1:3] * n_pop) / n_pop)
  )
  
  
  n_pop_scale <- n_pop / k_matrix_pop
  
  params <- list(
    keep_all_data = keep_all_data,
    print_detail = print_detail,
    n_sims = n_sims, 
    n_iterations = n_iterations,
    n_pop = n_pop, 
    n_initial_cases = n_initial_cases,
    group_props = group_props,
    dt_approx = dt_approx,
    recov_val = recov_val,
    params_timing = params_timing,
    params_symptom_timing = params_symptom_timing,
    params_serial = params_serial,
    test_delay_data = test_delay_data,
    # test_choice_delay_data = test_choice_delay_data,
    # test_results_delay_data = test_results_delay_data,
    test_sensitivity_data = test_sensitivity_data,
    ct_delay_data = ct_delay_data,
    probs_self_test_df = probs_self_test_df, 
    probs_isolate_symptoms_df = probs_isolate_symptoms_df, 
    probs_isolate_test_df = probs_isolate_test_df, 
    probs_isolate_ct_df = probs_isolate_ct_df,
    p_contact_if_isolated_home = p_contact_if_isolated_home, p_contact_traced = p_contact_traced, 
    p_hh_quarantine = p_hh_quarantine,
    k_matrix = k_matrix * n_pop_scale,
    contact_dispersion = contact_dispersion,
    p_immune = p_immune,
    hh_size_data = hh_size_data, 
    alpha = alpha, 
    sar_home = sar_home, sar_out = sar_out,
    infectiousness_by_symptoms = infectiousness_by_symptoms
  )
  
  null_params <- params %>% map_lgl(~ is.null(.x)) %>% any()
  
  if (null_params) {
    # print(detect(params, ~ is.null(.x)))
    print(params)
    stop("one of the params is null")
  }
  
  if (n_cores == 1) {
    out <- outbreak_sims(keep_all_data = keep_all_data,
                         print_detail = print_detail,
                         n_sims = n_sims, 
                         n_iterations = n_iterations,
                         n_pop = n_pop, 
                         n_initial_cases = n_initial_cases,
                         group_props = group_props,
                         dt_approx = dt_approx,
                         recov_val = recov_val,
                         params_timing = params_timing,
                         params_symptom_timing = params_symptom_timing,
                         params_serial = params_serial,
                         test_delay_data = test_delay_data,
                         # test_choice_delay_data = test_choice_delay_data,
                         # test_results_delay_data = test_results_delay_data,
                         test_sensitivity_data = test_sensitivity_data,
                         ct_delay_data = ct_delay_data,
                         probs_self_test_df = probs_self_test_df, 
                         probs_isolate_symptoms_df = probs_isolate_symptoms_df, 
                         probs_isolate_test_df = probs_isolate_test_df, 
                         probs_isolate_ct_df = probs_isolate_ct_df,
                         p_contact_if_isolated_home = p_contact_if_isolated_home, p_contact_traced = p_contact_traced, 
                         p_hh_quarantine = p_hh_quarantine,
                         k_matrix = k_matrix * n_pop_scale,
                         contact_dispersion = contact_dispersion,
                         hh_size_data = hh_size_data, 
                         p_immune = p_immune,
                         alpha = alpha, 
                         sar_home = sar_home, sar_out = sar_out,
                         infectiousness_by_symptoms = infectiousness_by_symptoms
    )
  } else if (n_cores > 1)
  
  out <- outbreak_sims_parallel(keep_all_data = keep_all_data,
                                n_cores = n_cores,
                                print_detail = print_detail,
                                n_sims = n_sims, 
                                n_iterations = n_iterations,
                                n_pop = n_pop, 
                                n_initial_cases = n_initial_cases,
                                group_props = group_props,
                                dt_approx = dt_approx,
                                recov_val = recov_val,
                                params_timing = params_timing,
                                params_symptom_timing = params_symptom_timing,
                                params_serial = params_serial,
                                test_delay_data = test_delay_data,
                                # test_choice_delay_data = test_choice_delay_data,
                                # test_results_delay_data = test_results_delay_data,
                                test_sensitivity_data = test_sensitivity_data,
                                ct_delay_data = ct_delay_data,
                                probs_self_test_df = probs_self_test_df, 
                                probs_isolate_symptoms_df = probs_isolate_symptoms_df, 
                                probs_isolate_test_df = probs_isolate_test_df, 
                                probs_isolate_ct_df = probs_isolate_ct_df,
                                p_contact_if_isolated_home = p_contact_if_isolated_home, p_contact_traced = p_contact_traced, 
                                p_hh_quarantine = p_hh_quarantine,
                                k_matrix = k_matrix * n_pop_scale,
                                contact_dispersion = contact_dispersion,
                                hh_size_data = hh_size_data, 
                                p_immune = p_immune,
                                alpha = alpha, 
                                sar_home = sar_home, sar_out = sar_out,
                                infectiousness_by_symptoms = infectiousness_by_symptoms
  )
  
  return(out)
}








# MOBILITY CHANGES v3 -----------------------------------------------------


group_diff_mobility_change <- function(n_sims = 4,
                                       n_cores,
                                       n_iterations = 5,
                                       keep_all_data = FALSE,
                                       print_detail = TRUE,
                                       n_pop = 100000,
                                       hh_size_data = data_save$hh_data_bogota,
                                       n_initial_cases = NULL, # calculated within function if remains NULL
                                       group_props = data_save$group_props,     # ditto
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
                                       p_immune = c(0, 0, 0, 0),
                                       p_contact_traced = data_save$p_contact_traced,
                                       p_hh_quarantine = data_save$p_hh_quarantine,
                                       k_matrix = data_save$k_matrix, # multiply by new n_pop here,
                                       # k_scale_factor = data_save$k_scale_factor,
                                       k_matrix_pop = data_save$k_matrix_pop,
                                       contact_dispersion = data_save$contact_dispersion,
                                       start_date = data_save$start_date,
                                       lockdown_end_date = data_save$lockdown_end_date,
                                       lockdown_end_mobility_increase = data_save$lockdown_end_mobility_increase,
                                       # beta_matrix = beta_matrix_6,
                                       # r0_group_1 = 5, 
                                       # contacts_mean = contacts_mean_6,
                                       sar_out = data_save$sar_out_input,
                                       sar_home = data_save$sar_home_input, 
                                       infectiousness_by_symptoms = data_save$infectiouness_by_symptoms,
                                       alpha = c(0, 0, 0, 0),
                                       policy_start = data_save$policy_start_times,
                                       mobility_factors = data_save$mobility_factors) {
  
  # Calculate number of initial cases and group proportions
  if (is.null(n_initial_cases)) {
    n_initial_cases <- round(group_props * n_pop / 5000) # 1000th of the pop is infected at the model start
    
    if (sum(n_initial_cases == 0) > 0) stop("n_initial_cases contains 0s, this causes bugs")
  }
  
  group_props <- c(
    round(group_props[1:3] * n_pop) / n_pop,
    1 - sum(round(group_props[1:3] * n_pop) / n_pop)
  )
  
  
  # Generate K matrix that changes over time based on mobility
  n_pop_scale <- n_pop / k_matrix_pop
  k_matrix_basic <- k_matrix * n_pop_scale
  
  k_change <- mobility_factors %>% 
    map(~ list(k_matrix = k_matrix_basic * .x)) %>% 
    set_names(paste0("phase_", as.character(1:length(.))))
  
  params <- list(
    keep_all_data = keep_all_data,
    print_detail = print_detail,
    n_sims = n_sims, 
    n_iterations = n_iterations,
    n_pop = n_pop, 
    n_initial_cases = n_initial_cases,
    group_props = group_props,
    dt_approx = dt_approx,
    recov_val = recov_val,
    params_timing = params_timing,
    params_symptom_timing = params_symptom_timing,
    params_serial = params_serial,
    test_delay_data = test_delay_data,
    # test_choice_delay_data = test_choice_delay_data,
    # test_results_delay_data = test_results_delay_data,
    test_sensitivity_data = test_sensitivity_data,
    ct_delay_data = ct_delay_data,
    probs_self_test_df = probs_self_test_df, 
    probs_isolate_symptoms_df = probs_isolate_symptoms_df, 
    probs_isolate_test_df = probs_isolate_test_df, 
    probs_isolate_ct_df = probs_isolate_ct_df,
    p_contact_if_isolated_home = p_contact_if_isolated_home, p_contact_traced = p_contact_traced, 
    p_hh_quarantine = p_hh_quarantine,
    k_matrix = k_matrix * n_pop_scale,
    contact_dispersion = contact_dispersion,
    hh_size_data = hh_size_data, 
    alpha = alpha, 
    sar_home = sar_home, sar_out = sar_out,
    infectiousness_by_symptoms = infectiousness_by_symptoms
  )
  
  null_params <- params %>% map_lgl(~ is.null(.x)) %>% any()
  
  if (null_params) {
    print(detect(params, ~ is.null(.x)))
    stop("one of the params is null")
  }
  
  out <- outbreak_sims_policy_change_parallel(
    n_sims = n_sims,
    n_cores = n_cores,
    n_iterations = n_iterations,
    keep_all_data = keep_all_data,
    print_detail = print_detail,
    constant_params = list(
      n_pop = n_pop, 
      n_initial_cases = n_initial_cases,
      group_props = group_props,
      dt_approx = dt_approx,
      recov_val = recov_val,
      params_timing = params_timing,
      params_symptom_timing = params_symptom_timing,
      params_serial = params_serial,
      test_delay_data = test_delay_data,
      test_sensitivity_data = test_sensitivity_data,
      ct_delay_data = ct_delay_data,
      probs_self_test_df = probs_self_test_df, 
      probs_isolate_symptoms_df = probs_isolate_symptoms_df, 
      probs_isolate_test_df = probs_isolate_test_df, 
      probs_isolate_ct_df = probs_isolate_ct_df,
      p_contact_if_isolated_home = p_contact_if_isolated_home, 
      p_contact_traced = p_contact_traced, 
      p_hh_quarantine = p_hh_quarantine,
      p_immune = p_immune,
      contact_dispersion = contact_dispersion,
      hh_size_data = hh_size_data, 
      alpha = alpha, 
      sar_home = sar_home, sar_out = sar_out,
      infectiousness_by_symptoms = infectiousness_by_symptoms
    ),
    policy_start_times = policy_start,
    time_varying_params = k_change
  )
  
  
  return(out)
  
}



# TEST
# library("pbmcapply")
# tic()
# test_group_diff <- group_diff_mobility_change(n_cores = 4)
# toc() # new (parallel): 70s, then 69, then 40 [only with time series], then 29 (using 4 cores)
# 
# 
# 
# 
# tic()
# test_group_diff <- group_diff_mobility_change()
# toc() # 39s, then 29s


# How many cores?
