


# .................... ---------------------------------------------------------------------
# MAIN RESULTS - GROUP DIFFS - 5+ groups --------------------------------------------------
  


# Group diff default function ---------------------------------------------



    group_diff_default <- function(n_sims = 1,
                                   n_iterations = 1000,
                                   keep_all_data = FALSE,
                                   print_detail = TRUE,
                                   n_pop = 10000, 
                                   hh_size_data = data_save$hh_data_bogota,
                                   n_initial_cases = NULL, # calculated within function if remains NULL
                                   group_props = NULL,     # ditto
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
                                   p_hh_quarantine = data_save$params_data$p_hh_quarantine,
                                   k_matrix = data_save$k_matrix * data_save$k_scale_factor, # multiply by new n_pop here,
                                   contact_dispersion = data_save$contact_dispersion,
                                   # beta_matrix = beta_matrix_6,
                                   # r0_group_1 = 5, 
                                   # contacts_mean = contacts_mean_6,
                                   p_immune = c(0, 0, 0, 0),
                                   sar_out = data_save$params_data$sar_out_input,
                                   sar_home = data_save$params_data$sar_home_input, 
                                   infectiousness_by_symptoms =  data_save$infectiouness_by_symptoms,
                                   alpha = c(0, 0, 0, 0)) {
      
      # Calculate number of initial cases and group proportions
      if (is.null(n_initial_cases)) {
        n_initial_cases <- round(data_save$group_props * n_pop / 5000) # 1000th of the pop is infected at the model start
        
        if (sum(n_initial_cases == 0) > 0) stop("n_initial_cases contains 0s, this causes bugs")
      }
      
      
      
      group_props <- c(
        round(data_save$group_props[1:3] * n_pop) / n_pop,
        1 - sum(round(data_save$group_props[1:3] * n_pop) / n_pop)
      )
      
      
      n_pop_scale <- n_pop / data_save$k_matrix_pop
      
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
      
      outbreak_sims(keep_all_data = keep_all_data,
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
    }
    
    # test_group_diff <- group_diff_default(n_iterations = 30, n_initial_cases = c(10, 10, 10, 10))
    # test_group_diff <- group_diff_default(n_iterations = 100, n_initial_cases = c(10, 10, 10, 10), p_immune = c(0.5, 0.5, 0.5, 0.5))
  
    # test_group_diff$params$next_gen_matrix
    # test_group_diff$params$contact_means

    
    
    

# Run baseline in detail --------------------------------------------------

    n_pop <- 100000
    
    set.seed(12345)
    model_baseline <- group_diff_default(n_pop = n_pop, n_sims = 1, print_detail = FALSE)
    
    model_baseline_time_series <- model_baseline$time_series
    
    # save(model_baseline_time_series, file = "data/processed/model_baseline_final.RData")
    
    

    save(model_baseline_time_series, file = "data/processed/model_baseline_final_v2.RData")
    
    

    # save(model_baseline_time_series, file = "model_baseline_100kpop.RData")
    # save(model_baseline_time_series, file = "model_baseline_100kpop_2xtest.RData")
    # save(model_baseline, file = "model_baseline_7.6pop.RData")
    # save(model_baseline_time_series, file = "model_baseline_1m.RData")
           
    # load("OLD/model_baseline_100kpop.RData")
    
    # load("OLD/model_baseline_1m.RData", verbose = TRUE)
  
    # sum(best_guess_initial_cases$n)
    # # 1000 people = 0.
    # format(500 / 7600000, scientific=F)   # 0.00006578947  of the population
    # 
    # 
    # round(best_guess_initial_cases$n * (1e5L / 7.6e6))
    
    # 1000 people actually detected
    
    model_baseline %>% graph_prevalence()
    model_baseline %>% graph_prop_infected()
    
    model_baseline$outbreak_t[[1]]$params$k_matrix
    
    # model_baseline_time_series <- model_basline_time_series
    


    
    
    
    
    

# .............................. ------------------------------------------



# MOBILITY CHANGES with SPLIT -----------------------------------------------------


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
                                           # probs_self_test_df = data_save$probs_self_test,
                                           
                                           probs_isolate_symptoms_df = data_save$probs_isolate_symptoms,
                                           probs_isolate_test_df = data_save$probs_isolate_test,
                                           probs_isolate_ct_df =   data_save$probs_isolate_ct,
                                           p_contact_if_isolated_home = rep(1, 4), 
                                           p_immune = c(0, 0, 0, 0),
                                           p_contact_traced = data_save$p_contact_traced,
                                           p_hh_quarantine = data_save$params_data$p_hh_quarantine,
                                           k_matrix = data_save$k_matrix, # multiply by new n_pop here,
                                           # k_scale_factor = data_save$k_scale_factor,
                                           k_matrix_pop = data_save$k_matrix_pop,
                                           contact_dispersion = data_save$contact_dispersion,
                                           start_date = data_save$start_date,
                                           # lockdown_end_date = data_save$lockdown_end_date,
                                           # lockdown_end_mobility_increase = data_save$lockdown_end_mobility_increase,
                                           # beta_matrix = beta_matrix_6,
                                           # r0_group_1 = 5, 
                                           # contacts_mean = contacts_mean_6,
                                           sar_out = data_save$params_data$sar_out_input,
                                           sar_home = data_save$params_data$sar_home_input, 
                                           infectiousness_by_symptoms = data_save$infectiouness_by_symptoms,
                                           alpha = c(0, 0, 0, 0),
                                           policy_start = data_save$policy_start_times,
                                           mobility_factors = data_save$mobility_factors,
                                           t_param_split,
                                           probs_self_test_df_1 = data_save$probs_self_test,
                                           probs_self_test_df_2 = data_save$probs_self_test %>% mutate(prob = prob * rep(c(1.2, 1.3, 1.3, 1.4), each = 2))
                                           # probs_self_test_df_split = list(
                                           #   data_save$probs_self_test,
                                           #   data_save$probs_self_test %>% mutate(prob = prob * rep(c(1.2, 1.3, 1.3, 1.4), each = 2))
                                           # )
    ) {
      
      
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
      
      # Generate split parameters into a time_varying_params type form
      # t_param_split <- 100
      # probs_self_test_df_split
      
      time_varying_params <- k_change
      time_varying_params[which(policy_start < t_param_split)] <- time_varying_params[which(policy_start < t_param_split)] %>% map(~ append(.x, list(probs_self_test_df = probs_self_test_df_1)))
      time_varying_params[ which(policy_start >= t_param_split)] <- time_varying_params[ which(policy_start >= t_param_split)] %>% map(~ append(.x, list(probs_self_test_df = probs_self_test_df_2)))
      
      
      
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
        # probs_self_test_df = probs_self_test_df, 
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
          # probs_self_test_df = probs_self_test_df, 
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
        time_varying_params = time_varying_params
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
      
      
    
    
    

# ....... -----------------------------------------------------------------



# COUNTERFACTUALS ---------------------------------------------------------

    
    

# Run counterfactuals using function to equalise params --------------------

    # BASELINE PARAMETERS
    params_baseline <- list(
      group_props = data_save$group_props,
      n_pop = 100000,
      
      ct_delay_data             = data_save$ct_delay_data,
      probs_self_test_df_1        = data_save$probs_self_test,
      probs_self_test_df_2        = data_save$probs_self_test %>% mutate(prob = prob * rep(c(1.2, 1.3, 1.3, 1.4), each = 2)),
      probs_isolate_symptoms_df = data_save$probs_isolate_symptoms,
      probs_isolate_test_df     = data_save$probs_isolate_test,
      probs_isolate_ct_df       = data_save$probs_isolate_ct,
      p_contact_if_isolated_home= rep(1, 4),
      p_contact_traced          = data_save$p_contact_traced,
      p_hh_quarantine           = data_save$params_data$p_hh_quarantine,
      hh_size_data              = data_save$hh_data_bogota,
      test_delay_data           = data_save$test_delay_data,
      k_matrix                  = data_save$k_matrix,
      sar_out                   = data_save$params_data$sar_out_input,
      sar_home                  = data_save$params_data$sar_home_input,
      alpha                     = rep(0, 4)
    )
    
    # EQUALISED PARAMETERS
    params_equalised <- equalise_params_to(params_baseline, equalise_to = 4)
    
    # Calculate new list of equalised parameters
    # params_to_equalise_list <- list(
    #   "baseline"             = c(),
    #   "out_contacts"         = c("k_matrix", "sar_out"),
    #   "home_contacts"        = c("p_contact_if_isolated_home", "sar_home", "hh_size_data"),
    #   "isolation_behaviour"  = c("probs_isolate_symptoms_df", "probs_isolate_test_df", "probs_isolate_ct_df", "p_hh_quarantine"),
    #   "testing_tracing"      = c("probs_self_test_df", "test_choice_delay_data", "test_results_delay_data", "ct_delay_data", "p_contact_traced")
    # )
    params_to_equalise_list <- list(
      "baseline"             = c(),
      "out_contacts"         = c("k_matrix"),
      "home_contacts"        = c("hh_size_data"),
      "isolation_behaviour"  = c("probs_isolate_symptoms_df", "probs_isolate_test_df", "probs_isolate_ct_df", "p_hh_quarantine"),
      "testing_tracing"      = c("probs_self_test_df_1", "probs_self_test_df_2", "test_choice_delay_data", "test_results_delay_data", "p_contact_traced")
    )
    
    params_to_equalise_list$all <- flatten_chr(params_to_equalise_list)
    
    # Combine params from old and new
    params_combined <- tibble(
      parameter_set = names(params_to_equalise_list),
      params_to_equalise = params_to_equalise_list
    ) %>% 
      mutate(
        baseline_to_use = map(params_to_equalise, ~ params_baseline[! names(params_baseline) %in% .x]),
        equalised_to_use = map(params_to_equalise, ~ params_equalised[names(params_equalised) %in% .x]),
        params_to_use = map2(baseline_to_use, equalised_to_use, ~ c(.x, .y))
      ) %>% 
      select(-baseline_to_use, equalised_to_use)
    
    
    # Convert to tibble to be pmapped
    params_combined_df <- bind_cols(
      tibble(parameter_set = params_combined$parameter_set,
             params_to_equalise = params_combined$params_to_equalise),
      params_combined$params_to_use %>% purrr::transpose() %>% as_tibble()
    )
    
    
    # params_combined_df %>% head(2) %>% glimpse()
    
    # RUN THE MODELS
    RNGkind("L'Ecuyer-CMRG"); set.seed(12345)
    group_diffs <- params_combined_df %>% 
      mutate(model_output = pmap(select(., -parameter_set, -params_to_equalise),
                                 group_diff_mobility_change, 
                                 print_detail = FALSE,
                                 n_cores = 4,
                                 n_sims = 8, 
                                 n_iterations = 100000))
    
    system("say R has finished running")
    
  
    # Slim the data to just time series to be used in graphs
    group_diff_slim <- group_diffs %>% 
      select(parameter_set, any_of("params_to_equalise"), model_output) %>% 
      mutate(time_series = map(model_output, "time_series")) %>% 
      select(-model_output) %>% 
      unnest(time_series) %>% 
      mutate(parameter_set = factor(parameter_set, levels = names(params_to_equalise_list)))


    
    save(group_diff_slim, file = "data/temp/group_diff_mobility_8thFeb.RData")
    # save(group_diff_slim, file = "data/temp/group_diff_mobility_manual_v2.RData")
    # save(group_diff_slim, file = "data/temp/group_diff_mobility_manual.RData")
    # save(group_diff_slim, file = "data/processed/group_diff_slim_v2.RData")
    # load("data/processed/group_diff_slim.RData")
    
    
    
    
    

