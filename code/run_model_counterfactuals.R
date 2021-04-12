

# Run baseline in detail --------------------------------------------------

    # n_pop <- 100000
    # 
    # set.seed(12345)
    # model_baseline <- group_diff_default(n_pop = n_pop, n_sims = 1, print_detail = FALSE)
    # 
    # model_baseline_time_series <- model_baseline$time_series
    # 
    # # save(model_baseline_time_series, file = "data/processed/model_baseline_final.RData")
    # 
    # 
    # 
    # save(model_baseline_time_series, file = "data/processed/model_baseline_final_v2.RData")
    
    

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
    # 
    # model_baseline %>% graph_prevalence()
    # model_baseline %>% graph_prop_infected()
    # 
    # model_baseline$outbreak_t[[1]]$params$k_matrix
    
    # model_baseline_time_series <- model_basline_time_series
    


    
    
    
    
    

# .............................. ------------------------------------------




      
    
    
    

# ....... -----------------------------------------------------------------



# COUNTERFACTUALS ---------------------------------------------------------


# Run counterfactuals using function to equalise params --------------------


  
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
      "home_contacts"        = c("hh_size_data", "sar_home"),
      "isolation_behaviour"  = c("probs_isolate_symptoms_df", "probs_isolate_test_df", "probs_isolate_ct_df", "p_hh_quarantine"),
      "testing_tracing"      = c("probs_self_test_df", "test_choice_delay_data", "test_results_delay_data", "p_contact_traced")
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
    ) %>% 
      select(-starts_with("days_work"), -days_of_work)


    params_combined_df %>% head(2) %>% glimpse()

    # RUN THE MODELS
    RNGkind("L'Ecuyer-CMRG"); set.seed(12345)
    group_diffs <- params_combined_df %>%
      mutate(model_output = pmap(select(., -parameter_set, -params_to_equalise),
                                 group_diff,
                                 print_detail = FALSE,
                                 n_cores = n_cores,
                                 n_sims = n_sims,
                                 n_iterations = n_iterations))

    # system("say R has finished running")


    # Slim the data to just time series to be used in graphs
    group_diff_slim <- group_diffs %>%
      select(parameter_set, any_of("params_to_equalise"), model_output) %>%
      # mutate(time_series = map(model_output, "time_series")) %>%
      # select(-model_output) %>%
      rename(time_series = model_output) %>%
      unnest(time_series) %>%
      mutate(parameter_set = factor(parameter_set, levels = names(params_to_equalise_list)))



    save(group_diff_slim, file = "data/processed/group_diff_slim.RData")
    # save(group_diff_slim, file = "data/temp/group_diff_mobility_manual_v2.RData")
    # save(group_diff_slim, file = "data/temp/group_diff_mobility_manual.RData")
    # save(group_diff_slim, file = "data/processed/group_diff_slim_v2.RData")
    # load("data/processed/group_diff_slim.RData")
    
    
    
    
    


    
    
    

