# INEQUALITY IS BAD  ------------------------------------------------------
    
    # params_baseline -- comes from "run_model_main.R"
    # also need group_diff
  
    # CALCUALTE NEW TIGHTENED LIST OF PARAMETERS (takes a while to find rescale factors for the beta matrix, hh_size_data etc.)
    RNGkind("L'Ecuyer-CMRG"); set.seed(12345)
    params_tightened_list <- tibble(
      tighten_factor = c(0, 0.25, 0.5, 0.75, 1),
      # tighten_factor = c(0, 1),
      params_to_scale = list(params_baseline) 
    ) %>% 
      rowwise() %>% 
      mutate(params_tightened = list(tighten_params(params_to_scale, tighten_factor = tighten_factor, all_at_once = TRUE)))
    
    # Transpose so it's amenable to being pmapped 
    params_tightened_transposed <- bind_cols(
      tibble(tighten_factor = params_tightened_list$tighten_factor),
      params_tightened_list$params_tightened %>% 
        purrr::transpose() %>% 
        as_tibble()
    )
    
    
    
    
    rm(params_tightened_list)
    
    # CHECK NEW HH SIZE MEANS
    
    # params_tightened_transposed$contact_means
    
    # params_tightened_transposed$test_results_delay_data %>%
    #   set_names(as.character(c(0, 0.25, 0.5, 0.75, 1))) %>%
    #   bind_rows(.id = "tighten_factor") %>%
    #   filter(test_results_delay <= 10) %>%
    #   ggplot(aes(x = tighten_factor, y = test_results_delay)) +
    #   geom_boxplot()
    
    
    # CHECK ISOLATION PARAMS
    # params_tightened_transposed %>% 
    #   mutate(p_isolate_ct = map(probs_isolate_symptoms_df, ~ filter(.x, symptom_severity == 2))) %>% 
    #   select(tighten_factor, p_isolate_ct) %>% 
    #   unnest(p_isolate_ct) %>% 
    #   ggplot(aes(x = tighten_factor, y = prob, colour = factor(i_group))) + 
    #   geom_line() + geom_point()
  

    
    
    # Run all simulations with new models:
    RNGkind("L'Ecuyer-CMRG"); set.seed(12345)
    inequality_is_bad <- params_tightened_transposed %>% 
      mutate(model_output = pmap(select(., -tighten_factor, -contact_means, -starts_with("days_work"), -days_of_work), 
                                 group_diff, 
                                 # keep_all_data = TRUE, 
                                 print_detail = FALSE,
                                 n_cores = n_cores,
                                 n_sims = n_sims, 
                                 # n_pop = 100000,
                                 n_iterations = n_iterations)) # if you keep all data it's 40GB !! bloody hell
    
    
    # Remove all the parameters (for ease)
    inequality_is_bad_slim <- inequality_is_bad %>% 
      # mutate(time_series = map(model_output, "time_series")) %>% 
      mutate(time_series = model_output) %>%
      select(-model_output) %>%
      select(tighten_factor, time_series) %>% 
      unnest(time_series) 
    
    save(inequality_is_bad_slim, file = "data/processed/inequality_is_bad.RData")
    # save(inequality_is_bad_slim, file = "inequality_is_bad.RData")
    # load("data/processed/inequality_is_bad.RData", verbose = TRUE)
    
   
    
    
    
# Plot parameters to check ------------------------------------------------
    
    # # HH size
    # hh_size_bound <- params_tightened_transposed$hh_size_data %>%
    #   set_names(params_tightened_transposed$tighten_factor) %>%
    #   bind_rows(.id = "tighten_factor")
    # 
    # hh_size_bound %>%
    #   group_by(tighten_factor, i_group) %>%
    #   summarise(mean_se(hh_size)) %>%
    #   ggplot(aes(x = factor(i_group), y = y)) +
    #   geom_point() +
    #   geom_errorbar(aes(ymin = ymin, ymax = ymax), width = 0.2) +
    #   facet_wrap(~ tighten_factor)
    
    
    # named_list_to_tibble <- function(.list) {
    #   tib <- tibble(name = names(.list))
    #   tib$data <- unname(.list)
    #   return(tib)
    # }
    # 
    # params_tightened_long <- params_tightened_list %>% 
    #   select(-params_to_scale) %>% 
    #   ungroup %>% 
    #   mutate(params_tightened = map(params_tightened, named_list_to_tibble)) %>% 
    #   unnest(params_tightened) %>% 
    #   arrange(name)
    # 
    # params_tightened_long_no_data <- params_tightened_long %>% 
    #   filter(!str_detect(name, "data"), map_lgl(data, ~ is_atomic(.x) & !is.matrix(.x))) %>% 
    #   mutate(i_group = list(1:4)) %>% 
    #   unnest(c(data, i_group)) %>% 
    #   filter(name != "alpha") %>% 
    #   print_all
    # 
    # 
    # # PLOT ALL VECTOR BASED
    # params_tightened_long_no_data %>% 
    #   # filter(name == "contact_means") %>% 
    #   
    #   ggplot(aes(x = tighten_factor, y = data, colour = factor(i_group))) + 
    #   geom_line(alpha = 0.3) + geom_point() + 
    #   facet_wrap(~ name, scales = "free_y")
    # 
    # # PLOT ALL DATA BASED:
    # 
    # # (i) ct_delay_data
    # params_tightened_long %>% filter(name == "ct_delay_data") %>% #.$data %>% .[[1]] %>% 
    #   rowwise() %>% 
    #   mutate(
    #     summ_data = list(data %>% group_by(i_group) %>% summarise(mean_se(ct_delay)))
    #   ) %>% 
    #   unnest(summ_data) %>% 
    #   
    #   # Plot
    #   ggplot(aes(x = tighten_factor, y = y, colour = factor(i_group))) + 
    #   geom_line(alpha = 0.3) + geom_point() + 
    #   facet_wrap(~ name, scales = "free_y")
    # 
    # # (ii) test delay data
    # params_tightened_long %>% filter(name == "test_delay_data") %>% #.$data %>% .[[1]] %>% 
    #   rowwise() %>% 
    #   mutate(
    #     summ_data = list(data %>% group_by(i_group) %>% summarise(mean_se(test_results_delay)))
    #   ) %>% 
    #   unnest(summ_data) %>% 
    #   
    #   # Plot
    #   ggplot(aes(x = tighten_factor, y = y, colour = factor(i_group))) + 
    #   geom_line(alpha = 0.3) + geom_point() + 
    #   facet_wrap(~ name, scales = "free_y")
    # 
    # params_tightened_long %>% filter(name == "test_delay_data") %>% #.$data %>% .[[1]] %>% 
    #   rowwise() %>% 
    #   mutate(
    #     summ_data = list(data %>% group_by(i_group) %>% summarise(mean_se(test_choice_delay)))
    #   ) %>% 
    #   unnest(summ_data) %>% 
    #   
    #   # Plot
    #   ggplot(aes(x = tighten_factor, y = y, colour = factor(i_group))) + 
    #   geom_line(alpha = 0.3) + geom_point() + 
    #   facet_wrap(~ name, scales = "free_y")
    # 
    # 
    # # (iii) HH size data
    # params_tightened_long %>% filter(name == "hh_size_data") %>% #.$data %>% .[[1]] %>% 
    #   rowwise() %>% 
    #   mutate(
    #     summ_data = list(data %>% group_by(i_group) %>% summarise(mean_se(hh_size)))
    #   ) %>%  
    #   unnest(summ_data) %>% 
    #   
    #   # Plot
    #   ggplot(aes(x = tighten_factor, y = y, colour = factor(i_group))) + 
    #   geom_line(alpha = 0.3) + geom_point() + 
    #   facet_wrap(~ name, scales = "free_y")
    
    
    
    
# Trouble shooting inequality is bad --------------------------------------
    
    # inequality_is_bad_outbreak <- inequality_is_bad$model_output[[1]]
    # 
    # inequality_is_bad_outbreak$params$contact_weights %>% print_all
    # inequality_is_bad_outbreak$time_series %>% 
    #   mutate(p_infected = n_cases_live / n_pop, .before = n_susceptible) %>% 
    #   print(n = 200)
    # 
    # inequality_is_bad_live <- inequality_is_bad_outbreak$outbreak_t_record[[1]] %>% map("live_cases") %>% 
    #   set_names(as.character(0:50)) %>% 
    #   bind_rows(.id = "t") %>% 
    #   mutate(t = as.numeric(t), i_group = factor(i_group)) %>% 
    #   arrange(case_id, t) %>% 
    #   relocate(case_id, t)
    # 
    # inequality_is_bad_secondary <- inequality_is_bad_outbreak$outbreak_t_record[[1]] %>% map("secondary_cases") %>% 
    #   set_names(as.character(0:50)) %>% 
    #   bind_rows(.id = "t") %>% 
    #   mutate(t = as.numeric(t), i_group = factor(i_group)) %>% 
    #   arrange(case_id, secondary_case_id, t) %>% 
    #   relocate(case_id, secondary_case_id, t)
    # 
    # 
    # # Number of secondary cases 
    # inequality_is_bad_secondary %>% 
    #   select(case_id, secondary_case_id, i_group, secondary_i_group) %>% 
    #   distinct() %>% 
    #   group_by(case_id, i_group) %>% 
    #   summarise(n = n()) %>% 
    #   
    #   ggplot(aes(x = log(n + 1), colour = factor(i_group))) + 
    #   geom_density(adjust = 2)
    # 
    # 
    # # Avg number of seconray cases
    # inequality_is_bad_secondary %>% 
    #   select(case_id, secondary_case_id, i_group, secondary_i_group) %>% 
    #   distinct() %>% 
    #   group_by(i_group, case_id) %>% 
    #   summarise(n = n()) %>% 
    #   summarise(avg_g = mean(n))
    # 
    # # Proportion from each group
    # inequality_is_bad_secondary %>% 
    #   select(case_id, secondary_case_id, i_group, secondary_i_group) %>% 
    #   distinct() %>% 
    #   group_by(i_group, secondary_i_group) %>% 
    #   summarise(n = n()) %>% 
    #   mutate(p = n / sum(n))
    # 
    # 
    # inequality_is_bad_outbreak %>% plot_detected_by_i_group()  # group 6 is much slower to get detected for some reason !!
    # 
    # inequality_is_bad_outbreak %>% plot_by_i_group()
    
    
    
    
    
