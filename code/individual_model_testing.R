    source("individual_model_functions.R")    # imports all functions
    
    # source("import_data.R")
    load("data_save")                         # loads all cleaned data


# TESTING -----------------------------------------------------------------

# Test draw_ functions ----------------------------------------------------

    # (1) SYMPTOMS/RECOVERY
    probs_default <- c(0.3, 0.6, 0.8, 0.9, 1)
    probs_df_basic <- tibble(
      i_group = c(rep(1, 5), rep(2, 5), rep(3, 5)),
      symptom_severity = as.integer(rep(1:5, 3))
    )
    
    k_matrix_trial <- Matrix::forceSymmetric(
      matrix(
        c(4, 1, 1, 1, 1, 1,
          4,  6, 1, 1, 1, 1,
          3,  4, 7, 1, 1, 1,
          1,  3, 3, 4, 1, 1,
          1,  1, 2, 2, 3, 1, 
          0,  0, 1, 1, 2, 3),
        nrow = 6
      )
    ) %>% as.matrix()
    

    
    # Geenrate hh data
    hh_data_example <- gen_hh_data(n_pop = 10000, group_props = c(0.2, 0.5, 0.3), 
                              hh_size_data = list(
                                `1` = rbinom(10000, size = 10, prob = 0.3) + 1,
                                `2` = rbinom(10000, size = 5, prob = 0.5) + 1,
                                `3` = rbinom(10000, size = 5, prob = 0.5) + 1
                              ))
    
    hh_data_sampled <- hh_data_example %>% 
      rowwise() %>% 
      mutate(hh_ind_id = list(1:hh_size)) %>% 
      unnest(hh_ind_id) %>% 
      mutate(hh_ind_id = paste0(hh_id, "_", hh_ind_id)) %>% 
      slice_sample(n = 10000)
    
    
    draw_symptoms_debug <- draw_symptoms_recovery(ids = 1:10000,
                                                  hh_id = hh_data_sampled$hh_id,
                                                  i_group = hh_data_sampled$i_group,
                                                  hh_size = hh_data_sampled$hh_size,
                                                  infection_timings = 0,          
                                                  contact_tracing_test_timing = NULL,
                                                  contact_tracing_results_timing = NULL,
                                                  # random_test_yn = NULL,
                                                  # random_test_timing = NULL,
                                                  # isolate_random_test_timing = NULL,
                                                  test_results_delay_data = c(1.5, 1, 1),
                                                  recov_val = 5,
                                                  test_sensitivity = 0.85,
                                                  probs_self_test_df = probs_df_basic %>% mutate(prob = c(probs_default * 0.1, probs_default * 0.5, probs_default * 1)),
                                                  probs_isolate_ct_df = probs_df_basic %>% mutate(prob = c(probs_default * 0.1, probs_default * 0.5, probs_default * 1)),
                                                  probs_isolate_symptoms_df = probs_df_basic %>% mutate(prob = c(probs_default * 0.1, probs_default * 0.3, probs_default * 0.5)),
                                                  probs_isolate_test_df = probs_df_basic %>% mutate(prob = c(probs_default * 0.2, probs_default * 0.2, probs_default * 1))) %>%
      mutate(hh_ind_id = hh_data_sampled$hh_ind_id)
    
    # draw_symptoms_debug %>% 
      # ggplot(aes(x = isolate_after_ct, fill = factor(i_group))) + geom_bar(position = "dodge") # plot isolation decisions

    
    # (2) RANDOM TESTING
    draw_random_testing_debug <- draw_symptoms_debug %>% 
      draw_random_testing(alpha = c(0.1, 0.2, 0.3), dt_approx = 1, test_sensitivity = 0.5) %>% 
      update_timing(dt = 1) %>% 
      draw_random_testing(alpha = c(0.1, 0.2, 0.3), dt_approx = 1, test_sensitivity = 0.5) %>% 
      update_timing(dt = 1) %>% 
      draw_random_testing(alpha = c(0.1, 0.2, 0.3), dt_approx = 1, test_sensitivity = 0.5) %>% 
      count_nas() %>% 
      count_prop(random_test_yn, is.na(random_test_timing))
    
    draw_random_testing_debug %>% 
      ggplot(aes(x = random_test_yn, y = ..prop.., group = factor(i_group), fill = factor(i_group))) + geom_bar(position = "dodge") # plot isolation decisions
    
    
    
    # (3) TEST ALL THE DRAWS WORK TOGETHER
    
    draws_debug <- draw_symptoms_recovery(ids = 1:10000,
                                          hh_id = hh_data_sampled$hh_id,
                                          i_group = hh_data_sampled$i_group,
                                          hh_size = hh_data_sampled$hh_size,
                                          infection_timings = 0,          
                                          contact_tracing_test_timing = NULL,
                                          contact_tracing_results_timing = NULL,
                                          # random_test_yn = NULL,
                                          # random_test_timing = NULL,
                                          # isolate_random_test_timing = NULL,
                                          test_results_delay_data = c(1.5, 1, 1),
                                          recov_val = 5,
                                          test_sensitivity = 0.85,
                                          probs_self_test_df = probs_df_basic %>% mutate(prob = c(probs_default * 0.1, probs_default * 0.5, probs_default * 1)),
                                          probs_isolate_ct_df = probs_df_basic %>% mutate(prob = c(probs_default * 0.1, probs_default * 0.5, probs_default * 1)),
                                          probs_isolate_symptoms_df = probs_df_basic %>% mutate(prob = c(probs_default * 0.1, probs_default * 0.3, probs_default * 0.5)),
                                          probs_isolate_test_df = probs_df_basic %>% mutate(prob = c(probs_default * 0.2, probs_default * 0.2, probs_default * 1))) %>% 
      # group_by(hh_id) %>% 
      mutate(hh_ind_id = hh_data_sampled$hh_ind_id) %>% 
      # ungroup %>% 
      # mutate(hh_ind_id = paste0(hh_id, "_", sample(1:5, nrow(.), replace = TRUE))) %>% 
      # draw_random_testing(alpha = c(0.1, 0.2, 0.3), dt_approx = 1, test_sensitivity = 0.85) %>% 
      # update_timing(dt = 1) %>% 
      # draw_random_testing(alpha = c(0.1, 0.2, 0.3), dt_approx = 1, test_sensitivity = 0.85) %>% 
      # update_timing(dt = 1) %>% 
      # draw_random_testing(alpha = c(0.1, 0.2, 0.3), dt_approx = 1, test_sensitivity = 0.85) %>% 
      draw_secondary_cases(p_contact_if_isolated_home = c(0.05, 0.1, 0.2), p_contact_traced = c(0.1, 0.2, 0.3),
                           beta_matrix = crossing(to = 1:3, from = 1:3) %>% mutate(beta_val = c(10, 5, 2, 5, 7, 1, 2, 1, 4)),
                           test_results_delay_data = c(1.5, 1, 1),
                           group_props = c(0.25, 0.5, 0.25), 
                           contacts_mean = c(7, 3, 2),
                           sar_home = 0.2, sar_out = 0.05,
                           hh_data = hh_data_sampled)
      # count_prop(is.na(secondary_case_id))
    # dups_report(secondary_case_id) %>% 
    # relocate(secondary_case_id, potential_secondary_cases) %>% 
    # dups_view(secondary_case_id)
    
    
    # TIMING - was 78 seconds with 10k obs and secondary_ids
    # With secondary_ids_vec, it's 4 seconds
    
    draws_debug %>% view()
    
    
    
    
    

# Test outbreak functions -------------------------------------------------

    
    # (1) OUTBREAK SETUP
    probs_default <- c(0.3, 0.6, 0.8, 0.9, 1)
    probs_df_basic <- tibble(
      i_group = c(rep(1, 5), rep(2, 5), rep(3, 5)),
      symptom_severity = as.integer(rep(1:5, 3))
    )
    
    outbreak_setup_test <- outbreak_setup(
      n_pop_total = 100000, 
      n_initial_cases = c(10, 50, 40), group_props = c(0.3, 0.5, 0.2), dt_approx = 1,
      recov_val = 10, # people recover 10 days after symptoms show (cuts off only ~ 1% of infections at the upper end of timing distribution)
      params_timing = data_save$params_timing,
      params_symptom_timing = data_save$params_symptom_timing,
      params_serial = data_save$params_serial,
      hh_size_data = data_save$hh_data_bogota %>% filter(i_group %in% 1:3),
      test_delay_data = data_save$test_delay_data %>% filter(i_group %in% 1:3),
      # test_choice_delay_data = data_save$test_choice_delay_data %>% filter(i_group %in% 1:3), 
      # test_results_delay_data = data_save$test_results_delay_data %>% filter(i_group %in% 1:3), 
      test_sensitivity_data = data_save$test_sensitivity_data, 
      ct_delay_data = ct_delay_data_real,
      probs_self_test_df = probs_df_basic %>% mutate(prob = c(probs_default * 0.1, probs_default * 0.5, probs_default * 1)),
      probs_isolate_symptoms_df = probs_df_basic %>% mutate(prob = c(probs_default * 0.2, probs_default * 0.4, probs_default * 0.8)),
      probs_isolate_test_df = probs_df_basic %>% mutate(prob = c(probs_default * 0.5, probs_default * 0.8, probs_default * 1)),
      probs_isolate_ct_df = probs_df_basic %>% mutate(prob = c(probs_default * 0.6, probs_default * 0.6, probs_default * 1)),
      p_contact_if_isolated_home = c(0.3, 0.2, 0.1), 
      p_contact_traced = c(0.1, 0.2, 0.2),
      p_hh_quarantine = c(0.4, 0.5, 0.5),
      k_matrix = k_matrix_3,
      contact_dispersion = 0.2,
      sar_home = c(0.4, 0.3, 0.3), sar_out = c(0.2, 0.2, 0.2), 
      infectiousness_by_symptoms = c(0.7, 1, 1, 1, 1),
      alpha = c(0.1, 0.1, 0.1)
    )
    
    
    outbreak_setup_test$secondary_cases %>% count_nas()
    
    # (2) SINGLE OUTBREAK SIM
    # single_sim_test <- outbreak_sim(
    #   n_iterations = 10,
    #   n_initial_cases = c(10, 50, 40), group_props = c(0.3, 0.5, 0.2), n_pop = 10000, dt_approx = 1,
    #   test_results_delay_data = 2, recov_val = 7,
    #   probs_self_test_df = probs_df_basic %>% mutate(prob = c(probs_default * 0.1, probs_default * 0.5, probs_default * 1)),
    #   probs_isolate_symptoms_df = probs_df_basic %>% mutate(prob = c(probs_default * 0.1, probs_default * 0.3, probs_default * 0.5)),
    #   probs_isolate_test_df = probs_df_basic %>% mutate(prob = c(probs_default * 0.2, probs_default * 0.2, probs_default * 1)),
    #   probs_isolate_ct_df = probs_df_basic %>% mutate(prob = c(probs_default * 0.2, probs_default * 0.2, probs_default * 1)),
    #   test_sensitivity = 0.85,
    #   p_contact_if_isolated_home = c(0.1, 0.2, 0.3), p_contact_traced = c(0.1, 0.2, 0.4),
    #   beta_matrix = crossing(to = 1:3, from = 1:3) %>% mutate(beta_val = c(10, 5, 2, 5, 7, 1, 2, 1, 4)),
    #   r0_group_1 = 2.5, dispersion = 0.16, alpha =  c(0.05, 0.1, 0.15)
    # )
    
    
    
    
    
    
    
    
    
  
    
# ** Test all functions -------------------------------------------------------------
    
    # QUICK TEST
    # probs_df_basic_4 <- crossing(i_group = 1:4, symptom_severity = as.integer(1:5))
    # # probs_df_basic <- crossing(i_group = 1:3, symptom_severity = as.integer(1:5))
    # probs_default <- c(0.3, 0.6, 0.8, 0.9, 1)
  
    
    # k_matrix_3 <- Matrix::forceSymmetric(
    #   matrix(
    #     c(10, 1, 1,
    #       5,  9, 1,
    #       3,  4, 7),
    #     nrow = 3
    #   )
    # ) %>% as.matrix()
    

    
    
    speedy_test <- outbreak_sims(
      n_sims = 2,
      n_iterations = 300,
      n_pop_total = 10000, 
      n_initial_cases = c(2, 3, 4, 1), group_props = data_save$group_props, dt_approx = 1,
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
      p_contact_traced = data_save$params_data$p_contact_traced,
      p_hh_quarantine = data_save$params_data$p_hh_quarantine,
      k_matrix = data_save$k_matrix * 5,
      contact_dispersion = data_save$contact_dispersion,
      sar_home = data_save$params_data$sar_home_input, 
      sar_out = data_save$params_data$sar_out_input, 
      infectiousness_by_symptoms = c(0.7, 1),
      # alpha = c(0.1, 0.1, 0.1)
      alpha = c(0, 0, 0, 0)
    )
    
    plot(speedy_test)
    speedy_test$time_series %>% 
      group_by(sim_id, t) %>% 
      summarise(n_confirmed_cases = sum(n_confirmed_cases)) %>% 
      ggplot(aes(colour = factor(sim_id), x = t, y = n_confirmed_cases)) + 
      geom_line()
    
    plot_detected_by_i_group(speedy_test)
    
    plot_by_i_group(speedy_test)
    plot(speedy_test)
    
    speedy_test$time_series %>% 
      group_by(sim_id, t) %>% 
      summarise(n_confirmed_cases = sum(n_confirmed_cases),
                n_cases_cum = sum(n_cases_cum)) %>% 
      mutate(share_confirmed = n_confirmed_cases / n_cases_cum) %>% 
      print_all

    
    fspeedy_test$time_series %>% group_by(i_group) %>% quantile_summarise(n_susceptible, n_detected)
    speedy_test$outbreak_t[[1]]$ind_case_counter
    
    plot(speedy_test, indiv_sims = TRUE)
    plot(speedy_test, indiv_sims = FALSE)
    plot_detected(speedy_test, prop = FALSE)
    plot_by_i_group(speedy_test, stacked = FALSE, prop = FALSE)
    plot_detected_by_i_group(speedy_test, prop = TRUE)
    
    
    # SLOWER TEST
    epidemic_1_multi <- outbreak_sims(
      n_sims = 3,
      n_iterations = 100,
      n_pop_total = 100000, 
      hh_size_data = data_save$hh_data_bogota %>% filter(i_group %in% 1:3),
      # hh_size_data = data_save$hh_data_bogota %>% filter(i_group == 3) %>% rep_tibble() %>%  count_prop(i_group),
      n_initial_cases = c(10, 50, 40), group_props = c(0.3, 0.5, 0.2), dt_approx = 1,
      recov_val = 10, # people recover 10 days after symptoms show (cuts off only ~ 1% of infections at the upper end of timing distribution)
      params_timing = data_save$params_timing,
      params_symptom_timing = data_save$params_symptom_timing,
      params_serial = data_save$params_serial,
      test_choice_delay_data = data_save$test_choice_delay_data %>% filter(i_group %in% 1:3), 
      test_results_delay_data = data_save$test_results_delay_data %>% filter(i_group %in% 1:3), 
      test_sensitivity_data = data_save$test_sensitivity_data, 
      ct_delay_data = ct_delay_data_real,
      probs_self_test_df = probs_df_basic %>% mutate(prob = c(probs_default * 0.1, probs_default * 0.5, probs_default * 1)),
      probs_isolate_symptoms_df = probs_df_basic %>% mutate(prob = c(probs_default * 0.2, probs_default * 0.4, probs_default * 0.8)),
      probs_isolate_test_df = probs_df_basic %>% mutate(prob = c(probs_default * 0.5, probs_default * 0.8, probs_default * 1)),
      probs_isolate_ct_df = probs_df_basic %>% mutate(prob = c(probs_default * 0.6, probs_default * 0.6, probs_default * 1)),
      p_contact_if_isolated_home = c(0.3, 0.2, 0.1), 
      p_contact_traced = c(0.1, 0.2, 0.2),
      p_hh_quarantine = c(0.4, 0.5, 0.5),
      beta_matrix = crossing(to = 1:3, from = 1:3) %>% mutate(beta_val = c(10, 5, 2, 5, 7, 1, 2, 1, 4)),
      # r0_group_1 = 5, 
      contacts_mean = c(30, 20, 20),
      sar_home = c(0.4, 0.3, 0.3), sar_out = c(0.2, 0.2, 0.2), 
      infectiousness_by_symptoms = c(0.7, 1, 1, 1, 1),
      alpha = c(0, 0, 0)
    )
    
    
    epidemic_1_multi$time_series %>% tail(20) %>% print_all
    plot(epidemic_1_multi, indiv_sims = FALSE, conf_level = 0.99)
    plot(epidemic_1_multi, indiv_sims = TRUE)
    plot_by_i_group(epidemic_1_multi, stacked = FALSE)
    plot_detected(epidemic_1_multi, prop = TRUE)
    plot_detected_by_i_group(epidemic_1_multi, prop = FALSE)
    
    # PLOT CONFIRMED CASES
    speedy_test %>% 
      .$time_series %>% 
      # group_by(t) %>% summarise(n_confirmed_cases = sum(n_confirmed_cases)) %>% 
      ggplot(aes(x = t, y = n_confirmed_cases, colour = factor(i_group))) + 
      geom_line()
    
    plot_confirmed_cases <- function(outbreak) {
      outbreak$time_series %>% 
        group_by(i_group) %>% arrange(i_group, t) %>% 
        mutate(weekly_confirmed_cases = n_confirmed_cases - lag(n_confirmed_cases, 7)) %>% 
        # summarise(n_confirmed_cases = sum(n_confirmed_cases)) %>%
        group_by(t, i_group) %>% 
        summarise(weekly_confirmed_cases = median(weekly_confirmed_cases)) %>% 
        ggplot(aes(x = t, y = weekly_confirmed_cases, colour = factor(i_group))) + 
        geom_line(show.legend = FALSE)
    }
    
    plot_new_cases <- function(outbreak) {
      epidemic_1_multi$time_series %>% 
        group_by(i_group) %>% arrange(i_group, t) %>% 
        mutate(weekly_new_cases = n_cases_cum - lag(n_cases_cum, 7)) %>% 
        # summarise(n_confirmed_cases = sum(n_confirmed_cases)) %>%
        group_by(t, i_group) %>% 
        summarise(weekly_new_cases = median(weekly_new_cases)) %>% 
        ggplot(aes(x = t, y = weekly_new_cases, colour = factor(i_group))) + 
        geom_line(show.legend = FALSE)
    }


    ggpubr::ggarrange(
      plot_confirmed_cases(epidemic_1_multi) + coord_cartesian(xlim = c(0, 80)),
      plot_new_cases(epidemic_1_multi) + coord_cartesian(xlim = c(0, 80)),
      nrow = 2,
      align = "v"
    )
    
    
  
    
    
    
    # TEST FUNCTION TO LOOK AT EFFECT OF PARAMETERS
    # PROBABILITY OF TESTING (for given symptoms) - lower values means more cases
    probs_self_test_default <- c(0, 0.04, 0.1, 0.15, 0.4)
    
    # SPEEDY 
    speedy_test_params <- outbreak_sims_by_params(
      n_sims = 3,
      n_iterations = 20,
      n_initial_cases = list(c(10, 50, 40)), 
      group_props = list(c(0.3, 0.5, 0.2)),
      n_pop = 10000, 
      dt_approx = 1,
      test_results_delay_data = 2, 
      recov_val = 7,
      test_sensitivity = c(0.5, 1),
      probs_self_test_df = list(probs_df_basic %>% mutate(prob = c(probs_default * 0.8, probs_default * 0.9, probs_default * 1))),
      # probs_self_test_df = list(
      #   "low" = probs_df_basic %>% mutate(prob = c(probs_default * 0.1, probs_default * 0.1, probs_default * 0.1)),
      #   "medium" = probs_df_basic %>% mutate(prob = c(probs_default * 0.5, probs_default * 0.5, probs_default * 0.5)),
      #   "high" = probs_df_basic %>% mutate(prob = c(probs_default * 1, probs_default * 1, probs_default * 1))
      # ),
      probs_isolate_symptoms_df = list(probs_df_basic %>% mutate(prob = c(probs_default * 0.1, probs_default * 0.3, probs_default * 0.5))),
      probs_isolate_test_df = list(probs_df_basic %>% mutate(prob = c(probs_default * 0.2, probs_default * 0.2, probs_default * 1))),
      probs_isolate_ct_df = list(probs_df_basic %>% mutate(prob = c(probs_default * 0.2, probs_default * 0.2, probs_default * 1))),
      p_contact_if_isolated_home = list(c(0.3, 0.2, 0.1)), 
      p_contact_traced = list(c(0.1, 0.2, 0.4)),
      beta_matrix = list(crossing(to = 1:3, from = 1:3) %>% mutate(beta_val = c(10, 5, 2, 5, 7, 1, 2, 1, 4))),
      r0_group_1 = 3, 
      dispersion = 0.16, 
      alpha =  list(c(0.05, 0.1, 0.15))
    )
    
    
    
    
    
    
    
    plot(data = speedy_test_params, param_groups = c("test_sensitivity"), indiv_sims = FALSE, conf_level = 0.9)
    
    plot_detected(data = speedy_test_params, param_groups = c("test_sensitivity"), prop = TRUE)
    
    
    
    
 
    
    
    

# Test policy change functions --------------------------------------------

    
    
    
    policy_change_test <- outbreak_sim_policy_change(
      n_iterations = 100, keep_all_data = FALSE, record_times = FALSE, 
      constant_params = list(
        n_pop_total = 500, 
        hh_size_data = list(
          `1` = rbinom(1000, size = 10, prob = 0.3) + 1,
          `2` = rbinom(1000, size = 5, prob = 0.5) + 1,
          `3` = rbinom(1000, size = 5, prob = 0.5) + 1
        ),
        n_initial_cases = c(10, 50, 40), group_props = c(0.3, 0.2, 0.5), dt_approx = 1,
        recov_val = 7,
        test_sensitivity = 0.85, 
        infectiousness_by_symptoms = c(0.7, 1, 1, 1, 1)
      ),
      policy_start_times = c(0, 10, 40),
      time_varying_params = list(
        "prelockdown" = list(
          test_results_delay_data = c(1.5, 1.5, 1.5),
          probs_self_test_df = probs_df_basic %>% mutate(prob = c(probs_default * 0.1, probs_default * 0.5, probs_default * 1)),
          probs_isolate_symptoms_df = probs_df_basic %>% mutate(prob = c(probs_default * 0.2, probs_default * 0.4, probs_default * 0.6)),
          probs_isolate_test_df = probs_df_basic %>% mutate(prob = c(probs_default * 0.5, probs_default * 0.8, probs_default * 1)),
          probs_isolate_ct_df = probs_df_basic %>% mutate(prob = c(probs_default * 0.6, probs_default * 0.6, probs_default * 1)),
          p_contact_if_isolated_home = c(0.3, 0.2, 0.1), 
          p_contact_traced = c(0.1, 0.2, 0.2),
          p_hh_quarantine = c(0.4, 0.5, 0.5),
          beta_matrix = crossing(to = 1:3, from = 1:3) %>% mutate(beta_val = c(10, 5, 2, 5, 7, 1, 2, 1, 4)),
          r0_group_1 = 5, contacts_mean = c(30, 20, 20),
          secondary_timing_params = inf_profile_params,       sar_home = c(0.4, 0.3, 0.3), sar_out = c(0.2, 0.2, 0.2), 
          alpha =  c(0, 0, 0)
        ),
        "lockdown" = list(
          test_results_delay_data = c(1.5, 1.5, 1.5),
          probs_self_test_df = probs_df_basic %>% mutate(prob = c(probs_default * 0.1, probs_default * 0.5, probs_default * 1)),
          probs_isolate_symptoms_df = probs_df_basic %>% mutate(prob = c(probs_default * 0.2, probs_default * 0.4, probs_default * 0.6)),
          probs_isolate_test_df = probs_df_basic %>% mutate(prob = c(probs_default * 0.5, probs_default * 0.8, probs_default * 1)),
          probs_isolate_ct_df = probs_df_basic %>% mutate(prob = c(probs_default * 0.6, probs_default * 0.6, probs_default * 1)),
          p_contact_if_isolated_home = c(0.3, 0.2, 0.1), 
          p_contact_traced = c(0.1, 0.2, 0.2),
          p_hh_quarantine = c(0.4, 0.5, 0.5),
          beta_matrix = crossing(to = 1:3, from = 1:3) %>% mutate(beta_val = c(10, 5, 2, 5, 7, 1, 2, 1, 4)),
          r0_group_1 = 5, contacts_mean = c(2, 2, 2),
          secondary_timing_params = inf_profile_params,       sar_home = c(0.4, 0.3, 0.3), sar_out = c(0.2, 0.2, 0.2), 
          alpha =  c(0, 0, 0)
        ),
        "postlockdown" = list(
          test_results_delay_data = c(1.5, 1.5, 1.5),
          probs_self_test_df = probs_df_basic %>% mutate(prob = c(probs_default * 0.1, probs_default * 0.5, probs_default * 1)),
          probs_isolate_symptoms_df = probs_df_basic %>% mutate(prob = c(probs_default * 0.2, probs_default * 0.4, probs_default * 0.6)),
          probs_isolate_test_df = probs_df_basic %>% mutate(prob = c(probs_default * 0.5, probs_default * 0.8, probs_default * 1)),
          probs_isolate_ct_df = probs_df_basic %>% mutate(prob = c(probs_default * 0.6, probs_default * 0.6, probs_default * 1)),
          p_contact_if_isolated_home = c(0.3, 0.2, 0.1), 
          p_contact_traced = c(0.1, 0.2, 0.2),
          p_hh_quarantine = c(0.4, 0.5, 0.5),
          beta_matrix = crossing(to = 1:3, from = 1:3) %>% mutate(beta_val = c(10, 5, 2, 5, 7, 1, 2, 1, 4)),
          r0_group_1 = 5, contacts_mean = c(30, 20, 20),
          secondary_timing_params = inf_profile_params,       sar_home = c(0.4, 0.3, 0.3), sar_out = c(0.2, 0.2, 0.2), 
          alpha =  c(0, 0, 0)
        )
      )
    )
    
    
    
    
    # FOR PLOTTING the policy change graphs
    lockdown_rect <- policy_change_test$time_series %>% 
      summarise(xmin = min(t[policy == "lockdown"], na.rm = TRUE), 
                xmax = max(t[policy == "lockdown"], na.rm = TRUE))
    
    policy_change_test$time_series %>% 
      group_by(t, policy) %>% 
      summarise(n_cases_live = sum(n_cases_live)) %>% 
      ggplot() + 
      geom_rect(data = lockdown_rect,
                aes(xmin = xmin,
                    xmax = xmax,
                    ymin=0, ymax=Inf), alpha = 0.1, fill = "indianred") +
      geom_vline(aes(xintercept = min(t[policy == "lockdown"], na.rm = TRUE)), linetype = "dashed") +
      geom_vline(aes(xintercept = max(t[policy == "lockdown"], na.rm = TRUE)), linetype = "dashed") + 
      geom_line(aes(x = t, y = n_cases_live))
    
    
    
    policy_change_multi_test <- outbreak_sims_policy_change(
      n_sims = 5, keep_all_data = FALSE, record_times = FALSE, 
      n_iterations = 70, 
      constant_params = list(
        n_pop_total = 10000, 
        hh_size_data = list(
          `1` = rbinom(10000, size = 10, prob = 0.3) + 1,
          `2` = rbinom(10000, size = 5, prob = 0.5) + 1,
          `3` = rbinom(10000, size = 5, prob = 0.5) + 1
        ),
        n_initial_cases = c(10, 50, 40), group_props = c(0.3, 0.2, 0.5), dt_approx = 1,
        recov_val = 7,
        test_sensitivity = 0.85, 
        infectiousness_by_symptoms = c(0.7, 1, 1, 1, 1)
      ),
      policy_start_times = c(0, 10, 30),
      time_varying_params = list(
        "prelockdown" = list(
          test_results_delay_data = c(1.5, 1.5, 1.5),
          probs_self_test_df = probs_df_basic %>% mutate(prob = c(probs_default * 0.1, probs_default * 0.5, probs_default * 1)),
          probs_isolate_symptoms_df = probs_df_basic %>% mutate(prob = c(probs_default * 0.2, probs_default * 0.4, probs_default * 0.6)),
          probs_isolate_test_df = probs_df_basic %>% mutate(prob = c(probs_default * 0.5, probs_default * 0.8, probs_default * 1)),
          probs_isolate_ct_df = probs_df_basic %>% mutate(prob = c(probs_default * 0.6, probs_default * 0.6, probs_default * 1)),
          p_contact_if_isolated_home = c(0.3, 0.2, 0.1), 
          p_contact_traced = c(0.1, 0.2, 0.2),
          p_hh_quarantine = c(0.4, 0.5, 0.5),
          beta_matrix = crossing(to = 1:3, from = 1:3) %>% mutate(beta_val = c(10, 5, 2, 5, 7, 1, 2, 1, 4)),
          r0_group_1 = 5, contacts_mean = c(30, 20, 20),
          secondary_timing_params = inf_profile_params,       sar_home = c(0.4, 0.3, 0.3), sar_out = c(0.2, 0.2, 0.2), 
          alpha =  c(0, 0, 0)
        ),
        "lockdown" = list(
          test_results_delay_data = c(1.5, 1.5, 1.5),
          probs_self_test_df = probs_df_basic %>% mutate(prob = c(probs_default * 0.1, probs_default * 0.5, probs_default * 1)),
          probs_isolate_symptoms_df = probs_df_basic %>% mutate(prob = c(probs_default * 0.2, probs_default * 0.4, probs_default * 0.6)),
          probs_isolate_test_df = probs_df_basic %>% mutate(prob = c(probs_default * 0.5, probs_default * 0.8, probs_default * 1)),
          probs_isolate_ct_df = probs_df_basic %>% mutate(prob = c(probs_default * 0.6, probs_default * 0.6, probs_default * 1)),
          p_contact_if_isolated_home = c(0.3, 0.2, 0.1), 
          p_contact_traced = c(0.1, 0.2, 0.2),
          p_hh_quarantine = c(0.4, 0.5, 0.5),
          beta_matrix = crossing(to = 1:3, from = 1:3) %>% mutate(beta_val = c(10, 5, 2, 5, 7, 1, 2, 1, 4)),
          r0_group_1 = 5, contacts_mean = c(2, 2, 2),
          secondary_timing_params = inf_profile_params,       sar_home = c(0.4, 0.3, 0.3), sar_out = c(0.2, 0.2, 0.2), 
          alpha =  c(0, 0, 0)
        ),
        "postlockdown" = list(
          test_results_delay_data = c(1.5, 1.5, 1.5),
          probs_self_test_df = probs_df_basic %>% mutate(prob = c(probs_default * 0.1, probs_default * 0.5, probs_default * 1)),
          probs_isolate_symptoms_df = probs_df_basic %>% mutate(prob = c(probs_default * 0.2, probs_default * 0.4, probs_default * 0.6)),
          probs_isolate_test_df = probs_df_basic %>% mutate(prob = c(probs_default * 0.5, probs_default * 0.8, probs_default * 1)),
          probs_isolate_ct_df = probs_df_basic %>% mutate(prob = c(probs_default * 0.6, probs_default * 0.6, probs_default * 1)),
          p_contact_if_isolated_home = c(0.3, 0.2, 0.1), 
          p_contact_traced = c(0.1, 0.2, 0.2),
          p_hh_quarantine = c(0.4, 0.5, 0.5),
          beta_matrix = crossing(to = 1:3, from = 1:3) %>% mutate(beta_val = c(10, 5, 2, 5, 7, 1, 2, 1, 4)),
          r0_group_1 = 5, contacts_mean = c(30, 20, 20),
          secondary_timing_params = inf_profile_params,       sar_home = c(0.4, 0.3, 0.3), sar_out = c(0.2, 0.2, 0.2), 
          alpha =  c(0, 0, 0)
        )
      )
    )
    
    
    
    # FOR PLOTTING the policy change graphs
    lockdown_rect <- policy_change_multi_test$time_series %>% 
      summarise(xmin = min(t[policy == "lockdown"], na.rm = TRUE), 
                xmax = max(t[policy == "lockdown"], na.rm = TRUE))
    
    policy_change_multi_test$time_series %>% 
      group_by(sim_id, t, policy) %>% 
      summarise(n_cases_live = sum(n_cases_live)) %>% 
      group_by(t, policy) %>% 
      quantile_summarise(n_cases_live, conf_level = 0.8) %>% 
      ggplot() + 
      geom_rect(data = lockdown_rect,
                aes(xmin = xmin,
                    xmax = xmax,
                    ymin=0, ymax=Inf), alpha = 0.1, fill = "indianred") +
      geom_vline(aes(xintercept = min(t[policy == "lockdown"], na.rm = TRUE)), linetype = "dashed") +
      geom_vline(aes(xintercept = max(t[policy == "lockdown"], na.rm = TRUE)), linetype = "dashed") + 
      geom_line(aes(x = t, y = n_cases_live_median)) + 
      geom_ribbon(aes(x = t, ymax = n_cases_live_upper, ymin = n_cases_live_lower), alpha = 0.2)
    
    
    
    
    
    
    
    
    

# Test with actual hh size data (with 6 groups) -------------------------------------------

    
    beta_matrix_6 <- Matrix::forceSymmetric(
      matrix(
        c(10, 1, 1, 1, 1, 1,
          4,  6, 1, 1, 1, 1,
          3,  4, 7, 1, 1, 1,
          1,  3, 3, 4, 1, 1,
          1,  1, 2, 2, 3, 1, 
          0,  0, 1, 1, 2, 3),
        nrow = 6
      )
    ) %>% 
      as.vector() %>% 
      {mutate(crossing(to = 1:6, from = 1:6), beta_val = .)}
    
    
    probs_df_basic_6 <- tibble(
      i_group = c(rep(1, 5), rep(2, 5), rep(3, 5), rep(4, 5), rep(5, 5), rep(6, 5)),
      symptom_severity = as.integer(rep(1:5, 6))
    ) %>% 
      mutate(
        prob = rep(probs_default, 6)
      )
    
    test_bogota_pop <- outbreak_sims(
      n_sims = 1,
      n_iterations = 20,
      n_pop_total = n_pop_bogota, 
      hh_size_data = hh_data_bogota,
      sample_hh_data = FALSE, # use the full hh data, rather than sampling
      n_initial_cases = c(10, 50, 40, 20, 20, 20), 
      group_props = group_props_bogota, 
      dt_approx = 1,
      recov_val = 7,
      test_results_delay_data = test_results_delay_data_6, 
      test_sensitivity_data = test_sensitivity_data, 
      ct_delay_data = ct_delay_data_real_6,
      probs_self_test_df =        probs_df_basic_6 %>% mutate(prob = prob * 0.4),
      probs_isolate_symptoms_df = probs_df_basic_6 %>% mutate(prob = prob * 0.6),
      probs_isolate_test_df =     probs_df_basic_6 %>% mutate(prob = prob * 0.9),
      probs_isolate_ct_df =       probs_df_basic_6 %>% mutate(prob = prob * 0.5),
      p_contact_if_isolated_home = c(0.7, 0.5, 0.5, 0.5, 0.5, 0.5), 
      p_contact_traced = c(0.1, 0.2, 0.2, 0.2, 0.3, 0.3),
      p_hh_quarantine = c(0.4, 0.5, 0.5, 0.5, 0.5, 0.5),
      beta_matrix = beta_matrix_6,
      # r0_group_1 = 5, 
      contacts_mean = c(30, 20, 20, 15, 15, 15),
      secondary_timing_params = inf_profile_params,
      sar_home = c(0.4, 0.3, 0.3, 0.2, 0.2, 0.2), sar_out = c(0.2, 0.2, 0.2, 0.2, 0.2, 0.2), 
      infectiousness_by_symptoms = c(0.7, 1, 1, 1, 1),
      # alpha =  c(0.05, 0.1, 0.1, 0.1)
      alpha = c(0, 0, 0, 0, 0, 0)
    )
    
    
    test_bogota_pop %>% plot(indiv_sims = TRUE)
    test_bogota_pop %>% plot_by_i_group()
    test_bogota_pop %>% plot_incidence()
    
    
    
    
# Look at parameter effects 1 by 1 ----------------------------------------
    
    
    
    # PROBABILITY OF SELF-TESTING - if people are more likely to get tested, the virus spreads less
    effect_probs_self_test <- outbreak_sims_by_params(
      n_sims = 10,
      n_iterations = 70,
      n_initial_cases = list(c(10, 50, 40)), 
      group_props = list(c(0.3, 0.5, 0.2)), n_pop = 100000, dt_approx = 1,
      test_results_delay_data = 1, recov_val = 7,
      probs_self_test_df = list(
        "1_low" = probs_df_basic %>% mutate(prob = 0.2),
        "2_medium" = probs_df_basic %>% mutate(prob = 0.5),
        "3_high" = probs_df_basic %>% mutate(prob = 0.9)
      ),
      probs_isolate_symptoms_df = list(
        probs_df_basic %>% mutate(prob = c(probs_default * 0, probs_default * 0, probs_default * 0))
      ),
      probs_isolate_test_df = list(
        probs_df_basic %>% mutate(prob = c(probs_default * 1, probs_default * 1, probs_default * 1))
      ),
      probs_isolate_ct_df = list(
        probs_df_basic %>% mutate(prob = c(probs_default * 0.1, probs_default * 0.1, probs_default * 0.1))
      ),
      test_sensitivity = 0.9,
      p_contact_if_isolated_home = list(c(0.1, 0.1, 0.1)), p_contact_traced = list(c(0.2, 0.3, 0.4)),
      beta_matrix = list(
        crossing(to = 1:3, from = 1:3) %>% mutate(beta_val = c(10, 5, 2, 5, 7, 1, 2, 1, 4))
      ),
      r0_group_1 = 2.5, dispersion = 0.16, alpha =  list(c(0.05, 0.1, 0.1))
    )
    
    
    plot(effect_probs_self_test,  
         param_groups = "probs_self_test_df", indiv_sims = FALSE,
         conf_level = 0.9)
    
    plot_detected(effect_probs_self_test, param_groups = "probs_self_test_df", prop = TRUE)
    
    
    # PROBABILITY OF ISOLATING - if people are less likely to isolate, the virus spreads more
    effect_probs_isolate_symptoms <- outbreak_sims_by_params(
      n_sims = 10, dt_approx = 1, n_iterations = 10,
      probs_self_test = list(probs_self_test_default), 
      probs_isolate_symptoms = list(probs_self_test_default, c(0.7, 0.8, 0.9, 1, 1)),
      probs_isolate_test = list(c(0.4, 0.6, 0.9, 0.95, 0.99)),
      test_results_delay_data = 0,
      recov_val = 5, r0 = 2, dispersion = 0.16,
      p_contact_if_isolated_home = c(0.1)
    )
    
    plot(effect_probs_isolate_symptoms, probs_isolate_symptoms)
    
    # PROBABILITY OF ISOLATING AFTER TEST - if people are less likely to isolate, the virus spreads more
    effect_probs_isolate_test <- outbreak_sims_by_params(
      n_sims = 10, dt_approx = 1, n_iterations = 10,
      probs_self_test = list(0.5 + (1:5) / 10), 
      probs_isolate_symptoms = list(probs_self_test_default),
      probs_isolate_test = list(1:5 / 10, c(1, 1, 1, 1, 1)),
      test_results_delay_data = 0,
      recov_val = 5, r0 = 2, dispersion = 0.16,
      p_contact_if_isolated_home = c(0.1)
    )
    
    plot(effect_probs_isolate_test, probs_isolate_test)
    
    
    
    
    # test_results_delay - higher means worse outbreak
    effect_test_results_delay <- outbreak_sims_by_params(
      n_sims = 10,
      n_iterations = 50,
      n_initial_cases = list(c(10, 50, 40)), 
      group_props = list(c(0.3, 0.5, 0.2)), n_pop = 10000, dt_approx = 1,
      test_results_delay_data = c(0, 0.5, 1, 2, 5), recov_val = 7,
      probs_self_test_df = list(
        probs_df_basic %>% mutate(prob = c(probs_default * 0.4, probs_default * 0.7, probs_default * 0.9))
      ),
      probs_isolate_symptoms_df = list(
        probs_df_basic %>% mutate(prob = c(probs_default * 0.1, probs_default * 0.3, probs_default * 0.5))
      ),
      probs_isolate_test_df = list(
        probs_df_basic %>% mutate(prob = c(probs_default * 1, probs_default * 1, probs_default * 1))
      ),
      p_contact_if_isolated_home = list(c(0.3, 0.2, 0.1)), p_contact_traced = list(c(0.1, 0.2, 0.4)),
      beta_matrix = list(
        crossing(to = 1:3, from = 1:3) %>% mutate(beta_val = c(10, 5, 2, 5, 7, 1, 2, 1, 4))
      ),
      r0_group_1 = 2.5, dispersion = 0.16, alpha =  list(c(0.05, 0.1, 0.1))
    )
    
    plot(effect_test_results_delay, param_group = test_results_delay_data)
    
    
    # p_contact_if_isolated_home - higher should mean worse outbreak
    effect_isolated <- outbreak_sims_by_params(
      n_sims = 20, dt_approx = 1, n_iterations = 40, 
      probs_self_test = list(probs_self_test_default) ,    # << high probability of testing 
      probs_isolate_symptoms = list(c(0, 0.04, 0.1, 0.15, 0.4)),
      probs_isolate_test = list(c(0.6, 0.6, 0.9, 0.95, 0.99)), 
      p_contact_if_isolated_home = 0.1, 
      p_contact_traced = 0.2,
      test_results_delay_data = c(0, 1, 3, 5, 10),
      recov_val = 5, r0 = 2, dispersion = 0.16
    )
    
    plot(effect_isolated, p_contact_if_isolated_home)
    
    
    # R0 - higher R0 is worse
    effect_r0 <- outbreak_sims_by_params(
      n_sims = 20, dt_approx = 1, n_iterations = 50, test_thresh = 0.35, p_contact_if_isolated_home = 0.35,
      test_results_delay_data = 2, recov_val = 7, dispersion = 0.16, r0 = c(0.8, 1.5, 2.5, 5)
    )
    
    plot(effect_r0, r0)
    
    
    
    # Recov val - higher recov_val is worse (slower recovery, more infectiousness)
    effect_recov_val <- outbreak_sims_by_params(
      n_sims = 20, dt_approx = 1, n_iterations = 50, test_thresh = 0.35, p_contact_if_isolated_home = 0.35,
      test_results_delay_data = 2, recov_val = c(3, 6.5, 10), dispersion = 0.16, r0 = 2
    )
    
    plot(effect_recov_val, recov_val)
    
    
    
    # Contact tracing effectiveness - higher p means fewer cases. 
    # You need to have a high probability of getting tested otherwise contact tracing is basically useless
    effect_contact_tracing_prob <- outbreak_sims_by_params(
      n_sims = 30, dt_approx = 1, n_iterations = 30,
      probs_self_test = list(c(0.2, 0.5, 0.8, 0.9, 1)) ,    # << high probability of testing 
      probs_isolate_symptoms = list(c(0, 0.04, 0.1, 0.15, 0.4)),
      probs_isolate_test = list(c(0.4, 0.6, 0.9, 0.95, 0.99)), 
      p_contact_if_isolated_home = 0.1, 
      p_contact_traced = c(0.1, 0.4, 0.9),
      test_results_delay_data = 2,
      recov_val = 5, r0 = 2, dispersion = 0.16
    )
    
    plot(effect_contact_tracing_prob, p_contact_traced)
    
    effect_contact_tracing_prob$epidemic_results[[1]]$outbreak_t[[20]]$live_cases %>% 
      view()
    
    
    
    # Effect of alpha
    effect_random_testing <- outbreak_sims_by_params(
      n_sims = 30, dt_approx = 1, n_iterations = 30,
      probs_self_test = list(c(0.2, 0.5, 0.8, 0.9, 1)) ,    # << high probability of testing 
      probs_isolate_symptoms = list(c(0, 0.04, 0.1, 0.15, 0.4)),
      probs_isolate_test = list(c(0.7, 0.8, 0.9, 0.95, 0.99)), 
      p_contact_if_isolated_home = 0.1, 
      p_contact_traced = 0.3,
      test_results_delay_data = 2,
      alpha = c(0, 0.3, 0.6),
      recov_val = 5, r0 = 3, dispersion = 0.16
    )
    
    plot(effect_random_testing, param_group = alpha)
    
    
    
    probs_default_high <- c(0.5, 0.7, 0.85, 0.9, 1)
    
    # Effect of test_sensitivity
    effect_test_sensitivity <- outbreak_sims_by_params(
      n_sims = 20,
      n_iterations = 60,
      n_initial_cases = list(c(10, 50, 40)), 
      group_props = list(c(0.3, 0.5, 0.2)), n_pop = 10000, dt_approx = 1,
      test_results_delay_data = 1.5, recov_val = 7,
      test_sensitivity = c(0.5, 0.8, 1),
      probs_self_test_df = list(
        probs_df_basic %>% mutate(prob = c(probs_default_high * 0.5, probs_default_high * 0.8, probs_default_high * 1))
      ),
      probs_isolate_symptoms_df = list(
        probs_df_basic %>% mutate(prob = c(probs_default_high * 0.1, probs_default_high * 0.2, probs_default_high * 0.3))
      ),
      probs_isolate_test_df = list(
        probs_df_basic %>% mutate(prob = c(probs_default_high * 0.8, probs_default_high * 0.9, probs_default_high * 1))
      ),
      probs_isolate_ct_df = list(
        probs_df_basic %>% mutate(prob = c(probs_default_high * 0.1, probs_default_high * 0.2, probs_default_high * 0.3))
      ),
      p_contact_if_isolated_home = list(c(0.2, 0.1, 0)), p_contact_traced = list(c(0.1, 0.2, 0.4)),
      beta_matrix = list(
        crossing(to = 1:3, from = 1:3) %>% mutate(beta_val = c(10, 5, 2, 5, 7, 1, 2, 1, 4))
      ),
      r0_group_1 = 2.5, dispersion = 0.16, alpha =  list(c(0.05, 0.1, 0.1))
    )
    
    
    plot(effect_test_sensitivity, param_groups = "test_sensitivity", indiv_sims = FALSE, conf_level = 0.5)
    
    plot_detected(effect_test_sensitivity, param_groups = "test_sensitivity", prop = TRUE)
    
    

# .............. ---------------------------------------------------------------------
# CONSISTENCY CHECKS ------------------------------------------------------
    

    consistency_check_data <- outbreak_sims(
      keep_all_data = TRUE,
      n_sims = 1,
      n_iterations = 30,
      n_pop_total = 1000, 
      hh_size_data = data_save$hh_data_bogota %>% filter(i_group %in% 1:3),
      n_initial_cases = c(10, 50, 40), group_props = c(0.3, 0.5, 0.2), dt_approx = 1,
      recov_val = 10, # people recover 10 days after symptoms show (cuts off only ~ 1% of infections at the upper end of timing distribution)
      params_timing = data_save$params_timing,
      params_symptom_timing = data_save$params_symptom_timing,
      params_serial = data_save$params_serial,
      test_delay_data = data_save$test_delay_data %>% filter(i_group %in% 1:3),
      # test_choice_delay_data = data_save$test_choice_delay_data %>% filter(i_group %in% 1:3), 
      # test_results_delay_data = data_save$test_results_delay_data %>% filter(i_group %in% 1:3), 
      test_sensitivity_data = data_save$test_sensitivity_data, 
      ct_delay_data = ct_delay_data_real,
      probs_self_test_df = probs_df_basic %>% mutate(prob = c(probs_default * 0.1, probs_default * 0.5, probs_default * 1)),
      probs_isolate_symptoms_df = probs_df_basic %>% mutate(prob = c(probs_default * 0.2, probs_default * 0.4, probs_default * 0.8)),
      probs_isolate_test_df = probs_df_basic %>% mutate(prob = c(probs_default * 0.5, probs_default * 0.8, probs_default * 1)),
      probs_isolate_ct_df = probs_df_basic %>% mutate(prob = c(probs_default * 0.6, probs_default * 0.6, probs_default * 1)),
      p_contact_if_isolated_home = c(0.3, 0.2, 0.1), 
      p_contact_traced = c(0.1, 0.2, 0.2),
      p_hh_quarantine = c(0.4, 0.5, 0.5),
      k_matrix = data_save$k_matrix_5[1:3, 1:3] / 10,
      contact_dispersion = 0.2,
      sar_home = c(0.4, 0.3, 0.3), sar_out = c(0.2, 0.2, 0.2), 
      infectiousness_by_symptoms = c(0.7, 1, 1, 1, 1),
      alpha = c(0.1, 0.1, 0.1)
    )
    
    value <- function(x) {!is.na(x)}
    
    # Get a dataset including live cases at each moment t
    t_vec <- consistency_check_data$outbreak_t_record[[1]] %>% map_dbl("t")
    
    live_cases_list <- consistency_check_data$outbreak_t_record[[1]] %>% map("live_cases") %>% 
      set_names(t_vec) %>% 
      bind_rows(.id = "t") %>% 
      mutate(t = as.numeric(t), i_group = factor(i_group)) %>% 
      arrange(case_id, t)
    
    secondary_cases_list <- consistency_check_data$outbreak_t_record[[1]] %>% map("secondary_cases") %>% 
      set_names(t_vec) %>% 
      bind_rows(.id = "t") %>% 
      mutate(t = as.numeric(t), i_group = factor(i_group)) %>% 
      arrange(case_id, secondary_case_id, t) %>% 
      relocate(case_id, secondary_case_id, t)
    
    # secondary_cases_list %>% count_nas(sort = FALSE)
    
    # secondary_old_list <- consistency_check_data$outbreak_t_record[[1]] %>% map("secondary_cases_old") %>% 
    #   set_names(t_vec) %>% 
    #   bind_rows(.id = "t") %>% 
    #   mutate(t = as.numeric(t), i_group = factor(i_group)) %>% 
    #   arrange(case_id, secondary_case_id, t) %>% 
    #   relocate(case_id, secondary_case_id, t)
    
    hh_status_list <- consistency_check_data$outbreak_t_record[[1]] %>% map("hh_status") %>% 
      set_names(t_vec) %>% 
      bind_rows(.id = "t") %>% 
      mutate(t = as.numeric(t), i_group = factor(i_group))
    
    # cases_secondary_all_list <- consistency_check_data$outbreak_t_record[[1]] %>% map("cases_secondary_all") %>% 
    #   set_names(t_vec) %>% 
    #   bind_rows(.id = "t") %>% 
    #   mutate(t = as.numeric(t), i_group = factor(i_group)) %>% 
    #   arrange(case_id, secondary_case_id, t)
    
    
    # secondary_cases_list %>% group_by(case_id) %>% 
    #   select(case_id, secondary_case_id, secondary_case_timing, t, new_potential_case, new_actual_case) %>% 
    #   arrange(case_id, secondary_case_id, t) %>% 
    #   print(n = 100)
    #   summarise(n_potential_cases = if_else(
    #     sum(!is.na(secondary_case_id)) == 0, 0, n_distinct(secondary_case_id)
    #     n_actual_cases = sum(new_actual_case == TRUE, na.rm = TRUE)
    #   )
    
    
    # TAKES A WHILE - match variables on live / secondary cases
    all_identical <- function(x, include_NAs = FALSE) {
      if (!include_NAs) x <- x[!is.na(x)]
      length(unique(x)) == 1
    }
    
    # (1) Are there people in secondary cases who aren't in live_cases? - there shouldn't be
    cases_merged <- trackr::full_join_track(live_cases_list, secondary_cases_list,
                                            by = c("case_id", "t"),
                                            suffix = c("_L", "_S"),
                                            .merge = TRUE) # should be none in the y_only section
    # Previous bug - some people were in secondary cases but not live cases when they RECOVER before the secondary case is meant to happen
    
    
    
    # (TODO: I can create this dataset much more easily using bind_rows)
    # variable_match <- cases_merged %>% 
    #   select(t, case_id, secondary_case_id, ends_with("_L"), ends_with("_S")) %>% 
    #   pivot_longer(cols = ends_with("_L") | ends_with("_S"),
    #                names_to = c(".value", "dataset"),
    #                names_pattern = "(.*)_(L|S)") %>% 
    #   relocate(t, case_id, dataset) %>% 
    #   group_by(t, case_id, secondary_case_id) %>% 
    #   # sample_n_groups(100) %>% 
    #   filter(!is.na(secondary_case_id)) %>% 
    #   mutate(
    #     across(
    #       -c(dataset),
    #       list(conflict = ~ !all_identical(., include_NAs = TRUE))
    #     )
    #   )
    
    
    # RUN UNIT TESTS
    library("testthat")
    test_file("unit_tests.R")
    
    
  

# Visual checks -----------------------------------------------------------

    # Histograms of variables
    live_cases_list %>% ggplot(aes(x = symptom_severity)) + geom_bar()
    live_cases_list %>% ggplot(aes(x = self_test_yn)) + geom_bar()
    live_cases_list %>% filter(!is.na(ct_test_negative)) %>% ggplot(aes(x = ct_test_negative)) + geom_bar()
    live_cases_list %>% filter(!is.na(ct_false_negative)) %>% ggplot(aes(x = ct_false_negative)) + geom_bar()
    
    
    # (1) Check testing decision based on symptom severity
    live_cases_list %>% 
      group_by(i_group, symptom_severity) %>% 
      count_prop(self_test_yn, return_count = TRUE) %>% 
      filter(self_test_yn) %>% 
      ggplot(aes(x = factor(i_group), y = prop, fill = factor(symptom_severity)), stat = "identity") + geom_col(position = "dodge")
    
    
    # (2) Check proportion of people who are randomly tested in each period - this goes way above 0.2 because it's "cumulative"
    live_cases_list %>% 
      group_by(i_group, t) %>% 
      summarise(random_test_yn = mean(random_test_yn, na.rm = TRUE)) %>% 
      ggplot(aes(x = t, y = random_test_yn, colour = factor(i_group))) + geom_line()
    
    # (3) Check proportion of "new" random tests
    live_cases_list %>% 
      arrange(case_id, t) %>% 
      mutate(new_random_test = random_test_yn_lag1 == FALSE & random_test_yn == TRUE) %>% 
      group_by(i_group, t) %>% 
      summarise(new_random_test = mean(new_random_test, na.rm = TRUE)) %>% 
      mutate(new_random_test_roll3 = rollmean(new_random_test, k = 3, na.pad = TRUE, align = "center")) %>% 
      ggplot(aes(x = t, y = new_random_test_roll3, colour = factor(i_group))) + geom_line()
    
    # (12) Infection timing always before symptom timing
    live_cases_list %>% hist_basic(x = symptom_timing - infection_timing) # should have none < 0
    
    
    # (13) Check test sensitivity
    live_cases_list %>% group_by(self_test_days_exposure) %>% 
      summarise(self_test_false_negative = mean(self_test_false_negative),
                n = n()) %>% 
      print_all %>% 
      ggplot(aes(x = self_test_days_exposure, y = 1 - self_test_false_negative)) + geom_line() + 
      geom_line(data = data_save$test_sensitivity_data, aes(x = days_exposure, y = sensitivity), linetype = "dashed", colour = "indianred")
    
    live_cases_list %>% group_by(ct_test_days_exposure) %>% summarise(ct_false_negative = mean(ct_false_negative), n = n()) %>% print_all %>% 
      ggplot(aes(x = ct_test_days_exposure, y = 1 - ct_false_negative)) + geom_line() + 
      geom_line(data = data_save$test_sensitivity_data, aes(x = days_exposure, y = sensitivity), linetype = "dashed", colour = "indianred")
    
    
    
    # (13) How many people isolate presymptomatically?
    secondary_cases_list %>% filter(value(isolate_ct_test_timing)) %>% 
      count_prop(isolate_ct_test_timing < symptom_timing) # mostly pre-symptom
    
    secondary_cases_list %>% filter(value(isolate_ct_results_timing)) %>% 
      count_prop(isolate_ct_results_timing < symptom_timing) # about half and half
    
    secondary_cases_list %>% filter(value(isolate_random_test_timing)) %>% 
      count_prop(isolate_random_test_timing < symptom_timing) # VERY often presymptom!
    
    
    
    

# UNIT TESTS [now moved to unit_tests.R] ----------------------------------


# 1. Live case consistency checks --------------------------------------------
    
    # MOVED TO UNIT_TESTS.R
    # # (0) Dups
    # live_cases_list %>% dups_report(case_id, t)
    # 
    # # (4) Check we only have isolate_symptoms_timing for people who actually isolate after symptoms
    # live_cases_list %>% 
    #   count_prop(isolate_after_symptoms, value(isolate_symptoms_timing))
    # 
    # # Check that contact tracing times never stay the same for a given person
    # # Former bug - contact_tracing_timing and contact_tracing_test_timing don't update correctly over periods [should increment by 1 each time]
    # live_cases_list %>% 
    #   group_by(case_id) %>% 
    #   arrange(t) %>% 
    #   mutate(ct_test_timing_bug = contact_tracing_test_timing == lag(contact_tracing_test_timing),
    #          ct_results_timing_bug = contact_tracing_results_timing == lag(contact_tracing_results_timing)) %>% 
    #   ungroup %>% 
    #   count_prop(ct_test_timing_bug, ct_results_timing_bug) # should always be false / NA
    # 
    # # (6) Check that symptom severity never changes for an individual
    # live_cases_list %>% group_by(case_id) %>% summarise(n_distinct = n_distinct(symptom_severity)) %>% 
    #   count_prop(n_distinct) # should only have n_distinct = 1
    # 
    # # (7) Check i_group never changes
    # live_cases_list %>% group_by(case_id) %>% summarise(n_distinct = n_distinct(i_group)) %>% 
    #   count_prop(n_distinct) # should only have n_distinct = 1
    # 
    # # (8) Check currently testing / previously tested 
    # 
    # # No one is both currently_testing and previously_tested? Actually that's OK if they are - people can still self-test/contact traced after being randomly tested
    # sum(live_cases_list$previously_tested & live_cases_list$currently_testing)
    # 
    # # No one is previously tested having NEVER been currently tested before?
    # # YES - people come up as previously tested when they were contact tested BEFORE they were infected...
    # live_cases_list %>% 
    #   # select(case_id, t, currently_testing, previously_tested) %>% 
    #   group_by(case_id) %>% 
    #   mutate(across(c(currently_testing, previously_tested), as.numeric)) %>% 
    #   mutate(ever_tested_debug = cummax(currently_testing)) %>% 
    #   # mutate(testing_offswitch = currently_testing == 0 & lag(currently_testing == 1)) %>% 
    #   mutate(testing_bug = previously_tested == TRUE & ever_tested_debug == FALSE) %>% 
    #   # ungroup %>% 
    #   # count_prop(testing_bug) %>% 
    #   mutate(traced_before_infection = sum(infection_timing > contact_tracing_test_timing, na.rm = TRUE) > 0,
    #          id_has_testing_bug = sum(testing_bug) > 0) %>% 
    #   ungroup %>% 
    #   count_prop(traced_before_infection, id_has_testing_bug) 
    # # there should be no people with traced_before_infection = FALSE and id_has_testing_bug = TRUE
    # 
    # 
    # # (9) [Previous bug] random_test_yn_ever should never be NA
    # live_cases_list %$% sum(is.na(random_test_yn_ever) | is.na(random_test_yn))
    # 
    # 
    # # (10) Check isolation behaviour
    # 
    # 
    # 
    # 
    #     # No cases that have never been randomly tested and have an isolation time due to random testing
    #     live_cases_list %$% sum(!random_test_yn_ever & value(isolate_random_test_timing)) # should be 0
    #     
    #     # No cases that isolate due to random testing and have isolate_after_test is F
    #     live_cases_list %$% sum(!isolate_after_test_presymp & !isolate_after_test_symp & value(isolate_random_test_timing)) # should be 0
    #     
    #     # No cases for whom test is finished, and random_test_yn is still TRUE
    #     # This was not the case when random testing was still done on people who had been tested negative
    #     live_cases_list %$% sum(random_test_timing < 0 & random_test_yn) # should be 0
    #     
    #     live_cases_list %>% count_prop(value(symptom_timing), value(isolate_symptoms_timing), 
    #                                    symptom_timing == isolate_symptoms_timing) # should have no people who are symptom_timing != isolate_symptoms_timing
    #     
    #     
    # # (11) Check testing delays
    # live_cases_list %>% print_names
    # 
    #     # test_results_delay should always be the same as the difference between contact_tracing_test and results
    #     # live_cases_list  %>% mutate(
    #     #   test_delay_contact = contact_tracing_results_timing - contact_tracing_test_timing
    #     # ) %$%
    #     #   sum(round(test_delay_contact, 3) != test_results_delay, na.rm = TRUE) # should be 0
    #     
    #     
    #     # Make sure first random test timing is the same as test_results_delay as well
    #     live_cases_list %>% 
    #       group_by(case_id, test_results_delay, i_group) %>% 
    #       summarise(random_test_timing_max = max_na(random_test_timing), .groups = "drop") %$% 
    #       # count_prop(i_group, test_results_delay, random_test_timing_max) %$%  # make sure these match the test_delay
    #       sum(round(random_test_timing_max, 3) != test_results_delay, na.rm = TRUE) # should be 0
    #     
    # 
    #     
    # # (12) Infection timing always before symptom timing
    # # live_cases_list %>% hist_basic(x = symptom_timing - infection_timing) # should have none < 0
    

  
    
# 2. Secondary case consistency checks ---------------------------------------
    
    
    
    # (1) Number of people who have an ID but don't have a secondary_case_timing - should be 0
    secondary_cases_list %$%
      sum(!is.na(secondary_case_id) & is.na(secondary_case_timing)) # should be 0
    
    # (2) Secondary case timing should always be greater than infection timing
    secondary_cases_list %$%
      sum(secondary_case_timing - infection_timing < 0, na.rm = TRUE) # should be 0
    
    # (2b ) Symptom timing always after infection
    secondary_cases_list %$%
      sum(symptom_timing - infection_timing < 0, na.rm = TRUE) # should be 0
    
    
    # (3) make sure isolation happens only if isolation_timing is less than secondary_case_timing
    # secondary_isolated == TRUE only if isolation_timing < secondary_case_timing == TRUE
    secondary_cases_list %>% 
      count_prop(transmission_isolated, isolation_timing < secondary_case_timing, is.na(isolation_timing) | is.na(secondary_case_timing)) # should always be TRUE at the same time
      # OBSOLETE - you can now be isolated because of household level quarantining **
    
    
    # (4) Check relationship between secondary_isolated and contact_if_isolated
    secondary_cases_list %>% 
      count_prop(transmission_isolated, contact_if_isolated) # contact if isolated should be NA if secondary_isolated is FALSE, and T/F when secondary_isolated is TRUE
    
    # (5) Make sure only people who don't have a secondary_case_id [i.e. "non cases"] are NA for secondary_isolated
    secondary_cases_list %$%
      sum(!is.na(secondary_case_id) & is.na(transmission_isolated))
    
    # (6) Prop susceptible never over 1 or less than 0
    # secondary_cases_list %>% hist_basic(x = prop_susceptible)
    secondary_cases_list %>% filter(is.na(transmission_isolated)) %>% nrow() # should be 0
  
    # (7) [PREVIOUS BUG] Sometimes new_potential_case is NA in cases_secondary_all because
    # all details of the original case (including, importantly, recovery_timing) are missing from this dataset
    # This occurs when the secondary_case_timing is *just* before the recovery but in the same period
    # cases_secondary_all_list %>% 
    #   group_by(case_id, secondary_case_id) %>% 
    #   filter(sum(is.na(i_group)) > 0) %>% 
    #   relocate(case_id, secondary_case_id, t, secondary_case_timing, recovery_timing) %>% 
    #   view()
    
    # (8) Make sure (some) people isolate at contact_tracing_test_timing if they ahve isolate_after_ct
    secondary_cases_list %$% sum((isolate_after_ct_presymp | isolate_after_ct_symp) & isolation_timing == contact_tracing_test_timing, na.rm = TRUE) # should be >0
    
    # (9) Make sure deisolate_timing is not only NA
    secondary_cases_list %$% sum(value(deisolate_timing)) # should be greater than 0
        
        # (9a) Make sure isolation_timing_2 only exists if deisolate_timing exists
        secondary_cases_list %$% sum(value(isolation_timing_2) == TRUE & value(deisolate_timing) == FALSE) # should be 0
        secondary_cases_list %$% sum(value(isolation_timing_2) == TRUE & value(deisolate_timing) == TRUE & deisolate_timing >= isolation_timing_2) # (deisolate timing is always less than isolaton_timing_2) should be 0
        
        # (9b) Make sure secondary case is isolated if it's between isolation and deisolate, or after isolation_timing_2
        secondary_cases_list %$% sum(transmission_isolated == FALSE & infection_type == "out" & secondary_case_timing >= isolation_timing  & secondary_case_timing < deisolate_timing, na.rm = TRUE)  #should be 0
        secondary_cases_list %$% sum(transmission_isolated == FALSE & infection_type == "out" & secondary_case_timing >= isolation_timing_2, na.rm = TRUE)  #should be 0
        secondary_cases_list %$% sum(transmission_isolated == TRUE & secondary_case_timing < isolation_timing, na.rm = TRUE)  #should be 0 - OBSOLETE because of hh quarantining
        secondary_cases_list %$% sum(transmission_isolated == TRUE & secondary_case_timing > deisolate_timing & (!value(isolation_timing_2) | secondary_case_timing < isolation_timing_2), na.rm = TRUE)  #should be 0 - OBSOLETE because of hh quarantining
        
        # (9c) Make sure deisolate timing only exists when isolate_ct_test_timing exists
        secondary_cases_list %$% sum(value(deisolate_timing) & value(isolate_ct_test_timing) == FALSE) # should be 0
        
        # No transmission isolated at home
        secondary_cases_list %$% sum(transmission_isolated == TRUE & hh_id == secondary_hh_id, na.rm = TRUE)  #should be 0
        
    
    # (10) Make sure isolate_ct_test_timing when people are isolated and ct tested
    # (FORMER BUG - solved by improving update_isolation to make sure values for isolation are kept up to date)
    secondary_cases_list %>% count_prop(isolate_after_ct_symp, value(contact_tracing_test_timing), value(isolate_ct_test_timing)) %$%
      sum(isolate_after_ct_symp & value(contact_tracing_test_timing) & !value(isolate_ct_test_timing)) # should be 0
    
    
    # (11) Check contact tracing consistency
    
        # - Only isolate if you are tested
        sum(
          is.na(secondary_cases_list$contact_tracing_test_timing) &
            !is.na(secondary_cases_list$isolate_ct_test_timing)
        ) # should be 0
        
        # Only isolate if you have results
        sum(
          is.na(secondary_cases_list$contact_tracing_results_timing) &
            !is.na(secondary_cases_list$isolate_ct_results_timing)
        ) # should be 0
    
        
        
    # (12) Check people don't isolate if they have a false negative
        
        # live_cases_list %$% sum(ct_test_negative & value(isolate_ct_test_timing)) 
        # actually this can be above 0 - people isolate just from CT not based on results
        
        secondary_cases_list %$% sum(ct_false_negative & value(isolate_ct_results_timing)) # should be 0
        
        secondary_cases_list %$% sum(random_test_false_negative & value(isolate_random_test_timing)) # should be 0
        
        secondary_cases_list %$% sum(self_test_false_negative & value(isolate_self_test_timing)) # should be 0
        
        # Previous bug - involving people being randomly tested twice
        secondary_cases_list %>%  
          group_by(case_id) %>% 
          select(case_id, contains("random_test")) %>% 
          view_filter(sum(random_test_timing < 0 & random_test_yn & random_test_yn_ever) > 0, n = 1) 
        # BUG : some people for whom random test is finished (random_test_timing < 0) don't have random_test_yn
        # this was people who are randomly tested in the same period 
        # SOLVED by changing previously_tested_positive condition to previously_tested
        
    
        
    # (13) How many people isolate presymptomatically?
    secondary_cases_list %>% filter(value(isolate_ct_test_timing)) %>% 
      count_prop(isolate_ct_test_timing < symptom_timing) # mostly pre-symptom
    
    secondary_cases_list %>% filter(value(isolate_ct_results_timing)) %>% 
      count_prop(isolate_ct_results_timing < symptom_timing) # about half and half
    
    secondary_cases_list %>% filter(value(isolate_random_test_timing)) %>% 
      count_prop(isolate_random_test_timing < symptom_timing) # VERY often presymptom!
    
    # people can never isolate presymptoms for self test 
    secondary_cases_list %$% sum(value(isolate_self_test_timing) & isolate_self_test_timing < symptom_timing) # NEVER pre-symptoms - should be 0

    
    # (14) Make sure secondary case group is always the same for every row secondary_case_id
    secondary_cases_list %>% select(secondary_case_id, secondary_i_group) %>% 
      dups_drop() %>% 
      dups_report(secondary_case_id) # should be no duplicates
    
  
    # VIEWING PORTAL
    secondary_cases_list %>% group_by(case_id, secondary_case_id) %>% arrange(case_id, t) %>%
      filter(!is.na(secondary_case_id)) %>% 
      relocate(case_id, t, secondary_case_id, secondary_case_timing, isolate_after_ct, deisolate_timing,
               contact_tracing_test_timing, contact_tracing_results_timing, isolation_timing, ct_false_negative, ct_test_negative, starts_with("isolate_"),
               isolation_timing_2) %>% 
      filter(value(deisolate_timing)) %>% 
      filter(case_id == "29_2_3_18_2_3_7") %>% view()
      # ungroup %>% count_prop(ct_false_negative & value(contact_tracing_results_timing))
      # filter(sum(value(deisolate_timing)) > 0) %>%
      # filter(isolate_after_ct & (ct_false_negative | ct_test_negative) & value(contact_tracing_results_timing)) %>% 
      # filter((ct_false_negative | ct_test_negative), value(contact_tracing_results_timing)) %>%
      view_n(1)
      
      secondary_cases_list$secondary_case_group %>% typeof()
      
      secondary_cases_list %>% group_by(case_id) %>% 
        filter(!is.na(secondary_case_id)) %>% 
        arrange(case_id, t) %>%
        view_n(1)
    
  
        
      secondary_cases_list %>% count_nas()
    
    
# 3. Live/secondary matching consistency checks ---------------------------
    
    
    
    # (2) No duplicates
    cases_merged %>% dups_report(t, case_id, secondary_case_id) 
    
    # (3) Check no duplicates for the x_only
    cases_merged %>% 
      filter(.merge == "x_only") %>% 
      dups_report(t, case_id) # there should be only one row per obs for the people only in live_cases
    
    # (4) Check whether the _L and _S variables are the same when they exist in both datasets
    

    
    # Check there are no conflicts between values across datasets
    variable_match %>% ungroup %>% summarise(across(ends_with("_conflict"), list(total = sum, prop = mean))) %>% 
      pivot_longer(everything(), names_to = c("variable", ".value"), names_pattern = "(.*)_(total|prop)") %>% 
      print_all   # there should be no conflicting variables
    
    
    # Previously, there were conflicts for contact tracing timings...
    # This is due to: the first period in which contact_tracing is "added" to a case, 
    #  it only gets added to live, although it subsequently gets added to secondary afterwards
    # Do they always get updated 1 period AFTER being infected? i.e. first period of infection it's not there, and then it's probably
    # update through the update_live_contact_tracing thing
    
    # SOLVED - by changing the order of (2) RANDOM TESTING in the outbreak_step function
    # so that secondary cases are updated AFTER live cases had their CT times updated
    
    # variable_match %>% 
    #   group_by(case_id) %>% 
    #   filter(sum(contact_tracing_test_timing_conflict > 0) | sum(contact_tracing_results_timing_conflict > 0)) %>% 
    #   select(sort(names(.))) %>% relocate(t, case_id, secondary_case_id, dataset) %>% 
    #   view_n(1)  
    # 
    # # Look at an individual case, and look at parent
    # ind_case <- cases_merged %>% 
    #   filter(case_id == "39_3_5") %>% 
    #   select(sort(names(.))) %>% relocate(t, case_id, secondary_case_id) %>% 
    #   view()
    # 
    # parent_case <- cases_merged %>% 
    #   select(sort(names(.))) %>% relocate(t, case_id, secondary_case_id) %>% 
    #   filter(case_id == "39_3") %>% 
    #   view()
    
    
    # (5) Why is contact_tracing stuff often added when negative? 
    # i.e. the "first" time we see a contact tracing timing is when it's negative?
    # - SOLVED - edited update_live_contact_tracing so that contact tracing timings are removed (set to NA) when they occur before infection time
    
    # Check whether timing conflict is only the people for whom it starts off negative  (FALSE)
    variable_match %>% 
      group_by(case_id) %>% 
      arrange(case_id, t) %>% 
      summarise(
        contact_tracing_test_timing = first_non_na(contact_tracing_test_timing),
        contact_tracing_results_timing = first_non_na(contact_tracing_results_timing),
        conflict = sum(contact_tracing_test_timing_conflict > 0) | sum(contact_tracing_results_timing_conflict > 0)
      ) %>% 
      view()
    
    # Look at characteristics of the negative-starters
    variable_match %>% 
      group_by(case_id) %>% 
      arrange(case_id, t) %>% 
      mutate(
        contact_tracing_test_timing_first = first_non_na(contact_tracing_test_timing),
        contact_tracing_results_timing_first = first_non_na(contact_tracing_results_timing),
        conflict = sum(contact_tracing_test_timing_conflict > 0) | sum(contact_tracing_results_timing_conflict > 0)
      ) %>% 
      filter(contact_tracing_test_timing_first < 0) %>% 
      view()
    # They are people who get tested before they are infected
    
    
    # Look at one case, and look at when parent is "detected"
    variable_match %>% 
      group_by(case_id) %>% 
      arrange(case_id, t) %>% 
      mutate(
        contact_tracing_test_timing_first = first_non_na(contact_tracing_test_timing),
        contact_tracing_results_timing_first = first_non_na(contact_tracing_results_timing),
        conflict = sum(contact_tracing_test_timing_conflict > 0) | sum(contact_tracing_results_timing_conflict > 0)
      ) %>% 
      filter(contact_tracing_test_timing_first < 0) %>% 
      view()
    
    
    
    # (6) 
    # cases_merged %>% names %>% enframe %>% print_all
    # cases_merged %>% group_by(case_id, secondary_case_id) %>% arrange(case_id, t) %>% 
    #   relocate(case_id, t, secondary_case_id, secondary_case_timing) %>% 
    #   filter(sum(value(contact_tracing_results_timing_L)) > 0) %>% 
    #   view_n(1)
    #   view_filter(case_id == "30_6_18", n = 3)
    # view_n(1)
      
      
    cases_merged %>% 
      count_prop(ct_test_negative_L)
    
      
    cases_merged %>% 
      mutate(asymptomatic_transmission = secondary_case_timing < symptom_timing_L) %>% 
      count_prop(asymptomatic_transmission, is.na(secondary_case_id))
      relocate(case_id, t, secondary_case_id, secondary_case_timing, symptom_timing_L, asymptomatic_transmission) %>% 
      group_by(case_id, secondary_case_id) %>% 
      view_filter(is.na(asymptomatic_transmission), n = 10)
    
      
      
      
      

# 4. HH quarantining consistency checks -----------------------------------

    # (1) if in the same household, primary_quarantine should always be the same as secondary_quarnatined
    secondary_cases_list %$% 
      sum(hh_id == secondary_hh_id & primary_quarantined != secondary_quarantined, na.rm = TRUE) # should be 0
      
    # (2) All NAs for new_potential_case, susceptible, new_actual_case, pipped should come from NEW ENTRIES
    secondary_cases_list %>% select(new_potential_case, susceptible, new_actual_case, pipped) %>% count_nas()
    new_entries <- secondary_cases_list %>% group_by(case_id, secondary_case_id) %>% filter(t == min(t))
    new_entries %>% select(new_potential_case, susceptible, new_actual_case, pipped) %>% count_nas()
      
    
    
    primary_quarantine_check <- function(i) {
      
      print(paste0("rep ", i))
      
      # Randomly select a case_id / secondary_id combination that's primary quarantined
      random_primary_quarantined <- secondary_cases_list %>% filter(primary_quarantined) %>% 
        mutate(secondary_case_t = t + secondary_case_timing, .after = secondary_case_id) %>% 
        select(case_id, secondary_case_id, secondary_case_t, primary_quarantined, hh_id) %>% 
        dups_drop(warn = FALSE) %>% 
        slice_sample(n = 1) #%>% 
      # {print(.$secondary_case_t); invisible(.)}
      
      # Look at the corresponding hh status data to check interval
      random_primary_hh <- hh_status_list %>% 
        semi_join(random_primary_quarantined, by = "hh_id") %>% 
        mutate(across(matches("isolate_interval"),
                      ~ t + .x)) %>% 
        select(matches("isolate_interval")) %>% 
        dups_drop(warn = FALSE) %>% 
        janitor::remove_empty(which = c("rows"))
      
      
      # Check whether the value's in the interval
      check_tf <- any(data.table::between(random_primary_quarantined$secondary_case_t, random_primary_hh$isolate_interval_1_l, random_primary_hh$isolate_interval_1_r)) | 
        any(between(random_primary_quarantined$secondary_case_t, random_primary_hh$isolate_interval_2_l, random_primary_hh$isolate_interval_2_r))
      
      return(check_tf)
      
    }
    
    # Check 100 times
    tibble(i = 1:100) %>% 
      rowwise() %>% 
      mutate(primary_quarantine_check = primary_quarantine_check(i)) %>%
      count_prop(primary_quarantine_check)
    
    
    
    # Secondary
    secondary_quarantine_check <- function(i) {
      
      print(paste0("rep ", i))
      
      # Randomly select a case_id / secondary_id combination that's secondary quarantined
      random_secondary_quarantined <- secondary_cases_list %>% filter(secondary_quarantined) %>% 
        mutate(secondary_case_t = t + secondary_case_timing, .after = secondary_case_id) %>% 
        select(case_id, secondary_case_id, secondary_case_t, secondary_quarantined, secondary_hh_id) %>% 
        dups_drop(warn = FALSE) %>% 
        slice_sample(n = 1) #%>% 
      # {print(.$secondary_case_t); invisible(.)}
      
      # Look at the corresponding hh status data to check interval
      random_secondary_hh <- hh_status_list %>% 
        semi_join(random_secondary_quarantined, by = c("hh_id" = "secondary_hh_id")) %>% 
        mutate(across(matches("isolate_interval"),
                      ~ t + .x)) %>% 
        select(matches("isolate_interval")) %>% 
        dups_drop(warn = FALSE) %>% 
        janitor::remove_empty(which = c("rows"))
      
      
      # Check whether the value's in at least one of the intervals
      check_tf <- any(data.table::between(random_secondary_quarantined$secondary_case_t, random_secondary_hh$isolate_interval_1_l, random_secondary_hh$isolate_interval_1_r)) | 
        any(between(random_secondary_quarantined$secondary_case_t, random_secondary_hh$isolate_interval_2_l, random_secondary_hh$isolate_interval_2_r))
      
      return(check_tf)
      
    }
    
    # Check 100 times
    tibble(i = 1:100) %>% 
      rowwise() %>% 
      mutate(secondary_quarantine_check = secondary_quarantine_check(i)) %>%
      count_prop(secondary_quarantine_check)
    
    
    
    
      
    
      

# ............. -----------------------------------------------------------


# Timing tests -------------------------------------------------------------
    
    
    tic.clearlog()
    
    timing_test <- outbreak_sims(
      record_times = TRUE,
      n_sims = 1,
      n_iterations = 20,
      n_pop_total = 100000, 
      hh_size_data = data_save$hh_data_bogota %>% filter(i_group %in% 1:3),
      n_initial_cases = c(10, 50, 40), group_props = c(0.3, 0.5, 0.2), dt_approx = 1,
      recov_val = 10, # people recover 10 days after symptoms show (cuts off only ~ 1% of infections at the upper end of timing distribution)
      params_timing = data_save$params_timing,
      params_symptom_timing = data_save$params_symptom_timing,
      params_serial = data_save$params_serial,
      test_choice_delay_data = data_save$test_choice_delay_data %>% filter(i_group %in% 1:3), 
      test_results_delay_data = data_save$test_results_delay_data %>% filter(i_group %in% 1:3), 
      test_sensitivity_data = data_save$test_sensitivity_data, 
      ct_delay_data = ct_delay_data_real,
      probs_self_test_df = probs_df_basic %>% mutate(prob = c(probs_default * 0.1, probs_default * 0.5, probs_default * 1)),
      probs_isolate_symptoms_df = probs_df_basic %>% mutate(prob = c(probs_default * 0.2, probs_default * 0.4, probs_default * 0.8)),
      probs_isolate_test_df = probs_df_basic %>% mutate(prob = c(probs_default * 0.5, probs_default * 0.8, probs_default * 1)),
      probs_isolate_ct_df = probs_df_basic %>% mutate(prob = c(probs_default * 0.6, probs_default * 0.6, probs_default * 1)),
      p_contact_if_isolated_home = c(0.3, 0.2, 0.1), 
      p_contact_traced = c(0.1, 0.2, 0.2),
      p_hh_quarantine = c(0.4, 0.5, 0.5),
      k_matrix = k_matrix_trial,
      contact_dispersion = 0.2,
      sar_home = c(0.4, 0.3, 0.3), sar_out = c(0.2, 0.2, 0.2), 
      infectiousness_by_symptoms = c(0.7, 1, 1, 1, 1),
      alpha = c(0, 0, 0)
    )
    
    plot(timing_test, indiv_sims = TRUE)
    
    
    # Store the tictoc log
    tic_log_basic <- tic.log(format = FALSE)
    
    # Make timings into a readable format
    tic_log_df <- tibble(
      type = map_chr(tic_log_basic, "msg"),
      tic = map_dbl(tic_log_basic, "tic"),
      toc = map_dbl(tic_log_basic, "toc"),
      time_elapsed = toc - tic
    ) %>% 
      separate(type, into = c("type", "t")) %>% 
      mutate(t = as.integer(t)) %>% 
      arrange(type, t) %>% 
      group_by(type, t) %>% 
      select(-tic, -toc) %>% 
      mutate(sim = row_number()) %>% 
      summarise(time_elapsed = mean(time_elapsed))
    
    # SAVE CURRENT TO USE AS COMPARISON
    # tic_log_previous <- tic_log_df
    

    # Plot timings 
    tic_log_df %>% 
      filter(type %in% c("findpotentialcases", "newcases", "bindcases") | str_detect(type, "newcases")) %>%
      ggplot(aes(x = t, y = time_elapsed, colour = type)) + 
      coord_cartesian(ylim = c(0, NA)) + 
      geom_hline(yintercept = 0) + 
      geom_line()
    # facet_wrap(~ version)
    
    
    
    # Or plot comparison
    tic_log_comparison <- bind_rows(
      previous = tic_log_previous, new = tic_log_df,
      .id = "version"
    )
    
    tic_log_comparison %>% 
      # filter(type %in% c("findpotentialcases", "newcases", "bindcases")) %>%
      ggplot(aes(x = t, y = time_elapsed, colour = type)) + 
      coord_cartesian(ylim = c(0, NA)) + 
      geom_hline(yintercept = 0) + 
      geom_line() + 
      facet_wrap(~ version)
    
    tic_log_previous$time_elapsed %>% sum() # used to be 1024! hooray
    tic_log_df$time_elapsed %>% sum()       # now is 144!
    
    # Slowness is driven by draw_secondary_cases - test this individually
    primary_cases_for_timing <- draw_symptoms_recovery(ids = 1:100,
                                                       i_group = rep(1:2, 50) ,
                                                       hh_id = NULL,
                                                       hh_ind_id = NULL,
                                                       hh_size = NULL,
                                                       would_quarantine = NULL,
                                                       marg_dist_primary = marg_dist_primary,
                                                       contact_tracing_test_timing = NULL,    
                                                       contact_tracing_results_timing = NULL,
                                                       test_sensitivity_data = tibble(days_exposure = as.integer(1:40), sensitivity = 0.8),
                                                       infection_timings = 0,              # inputs from the other draws
                                                       test_choice_delay_data = tibble(i_group = 1:2, test_choice_delay = as.double(1:2)),
                                                       test_results_delay_data = tibble(i_group = 1:2, test_results_delay = as.double(3:4)), recov_val = c(10, 10),
                                                       probs_self_test_df = probs_df_dummy, probs_isolate_symptoms_df = probs_df_dummy, 
                                                       probs_isolate_test_df = probs_df_dummy, probs_isolate_ct_df = probs_df_dummy,
                                                       infectiousness_by_symptoms = c(0.7, 1, 1, 1, 1))
    
    # %>% 
      draw_random_testing(alpha = c(0.1, 0.2, 0.3), dt_approx = 1) %>% 
      update_timing(dt = 1) %>% 
      draw_random_testing(alpha = c(0.1, 0.2, 0.3), dt_approx = 1) %>% 
      update_timing(dt = 1) %>% 
      draw_random_testing(alpha = c(0.1, 0.2, 0.3), dt_approx = 1) 
    
    
    tic()
    
    primary_cases_for_timing %>% draw_secondary_cases(p_contact_if_isolated_home = c(0.05, 0.1, 0.2), p_contact_traced = c(0.1, 0.2, 0.3),
                                                      beta_matrix = crossing(to = 1:3, from = 1:3) %>% mutate(beta_val = c(10, 5, 2, 5, 7, 1, 2, 1, 4)),
                                                      group_props = c(0.25, 0.5, 0.25), 
                                                      r0 = c(1, 2, 3),
                                                      dispersion = 0.16)
    
    toc() # was approx 2.5-3 sec, now 1.15 sec
    
    
    
    
    
    
    # SPEED UP SECONDDARY_ID_OUT
    hh_data <- timing_test$outbreak_t[[1]]$hh_status
    hh_ind_id <- sample(hh_data$hh_ind_id, size = 10000, replace = FALSE)
    hh_id <- tibble(hh_ind_id = hh_ind_id) %>% left_join(hh_data) %>% .$hh_id
    n_secondary_cases_home <- tibble(hh_ind_id = hh_ind_id) %>% left_join(hh_data) %>%
      mutate(n_secondary_cases_home = if_else(hh_size > 3, hh_size - 2, hh_size - 1)) %>%
      .$n_secondary_cases_home
    i_group <- tibble(hh_ind_id = hh_ind_id) %>% left_join(hh_data) %>% .$i_group
    n_secondary_cases_out <- rep(c(10, 5), 5000)
    beta_matrix <- crossing(
      to = 1:3, from = 1:3
    ) %>%
      mutate(beta_val = c(10, 5, 2, 5, 7, 1, 2, 1, 4))
    
    
    hh_data_weighted <- map(
      .x = unique(beta_matrix$from), 
      function(x) {
        left_join(hh_data, beta_matrix[beta_matrix$from == x, ], by = c("i_group" = "to"))
      }
    )
    
    
    hh_data_weighted_slim <- hh_data_weighted %>% 
      map(select, beta_val, hh_id, hh_ind_id)
    
    f1 <- function() {
      pmap(
        .l = list(n_secondary_cases_out, i_group, hh_id),
        function(x, y, z) {
          hh_data_w <- hh_data_weighted[[y]] # choose the hh_data_weighted df for that i_group
          probs <- hh_data_w$beta_val      # probs are set to values from beta matrix
          probs[hh_data_w$hh_id == z] <- 0 # set probability to 0 if hh_id is the same
          if (x == 0) NA_character_
          else sample(x = hh_data_w$hh_ind_id, size = x, prob = probs, replace = FALSE)
        }
      )
    }
    
    f1_slim <- function() {
      pmap(
        .l = list(n_secondary_cases_out, i_group, hh_id),
        function(x, y, z) {
          hh_data_w <- hh_data_weighted_slim[[y]] # choose the hh_data_weighted df for that i_group
          probs <- hh_data_w$beta_val      # probs are set to values from beta matrix
          probs[hh_data_w$hh_id == z] <- 0 # set probability to 0 if hh_id is the same
          if (x == 0) NA_character_
          else sample(x = hh_data_w$hh_ind_id, size = x, prob = probs, replace = FALSE)
        }
      )
    }
    
    f2 <- function() {
      l <- list()
      for (i in seq_along(n_secondary_cases_out)) {
        # i <- 3
        hh_data_w <- hh_data_weighted[[ i_group[[i]] ]] # choose the hh_data_weighted df for that i_group
        probs <- hh_data_w$beta_val      # probs are set to values from beta matrix
        probs[hh_data_w$hh_id == hh_id[[i]] ] <- 0 # set probability to 0 if hh_id is the same
        if (n_secondary_cases_out[[i]] == 0) l[[i]] <- NA_character_
        else l[[i]] <- sample(x = hh_data_w$hh_ind_id, size = n_secondary_cases_out[[i]], prob = probs, replace = FALSE)
      }
    }
    
    microbenchmark::microbenchmark(
      f1(), f1_slim(), f2(),
      times = 5
    ) # starting point - 
    
    
    # OUTCOME - very little difference between them.. stick to f1_slim()
    
    
    
    
    
    
    
    
    
    
# Try calculating Re ------------------------------------------------------
    
    all_case_ids <- secondary_cases_list %>% 
      mutate(infection_t = t + infection_timing) %>% 
      select(case_id, potential_secondary_cases, i_group, infection_t) %>% 
      dups_drop()
    
    ind_case_counter <- consistency_check_data$outbreak_t_record[[1]] %>% last() %>% .$ind_case_counter %>% 
      group_by(case_id) %>% 
      summarise(new_potential_case = sum(new_potential_case),
                new_actual_case = sum(new_actual_case), .groups = "drop")
    
    all_case_counter <- left_join(all_case_ids, ind_case_counter, by = "case_id") %>% 
      mutate(
        across(c(new_potential_case, new_actual_case),
               ~ if_else(is.na(.), 0L, .))
      ) %>% 
      mutate(unfinished_infection_cycle = new_potential_case < potential_secondary_cases) %>% 
      count_prop(unfinished_infection_cycle)
    
    re_estimate <- all_case_counter #%>% 
    # filter(!unfinished_infection_cycle)
    
    re_estimate %>% ggplot(aes(x = infection_t)) + 
      geom_smooth(aes(y = potential_secondary_cases), colour = "indianred") + 
      geom_smooth(aes(y = new_actual_case), colour = "darkgreen")
    
    re_estimate %>% summarise(potential_secondary = mean(potential_secondary_cases),
                              actual_secondary = mean(new_actual_case))
    
    # At the moment, the gap between potential and actual is due to (1) immnunity, and (2) control measures
    # It would be good to be able to isolate the difference just to control measures
    # E.g. possibly track the "reason" that case doesn't become actual.
    
    
    
    

# Calculate Re v2 ---------------------------------------------------------

    
    # Get all live cases by infection time (for denominator of Re calc)
    live_cases_by_infection_time <- live_cases_list %>% 
      mutate(infection_t = t + infection_timing) %>% 
      select(case_id, infection_t, i_group) %>% 
      dups_drop() %>% 
      mutate(infection_day = infection_t %/% 1)
    

    
    # Get each secondary infection and outcome
    ind_case_counter_list <- consistency_check_data$outbreak_t_record[[1]] %>% map("ind_case_counter") %>% 
      set_names(t_vec) %>% 
      bind_rows(.id = "t") %>% 
      mutate(t = as.numeric(t),
             infection_t = t + infection_timing) %>% 
      select(-t, -infection_timing) %>% 
      relocate(infection_t)
    
    dups_report(ind_case_counter_list, case_id, secondary_case_id) # no dups
    
    
    # Get OVERALL Re first (not by i_group)
    n_secondary_by_week <- ind_case_counter_list %>% 
      # sample_n(20) %>%
      mutate(infection_day = infection_t %/% 1) %>% 
      group_by(infection_day, case_type) %>% 
      summarise(n_secondary = n())
    
    n_live_by_week <- live_cases_by_infection_time %>% 
      group_by(infection_day) %>% 
      summarise(n_live = n())
      
    
    
    
    re_calc <- left_join(
      n_secondary_by_week, n_live_by_week, by = "infection_day"
    ) %>% 
      filter(case_type != "immune") %>% 
      group_by(case_type, infection_day) %>% arrange(case_type, infection_day) %>% 
      mutate(re = n_secondary / n_live,
             across(where(is.numeric), as.double))
    
    weighted_re_vec_actual <- re_calc %>% 
      filter(case_type == "actual_case") %>% 
      rollapply(
        width = 7, fill = NA_real_, align = "right",
        function(z){
          #uncomment if you want to see the structure
          return(
            weighted_mean = weighted.mean(x = as.double(z[, "re"]), w = as.double(z[, "n_live"]))
          )
        },
        by.column = FALSE
      )
    
    weighted_re_vec_isolated <- re_calc %>% 
      filter(case_type == "isolated_no_contact") %>% 
      rollapply(
        width = 7, fill = NA_real_, align = "right",
        function(z){
          #uncomment if you want to see the structure
          return(
            weighted_mean = weighted.mean(x = as.double(z[, "re"]), w = as.double(z[, "n_live"]))
          )
        },
        by.column = FALSE
      )
    
    re_calc_with_roll <- re_calc %>% 
      ungroup %>% 
      mutate(
        re_roll = c(weighted_re_vec_actual, weighted_re_vec_isolated)
      )
    
    
    max_second_case_timing <- secondary_cases_list %>% .$secondary_case_timing %>% max(na.rm = TRUE)
    
    # re_calc_with_roll %>% view()
    
    # PLOT 
    re_calc_with_roll %>% 
      complete(infection_day, case_type) %>% 
      group_by(infection_day) %>% 
      filter(sum(is.na(re_roll)) == 0) %>% 
      ungroup %>% 
      filter(infection_day < max(infection_day) - max_second_case_timing) %>%  # conservative cut off point to make sure we're not underestimtaing at the end of simulation
      mutate(case_type = fct_rev(fct_recode(case_type, 
                                           "no isolation" = "isolated_no_contact",
                                           "with_isolation" = "actual_case"))) %>% 
               ggplot(aes(x = infection_day, y = re_roll, colour = case_type)) + 
               geom_line(position = "stack", size = 1.2) + 
      geom_hline(yintercept = 1, linetype = "dashed", colour = "darkgrey") + 
      labs(y = expression(R[e]), x = "Day infected", colour = element_blank())
    
      # mutate(
      #   case_type = if_else(case_type == "isolated_no_contact"),sk
      # )
      # mutate(n_secondary_no_isolation = sum(n_secondary)) %>% 
      # filter(case_type != "isolated_no_contact") %>% 
      # select(-case_type) %>% 
      # mutate(re_no_isolation = n_secondary_no_isolation / n_live) %>% 
      # relocate(infection_day, n_secondary, n_secondary_no_isolation, n_live) %>% 
      # pivot_longer(
      #   -c(infection_day, n_live),
      #   names_to = c("n_secondary", "re"),
      #   names_pattern = "(.*secondary.*)(_no_isolation)"
      # )
    
      
      
    
    # ind_case_counter_list %>% view_n()
    

# Correlation between infectiousness and infectious period ----------------

    secondary_cases_list %>% 
      mutate(
        infectious_period = recovery_timing - infection_timing
      ) %>% 
      select(case_id, infectious_period, potential_secondary_cases) %>% 
      dups_drop(case_id) %>% 
      
      ggplot(aes(x = infectious_period, y = potential_secondary_cases)) + 
      geom_smooth()
    
    
    
    


# Correlations on timing variables ----------------------------------------

    k_matrix_trial <- Matrix::forceSymmetric(
      matrix(
        c(20, 1, 1, 
          4,  15, 1,
          1,  2, 10),
        nrow = 3
      )
    ) %>% as.matrix()
    
    
    
    outbreak_setup_test <- outbreak_setup(
      n_pop_total = 10000, 
      n_initial_cases = c(2, 3, 4, 1), group_props = data_save$group_props, dt_approx = 1,
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
      p_contact_traced = data_save$params_data$p_contact_traced,
      p_hh_quarantine = data_save$params_data$p_hh_quarantine,
      k_matrix = data_save$k_matrix * 1.3,
      contact_dispersion = data_save$contact_dispersion,
      sar_home = data_save$params_data$sar_home_input, 
      sar_out = data_save$params_data$sar_out_input, 
      infectiousness_by_symptoms = c(0.7, 1),
      # alpha = c(0.1, 0.1, 0.1)
      alpha = c(0, 0, 0, 0)
    )
    
    
    # CHECK PROPORTION of infections from each group
    outbreak_setup_test$secondary_cases %>% print_names %>% 
      filter(infection_type == "out") %>% 
      select(i_group, secondary_i_group) %>% 
      group_by(i_group, secondary_i_group) %>% 
      summarise(n = n()) %>% 
      summarise(p = n / sum(n))
    
    k_matrix_trial / rowSums(k_matrix_trial)
    
    
    dist_data <- outbreak_setup_test$secondary_cases
    


    ggplot(dist_data, aes(x = symptom_timing, y = serial_interval)) + 
      geom_point(alpha = 0.05, size = 0.3) + 
      stat_density_2d(aes(fill = ..level..), geom = "polygon")
    
    ggplot(dist_data, aes(x = secondary_symptoms_delay, y = serial_interval)) + 
      geom_point(alpha = 0.1, size = 0.3) + 
      stat_density_2d(aes(fill = ..level..), geom = "polygon") 
    
    ggplot(dist_data, aes(x = symptom_timing, y = secondary_symptoms_delay)) + 
      geom_point(alpha = 0.05, size = 0.3) + 
      stat_density_2d(aes(fill = ..level..), geom = "polygon")
    

    ggplot(dist_data, aes(x = symptom_timing, y = secondary_case_timing)) + 
      geom_abline(slope = 1, intercept = 0, linetype = "dashed", colour = "indianred") + 
      geom_point(alpha = 0.05, size = 0.3) + 
      coord_cartesian(xlim = c(0, 20), ylim = c(0, 20)) + 
      stat_density_2d(aes(fill = ..level..), geom = "polygon") + 
      geom_smooth(colour = "darkgreen", method = "lm")
    
    ggplot(dist_data, aes(x = secondary_symptoms_delay, y = secondary_case_timing)) + 
      geom_point(alpha = 0.05) + 
      stat_density_2d(aes(fill = ..level..), geom = "polygon")
    
    
    ggplot(dist_data, aes(x = symptom_timing, y = secondary_case_timing - symptom_timing)) + 
      geom_point(alpha = 0.1, size = 0.3) + 
      stat_density_2d(aes(fill = ..level..), geom = "polygon")
    
    
    dist_data %>% hist_basic(x = secondary_case_timing - symptom_timing)
    dist_data %$% max(secondary_case_timing - symptom_timing) # never above 7
    # Count proportion above 15
    dist_data %>% count_prop(secondary_case_timing - symptom_timing >= 10)
    
    # library("data.table")
    # ggplot(dist_data %>% filter(primary_symptoms %between% c(10, 11)), 
    #        aes(x = secondary_symptoms_delay, y = secondary_timing)) + 
    #   coord_cartesian(xlim = c(0, 20), ylim = c(0, 20)) + 
    #   geom_point(alpha = 0.1) + 
    #   stat_density_2d(aes(fill = ..level..), geom = "polygon")
    
    
    # secondary_timing_data <- rweibull(n = 100000,
    #                                   scale = 6 / (log(2))^(1/2.8),
    #                                   shape = 2.8)
    
    
    
    
    
    
    
    # PLOT COMPARISONS OF DISTRIBUTIONS
    gen_data_clean %>% hist_basic(x = secondary_timing - primary_symptoms)
    gen_data_clean %>% hist_basic(x = secondary_timing)
    
    
    
    compare_dist <- function(data, generated) {
      ggplot() + 
        geom_density(data = tibble(x = data), aes(x = x), colour = "indianred", linetype = "dashed", show.legend = TRUE) + 
        geom_density(data = tibble(x = generated), aes(x = x), colour = "skyblue")
    }
    
    secondary_timing_data <- rgamma(n = 10000,  
                                    shape = params_test[[8]],
                                    scale = 6/params_test[[8]])
    
    
    compare_dist(data = secondary_timing_data,
                 generated = gen_data_clean$secondary_timing)
    
    
    
    
    compare_dist(data = rlnorm(10000, meanlog = 1.63, sdlog = 0.5),
                 generated = gen_data_clean$primary_symptoms)   
    
    compare_dist(data = rgamma(10000, shape = serial_parameters$shape, rate = serial_parameters$rate) + min.serial - 0.5,
                 generated = gen_data_clean$serial_interval)   
    
    
    
    
        
    
    

    
    


# ............ ------------------------------------------------------------


# OLD - group differences---------------------------------------------------------------------


    
    
    # Calculate the differences each parameter set make
    group_diff_quantiles <- group_diff_df %>%
      .[names(.) != "all"] %>%
      map("time_series") %>%
      bind_rows(.id = "parameter_set") %>%
      group_by(parameter_set, i_group, sim_id) %>% summarise(n_cases_cum = max(n_cases_cum)) %>%
      group_by(parameter_set, sim_id) %>%
      mutate(n_cases_cum_diff = n_cases_cum[i_group == 1] - n_cases_cum[i_group == 2]) %>%
      group_by(parameter_set) %>%
      summarise(quibble(n_cases_cum_diff, q = c(0.1, 0.25, 0.5, 0.75, 0.9)),
                quibble(n_cases_cum, q = c(0.1, 0.25, 0.5, 0.75, 0.9)))
    
    
    
    
    # Plot the contribution of each factor (at the medians)
    total_diff_median <- group_diff_all_summ$n_cases_cum_diff[group_diff_all_summ$q == 0.5]
    total_cases_median <- group_diff_all_summ$n_cases_cum[group_diff_all_summ$q == 0.5]
    
    # THIS NEEDS TO BE UPDATED WHEN WE HAVE REAL RESULTS
    group_diff_effect <- group_diff_quantiles %>% 
      filter(parameter_set != "baseline") %>% 
     
      ungroup %>% 
      mutate(effect = total_diff_median - n_cases_cum_diff, 
             effect_on_total_cases = total_cases_median - n_cases_cum)
    
    
    group_diff_effect %>% 
      mutate(type = case_when(q == 0.5 ~ "median",
                              q == 0.1 ~ "upper",
                              q == 0.9 ~ "lower")) %>% 
      drop_na() %>% 
      select(parameter_set, effect, type) %>% 
      pivot_wider(names_from = type, values_from = effect) %>% 
      ggplot(aes(x = parameter_set, fill = parameter_set)) + 
      geom_col(aes(y = median), width = 0.8) + 
      geom_errorbar(aes(ymax = upper, ymin = lower), width = 0.2) + 
      coord_flip() + 
      labs(y = "Reduction in rich-poor gap in total cases", x = element_blank())
    
    group_diff_effect %>% 
      mutate(type = case_when(q == 0.5 ~ "median",
                              q == 0.1 ~ "upper",
                              q == 0.9 ~ "lower")) %>% 
      drop_na() %>% 
      select(parameter_set, effect_on_total_cases, type) %>% 
      pivot_wider(names_from = type, values_from = effect_on_total_cases) %>% 
      ggplot(aes(x = parameter_set, fill = parameter_set)) + 
      geom_col(aes(y = median), width = 0.8) + 
      geom_errorbar(aes(ymax = upper, ymin = lower), width = 0.2) + 
      coord_flip() + 
      labs(y = "Reduction in total cases", x = element_blank())
    
    
    
      
    
    # PERCENTAGE CONTRIBUTION - manually edit for presentation purposes
    group_diff_contrib <- group_diff_effect %>% 
      filter(q == 0.5) %>% # only look at the median
    # bind_rows(tibble(parameter_set = "unidentified", n_cases_cum_diff = total_diff_median - sum(.$n_cases_cum_diff))) %>%
      # MANUALLY SET TO 100 for illustration purposes
      mutate(effect = if_else(effect < 0, -effect, effect)) %>% 
      bind_rows(tibble(parameter_set = "interaction effects", effect = 100)) %>%
      mutate(
        prop_contrib = effect / sum(effect),
        perc_contrib = paste0(round(prop_contrib, 2) * 100, "%")
      ) %>% 
      mutate(
        parameter_set = case_when(
          parameter_set == "r0" ~ "contact rates",
          parameter_set == "healthcare" ~ "health care / testing delays", 
          TRUE ~ parameter_set
        ),
        parameter_set = fct_rev(factor(parameter_set, 
                                       levels = c("contact rates",
                                                  "isolation_behaviour",
                                                  "testing_decision",
                                                  "health care / testing delays",
                                                  "interaction effects")))
      )
             
    # library("ggrepel")
    
    
      
      
      
    # PERCENTAGE CONTRIBUTION
    group_diff_contrib %>% 
      ggplot(aes(x = TRUE)) + 
      # geom_col(aes(y = total_diff_median), position = "identity") +
      geom_col(aes(y = effect, fill = parameter_set), colour = "white", position = "stack", width = 0.3) +
      geom_label_repel(
        aes(y = effect, group = parameter_set, label = perc_contrib), 
        position = position_stack(vjust = 0.5), size = 2.5, label.padding = 0.35
      ) + 
      theme_custom() + 
      scale_fill_manual(values = c("lightgrey", "#E17C6F", "#5BAFF3", "#D175F0", "#64BC80")) +
      scale_x_discrete(labels = element_blank()) + 
      # geom_bar(aes(y = total_diff_median), stat = "identity") + 
      coord_flip() + 
      labs(x = element_blank(), y = "Differences in total cases (group 2 - group 1)")
    
    
    
    ?geom_label_repel
    
    
    
    # --- OLD 
    
    
    # JUST R0 / transmission matrix
    group_diff_r0 <- outbreak_sims(
      n_sims = 100,
      n_iterations = 100,
      n_pop = 10000, 
      n_initial_cases = c(100, 100), group_props = c(0.5, 0.5), dt_approx = 1,
      test_results_delay_data = 1.5, recov_val = 7,
      test_sensitivity = 0.85, 
      probs_self_test_df =        probs_df_basic_2 %>% mutate(prob = rep(c(0, 0.4, 0.7, 0.8, 1), 2)),
      probs_isolate_symptoms_df = probs_df_basic_2 %>% mutate(prob = rep(c(0, 0.4, 0.7, 0.8, 1), 2)),
      probs_isolate_test_df =     probs_df_basic_2 %>% mutate(prob = rep(c(0.6, 0.7, 0.8, 1, 1), 2)),
      probs_isolate_ct_df =       probs_df_basic_2 %>% mutate(prob = rep(c(0.3, 0.4, 0.6, 1, 1), 2)),
      p_contact_if_isolated_home = c(0.2, 0.2), 
      p_contact_traced =      c(0.3, 0.3),
      beta_matrix = crossing(to = 1:2, from = 1:2) %>% mutate(beta_val = c(8, 2, 2, 4)), # implies R0 = (3.5, 2.1)
      r0_group_1 = 3.5, 
      dispersion = 0.16, 
      alpha = c(0, 0)
    )
    
    # calculate_r0s(beta_matrix = crossing(to = 1:2, from = 1:2) %>% mutate(beta_val = c(8, 2, 2, 6)),
    #               group_props = c(0.5, 0.5),
    #               r0_group_1 = 3.5)
    
    
    plot(group_diff_r0, indiv_sims = FALSE, conf_level = 0.8)
    plot_by_i_group(group_diff_r0, conf_level = 0.9)
    plot_by_i_group_policymaker_mode(group_diff_r0, conf_level = 0.8)
    
    
    group_diff_r0$time_series %>% group_by(i_group, sim_id) %>% summarise(n_cases_cum = max(n_cases_cum)) %>% 
      ungroup %>% 
      ggplot(aes(x = n_cases_cum, fill = factor(i_group))) + 
      geom_histogram(boundary = 0, position = "identity", colour = "white")
      # facet_wrap(~ factor(i_group))
      # quantile_summarise(n_cases_cum, conf_level = 0.5)
    
    
  
    
    group_diff_r0$time_series %>% group_by(i_group, sim_id) %>% summarise(n_cases_cum = max(n_cases_cum)) %>% 
      group_by(sim_id) %>% 
      mutate(n_cases_cum_diff = n_cases_cum[i_group == 1] - n_cases_cum[i_group == 2]) %>% 
      ungroup %>% 
      summarise(quibble(n_cases_cum_diff, q = c(0.1, 0.25, 0.5, 0.75, 0.9)))
    
    
    
    
    # JUST SELF TESTING BEHAVIOUR 
    group_diff_testing_behaviour <- outbreak_sims(
      n_sims = 20,
      n_iterations = 100,
      n_pop = 10000, 
      n_initial_cases = c(100, 100), group_props = c(0.5, 0.5), dt_approx = 1,
      test_results_delay_data = c(1.5, 1.5), recov_val = 7,
      test_sensitivity = 0.85, 
      probs_self_test_df =        probs_df_basic_2 %>% mutate(prob = c(c(0, 0.2, 0.35, 0.4, 0.5), c(0, 0.4, 0.7, 0.8, 1))),
      probs_isolate_symptoms_df = probs_df_basic_2 %>% mutate(prob = rep(c(0, 0.4, 0.7, 0.8, 1), 2)),
      probs_isolate_test_df =     probs_df_basic_2 %>% mutate(prob = rep(c(0.6, 0.7, 0.8, 1, 1), 2)),
      probs_isolate_ct_df =       probs_df_basic_2 %>% mutate(prob = rep(c(0.3, 0.4, 0.6, 1, 1), 2)),
      p_contact_if_isolated_home = c(0.2, 0.2), 
      p_contact_traced =      c(0.3, 0.3),
      beta_matrix = crossing(to = 1:2, from = 1:2) %>% mutate(beta_val = c(4, 2, 2, 4)), # implies R0 = (3.5, 2.1)
      r0_group_1 = 3.5, 
      dispersion = 0.16, 
      alpha = c(0, 0)
    )
    
    group_diff_testing_behaviour$time_series %>% group_by(i_group, sim_id) %>% summarise(n_cases_cum = max(n_cases_cum)) %>% 
      ungroup %>% 
      ggplot(aes(x = n_cases_cum, fill = factor(i_group))) + 
      geom_histogram(boundary = 0, position = "identity", colour = "white")
    
    
    
    

# OLD - Work out how to rescale the beta/r0 -------------------------------------

    
    # Check the R0s
    calculate_r0s(beta_matrix = crossing(to = 1:2, from = 1:2) %>% mutate(beta_val = c(8, 2, 2, 6)),
                  group_props = c(0.5, 0.5), r0_group_1 = 3.5)
    
    calculate_r0s(beta_matrix = crossing(to = 1:2, from = 1:2) %>% mutate(beta_val = c(6, 2, 2, 6)),
                  group_props = c(0.5, 0.5), r0_group_1 = 2.8)
    
    # PROBLEM - at the moment I've made the counterfactual a case where not only R0 is lower,
    # but also group 1 is less assortative as well! How to change beta to reflect this?
    
    beta_matrix_test <- matrix(c(8, 4, 4, 6), nrow = 2)
    beta_eigen <- eigen(beta_matrix_test) # eigen values are 8 and 4
    beta_eigen$values[2] / beta_eigen$values[1]
    
    
    
    
    
    calc_assortative_alpha <- function(beta_matrix, group_props) {
      # beta_matrix <- beta_matrix_test; group_props <- c(0.5, 0.5)
      beta_ii <- beta_matrix[1, 1]
      beta_ij <- beta_matrix[1, 2]
      beta_jj <- beta_matrix[2, 2]
      p_i <- group_props[1]
      p_j <- group_props[2]
      
      num <- (beta_ii * beta_jj - beta_ij^2) * (p_i * p_j)
      den <- (beta_ii * p_i + beta_ij * p_j) * (beta_ij * p_i + beta_jj * p_j)
      alpha <- num/den
      return(alpha)
    }
    
    calc_contact_tau <- function(beta_matrix, group_props, group) {
      # beta_matrix <- beta_matrix_test; group_props <- c(0.5, 0.5); group <- 2
      other_group <- if_else(group == 1, 2, 1)
      beta_self <- beta_matrix[group, group]
      beta_other <- beta_matrix[group, other_group]
      p_self <- group_props[group]
      p_other <- group_props[other_group]
      
      tau <- (beta_self * p_self) + (beta_other * p_other)
      return(tau)
    }
    
    calc_assortative_alpha(matrix(c(10, 5, 5, 8), nrow = 2), group_props = c(0.5, 0.5))    
    calc_contact_tau(matrix(c(10, 2, 2, 5), nrow = 2), group_props = c(0.5, 0.5), group = 1)
    
    
    recover_beta <- function(beta_matrix_original, group_props) {
      alpha <- calc_assortative_alpha(beta_matrix_original, group_props = group_props)
      tau_i <- calc_contact_tau(beta_matrix_original, group_props = group_props, group = 1)
      tau_j <- calc_contact_tau(beta_matrix_original, group_props = group_props, group = 2)
      
      
      return(beta_matrix)
    }
    
    # recover_beta(matrix(c(10, 5, 5, 20), nrow = 2), group_props = c(0.5, 0.5))
    
    
    
    
    
    rescale_beta <- function(beta_matrix_original, group_props, scale_type, spread_val = NULL) {
      # beta_matrix_original <- matrix(c(8, 4, 4, 6), nrow = 2); group_props <- c(0.5, 0.5)
      alpha <- calc_assortative_alpha(beta_matrix_original, group_props = group_props)
      tau_i <- calc_contact_tau(beta_matrix_original, group_props = group_props, group = 1)
      tau_j <- calc_contact_tau(beta_matrix_original, group_props = group_props, group = 2)
      
      # CHECK it matches
      beta_ii_old <- (alpha * tau_i / group_props[1]) + ((1 - alpha) * tau_i^2 / (tau_i * group_props[1] + tau_j * group_props[2]))
      beta_ij_old <- (1 - alpha) * (tau_i * tau_j) / (tau_i * group_props[1] + tau_j * group_props[2])
      beta_jj_old <- (alpha * tau_j / group_props[2]) + (1 - alpha) * (tau_j^2 / (tau_i * group_props[1] + tau_j * group_props[2]))
      
      beta_matrix_recalc <- matrix(c(beta_ii_old, beta_ij_old, beta_ij_old, beta_jj_old), nrow = 2)
      if (!all.equal(beta_matrix_original, beta_matrix_recalc)) stop("something is going wrong with beta calculation")
      
      
      # Rescale group 1 to group 2
      if (scale_type == "rescale_1") {
        # Update group 1 to the same tau as group 2
        beta_ii_new <- (alpha * tau_j / group_props[1]) + ((1 - alpha) * tau_j^2 / (tau_j * group_props[1] + tau_j * group_props[2]))
        beta_ij_new <- (1 - alpha) * (tau_j * tau_j) / (tau_j * group_props[1] + tau_j * group_props[2])
        beta_jj_new <- (alpha * tau_j / group_props[2]) + (1 - alpha) * (tau_j^2 / (tau_j * group_props[1] + tau_j * group_props[2]))
        
      } else if (scale_type == "rescale_2") {
        # Update group 2 to the same tau as group 1
        beta_ii_new <- (alpha * tau_i / group_props[1]) + ((1 - alpha) * tau_i^2 / (tau_i * group_props[1] + tau_i * group_props[2]))
        beta_ij_new <- (1 - alpha) * (tau_i * tau_i) / (tau_i * group_props[1] + tau_i * group_props[2])
        beta_jj_new <- (alpha * tau_i / group_props[2]) + (1 - alpha) * (tau_i^2 / (tau_i * group_props[1] + tau_i * group_props[2]))
        
      } else if (scale_type == "spread") {
        if (is.null(spread_val)) stop("spread_type is spread, but spread_val is null")
        tau_j <- tau_j - spread_val
        tau_i <- tau_i + spread_val
      
        beta_ii_new <- (alpha * tau_i / group_props[1]) + ((1 - alpha) * tau_i^2 / (tau_i * group_props[1] + tau_j * group_props[2]))
        beta_ij_new <- (1 - alpha) * (tau_i * tau_j) / (tau_i * group_props[1] + tau_j * group_props[2])
        beta_jj_new <- (alpha * tau_j / group_props[2]) + (1 - alpha) * (tau_j^2 / (tau_i * group_props[1] + tau_j * group_props[2]))
        
      }
      
      new_beta_matrix <- matrix(c(beta_ii_new, beta_ij_new, beta_ij_new, beta_jj_new), nrow = 2)
    
      return(new_beta_matrix)
    }
    
    
    
    # WHY DOES RESCALE BETA INCREASE THE Bjj term?
    # intuition - when forming random partnerships, you are less likely to form a partnership with 
    # the OTHER group randomly (because their contact rate has dropped), so you are therefore slightly MORE
    # likely to form a partnership with your OWN group randomly.
    
    
    rescale_beta(matrix(c(10, 4, 4, 8), nrow = 2), group_props = c(0.5, 0.5),
                 scale_type = "spread", spread_val = 0)
    
    
    
    # 
    # rescale_r0 <- function(beta_matrix, group_props = c(0.5, 0.5)) {
    #   
    #   # beta_matrix <- matrix(c(10, 5, 5, 8), nrow = 2)
    #   
    #   
    #   rescaled_beta <- rescale_beta(beta_matrix, group_props = c(0.5, 0.5))
    #   
    #   
    #   p_conditional_rescaled <- combine_beta_group_props(beta_matrix = crossing(to = 1:2, from = 1:2) %>% mutate(beta_val = as.vector(rescaled_beta)),
    #                                                      group_props = c(0.5, 0.5)) %>% 
    #     group_by(from) %>% 
    #     mutate(p_conditional = beta_by_n / sum(beta_by_n)) %>% 
    #     arrange(from) %>% 
    #     summarise(total_beta_by_n = sum(beta_by_n), .groups = "drop") %>% 
    #     mutate(r0_ratio = total_beta_by_n / total_beta_by_n[from == 1])
    #   
    #   p_conditional_old <- combine_beta_group_props(beta_matrix = crossing(to = 1:2, from = 1:2) %>% mutate(beta_val = as.vector(beta_matrix)),
    #                                                 group_props = c(0.5, 0.5)) %>% 
    #     group_by(from) %>% 
    #     mutate(p_conditional = beta_by_n / sum(beta_by_n)) %>% 
    #     arrange(from) %>% 
    #     summarise(total_beta_by_n = sum(beta_by_n), .groups = "drop") %>% 
    #     mutate(r0_ratio = total_beta_by_n / total_beta_by_n[from == 1])
    #   
    #   
    # }
    
    
    
    
    # FIND SCALE from BETA TO R0
    # TEST: 
    
    rescale_r0 <- function(beta_matrix, r0_group_1, group_props, scale_type, spread_val = NULL) {
      # r0_group_1 <- 2; group_props = c(0.5, 0.5)
      # beta_matrix = crossing(to = 1:2, from = 1:2) %>% mutate(beta_val = as.vector(matrix(c(10, 5, 5, 8), nrow = 2)))
      
      beta_by_n <- combine_beta_group_props(beta_matrix = beta_matrix,
                                            group_props = c(0.5, 0.5)) %>% 
        group_by(from) %>% 
        mutate(p_conditional = beta_by_n / sum(beta_by_n)) %>% 
        arrange(from) %>% 
        summarise(total_beta_by_n = sum(beta_by_n), .groups = "drop") %>% 
        .$total_beta_by_n
      
      r0s <- calculate_r0s(beta_matrix = beta_matrix,
                           group_props = group_props,
                           r0_group_1 = r0_group_1)
      
      r0_1 <- r0s[1]; r0_2 <- r0s[2]
      
      r0_scaling_1 <- r0_1 / beta_by_n[1]
      r0_scaling_2 <- r0_2 / beta_by_n[2]
      
      if (!all.equal(r0_scaling_1, r0_scaling_2)) stop("r0_scaling factors should be the same across both groups")
      
      # r0_scaling_1
      
      rescaled_beta <- rescale_beta(matrix(beta_matrix$beta_val, nrow = 2), group_props = group_props,
                                    scale_type = scale_type, spread_val = spread_val)
      beta_matrix_new = crossing(to = 1:2, from = 1:2) %>% mutate(beta_val = as.vector(rescaled_beta))
      
      beta_by_n <- combine_beta_group_props(beta_matrix = beta_matrix_new,
                                            group_props = group_props) %>% 
        group_by(from) %>% 
        mutate(p_conditional = beta_by_n / sum(beta_by_n)) %>% 
        arrange(from) %>% 
        summarise(total_beta_by_n = sum(beta_by_n), .groups = "drop") %>%
        .$total_beta_by_n
      
      r0_new <- beta_by_n * r0_scaling_1
      
      return(r0_new)
      
    }
    
    
    # NEXT - redo rescale_r0 function to make sure it works with new rescale_beta
    
    
    # rescale_beta(beta_matrix_original = matrix(c(10, 5, 5, 8), nrow = 2),
    #              group_props = c(0.5, 0.5))
    
    
    
    
    rescale_r0(beta_matrix = crossing(to = 1:2, from = 1:2) %>% mutate(beta_val = as.vector(matrix(c(8, 5, 5, 8), nrow = 2))),
               r0_group_1 = 3.5, 
               group_props = c(0.5, 0.5),
               scale_type = "spread", spread_val = 0.5)
    
    calculate_r0s(
      beta_matrix = crossing(to = 1:2, from = 1:2) %>% mutate(beta_val = as.vector(matrix(c(10, 5, 5, 8), nrow = 2))),
      r0_group_1 = 3.5,
      group_props = c(0.5, 0.5)
    )

    
       
    


    

# Inequality is bad -------------------------------------------------------

    n_sims_inequality <- 5
    n_iterations_inequality <- 80
    probs_df_basic_2 <- crossing(i_group = 1:2, symptom_severity = (1L:5L))
    n_pop_inequality <- 10000
    contacts_baseline <- 50
    
    calc_contacts_mean(
      beta_matrix = crossing(to = 1:2, from = 1:2) %>% mutate(beta_val = c(6, 2, 2, 6)),
      sar_out = c(0.2, 0.2),
      group_props = c(0.5, 0.5),
      contacts_mean_scale = 15,
      scale_on = 1
    )
    
    
    
    group_symmetric <- outbreak_sims(
      n_sims = n_sims_inequality,
      n_iterations = n_iterations_inequality,
      n_pop = n_pop_inequality, 
      hh_size_data = list(
        `1` = rbinom(n_pop_inequality, size = 5, prob = 0.5) + 1,
        `2` = rbinom(n_pop_inequality, size = 5, prob = 0.5) + 1
      ),
      n_initial_cases = c(30, 30), group_props = c(0.5, 0.5), dt_approx = 1,
      test_results_delay_data = c(2, 2), recov_val = 7,
      test_sensitivity = 0.9, 
      probs_self_test_df =        probs_df_basic_2 %>% mutate(prob = c(c(0, 0.3, 0.525, 0.5, 0.75), c(0, 0.3, 0.525, 0.5, 0.75))),
      probs_isolate_symptoms_df = probs_df_basic_2 %>% mutate(prob = c(c(0, 0.3, 0.525, 0.5, 0.75), c(0, 0.3, 0.525, 0.5, 0.75))),
      probs_isolate_test_df =     probs_df_basic_2 %>% mutate(prob = c(c(0.45, 0.525, 0.6, 0.75, 0.75), c(0.45, 0.525, 0.6, 0.75, 0.75))),
      probs_isolate_ct_df =       probs_df_basic_2 %>% mutate(prob = c(c(0.225, 0.3, 0.45, 0.75, 0.75), c(0.225, 0.3, 0.45, 0.75, 0.75))),
      p_hh_quarantine = c(0.5, 0.5),
      p_contact_if_isolated_home = c(0.3, 0.3), 
      p_contact_traced =      c(0.2, 0.2),
      beta_matrix = crossing(to = 1:2, from = 1:2) %>% mutate(beta_val = c(6, 2, 2, 6)), 
      contacts_mean = c(contacts_baseline, contacts_baseline),
      sar_out = c(0.05, 0.05),
      secondary_timing_params = inf_profile_params,       sar_home = c(0.2, 0.2),
      alpha = c(0, 0),
      infectiousness_by_symptoms = c(0.7, 1, 1, 1, 1)
    )
    
    
    
    

    
    plot_by_i_group(group_symmetric)
    
    
    # Calculate rescaled betas & contact means
    
        beta_2 <- crossing(to = 1:2, from = 1:2) %>% mutate(beta_val = c(6, 2, 2, 6))
    
        # MILD DIFF
        beta_mild <- rescale_beta(
          beta_matrix = beta_2,
          rescale_factors = c(1.2, 0.778)
        )
        
        contacts_mild <- rescale_contacts_mean(
          beta_matrix = beta_2,
          rescale_factors = c(1.2, 0.778),
          group_props = c(0.5, 0.5),
          sar_out = c(0.05, 0.05),
          contacts_mean_scale = contacts_baseline,
          scale_on = 2
        ) #%>% 
          # {print(sum(.))}
      
        
        # EXTREME DIFF
        beta_extreme <- rescale_beta(
          beta_matrix = beta_2,
          rescale_factors = c(1.3, 0.646)
        )
        
        contacts_extreme <- rescale_contacts_mean(
          beta_matrix = beta_2,
          rescale_factors = c(1.3, 0.646),
          group_props = c(0.5, 0.5),
          sar_out = c(0.05, 0.05),
          contacts_mean_scale = contacts_baseline,
          scale_on = 1
        )
    

    

    group_mild_diff <- outbreak_sims(
      n_sims = n_sims_inequality,
      n_iterations = n_iterations_inequality,
      n_pop = n_pop_inequality,
      hh_size_data = list(
        `1` = rbinom(n_pop_inequality, size = 5, prob = 0.5) + 1,
        `2` = rbinom(n_pop_inequality, size = 5, prob = 0.5) + 1
      ),
      n_initial_cases = c(30, 30), group_props = c(0.5, 0.5), dt_approx = 1,
      test_results_delay_data = c(2.25, 1.75), recov_val = 7,
      test_sensitivity = 0.9, 
      probs_self_test_df =        probs_df_basic_2 %>% mutate(prob = c(c(0, 0.2, 0.35, 0.4, 0.5) * 1.25, c(0, 0.4, 0.7, 0.8, 1) * 7/8)),
      probs_isolate_symptoms_df = probs_df_basic_2 %>% mutate(prob = c(c(0, 0.2, 0.35, 0.4, 0.5) * 1.25, c(0, 0.4, 0.7, 0.8, 1) * 7/8)),
      probs_isolate_test_df =     probs_df_basic_2 %>% mutate(prob = c(c(0.3, 0.35, 0.4, 0.5, 0.5) * 1.25, c(0.6, 0.7, 0.8, 1, 1) * 7/8)),
      probs_isolate_ct_df =       probs_df_basic_2 %>% mutate(prob = c(c(0.15, 0.2, 0.3, 0.5, 0.5) * 1.25, c(0.3, 0.4, 0.6, 1, 1) * 7/8)),
      p_contact_if_isolated_home = c(0.35, 0.25), 
      p_contact_traced =      c(0.125, 0.275),
      p_hh_quarantine = c(0.4, 0.6),
      beta_matrix = beta_mild, 
      contacts_mean = contacts_mild,
      sar_out = c(0.05, 0.05),
      secondary_timing_params = inf_profile_params,       sar_home = c(0.2, 0.2),
      alpha = c(0, 0),
      infectiousness_by_symptoms = c(0.7, 1, 1, 1, 1)
    )
    
    
    group_extreme_diff <- outbreak_sims(
      n_sims = n_sims_inequality,
      n_iterations = n_iterations_inequality,
      n_pop = n_pop_inequality, 
      hh_size_data = list(
        `1` = rbinom(n_pop_inequality, size = 5, prob = 0.5) + 1,
        `2` = rbinom(n_pop_inequality, size = 5, prob = 0.5) + 1
      ),
      n_initial_cases = c(30, 30), group_props = c(0.5, 0.5), dt_approx = 1,
      test_results_delay_data = c(2.5, 1.5), recov_val = 7,
      test_sensitivity = 0.9, 
      probs_self_test_df =        probs_df_basic_2 %>% mutate(prob = c(c(0, 0.2, 0.35, 0.4, 0.5), c(0, 0.4, 0.7, 0.8, 1))),
      probs_isolate_symptoms_df = probs_df_basic_2 %>% mutate(prob = c(c(0, 0.2, 0.35, 0.4, 0.5), c(0, 0.4, 0.7, 0.8, 1))),
      probs_isolate_test_df =     probs_df_basic_2 %>% mutate(prob = c(c(0.3, 0.35, 0.4, 0.5, 0.5), c(0.6, 0.7, 0.8, 1, 1))),
      probs_isolate_ct_df =       probs_df_basic_2 %>% mutate(prob = c(c(0.15, 0.2, 0.3, 0.5, 0.5), c(0.3, 0.4, 0.6, 1, 1))),
      p_contact_if_isolated_home = c(0.4, 0.2), 
      p_contact_traced =      c(0.05, 0.35),
      p_hh_quarantine = c(0.3, 0.7),
      beta_matrix = beta_extreme, 
      contacts_mean = contacts_extreme,
      sar_out = c(0.05, 0.05),
      secondary_timing_params = inf_profile_params,       sar_home = c(0.2, 0.2),
      alpha = c(0, 0),
      infectiousness_by_symptoms = c(0.7, 1, 1, 1, 1)
    )
    
    
    plot_by_i_group(group_mild_diff)
    plot_by_i_group(group_extreme_diff)
    
    
    # PLOT ALL THE MODELS TOGETHER
    lst(
      group_symmetric,
      group_mild_diff,
      group_extreme_diff
    ) %>% 
      map("time_series") %>% 
      bind_rows(.id = "model") %>% 
      group_by(model, sim_id, t) %>% 
      summarise(across(starts_with("n_") & !where(is.list), sum), .groups = "drop_last") %>% 
      group_by(model, sim_id) %>%
      # mutate(new_cases = n_cases_cum - lag(n_cases_cum),
      #        new_cases_week = rollsum(new_cases, 7, align = "right", na.pad = TRUE)) %>%
      group_by(model, t) %>%
      select(t, sim_id, n_cases_live) %>%
      filter(!is.na(n_cases_live)) %>%
      quantile_summarise(n_cases_live, conf_level = 0.8) %>% 
      # summarise(
      #   across(n_cases_live,
      #          list(median = ~ median(.x),
      #               upper = ~ quantile(.x, conf_upper),
      #               lower = ~ quantile(.x, conf_lower)))
      # ) %>%
      ungroup %>%
      ggplot(aes(x = t, group = model)) +
      geom_line(aes(y = n_cases_live_median,
                    colour = model), size = 1.3) +
      geom_ribbon(aes(ymax = n_cases_live_upper, ymin = n_cases_live_lower,
                      fill = model), alpha = 0.3) +
      theme_custom() + theme(legend.position = "top")
   
    
    
    
    

# Effect of testing delays ------------------------------------------------

    
    
    # DO A PLOT TO SEE HOW LOW TESTING DELAYS HAVE TO GO TO HELP
    
    # testing_delays <- outbreak_sims_by_params(
    #   n_sims = 5,
    #   n_iterations = 60,
    #   n_pop = 100000, 
    #   hh_size_data = list(
    #     `1` = rbinom(100000, size = 5, prob = 0.5) + 1
    #   ),
    #   n_initial_cases = c(100), group_props = c(1), dt_approx = 1,
    #   test_results_delay_data = c(2, 2), recov_val = 7,
    #   test_sensitivity = 0.9, 
    #   probs_self_test_df =        probs_df_basic_2 %>% mutate(prob = c(c(0, 0.3, 0.525, 0.5, 0.75), c(0, 0.3, 0.525, 0.5, 0.75))),
    #   probs_isolate_symptoms_df = probs_df_basic_2 %>% mutate(prob = c(c(0, 0.3, 0.525, 0.5, 0.75), c(0, 0.3, 0.525, 0.5, 0.75))),
    #   probs_isolate_test_df =     probs_df_basic_2 %>% mutate(prob = c(c(0.45, 0.525, 0.6, 0.75, 0.75), c(0.45, 0.525, 0.6, 0.75, 0.75))),
    #   probs_isolate_ct_df =       probs_df_basic_2 %>% mutate(prob = c(c(0.225, 0.3, 0.45, 0.75, 0.75), c(0.225, 0.3, 0.45, 0.75, 0.75))),
    #   p_hh_quarantine = c(0.5, 0.5),
    #   p_contact_if_isolated_home = c(0.3, 0.3), 
    #   p_contact_traced =      c(0.2, 0.2),
    #   beta_matrix = crossing(to = 1:2, from = 1:2) %>% mutate(beta_val = c(6, 2, 2, 6)), 
    #   contacts_mean = c(30, 30),
    #   sar_out = c(0.05, 0.05),
    #   secondary_timing_params = inf_profile_params,       sar_home = c(0.2, 0.2),
    #   alpha = c(0, 0),
    #   infectiousness_by_symptoms = c(0.7, 1, 1, 1, 1)
    # )
    
    
    
    

# ..................... ---------------------------------------------------



# Testing rescaling functions ---------------------------------------------

    beta_default_3 <- crossing(to = 1:3, from = 1:3) %>% mutate(beta_val = c(10, 5, 2, 5, 7, 1, 2, 1, 2))
    
    b_rescaled <- rescale_beta(
      beta_matrix = beta_default_3,
      rescale_factors = c(1.5, 1, 1)
    )
    
    calc_contacts_mean(
      beta_matrix = beta_default_3,
      group_props = c(0.3, 0.35, 0.35),
      sar_out = c(0.2, 0.15, 0.15),
      contacts_mean_scale = 3,
      scale_on = 3
    )
    
    
    calc_contacts_mean(
      beta_matrix = b_rescaled,
      group_props = c(0.3, 0.35, 0.35),
      sar_out = c(0.2, 0.15, 0.15),
      contacts_mean_scale = 3,
      scale_on = 3
    )     # but actually even the thing being scaled on should change!! (hence next function)
    
    
    
    rescale_contacts_mean(
      beta_matrix = beta_default_3,
      rescale_factors = c(1.5, 1, 1),
      sar_out = c(0.2, 0.1, 0.15),
      group_props = c(0.3, 0.4, 0.3),
      contacts_mean_scale = 3,
      scale_on = 3
    )
    
    
    
    contacts_tightened <- calc_contacts_mean(
      beta_matrix = beta_matrix_3,
      sar_out = c(0.1, 0.1, 0.1),
      group_props = data_save$group_props[1:3],
      contacts_mean_scale = 10,
      scale_on = 3
    ) %>% print %>% 
      rescale_tighten(group_props = data_save$group_props[1:3], tighten_factor = 0.3) %>% print
    
    
    rescale_factor_tighten <- nloptr(x0 = c(1, 1, 1), 
                                     eval_f = find_the_rescale_factor_vec, 
                                     lb = c(0, 0, 0),
                                     contacts_mean_aim = contacts_tightened,
                                     opts = list(
                                       "algorithm" = "NLOPT_LN_SBPLX",
                                       "xtol_abs"=1.0e-10,
                                       "maxeval" = 200,
                                       "print_level" = 3
                                     ),
                                     beta_matrix = beta_matrix_3,
                                     sar_out = c(0.1, 0.1, 0.1),
                                     group_props = data_save$group_props[1:3],
                                     contacts_mean_scale = 10,
                                     scale_on = 3)   
    
    rescale_contacts_mean(
      rescale_factors = rescale_factor_tighten$solution,
      beta_matrix = beta_matrix_3,
      sar_out = c(0.1, 0.1, 0.1),
      group_props = data_save$group_props[1:3],
      contacts_mean_scale = 10,
      scale_on = 3
    )
    
    
    
    
    
    
    


    

    