    
    load("data/processed/tweak_final_sim_v2.RData")
    # load("data/processed/tweak_final_sim.RData")
    load("data/processed/tweaks_mean.RData", verbose = TRUE)
    # Gives you final_sim
    
    load("data/processed/group_diff_slim.RData")
    
    
    
# PLOTS to check ----------------------------------------------------------
    
    sim_to_plot <- final_sim
    
    roll_smooth <- 14
    
    # lower_bound_t_plot <- updated_params$lower_bound_t
    lower_bound_date <- ymd("2020-07-01")
    

# Mobility adjustments ----------------------------------------------------

    
    
    # 1. Plot mobility adjustments
    policy_vec_tweaks <- tweaks_df_slim %>% 
      rowwise() %>% 
      mutate(diff_times = list(lead(policy_start_times) - policy_start_times),
             diff_times = list(if_else(is.na(diff_times), 500, diff_times)),
             policy_vec = list(rep(mobility_factors, times = diff_times)),
             t = list(0:(length(policy_vec) -  1))) %>% 
      # filter(iteration %in% c(1, 5, 10, 15, 19)) %>%
      tail() %>% 
      # filter(iteration %in% c(0, 10, 15, 16, 17)) %>%
      unnest(c(policy_vec, t))
  
    policy_vec_tweaks %>% 
      ggplot(aes(x = t, y = policy_vec, colour = factor(iteration))) + 
      geom_line() + 
      coord_cartesian(xlim = c(0, 250)) + 
      scale_y_continuous(breaks = seq(0, 3, 0.5), minor_breaks = seq(0, 3, 0.1))
    
    
  
 
# Merge real/model totals  -----------------------------------------------------------

    # REAL DATA
    real_cases_total <- confirmed_cases_total
    
    # Get the rolling average pc in real data
    real_cases_slim <- real_cases_total %>%
      filter(n_confirmed_new != 0) %>% 
      mutate(across(c(n_confirmed_new_pc), list(roll = ~ zoo::rollmean(.x, k = roll_smooth, na.pad = TRUE, align = "center")))) %>% 
      mutate(real_pc_roll = n_confirmed_new_pc_roll)
    
    # What's the value of pc_roll on lower bound date
    lower_bound_date_value <- real_cases_slim %>% filter(date_results == lower_bound_date) %>% .$real_pc_roll
    
    
    # MODEL DATA
    model_cases_total <- time_series_to_case_count(sim_to_plot)
    
    model_cases_slim <- model_cases_total %>% 
      mutate(across(c(n_confirmed_new_model_pc), list(roll = ~ zoo::rollmean(.x, k = roll_smooth, na.pad = TRUE, align = "center")))) %>% 
      mutate(model_pc_roll = n_confirmed_new_model_pc_roll) %>% 
      select(sim_id, t, model_pc_roll) %>% 
      group_by(t) %>% 
      quantile_summarise(model_pc_roll, conf_level = 0.95)
    
    # lower bound date is equivalent to which t?
    lower_bound_t <- model_cases_slim %>% filter(model_pc_roll_mean >= lower_bound_date_value) %>% .$t %>% first()
    
    model_data_with_date <- model_cases_slim %>% 
      mutate(date_model = days(t) + lower_bound_date - days(lower_bound_t))
    
    
    # MERGE REAL AND MODEL
    match_model <- real_cases_slim %>% 
      left_join(model_data_with_date, by = c("date_results" = "date_model"))
    
    max_date <- max(match_model$date_results, na.rm = TRUE) - days(detection_delay)
    max_t <- match_model$t[match_model$date_results == max_date & !is.na(match_model$t)]
    

    
    
    

# Merge real / model by group ---------------------------------------------
    
    # All data
    model_by_group <- sim_to_plot %>% 
      select(sim_id, t, i_group, policy, n_pop, n_confirmed_cases) %>% 
      arrange(sim_id, i_group, t) %>% 
      group_by(sim_id, i_group) %>% 
      mutate(
        model_new_confirmed = n_confirmed_cases - lag(n_confirmed_cases),
        model_new_confirmed_pc = model_new_confirmed / n_pop,
        # model_new_confirmed_week_roll = zoo::rollsum(model_new_confirmed, k = 20/7, fill = NA, align = "center") / (20/7),
        # model_new_confirmed_week_roll_pc = model_new_confirmed_week_roll / n_pop,
        model_cum_confirmed = n_confirmed_cases,
        model_cum_confirmed_pc = n_confirmed_cases / n_pop
        # new_cases_per_cap = (n_new_cases_week_roll / n_pop)
      ) %>% 
      select(-n_confirmed_cases) %>% 
      mutate(date_model = days(t) + lower_bound_date - days(lower_bound_t))
    
    
    # BY GROUP
    model_by_group_summ <- model_by_group %>% 
      group_by(t, date_model, i_group) %>% 
      quantile_summarise(c(model_new_confirmed_pc, model_cum_confirmed_pc), conf_level = 0.9) %>% 
      group_by(i_group) %>% 
      mutate(
        across(starts_with("model_new_confirmed_pc_"), 
               list(roll = ~ zoo::rollmean(.x, k = 7, fill = NA_real_, align = "center")))
      ) %>% 
      ungroup
    
    # model_by_group_summ %>% view_n()
    
    model_by_group_with_data <- full_join(
      confirmed_cases_by_group %>% rename(
        data_new_confirmed = n_confirmed_new,
        data_new_confirmed_pc = n_confirmed_new_pc,
        data_cum_confirmed_pc = n_confirmed_cum_pc,
        data_cum_confirmed = n_confirmed_cum),
      model_by_group_summ,
      by = c("date_results" = "date_model", "i_group")
    ) %>% 
      mutate(
        across(c(data_new_confirmed_pc), 
               list(roll = ~ zoo::rollmean(.x, k = 7, fill = NA_real_, align = "center")))
      ) %>% 
      rename(date = date_results) %>% 
      mutate(i_group = recode_i_group(i_group))
    
    
    

# All groups (for last graph) ---------------------------------------------

    model_data_with_date <- model_cases_slim %>% 
      mutate(date_model = days(t) + lower_bound_date - days(lower_bound_t))
    
    
    # MERGE REAL AND MODEL
    cum_all_groups <- full_join(
      real_cases_slim %>% select(date_results, n_confirmed_cum_pc),
      model_cases_total %>% mutate(date_model = days(t) + lower_bound_date - days(lower_bound_t)) %>% select(sim_id, date_model, n_confirmed_cum_model_pc),
      by = c("date_results" = "date_model")
    ) 
    

# PLOT COMPARISONS --------------------------------------------------------

    # NEW CONFIRMED CASES (TOTAL)
    ggplot(match_model, aes(x = date_results)) +
      geom_line(aes(y = real_pc_roll), size = 1.5, colour = "skyblue") +
      geom_ribbon(aes(ymin = model_pc_roll_lower, ymax = model_pc_roll_upper), fill = "indianred", alpha = 0.2) +
      geom_line(aes(y = model_pc_roll_mean), size = 1.5, colour = "indianred") +
      # geom_point(aes(y = model_pc_roll_median)) +
      # geom_vline(xintercept = model_deviations$date_results[model_deviations$t == t_deviate]) +
      theme_custom() +
      scale_x_date(breaks = "months", date_labels = "%b %d")
    
    ggsave("figures/tweak_match_new_total.png", width = 6, height = 4, scale = 1.5)
    
    
    # NEW BY GROUP v2 (separate with confidence intervals)
    model_by_group_with_data %>% 
      pivot_longer(c(data_new_confirmed_pc_roll, model_new_confirmed_pc_mean_roll)) %>% 
      filter(t <= max_t) %>%
      mutate(name = fct_relevel(name, "model_new_confirmed_pc_mean_roll")) %>% 
      dups_report(date, name, i_group) %>% 
      ggplot(aes(x = date)) + 
      geom_ribbon(aes(ymin = model_new_confirmed_pc_lower_roll,
                      ymax = model_new_confirmed_pc_upper_roll,
                      fill = i_group,
                      name = "model_new_confirmed_pc_mean_roll"),
                  alpha = 0.1, size = 0.3) +
      geom_line(aes(y = value, colour = i_group, linetype = name), size = 1) + 
      facet_wrap(i_group ~ .) + 
      scale_colour_viridis(option = "viridis",
                           begin = 0.3, 
                           discrete = TRUE) + 
      scale_fill_viridis(option = "viridis",
                         begin = 0.3, 
                         discrete = TRUE) + 
      theme_custom()
    
    ggsave("figures/tweak_match_new_by_group_conf_int.png", width = 6, height = 4, scale = 1.5)
    
    
    # NEW BY GROUP 
    model_by_group_with_data %>% 
      pivot_longer(c(data_new_confirmed_pc_roll, model_new_confirmed_pc_mean_roll)) %>% 
      filter(t <= max_t) %>%
      # mutate(name = fct_relevel(name, "model_new_confirmed_pc_mean_roll")) %>% 
      dups_report(date, name, i_group) %>% 
      ggplot(aes(x = date)) + 
      # geom_ribbon(aes(ymin = model_new_confirmed_pc_lower_roll,
      #                 ymax = model_new_confirmed_pc_upper_roll,
      #                 fill = i_group,
      #                 name = "model_new_confirmed_pc_mean_roll"),
      #             alpha = 0.1, size = 0.3) +
      geom_line(aes(y = value, colour = i_group), size = 1) + 
      facet_wrap( ~ name) + 
      scale_colour_viridis(option = "viridis",
                           begin = 0.3, 
                           discrete = TRUE) + 
      scale_fill_viridis(option = "viridis",
                         begin = 0.3, 
                         discrete = TRUE) + 
      theme_custom()
    
    ggsave("figures/tweak_match_new_by_group.png", width = 6, height = 4, scale = 1.5)
    
    
    # CUMULATIVE BY GROUP
    model_by_group_with_data %>% 
      # group_by(t, date, i_group, data_cum_confirmed_pc) %>% 
      # quantile_summarise(c(model_cum_confirmed_pc), conf_level = 0.9) %>% 
      filter(t <= max_t) %>%
      ggplot(aes(x = date, colour = i_group)) + 
      geom_line(aes(y = model_cum_confirmed_pc_mean), size = 1.5) + 
      geom_ribbon(aes(ymin = model_cum_confirmed_pc_lower,
                      ymax = model_cum_confirmed_pc_upper,
                      fill = i_group),
                  alpha = 0.1, size = 0.3) +
      geom_line(aes(y = data_cum_confirmed_pc), linetype = "dashed", size = 1) + 
      scale_colour_viridis(option = "viridis",
                           begin = 0.3, 
                           discrete = TRUE) + 
      scale_fill_viridis(option = "viridis",
                         begin = 0.3, 
                         discrete = TRUE) + 
      theme_custom() + 
      facet_wrap(~ i_group) + 
      scale_x_date(breaks = "months", date_labels = "%b %d")
    
    ggsave("figures/tweak_match_cum_by_group.png", width = 6, height = 4, scale = 1.5)
    
    match_model %>% print_names
    
    # n_confirmed_cum_pc
    
    # CUMULATIVE ALL GROUPS
    ggplot(cum_all_groups, aes(x = date_results)) + 
      geom_line(aes(y = n_confirmed_cum_model_pc, group = factor(sim_id)), size = 0.1) + 
      # geom_line(data = cum_all_groups %>% filter(median_model), aes(y = model_cum_confirmed_pc, group = factor(sim_id)), colour = "skyblue", size = 1) +
      geom_line(aes(y = n_confirmed_cum_pc), colour = "indianred", size = 1.5) + 
      theme_custom() + 
      scale_x_date(breaks = "months", date_labels = "%b %d")
    
    ggsave("figures/tweak_match_cum_total.png", width = 6, height = 4, scale = 1.5)
    
    

    

    
    
    

# LEGACY ------------------------------------------------------------------


    # 
    # # >>>
    # 
    # 
    # # FIND THE "MEDIAN MODEL"
    # median_model <- model_by_group_with_data %>%
    #   group_by(sim_id, t, date) %>% 
    #   # filter(date <= ymd("2020-07-01")) %>% 
    #   summarise(across(c(model_cum_confirmed_pc, data_cum_confirmed_pc), sum)) %>% 
    #   group_by(sim_id) %>% 
    #   mutate(final_cum = max(model_cum_confirmed_pc)) %>% 
    #   ungroup %>% 
    #   mutate(median_model = final_cum == median(final_cum, na.rm = TRUE)) %>% 
    #   filter(median_model) %>% 
    #   .$sim_id %>% 
    #   unique()
    # 
    # # EXAMINE NEW CASES for the median model
    # median_model_data <- manual_adjustment_total_cases %>% 
    #   filter(sim_id == median_model) %>% 
    #   mutate(model_pc = new_cases_detected / n_pop_total,
    #          model_cum_pc = n_confirmed_cases / n_pop_total) %>% 
    #   mutate(across(c(model_pc), list(roll = ~ zoo::rollmean(.x, k = roll_smooth, na.pad = TRUE, align = "center")))) %>%
    #   select(sim_id, t, model_pc, model_cum_pc, model_pc_roll) %>%
    #   mutate(date_model = days(t) + lower_bound_date - days(lower_bound_t))
    # 
    # 
    # 
    # # manual_adjustment_total_cases %>% 
    # # ungroup %>% 
    # # view_n()
    # 
    # # MERGE REAL AND MODEL
    # real_data %>% 
    #   left_join(median_model_data, by = c("date_results" = "date_model")) %>% 
    #   
    #   ggplot(aes(x = date_results)) +
    #   # geom_line(aes(y = median_model_pc_roll - real_pc_roll )) +
    #   geom_point(aes(y = real_pc_roll), size = 1.5, colour = "skyblue") +
    #   # geom_line(aes(y = model_pc), size = 1.5, colour = "grey") +
    #   geom_point(aes(y = model_pc_roll)) +
    #   theme_custom() +
    #   scale_x_date(breaks = "months", date_labels = "%b %d")
    # 
    # 
    # 
    # # ALL GROUPS
    # confirmed_cases_all %>% 
    #   mutate(real_pc = n_confirmed_cum / monthly_positivity$n_pop_total[[1]]) %>% 
    #   left_join(median_model_data, by = c("date_results" = "date_model")) %>% 
    #   
    #   ggplot(aes(x = date_results)) + 
    #   geom_line(aes(y = real_pc), size = 1.5, colour = "skyblue") +
    #   geom_point(aes(y = model_cum_pc)) +
    #   theme_custom() +
    #   scale_x_date(breaks = "months", date_labels = "%b %d")
    # 
    # 
    # # DIFF CUMULATIVE
    # confirmed_cases_all %>% 
    #   mutate(real_pc = n_confirmed_cum / monthly_positivity$n_pop_total[[1]]) %>% 
    #   left_join(median_model_data, by = c("date_results" = "date_model")) %>% 
    #   
    #   ggplot(aes(x = date_results)) + 
    #   geom_line(aes(y = model_cum_pc - real_pc), size = 1.5, colour = "skyblue") +
    #   # geom_point(aes(y = model_cum_pc)) +
    #   theme_custom() +
    #   scale_x_date(breaks = "months", date_labels = "%b %d")
    # 
    # # DIFF NEW
    # real_data %>% 
    #   left_join(median_model_data, by = c("date_results" = "date_model")) %>% 
    #   
    #   ggplot(aes(x = date_results)) +
    #   # geom_line(aes(y = median_model_pc_roll - real_pc_roll )) +
    #   geom_point(aes(y = model_pc_roll - real_pc_roll), size = 1.5, colour = "skyblue") +
    #   # geom_line(aes(y = model_pc), size = 1.5, colour = "grey") +
    #   # geom_point(aes(y = model_pc_roll)) +
    #   theme_custom() +
    #   scale_x_date(breaks = "months", date_labels = "%b %d")
    # 
    # 
    # 
    # 
    # manual_adjustment_total_cases %>% 
    #   mutate(date_model = days(t) + lower_bound_date - days(lower_bound_t)) %>% 
    #   filter(sim_id == median_model)
    # 
    # 
    # 
    # 
    # median_model_data %>% 
    #   ungroup %>% 
    #   mutate(cum_check = cumsum(if_else(is.na(model_pc), 0, model_pc))) %>% 
    #   mutate(check_tf = abs(model_cum_pc -  cum_check) < 0.000001) %>% 
    #   print_all
    # 
    # 
    # ?cumsum
    # 
    # 
    # 
    # 
    # 
    # 
    # 
    # 
