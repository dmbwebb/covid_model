
# Set up ------------------------------------------------------------------

      


    
    library("viridis")



# ...... ------------------------------------------------------------------


# BASELINE ----------------------------------------------------------------


    # load("data/processed/model_baseline_final.RData")
    
    # load("data/processed/model_baseline_final_v2.RData")
    
# 1. Basic graph of main simulation ---------------------------------------------------------
    
    live_cases <- model_baseline_time_series %>% 
      mutate(n_cases_per_cap = n_cases_live / n_pop) %>% 
      group_by(t, i_group) %>% 
      summarise_quantiles(n_cases_per_cap) %>% 
      filter(q %in% c(0.05, 0.5, 0.95)) %>% 
      mutate(q_label = case_when(q == 0.05 ~ "_lower",
                                 q == 0.5   ~ "",
                                 q == 0.95 ~ "_upper")) %>% 
      select(-q) %>%
      ungroup %>% 
      pivot_wider(names_from = q_label,
                  values_from = n_cases_per_cap,
                  names_prefix = "n_cases_per_cap")
    
    live_cases %>% 
      mutate(i_group = recode_i_group(i_group)) %>% 
      ggplot(aes(x = t, y = n_cases_per_cap)) + 
      geom_line(aes(colour = factor(i_group)), size = 1) + 
      geom_ribbon(aes(ymax = n_cases_per_cap_upper, ymin = n_cases_per_cap_lower, fill = factor(i_group)), alpha = 0.15) + 
      # geom_label(aes(label = factor(stratum))) + 
      scale_colour_viridis(option = "viridis",
                           begin = 0.3, 
                           discrete = TRUE) + 
      scale_y_continuous(labels = scales::percent_format(accuracy = 1)) + 
      scale_fill_viridis(option = "viridis",
                         begin = 0.3, 
                         discrete = TRUE) + 
      theme_custom() + 
      theme(legend.position = c(0.12, 0.67)) + 
      labs(x = "t", y = "Prevalence", fill = "SES Group",
           colour = "SES Group") + # + 
      geom_hline(yintercept = 0, colour = "darkgrey")
    
    ggsave("figures/live_cases.png", scale = 0.8)
    
    
    

# 1b - Proportion infected ------------------------------------------------

    
    prop_infected <- model_baseline_time_series %>% 
      mutate(n_cases_per_cap = n_cases_cum / n_pop) %>% 
      group_by(t, i_group) %>% 
      summarise_quantiles(n_cases_per_cap) %>% 
      filter(q %in% c(0.05, 0.5, 0.95)) %>% 
      mutate(q_label = case_when(q == 0.05 ~ "_lower",
                                 q == 0.5   ~ "",
                                 q == 0.95 ~ "_upper")) %>% 
      select(-q) %>%
      ungroup %>% 
      pivot_wider(names_from = q_label,
                  values_from = n_cases_per_cap,
                  names_prefix = "n_cases_per_cap")
    
    prop_infected %>% 
      mutate(i_group = recode_i_group(i_group)) %>% 
      ggplot(aes(x = t, y = n_cases_per_cap)) + 
      geom_line(aes(colour = factor(i_group)), size = 1) + 
      geom_ribbon(aes(ymax = n_cases_per_cap_upper, ymin = n_cases_per_cap_lower, fill = factor(i_group)), alpha = 0.15) + 
      # geom_label(aes(label = factor(stratum))) + 
      scale_colour_viridis(option = "viridis",
                           begin = 0.3, 
                           discrete = TRUE) + 
      scale_y_continuous(labels = scales::percent_format(accuracy = 1)) + 
      scale_fill_viridis(option = "viridis",
                         begin = 0.3, 
                         discrete = TRUE) + 
      theme_custom() + 
      theme(legend.position = c(0.12, 0.67)) + 
      labs(x = "t", y = "Prop infected", fill = "SES Group",
           colour = "SES Group") + # + 
      geom_hline(yintercept = 0, colour = "darkgrey")
    
    # ggsave("figures/live_cases.png", scale = 0.8)
    
    
    
    
# 2. New confirmed cases (model) ------------------------------------------------
    
    
    
    
    confirmed_cases_model <- model_baseline_time_series %>% 
      select(sim_id, t, i_group, n_pop, n_confirmed_cases) %>% 
      arrange(sim_id, i_group, t) %>% 
      group_by(sim_id, i_group) %>% 
      mutate(
        new_confirmed_cases = n_confirmed_cases - lag(n_confirmed_cases),
        n_new_cases_week_roll = zoo::rollsum(new_confirmed_cases, k = 20, fill = NA, align = "right") / (20/7),
        n_confirmed_cum = n_confirmed_cases / n_pop,
        new_cases_per_cap = (n_new_cases_week_roll / n_pop)
      ) %>% 
      group_by(i_group, t) %>% 
      summarise_quantiles(c(new_cases_per_cap, n_confirmed_cum)) %>% 
      filter(q %in% c(0.05, 0.5, 0.95)) %>% 
      mutate(q_label = case_when(q == 0.05 ~ "_lower",
                                 q == 0.5   ~ "",
                                 q == 0.95 ~ "_upper")) %>% 
      select(-q) %>%
      ungroup %>% 
      pivot_wider(names_from = q_label,
                  values_from = c(new_cases_per_cap, n_confirmed_cum),
                  names_sep = "")
    
    
    ggplot(confirmed_cases_model, aes(x = t, y = new_cases_per_cap)) + 
      geom_line(aes(colour = factor(i_group)), size = 1) + 
      geom_ribbon(aes(ymax = new_cases_per_cap_upper, ymin = new_cases_per_cap_lower, fill = factor(i_group)), alpha = 0.1) + 
      # geom_label(aes(label = factor(stratum))) + 
      scale_colour_viridis(option = "viridis",
                           begin = 0.3, 
                           discrete = TRUE) + 
      scale_fill_viridis(option = "viridis",
                         begin = 0.3, 
                         discrete = TRUE) + 
      theme_custom() + 
      theme(legend.position = c(0.12, 0.67)) + 
      labs(x = "t", y = "Weekly new confirmed cases per capita in group [model]", fill = "SES Group",
           colour = "SES Group") + # + 
      geom_hline(yintercept = 0, colour = "darkgrey")
    
    
    
    
    

    
    

# ..... -------------------------------------------------------------------


# MOBILITY CHANGE ---------------------------------------------------------
    
    load("data/processed/policy_change_time_series.RData")

    # FOR PLOTTING the policy change graphs
    lockdown_rect_1 <- policy_change_time_series %>% 
      summarise(xmin = min(t[policy == "1_severe_lockdown"], na.rm = TRUE), 
                xmax = max(t[policy == "1_severe_lockdown"], na.rm = TRUE))
    
    # lockdown_rect_2 <- policy_change_time_series %>% 
    #   summarise(xmin = min(t[policy == "2_lighter_lockdown"], na.rm = TRUE), 
    #             xmax = max(t[policy == "2_lighter_lockdown"], na.rm = TRUE))
    
    
    policy_change_time_series %>% filter(i_group ==1) %>% 
      print(n = 200)
    
    

# 1. Basic graph ----------------------------------------------------------
    
    
    
    # PLOT POLICY CHANGE GRAPH
    policy_change_time_series %>% 
      group_by(t, i_group, n_pop, policy) %>% 
      mutate(prevalence_prop = n_cases_live / n_pop) %>% 
      quantile_summarise(c(prevalence_prop), conf_level = 0.9) %>% 
      mutate(i_group = recode_i_group(i_group)) %>% 
      # ungroup %>% 
      ggplot() + 
      
      # Shaded area
      geom_rect(data = lockdown_rect_1,
                aes(xmin = xmin,
                    xmax = xmax,
                    ymin=0, ymax=Inf), alpha = 0.08, fill = "indianred") +
      
      # geom_rect(data = lockdown_rect_2,
      #           aes(xmin = xmin,
      #               xmax = xmax,
      #               ymin=0, ymax=Inf), alpha = 0.08, fill = "skyblue") +
      
      # LINE
      geom_line(aes(x = t, y = prevalence_prop_median, 
                    colour = factor(i_group), group = factor(i_group)),
                show.legend = TRUE, size = 1) + 
      geom_ribbon(
        aes(
          x = t, 
          ymax = prevalence_prop_upper,
          ymin = prevalence_prop_lower,
          fill = factor(i_group),
          # colour = factor(i_group)
        ),
        linetype = 2,
        alpha = 0.1
      ) +
      # facet_wrap(~ parameter_set, labeller = as_labeller(counterfactual_labels)) + 
      
      # Formatting
      theme_custom() + 
      theme(legend.position = "right") + 
      scale_alpha_continuous(guide = FALSE) + 
      scale_y_continuous(labels = scales::percent_format(accuracy = 1), minor_breaks = NULL) + 
      labs(colour = "SES Group", fill = "SES Group", y = "% Infected at t") +
      geom_hline(yintercept = 0, colour = "darkgrey", size = 0.4) + 
      # coord_cartesian(xlim = c(0, 300)) +
      scale_colour_viridis(option = "viridis",
                           begin = 0.3, 
                           discrete = TRUE) + 
      scale_fill_viridis(option = "viridis",
                         begin = 0.3, 
                         discrete = TRUE) + 
      # scale_x_continuous(minor_breaks = c(0, 25, 50, 75), breaks = c(0, 50, 100, 150, 200)) + 
      geom_vline(aes(xintercept = min(t[policy == "lockdown"], na.rm = TRUE)), linetype = "dashed") +
      geom_vline(aes(xintercept = max(t[policy == "lockdown"], na.rm = TRUE)), linetype = "dashed")
    
    
    ggsave("figures/second_wave_example.png", width = 6, height = 4, scale = 0.85)    
    
    

    

# 2. Compare prop infected across both -------------------------------------------------------------

    
    prop_infected_policy_change <- policy_change_time_series %>% 
      mutate(n_cases_per_cap = n_cases_cum / n_pop) %>% 
      group_by(t, i_group) %>% 
      summarise_quantiles(n_cases_per_cap) %>% 
      filter(q %in% c(0.05, 0.5, 0.95)) %>% 
      mutate(q_label = case_when(q == 0.05 ~ "_lower",
                                 q == 0.5   ~ "",
                                 q == 0.95 ~ "_upper")) %>% 
      select(-q) %>%
      ungroup %>% 
      pivot_wider(names_from = q_label,
                  values_from = n_cases_per_cap,
                  names_prefix = "n_cases_per_cap")
    
    
    prop_infected_both <- bind_rows(
      "baseline" = prop_infected,
      "mobility_change" = prop_infected_policy_change,
      .id = "model_type"
    )
    
    
    prop_infected_both %>% 
      mutate(i_group = recode_i_group(i_group)) %>% 
      ggplot(aes(x = t, y = n_cases_per_cap)) + 
      geom_line(aes(colour = factor(i_group)), size = 1) + 
      geom_ribbon(aes(ymax = n_cases_per_cap_upper, ymin = n_cases_per_cap_lower, fill = factor(i_group)), alpha = 0.15) + 
      # geom_label(aes(label = factor(stratum))) + 
      scale_colour_viridis(option = "viridis",
                           begin = 0.3, 
                           discrete = TRUE) + 
      scale_y_continuous(labels = scales::percent_format(accuracy = 1)) + 
      scale_fill_viridis(option = "viridis",
                         begin = 0.3, 
                         discrete = TRUE) + 
      theme_custom() + 
      theme(legend.position = c(0.12, 0.67)) + 
      labs(x = "t", y = "Prop infected", fill = "SES Group",
           colour = "SES Group") + # + 
      geom_hline(yintercept = 0, colour = "darkgrey") + 
      facet_wrap(~ model_type)
  
    
    # CUMULATIVE - cases!!
    policy_change_time_series %>% 
      group_by(sim_id, t) %>% 
      summarise(
        n_confirmed_cases = sum(n_confirmed_cases),
        n_cases_cum = sum(n_cases_cum),
        n_pop = sum(n_pop)) %>% 
      mutate(n_cases_per_cap = n_cases_cum / n_pop,
             n_confirmed_cases_per_cap = n_confirmed_cases / n_pop) %>% 
      ungroup %>% 
      left_join(
        confirmed_cases_actual %>% group_by(t) %>% summarise(n_confirmed_cum = sum(n_confirmed_cum_abs)),
        by = "t"
      ) %>% 
      mutate(
        actual_per_cap = n_confirmed_cum / 7600000
      ) %>% 
      ggplot(aes(x = t, y = n_cases_per_cap)) + 
      geom_line(size = 1) + 
      geom_line(aes(y = n_confirmed_cases_per_cap), colour = "indianred") + 
      geom_line(aes(y = actual_per_cap), colour = "indianred", linetype = "dashed") +
      theme_custom() + 
      theme(legend.position = c(0.12, 0.67)) + 
      labs(x = "t", y = "Prop infected", fill = "SES Group",
           colour = "SES Group") + # + 
      geom_hline(yintercept = 0, colour = "darkgrey")
    

# 3. Confirmed cases in mobility change model -----------------------------

    confirmed_cases_model_mobility_change <- policy_change_time_series %>% 
      select(sim_id, t, i_group, n_pop, n_confirmed_cases) %>% 
      arrange(sim_id, i_group, t) %>% 
      group_by(sim_id, i_group) %>% 
      mutate(
        new_confirmed_cases = n_confirmed_cases - lag(n_confirmed_cases),
        n_new_cases_week_roll = zoo::rollsum(new_confirmed_cases, k = 20, fill = NA, align = "right") / (20/7),
        n_confirmed_cum = n_confirmed_cases / n_pop,
        new_cases_per_cap = (n_new_cases_week_roll / n_pop)
      ) %>% 
      group_by(i_group, t) %>% 
      summarise_quantiles(c(new_cases_per_cap, n_confirmed_cum)) %>% 
      filter(q %in% c(0.05, 0.5, 0.95)) %>% 
      mutate(q_label = case_when(q == 0.05 ~ "_lower",
                                 q == 0.5   ~ "",
                                 q == 0.95 ~ "_upper")) %>% 
      select(-q) %>%
      ungroup %>% 
      pivot_wider(names_from = q_label,
                  values_from = c(new_cases_per_cap, n_confirmed_cum),
                  names_sep = "")
    
    
    
    
# 3. Real and model new confirmed cases in same graph ----------------------------------
    
    
    # WHEN IS Start date on model - April 1st when there are 486 cases
    start_cases <- sum(data_save$best_guess_initial_cases$n) 
    start_cases_prop <- start_cases / sum(data_save$best_guess_initial_cases$stratum_pop)
    format(start_cases_prop, scientific=F)   # 0.00006578947  of the population
    
    
    # Calculate which t should align with April 1st
    start_t <- model_baseline_time_series %>% 
      group_by(sim_id, t) %>% # sum all i_groups
      summarise(across(c(n_pop, n_confirmed_cases), sum)) %>% 
      ungroup %>% 
      mutate(n_confirmed_prop = n_confirmed_cases / n_pop,
             over_start_thresh = n_confirmed_prop >= start_cases_prop) %>% 
      
      # Find first t that goes over starting threshold
      group_by(sim_id) %>% 
      summarise(start_t = first_non_na(t[over_start_thresh])) %>% 
      
      # Take rounded median of all as starting t
      summarise(start_t = floor(median_na(start_t))) %>% 
      .$start_t
    
    print(str_glue("Median starting time (April 1st) is t = {start_t}"))
    
    

    stratum_pops <- data_save$sds_case_data %>% 
      mutate(stratum = as.integer(stratum)) %>% 
      mutate(stratum = case_when(stratum %in% c(1, 2) ~ 1L,
                                 stratum %in% c(3, 4) ~ stratum - 1L,
                                 stratum %in% c(5, 6) ~ 4L)) %>% 
      group_by(stratum, stratum_pop) %>% 
      summarise() %>% 
      filter(!is.na(stratum)) %>% 
      summarise(stratum_pop = sum(stratum_pop))
    
    confirmed_cases_actual <- data_save$sds_case_data %>% 
      mutate(stratum = as.integer(stratum)) %>% 
      mutate(stratum = case_when(stratum %in% c(1, 2) ~ 1L,
                                 stratum %in% c(3, 4) ~ stratum - 1L,
                                 stratum %in% c(5, 6) ~ 4L)) %>% 
      filter(!is.na(stratum)) %>% 
      group_by(stratum, date_results) %>% 
      summarise(n_cases_day = n()) %>% 
      full_join(
        crossing(stratum = 1:4, date_results = ymd("2020-03-01") + days(1:400)), by = c("stratum", "date_results")
      ) %>% 
      mutate(n_cases_day = if_else(is.na(n_cases_day), 0L, n_cases_day)) %>% 
      left_join(
        stratum_pops, by = "stratum"
      ) %>% 
      # Calculate model t in real data
      mutate(t = interval_days(ymd("2020-04-01"), date_results) + start_t) %>% 
      arrange(stratum, date_results) %>% 
      group_by(stratum, stratum_pop) %>% 
      mutate(n_new_cases_week_roll = zoo::rollsum(n_cases_day, k = 7, fill = NA, align = "right")) %>% 
      mutate(new_cases_per_cap = (n_new_cases_week_roll / stratum_pop)) %>% 
      mutate(n_confirmed_cum = cumsum(n_cases_day) / stratum_pop,
             n_confirmed_cum_abs = cumsum(n_cases_day)) %>% 
      ungroup %>% 
      filter(
        # date_results <= ymd("2020-11-01")
        # date_results > ymd("2020-04-01")
      ) %>% 
      mutate(
        new_cases_per_cap = ifelse(date_results >= ymd("2020-12-01"), NA, new_cases_per_cap),
        n_confirmed_cum = ifelse(date_results >= ymd("2020-12-01"), NA, n_confirmed_cum)
      ) %>% 
      select(i_group = stratum, t, date_results, new_cases_per_cap, n_confirmed_cum, n_confirmed_cum_abs)
      # mutate(source = "data")
    
    t_adj <- -15
    
    date_to_t <- confirmed_cases_actual %>% 
      mutate(t = t - t_adj) %>% 
      select(t, date_results) %>% dups_drop()
    
    
    # combine both into one DF
    confirmed_cases_combined <- bind_rows(
      "data" = confirmed_cases_actual %>% select(-date_results) %>% mutate(t = t - t_adj),
      "model_baseline" = confirmed_cases_model,
      # "model_mobility_change" = confirmed_cases_model_mobility_change,
      .id = "source"
    ) %>% 
      left_join(date_to_t, by = "t")
    
    
    # PLOT - weekly new confirmed
    confirmed_cases_combined %>% 
      # filter(source != "model_mobility_change") %>% 
      filter(t > 0) %>% 
      ggplot(aes(x = date_results, y = new_cases_per_cap)) + 
      geom_line(aes(colour = factor(i_group)), size = 1) + 
      geom_ribbon(aes(ymax = new_cases_per_cap_upper, ymin = new_cases_per_cap_lower, fill = factor(i_group)), alpha = 0.2) +
      # geom_label(aes(label = factor(stratum))) + 
      scale_colour_viridis(option = "viridis",
                           begin = 0.3, 
                           discrete = TRUE) + 
      scale_fill_viridis(option = "viridis",
                         begin = 0.3, 
                         discrete = TRUE) + 
      theme_custom() + 
      # theme(legend.position = c(0.12, 0.67)) + 
      labs(x = "Date", y = "Weekly new confirmed cases per capita in group", fill = "SES Group",
           colour = "SES Group") + # + 
      geom_hline(yintercept = 0, colour = "darkgrey") + 
      facet_wrap(~ source)
    
    
    
    
    
    # PLOT
    confirmed_cases_combined %>% 
      # mutate(n_confirmed_cum = if_else(source == "model_baseline", n_confirmed_cum / 2, n_confirmed_cum)) %>%
      mutate(i_group = recode_i_group(i_group)) %>% 
      filter(t > 0, t < 300) %>% 
      ggplot(aes(x = date_results, y = n_confirmed_cum)) + 
      geom_line(aes(colour = factor(i_group), linetype = source), size = 1) + 
      geom_ribbon(
        aes(
          ymax = n_confirmed_cum_upper,
          ymin = n_confirmed_cum_lower,
          fill = factor(i_group),
          colour = factor(i_group)
        ),
        linetype = 2,
        alpha = 0.1
      ) +
      # geom_label(aes(label = factor(stratum))) + 
      scale_colour_viridis(option = "viridis",
                           begin = 0.3, 
                           discrete = TRUE) + 
      scale_fill_viridis(option = "viridis",
                         begin = 0.3, 
                         discrete = TRUE) + 
      theme_custom() + 
      # theme(legend.position = c(0.12, 0.67)) + 
      labs(x = "Date", y = "Cumulative confirmed cases per capita\nin each group", fill = "SES Group",
           colour = "SES Group") + # + 
      geom_hline(yintercept = 0, colour = "darkgrey")
      # facet_wrap(~ source, labeller = as_labeller(c("data" = "Data", "model_baseline" = "Model (no mobility change)", "model_mobility_change" = "Model (mobility change)")))
    
    
    # ggsave("figures/data_vs_model_confirmed_cum.png", scale = 0.7)
    
    
    
    
    
    

# ..... -------------------------------------------------------------------


