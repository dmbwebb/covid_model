    load("data/processed/group_diff_tightened_half.RData", verbose = TRUE)
    # load("data/processed/tweak_final_sim.RData", verbose = TRUE) # or could use tweak_final_sim_v2
    load("data/processed/tweak_final_sim_v2.RData", verbose = TRUE)
    # load("data/processed/tweaks_mean.RData", verbose = TRUE)
    # tweaks_df_slim %>% print_all
    
    library("ggrepel")
    
    # Parameters
    roll_smooth <- 14
    
    lag_period <- 14
    # conf_level <- 0.95
    
    detection_delay <- 21
    lower_bound_date <- ymd("2020-07-01") # date to match to
    
    
    # Isolate no mobility change
    baseline <- group_diff_slim_tightened %>% filter(parameter_set == "baseline")
    
    mobility_change <- final_sim
    
    model_results <- bind_rows(
      "baseline" = baseline,
      "mobility_change" = mobility_change,
      .id = "model_type"
    )
    
    
# Merge real/model totals to find the DATE to match with  -----------------------------------------------------------
    
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
    n_pop_total <- 100000
    
    model_cases_total <- model_results %>% 
      group_by(model_type, sim_id, t, policy) %>% 
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
      relocate(model_type, sim_id, t, policy, n_new_model, n_cum_model, n_confirmed_cum_model, n_confirmed_new_model) %>% 
      mutate(
        across(starts_with("n_"), list(pc = ~ . / n_pop_total)),
      )
    
    model_cases_slim <- model_cases_total %>% 
      mutate(across(c(n_confirmed_new_model_pc), list(roll = ~ zoo::rollmean(.x, k = roll_smooth, na.pad = TRUE, align = "center")))) %>% 
      mutate(model_pc_roll = n_confirmed_new_model_pc_roll) %>% 
      select(sim_id, t, model_pc_roll) %>% 
      group_by(t) %>% 
      quantile_summarise(c(model_pc_roll), conf_level = conf_level)
    
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
    
    
    
    # All data from model
    model_by_group <- model_results %>% 
      complete(t, nesting(model_type, sim_id, i_group)) %>%           # NEED TO FILL OUT MISSING VALUES AT END OF EPIDEMIC
      select(model_type, sim_id, t, i_group, policy, n_pop, n_confirmed_cases, n_cases_cum) %>% 
      arrange(model_type, sim_id, i_group, t) %>% 
      group_by(model_type, sim_id, i_group) %>% 
      fill(policy, n_pop, n_confirmed_cases, n_cases_cum, .direction = "down") %>% # NEED TO FILL OUT MISSING VALUES AT END OF EPIDEMIC
      mutate(
        model_cum_confirmed = n_confirmed_cases,
        model_cum_confirmed_pc = n_confirmed_cases / n_pop,

        model_new_confirmed = n_confirmed_cases - lag(n_confirmed_cases),
        model_new_confirmed_pc = model_new_confirmed / n_pop,
        model_new_confirmed_week_pc = model_cum_confirmed_pc - lag(model_cum_confirmed_pc, lag_period),
        
        model_cum = n_cases_cum,
        model_cum_pc = n_cases_cum / n_pop,
        model_new_week_pc = (model_cum - lag(model_cum, lag_period)) / n_pop
      ) %>% 
      # Calculate share of cases from each group
      group_by(model_type, sim_id, t) %>% 
      mutate(
        detected_ratio = n_confirmed_cases / model_cum,
        cum_group_share = model_cum_confirmed / sum(model_cum_confirmed),
        cum_group_share_actual = model_cum / sum(model_cum)
      ) %>% 
      ungroup %>% 
      select(-n_confirmed_cases, -n_cases_cum) %>% 
      mutate(date_model = days(t) + lower_bound_date - days(lower_bound_t))
    
    
    # Summarise over multiple sims [still by group]
    model_by_group_summ <- model_by_group %>% 
      group_by(model_type, t, date_model, i_group) %>% 
      
      quantile_summarise(c(model_new_confirmed_week_pc, model_cum_confirmed_pc, model_new_week_pc, model_cum_pc, cum_group_share, cum_group_share_actual, detected_ratio), conf_level = conf_level) %>% 
      group_by(model_type, i_group) %>% 
      # mutate(
      #   across(starts_with("model_new_confirmed_pc_"), 
      #          list(roll = ~ zoo::rollmean(.x, k = 7, fill = NA_real_, align = "center")))
      # ) %>% 
      ungroup
    
    
    # Rename vars to be consistent with data
    model_renamed <- model_by_group_summ %>% 
      rename_with(~ str_replace_all(., "^model_", "n_")) %>% 
      rename(model_type = n_type) %>% 
      rename_with(~ str_replace_all(.x, "_mean", "")) %>% 
      rename(date = date_model) %>% 
      glimpse()
    
    
    # Rename data vars to be consistent with model
    data_renamed <- confirmed_cases_by_group %>% 
      group_by(i_group) %>% 
      mutate(n_new_confirmed_week_pc = n_confirmed_cum_pc - lag(n_confirmed_cum_pc, lag_period)) %>% 
      ungroup %>% 
      rename(
        n_new_confirmed = n_confirmed_new,
        n_new_confirmed_pc = n_confirmed_new_pc,
        n_cum_confirmed_pc = n_confirmed_cum_pc,
        n_cum_confirmed = n_confirmed_cum
      ) %>% 
      group_by(date_results) %>% 
      mutate(
        cum_group_share = n_cum_confirmed / sum(n_cum_confirmed)
      ) %>% 
      ungroup %>% 
      rename(date = date_results) %>% 
      mutate(model_type = "data")

    
    # Combine model and data
    model_by_group_with_data <- bind_rows(
      model_renamed, data_renamed
    ) %>% 
      group_by(model_type, i_group) %>% 
      arrange(model_type, i_group, date) %>% 
      # mutate(
      #   across(c(n_new_confirmed_pc, n_new_confirmed_pc_upper, n_new_confirmed_pc_lower),
      #          list(roll = ~ .x - lag(.x, )))
      # ) %>% 
      mutate(i_group = recode_i_group(i_group)) %>% 
      count_prop(i_group) %>% 
      ungroup %>% 
      glimpse()
      
      
    # Reformat to allow for comparison of both with the data
    model_vs_data_comp <- bind_rows(
      model_by_group_with_data %>% filter(model_type != "mobility_change") %>% 
        mutate(model_vs_data = case_when(model_type == "data" ~ "data", model_type == "baseline" ~ "model"),
               model_type = "baseline"),
      model_by_group_with_data %>% filter(model_type != "baseline") %>% 
        mutate(model_vs_data = case_when(model_type == "data" ~ "data", model_type == "mobility_change" ~ "model"),
               model_type = "mobility_change")
    ) %>% 
      mutate(model_vs_data = factor(model_vs_data, levels = c("model", "data"), labels = c("Model", "Data"))) %>% 
      mutate(model_type = factor(model_type, levels = c("baseline", "mobility_change"),
                                 labels = c("Model (No mobility change)", "Model (Mobility change)")))
    
    
    
    
    

# Merge real/model to get totals by model_type ----------------------------

    # Amalgamate across all groups
    model_all_groups <- model_by_group %>% 
      group_by(model_type, sim_id, t, date_model) %>% 
      summarise(cum_confirmed = sum(model_cum_confirmed),
                cum = sum(model_cum),
                n_pop = sum(n_pop)) %>% 
      group_by(model_type, sim_id) %>% 
      mutate(new_confirmed_week = cum_confirmed - lag(cum_confirmed, lag_period),
             new_confirmed_week_pc = new_confirmed_week / n_pop,
             cum_confirmed_pc = cum_confirmed / n_pop,
             new_week_pc = (cum - lag(cum, lag_period)) / n_pop) %>% 
      group_by(model_type, date_model, t) %>% 
      quantile_summarise(c(new_confirmed_week_pc, cum_confirmed_pc, new_week_pc), conf_level = conf_level) %>% 
      rename(date = date_model) %>% 
      mutate(model_vs_data = "model") %>% 
      rename_with(~ str_replace_all(.x, "_mean", ""))
      
    data_all_groups <- real_cases_slim %>% 
      select(date_results, total_pop, n_confirmed_cum) %>% 
      mutate(
        new_confirmed_week = n_confirmed_cum - lag(n_confirmed_cum, lag_period),
        new_confirmed_week_pc = new_confirmed_week / total_pop,
        cum_confirmed_pc = n_confirmed_cum / total_pop
      ) %>% 
      select(date = date_results, new_confirmed_week_pc, cum_confirmed_pc) %>% 
      mutate(model_type = "data", model_vs_data = "data")
    
    
    all_groups <- bind_rows(
      model_all_groups,
      data_all_groups
      # data_all_groups %>% mutate(model_type = "baseline"),
      # data_all_groups %>% mutate(model_type = "mobility_change")
    ) %>% 
      ungroup %>% 
      mutate(model_vs_data = factor(model_vs_data, levels = c("model", "data"), labels = c("Model", "Data"))) %>% 
      mutate(model_type = factor(model_type, levels = c("data", "baseline", "mobility_change"),
                                 labels = c("Data", "Model (No mobility change)", "Model (Mobility change)")))
    
    
    all_groups_comp <- bind_rows(
      model_all_groups,
      data_all_groups %>% mutate(model_type = "baseline"),
      data_all_groups %>% mutate(model_type = "mobility_change")
    )
    
    

# PLOTS -------------------------------------------------------------------

    date_lim_plot <- as.Date("2021-06-01")
    
    
    # (1) New by group 
    # NEW BY GROUP 
    plot_model_match_detected <- model_by_group_with_data %>% 
      mutate(model_type = factor(model_type, levels = c("data", "baseline", "mobility_change"),
                                 labels = c("(a) Data (HSB)", "(b) Model (No mobility change)", "(c) Model (Mobility change)"))) %>% 
      # mutate(i_group = recode_i_group(i_group)) %>% 
      # filter(t <= max_t) %>%
      ggplot(aes(x = date)) + 
      # geom_ribbon(aes(ymin = model_new_confirmed_pc_lower_roll,
      #                 ymax = model_new_confirmed_pc_upper_roll,
      #                 fill = i_group,
      #                 name = "model_new_confirmed_pc_mean_roll"),
      #             alpha = 0.1, size = 0.3) +
      geom_line(aes(y = n_new_confirmed_week_pc, colour = i_group), size = 1) + 
      facet_wrap( ~ model_type) + 
      scale_colour_viridis(option = "viridis",
                           begin = 0.3, 
                           discrete = TRUE) + 
      scale_fill_viridis(option = "viridis",
                         begin = 0.3, 
                         discrete = TRUE) + 
      theme_custom() + 
      scale_y_continuous(labels = scales::percent) + 
      scale_x_date(date_labels = "%e/%m/%y", breaks = "4 months", minor_breaks = "months", limits = c(as.Date(NA), date_lim_plot)) + 
      # coord_cartesian(xlim = c(NA, date_lim_plot)) + 
      labs(y = "Per capita incidence\n(confirmed cases only,\nprevious 2 weeks)", colour = "SES group", x = "Date of Confirmed Case")
    
    plot_model_match_detected
    
    ggsave("figures/model_match.pdf",   width = 9, height = 3, scale = 1.5)
    
    
    
    model_by_group_with_data %>% 
      filter(model_type == "data") %>% 
      filter(i_group == "5&6") %>% 
      select(date, n_new_confirmed_pc) %>% 
      filter(date <= ymd("2020-08-01")) %>% 
      filter(n_new_confirmed_pc == max(n_new_confirmed_pc, na.rm = TRUE))
    
    
    
    # model_by_group_with_data %>% 
    #   count(value(n_new_confirmed_week_pc_upper), model_type)
    
    # (2) New by group with confidence intervals
    # Transform to get data and model 2x2
    ggplot(model_vs_data_comp, aes(x = date)) + 
      geom_ribbon(aes(ymin = n_new_confirmed_week_pc_lower,
                      ymax = n_new_confirmed_week_pc_upper,
                      fill = i_group,
                      group = model_vs_data),
                  alpha = 0.2, size = 0.3) +
      geom_line(aes(y = n_new_confirmed_week_pc, colour = i_group, linetype = model_vs_data), size = 1) + 
      facet_grid(i_group ~ model_type) + 
      scale_x_date(date_labels = "%e/%m/%y", breaks = "3 months", limits = c(as.Date(NA), date_lim_plot)) + 
      scale_y_continuous(labels = scales::percent) + 
      labs(x = "Date of Confirmed Case", y = "Per capita incidence (confirmed cases only, previous 2 weeks)", linetype = element_blank(), colour = "SES Group", fill = "SES Group") + 
      scale_colour_viridis(option = "viridis",
                           begin = 0.3, 
                           discrete = TRUE) + 
      scale_fill_viridis(option = "viridis",
                         begin = 0.3, 
                         discrete = TRUE) + 
      theme_custom()
    
    ggsave("figures/model_match_conf_int.pdf",   width = 6, height = 6, scale = 1.5)
    
    
    # (3) Cumulative cases
    ggplot(model_vs_data_comp, aes(x = date)) + 
      geom_ribbon(aes(ymin = n_cum_confirmed_pc_lower,
                      ymax = n_cum_confirmed_pc_upper,
                      fill = i_group,
                      group = model_vs_data),
                  alpha = 0.2, size = 0.3) +
      geom_line(aes(y = n_cum_confirmed_pc, colour = i_group, linetype = model_vs_data), size = 1) + 
      facet_grid(i_group ~ model_type) + 
      scale_x_date(date_labels = "%e/%m/%y", breaks = "3 months", limits = c(as.Date(NA), date_lim_plot)) + 
      scale_y_continuous(labels = scales::percent) + 
      labs(x = "Date of Confirmed Case", y = "Cumulative per capita incidence (confirmed cases only)", linetype = element_blank(), colour = "SES Group", fill = "SES Group") + 
      scale_colour_viridis(option = "viridis",
                           begin = 0.3, 
                           discrete = TRUE) + 
      scale_fill_viridis(option = "viridis",
                         begin = 0.3, 
                         discrete = TRUE) + 
      theme_custom()
    
    ggsave("figures/model_match_cumulative.pdf",   width = 6, height = 6, scale = 1.5)
    
    
    
    # (4) Comparing the SHARES of infected in each group
    model_vs_data_comp %>% 
      filter(date <= ymd("2021-02-01") & date >= ymd("2020-05-01")) %>% 
      ggplot(aes(x = date)) + 
      geom_line(aes(y = cum_group_share, colour = i_group, linetype = model_vs_data), size = 1) + 
      geom_ribbon(aes(ymin = cum_group_share_lower,
                      ymax = cum_group_share_upper,
                      fill = i_group,
                      group = model_vs_data),
                  alpha = 0.2, size = 0.3) +
      facet_grid(i_group ~ model_type, scales = "free_y") + 
      scale_colour_viridis(option = "viridis",
                           begin = 0.3, 
                           discrete = TRUE) + 
      scale_fill_viridis(option = "viridis",
                         begin = 0.3, 
                         discrete = TRUE) + 
      theme_custom() + 
      scale_x_date(date_labels = "%e/%m/%y", breaks = "3 months") + 
      scale_y_continuous(labels = scales::percent_format(accuracy = 1)) + 
      labs(y = "Share of cumulative confirmed cases from SES group", x = "Date of Confirmed Case", colour = "SES Group", fill = "SES Group", linetype = element_blank())
    
    
    ggsave("figures/model_match_shares.pdf",   width = 6, height = 6, scale = 1.5)
    
    
    
    # (5) NEW (all groups) (TOTAL)
    ggplot(all_groups, aes(x = date)) +
      geom_line(aes(y = new_confirmed_week_pc, colour = model_vs_data), size = 1) + 
      facet_wrap(~ model_type) + 
      theme_custom(legend.position = "bottom") + 
      geom_ribbon(aes(ymin = new_confirmed_week_pc_lower, ymax = new_confirmed_week_pc_upper, fill = model_vs_data), alpha = 0.2) +
      scale_x_date(date_labels = "%e/%m/%y", minor_breaks = "months", limits = c(as.Date(NA), date_lim_plot)) + 
      labs(y = "Per capita incidence\n(confirmed cases only, previous 2 weeks)", x = "Date of Confirmed Case", colour = element_blank(), fill = element_blank()) + 
      scale_y_continuous(labels = scales::percent) + 
      coord_cartesian(ylim = c(NA, 0.012))
      
    ggsave("figures/model_match_total.pdf",   width = 6, height = 2.5, scale = 1.5)
    
    
    # (6) CUMULATIVE (all groups)
    ggplot(all_groups, aes(x = date)) +
      geom_line(aes(y = cum_confirmed_pc, colour = model_vs_data), size = 1) + 
      facet_wrap(~ model_type) + 
      theme_custom(legend.position = "bottom") + 
      geom_ribbon(aes(ymin = cum_confirmed_pc_lower, ymax = cum_confirmed_pc_upper, fill = model_vs_data), alpha = 0.2) +
      scale_x_date(date_labels = "%e/%m/%y", minor_breaks = "months", limits = c(as.Date(NA), date_lim_plot)) + 
      labs(y = "Cumulative per capita confirmed incidence\n(confirmed cases only)", x = "Date of Confirmed Case", colour = element_blank(), fill = element_blank()) + 
      scale_y_continuous(labels = scales::percent)
      # coord_cartesian(ylim = c(NA, 0.012))
    
    
    ggsave("figures/model_match_total_cum.pdf",   width = 6, height = 2.5, scale = 1.5)
    
    
    
    
    # (7) COMPARING detected vs "true" cases [full curve]
    all_groups %>% print_names
    
    all_groups %>% 
      pivot_longer(c(new_week_pc, new_confirmed_week_pc), names_to = "case_type", values_to = "y") %>% 
      mutate(
        case_type = fct_recode(case_type,
                               "Confirmed Cases Only" = "new_confirmed_week_pc",
                               "All Cases" = "new_week_pc"
        )
      ) %>% 
      ggplot(aes(x = date)) + 
      geom_line(aes(y = y, colour = case_type, linetype = case_type), size = 1) + 
      # geom_line(aes(y = new_week_pc), colour = "seagreen4", linetype = "dashed") + 
      # geom_line(aes(y = new_confirmed_week_pc), colour = "indianred") + 
      facet_wrap(~ model_type) + 
      theme_custom(legend.position = "top") + 
      scale_colour_manual(values = c("orange", "#56B4E9")) + 
      # geom_ribbon(aes(ymin = new_confirmed_week_pc_lower, ymax = new_confirmed_week_pc_upper, fill = model_vs_data), alpha = 0.2) +
      scale_x_date(date_labels = "%e/%m/%y", minor_breaks = "months", limits = c(as.Date(NA), date_lim_plot)) + 
      labs(y = "Per capita incidence\n(all groups, previous 2 weeks)", x = "Date of Case", colour = "Case Type", linetype = "Case Type") + 
      scale_y_continuous(labels = scales::percent)
    
    ggsave("figures/model_match_real_vs_detected.pdf",   width = 6, height = 2.5, scale = 1.5)
    
    
    # (8) Comparing detected vs true cases (by group)
    model_vs_data_comp %>% 
      complete(model_vs_data, model_type, date, i_group) %>% 
      filter(model_vs_data != "Data") %>% 
      select(model_type, t, date, i_group, n_new_week_pc, n_new_confirmed_week_pc) %>% 
      pivot_longer(c(n_new_week_pc, n_new_confirmed_week_pc)) %>% 
      mutate(name = fct_recode(name, "Confirmed Cases Only" = "n_new_confirmed_week_pc",
                               "All Cases" = "n_new_week_pc")) %>% 
      
      ggplot(aes(x = date)) + 
      geom_line(aes(y = value, colour = i_group)) + 
      # geom_line(aes(y = n_new_confirmed_week_pc, colour = i_group)) + 
      facet_grid(name ~ model_type, scales = "free_y") + 
      theme_custom(legend.position = "bottom") + 
      scale_colour_viridis(option = "viridis",
                           begin = 0.3, 
                           discrete = TRUE) + 
      scale_fill_viridis(option = "viridis",
                         begin = 0.3, 
                         discrete = TRUE) + 
      # geom_ribbon(aes(ymin = new_confirmed_week_pc_lower, ymax = new_confirmed_week_pc_upper, fill = model_vs_data), alpha = 0.2) +
      scale_x_date(date_labels = "%e/%m/%y", minor_breaks = "months", limits = c(as.Date(NA), date_lim_plot)) +
      labs(y = "Per capita incidence (previous 2 weeks)", x = "Date of Case", colour = "SES Group") + 
      scale_y_continuous(labels = scales::percent)
    
    ggsave("figures/model_match_real_vs_detected_by_group.pdf",   width = 6, height = 3.5, scale = 1.5)
    
    
    # (9) Summary - detected ratio
    plot_model_match_prop_detected <- model_vs_data_comp %>% 
      complete(model_vs_data, model_type, date, i_group) %>% 
      filter(model_vs_data != "Data") %>% 
      filter(date >= ymd("2020-08-01")) %>% 
      
      ggplot(aes(x = date)) + 
      geom_line(aes(y = detected_ratio, colour = i_group), size = 1) + 
      geom_ribbon(aes(ymin = detected_ratio_lower, ymax = detected_ratio_upper, fill = i_group, colour = i_group), alpha = 0.2, linetype = "dotted") +
      facet_wrap(~ model_type) + 
      theme_custom() + 
      scale_colour_viridis(option = "viridis",
                           begin = 0.3, 
                           discrete = TRUE) + 
      scale_fill_viridis(option = "viridis",
                         begin = 0.3, 
                         discrete = TRUE) + 
      # geom_ribbon(aes(ymin = new_confirmed_week_pc_lower, ymax = new_confirmed_week_pc_upper, fill = model_vs_data), alpha = 0.2) +
      scale_x_date(date_labels = "%e/%m/%y", minor_breaks = "months", limits = c(as.Date(NA), date_lim_plot)) + 
      labs(y = "Proportion of total infections detected", x = "Date", colour = "SES Group", fill = "SES Group") + 
      scale_y_continuous(labels = scales::percent)
    
    plot_model_match_prop_detected

    ggsave("figures/model_match_prop_detected.pdf",   width = 6, height = 2.5, scale = 1.5)    
    
    
    plot_model_match_prop_detected + geom_label_repel(data = model_vs_data_comp %>% 
                                                  complete(model_vs_data, model_type, date, i_group) %>% 
                                                  filter(model_vs_data != "Data") %>% 
                                                  filter(date == date_lim_plot),
                                                aes(y = detected_ratio, label = round(detected_ratio, 3)))
    
    
    ggsave("figures/model_match_prop_detected_vallabels.pdf",   width = 6, height = 2.5, scale = 1.5)    
    
    
    
    
    
    
    
    
    
    
    
    

# DEATHS (temp) -----------------------------------------------------------

    
    model_renamed %>% print_names
    
    # Rename data vars to be consistent with model
    death_data_renamed <- deaths_by_group %>% 
      group_by(i_group) %>% 
      mutate(n_new_week_pc = n_inferred_cum_pc - lag(n_inferred_cum_pc, lag_period)) %>% 
      ungroup %>% 
      # rename(
      #   # n_new_confirmed = n_confirmed_new,
      #   # n_new_confirmed_pc = n_confirmed_new_pc,
      #   # n_cum_confirmed_pc = n_confirmed_cum_pc,
      #   # n_cum_confirmed = n_confirmed_cum
      # ) %>% 
      # group_by(date_results) %>% 
      # mutate(
      #   cum_group_share = n_cum_confirmed / sum(n_cum_confirmed)
      # ) %>% 
      # ungroup %>% 
      rename(date = date_death) %>% 
      mutate(model_type = "data")
    
    # Combine model and data
    model_by_group_with_death_data <- bind_rows(
      model_renamed, death_data_renamed
    ) %>% 
      group_by(model_type, i_group) %>% 
      arrange(model_type, i_group, date) %>% 
      # mutate(
      #   across(c(n_new_confirmed_pc, n_new_confirmed_pc_upper, n_new_confirmed_pc_lower),
      #          list(roll = ~ .x - lag(.x, )))
      # ) %>% 
      mutate(i_group = recode_i_group(i_group)) %>% 
      count_prop(i_group) %>% 
      ungroup %>% 
      glimpse()
    
    
    # Reformat to allow for comparison of both with the data
    model_vs_data_with_deaths_comp <- bind_rows(
      model_by_group_with_death_data %>% filter(model_type != "mobility_change") %>% 
        mutate(model_vs_data = case_when(model_type == "data" ~ "data", model_type == "baseline" ~ "model"),
               model_type = "baseline"),
      model_by_group_with_death_data %>% filter(model_type != "baseline") %>% 
        mutate(model_vs_data = case_when(model_type == "data" ~ "data", model_type == "mobility_change" ~ "model"),
               model_type = "mobility_change")
    ) %>% 
      mutate(model_vs_data = factor(model_vs_data, levels = c("model", "data"), labels = c("Model", "Data"))) %>% 
      mutate(model_type = factor(model_type, levels = c("baseline", "mobility_change"),
                                 labels = c("Model (No mobility change)", "Model (Mobility change)")))
    
    
    
    
    # (1) No confidence intervals
    plot_model_match_death <- model_by_group_with_death_data %>% 
      mutate(model_type = factor(model_type, levels = c("data", "baseline", "mobility_change"),
                                 labels = c("(a) Data (HSB)", "(b) Model (No mobility change)", "(c) Model (Mobility change)"))) %>% 
      # mutate(i_group = recode_i_group(i_group)) %>% 
      # filter(t <= max_t) %>%
      ggplot(aes(x = date)) + 
      # geom_ribbon(aes(ymin = model_new_confirmed_pc_lower_roll,
      #                 ymax = model_new_confirmed_pc_upper_roll,
      #                 fill = i_group,
      #                 name = "model_new_confirmed_pc_mean_roll"),
      #             alpha = 0.1, size = 0.3) +
      geom_line(aes(y = n_new_week_pc, colour = i_group), size = 1) + 
      facet_wrap( ~ model_type) + 
      scale_colour_viridis(option = "viridis",
                           begin = 0.3, 
                           discrete = TRUE) + 
      scale_fill_viridis(option = "viridis",
                         begin = 0.3, 
                         discrete = TRUE) + 
      theme_custom() + 
      theme(legend.position = "top") + 
      scale_y_continuous(labels = scales::percent) + 
      scale_x_date(date_labels = "%e/%m/%y", breaks = "4 months", minor_breaks = "months", limits = c(as.Date(NA), date_lim_plot)) + 
      # coord_cartesian(xlim = c(NA, date_lim_plot)) + 
      labs(y = "Per capita incidence\n(previous 2 weeks)", x = "Date", colour = "SES group")
    
    plot_model_match_death
    
    ggsave("figures/model_match_deaths.pdf",   width = 9, height = 3, scale = 1.5)
    
    # (2) New by group with confidence intervals
    # Transform to get data and model 2x2
    ggplot(model_vs_data_with_deaths_comp, aes(x = date)) + 
      geom_ribbon(aes(ymin = n_new_week_pc_lower,
                      ymax = n_new_week_pc_upper,
                      fill = i_group,
                      group = model_vs_data),
                  alpha = 0.2, size = 0.3) +
      geom_line(aes(y = n_new_week_pc, colour = i_group, linetype = model_vs_data), size = 1) + 
      facet_grid(i_group ~ model_type) + 
      scale_x_date(date_labels = "%e/%m/%y", breaks = "3 months", limits = c(as.Date(NA), date_lim_plot)) + 
      scale_y_continuous(labels = scales::percent_format(accuracy = 1)) + 
      labs(x = "Date", y = "Per capita incidence\n(previous 2 weeks)", linetype = element_blank(), colour = "SES Group", fill = "SES Group") + 
      scale_colour_viridis(option = "viridis",
                           begin = 0.3, 
                           discrete = TRUE) + 
      scale_fill_viridis(option = "viridis",
                         begin = 0.3, 
                         discrete = TRUE) + 
      theme_custom()
    
    ggsave("figures/model_match_deaths_conf_int.pdf",   width = 6, height = 6, scale = 1.5)
    
  
    
    

# Match with COVIDA data --------------------------------------------------

    # FROM import_covida_data.R
    covida_rates_renamed <- covida_rates %>% 
      mutate(
        date = case_when(grp == 1 ~ ymd("2020-11-30"),
                         grp == 2 ~ ymd("2021-03-03"))
      ) %>% 
      print_names %>% 
      rename(
        n_cum_pc = acumm_covid_covida,
        n_cum_pc_lower = q025_acumm_covid_covida,
        n_cum_pc_upper = q975_acumm_covid_covida
      ) %>% 
      mutate(i_group = recode_i_group(stratum)) %>% 
      mutate(model_type = "covida")
    
    # RATIO OF INFECTIONS
    covida_rates_renamed %>% 
      filter(date == ymd("2021-03-3")) %>% 
      select(stratum, n_cum_pc) %>% 
      mutate(n_cum_pc_ratio = n_cum_pc / n_cum_pc[stratum == 4])
    
      
    # Incorporate model data
    model_match_covida <- model_by_group_with_data %>% 
      glimpse() %>% 
      filter(date %in% c(ymd("2020-11-30"), ymd("2021-03-03"))) %>% 
      count_prop(model_type) %>% 
      filter(model_type != "data") %>% 
      select(model_type, t, date, i_group, matches("n_cum_pc")) %>% 
      bind_rows(covida_rates_renamed) %>% 
      arrange(date)
  
      
    
    plot_model_match_covida <- model_match_covida %>% 
      # filter(! (model_type == "baseline" & date == ymd("2021-03-03"))) %>% 
      filter(date == ymd("2021-03-3")) %>% 
      mutate(
        model_type = fct_recode(
          factor(model_type, levels = c("covida", "baseline", "mobility_change")),
          "(d) Data (CoVIDA)" = "covida",
          "(e) Model (No mobility change)" = "baseline",
          "(f) Model (Mobility change)" = "mobility_change"
        )
      ) %>% 
      ggplot(aes(x = i_group, colour = i_group)) + 
      geom_point(aes(y = n_cum_pc), show.legend = FALSE, position = position_dodge(width = 0.2)) + 
      geom_errorbar(aes(ymin = n_cum_pc_lower, ymax = n_cum_pc_upper), show.legend = FALSE,  width = 0.2, position = position_dodge(width = 0.2)) + 
      facet_grid(~ model_type) +
      scale_y_continuous(labels = scales::percent) + 
      coord_cartesian(ylim = c(0, 1)) + 
      scale_colour_viridis(option = "viridis",
                           begin = 0.3, end = 0.98, 
                           discrete = TRUE) + 
      labs(x = "SES Group", y = "Cumulative per capita\nincidence\n(by 03/03/2021)", colour = element_blank()) + 
      theme_custom(panel.grid.major = element_blank())
    
    plot_model_match_covida
      
    ggsave("figures/model_match_covida.pdf",   width = 6, height = 3, scale = 0.8)

        
    plot_model_match_covida + geom_label_repel(aes(y = n_cum_pc, label = round(n_cum_pc, 3)))
    ggsave("figures/model_match_covida_vallabels.pdf",   width = 6, height = 3, scale = 1.5)
    
      
    
    

# Deaths and COVIDA -------------------------------------------------------

    library("ggpubr")
    
    
    
    plot_final <- ggarrange(
      plot_model_match_death,
      plot_model_match_covida,
      nrow = 2,
      heights = c(2, 2),
      align = "v",
      common.legend = TRUE
    )
    # ?ggarrange
    plot_final
    
    ggsave("figures/model_match_covida_and_deaths.pdf",   width = 7, height = 4.5, scale = 1.2)
     
    
    library("ggpubr")
    
    
    
    plot_final_detected <- ggarrange(
      plot_model_match_detected,
      plot_model_match_covida,
      nrow = 2,
      heights = c(2, 2),
      align = "v",
      common.legend = TRUE
    )
    # ?ggarrange
    plot_final_detected
    
    ggsave("figures/model_match_covida_detected.pdf", width = 7, height = 4.5, scale = 1.2)
       
    
    
    

# PLOT DATA ---------------------------------------------------------------


    
    # PLOT DATA FOR THIS FIGURE
    
    # TOP HALF 
    model_by_group_with_data %>% count_prop(model_type) %>% 
      mutate(model_type = factor(model_type, levels = c("data", "baseline", "mobility_change"),
                                 labels = c("(a) Data (HSB)", "(b) Model (No mobility change)", "(c) Model (Mobility change)"))) %>% 
      select(model_type, t, date, i_group, n_new_confirmed_week_pc) %>% 
      write_excel_csv("data/processed/plot_data_fig2a.csv")
    
    # BOTTOM HALF
    model_match_covida %>% 
      filter(date == ymd("2021-03-3")) %>% 
      mutate(
        model_type = fct_recode(
          factor(model_type, levels = c("covida", "baseline", "mobility_change")),
          "(d) Data (CoVIDA)" = "covida",
          "(e) Model (No mobility change)" = "baseline",
          "(f) Model (Mobility change)" = "mobility_change"
        )
      ) %>% 
      select(model_type, ses_group = i_group, n_cum_pc, n_cum_pc_lower, n_cum_pc_upper)
    
    model_by_group_with_data %>% count_prop(model_type) %>% 
      mutate(model_type = factor(model_type, levels = c("data", "baseline", "mobility_change"),
                                 labels = c("(a) Data (HSB)", "(b) Model (No mobility change)", "(c) Model (Mobility change)"))) %>% 
      select(model_type, t, date, i_group, n_new_confirmed_week_pc) %>% 
      write_excel_csv("data/processed/plot_data_fig2b.csv")
    

# Share vs COVIDA ---------------------------------------------------------

    
    monte_carlo_conf_int <- function(means, sds, n_sims = 1000) {
      
      # means <- c(1, 2, 3, 4)
      # sds <- c(0.3, 0.2, 0.1, 0.5)
      # n_sims = 1000
      
      tibble(
        group = 1:length(means),
        mean = means,
        sd = sds
      ) %>% 
        rowwise() %>% 
        mutate(sim_val = list(rnorm(n = n_sims, mean = mean, sd = sd)),
               id = list(1:n_sims)) %>% 
        unnest(c(sim_val, id)) %>% 
        group_by(id) %>%
        mutate(ratio = sim_val / sum(sim_val)) %>% 
        ungroup %>% 
        group_by(group, mean, sd) %>% 
        summarise(
          conf_lower = quantile(ratio, 0.025),
          conf_upper = quantile(ratio, 0.975)
        ) %>% 
        select(-group)
      
    }
    
    covida_group_share <- covida_rates_renamed %>% 
      select(date, i_group, matches("tot_day")) %>% 
      # mutate(var_upper = (q975_tot_day_cases_covida - tot_day_cases_covida) / 1.96,
             # var_lower = abs(q025_tot_day_cases_covida - tot_day_cases_covida) / 1.96) %>% 
      mutate(sd = (q975_tot_day_cases_covida - q025_tot_day_cases_covida) / (2 * qnorm(0.975))) %>% 
      group_by(date) %>% 
      summarise(monte_carlo_conf_int(means = tot_day_cases_covida, sds = sd, n_sims = 10000)) %>% 
      mutate(i_group = recode_i_group(group)) %>% 
      mutate(val = mean / sum(mean)) %>% 
      print
      # group_by(date) %>% 
      # mutate(cum_group_share_covida = tot_day_cases_covida / sum(tot_day_cases_covida))
      
      
      

    
      
      
      
    
    
  
    # (4) Comparing the SHARES of infected in each group
    model_by_group_with_data %>% 
      filter(model_type != "data") %>% 
      filter(date <= ymd("2021-03-01") & date >= ymd("2020-05-01")) %>% 
      mutate(model_type = factor(model_type, levels = c("baseline", "mobility_change"),
                                 labels = c("Model (No mobility change)", "Model (Mobility change)"))) %>% 
      ggplot(aes(x = date)) + 
      geom_line(aes(y = cum_group_share_actual, colour = i_group)) + 
      geom_ribbon(aes(ymin = cum_group_share_actual_lower,
                      ymax = cum_group_share_actual_upper,
                      fill = i_group),
                  alpha = 0.2, size = 0.3) +
      geom_point(data = covida_group_share, 
                 aes(x = date, y = val, colour = i_group), size = 0.8) + 
      geom_errorbar(data = covida_group_share, 
                    aes(x = date, ymin = conf_lower, ymax = conf_upper, colour = i_group),
                    width = 6) + 
      facet_grid(i_group ~ model_type, scales = "free_y") + 
      scale_colour_viridis(option = "viridis",
                           begin = 0.3, 
                           discrete = TRUE) + 
      scale_fill_viridis(option = "viridis",
                         begin = 0.3, 
                         discrete = TRUE) + 
      theme_custom(legend.position = "top") + 
      scale_x_date(date_labels = "%e/%m/%y", breaks = "3 months") + 
      scale_y_continuous(labels = scales::percent_format(accuracy = 1)) + 
      labs(y = "Share of cumulative cases from SES group", x = "Date of case", colour = "SES Group", fill = "SES Group", linetype = element_blank())
    
    
    ggsave("figures/model_match_shares_covida.pdf",   width = 6, height = 6, scale = 0.8)
    
    
    
    