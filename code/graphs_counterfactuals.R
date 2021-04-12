# COUNTERFACTUALs ---------------------------------------------------------
    
    load("data/processed/group_diff_slim.RData", verbose = TRUE)
    
    # Prepare plot DF
    group_diff_at_t <- group_diff_slim %>% 
      # mutate(parameter_set = paste0(parameter_set, "_", tighten_factor)) %>% 
      group_by(parameter_set, sim_id, i_group) %>% 
      arrange(parameter_set, sim_id, i_group, t) %>% 
      mutate(new_cases_2_week = n_cases_cum - lag(n_cases_cum, 14)) %>% 
      group_by(parameter_set, t, i_group, n_pop) %>% 
      quantile_summarise(c(n_cases_live, new_cases_2_week), conf_level = 0) %>% 
      mutate(incidence_prop = new_cases_2_week_median / n_pop) %>% 
      ungroup
      # mutate(prevalence_prop = n_cases_live_median / n_pop)
    
    # Labels for each counterfacutal type
    counterfactual_labels <- c(
      "baseline" = "Baseline",
      "out_contacts" = "Out of home",
      "home_contacts" = "Within home",
      "isolation_behaviour" = "Isolation behavior",
      "testing_tracing" = "Testing & Tracing",
      "all" = "All"
    )
    
    
# 1. Counterfactual epidemic curve -----------------------------------------------------------
    
    
    
    # PLOT epidemic curves
    # included_curves <- c("baseline", "out_contacts", "home_contacts", "isolation_behaviour", "testing_tracing", "all")
    # group_diff_at_t %>% ungroup %>% count_prop(parameter_set)
    
    group_diff_at_t %>% 
      ungroup %>% 
      mutate(i_group = recode_i_group(i_group)) %>% 
      # mutate(alpha = if_else(parameter_set %in% included_curves, 1, 0)) %>% 
      ggplot(aes(x = t, y = incidence_prop, 
                 colour = factor(i_group), group = factor(i_group)
                 # alpha = 1
      )) + 
      geom_line(show.legend = TRUE) + 
      facet_wrap( ~ parameter_set) + 
      facet_wrap(~ parameter_set, labeller = as_labeller(counterfactual_labels)) +
      
      # Formatting
      theme_custom() + 
      theme(legend.position = "right") + 
      # scale_alpha_continuous(guide = FALSE) + 
      scale_y_continuous(labels = scales::percent_format(accuracy = 1), minor_breaks = NULL) + 
      labs(colour = "SES Group", fill = "SES Group", y = "Per capita incidence (previous 2 weeks)") +
      geom_hline(yintercept = 0, colour = "darkgrey", size = 0.4) + 
      # coord_cartesian(xlim = c(0, 500)) +
      scale_colour_viridis(option = "viridis",
                           begin = 0.3, 
                           discrete = TRUE) + 
      coord_cartesian(xlim = c(0, 300))
    # scale_x_continuous(minor_breaks = c(0, 25, 50, 75), breaks = c(0, 50, 100))
    
    
    ggsave("figures/counterfactuals_epidemic_curve.png", width = 8, height = 4.5, scale = 0.8)
    # ggsave("epidemic_curve_2.png", width = 8, height = 4.5, scale = 0.8)
    # ggsave("epidemic_curve_3.png", width = 8, height = 4.5, scale = 0.8)
    
    
    
    
    
    
    
    # 2. Counterfactual - summary results -------------------------------------
    
    # group_diff_slim %>% count_prop(parameter_set, tighten_factor)
    
    # SUMMARISE ALL INFO
    group_diff_all_summ <- group_diff_slim %>% 
      # mutate(parameter_set = paste0(parameter_set, "_", tighten_factor)) %>% 
      
      # (1) Calculate total cases
      group_by(parameter_set, i_group, sim_id, n_pop) %>% 
      summarise(n_cases = max(n_cases_cum)) %>% 
      
      # Calculate proportion infected, and gaps with group 5
      group_by(parameter_set, sim_id) %>% arrange(sim_id) %>% 
      mutate(prop_infected = n_cases / n_pop) %>% 
      mutate(diff_prop_infected = prop_infected - prop_infected[i_group == 4],
             diff_n_cases = n_cases - n_cases[i_group == 4]) %>% 
      
      
      # Calculate "effect" relative to baseline
      group_by(sim_id, i_group) %>% 
      mutate(across(-c(parameter_set), list(effect = ~ .x - .x[parameter_set == "baseline"]))) %>% 
      # Calculte implied contribution for each group, for each sim
      mutate(contribution = diff_prop_infected_effect / diff_prop_infected_effect[parameter_set == "all"]) %>% 
      group_by(parameter_set, i_group) %>% 
      arrange(parameter_set, i_group)
    # summarise_quantiles(-c(sim_id)) %>% print_all
    
    
    
    # Make into function
    plot_group_diffs <- function(.data, y, ylab = "Y") {
      
      # Labels for each counterfacutal type
      counterfactual_labels <- c(
        "baseline" = "Baseline",
        "out_contacts" = "Out of home",
        "home_contacts" = "Within home",
        "isolation_behaviour" = "Isolation behavior",
        "testing_tracing" = "Testing & Tracing",
        "all" = "All"
      )
      
      .data %>% 
        mutate(i_group = recode_i_group(i_group)) %>% 
        summarise_quantiles({{y}}) %>% 
        select(parameter_set, i_group, q, {{y}}) %>% 
        filter(q %in% c(0.025, 0.5, 0.975)) %>% 
        pivot_wider(names_from = q, names_prefix = "y_", values_from = {{y}}) %>% 
        ggplot(aes(x = i_group, fill = factor(i_group))) + 
        geom_col(aes(y = y_0.5)) + 
        geom_errorbar(aes(ymin = y_0.025, ymax = y_0.975), width = 0.2, size = 0.5, colour = "#374561") + 
        facet_wrap(~ parameter_set, labeller = as_labeller(counterfactual_labels)) +
        # facet_wrap(~ parameter_set) + 
        labs(y = ylab, x = "SES Group", fill = "SES Group") + 
        theme_custom() + 
        theme(legend.position = "right", panel.grid.minor = element_blank(), panel.grid.major.x = element_blank()) + 
        geom_hline(yintercept = 0, colour = "darkgrey", size = 0.4) + 
        scale_fill_viridis(option = "viridis",
                           begin = 0.3, 
                           discrete = TRUE)
      # scale_x_continuous(breaks = 1:5)
      
    }
    
    
    # PLOT
    # Prop infected in each scenario
    plot_group_diffs(group_diff_all_summ, y = prop_infected, "Cumulative per capita incidence\n(entire epidemic)")
    ggsave("figures/prop_infected.png", width = 8, height = 4.5, scale = 0.8)
    
    plot_group_diffs(group_diff_all_summ, y = prop_infected, "Cumulative per capita incidence\n(entire epidemic)") + geom_label_repel(aes(y = y_0.5, label = round(y_0.5, 3)))
    ggsave("figures/prop_infected_vallabels.png", width = 8, height = 4.5, scale = 0.8)
    
    
    # Gap (relative to group 6) in each scenario
    plot_group_diffs(group_diff_all_summ, y = diff_prop_infected, "Gap relative to group 6")
    ggsave("figures/inequality.png", width = 8, height = 4.5, scale = 0.8)
    
    plot_group_diffs(group_diff_all_summ, y = diff_prop_infected, "Gap relative to group 6") + geom_label_repel(aes(y = y_0.5, label = round(y_0.5, 3)))
    ggsave("figures/inequality_valllabels.png", width = 8, height = 4.5, scale = 0.8)
    
    # TO DO !!!! SHOULD CALCULATE THE QUANTILES AFTER CALCULATING DESIRED QUANTITY
    
    # CHANGE RELATIVE TO BASELINE of prop infected
    # group_diff_all_summ %>% 
    #   group_by(i_group, q) %>% 
    #   mutate(across(-c(parameter_set), ~ . - .[parameter_set == "baseline"])) %>% 
    #   ungroup %>% 
    #   
    #   plot_group_diffs(
    #     y = prop_infected,
    #     ylab = "Change in Total Proportion Infected relative to Baseline"
    #   )
    
    # CHANGE RELATIVE TO BASELINE of group 1 gap (inequality)
    # group_diff_all_summ %>% 
    #   group_by(i_group, q) %>% 
    #   mutate(across(-c(parameter_set), ~ . - .[parameter_set == "baseline"])) %>% 
    #   ungroup %>% 
    #   
    #   plot_group_diffs(
    #     y = diff_prop_infected,
    #     ylab = "Change relative to baseline case in inequality relative to group 6"
    #   )
    
    
    
    # 3. Decomposition (v1) ----------------------------------------------------
    
    
    
    # colours <- c(
    #   "darkgrey", "#E87D72", "#55BA70", "#4FB7DF", "#D378EF" 
    # )
    
    
    # PLOT implied "proportion" contribution of each parameter group
    # group_diff_all_summ %>% 
    #   group_by(i_group, q) %>% 
    #   mutate(across(-c(parameter_set), ~ . - .[parameter_set == "baseline"])) %>% 
    #   ungroup %>% 
    #   filter(q == 0.5) %>% 
    #   filter(parameter_set != "baseline") %>% 
    #   group_by(i_group) %>% 
    #   filter(i_group != 4) %>% 
    #   arrange(i_group) %>% 
    #   mutate(diff_prop_infected = diff_prop_infected / diff_prop_infected[parameter_set == "all"]) %>% 
    #   select(parameter_set, i_group, diff_prop_infected) %>% 
    #   mutate(parameter_set = as.character(parameter_set)) %>%
    #   mutate(parameter_set = if_else(parameter_set == "all", "Interaction", parameter_set)) %>%
    #   mutate(parameter_set = factor(parameter_set)) %>%
    #   group_by(i_group) %>%
    #   mutate(diff_prop_infected = if_else(parameter_set == "Interaction",
    #                                       1 - sum(diff_prop_infected[parameter_set != "Interaction"]),
    #                                       diff_prop_infected)) %>%
    #   # filter(!(parameter_set == "Interaction" & diff_prop_infected < 0)) %>% 
    #   mutate(parameter_set = forcats::fct_relevel(parameter_set, "Interaction"),
    #          parameter_set = fct_relevel(parameter_set, "home_contacts", after = 6)) %>% 
    #   
    #   ggplot(aes(x = fct_rev(factor(i_group)), y = diff_prop_infected)) + 
    #   # ggpattern::geom_col_pattern(
    #   #   aes(fill = parameter_set, pattern = parameter_set, pattern_angle = class), 
    #   #   colour = "white"
    #   # ) + 
    #   geom_col(aes(fill = parameter_set), colour = "white", width = 0.4) + 
    #   # geom_label_repel(aes(label = parameter_set), size = 0.
    #   #                  position = position_stack(vjust = .5)) + 
    #   scale_fill_manual(values = colours, labels = counterfactual_labels, name = element_blank()) + 
    #   theme_custom() + 
    #   geom_hline(yintercept = 1, size = 0.3, linetype = "longdash") + 
    #   geom_hline(yintercept = 0, size = 0.3, linetype = "longdash") + 
    #   geom_label(data = . %>% mutate(label = if_else(diff_prop_infected > 0.1, 
    #                                                  paste0(round(diff_prop_infected * 100), "%"),
    #                                                  NA_character_)),
    #              aes(label = label, group = parameter_set),
    #              position = position_stack(vjust = 0.5),
    #              alpha = 0.5,
    #              size = 2) + 
    #   coord_flip() + 
    #   scale_y_continuous(labels = scales::percent_format(accuracy = 1), breaks = 0:5 / 5) + 
    #   # scale_fill_wsj() + 
    #   labs(x = "SES Group", y = "Contribution to reduction in inequality")
    
    
    
    # ggsave("inequality_decomposition.png",  width = 6, height = 4.5, scale = 0.8)
    
    
    
    
    
    
# 4. Decoposition (v2) ----------------------------------------------------
    
    
    
    
    # V2 - no pie chart, just bars
    # Relative to total inequality, how much does each parameter change reduce
    
    
    contributions_df <- group_diff_all_summ %>% 
      # mutate(parameter_set = paste0(parameter_set, "_", tighten_factor)) %>% 
      # group_by(sim_id, parameter_set,  i_group) %>% 
      arrange(sim_id, parameter_set,  i_group) %>% 
      select(sim_id, parameter_set,  i_group, diff_n_cases) %>% 
      
      # For each sim_id, i_group, what's the prop reduction
      group_by(sim_id, i_group) %>% 
      filter(i_group != 4) %>%  # remove 4 because we're discussing related to top group
      mutate(prop_reduction = (diff_n_cases[parameter_set == "baseline"] - diff_n_cases) / diff_n_cases[parameter_set == "baseline"]) %>% 
      
      # Calculate quantiles across sim_ids
      group_by(i_group, parameter_set) %>% 
      summarise_quantiles(prop_reduction) %>% 
      filter(q %in% c(0.05, 0.5, 0.95)) %>% 
      mutate(q_label = case_when(q == 0.05 ~ "_lower",
                                 q == 0.5   ~ "",
                                 q == 0.95 ~ "_upper")) %>% 
      select(-q) %>%
      ungroup %>% 
      pivot_wider(names_from = q_label,
                  values_from = prop_reduction,
                  names_prefix = "prop_reduction")
    
    # Plot
    plot_inequality_decomposition <- contributions_df %>% 
      mutate(i_group = recode_i_group(i_group)) %>% 
      ggplot(aes(x = parameter_set, y = prop_reduction, fill = i_group)) + 
      geom_col(position = "dodge", width = 0.7) + 
      geom_errorbar(aes(ymin = prop_reduction_lower, ymax = prop_reduction_upper), position = position_dodge(0.7), width = 0.2, size = 0.5, colour = "#374561") + 
      labs(x = "Scenario", y = "% reduction in inequality\nrelative to group 5&6") + 
      scale_y_continuous(labels = scales::percent_format(accuracy = 1)) + 
      scale_fill_viridis(option = "viridis",
                         begin = 0.3, 
                         end = 0.766666667,
                         discrete = TRUE,
                         name = "SES Group") + 
      scale_x_discrete(labels = counterfactual_labels) +
      theme_custom() +
      theme(axis.text.x = element_text(angle = 45, hjust = 1),
            panel.grid.minor = element_blank(),
            panel.grid.major.x = element_blank()) + 
      coord_cartesian(ylim = c(-0.1, 1)) + 
      geom_hline(yintercept = 0)
    
    
    plot_inequality_decomposition
    ggsave("figures/inequality_decomposition.png",  width = 6, height = 4.5, scale = 0.8)
    
    
    plot_inequality_decomposition + geom_label_repel(aes(y = prop_reduction, label = round(prop_reduction, 3)), position = position_dodge(0.7))
    ggsave("figures/inequality_decomposition_vallabels.png",  width = 6, height = 4.5, scale = 0.8)
        
# Output tables ------------------------------------------------------------------
    
    # Column 0 - N infected in each one
    # Column 1 - what proportion of population in each group are infected in each one
    # Column 2 - implied effect relative to baseline on inequality
    # Column 3 - implied contribution of effect on inequality
    
    # Should report the median, the 95th centiles
    # 
    
    
    # group_diff_table <- group_diff_all_summ %>% 
    #   pivot_longer(-c(q, n_pop, parameter_set, i_group), names_to = "outcome_name", values_to = "outcome_val") %>% 
    #   filter(q %in% c(0.01, 0.5, 0.99)) %>% 
    #   
    #   # Combine confidence intervals
    #   pivot_wider(names_from = q, values_from = outcome_val) %>% 
    #   rename(lower =  `0.01`, upper = `0.99`, median = `0.5`) %>% 
    #   
    #   # Keep only outcomes we want in the table
    #   filter(outcome_name %in% c("prop_infected", "diff_prop_infected_effect", "contribution")) %>% 
    #   
    #   # Clean up vars
    #   mutate(
    #     across(c(lower, upper), ~ replace_with_na(.x, c(0, 1))),
    #     across(c(lower, upper), ~ replace_with_na(.x, c(0))),
    #     across(c(lower, median, upper), ~ round(.x, 2)),  # Round all to 3 d.p.
    #     across(c(lower, upper), ~ format(.x, digits = 2, nsmall = 2, trim = TRUE)), 
    #     median = format(median, digits = 2, nsmall = 2),
    #     across(c(median), ~ if_else(str_detect(.x, "NA"), "", .x))
    #   ) %>% # remove 0s/1s
    #   
    #   # Combine lower / upper
    #   mutate(conf = if_else(str_detect(lower, "NA"), "", as.character(str_glue("[{lower}, {upper}]")))) %>% 
    #   mutate(output = str_glue("{median} {conf}")) %>% 
    #   select(-lower, -upper, -median, -conf) %>%
    #   select(-n_pop) %>% 
    #   
    #   # Get each variable in its own column 
    #   pivot_wider(names_from = outcome_name, values_from = output) %>% 
    #   mutate(
    #     parameter_set = fct_recode(parameter_set,
    #                                "Baseline" = "baseline",
    #                                "Interactions (Out)" = "out_contacts",
    #                                "Interactions (Home)" = "home_contacts",
    #                                "Isolation" = "isolation_behaviour",
    #                                "Testing & Tracing" = "testing_tracing",
    #                                "All" = "all"
    #     )
    #   ) %>% 
    #   
    #   # Only write parameter set once each time
    #   mutate(parameter_set = as.character(parameter_set)) %>% 
    #   group_by(parameter_set) %>% 
    #   mutate(parameter_set = if_else(parameter_set == lag(parameter_set) & !is.na(lag(parameter_set)), "", parameter_set)) %>% 
    #   
    #   # Rename columns
    #   set_names(
    #     c("Parameter Set", "SES Group", "Proportion Infected", "Effect on group 5 gap", "Contribution")
    #   ) %>% 
    #   
    #   print_all
    # 
    # library("xtable")
    # 
    # print.xtable(
    #   xtable(group_diff_table, caption = "Main Output"),
    #   # type = "html",
    #   caption.placement = "top",
    #   table.placement = "!htbp",
    #   file = "group_diffs.tex",
    #   include.rownames = FALSE,
    #   hline.after = 0:5 * 5
    # )
    # 
    # 
    # 
    # 
    # 
    # ?as.character
    # 
    # format_dbl <- function(x, n_digits) {
    #   as.character(format(round(x, n_digits), nsmall = n_digits))
    # }
    # 
    # format_chr <- function(x) {
    #   str_replace_all(x, "_", " ") %>% 
    #     str_to_title(locale = "en")
    # }
    # 
    # 
    # 
    # out_table <- sds_means %>% 
    #   ungroup %>% 
    #   mutate(p.value = sds_mean_diffs$p.value) %>% 
    #   select(-conf.low, -conf.high) %>% 
    #   mutate(
    #     stratum = as.integer(stratum),
    #     p.value = format_dbl(p.value, 3),
    #     p.value = if_else(stratum == 1, "-", p.value),
    #     std.error = format_dbl(std.error, 2),
    #     estimate = format_dbl(estimate, 2),
    #     delay_var = format_chr(delay_var)
    #   ) %>% 
    #   mutate(mean_se = str_glue("{estimate} ({std.error})")) %>% 
    #   mutate(delay_var = if_else(delay_var == lag(delay_var) & !is.na(lag(delay_var)), "", delay_var)) %>% 
    #   select("Delay Type" = delay_var, 
    #          "SES stratum" = stratum, 
    #          "Mean (SE)" = mean_se,
    #          "p val diff (to stratum 1)" = p.value,
    #          "N" = N)
    # 
    # 
    # 
    # 
    # print.xtable(
    #   xtable(out_table, caption = "Testing Delay Means"),
    #   # type = "html",
    #   caption.placement = "top",
    #   table.placement = "!htbp",
    #   file = "test_delay_means.tex",
    #   include.rownames = FALSE
    # )
    
    
    
    # 5. Share of detected-------------------------------------------------------
    
    
    
    # 2 options:
    # (1) What share of people who are infected are detected at at least one time during their infection?
    # (2) What's the average proportion of infected people at a given t who are detected [weighted by # infected at that t?]
    
    # Number 2 I can already calculate:
    
    # share_detected_instantaneous <- group_diff_slim %>% 
    #   filter(parameter_set != "out_plus_home") %>% 
    #   filter(parameter_set %in% c("baseline", "testing_tracing", "all")) %>%
    #   
    #   mutate(share_detected_inst = n_detected / n_cases_live) %>% 
    #   group_by(parameter_set, sim_id, i_group) %>% 
    #   select(t, share_detected_inst, n_cases_live) %>% 
    #   summarise(share_detected_w = weighted.mean(x = share_detected_inst, w = n_cases_live)) %>% 
    #   group_by(parameter_set, i_group) %>% 
    #   summarise(
    #     q = c(0.025, 0.25, 0.5, 0.75, 0.975),
    #     across(c(share_detected_w),
    #            ~ quantile(.x, q = c(0.025, 0.25, 0.5, 0.75, 0.975)))
    #     # quibble
    #     # quibble(diff_prop_infected, q = c(0.1, 0.25, 0.5, 0.75, 0.9)),
    #     # quibble(diff_tot_infected, q = c(0.1, 0.25, 0.5, 0.75, 0.9))
    #   ) %>% 
    #   mutate(
    #     parameter_set = factor(parameter_set, levels = unique(group_diff_slim$parameter_set))
    #   )
    # 
    # share_detected_instantaneous %>% plot_group_diffs(y = share_detected_w, 
    #                                                   # ylab = "Average share of infected who have tested positive"
    #                                                   ylab = "Avg. share detected"
    # )
    # 
    # 
    # # ggsave("share_detected.png",width = 8, height = 4.5, scale = 0.8)
    # ggsave("share_detected_testing_tracing.png", width = 8, height = 2.5, scale = 0.8)
    
    
    
    
    
    
    # V2 - ratio of confirmed cases / total cases
    
    share_confirmed <- group_diff_slim %>% 
      filter(parameter_set != "out_plus_home") %>% 
      filter(parameter_set %in% c("baseline", "testing_tracing", "all")) %>%
      
      group_by(parameter_set, sim_id, i_group) %>% 
      summarise(n_confirmed_cases = max(n_confirmed_cases),
                n_cases_cum = max(n_cases_cum)) %>% 
      mutate(share_confirmed = n_confirmed_cases / n_cases_cum) %>% 
      
      group_by(parameter_set, i_group) %>% 
      summarise_quantiles(
        share_confirmed, q = c(0.025, 0.25, 0.5, 0.75, 0.975)
      )
    
    
    share_confirmed %>% plot_group_diffs(y = share_confirmed, 
                                         # ylab = "Average share of infected who have tested positive"
                                         ylab = "Avg. share detected")
    
    
    ggsave("figures/share_detected_testing_tracing.png", width = 8, height = 2.5, scale = 0.8)
    
    
# 6. Share isolated ----------------------------------------------------------
    
    
    
    
    share_isolated_instantaneous <- group_diff_slim %>% 
      filter(parameter_set != "out_plus_home") %>% 
      filter(parameter_set %in% c("baseline", "testing_tracing", "all")) %>%
      
      # filter(sim_id == 1, i_group == 1) %>% 
      # filter(n_cases_live > 200) %>% 
      # print(n = 200) %>% 
      # ggplot(aes(x = t, y = n_isolators / n_cases_live, colour = factor(parameter_set))) + 
      # geom_line()
      mutate(isolation_ratio_inst = n_isolators / n_cases_live) %>% 
      group_by(parameter_set, sim_id, i_group) %>% 
      select(t, isolation_ratio_inst, n_cases_live) %>% 
      summarise(isolation_ratio_w = weighted.mean(x = isolation_ratio_inst, w = n_cases_live)) %>% 
      group_by(parameter_set, i_group) %>% 
      summarise(
        q = c(0.025, 0.25, 0.5, 0.75, 0.975),
        across(c(isolation_ratio_w),
               ~ quantile(.x, q = c(0.025, 0.25, 0.5, 0.75, 0.975)))
        # quibble
        # quibble(diff_prop_infected, q = c(0.1, 0.25, 0.5, 0.75, 0.9)),
        # quibble(diff_tot_infected, q = c(0.1, 0.25, 0.5, 0.75, 0.9))
      ) %>% 
      mutate(
        parameter_set = factor(parameter_set, levels = unique(group_diff_slim$parameter_set))
      )
    
    share_isolated_instantaneous %>% plot_group_diffs(y = "isolation_ratio_w", ylab = "Avg. # Isolating / # Infected")
    
    
    # ggsave("share_isolated.png",width = 8, height = 4.5, scale = 0.8)
    ggsave("figures/share_isolated_testing_tracing.png", width = 8, height = 2.5, scale = 0.8)
    
    
    
    
    
# .... --------------------------------------------------------------------
    
    
