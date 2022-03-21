# COUNTERFACTUALs ---------------------------------------------------------



# load("data/processed/group_diff_slim.RData", verbose = TRUE)
load("data/processed/group_diff_tightened_half.RData", verbose = TRUE)
load("data/processed/group_diff_split.RData", verbose = TRUE)


# group_diff_slim_split

# Labels for each counterfacutal type
counterfactual_labels <- c(
  "baseline" = "Baseline",
  "out_contacts" = "Out of home",
  "home_contacts" = "Within home",
  "isolation_behaviour" = "Isolation behavior",
  "testing_tracing" = "Testing & Tracing",
  "all" = "All",
  "hh_size" = "HH size",
  "sar_home" = "SAR Home"
)


# Prepare plot DF
group_diff_at_t <- bind_rows(group_diff_slim_tightened, group_diff_slim_split) %>% 
  mutate(
    parameter_set_label = str_replace_all(parameter_set, !!counterfactual_labels),
    parameter_set_label = as.character(case_when(
      tighten_factor == 0.5 ~ str_glue("{parameter_set_label} (50%)"),
      TRUE                  ~ parameter_set_label
    )),
    parameter_set_label = fct_relevel(factor(parameter_set_label), "Baseline")
  ) %>% 
  count_prop(parameter_set_label) %>% 
  # mutate(parameter_set = paste0(parameter_set, "_", tighten_factor)) %>% 
  # mutate(parameter_set = paste0(parameter_set, "_", tighten_factor)) %>% 
  # group_by(parameter_set_label, sim_id, i_group) %>% 
  arrange(parameter_set_label, sim_id, i_group, t) %>% 
  complete(parameter_set_label, sim_id, i_group, t) %>% 
  group_by(parameter_set_label, sim_id, i_group) %>% 
  fill(n_pop, n_cases_cum, .direction = "down") %>% 
  # select(parameter_set_label, sim_id, i_group, t, n_cases_cum, n_pop) %>% 
  group_by(parameter_set_label, sim_id, i_group) %>% 
  mutate(new_cases_2_week = n_cases_cum - lag(n_cases_cum, 14)) %>% 
  
  # group_diff_at_t %>% sample_n_groups(2) %>% head() %>%  print
  group_by(parameter_set_label, t, i_group, n_pop) %>% 
  quantile_summarise(c(n_cases_live, new_cases_2_week), conf_level = 0) %>% 
  # print
  mutate(incidence_prop = if_else(is.nan(new_cases_2_week_median), NA_real_, new_cases_2_week_median / n_pop)) %>% 
  ungroup %>% 
  mutate(
    parameter_set_label = factor(parameter_set_label,
                                 levels = c(
                                   "Baseline", 
                                   "Out of home", "Out of home (50%)", "Within home",
                                   "Within home (50%)", "Isolation behavior", "Isolation behavior (50%)",
                                   "Testing & Tracing", "Testing & Tracing (50%)", "HH size", "HH size (50%)",
                                   "SAR Home", "SAR Home (50%)", "All", "All (50%)"
                                 )
    )
  ) %>% 
  filter(!str_detect(parameter_set_label, "All"))



# Check how many people are false negative
bind_rows(group_diff_slim_tightened, group_diff_slim_split) %>% 
  mutate(
    parameter_set_label = str_replace_all(parameter_set, !!counterfactual_labels),
    parameter_set_label = as.character(case_when(
      tighten_factor == 0.5 ~ str_glue("{parameter_set_label} (50%)"),
      TRUE                  ~ parameter_set_label
    )),
    parameter_set_label = fct_relevel(factor(parameter_set_label), "Baseline")
  ) %>% 
  count_prop(parameter_set_label) %>% 
  select(parameter_set_label, i_group, sim_id, t, eventually_confirmed_self, eventually_false_neg_self) %>% 
  group_by(parameter_set_label, sim_id) %>% 
  filter(t == max(t)) %>% 
  group_by(parameter_set_label, sim_id) %>% 
  summarise(across(c(eventually_confirmed_self, eventually_false_neg_self), sum)) %>% 
  mutate(prop_confirmed = eventually_confirmed_self / (eventually_false_neg_self + eventually_confirmed_self)) %>% 
  
  summarise(mean_se(prop_confirmed)) %>% 
  
  ggplot(aes(x = parameter_set_label, y = y)) + 
  geom_point() + 
  geom_errorbar(aes(ymin = ymin, ymax = ymax)) + 
  coord_flip()

# group_diff_slim_tightened %>% print_names


# check_false_negs <- sim_for_r0 %>%
#   group_by(sim_id, t) %>%
#   # sum up across i_groups
#   summarise(across(-c(i_group, policy, prop_susceptible), sum)) %>%
#   mutate(prop_confirmed = eventually_confirmed_self / (eventually_false_neg_self + eventually_confirmed_self)) %>%
#   select(-c(n_isolators:n_confirmed_cases), n_cases_cum, n_confirmed_cases) %>%
#   mutate(prop_self_test_yn = all_self_test_yn / n_cases_cum) %>%
#   tail()



# 1. Counterfactual epidemic curve -----------------------------------------------------------



# PLOT epidemic curves
# included_curves <- c("baseline", "out_contacts", "home_contacts", "isolation_behaviour", "testing_tracing", "all")
# group_diff_at_t %>% ungroup %>% count_prop(parameter_set)

plot_epidemic_curve <- group_diff_at_t %>% 
  filter(!str_detect(parameter_set_label, "SAR|HH")) %>% 
  mutate(fifty = fct_recode(factor(str_detect(parameter_set_label, "50%")), "50%" = "TRUE", "100%" = "FALSE"),
         parameter_set = str_replace_all(parameter_set_label, " \\(50%\\)", "")) %>% 
  mutate(parameter_set = factor(parameter_set, levels = c("Baseline", "Out of home", "Within home", "Isolation behavior", "Testing & Tracing"))) %>% 
  ungroup %>% 
  ungroup %>% 
  mutate(i_group = recode_i_group(i_group)) %>% 
  ggplot(aes(x = t, y = incidence_prop, 
             colour = factor(i_group), group = factor(i_group)
             # alpha = 1
  )) + 
  geom_line(show.legend = TRUE) + 
  facet_grid(fifty ~ parameter_set) + 
  # facet_wrap(~ parameter_set, labeller = as_labeller(counterfactual_labels)) + 
  
  # Formatting
  theme_custom() + 
  theme(legend.position = "bottom") + 
  # scale_alpha_continuous(guide = FALSE) + 
  scale_y_continuous(labels = scales::percent_format(accuracy = 1), minor_breaks = NULL) + 
  labs(colour = "SES Group", fill = "SES Group", y = "Per capita incidence\n(previous 2 weeks)") +
  geom_hline(yintercept = 0, colour = "darkgrey", size = 0.4) + 
  # coord_cartesian(xlim = c(0, 500)) +
  scale_colour_viridis(option = "viridis",
                       begin = 0.3, 
                       discrete = TRUE) + 
  coord_cartesian(xlim = c(0, 400))
# scale_x_continuous(minor_breaks = c(0, 25, 50, 75), breaks = c(0, 50, 100))

plot_epidemic_curve <- graph_remove_bottom_left(plot_epidemic_curve)


ggsave("figures/counterfactuals_epidemic_curve_appendix.pdf", width = 10, height = 5.5, scale = 0.8, plot = plot_epidemic_curve)


plot_home_epidemic_curve <- group_diff_at_t %>% 
  filter(str_detect(parameter_set_label, "Baseline|Within home|SAR|HH")) %>% 
  mutate(fifty = fct_recode(factor(str_detect(parameter_set_label, "50%")), "50%" = "TRUE", "100%" = "FALSE"),
         parameter_set = str_replace_all(parameter_set_label, " \\(50%\\)", "")) %>% 
  mutate(parameter_set = factor(parameter_set,levels = c(
    "Baseline", "Within home", "HH size", "SAR Home"
  ))) %>% 
  ungroup %>% 
  mutate(i_group = recode_i_group(i_group)) %>% 
  ggplot(aes(x = t, y = incidence_prop, 
             colour = factor(i_group), group = factor(i_group)
             # alpha = 1
  )) + 
  geom_line(show.legend = TRUE) + 
  facet_grid(fifty ~ parameter_set) + 
  # facet_wrap(~ parameter_set, labeller = as_labeller(counterfactual_labels)) + 
  
  # Formatting
  theme_custom() + 
  theme(legend.position = "top") + 
  # scale_alpha_continuous(guide = FALSE) + 
  scale_y_continuous(labels = scales::percent_format(accuracy = 1), minor_breaks = NULL) + 
  labs(colour = "SES Group", fill = "SES Group", y = "Per capita incidence (previous 2 weeks)") +
  geom_hline(yintercept = 0, colour = "darkgrey", size = 0.4) + 
  # coord_cartesian(xlim = c(0, 500)) +
  scale_colour_viridis(option = "viridis",
                       begin = 0.3, 
                       discrete = TRUE) + 
  coord_cartesian(xlim = c(0, 400))
# scale_x_continuous(minor_breaks = c(0, 25, 50, 75), breaks = c(0, 50, 100))



plot_home_epidemic_curve <- graph_remove_bottom_left(plot_home_epidemic_curve)

# plot_home_epidemic_curve
# library(grid)
# library(gridExtra)
# g <- ggplotGrob(plot_home_epidemic_curve)
# # get the grobs that must be removed
# # rm_grobs <- g$layout$name %in% c("panel-2-1", "panel-3-3", "strip-t-3-1", "strip-t-3-3")
# rm_grobs <- g$layout$name %in% c("panel-2-1", "axis-b-1")
# # remove grobs
# g$grobs[rm_grobs] <- NULL
# g$layout <- g$layout[!rm_grobs, ]
# ## move axis closer to panel
# 
# g <- gtable::gtable_add_cols(g, unit(0.3, "null"), pos = 5)
# 
# # g$layout[g$layout$name == "panel-1-1", c("t", "b")] <- c(, 8)
# grid.newpage()
# grid_plot <- grid.arrange(g)

ggsave("figures/counterfactuals_home_epidemic_curve.pdf", width = 8, height = 4.5, scale = 0.9, plot = plot_home_epidemic_curve)


# ggsave("figures/prop_infected_main.pdf", width = 9.5, height = 4, scale = 0.8, plot = grid_plot)




# 2. Counterfactual - summary results -------------------------------------

# group_diff_slim %>% count_prop(parameter_set, tighten_factor)

# SUMMARISE ALL INFO
group_diff_all_summ <- bind_rows(group_diff_slim_tightened, group_diff_slim_split) %>% 
  mutate(
    parameter_set_label = str_replace_all(parameter_set, !!counterfactual_labels),
    parameter_set_label = as.character(case_when(
      tighten_factor == 0.5 ~ str_glue("{parameter_set_label} (50%)"),
      TRUE                  ~ parameter_set_label
    )),
    parameter_set_label = fct_relevel(factor(parameter_set_label), "Baseline")
  ) %>% 
  count_prop(parameter_set_label) %>%
  select(-parameter_set) %>% 
  
  # (1) Calculate total cases
  group_by(parameter_set_label, i_group, sim_id, n_pop) %>% 
  summarise(n_cases = max(n_cases_cum)) %>% 
  
  # Calculate proportion infected, and gaps with group 5
  group_by(parameter_set_label, sim_id) %>% arrange(sim_id) %>% 
  mutate(prop_infected = n_cases / n_pop) %>% 
  mutate(diff_prop_infected = prop_infected - prop_infected[i_group == 4],
         diff_n_cases = n_cases - n_cases[i_group == 4]) %>% 
  
  
  # Calculate "effect" relative to baseline
  group_by(sim_id, i_group) %>% 
  mutate(across(-c(parameter_set_label), list(effect = ~ .x - .x[parameter_set_label == "Baseline"]))) %>% 
  # Calculte implied contribution for each group, for each sim
  # mutate(contribution = diff_prop_infected_effect / diff_prop_infected_effect[parameter_set_label == "All"]) %>% 
  group_by(parameter_set_label, i_group) %>% 
  arrange(parameter_set_label, i_group) %>% 
  mutate(
    parameter_set_label = factor(parameter_set_label,
                                 levels = c(
                                   "Out of home", "Out of home (50%)", "Within home",
                                   "Within home (50%)", "Isolation behavior", "Isolation behavior (50%)",
                                   "Testing & Tracing", "Testing & Tracing (50%)", "HH size", "HH size (50%)",
                                   "SAR Home", "SAR Home (50%)", "All", "All (50%)", "Baseline"
                                 )
    )
  ) %>% 
  filter(!str_detect(parameter_set_label, "All"))
# summarise_quantiles(-c(sim_id)) %>% print_all



# Make into function
plot_group_diffs <- function(.data, y, ylab = "Y", text_label = FALSE, text_label_thresh = 0, label_nudge = 0.1) {
  
  # .data <- group_diff_all_summ; y = 
  # plot_group_diffs(group_diff_all_summ, y = prop_infected, "Cumulative per capita incidence (entire epidemic)", text_label = TRUE, text_label_thresh = 10)
  
  
  data_processed <- .data %>% 
    mutate(i_group = recode_i_group(i_group)) %>% 
    summarise_quantiles({{y}}) %>% 
    select(parameter_set_label, i_group, q, {{y}}) %>% 
    filter(q %in% c(0.025, 0.5, 0.975)) %>% 
    pivot_wider(names_from = q, names_prefix = "y_", values_from = {{y}})
  
  # data_processed <- group_diff_all_summ %>% 
  #   mutate(i_group = recode_i_group(i_group)) %>% 
  #   summarise_quantiles(prop_infected) %>% 
  #   select(parameter_set_label, i_group, q, prop_infected) %>% 
  #   filter(q %in% c(0.025, 0.5, 0.975)) %>% 
  #   pivot_wider(names_from = q, names_prefix = "y_", values_from = prop_infected)
  
  p <- ggplot(data_processed, aes(x = i_group, fill = factor(i_group))) + 
    geom_col(aes(y = y_0.5)) + 
    geom_errorbar(aes(ymin = y_0.025, ymax = y_0.975), width = 0.2, size = 0.5, colour = "#374561") + 
    facet_wrap(~ parameter_set_label) +
    labs(y = ylab, x = "SES Group", fill = "SES Group") + 
    theme_custom() + 
    theme(legend.position = "bottom", panel.grid.minor = element_blank(), panel.grid.major.x = element_blank()) + 
    geom_hline(yintercept = 0, colour = "darkgrey", size = 0.4) + 
    scale_y_continuous(labels = scales::percent_format(accuracy = 1)) +
    scale_fill_viridis(option = "viridis",
                       begin = 0.3,
                       end = 1,
                       discrete = TRUE)
  # scale_x_continuous(breaks = 1:5)
  
  if (text_label) {
    p <- p + 
      geom_label(data = data_processed %>% filter(y_0.5 * 100 >= text_label_thresh),
                aes(y = y_0.5, label = paste0(round(y_0.5, 2) * 100)),
                size = 2,
                label.size = 0.2,
                label.padding = unit(0.15, "lines"),
                colour = "black", fill = "white", alpha = 0.5, position = position_stack(vjust = 0.5)) +
      geom_label(data = data_processed %>% filter(y_0.5 * 100 < text_label_thresh & round(y_0.5, 2) * 100 >= 1),
                aes(y = y_0.5, label = paste0(round(y_0.5, 2) * 100)),
                size = 2,
                label.size = 0.2,
                label.padding = unit(0.15, "lines"),
                colour = "black",
                nudge_y = label_nudge,
                fill = "white",
                alpha = 0.5
                # position = position_stack(vjust = 0.5)
                )
  }
  
  return(p)
  
}


# PLOT
# Prop infected in each scenario
plot_group_diffs(group_diff_all_summ, y = prop_infected, "Cumulative per capita incidence (entire epidemic)")
ggsave("figures/prop_infected_appendix.pdf", width = 8, height = 4.5, scale = 1)

# plot_group_diffs(group_diff_all_summ, y = prop_infected, "Cumulative per capita incidence (entire epidemic)") + geom_label_repel(aes(y = y_0.5, label = round(y_0.5, 3)))
plot_group_diffs(group_diff_all_summ, y = prop_infected, "Cumulative per capita incidence (entire epidemic)", text_label = TRUE, text_label_thresh = 10)
  # geom_text(aes(y = y_0.5, label = paste0(round(y_0.5, 2) * 100)), colour = "white", position = position_stack(vjust = 0.5), show.legend = FALSE)
ggsave("figures/prop_infected_appendix_vallabels.pdf", width = 8, height = 4.5, scale = 1.3)


 
    # PLOT NEW VERSION FOR MAIN
    group_diff_main <- group_diff_all_summ %>% 
      filter(!str_detect(parameter_set_label, "SAR Home|HH size")) %>% 
      mutate(fifty = fct_recode(factor(str_detect(parameter_set_label, "50%")), "50%" = "TRUE", "100%" = "FALSE"),
             parameter_set = str_replace_all(parameter_set_label, " \\(50%\\)", "")) %>% 
      # count_prop(parameter_set, fifty) %>% 
      mutate(parameter_set = factor(parameter_set,levels = c(
        "Baseline", "Out of home", "Within home", "Isolation behavior",
        "Testing & Tracing", "All"
      ))) %>% 
      group_by(parameter_set, fifty, .add = TRUE)
    
    # NO VAL LABELS
    plot_group_diff_main <- plot_group_diffs(group_diff_main, y = prop_infected, "Cumulative per capita incidence\n(entire epidemic)", text_label = FALSE) + 
      facet_grid(fifty ~ parameter_set, scales = "free_x") + 
      theme(legend.position = c(0.08, 0.24))
    
    plot_group_diff_main
    library(grid)
    library(gridExtra)
    g <- ggplotGrob(plot_group_diff_main)
    # get the grobs that must be removed
    # rm_grobs <- g$layout$name %in% c("panel-2-1", "panel-3-3", "strip-t-3-1", "strip-t-3-3")
    rm_grobs <- g$layout$name %in% c("panel-2-1", "axis-b-1")
    # remove grobs
    g$grobs[rm_grobs] <- NULL
    g$layout <- g$layout[!rm_grobs, ]
    ## move axis closer to panel
    
    g <- gtable::gtable_add_cols(g, unit(0.3, "null"), pos = 5)
    
    # g$layout[g$layout$name == "panel-1-1", c("t", "b")] <- c(, 8)
    grid.newpage()
    grid_plot <- grid.arrange(g)
    
    ggsave("figures/prop_infected_main.pdf", width = 9.5, height = 4, scale = 0.8, plot = grid_plot)
    
    
    # WITH VAL LABELS
    plot_group_diff_main <- plot_group_diffs(
      group_diff_main,
      y = prop_infected,
      "Cumulative per capita incidence\n(entire epidemic)",
      text_label = TRUE,
      text_label_thresh = 7,
      label_nudge = 0.06
    ) +
      facet_grid(fifty ~ parameter_set, scales = "free_x") +
      theme(legend.position = c(0.08, 0.24))
    
    plot_group_diff_main
    library(grid)
    library(gridExtra)
    g <- ggplotGrob(plot_group_diff_main)
    # get the grobs that must be removed
    # rm_grobs <- g$layout$name %in% c("panel-2-1", "panel-3-3", "strip-t-3-1", "strip-t-3-3")
    rm_grobs <- g$layout$name %in% c("panel-2-1", "axis-b-1")
    # remove grobs
    g$grobs[rm_grobs] <- NULL
    g$layout <- g$layout[!rm_grobs, ]
    ## move axis closer to panel
    
    g <- gtable::gtable_add_cols(g, unit(0.3, "null"), pos = 5)
    
    # g$layout[g$layout$name == "panel-1-1", c("t", "b")] <- c(, 8)
    grid.newpage()
    grid_plot <- grid.arrange(g)
    
    ggsave("figures/prop_infected_main_vallabels.pdf", width = 9.5, height = 4, scale = 0.8, plot = grid_plot)






# PLOT SAR HOME / HH size
group_diff_split_params <- group_diff_all_summ %>% 
  filter(str_detect(parameter_set_label, "Baseline|Within home|SAR Home|HH size")) %>% 
  mutate(fifty = fct_recode(factor(str_detect(parameter_set_label, "50%")), "50%" = "TRUE", "100%" = "FALSE"),
         parameter_set = str_replace_all(parameter_set_label, " \\(50%\\)", "")) %>% 
  # count_prop(parameter_set, fifty) %>% 
  mutate(parameter_set = factor(parameter_set,levels = c(
    "Baseline", "Within home", "HH size", "SAR Home"
  ))) %>% 
  group_by(parameter_set, fifty, .add = TRUE)

plot_group_diff_split_params <- plot_group_diffs(group_diff_split_params, 
                                                 y = prop_infected, 
                                                 "Cumulative per capita incidence (entire epidemic)",
                                                 text_label = TRUE,
                                                 text_label_thresh = 7,
                                                 label_nudge = 0.06
                                                 ) + 
  facet_grid(fifty ~ parameter_set)

plot_group_diff_split_params <- graph_remove_bottom_left(plot_group_diff_split_params)

ggsave("figures/prop_infected_within_home_split.pdf", width = 8, height = 4.5, scale = 0.9, plot = plot_group_diff_split_params)


 # Gap (relative to group 6) in each scenario
plot_group_diffs(group_diff_all_summ, y = diff_prop_infected, "Gap relative to group 6")
ggsave("figures/inequality_appendix.pdf", width = 8, height = 4.5, scale = 1)

plot_group_diffs(group_diff_all_summ, y = diff_prop_infected, "Gap relative to group 6") + geom_label_repel(aes(y = y_0.5, label = round(y_0.5, 3)))
ggsave("figures/inequality_appendix_vallabels.pdf", width = 8, height = 4.5, scale = 1.3)


# 3. Decoposition (v2) ----------------------------------------------------

# DOESN'T WORK ANY MORE BECAUSE WE DON'T HAVE "ALL"


# V2 - no pie chart, just bars
# Relative to total inequality, how much does each parameter change reduce


# contributions_df <- group_diff_all_summ %>% 
#   # mutate(parameter_set = paste0(parameter_set, "_", tighten_factor)) %>% 
#   # group_by(sim_id, parameter_set,  i_group) %>% 
#   arrange(sim_id, parameter_set_label,  i_group) %>% 
#   select(sim_id, parameter_set_label,  i_group, diff_n_cases) %>% 
#   
#   # For each sim_id, i_group, what's the prop reduction
#   group_by(sim_id, i_group) %>% 
#   filter(i_group != 4) %>%  # remove 4 because we're discussing related to top group
#   mutate(prop_reduction = (diff_n_cases[parameter_set_label == "Baseline"] - diff_n_cases) / diff_n_cases[parameter_set_label == "Baseline"]) %>% 
#   
#   # Calculate quantiles across sim_ids
#   group_by(i_group, parameter_set_label) %>% 
#   summarise_quantiles(prop_reduction) %>% 
#   filter(q %in% c(0.05, 0.5, 0.95)) %>% 
#   mutate(q_label = case_when(q == 0.05 ~ "_lower",
#                              q == 0.5   ~ "",
#                              q == 0.95 ~ "_upper")) %>% 
#   select(-q) %>%
#   ungroup %>% 
#   pivot_wider(names_from = q_label,
#               values_from = prop_reduction,
#               names_prefix = "prop_reduction")
# 
# # Plot
# plot_inequality_decomposition <- contributions_df %>% 
#   mutate(i_group = recode_i_group(i_group)) %>% 
#   ggplot(aes(x = parameter_set_label, y = prop_reduction, fill = i_group)) + 
#   geom_col(position = "dodge", width = 0.7) + 
#   geom_errorbar(aes(ymin = prop_reduction_lower, ymax = prop_reduction_upper), position = position_dodge(0.7), width = 0.2, size = 0.5, colour = "#374561") + 
#   labs(x = "Scenario", y = "% reduction in inequality\nrelative to group 5&6") + 
#   scale_y_continuous(labels = scales::percent_format(accuracy = 1)) + 
#   scale_fill_viridis(option = "viridis",
#                      begin = 0.3, 
#                      end = 0.766666667,
#                      discrete = TRUE,
#                      name = "SES Group") + 
#   # scale_x_discrete(labels = counterfactual_labels) + 
#   theme_custom() +
#   theme(axis.text.x = element_text(angle = 45, hjust = 1),
#         panel.grid.minor = element_blank(),
#         panel.grid.major.x = element_blank()) + 
#   coord_cartesian(ylim = c(-0.1, 1)) + 
#   geom_hline(yintercept = 0)
# 
# plot_inequality_decomposition
# 
# ggsave("figures/inequality_decomposition_appendix.pdf",  width = 6, height = 4.5, scale = 1)




# plot_inequality_decomposition + geom_label_repel(aes(y = prop_reduction, label = round(prop_reduction, 3)), position = position_dodge(0.7))
# ggsave("figures/inequality_decomposition_appendix_vallabels.pdf",  width = 6, height = 4.5, scale = 3)





# Plot data ---------------------------------------------------------------



plot_group_diff_main <- plot_group_diffs(
  group_diff_main,
  y = prop_infected,
  "Cumulative per capita incidence\n(entire epidemic)",
  text_label = TRUE,
  text_label_thresh = 7,
  label_nudge = 0.06
) +
  facet_grid(fifty ~ parameter_set, scales = "free_x") +
  theme(legend.position = c(0.08, 0.24))


group_diff_main %>% 
  mutate(i_group = recode_i_group(i_group)) %>% 
  summarise_quantiles(prop_infected) %>% 
  select(parameter_set_label, i_group, q, prop_infected) %>% 
  filter(q %in% c(0.025, 0.5, 0.975)) %>% 
  pivot_wider(names_from = q, names_prefix = "y_", values_from = prop_infected) %>% 
  select(parameter_set, prop_reduction_in_inequality = fifty, parameter_set_label, ses_group = i_group, y_lower = y_0.025, y = y_0.5, y_upper = y_0.975) %>% 
  write_excel_csv("data/processed/plot_data_fig3.csv")





# Make into function
plot_group_diffs <- function(.data, y, ylab = "Y", text_label = FALSE, text_label_thresh = 0, label_nudge = 0.1) {
  
  # .data <- group_diff_all_summ; y = 
  # plot_group_diffs(group_diff_all_summ, y = prop_infected, "Cumulative per capita incidence (entire epidemic)", text_label = TRUE, text_label_thresh = 10)
  
  
  data_processed <- .data %>% 
    mutate(i_group = recode_i_group(i_group)) %>% 
    summarise_quantiles(prop_infected) %>% 
    select(parameter_set_label, i_group, q, prop_infected) %>% 
    filter(q %in% c(0.025, 0.5, 0.975)) %>% 
    pivot_wider(names_from = q, names_prefix = "y_", values_from = prop_infected)
  
  # data_processed <- group_diff_all_summ %>% 
  #   mutate(i_group = recode_i_group(i_group)) %>% 
  #   summarise_quantiles(prop_infected) %>% 
  #   select(parameter_set_label, i_group, q, prop_infected) %>% 
  #   filter(q %in% c(0.025, 0.5, 0.975)) %>% 
  #   pivot_wider(names_from = q, names_prefix = "y_", values_from = prop_infected)
  
  p <- ggplot(data_processed, aes(x = i_group, fill = factor(i_group))) + 
    geom_col(aes(y = y_0.5)) + 
    geom_errorbar(aes(ymin = y_0.025, ymax = y_0.975), width = 0.2, size = 0.5, colour = "#374561") + 
    facet_wrap(~ parameter_set_label) +
    labs(y = ylab, x = "SES Group", fill = "SES Group") + 
    theme_custom() + 
    theme(legend.position = "bottom", panel.grid.minor = element_blank(), panel.grid.major.x = element_blank()) + 
    geom_hline(yintercept = 0, colour = "darkgrey", size = 0.4) + 
    scale_y_continuous(labels = scales::percent_format(accuracy = 1)) +
    scale_fill_viridis(option = "viridis",
                       begin = 0.3,
                       end = 1,
                       discrete = TRUE)
  # scale_x_continuous(breaks = 1:5)
  
  if (text_label) {
    p <- p + 
      geom_label(data = data_processed %>% filter(y_0.5 * 100 >= text_label_thresh),
                 aes(y = y_0.5, label = paste0(round(y_0.5, 2) * 100)),
                 size = 2,
                 label.size = 0.2,
                 label.padding = unit(0.15, "lines"),
                 colour = "black", fill = "white", alpha = 0.5, position = position_stack(vjust = 0.5)) +
      geom_label(data = data_processed %>% filter(y_0.5 * 100 < text_label_thresh & round(y_0.5, 2) * 100 >= 1),
                 aes(y = y_0.5, label = paste0(round(y_0.5, 2) * 100)),
                 size = 2,
                 label.size = 0.2,
                 label.padding = unit(0.15, "lines"),
                 colour = "black",
                 nudge_y = label_nudge,
                 fill = "white",
                 alpha = 0.5
                 # position = position_stack(vjust = 0.5)
      )
  }
  
  return(p)
  
}