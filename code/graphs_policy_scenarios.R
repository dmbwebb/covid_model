load("data/processed/group_diff_tightened_half.RData", verbose = TRUE)
load("data/processed/target_slim.RData", verbose = TRUE)
load("data/processed/target_isolation_slim.RData", verbose = TRUE)
load("data/processed/target_immunity_slim.RData", verbose = TRUE)
# load("data/processed/rapid_testing_slim_v2.RData", verbose = TRUE)
load("data/processed/fast_testing_slim.RData", verbose = TRUE)
load("data/processed/0_testing_slim.RData", verbose = TRUE)

equalish <- function(x, y, tol = 0.001) {
  abs(x - y) < tol
}



# Collect scenarios -------------------------------------------------------




baseline <- group_diff_slim_tightened %>% 
  mutate(scenario = if_else(parameter_set == "baseline", "Baseline", NA_character_)) %>% 
  relocate(scenario) %>% 
  filter(!is.na(scenario))

contacts <- target_slim %>% 
  mutate(scenario = if_else(equalish(mean_reduc_per_person, 1), "Reduce outside-home contacts by 1", NA_character_)) %>% 
  filter(!is.na(scenario))


# target_slim %>% group_by(mean_reduc_per_person) %>% summarise(n = n_distinct(sim_id))


isolation <- target_isolation_slim %>% 
  mutate(scenario = if_else(equalish(mean_increase_per_person, 0.2), "Increase isolation by 20 p.p.", NA_character_)) %>% 
  relocate(scenario) %>% 
  filter(!is.na(scenario))




immunity <- target_immunity_slim %>% 
  mutate(scenario = if_else(equalish(prop_immune, 0.1), "10% initially vaccinated", NA_character_)) %>% 
  relocate(scenario) %>% 
  filter(!is.na(scenario))

testing <- rapid_testing_slim %>%
  select(-p_self_test, -probs_self_test) %>% 
  mutate(
    scenario = case_when(
      equalish(prop_rapid_testing, 1) & mean_increase_self_test == 0 & target_type == "untargeted" ~ "Fast testing",
      equalish(prop_rapid_testing, 0) & equalish(mean_increase_self_test, 0.3)                     ~ "Increase self-testing by 30 p.p.",
      equalish(prop_rapid_testing, 1) & equalish(mean_increase_self_test, 0.3)                     ~ "Fast testing + Increase self-testing by 30 p.p."
    ),
    target_type = if_else(scenario == "Fast testing", NA_character_, target_type) # remove target type for Fast testing scenario (not relevant)
  ) %>% 
  filter(!is.na(scenario))

zero_testing <- model_0_testing %>% 
  mutate(scenario = "No testing")

# zero_testing %>% summarise(n = n_distinct(sim_id))

scenarios <- bind_rows(
  baseline, contacts, isolation, immunity, testing, zero_testing
)

scenarios_summ <- scenarios %>% 
  # Amalgamate all groups
  group_by(scenario, target_type, sim_id, t) %>% 
  summarise(n_pop = sum(n_pop), 
            n_cases_cum = sum(n_cases_cum)) %>% 
  
  # Look at maximum cumulative cases (for end of epidemic)
  group_by(scenario, target_type, sim_id, n_pop) %>% 
  summarise(n_cases_cum = max(n_cases_cum)) %>% 
  
  mutate(prop_infected = n_cases_cum / n_pop) %>% 
  
  # Summarise for quantiles
  group_by(scenario, target_type) %>% 
  quantile_summarise(prop_infected, conf_level = conf_level) %>% 
  ungroup



# P values ----------------------------------------------------------------


# P values (targeted vs non-targeted)
p_target <- scenarios %>% 
  # Amalgamate all groups
  group_by(scenario, target_type, sim_id, t) %>% 
  summarise(n_pop = sum(n_pop), 
            n_cases_cum = sum(n_cases_cum)) %>% 
  
  # Look at maximum cumulative cases (for end of epidemic)
  group_by(scenario, target_type, sim_id, n_pop) %>% 
  summarise(n_cases_cum = max(n_cases_cum)) %>% 
  mutate(prop_infected = n_cases_cum / n_pop) %>% 
  
  # For each scenario, count how many sim_ids are lower on target than on nontargeted
  filter(!is.na(target_type)) %>% 
  group_by(scenario, sim_id) %>% 
  summarise(prop_infected_targeted = prop_infected[target_type == "targeted"],
            prop_infected_untargeted = prop_infected[target_type == "untargeted"]) %>% 
  mutate(targeted_lower = prop_infected_targeted < prop_infected_untargeted) %>% 
  group_by(scenario) %>% 
  summarise(n_targeted_lower = sum(targeted_lower),
            total = n(),
            p_targeted_lower = mean(targeted_lower))

write_excel_csv(p_target, path = "figures/policy_scenarios_n_below_target.csv")


p_target_ttest <- scenarios %>% 
  # Amalgamate all groups
  group_by(scenario, target_type, sim_id, t) %>% 
  summarise(n_pop = sum(n_pop), 
            n_cases_cum = sum(n_cases_cum)) %>% 
  
  # Look at maximum cumulative cases (for end of epidemic)
  group_by(scenario, target_type, sim_id, n_pop) %>% 
  summarise(n_cases_cum = max(n_cases_cum)) %>% 
  mutate(prop_infected = n_cases_cum / n_pop) %>% 
  
  # For each scenario, count how many sim_ids are lower on target than on nontargeted
  filter(!is.na(target_type)) %>% 
  group_by(scenario) %>% 
  summarise(
    t_test = list(stats::t.test(x = prop_infected[target_type == "targeted"],
                                y = prop_infected[target_type == "untargeted"]))
  ) %>% 
  rowwise() %>% 
  mutate(
    p = t_test$p.value
  ) %>% 
  mutate(p_round = round(p, digits = 5)) %>% 
  select(scenario, p, p_round)

write_excel_csv(p_target_ttest, path = "figures/policy_scenarios_p_target.csv")














# P values (compared to baseline)
p_baseline <- scenarios %>% 
  # Amalgamate all groups
  group_by(scenario, target_type, sim_id, t) %>% 
  summarise(n_pop = sum(n_pop), 
            n_cases_cum = sum(n_cases_cum)) %>% 
  
  # Look at maximum cumulative cases (for end of epidemic)
  group_by(scenario, target_type, sim_id, n_pop) %>% 
  summarise(n_cases_cum = max(n_cases_cum)) %>% 
  mutate(prop_infected = n_cases_cum / n_pop) %>% 
  
  ungroup %>% 
  group_by(sim_id) %>% 
  mutate(prop_infected_baseline = prop_infected[scenario == "Baseline"]) %>% 
  mutate(lower_than_baseline = prop_infected < prop_infected_baseline) %>% 
  group_by(scenario, target_type) %>% 
  summarise(
    n_lower_than_baseline = sum(lower_than_baseline),
    p_lower_than_baseline = mean(lower_than_baseline),
    total = n()
  ) %>% 
  count_prop(scenario)

write_excel_csv(p_baseline, path = "figures/policy_scenarios_n_below_baseline.csv")

# P VALUES v2 (t test)
p_scenarios <- scenarios %>% 
  # Amalgamate all groups
  group_by(scenario, target_type, sim_id, t) %>% 
  summarise(n_pop = sum(n_pop), 
            n_cases_cum = sum(n_cases_cum)) %>% 
  
  # Look at maximum cumulative cases (for end of epidemic)
  group_by(scenario, target_type, sim_id, n_pop) %>% 
  summarise(n_cases_cum = max(n_cases_cum)) %>% 
  mutate(prop_infected = n_cases_cum / n_pop) %>% 
  
  ungroup

outcome_baseline <- p_scenarios %>% filter(scenario == "Baseline")

# data_save$probs_self_test
p_inequality_t_test <- p_scenarios %>% 
  filter(scenario != "Baseline") %>% 
  group_by(scenario, target_type) %>% 
  select(scenario, target_type, prop_infected) %>% 
  nest(data = prop_infected) %>% 
  rowwise() %>% 
  mutate(t_test = list(stats::t.test(x = data, y = outcome_baseline$prop_infected)),
         p = t_test$p.value) %>% 
  select(-data, -t_test) %>% 
  mutate(p_round = round(p, digits = 5))


write_excel_csv(p_inequality_t_test, path = "figures/policy_scenarios_p_baseline.csv")










# PLOT --------------------------------------------------------------------


# PLOT SCENARIOS
plot_policy_scenarios <- scenarios_summ %>% 
  mutate(scenario = factor(scenario,
                           levels = c(
                             "Baseline",
                             "10% initially vaccinated",
                             "Reduce outside-home contacts by 1",
                             "Increase isolation by 20 p.p.",
                             "Increase self-testing by 30 p.p.",
                             "No testing",
                             "Fast testing",
                             "Fast testing + Increase self-testing by 30 p.p."
                           )),
         target_type = fct_relevel(factor(target_type), "untargeted"),
         target_type = fct_explicit_na(target_type),
         target_type = factor(target_type, labels = c("Untargeted", "Targeted (SES 1&2)", "NA"))) %>% 
  ggplot(aes(x = scenario, colour = target_type)) + 
  geom_point(aes(y = prop_infected_median), position = position_dodge(width = 0.2)) + 
  geom_errorbar(aes(ymin = prop_infected_lower, ymax = prop_infected_upper), width = 0.1, position = position_dodge(width = 0.2)) + 
  theme_custom() +
  scale_y_continuous(labels = scales::percent_format(accuracy = 1)) + 
  scale_color_manual(breaks = c("Untargeted", "Targeted (SES 1&2)"),
                     values = c("#619CFF", "#F8766D", "#2e4057")) + 
  theme(axis.text.x = element_text(angle = 45, hjust = 1),
        panel.grid.minor = element_blank(),
        panel.grid.major.x = element_blank(),
        legend.position = "bottom") + 
  labs(y = "Cumulative per capita incidence\n(all SES groups)", x = element_blank(), colour = element_blank())
  
plot_policy_scenarios
ggsave("figures/policy_scenarios.pdf", device = cairo_pdf, width = 8, height = 6, scale = 0.8)

plot_policy_scenarios + geom_label_repel(aes(y = prop_infected_median, label = round(prop_infected_median, 3)), show.legend = FALSE)
ggsave("figures/policy_scenarios_vallabels.pdf", device = cairo_pdf, width = 8, height = 6, scale = 0.8)
   







# Fast testing by SES -----------------------------------------------------




testing_by_group <- testing %>% 
  bind_rows(baseline) %>% 
  bind_rows(zero_testing) %>% 
  # Amalgamate all groups
  # group_by(scenario, target_type, sim_id, t) %>% 
  # summarise(n_pop = sum(n_pop), 
  #           n_cases_cum = sum(n_cases_cum)) %>% 
  
  # Look at maximum cumulative cases (for end of epidemic)
  group_by(scenario, target_type, i_group, sim_id, n_pop) %>% 
  summarise(n_cases_cum = max(n_cases_cum)) %>% 
  
  mutate(prop_infected = n_cases_cum / n_pop) %>% 
  
  # Summarise for quantiles
  group_by(scenario, target_type, i_group) %>% 
  quantile_summarise(prop_infected, conf_level = conf_level) %>% 
  mutate(scenario_label = paste0(scenario, " ", target_type)) %>% 
  ungroup




## Reordering testing_by_group$scenario_label
testing_by_group %>% 
  mutate(scenario_label = factor(scenario_label,  levels = c(
    "Baseline NA",
    "Increase self-testing by 30 p.p. untargeted",
    "Increase self-testing by 30 p.p. targeted",
    "No testing NA",
    "Fast testing NA",
    "Fast testing + Increase self-testing by 30 p.p. untargeted",
    "Fast testing + Increase self-testing by 30 p.p. targeted"
  ))) %>% 
  mutate(scenario_label = fct_recode(scenario_label,
                                     "Baseline" = "Baseline NA",
                                     "No testing" = "No testing NA",
                                     "Fast testing +\nIncrease self-testing by 30 p.p.\n(Targeted)" = "Fast testing + Increase self-testing by 30 p.p. targeted",
                                     "Fast testing +\nIncrease self-testing by 30 p.p.\n(Untargeted)" = "Fast testing + Increase self-testing by 30 p.p. untargeted",
                                     "Fast testing" = "Fast testing NA",
                                     "Increase self-testing by 30 p.p.\n(Targeted)" = "Increase self-testing by 30 p.p. targeted",
                                     "Increase self-testing by 30 p.p.\n(Untargeted)" = "Increase self-testing by 30 p.p. untargeted")) %>% 
  ggplot(aes(x = i_group, fill = factor(i_group))) + 
  geom_col(aes(y = prop_infected_mean)) + 
  geom_errorbar(aes(ymin = prop_infected_lower, ymax = prop_infected_upper), width = 0.2, size = 0.5, colour = "#374561") + 
  facet_wrap(~ scenario_label, ncol = 4) +
  labs(y = "Cumulative per capita incidence (entire epidemic)", x = "SES Group", fill = "SES Group") + 
  theme_custom() + 
  theme(legend.position = c(0.9, 0.2), panel.grid.minor = element_blank(), panel.grid.major.x = element_blank()) + 
  geom_hline(yintercept = 0, colour = "darkgrey", size = 0.4) + 
  scale_y_continuous(labels = scales::percent_format(accuracy = 1)) +
  scale_fill_viridis(option = "viridis",
                     begin = 0.3,
                     end = 1,
                     discrete = TRUE)


ggsave("figures/fast_testing_by_group.pdf", device = cairo_pdf, width = 11, height = 6, scale = 0.8)
  

