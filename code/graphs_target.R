save_pdf <- function(file_path, width = 6) {
  ggsave(str_glue("figures/{file_path}.pdf"), device = cairo_pdf, width = width, height = 4.5, scale = 0.8)
}

load("data/processed/target_slim.RData")
load("data/processed/target_isolation_slim.RData")
load("data/processed/target_immunity_slim.RData")
# load("data/processed/rapid_testing_slim_v2.RData")
load("data/processed/fast_testing_slim.RData")

target_slim %>% group_by(mean_reduc_per_person, target_type) %>% summarise(n = n_distinct(sim_id))
target_isolation_slim %>% group_by(mean_increase_per_person, target_type) %>% summarise(n = n_distinct(sim_id))
target_isolation_slim %>% group_by(mean_increase_per_person, target_type) %>% summarise(n = n_distinct(sim_id))
rapid_testing_slim %>% group_by(prop_rapid_testing, mean_increase_self_test, target_type) %>% summarise(n = n_distinct(sim_id)) %>% print_all

# Contacts ----------------------------------------------------------------

n_pop_total <- 100000


target_plot <- function(data, x, x_lab, show.legend = TRUE) {
  data %>% 
    mutate(target_type = fct_recode(factor(target_type, levels = c("untargeted", "targeted")), "Untargeted" = "untargeted", "Targeted" = "targeted")) %>% 
    ggplot(aes(x = {{x}}, y = prop_infected_mean, colour = target_type)) + 
    # geom_errorbar(aes(ymin = prop_infected_lower, ymax = prop_infected_upper), width = 0.1) + 
    geom_ribbon(aes(ymin = prop_infected_lower, ymax = prop_infected_upper, fill = target_type), linetype = "dashed", size = 0.5, alpha = 0.05, show.legend = show.legend) +
    geom_line(size = 1.5, show.legend = show.legend) + 
    geom_point(size = 2, show.legend = show.legend) + 
    theme_custom(legend.position = "top") + 
    scale_colour_manual(values = c("#00A5FF", "#F8766D")) + 
    scale_fill_manual(values = c("#00A5FF", "#F8766D")) + 
    scale_y_continuous(labels = scales::percent_format(accuracy = 1)) + 
    labs(y = "Cumulative Per Capita Incidence\n(all SES groups, entire epidemic)", x = x_lab, colour = element_blank(), fill = element_blank())
}


target_slim %>% 
  
  group_by(mean_reduc_per_person, target_type, sim_id, t) %>% 
  mutate(n_pop_total = sum(n_pop)) %>% 
  # mutate(budget = budget / n_pop_total) %>% 
  ungroup %>% 
  
  # Calculate total infected in each run
  group_by(mean_reduc_per_person, n_pop_total, target_type, sim_id, i_group) %>% 
  summarise(n_cases_cum = max(n_cases_cum)) %>% 
  
  # Sum up each i_group
  summarise(n_cases_cum = sum(n_cases_cum)) %>% 
  # print_all
  
  # Calculate median and confidence intervals
  group_by(mean_reduc_per_person, target_type) %>% 
  mutate(prop_infected = n_cases_cum / n_pop_total) %>% 
  quantile_summarise(prop_infected, conf_level = conf_level) %>% 
  {bind_rows(., .[1, ] %>% mutate(target_type = "targeted"))} %>%  # duplicate first row
  
  target_plot(x = mean_reduc_per_person, x_lab = "Mean reduction in contacts over whole population", show.legend = FALSE)


save_pdf("target_contacts")






# Isolation ---------------------------------------------------------------


target_isolation_slim %>% 
  
  filter(mean_increase_per_person <= 0.2) %>% 
  group_by(mean_increase_per_person, target_type, sim_id, t) %>% 
  mutate(n_pop_total = sum(n_pop)) %>% 
  # mutate(budget = budget / n_pop_total) %>% 
  ungroup %>% 
  
  # Calculate total infected in each run
  group_by(mean_increase_per_person, n_pop_total, target_type, sim_id, i_group) %>% 
  summarise(n_cases_cum = max(n_cases_cum)) %>% 
  
  # Sum up each i_group
  summarise(n_cases_cum = sum(n_cases_cum)) %>% 
  # print_all
  
  # Calculate median and confidence intervals
  group_by(mean_increase_per_person, target_type) %>% 
  mutate(prop_infected = n_cases_cum / n_pop_total) %>% 
  quantile_summarise(prop_infected, conf_level = conf_level) %>% 
  # { d <- bind_rows(., .[1, ] %>% mutate(target_type = "targeted"))
  # return(d)} %>%  # duplicate first row
  
  # PLOT
  target_plot(x = mean_increase_per_person, x_lab = "Mean increase in isolation probability over whole population", show.legend = FALSE)



save_pdf("target_isolation")




# Immunity ----------------------------------------------------------------






target_immunity_slim %>% 
  
  filter(prop_immune <= 0.3) %>%
  group_by(prop_immune, target_type, sim_id, t) %>% 
  mutate(n_pop_total = sum(n_pop)) %>% 
  # mutate(budget = budget / n_pop_total) %>% 
  ungroup %>% 
  
  # Calculate total infected in each run
  group_by(prop_immune, n_pop_total, target_type, sim_id, i_group) %>% 
  summarise(n_cases_cum = max(n_cases_cum)) %>% 
  
  # Sum up each i_group
  summarise(n_cases_cum = sum(n_cases_cum)) %>% 
  # print_all
  
  # Calculate median and confidence intervals
  group_by(prop_immune, target_type) %>% 
  mutate(prop_infected = n_cases_cum / n_pop_total) %>% 
  quantile_summarise(prop_infected, conf_level = conf_level) %>% 
  # { d <- bind_rows(., .[1, ] %>% mutate(target_type = "targeted"))
  # return(d)} %>%  # duplicate first row
  
  # PLOT
  target_plot(x = prop_immune, x_lab = "Proportion of entire population who are immune", show.legend = FALSE)


save_pdf("target_immunity")





# Testing (OLD) -----------------------------------------------------------------


# load("data/processed/rapid_testing_slim_v2.RData")
# 
# ## Recoding rapid_testing_slim$self_test_type
# 
# rapid_testing_slim %>% 
#   mutate(
#     self_test_type = fct_recode(factor(self_test_type),
#                "Baseline probability of test" = "baseline",
#                "High probability of test (+30p.p.)" = "high"
#     )
#   ) %>% 
#   group_by(total_prop_with_rapid_testing, target_type, self_test_type, sim_id, t) %>% 
#   mutate(n_pop_total = sum(n_pop)) %>% 
#   # mutate(budget_rapid = budget_rapid / n_pop_total) %>% 
#   ungroup %>% 
#   
#   # Calculate total infected in each run
#   group_by(total_prop_with_rapid_testing, target_type, self_test_type, sim_id, i_group) %>% 
#   summarise(n_cases_cum = max(n_cases_cum)) %>% 
#   
#   # Sum up each i_group
#   summarise(n_cases_cum = sum(n_cases_cum)) %>% 
#   # print_all
#   
#   # Calculate median and confidence intervals
#   group_by(total_prop_with_rapid_testing, target_type, self_test_type) %>% 
#   mutate(prop_infected = n_cases_cum / n_pop_total) %>%
#   quantile_summarise(prop_infected, conf_level = conf_level) %>% 
#   arrange(self_test_type) %>% 
#   {bind_rows(., 
#              .[1, ] %>% mutate(target_type = "targeted"),
#              .[nrow(.)/2 + 1, ] %>% mutate(target_type = "targeted"))} %>% 
#   
#   # PLOT
#   target_plot(x = total_prop_with_rapid_testing, x_lab = "Proportion of population\nwith access to rapid testing") + 
#   facet_wrap(~ self_test_type)
#   
# 
# save_pdf("target_rapid_testing")




# Testing (NEW) -----------------------------------------------------------






rapid_testing_slim %>% group_by(prop_rapid_testing) %>% summarise(n = n_distinct(sim_id))

## Recoding rapid_testing_slim$self_test_type

rapid_testing_slim %>% 
  
  group_by(prop_rapid_testing, mean_increase_self_test, target_type, sim_id, t) %>% 
  mutate(n_pop_total = sum(n_pop)) %>% 
  # mutate(budget = budget / n_pop_total) %>% 
  ungroup %>% 
  
  # Calculate total infected in each run
  group_by(prop_rapid_testing, mean_increase_self_test, target_type, n_pop_total, sim_id, i_group) %>% 
  summarise(n_cases_cum = max(n_cases_cum)) %>% 
  
  # Sum up each i_group
  summarise(n_cases_cum = sum(n_cases_cum)) %>% 
  # print_all
  
  # Calculate median and confidence intervals
  group_by(prop_rapid_testing, mean_increase_self_test, target_type) %>% 
  mutate(prop_infected = n_cases_cum / n_pop_total) %>% 
  quantile_summarise(prop_infected, conf_level = conf_level) %>% 
  
  mutate(
    rapid_testing_label = fct_recode(factor(prop_rapid_testing, levels = c(0, 0.5, 1)),
                                     "0% Fast Testing" = "0",
                                     "50% Fast Testing" = "0.5",
                                     "100% Fast Testing" = "1")
  ) %>% 
  
  filter(prop_rapid_testing != 0.5) %>% 
  
  # Duplicate the values with 0 self-testing (for targeted/untargeted)
  {
    bind_rows(., mutate(filter(., mean_increase_self_test == 0 & target_type == "untargeted"), target_type = "targeted"))
  } %>% 
  arrange(target_type) %>% 
  
  # PLOT
  target_plot(x = mean_increase_self_test, x_lab = "Mean increase in self-testing probability\nacross whole population", show.legend = FALSE) + 
  facet_wrap(~ rapid_testing_label)


# rapid_testing_slim %>% 
#   mutate(
#     self_test_type = fct_recode(factor(self_test_type),
#                                 "Baseline probability of test" = "baseline",
#                                 "High probability of test (+30p.p.)" = "high"
#     )
#   ) %>% 
#   group_by(total_prop_with_rapid_testing, target_type, self_test_type, sim_id, t) %>% 
#   mutate(n_pop_total = sum(n_pop)) %>% 
#   # mutate(budget_rapid = budget_rapid / n_pop_total) %>% 
#   ungroup %>% 
#   
#   # Calculate total infected in each run
#   group_by(total_prop_with_rapid_testing, target_type, self_test_type, sim_id, i_group) %>% 
#   summarise(n_cases_cum = max(n_cases_cum)) %>% 
#   
#   # Sum up each i_group
#   summarise(n_cases_cum = sum(n_cases_cum)) %>% 
#   # print_all
#   
#   # Calculate median and confidence intervals
#   group_by(total_prop_with_rapid_testing, target_type, self_test_type) %>% 
#   mutate(prop_infected = n_cases_cum / n_pop_total) %>%
#   quantile_summarise(prop_infected, conf_level = conf_level) %>% 
#   arrange(self_test_type) %>% 
#   {bind_rows(., 
#              .[1, ] %>% mutate(target_type = "targeted"),
#              .[nrow(.)/2 + 1, ] %>% mutate(target_type = "targeted"))} %>% 
#   
#   # PLOT
#   target_plot(x = total_prop_with_rapid_testing, x_lab = "Proportion of population\nwith access to rapid testing") + 
#   facet_wrap(~ self_test_type)


save_pdf("target_fast_testing", width = 8)




# CHECK PROP TESTED POSITIVE


# NEXT : work out how to calculate this for the testing scenarios

rapid_testing_slim %>% 
  # group_by(prop_rapid_testing, mean_increase_self_test, target_type, sim_id, t)
  # mutate(
  #   parameter_set_label = str_replace_all(parameter_set, !!counterfactual_labels),
  #   parameter_set_label = as.character(case_when(
  #     tighten_factor == 0.5 ~ str_glue("{parameter_set_label} (50%)"),
  #     TRUE                  ~ parameter_set_label
  #   )),
  #   parameter_set_label = fct_relevel(factor(parameter_set_label), "Baseline")
  # ) %>% 
  # count_prop(parameter_set_label) %>% 
  select(prop_rapid_testing, mean_increase_self_test, target_type, i_group, sim_id, t, eventually_confirmed_self, eventually_false_neg_self) %>% 
  group_by(prop_rapid_testing, mean_increase_self_test, target_type, sim_id) %>% 
  filter(t == max(t)) %>% 
  group_by(prop_rapid_testing, mean_increase_self_test, target_type, sim_id) %>% 
  summarise(across(c(eventually_confirmed_self, eventually_false_neg_self), sum)) %>% 
  mutate(prop_confirmed = eventually_confirmed_self / (eventually_false_neg_self + eventually_confirmed_self)) %>% 
  
  summarise(mean_sdl(prop_confirmed)) %>% 
  ggplot(aes(x = mean_increase_self_test, y = y, colour = target_type)) + 
  geom_point() + 
  facet_wrap(~ factor(prop_rapid_testing)) + 
  geom_errorbar(aes(ymin = ymin, ymax = ymax))
  # coord_flip()
