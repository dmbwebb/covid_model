library("viridis")
load("data/processed/inequality_is_bad.RData")
save_pdf <- function(file_path) {
  ggsave(str_glue("figures/{file_path}.pdf"), device = cairo_pdf, width = 8, height = 4.5, scale = 0.8)
}

inequality_is_bad_slim <- inequality_is_bad_slim %>% 
  mutate(tighten_factor = 1-tighten_factor)


# P values
p_baseline <- inequality_is_bad_slim %>% 
  # Amalgamate all groups
  group_by(tighten_factor, sim_id, t) %>% 
  summarise(n_pop = sum(n_pop), 
            n_cases_cum = sum(n_cases_cum)) %>% 
  
  # Look at maximum cumulative cases (for end of epidemic)
  group_by(tighten_factor, sim_id, n_pop) %>% 
  summarise(n_cases_cum = max(n_cases_cum)) %>% 
  mutate(prop_infected = n_cases_cum / n_pop) %>% 
  
  ungroup %>% 
  group_by(sim_id) %>% 
  mutate(prop_infected_baseline = prop_infected[abs(tighten_factor - 0) < 0.001]) %>% 
  mutate(lower_than_baseline = prop_infected < prop_infected_baseline) %>% 
  group_by(tighten_factor) %>% 
  summarise(
    n_lower_than_baseline = sum(lower_than_baseline),
    p_lower_than_baseline = mean(lower_than_baseline),
    total = n()
  )

write_excel_csv(p_baseline, path = "figures/inequality_is_bad_n_below_baseline.csv")


# P VALUES v2 (t test)
p_inequality <- inequality_is_bad_slim %>% 
  # Amalgamate all groups
  group_by(tighten_factor, sim_id, t) %>% 
  summarise(n_pop = sum(n_pop), 
            n_cases_cum = sum(n_cases_cum)) %>% 
  
  # Look at maximum cumulative cases (for end of epidemic)
  group_by(tighten_factor, sim_id, n_pop) %>% 
  summarise(n_cases_cum = max(n_cases_cum)) %>% 
  mutate(prop_infected = n_cases_cum / n_pop) %>% 
  ungroup

p_baseline <- p_inequality %>% filter(tighten_factor == 0)


p_inequality_t_test <- p_inequality %>% 
  filter(tighten_factor != 0) %>% 
  group_by(tighten_factor) %>% 
  select(tighten_factor, prop_infected) %>% 
  nest(data = prop_infected) %>% 
  rowwise() %>% 
  mutate(t_test = list(stats::t.test(x = data, y = p_baseline$prop_infected)),
         p = t_test$p.value) %>% 
  select(-data, -t_test) %>% 
  mutate(p_round = round(p, digits = 5))


write_excel_csv(p_inequality_t_test, path = "figures/inequality_is_bad_p_values.csv")


# PLOT EPIDEMIC CURVE
inequality_is_bad_slim %>% 
  group_by(tighten_factor, sim_id, i_group) %>% 
  arrange(tighten_factor, sim_id, i_group, t) %>% 
  mutate(new_cases_2_week = n_cases_cum - lag(n_cases_cum, 14)) %>% 
  group_by(tighten_factor, t, i_group, n_pop) %>% 
  quantile_summarise(c(n_cases_live, new_cases_2_week), conf_level = 0) %>% 
  ungroup %>% 
  mutate(incidence_prop = new_cases_2_week_median / n_pop) %>% 
  # mutate(tighten_factor = expression(lambda, "-"))
  mutate(i_group = recode_i_group(i_group)) %>% 
  
  ggplot(aes(x = t, y = incidence_prop, 
             colour = factor(i_group), group = factor(i_group))) + 
  geom_line(show.legend = TRUE) + 
  facet_wrap(. ~ tighten_factor, labeller = label_bquote(cols = lambda == .(tighten_factor))) + 
  
  # Formatting
  theme_custom() + 
  theme(legend.position = "right") + 
  scale_alpha_continuous(guide = FALSE) + 
  scale_y_continuous(labels = scales::percent_format(accuracy = 1), minor_breaks = NULL) + 
  labs(colour = "SES Group", fill = "SES Group", y = "Per capita incidence (previous 2 weeks)") +
  geom_hline(yintercept = 0, colour = "darkgrey", size = 0.4) + 
  coord_cartesian(xlim = c(0, 300)) + 
  scale_colour_viridis(option = "viridis",
                       begin = 0.3, 
                       discrete = TRUE)
  # scale_x_continuous(minor_breaks = c(0, 25, 50, 75), breaks = c(0, 50, 100))

save_pdf("inequality_is_bad_curve")



# PLOT TOTAL INFECTED in each group
plot_ineqality_is_bad_prop_infected <- inequality_is_bad_slim %>% 
  mutate(i_group = recode_i_group(i_group)) %>% 
  group_by(tighten_factor, i_group, sim_id, n_pop) %>% 
  summarise(n_cases_cum = max(n_cases_cum)) %>% 
  group_by(tighten_factor, i_group, n_pop) %>% 
  quantile_summarise(n_cases_cum, conf_level = conf_level) %>% 
  mutate(prop_infected_median = n_cases_cum_median / n_pop,
         prop_infected_upper = n_cases_cum_upper / n_pop,
         prop_infected_lower = n_cases_cum_lower / n_pop) %>% 
  # mutate(label = paste0("lambda == ", tighten_factor))) %>% 
  
  ggplot(aes(x = i_group, fill = factor(i_group))) + 
  geom_col(aes(y = prop_infected_median)) + 
  geom_errorbar(aes(ymin = prop_infected_lower, ymax = prop_infected_upper), width = 0.2, size = 0.5, colour = "#374561") + 
  facet_grid(. ~ tighten_factor, labeller = label_bquote(cols = lambda == .(tighten_factor))) + 
  labs(y = "Cumulative per capita incidence\n(entire epidemic)", x = "SES Group", fill = "SES Group") + 
  theme_custom() + 
  theme(legend.position = "right") + 
  geom_hline(yintercept = 0, colour = "darkgrey", size = 0.4) + 
  scale_fill_viridis(option = "viridis",
                     begin = 0.3, 
                     discrete = TRUE)
# scale_x_continuous(breaks = 1:6)

save_pdf("ineqality_is_bad_prop_infected")


plot_ineqality_is_bad_prop_infected + 
  geom_label_repel(aes(y = prop_infected_median, label = round(prop_infected_median, 3)), alpha = 0.8)

save_pdf("ineqality_is_bad_prop_infected_vallabels")

# PLOT TOTAL INFECTED OVER WHOLE EPIDEMIC
plot_inequality_is_bad <- inequality_is_bad_slim %>% 
  group_by(tighten_factor, i_group, sim_id, n_pop) %>% 
  summarise(n_cases_cum = max(n_cases_cum)) %>% print(n = 100) %>% 
  group_by(tighten_factor, sim_id) %>% 
  summarise(n_pop = sum(n_pop), n_cases_cum = sum(n_cases_cum)) %>% 
  mutate(prop_infected = n_cases_cum / n_pop) %>% 
  group_by(tighten_factor) %>% 
  quantile_summarise(prop_infected, conf_level = conf_level) %>% 
  # summarise_quantiles(prop_infected) %>% 
  # mutate(type = case_when(
  #   q == 0.025 ~ "lower", 
  #   q == 0.05 ~ "lower_90",
  #   q == 0.5 ~ "median",
  #   q == 0.95 ~ "upper_90",
  #   q == 0.975 ~ "upper"
  # )) %>% 
  # filter(!is.na(type)) %>% 
  # select(-q) %>% 
  # pivot_wider(names_from = type, values_from = prop_infected, names_prefix = "y_") %>% 
  
  ggplot(aes(x = tighten_factor)) + 
  # geom_col(aes(y = prop_infected_median), width = 0.1,show.legend = FALSE) + 
  geom_errorbar(aes(ymin = prop_infected_lower, ymax = prop_infected_upper), width = 0.02, colour = "#374561") + 
  # geom_errorbar(aes(ymin = y_lower_90, ymax = y_upper_90), width = 0, size = 0.9, colour = "#374561") + 
  geom_point(aes(y = prop_infected_median), size = 2.5, colour = "indianred") +
  # facet_wrap(~ tighten_factor) + 
  labs(y = "Cumulative per capita incidence\n(all SES groups, entire epidemic)", x = expression(paste(lambda, " (Proportional reduction in inequality)"))) + 
  theme_custom(panel.grid.major.x = element_blank()) + 
  theme(legend.position = "right") + 
  scale_y_continuous(labels = scales::percent)
# geom_hline(yintercept = 0, colour = "darkgrey", size = 0.4) + 
# scale_fill_viridis(option = "viridis",
#                    begin = 0.3, 
#                    discrete = TRUE) + 
# scale_x_continuous(breaks = 1:6)


plot_inequality_is_bad
save_pdf("ineqality_is_bad_prop_infected_point")



plot_inequality_is_bad + 
  geom_label(aes(y = prop_infected_median, label = round(prop_infected_median, 3)), alpha = 0.8) + 
  geom_label(aes(y = prop_infected_upper, label = round(prop_infected_upper, 3)), alpha = 0.8) + 
  geom_label(aes(y = prop_infected_lower, label = round(prop_infected_lower, 3)), alpha = 0.8)

save_pdf("ineqality_is_bad_prop_infected_point_vallabels")



# 
# inequality_is_bad_slim %>% 
#   group_by(tighten_factor, i_group, sim_id, n_pop) %>% 
#   summarise(n_cases_cum = max(n_cases_cum)) %>% print(n = 100) %>% 
#   group_by(tighten_factor, sim_id) %>% 
#   summarise(n_pop = sum(n_pop), n_cases_cum = sum(n_cases_cum)) %>% 
#   mutate(prop_infected = n_cases_cum / n_pop) %>% 
#   group_by(tighten_factor) %>% 
#   
#   ggplot(aes(x = factor(tighten_factor), y = prop_infected)) + 
#   geom_violin(width = 0.5)



