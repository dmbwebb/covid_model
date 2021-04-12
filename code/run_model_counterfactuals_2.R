# (1) Counterfactuals v2 with halfway tightening ---------------------------

# EQUALISED PARAMETERS
# params_equalised <- tighten_params(params_baseline, equalise_to = 4)

# Calculate new list of equalised parameters
params_to_equalise_list <- list(
  "baseline"             = c(),
  "out_contacts"         = c("k_matrix"),
  "home_contacts"        = c("hh_size_data", "sar_home"),
  "isolation_behaviour"  = c("probs_isolate_symptoms_df", "probs_isolate_test_df", "probs_isolate_ct_df", "p_hh_quarantine"),
  "testing_tracing"      = c("probs_self_test_df", "test_choice_delay_data", "test_results_delay_data", "p_contact_traced")
)

# REMOVE ALL - it's no longer correct because isolation and out_contacts are interdependent (need to expand on all_at_once argument if we want to make this work)
# params_to_equalise_list$all <- flatten_chr(params_to_equalise_list)

# CALCULATE TIGHTENED PARAMS [takes a while]
params_tightened <- tibble(
  # tighten_factor = c(0, 0.25, 0.5, 0.75)
  tighten_factor = c(0, 0.5)
) %>% 
  rowwise() %>% 
  mutate(params_tightened = list(tighten_params(params_baseline, tighten_factor = tighten_factor, tighten_to = 4, all_at_once = FALSE)))


# Combine params from old and new
params_tightened_df <- params_to_equalise_list %>% 
  enframe(name = "parameter_set", value = "params_to_equalise") %>% 
  crossing(params_tightened) %>% 
  rowwise() %>% 
  mutate(
    baseline_to_use = list(params_baseline[! (names(params_baseline) %in% params_to_equalise)]),
    tightened_to_use = list(params_tightened[names(params_tightened) %in% params_to_equalise]),
    params_to_use = list(c(baseline_to_use, tightened_to_use))
  ) %>% 
  select(-baseline_to_use, -tightened_to_use) %>%
  filter(!(parameter_set == "baseline" & tighten_factor != 0)) %>% 
  {
    bind_cols(
      tibble(parameter_set = .$parameter_set,
             tighten_factor = .$tighten_factor,
             params_to_equalise = .$params_to_equalise),
      .$params_to_use %>% purrr::transpose() %>% as_tibble()
    )
  }

probs_to_p <- function(x) {
  x %>% 
    map(~ filter(.x, symptom_severity == 2)) %>% 
    map("prob")
}

params_check <- params_tightened_df %>% 
  select(parameter_set, tighten_factor, matches("days_"), matches("isolate")) %>% 
  print_names %>% 
  mutate(
    i_group = list(recode_i_group(1:4)),
    probs_isolate_ct_df = probs_to_p(probs_isolate_ct_df),
    probs_isolate_symptoms_df = probs_to_p(probs_isolate_symptoms_df),
    probs_isolate_test_df = probs_to_p(probs_isolate_test_df)
  ) %>% 
  unnest(-c(parameter_set, tighten_factor)) %>% 
  relocate(parameter_set, i_group, tighten_factor)

params_check %>% 
  filter(parameter_set %in% c("baseline", "isolation_behaviour", "all")) %>% 
  pivot_longer(days_of_work:probs_isolate_test_df) %>% 
  filter(!str_detect(name, "days_")) %>% 
  ggplot(aes(x = i_group, y = value, colour = factor(tighten_factor))) + 
  geom_point() + 
  # geom_line() + 
  facet_grid(name ~ parameter_set, scales = "free_y")



# # Convert to tibble to be pmapped
# params_combined_df <- bind_cols(
#   tibble(parameter_set = params_combined$parameter_set,
#          tighten_factor = params_combined$tighten_factor,
#          params_to_equalise = params_combined$params_to_equalise),
#   params_combined$params_to_use %>% purrr::transpose() %>% as_tibble()
# ) %>% 
#   filter(tighten_factor %in% c(0, 0.5)) # fewer sims to start with 

# RUN THE MODELS
RNGkind("L'Ecuyer-CMRG"); set.seed(12345)
group_diffs_tightened <- params_tightened_df %>% 
  # filter(tighten_factor != 0) %>% 
  mutate(model_output = pmap(select(., -parameter_set, -params_to_equalise, -tighten_factor, -days_of_work, -starts_with("days_work")),
                             group_diff, 
                             print_detail = FALSE,
                             n_cores = n_cores,
                             n_sims = n_sims, 
                             n_iterations = n_iterations))

# system("say R has finished running")

# Slim the data to just time series to be used in graphs
group_diff_slim_tightened <- group_diffs_tightened %>% 
  select(parameter_set, tighten_factor, any_of("params_to_equalise"), model_output) %>% 
  # mutate(time_series = map(model_output, "time_series")) %>% 
  # select(-model_output) %>% 
  unnest(model_output)
# mutate(parameter_set = factor(parameter_set, levels = names(params_to_equalise_list)))


save(group_diff_slim_tightened, file = "data/processed/group_diff_tightened_half.RData")







# (2) Split out parameters ----------------------------------------------------

# SPLIT OUT HOME VS NOT HOME
params_split <- list(
  "hh_size" = c("hh_size_data"),
  "sar_home" = c("sar_home")
) %>% 
  enframe(name = "parameter_set", value = "params_to_equalise") %>% 
  crossing(params_tightened) %>% 
  rowwise() %>% 
  mutate(
    baseline_to_use = list(params_baseline[! (names(params_baseline) %in% params_to_equalise)]),
    tightened_to_use = list(params_tightened[names(params_tightened) %in% params_to_equalise]),
    params_to_use = list(c(baseline_to_use, tightened_to_use))
  ) %>% 
  select(-baseline_to_use, -tightened_to_use) %>%
  {bind_cols(
      tibble(parameter_set = .$parameter_set,
             tighten_factor = .$tighten_factor,
             params_to_equalise = .$params_to_equalise),
      .$params_to_use %>% purrr::transpose() %>% as_tibble()
    )}

# RUN THE MODELS
RNGkind("L'Ecuyer-CMRG"); set.seed(12345)
group_diffs_split <- params_split %>% 
  mutate(model_output = pmap(select(., -parameter_set, -params_to_equalise, -tighten_factor, -days_of_work, -starts_with("days_work")),
                             group_diff, 
                             print_detail = FALSE,
                             n_cores = n_cores,
                             n_sims = n_sims, 
                             n_iterations = n_iterations))

group_diff_slim_split <- group_diffs_split %>% 
  select(parameter_set, tighten_factor, any_of("params_to_equalise"), model_output) %>%
  unnest(model_output)

save(group_diff_slim_split, file = "data/processed/group_diff_split.RData")



