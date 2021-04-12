
# Get the data for consistency checks -------------------------------------

# .............. ---------------------------------------------------------------------
# CONSISTENCY CHECKS ------------------------------------------------------

# source("")
# load("data/processed/data_save.RData")

consistency_check_data <- outbreak_sims(
  keep_all_data = TRUE,
  n_sims = 1,
  n_iterations = 50,
  print_detail = TRUE,
  n_pop = 5000, 
  hh_size_data = data_save$hh_data_bogota,
  n_initial_cases = c(5, 5, 5, 5), # calculated within function if remains NULL
  group_props = data_save$group_props,     # ditto
  dt_approx = 1,
  recov_val = 10, # people recover 10 days after symptoms show (cuts off only ~ 1% of infections at the upper end of timing distribution)
  params_timing = data_save$params_timing,
  params_symptom_timing = data_save$params_symptom_timing,
  params_serial = data_save$params_serial,
  test_delay_data = data_save$test_delay_data,
  # test_choice_delay_data = data_save$test_choice_delay_data,
  # test_results_delay_data = data_save$test_results_delay_data, 
  test_sensitivity_data = data_save$test_sensitivity_data, 
  ct_delay_data = data_save$ct_delay_data,
  probs_self_test_df = data_save$probs_self_test,
  probs_isolate_symptoms_df = data_save$probs_isolate_symptoms,
  probs_isolate_test_df = data_save$probs_isolate_test,
  probs_isolate_ct_df =   data_save$probs_isolate_ct,
  p_contact_if_isolated_home = rep(1, 4), 
  p_contact_traced = data_save$params_data$p_contact_traced,
  p_hh_quarantine = data_save$params_data$p_hh_quarantine,
  k_matrix = data_save$k_matrix * (5000 / data_save$k_matrix_pop) * 2, # multiply by new n_pop here,
  contact_dispersion = data_save$contact_dispersion,
  # beta_matrix = beta_matrix_6,
  # r0_group_1 = 5, 
  # contacts_mean = contacts_mean_6,
  sar_out = data_save$params_data$sar_out_input,
  sar_home = data_save$params_data$sar_home_input, 
  infectiousness_by_symptoms =  data_save$infectiouness_by_symptoms,
  alpha = c(0, 0, 0, 0)
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
# library("testthat")
# test_file("unit_tests.R")


# Run this after generating the consistency checks stuff in individual_model_testing
# and then use 
# test_file("unit_tests.R")