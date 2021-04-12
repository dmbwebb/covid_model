
# ....................... -------------------------------------------------

# NEW TARGETED POLICY ---------------------------------------------------------------------

# NEW VERSION where reduction in params is uniform across targeted groups
target_params <- function(params_original, n_pops, budget, target_groups = 1,
                          direction = 1) {
  
  # TESTING
  # params_original <- rep(0.1, 6)
  # n_pops <- group_diff_default(n_sims = 1, n_iterations = 2)$time_series$n_pop[1:6]
  # budget <- 10000 * 0.04
  # target_groups <- 1:2
  # direction <- -1
  
  if (length(n_pops) != length(params_original)) stop("number of groups is not same across n_pops and params_original")
  
  (per_capita_untargeted <- budget / sum(n_pops)) # reduction across all groups
  (per_capita_targeted <- budget / sum(n_pops[target_groups])) # reduction if only on 1:2 groups
  
  targeted_indices <- if_else(seq_along(n_pops) %in% target_groups, 1, 0)
  
  out <-  tibble(
    untargeted = params_original + (per_capita_untargeted * direction),
    targeted = params_original + (targeted_indices * per_capita_targeted * direction)
  ) %>% print
  
  return(out)
  
}


target_groups <- 1

# Calculate contact means
n_pops <- (data_save$k_matrix_pop * data_save$group_props)


# Contacts outside home ---------------------------------------------------



contact_means_example <- rowSums(data_save$k_matrix * data_save$k_scale_factor)  / n_pops

target_contacts <- tibble(
  budget = c(0, 10, 20, 30) * 1000,
  mean_reduc_per_person = budget / data_save$k_matrix_pop
) %>% 
  rowwise() %>% 
  mutate(target_params = list(target_params(params_original = contact_means_example,
                                            n_pops = n_pops,
                                            budget = budget, 
                                            target_groups = target_groups,
                                            direction = -1))) %>% 
  unnest(target_params) %>% 
  group_by(budget, mean_reduc_per_person) %>% 
  summarise(untargeted = list(untargeted),
            targeted = list(targeted)) %>% 
  pivot_longer(-c(budget, mean_reduc_per_person), names_to = "target_type", values_to = "contact_means")


# target_contacts %>% unnest(contact_means) %>% print_all





target_k <- function(k_matrix, contact_means_aim, pops, which_groups) {
  
  if (length(contact_means_aim) == 0) {  
    print(k_matrix %>% round())
    return(k_matrix)
  }
  
  n_groups <- length(contact_means_aim)
  
  rescale_factor_target <- nloptr::nloptr(x0 = rep(1, n_groups), 
                                          eval_f = find_the_rescale_factor_vec, 
                                          lb = rep(0, n_groups) + 1e-5,
                                          opts = list(
                                            "algorithm" = "NLOPT_LN_SBPLX",
                                            "xtol_abs"=1.0e-10,
                                            "maxeval" = 1000,
                                            "print_level" = 0
                                          ),
                                          contact_means_aim = contact_means_aim,
                                          k_matrix = k_matrix, 
                                          pops = pops,
                                          which_groups = which_groups)  
  
  rescale_factors <- rescale_factor_target$solution
  # print(rescale_factors)
  
  rescale_factors_all <- rep(1, nrow(k_matrix))
  rescale_factors_all[which_groups] <- rescale_factors
  
  k_rescaled <- rescale_k(k_matrix, rescale_factors_all)
  
  return(k_rescaled)
}



target_k(k_matrix = data_save$k_matrix,
         contact_means_aim = c(20, 10, 8, 7),
         pops = n_pops,
         which_groups = 1:4)


# Generate the K matrix
RNGkind("L'Ecuyer-CMRG"); set.seed(12345)
target_k_list <- target_contacts %>% 
  # ONLY WANT VALUES TO OPTIMISE, remove contact_mean values that should be unchanged
  mutate(which_groups = case_when(target_type == "targeted" ~ list(target_groups),
                                  target_type == "untargeted" ~ list(1:4))) %>% 
  mutate(contact_means_original = list(contact_means_example)) %>% 
  filter(!(target_type == "targeted" & budget == 0)) %>% 
  rowwise() %>% 
  mutate(contact_means = list(contact_means[abs(contact_means - contact_means_original) > 0.01])) %>%  # only want to optimise over contact menas that change
  select(-contact_means_original) %>% 
  mutate(k_matrix = list(target_k(k_matrix = data_save$k_matrix * data_save$k_scale_factor,
                                  contact_means_aim = contact_means,
                                  pops = n_pops,
                                  which_groups = which_groups))) %>% 
  ungroup



# target_k_list$k_matrix[[4]] %>% round()
# (data_save$k_matrix * data_save$k_scale_factor) %>% round()

RNGkind("L'Ecuyer-CMRG"); set.seed(12345)
target_results <- target_k_list %>% 
  mutate(model_output = pmap(select(., k_matrix), 
                             group_diff, 
                             n_cores = n_cores,
                             n_sims = n_sims, print_detail = FALSE, n_iterations = n_iterations)) # k_scale_factor is already accounted for in process above

# *** NEED TO DO 50 OF EACH ****

target_slim <- target_results %>% 
  # mutate(model_output = map(model_output, "time_series")) %>% 
  unnest(model_output)

target_slim %>% group_by(mean_reduc_per_person) %>% summarise(n = n_distinct(sim_id))


save(target_slim, file = "data/processed/target_slim.RData")

# load("data/processed/target_slim.RData")





# NEED TO MAKE SURE THAT rescale factor is one when it's unchanged

# CREATE 2D by crossing
# target_crossed <- full_join_track(
#   target_params_list %>% rename(budget_sar_out = budget),
#   target_k_list %>% rename(budget_contacts = budget)
# ) %>% 
#   rowwise() %>% 
#   mutate(contact_means_check = list(rowSums(k_matrix / n_pops))) %>% 
#   ungroup
# 
# target_crossed %>% unnest(c(sar_out, contact_means_check)) %>% 
#   print(n = 100)
# 
# target_crossed_results <- target_crossed %>% 
#   mutate(model_output = pmap(select(., sar_out, k_matrix), group_diff_default, n_sims = 1, n_iterations = 10000))
# 
# target_crossed_slim <- target_crossed_results %>% 
#   mutate(model_output = map(model_output, "time_series")) %>% 
#   unnest(model_output)
# 
# save(target_crossed_slim, file = "data/processed/target_crossed_slim.RData")





# Isolation behaviour -----------------------------------------------------


p_ct <- data_save$probs_isolate_ct %>% select(i_group, prob) %>% dups_drop() %>% .$prob
p_symp <- data_save$probs_isolate_symptoms %>% select(i_group, prob) %>% filter(prob != 0) %>% .$prob
# data_save$probs_isolate_test

target_isolation <- tibble(
  budget = c(0, 10, 20, 25) * 100,
  mean_increase_per_person = budget / data_save$k_matrix_pop
) %>% 
  rowwise() %>% 
  mutate(p_isolate_ct = list(target_params(params_original = p_ct,
                                           n_pops = n_pops,
                                           budget = budget, 
                                           target_groups = target_groups,
                                           direction = 1)),
         p_isolate_symp = list(target_params(params_original = p_symp,
                                             n_pops = n_pops,
                                             budget = budget, 
                                             target_groups = target_groups,
                                             direction = 1))) %>% 
  unnest(c(p_isolate_ct, p_isolate_symp), names_repair = "universal") %>% 
  set_names(c("budget", "mean_increase_per_person", "p_isolate_ct_untargeted", "p_isolate_ct_targeted", "p_isolate_symp_untargeted", "p_isolate_symp_targeted")) %>% 
  mutate(i_group = rep(1:4, times = nrow(.) / 4)) %>%
  pivot_longer(-c(budget, mean_increase_per_person, i_group), names_to = c(".value", "target_type"), names_pattern = "(p_.*)_(u?n?targeted)$") %>% 
  group_by(budget, mean_increase_per_person, target_type) %>% 
  arrange(budget, mean_increase_per_person, target_type) %>% 
  summarise(probs_isolate_ct = list(data_save$probs_isolate_ct %>% mutate(prob = rep(p_isolate_ct, each = 2))),
            probs_isolate_symptoms = list(data_save$probs_isolate_symptoms %>% 
                                            mutate(prob = rep(p_isolate_symp, each = 2)) %>% 
                                            mutate(prob = if_else(symptom_severity == 1L, 0, prob)))) %>% 
  ungroup


RNGkind("L'Ecuyer-CMRG"); set.seed(12345)
target_isolation_results <- target_isolation %>% 
  mutate(model_output = pmap(select(., probs_isolate_ct, probs_isolate_symptoms), 
                             group_diff,
                             n_cores = n_cores,
                             n_sims = n_sims, print_detail = FALSE, n_iterations = n_iterations))


target_isolation_slim <- target_isolation_results %>% 
  # mutate(model_output = map(model_output, "time_series")) %>% 
  unnest(model_output)

target_isolation %>% unnest(c(probs_isolate_symptoms)) %>% print_all

save(target_isolation_slim, file = "data/processed/target_isolation_slim.RData")









# Immunity / vaccines -----------------------------------------------------




target_immunity <- tibble(
  budget = c(0, 10, 20, 30, 40) * 100,
  # budget = c(10, 30) * 100,
  prop_immune = budget / data_save$k_matrix_pop
) %>% 
  rowwise() %>% 
  mutate(p_immune = list(target_params(params_original = c(0, 0, 0, 0),
                                       n_pops = n_pops,
                                       budget = budget, 
                                       target_groups = target_groups,
                                       direction = 1))) %>% 
  unnest(p_immune) %>% 
  mutate(i_group = rep(1:4, times = nrow(.) / 4)) %>%
  pivot_longer(c(untargeted, targeted), names_to = "target_type", values_to = "p_immune") %>% 
  group_by(budget, prop_immune, target_type) %>% 
  arrange(budget, prop_immune, target_type) %>% 
  select(-i_group) %>% 
  summarise(p_immune =  list(p_immune)) %>% 
  ungroup

target_immunity %>% unnest(p_immune) %>% print_all



RNGkind("L'Ecuyer-CMRG"); set.seed(12345)
target_immunity_results <- target_immunity %>% 
  mutate(model_output = pmap(select(., p_immune), 
                             group_diff, 
                             n_cores = n_cores,
                             n_sims = n_sims, print_detail = FALSE, n_iterations = n_iterations))

target_immunity_slim <- target_immunity_results %>% 
  # mutate(model_output = map(model_output, "time_series")) %>% 
  unnest(model_output)



save(target_immunity_slim, file = "data/processed/target_immunity_slim.RData")





# Rapid testing counterfactual --------------------------------------------

# Fixed number of people who have their testing delays set to 0?
# Testing choice delay ≈ 1/2 days
# Testing results delay = 0.05 days

# Caveats - this won't properly account for false positives in the "budget"



# INCREASE P(SELF TEST) as well?

# data_save$probs_self_test


# # RAPID TESTING
# target_rapid_testing <- tibble(
#   budget_rapid = c(0, 0.2, 0.4) * 10000,
#   total_prop_with_rapid_testing = budget_rapid / data_save$k_matrix_pop
# ) %>% 
#   rowwise() %>% 
#   mutate(target_params = list(target_params(params_original = c(0, 0, 0, 0), # no-one has rapid testing initially
#                                             n_pops = n_pops,
#                                             budget = budget_rapid, 
#                                             target_groups = target_groups,
#                                             direction = 1))) %>% 
#   unnest(target_params) %>% 
#   group_by(budget_rapid, total_prop_with_rapid_testing) %>% 
#   summarise(untargeted = list(untargeted),
#             targeted = list(targeted)) %>% 
#   pivot_longer(-c(budget_rapid, total_prop_with_rapid_testing), names_to = "target_type", values_to = "prop_rapid_testing") %>% 
#   print
# 
# # INCREASE P(SELF TEST)
# # self_testing_increase <- tibble(
# #   budget_self_test = c(0, 0.3) * 10000,
# #   mean_increase_self_test = budget_self_test / data_save$k_matrix_pop
# # ) %>% 
# #   rowwise() %>% 
# #   mutate(target_params = list(target_params(params_original = data_save$probs_self_test %>% filter(symptom_severity == 2) %>% .$prob, # no-one has rapid testing initially
# #                                             n_pops = n_pops,
# #                                             budget = budget_self_test, 
# #                                             target_groups = target_groups,
# #                                             direction = 1))) %>% 
# #   unnest(target_params) %>% 
# #   group_by(budget_self_test, mean_increase_self_test) %>% 
# #   summarise(untargeted = list(untargeted),
# #             targeted = list(targeted)) %>% 
# #   pivot_longer(-c(budget_self_test, mean_increase_self_test), names_to = "target_type", values_to = "p_self_test") %>% 
# #   filter(target_type == "untargeted") %>% 
# #   select(-target_type) %>% 
# #   print
# self_testing_increase <- tibble(
#   self_test_type = c("baseline", "high"),
#   probs_self_test = list(
#     data_save$probs_self_test,
#     data_save$probs_self_test %>% mutate(prob = if_else(prob != 0, prob + 0.3, prob))
#   )
# )
# 
# 
# # COMBINE BOTH:
# rapid_plus_self_testing <- crossing(target_rapid_testing, self_testing_increase) %>% 
#   filter(!(budget_rapid == 0 & target_type == "targeted"))
# 
# 
# 
# 
# 
# # prop_rapid_testing <- rep(0.2, 5)
# # test_delay_data <- data_save$test_delay_data
# # rapid_choice_val <- 0.1
# # rapid_results_val <- 0.2
# 
# # Convert prop_rapid_testing into the testing_delays data
# gen_rapid_testing_data <- function(prop_rapid_testing, test_delay_data, rapid_choice_val, rapid_results_val) {
#   n_groups <- length(prop_rapid_testing)
#   
#   # link the i_groups to their prob of being rapid tested
#   rapid_testing_tibble <- tibble(
#     i_group = 1:n_groups,
#     prop_rapid_testing = prop_rapid_testing
#   )
#   
#   rapid_delay_data <- test_delay_data %>%
#     left_join(rapid_testing_tibble, by = "i_group") %>%
#     mutate(
#       rapid_test_yn = bern_vec(nrow(.), p = prop_rapid_testing),
#       test_choice_delay = if_else(rapid_test_yn, rapid_choice_val, test_choice_delay),
#       test_results_delay = if_else(rapid_test_yn, rapid_results_val, test_results_delay),
#     )
#   
#   return(rapid_delay_data)
#   
# }
# 
# 
# target_rapid_testing_data <- rapid_plus_self_testing %>% 
#   rowwise() %>% 
#   mutate(
#     test_delay_data = list(gen_rapid_testing_data(prop_rapid_testing = prop_rapid_testing,
#                                                   test_delay_data = data_save$test_delay_data,
#                                                   rapid_choice_val = 1,
#                                                   rapid_results = 0.05))
#   ) %>% 
#   ungroup
# 
# 
# 
# # n_cores <- 1; n_sims <- 1; n_iterations <- 1
# # RUN THE MODELS
# rapid_testing_results <- target_rapid_testing_data %>% 
#   mutate(model_output = pmap(select(., -c(budget_rapid, target_type, prop_rapid_testing, total_prop_with_rapid_testing, self_test_type)), 
#                              group_diff, 
#                              n_cores = n_cores,
#                              n_sims = n_sims, 
#                              n_iterations = n_iterations))
# 
# # rapid_testing_results_2 <- target_rapid_testing_data %>% 
# #   mutate(model_output = pmap(select(., -c(budget, target_type, prop_rapid_testing)), group_diff_default, n_sims = 5, n_iterations = 10000))
# 
# 
# # rapid_testing_slim <- rapid_testing_results %>% 
# #   mutate(model_output = map(model_output, "time_series")) %>% 
# #   unnest(model_output)
# 
# rapid_testing_slim <- rapid_testing_results %>% 
#   select(-test_delay_data) %>% 
#   # mutate(model_output = map(model_output, "time_series")) %>% 
#   unnest(model_output)
# 
# 
# 
# save(rapid_testing_slim, file = "data/processed/rapid_testing_slim_v2.RData")





# "FAST" testing counterfactuals --------------------------------------------

# Fixed number of people who have their testing delays set to 0?
# Testing choice delay ≈ 1/2 days
# Testing results delay = 0.05 days

# Caveats - this won't properly account for false positives in the "budget"

rapid_choice_val <-  1
rapid_results <-  1

# INCREASE P(SELF TEST) as well?

# data_save$probs_self_test


# RAPID TESTING
# target_rapid_testing <- tibble(
#   budget_rapid = c(0, 0.2, 0.4) * 10000,
#   total_prop_with_rapid_testing = budget_rapid / data_save$k_matrix_pop
# ) %>% 
#   rowwise() %>% 
#   mutate(target_params = list(target_params(params_original = c(0, 0, 0, 0), # no-one has rapid testing initially
#                                             n_pops = n_pops,
#                                             budget = budget_rapid, 
#                                             target_groups = target_groups,
#                                             direction = 1))) %>% 
#   unnest(target_params) %>% 
#   group_by(budget_rapid, total_prop_with_rapid_testing) %>% 
#   summarise(untargeted = list(untargeted),
#             targeted = list(targeted)) %>% 
#   pivot_longer(-c(budget_rapid, total_prop_with_rapid_testing), names_to = "target_type", values_to = "prop_rapid_testing")



# Convert prop_rapid_testing into the testing_delays data
gen_rapid_testing_data <- function(prop_rapid_testing, test_delay_data, rapid_choice_val, rapid_results_val) {
  
  # n_groups <- length(prop_rapid_testing)
  # 
  # # link the i_groups to their prob of being rapid tested
  # rapid_testing_tibble <- tibble(
  #   i_group = 1:n_groups,
  #   prop_rapid_testing = prop_rapid_testing
  # )
  
  rapid_delay_data <- test_delay_data %>%
    # left_join(rapid_testing_tibble, by = "i_group") %>%
    mutate(
      rapid_test_yn = bern_vec(nrow(.), p = prop_rapid_testing),
      test_choice_delay = if_else(rapid_test_yn, rapid_choice_val, test_choice_delay),
      test_results_delay = if_else(rapid_test_yn, rapid_results_val, test_results_delay),
    )
  
  return(rapid_delay_data)
  
}

RNGkind("L'Ecuyer-CMRG"); set.seed(12345)
rapid_testing_increase <- tibble(
  prop_rapid_testing = c(0, 0.5, 1)
) %>% 
  rowwise() %>% 
  mutate(
    test_delay_data = list(gen_rapid_testing_data(prop_rapid_testing = prop_rapid_testing,
                                                  test_delay_data = data_save$test_delay_data,
                                                  rapid_choice_val = rapid_choice_val,
                                                  rapid_results = rapid_results))
  ) %>% 
  ungroup


# INCREASE P(SELF TEST)
p_to_probs <- function(p, probs) {
  # probs <- data_save$probs_self_test; p <- c(0.1, 0.2, 0.3, 0.4)
  stopifnot(nrow(probs) == length(p) * 2)
  probs$prob[probs$symptom_severity == 2] <- p
  return(probs)
}


self_testing_increase <- tibble(
  budget_self_test = c(0, 0.1, 0.2, 0.3, 0.4) * 10000,
  mean_increase_self_test = budget_self_test / data_save$k_matrix_pop
) %>%
  rowwise() %>%
  mutate(target_params = list(target_params(params_original = data_save$probs_self_test %>% filter(symptom_severity == 2) %>% .$prob, # no-one has rapid testing initially
                                            n_pops = n_pops,
                                            budget = budget_self_test,
                                            target_groups = target_groups,
                                            direction = 1))) %>%
  unnest(target_params) %>%
  group_by(budget_self_test, mean_increase_self_test) %>%
  summarise(untargeted = list(untargeted),
            targeted = list(targeted)) %>%
  pivot_longer(-c(budget_self_test, mean_increase_self_test), names_to = "target_type", values_to = "p_self_test") %>%
  rowwise() %>% 
  mutate(
    probs_self_test = list(p_to_probs(p_self_test, data_save$probs_self_test))
  ) %>% 
  ungroup %>% 
  print

# self_testing_increase$p_self_test[[1]]
# self_testing_increase$probs_self_test[[1]]



# self_testing_increase <- tibble(
#   self_test_type = c("baseline", "high"),
#   probs_self_test = list(
#     data_save$probs_self_test,
#     data_save$probs_self_test %>% mutate(prob = if_else(prob != 0, prob + 0.3, prob))
#   )
# )






# COMBINE BOTH:
rapid_plus_self_testing <- crossing(rapid_testing_increase, self_testing_increase) %>% 
  filter(!(budget_self_test == 0 & target_type == "targeted"))





# prop_rapid_testing <- rep(0.2, 5)
# test_delay_data <- data_save$test_delay_data
# rapid_choice_val <- 0.1
# rapid_results_val <- 0.2








# n_cores <- 1; n_sims <- 1; n_iterations <- 1
# RUN THE MODELS
RNGkind("L'Ecuyer-CMRG"); set.seed(12345)
rapid_testing_results <- rapid_plus_self_testing %>% 
  mutate(model_output = pmap(select(., -c(prop_rapid_testing, budget_self_test, mean_increase_self_test, target_type, p_self_test)), 
                             group_diff, 
                             n_cores = n_cores,
                             n_sims = n_sims, 
                             n_iterations = n_iterations))

# rapid_testing_results_2 <- target_rapid_testing_data %>% 
#   mutate(model_output = pmap(select(., -c(budget, target_type, prop_rapid_testing)), group_diff_default, n_sims = 5, n_iterations = 10000))


# rapid_testing_slim <- rapid_testing_results %>% 
#   mutate(model_output = map(model_output, "time_series")) %>% 
#   unnest(model_output)

rapid_testing_slim <- rapid_testing_results %>% 
  select(-test_delay_data) %>% 
  # mutate(model_output = map(model_output, "time_series")) %>% 
  unnest(model_output)



save(rapid_testing_slim, file = "data/processed/fast_testing_slim.RData")




# ................ --------------------------------------------------------

# OLD ---------------------------------------------------------------------



# # TARGETED POLICIES -------------------------------------------------------
# 
#     
#     
#     
#     
# # Targeted policies 1D ----------------------------------------------------
#     
#   
#     
# 
#     
#     
#     # COULD REWRITE with X = untargeted reduction per person...
#     # i.e. untargeted we reduce sar_out by 0.05, what's the equivalent reduction if 
#     # we expend the same "effort" on reducing only [assuming constant cost per person...]
#     
#     n_pops <- group_diff(n_sims = 1, n_iterations = 2, n_initial_cases = c(2, 2, 2, 2))$time_series$n_pop[1:4]
#     
#     target_params(params_original = c(0.1, 0.1, 0.1, 0.1),
#                   n_pops = n_pops,
#                   budget = 25, 
#                   target_groups = 1,
#                   direction = -1)
#     
#     
#     # Make a dataframe with parameters
#     target_params_list <- tibble(
#       budget = c(0, 25, 45) * 10
#     ) %>% 
#       rowwise() %>% 
#       mutate(target_params = list(target_params(params_original = c(0.1, 0.1, 0.1, 0.1),
#                                                 n_pops = n_pops,
#                                                 budget = budget, 
#                                                 target_groups = 1,
#                                                 direction = -1))) %>% 
#       unnest(target_params) %>% 
#       group_by(budget) %>% 
#       summarise(untargeted = list(untargeted),
#                 targeted = list(targeted)) %>% 
#       pivot_longer(-budget, names_to = "target_type", values_to = "sar_out") %>% 
#       print
#     
#     target_params_list %>% unnest(sar_out) %>% print_all
#     
#     
#     target_results <- target_params_list %>% 
#       mutate(model_output = pmap(select(., -c(budget, target_type)), group_diff_default, n_sims = 1, print_detail = TRUE, n_iterations = 1000))
#     
#     target_results_2 <- target_params_list %>% 
#       mutate(model_output = pmap(select(., -c(budget, target_type)), group_diff_default, n_sims = 5, print_detail = FALSE, n_iterations = 1000))
#     
#     target_results_slim <- bind_rows(
#       target_results %>% 
#         mutate(model_output = map(model_output, "time_series")) %>% 
#         unnest(model_output) %>% 
#         mutate(sim_id = as.character(as.numeric(sim_id) + 5)),
#       target_results_2 %>% 
#         mutate(model_output = map(model_output, "time_series")) %>% 
#         unnest(model_output)
#     )
#     
#     save(target_results_slim, file = "data/processed/target_results_slim.RData")
#       
#     # load("OLD/target_results_slim.RData")
#     
#     # PLOT budget, targeted, etc. [PROP PEOPLE INFECTED!!]
#     target_results_slim %>% 
#       
#       # SHADY
#       # filter(!(target_type == "targeted" & sim_id == 5)) %>% 
#       group_by(budget, target_type, sim_id, t) %>% 
#       mutate(n_pop_total = sum(n_pop)) %>% 
#       mutate(budget = budget / n_pop_total) %>% 
#       ungroup %>% 
#       
#       # Calculate total infected in each run
#       group_by(budget, target_type, sim_id, i_group) %>% 
#       summarise(n_cases_cum = max(n_cases_cum)) %>% 
#       
#       # Sum up each i_group
#       summarise(n_cases_cum = sum(n_cases_cum)) %>% 
#       # print_all
#       
#       # Calculate median and confidence intervals
#       group_by(budget, target_type) %>% 
#       quantile_summarise(n_cases_cum, conf_level = 0.95) %>% 
#       
#       # PLOT
#       ggplot(aes(x = budget, y = n_cases_cum_median, colour = target_type)) + 
#       geom_line() + 
#       geom_errorbar(aes(ymin = n_cases_cum_lower, ymax = n_cases_cum_upper), width = 0.001, size = 0.5) + 
#       geom_point()
#     # labs(x = "Untargeted per capita reduction in SAR(out)",
#     #      y = "Total cases (Pop = 1000)")
#     
#     
#     
#     
#     
# # Target policies 2D ------------------------------------------------------
#     
#     target_groups <- 1
#     
#     
#     # Calculate contact means
#     n_pops <- (data_save$k_matrix_pop * data_save$group_props)
#     contact_means_example <- rowSums(data_save$k_matrix * data_save$k_scale_factor) / n_pops
#     
#     target_contacts <- tibble(
#       budget = c(0, 10, 20, 40) * 1000,
#     ) %>% 
#       rowwise() %>% 
#       mutate(target_params = list(target_params(params_original = contact_means_example,
#                                                 n_pops = n_pops,
#                                                 budget = budget, 
#                                                 target_groups = target_groups,
#                                                 direction = -1))) %>% 
#       unnest(target_params) %>% 
#       group_by(budget) %>% 
#       summarise(untargeted = list(untargeted),
#                 targeted = list(targeted)) %>% 
#       pivot_longer(-budget, names_to = "target_type", values_to = "contact_means")
#     
#     
#     # target_contacts %>% unnest(contact_means) %>% print_all
#     
#     
#     
#     
#     
#     target_k <- function(k_matrix, contact_means_aim, pops, which_groups) {
#       
#       if (length(contact_means_aim) == 0) return(k_matrix)
#       n_groups <- length(contact_means_aim)
#       
#       rescale_factor_target <- nloptr::nloptr(x0 = rep(1, n_groups), 
#                                               eval_f = find_the_rescale_factor_vec, 
#                                               lb = rep(0, n_groups) + 1e-5,
#                                               opts = list(
#                                                 "algorithm" = "NLOPT_LN_SBPLX",
#                                                 "xtol_abs"=1.0e-10,
#                                                 "maxeval" = 1000,
#                                                 "print_level" = 0
#                                               ),
#                                               contact_means_aim = contact_means_aim,
#                                               k_matrix = k_matrix, 
#                                               pops = pops,
#                                               which_groups = which_groups)  
#       
#       rescale_factors <- rescale_factor_target$solution
#       # print(rescale_factors)
#       
#       rescale_factors_all <- rep(1, nrow(k_matrix))
#       rescale_factors_all[which_groups] <- rescale_factors
#       
#       k_rescaled <- rescale_k(k_matrix, rescale_factors_all)
#       
#       return(k_rescaled)
#     }
#     
#     target_k(k_matrix = data_save$k_matrix,
#              contact_means_aim = c(20, 10, 8, 7),
#              pops = n_pops,
#              which_groups = 1:4)
#     
#     
#     # Generate the K matrix
#     target_k_list <- target_contacts %>% 
#       # ONLY WANT VALUES TO OPTIMISE, remove contact_mean values that should be unchanged
#       mutate(which_groups = case_when(target_type == "targeted" ~ list(target_groups),
#                                       target_type == "untargeted" ~ list(1:4))) %>% 
#       mutate(contact_means_original = list(contact_means_example)) %>% 
#       filter(!(target_type == "targeted" & budget == 0)) %>% 
#       rowwise() %>% 
#       mutate(contact_means = list(contact_means[abs(contact_means - contact_means_original) > 0.01])) %>%  # only want to optimise over contact menas that change
#       select(-contact_means_original) %>% 
#       mutate(k_matrix = list(target_k(k_matrix = (data_save$k_matrix * data_save$k_scale_factor),
#                                       contact_means_aim = contact_means,
#                                       pops = n_pops,
#                                       which_groups = which_groups)))
#     
#     
#     
#     # NEED TO MAKE SURE THAT rescale factor is one when it's unchanged
#     
#     # CREATE 2D by crossing
#     target_crossed <- full_join_track(
#       target_params_list %>% rename(budget_sar_out = budget),
#       target_k_list %>% rename(budget_contacts = budget)
#     ) %>% 
#       rowwise() %>% 
#       mutate(contact_means_check = list(rowSums(k_matrix / n_pops))) %>% 
#       ungroup
#     
#     target_crossed %>% unnest(c(sar_out, contact_means_check)) %>% 
#       print(n = 100)
#     
#     target_crossed_results <- target_crossed %>% 
#       mutate(model_output = pmap(select(., sar_out, k_matrix), group_diff_default, n_sims = 1, n_iterations = 10000))
#     
#     target_crossed_slim <- target_crossed_results %>% 
#       mutate(model_output = map(model_output, "time_series")) %>% 
#       unnest(model_output)
#     
#     save(target_crossed_slim, file = "data/processed/target_crossed_slim.RData")
#     
#     # NEXT - plot the graph
#     library("plotly")
#     
#     
#     
#     median_infected <- function(outbreak_sim_list) {
#       outbreak_sim_list$time_series %>% 
#         # Calculate total infected in each run
#         group_by(sim_id, i_group) %>% 
#         summarise(n_cases_cum = max(n_cases_cum)) %>%
#         # Sum up each i_group
#         summarise(n_cases_cum = sum(n_cases_cum)) %>% 
#         # Calculate median and confidence intervals
#         quantile_summarise(n_cases_cum, conf_level = 0.95) %>% 
#         .$n_cases_cum_median
#     }
#     
#     target_crossed_summ <- target_crossed_results %>% 
#       rowwise() %>% 
#       mutate(contact_means_check = paste0(round(contact_means_check), collapse = "-"),
#              sar_out = paste0(round(sar_out, 2), collapse = "-")) %>% 
#       select(-contact_means, -which_groups, -k_matrix) %>% 
#       mutate(
#         infected = median_infected(model_output)
#       ) %>% 
#       select(-model_output) %>% 
#       print_all
#     
#     
#     
#     
#     # 
#     # target_crossed_summ <- target_crossed_slim %>% 
#     #   group_by(budget_contacts, budget_sar_out, target_type, sim_id, t) %>% 
#     #   # mutate(n_pop_total = sum(n_pop)) %>% 
#     #   # mutate(budget = budget / n_pop_total) %>% 
#     #   ungroup %>% 
#     #   
#     #   # Calculate total infected in each run
#     #   group_by(budget_contacts, budget_sar_out, target_type, sim_id, i_group) %>% 
#     #   summarise(n_cases_cum = max(n_cases_cum)) %>% 
#     #   
#     #   # Sum up each i_group
#     #   summarise(n_cases_cum = sum(n_cases_cum)) %>% 
#     #   # print_all
#     #   
#     #   # Calculate median and confidence intervals
#     #   group_by(budget_contacts, budget_sar_out, target_type) %>% 
#     #   quantile_summarise(n_cases_cum, conf_level = 0.95) %>% 
#     #   arrange(budget_contacts, budget_sar_out)
#     
#     z_target <- target_crossed_summ %>% filter(target_type == "targeted") %>% 
#       .$infected %>% 
#       matrix(byrow = TRUE, nrow = length(unique(target_crossed_summ$budget_contacts)))
#     z_untarget <- target_crossed_summ %>% filter(target_type == "untargeted") %>% 
#       .$infected %>% 
#       matrix(byrow = TRUE, nrow = length(unique(target_crossed_summ$budget_contacts)))
#     
#     surface_plot <- plot_ly(
#       x = unique(target_crossed_summ$budget_contacts),
#       y = unique(target_crossed_summ$budget_sar_out)
#     ) %>%
#       add_surface(z = ~ z_target,
#                   colorscale = list(c(0, 1), c("rgb(300,184,214)", "rgb(300,90,124)")),
#                   contours = list(
#                     x = list(show = TRUE, highlight = TRUE), y= list(show = TRUE,  highlight = TRUE)
#                   )) %>%
#       add_surface(z = ~ z_untarget,
#                   contours = list(
#                     x = list(show = TRUE,  highlight = TRUE), y= list(show = TRUE,  highlight = TRUE)
#                   ))
#     
#     
#     surface_plot
#     
#     # PLOT (2) reduction in cases from targeting
#     target_crossed_diff <- target_crossed_summ %>% 
#       group_by(budget_sar_out, budget_contacts) %>% 
#       summarise(infected_diff = infected[target_type == "untargeted"] - infected[target_type == "targeted"])
#     
#     z_diff <- target_crossed_diff$infected_diff %>% 
#       matrix(byrow = TRUE, nrow = length(unique(target_crossed_diff$budget_contacts)))
#     
#     surface_plot_diff <- plot_ly(
#       x = unique(target_crossed_diff$budget_contacts),
#       y = unique(target_crossed_diff$budget_sar_out)
#     ) %>%
#       add_surface(z = ~ z_diff,
#                   # colorscale = list(c(0, 1), c("rgb(300,184,214)", "rgb(300,90,124)")),
#                   contours = list(
#                     x = list(show = TRUE, highlight = TRUE), y= list(show = TRUE,  highlight = TRUE)
#                   ))
#     
#     # CONCLUSION - very little additional benefit to doing targeting if you only do x (reduce contacts)
#     # but if you do SAR out interventions that's more useful, and if you do it's especially helpful to target
#     
#     
#     surface_plot_diff
#     
#     
#     
# 
#     
#     
#     
# # Random testing counterfactual -------------------------------------------
#     
#     
#     target_alpha <- tibble(
#       budget = c(0, 0.1, 0.2, 0.3) * 10000
#     ) %>% 
#       rowwise() %>% 
#       mutate(target_params = list(target_params(params_original = c(0, 0, 0, 0, 0),
#                                                 n_pops = n_pops,
#                                                 budget = budget, 
#                                                 target_groups = 1:2,
#                                                 direction = 1))) %>% 
#       unnest(target_params) %>% 
#       group_by(budget) %>% 
#       summarise(untargeted = list(untargeted),
#                 targeted = list(targeted)) %>% 
#       pivot_longer(-budget, names_to = "target_type", values_to = "alpha") %>% 
#       print
#     
#     
#     target_alpha %>% unnest(alpha) %>% print_all
#     
#     
#     # RUN THE MODELS
#     target_alpha_results <- target_alpha %>% 
#       mutate(model_output = pmap(select(., alpha), group_diff_default, n_sims = 2, n_iterations = 10000))
#     
#     
#     
#     # PLOT
#     target_alpha_results %>% 
#       mutate(model_output = map(model_output, "time_series")) %>%
#       unnest(model_output) %>% 
#       
#       group_by(budget, target_type, sim_id, t) %>% 
#       mutate(n_pop_total = sum(n_pop)) %>% 
#       mutate(budget = budget / n_pop_total) %>% 
#       ungroup %>% 
#       
#       # Calculate total infected in each run
#       group_by(budget, target_type, sim_id, i_group) %>% 
#       summarise(n_cases_cum = max(n_cases_cum)) %>% 
#       
#       # Sum up each i_group
#       summarise(n_cases_cum = sum(n_cases_cum)) %>% 
#       # print_all
#       
#       # Calculate median and confidence intervals
#       group_by(budget, target_type) %>% 
#       quantile_summarise(n_cases_cum, conf_level = 0.95) %>% 
#       
#       # PLOT
#       ggplot(aes(x = budget, y = n_cases_cum_median, colour = target_type)) + 
#       geom_line() + 
#       geom_errorbar(aes(ymin = n_cases_cum_lower, ymax = n_cases_cum_upper), width = 0.02, size = 0.5) + 
#       geom_point()
#     
#     
#     
#     
# # ........................ ------------------------------------------------
#     
#     
# # Sample selection --------------------------------------------------------
#     
#     
#     dat_covida <- read_dta("Datos_Salesforce_treated.dta")
#     
#     dat_covida %>% 
#       count_prop(contact, symptom) %>% 
#       count_prop(isolates) %>% 
#       count_prop(stratum) %>% 
#       count_prop(tested)
#     
#     # SEPARATE BY TESTED
#     dat_covida_symptoms <- dat_covida %>% filter(symptom == 1) %>% 
#       group_by(stratum, tested, symptom) %>% 
#       summarise(mean_se(isolates))
#     
#     dat_covida_contact <- dat_covida %>% 
#       filter(contact == 1) %>% 
#       group_by(stratum, tested, contact) %>% 
#       summarise(mean_se(isolates))
#     # mutate(contact_or_symptoms = contact | symptom) %>% 
#     
#     dat_covida_neither <- dat_covida %>% 
#       filter(symptom == 0 & contact == 0) %>% 
#       group_by(stratum, tested, symptom, contact) %>% 
#       summarise(mean_se(isolates))
#     
#     dat_by_tested <- bind_rows(dat_covida_symptoms, dat_covida_contact, dat_covida_neither) %>% 
#       mutate(
#         case_category = case_when(
#           symptom == 1 & is.na(contact) ~ "Symptoms",
#           is.na(symptom) & contact == 1 ~ "Contacted case",
#           symptom == 0 & contact == 0 ~ "No symptoms, no contact"
#         )
#       )
#     
#     
#     dat_by_tested %>% 
#       filter(stratum %in% 1:6) %>% 
#       ggplot(aes(x = factor(stratum), y = y, colour = factor(tested))) + 
#       geom_point(position = position_dodge(width = 0.2)) + 
#       theme(legend.position = "top") + 
#       geom_errorbar(aes(ymin = ymin, ymax = ymax), 
#                     width = 0.2,
#                     position = position_dodge(width = 0.2)) + 
#       facet_grid(case_category ~ ., scales = "free_y")
#     
#     
#     # ALL TOGETHER
#     dat_covida_symptoms <- dat_covida %>% filter(symptom == 1) %>% 
#       group_by(stratum, symptom) %>% 
#       summarise(mean_se(isolates))
#     
#     dat_covida_contact <- dat_covida %>% 
#       filter(contact == 1) %>% 
#       group_by(stratum, contact) %>% 
#       summarise(mean_se(isolates))
#     # mutate(contact_or_symptoms = contact | symptom) %>% 
#     
#     dat_covida_neither <- dat_covida %>% 
#       filter(symptom == 0 & contact == 0) %>% 
#       group_by(stratum, symptom, contact) %>% 
#       summarise(mean_se(isolates))
#     
#     
#     dat_covida %>% 
#       mutate(symptom_or_contact = symptom | contact) %>% 
#       group_by(tested, symptom_or_contact, stratum) %>% 
#       summarise(mean_se(isolates)) %>% 
#       filter(stratum %in% 1:6) %>% 
#       
#       ggplot(aes(x = factor(stratum), y = y, colour = factor(tested))) + 
#       geom_point(position = position_dodge(width = 0.2)) + 
#       # geom_line(aes(y = y), alpha = 0.4) +
#       theme(legend.position = "top") + 
#       scale_colour_discrete(labels = c("Survey no test", "Test")) + 
#       geom_errorbar(aes(ymin = ymin, ymax = ymax), width = 0.2, position = position_dodge(width = 0.2)) + 
#       geom_line(aes(x = stratum, y = y, colour = factor(tested)), alpha = 0.5) + 
#       facet_wrap(~ symptom_or_contact, nrow = 2,
#                  labeller = as_labeller(
#                    c("TRUE" = "Symptoms or contact", "FALSE" = "No symptoms, no contact")
#                  )) + 
#       labs(x = "SES Group", y = "P(Stayed at home over last 14 days)", colour = element_blank()) + 
#       theme_custom()
#     
#     ggsave("figures/sample_selection_isolates.png", width = 7, height = 5.5, scale = 0.7)
#     
#     
#     
#     dat_covida$isolates
#     
