RNGkind("L'Ecuyer-CMRG"); set.seed(12345)
model_0_testing <- group_diff(
  n_cores = n_cores, n_sims = n_sims, print_detail = FALSE, n_iterations = n_iterations,
  probs_self_test_df = params_baseline$probs_self_test_df %>% mutate(prob = 0),    # 0 prob of self test
  p_contact_traced = rep(0, 4)
)

save(model_0_testing, file = "data/processed/0_testing_slim.RData")