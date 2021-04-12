# # BASELINE PARAMETERS
params_baseline <- list(
  group_props = data_save$group_props,
  n_pop = 100000,
  
  ct_delay_data             = data_save$ct_delay_data,
  probs_self_test_df        = data_save$probs_self_test,
  probs_isolate_symptoms_df = data_save$probs_isolate_symptoms,
  probs_isolate_ct_df       = data_save$probs_isolate_ct,
  probs_isolate_test_df     = data_save$probs_isolate_test,
  days_of_work              = data_save$days_of_work,
  days_work_ct              = data_save$days_work_ct,
  days_work_hh_quarantine   = data_save$days_work_hh_quarantine,
  days_work_symptoms        = data_save$days_work_symptoms,
  p_contact_if_isolated_home= rep(1, 4),
  p_hh_quarantine           = data_save$p_hh_quarantine,
  p_contact_traced          = data_save$p_contact_traced,
  hh_size_data              = data_save$hh_data_bogota,
  test_delay_data           = data_save$test_delay_data,
  k_matrix                  = data_save$k_matrix,
  sar_out                   = data_save$sar_out_input,
  sar_home                  = data_save$sar_home_input,
  alpha                     = rep(0, 4)
)