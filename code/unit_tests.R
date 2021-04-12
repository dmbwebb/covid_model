

    library(testthat)
    library(magrittr) # so I can use %$%


# 1. Live cases --------------------------------------------------------------



test_that("no dups", {
  expect_true(
    live_cases_list %>% dups_report(case_id, t, print_only = FALSE) %>% nrow() == 1
  )
})

# (4) Check we only have isolate_symptoms_timing for people who actually isolate after symptoms
test_that("isolate_symptoms_timing", {
  n_bugs <- live_cases_list %$% sum((isolate_after_symptoms & !value(isolate_symptoms_timing)) | (!isolate_after_symptoms & value(isolate_symptoms_timing)))
  expect_equal(n_bugs, 0)
})


# Check that contact tracing times never stay the same for a given person
# Former bug - contact_tracing_timing and contact_tracing_test_timing don't update correctly over periods [should increment by 1 each time]
test_that("contact_tracing_timing and contact_tracing_test_timing", {
  
  n_safe <- live_cases_list %>% 
    group_by(case_id) %>% 
    arrange(t) %>% 
    mutate(ct_test_timing_bug = contact_tracing_test_timing == lag(contact_tracing_test_timing),
           ct_results_timing_bug = contact_tracing_results_timing == lag(contact_tracing_results_timing)) %>% 
    ungroup %$% 
    sum((ct_test_timing_bug == FALSE & ct_results_timing_bug == FALSE) | 
          is.na(ct_test_timing_bug) & is.na(ct_results_timing_bug))  # should always be false / NA
  
  expect_equal(n_safe, nrow(live_cases_list))
  
})

# (6) Check that symptom severity never changes for an individual
test_that("symptom severity doesn't change", {
  out <- live_cases_list %>% group_by(case_id) %>% summarise(n_distinct = n_distinct(symptom_severity)) %>% 
    count_prop(n_distinct, return_count = TRUE) # should only have n_distinct = 1
  
  expect_equal(out$n_distinct, 1L)
})



# (7) Check i_group never changes
test_that("i_group never changes", {
  out <- live_cases_list %>% group_by(case_id) %>% summarise(n_distinct = n_distinct(i_group)) %>% 
    count_prop(n_distinct, return_count = TRUE) # should only have n_distinct = 1
  
  expect_equal(out$n_distinct, 1L)
})

# No one is previously tested having NEVER been currently tested before?
# YES - people come up as previously tested when they were contact tested BEFORE they were infected...
test_that("No one is previously tested having NEVER been currently tested before?",
          {
            out <- live_cases_list %>% 
              # select(case_id, t, currently_testing, previously_tested) %>% 
              group_by(case_id) %>% 
              mutate(across(c(currently_testing, previously_tested), as.numeric)) %>% 
              mutate(ever_tested_debug = cummax(currently_testing)) %>% 
              # mutate(testing_offswitch = currently_testing == 0 & lag(currently_testing == 1)) %>% 
              mutate(testing_bug = previously_tested == TRUE & ever_tested_debug == FALSE) %>% 
              # ungroup %>% 
              # count_prop(testing_bug) %>% 
              mutate(traced_before_infection = sum(infection_timing > contact_tracing_test_timing, na.rm = TRUE) > 0,
                     id_has_testing_bug = sum(testing_bug) > 0) %>% 
              ungroup %$%
              sum(
                traced_before_infection == FALSE & id_has_testing_bug == TRUE
              )
            # there should be no people with traced_before_infection = FALSE and id_has_testing_bug = TRUE
            
            expect_equal(out, 0)
          })


# (9) [Previous bug] random_test_yn_ever should never be NA
test_that("random_test_yn_ever should never be NA", {
  out <- live_cases_list %$% sum(is.na(random_test_yn_ever) | is.na(random_test_yn))
  expect_equal(out, 0)
})


# (10) Check isolation behaviour

test_that("No cases that have never been randomly tested and have an isolation time due to random testing", {
  expect_equal(0, live_cases_list %$% sum(!random_test_yn_ever & value(isolate_random_test_timing))) # should be 0
})

test_that("No cases that isolate due to random testing and have isolate_after_test is F", {
  expect_equal(0, live_cases_list %$% sum(!isolate_after_test_presymp & !isolate_after_test_symp & value(isolate_random_test_timing))) # should be 0
})

test_that("No cases for whom test is finished, and random_test_yn is still TRUE", {
  # This was not the case when random testing was still done on people who had been tested negative
  expect_equal(0, live_cases_list %$% sum(random_test_timing < 0 & random_test_yn)) # should be 0
})

test_that("should have no people who are symptom_timing != isolate_symptoms_timing", {
  expect_equal(0, live_cases_list %$% sum(symptom_timing != isolate_symptoms_timing, na.rm = TRUE))
})



test_that("No symptoms before infections", {
  expect_equal(0, live_cases_list %$% sum(symptom_timing < infection_timing)) # should have none < 0
})




# 2. Secondary cases ------------------------------------------------------

test_that("Number of people who have an ID but don't have a secondary_case_timing - should be 0", {
  expect_equal(0, secondary_cases_list %$% sum(!is.na(secondary_case_id) & is.na(secondary_case_timing))) # should be 0
})

test_that("Secondary case timing should always be greater than infection timing", {
  expect_equal(0, 
               secondary_cases_list %$% sum(secondary_case_timing - infection_timing < 0, na.rm = TRUE)) # should be 0
})

test_that("Symptom timing always after infection", {
  expect_equal(0, 
               secondary_cases_list %$%
                 sum(symptom_timing - infection_timing < 0, na.rm = TRUE))
})


# (3) make sure isolation happens only if isolation_timing is less than secondary_case_timing
# secondary_isolated == TRUE only if isolation_timing < secondary_case_timing == TRUE
# secondary_cases_list %>% 
  # count_prop(transmission_isolated, isolation_timing < secondary_case_timing, is.na(isolation_timing) | is.na(secondary_case_timing)) # should always be TRUE at the same time
# OBSOLETE - you can now be isolated because of household level quarantining **

test_that("Check relationship between secondary_isolated and contact_if_isolated", {
  # count <- secondary_cases_list %>% count_prop(transmission_isolated, contact_if_isolated, return_count = TRUE) # contact if isolated should be NA if secondary_isolated is FALSE, and T/F when secondary_isolated is TRUE

  
  out <- secondary_cases_list %$%
    sum(
      (transmission_isolated == FALSE & is.na(contact_if_isolated)) |
        (transmission_isolated == TRUE & (contact_if_isolated == TRUE | contact_if_isolated == FALSE))
    )
    
  
  # out <- c(count$transmission_isolated, count$contact_if_isolated)
  
  expect_equal(out, nrow(secondary_cases_list))

})

test_that("Make sure only people who don't have a secondary_case_id [i.e. non cases] are NA for secondary_isolated", {
  expect_equal(0, secondary_cases_list %$%
                 sum(!is.na(secondary_case_id) & is.na(transmission_isolated)))
}) 


test_that("Prop susceptible never over 1 or less than 0", {
  expect_equal(0, secondary_cases_list %>% filter(is.na(transmission_isolated)) %>% nrow()) # should be 0
  
})

test_that("Make sure (some) people isolate at contact_tracing_test_timing if they ahve isolate_after_ct", {
  out <- secondary_cases_list %$% sum((isolate_after_ct_presymp | isolate_after_ct_symp) & isolation_timing == contact_tracing_test_timing, na.rm = TRUE) # should be >0
  
  expect_true(
    out > 0
  )
  
})

test_that("Make sure deisolate_timing is not only NA", {
  out <- secondary_cases_list %$% sum(value(deisolate_timing)) # should be greater than 0

  expect_true(out > 0)
})

test_that("Make sure isolation_timing_2 only exists if deisolate_timing exists", {
  out_1 <- secondary_cases_list %$% sum(value(isolation_timing_2) == TRUE & value(deisolate_timing) == FALSE) # should be 0
  out_2 <- secondary_cases_list %$% sum(value(isolation_timing_2) == TRUE & value(deisolate_timing) == TRUE & deisolate_timing >= isolation_timing_2) # (deisolate timing is always less than isolaton_timing_2) should be 0

  expect_true(out_1 == 0 & out_2 == 0)
})

test_that("Make sure secondary case is isolated if it's between isolation and deisolate, or after isolation_timing_2", {
  
  out_1 <- secondary_cases_list %$% sum(transmission_isolated == FALSE & infection_type == "out" & secondary_case_timing >= isolation_timing  & secondary_case_timing < deisolate_timing, na.rm = TRUE)  #should be 0
  out_2 <- secondary_cases_list %$% sum(transmission_isolated == FALSE & infection_type == "out" & secondary_case_timing >= isolation_timing_2, na.rm = TRUE)  #should be 0

  expect_true(out_1 == 0 & out_2 == 0)
  # secondary_cases_list %$% sum(transmission_isolated == TRUE & secondary_case_timing < isolation_timing, na.rm = TRUE)  #should be 0 - OBSOLETE because of hh quarantining
  # secondary_cases_list %$% sum(transmission_isolated == TRUE & secondary_case_timing > deisolate_timing & (!value(isolation_timing_2) | secondary_case_timing < isolation_timing_2), na.rm = TRUE)  #should be 0 - OBSOLETE because of hh quarantining
  
})

test_that("Make sure deisolate timing only exists when isolate_ct_test_timing exists", {
  out <- secondary_cases_list %$% sum(value(deisolate_timing) & value(isolate_ct_test_timing) == FALSE) # should be 0
  expect_equal(0, out)
})

# No transmission isolated at home [OBSOLETE, change of setting]
# secondary_cases_list %$% sum(transmission_isolated == TRUE & hh_id == secondary_hh_id, na.rm = TRUE)  #should be 0


test_that("isolate_ct_test_timing when people are isolated and ct tested", {
  # (FORMER BUG - solved by improving update_isolation to make sure values for isolation are kept up to date)
  out <- secondary_cases_list %$% 
    # count_prop(isolate_after_ct_symp, value(contact_tracing_test_timing), value(isolate_ct_test_timing)) %$%
    sum(isolate_after_ct_symp & value(contact_tracing_test_timing) & !value(isolate_ct_test_timing)) # should be 0
  expect_equal(0, out)
}) 

test_that("Only isolate if you are tested", {
  out <- sum(
    is.na(secondary_cases_list$contact_tracing_test_timing) &
      !is.na(secondary_cases_list$isolate_ct_test_timing)
  ) # should be 0
  
  expect_equal(0, out)
})

test_that("only isolate if you have results", {
  out <- sum(
    is.na(secondary_cases_list$contact_tracing_results_timing) &
      !is.na(secondary_cases_list$isolate_ct_results_timing)
  ) # should be 0
  
  expect_equal(0, out)
})


test_that("people don't isolate if they have a false negative", {
  out_1 <- secondary_cases_list %$% sum(ct_false_negative & value(isolate_ct_results_timing)) # should be 0
  
  out_2 <- secondary_cases_list %$% sum(random_test_false_negative & value(isolate_random_test_timing)) # should be 0
  
  out_3 <- secondary_cases_list %$% sum(self_test_false_negative & value(isolate_self_test_timing)) # should be 0
  
  
  # live_cases_list %$% sum(ct_test_negative & value(isolate_ct_test_timing)) 
  # actually this can be above 0 - people isolate just from CT not based on results
  
  expect_true(out_1 + out_2 + out_3 == 0)
})


test_that("Previous bug - involving people being randomly tested twice", {
  # Previous bug - involving people being randomly tested twice
  out <- live_cases_list %>%  
    group_by(case_id) %>% 
    select(case_id, t, contains("random_test"), currently_testing, previously_tested) %>% 
    filter(sum(random_test_yn & random_test_timing < 0) > 0)
  
  n_bugs <- nrow(out)
  # BUG : some people for whom random test is finished (random_test_timing < 0) don't have random_test_yn
  # this was people who are randomly tested in the same period 
  # SOLVED by changing previously_tested_positive condition to previously_tested
  
  expect_equal(nrow(out), 0)
})


test_that("people can never isolate presymptoms for self test ", {
  out <- secondary_cases_list %$% sum(value(isolate_self_test_timing) & isolate_self_test_timing < symptom_timing) # NEVER pre-symptoms - should be 0
  expect_equal(0, out)
})

test_that("secondary case group is always the same for every row secondary_case_id", {
  out <- secondary_cases_list %>% select(secondary_case_id, secondary_i_group) %>% 
    dups_drop(warn = FALSE) %>% 
    dups_count(secondary_case_id) # should be no duplicates
  
  expect_equal(0, out)
})

test_that("symptom timing matches secondary symptom timing", {
  
  primary <- secondary_cases_list %>% 
    mutate(symptoms_t = symptom_timing + t) %>%
    select(case_id, symptoms_t) %>% 
    dups_drop(warn = FALSE)
  
  secondary <- secondary_cases_list %>%
    mutate(
      secondary_symptoms_t = secondary_symptoms_timing + t
    ) %>%
    select(secondary_case_id, secondary_symptoms_t) %>%
    dups_drop(warn = FALSE)
  
  matched <- inner_join(primary, secondary, by = c("case_id" = "secondary_case_id"))
  
  # matched %>% filter(symptoms_t - secondary_symptoms_t > 0.001) %>% 
  #   print_all

  n_bugs <- sum(abs(matched$symptoms_t - matched$secondary_symptoms_t) > 0.001)
  
  expect_equal(n_bugs, 0)
  
})


test_that("secondary case timing matches infection timing", {
  
  primary <- secondary_cases_list %>% 
    mutate(infection_t = infection_timing + t) %>%
    select(case_id, infection_t) %>% 
    dups_drop(warn = FALSE)
  
  secondary <- secondary_cases_list %>%
    mutate(
      secondary_case_t = secondary_case_timing + t
    ) %>%
    select(secondary_case_id, secondary_case_t) %>%
    dups_drop(warn = FALSE)
  
  matched <- inner_join(primary, secondary, by = c("case_id" = "secondary_case_id"))
  
  n_bugs <- sum(abs(matched$secondary_case_t - matched$infection_t) > 0.001)
  
  expect_true(n_bugs == 0)
  
})


test_that("primary symptoms timing = symptom_timing", {
  n_bugs <- secondary_cases_list %$%
    sum(abs(symptom_timing - primary_symptoms_timing) > 0.001)
  
  # secondary_cases_list %>% 
  #   select(case_id, t, symptom_timing, primary_symptoms_timing) %>%
  #   filter(symptom_timing - primary_symptoms_timing > 0.001) %>% 
  #   print(n = 1000)
  
  expect_true(n_bugs == 0)
})



live_ids <- live_cases_list %>% pull(case_id)

# test_that("doesn't become live if contact_if_isolated is false", {
#   
#   isolated_ids <- secondary_cases_list %>% 
#     group_by(case_id, secondary_case_id) %>% 
#     filter(sum(contact_if_isolated == FALSE) > 0 & sum(contact_if_isolated == TRUE) > 0) %>% 
#     select(case_id, secondary_case_id, t, secondary_case_timing, 
#            transmission_isolated,
#            # matches("isolate_interval"),
#            contact_if_isolated, random_test_yn) %>% 
#     print(n = 100)
#   
#   isolated_ids <- secondary_cases_list %>% 
#     group_by(case_id, secondary_case_id) %>% 
#     filter(last(contact_if_isolated) == FALSE) %>% 
#     select(case_id, secondary_case_id, t, secondary_case_timing, 
#            transmission_isolated,
#            # matches("isolate_interval"),
#            contact_if_isolated, random_test_yn) %>% 
#     print(n = 100) %>% 
#   
#     
#     pull(secondary_case_id)
#   
#   n_bugs <- sum(isolated_ids %in% live_ids)
#   
#   expect_equal(n_bugs, 0)
#   
# })


test_that("secondary hh details match new hh details", {
  
  primary <- secondary_cases_list %>% 
    select(case_id, hh_ind_id, hh_id, hh_size, would_quarantine)
  
  secondary <- secondary_cases_list %>% 
    select(secondary_case_id, secondary_hh_ind_id, secondary_hh_id, secondary_hh_size, secondary_would_quarantine)
  
  matched <- inner_join(primary, secondary,
                        by = c("case_id" = "secondary_case_id")) %>% 
    mutate(bug_1 = hh_ind_id != secondary_hh_ind_id, 
           bug_2 = hh_id != secondary_hh_id,
           bug_3 = hh_size != secondary_hh_size,
           bug_4 = would_quarantine != secondary_would_quarantine)
  
  n_bugs <- matched %$% sum(bug_1 | bug_2 | bug_3 | bug_4)
  
  expect_equal(n_bugs, 0)
  
})


test_that("n_contacts_out > n_secondary_cases_out", {
  n_bugs <- secondary_cases_list %>% 
    mutate(bug = n_contacts_out < n_secondary_cases_out) %$%
    sum(bug)
  
  expect_equal(n_bugs, 0)
})

test_that("primary_symptoms_timing + serial_interval = secondary_symptoms_timing", {
  n_bugs <- secondary_cases_list %>% 
    mutate(bug = abs(primary_symptoms_timing + serial_interval - secondary_symptoms_timing) > 0.001) %$%
    sum(bug)
  
  expect_equal(n_bugs, 0)
})


test_that("deisolate_timing should not be less than isolation_timing", {
  n_bugs <- secondary_cases_list %>% 
    mutate(bug = deisolate_timing < isolation_timing) %$%
    sum(bug, na.rm = TRUE)
  
  expect_equal(n_bugs, 0)
})

# secondary_cases_list %>% count_prop(round(isolation_timing - deisolate_timing))



# test_that("")


# secondary_cases_list %>% print_names



# 3. Live / secondary match checks ----------------------------------------

cases_merged <- trackr::full_join_track(live_cases_list, secondary_cases_list,
                                        by = c("case_id", "t"),
                                        suffix = c("_L", "_S"),
                                        .merge = TRUE)

# (1) Are there people in secondary cases who aren't in live_cases? - there shouldn't be
test_that("Are there people in secondary cases who aren't in live_cases? - there shouldn't be", {
   # should be none in the y_only section
  # Previous bug - some people were in secondary cases but not live cases when they RECOVER before the secondary case is meant to happen
  expect_equal(0, sum(cases_merged$.merge == "y_only"))
})


# (2)
test_that("no dups when matched", {
  expect_equal(0, cases_merged %>% dups_count(t, case_id, secondary_case_id))
})

# test_that("no conflicting variables", {
#   # Check there are no conflicts between values across datasets
# 
#   
#     out <- variable_match %>% ungroup %>% summarise(across(ends_with("_conflict"), list(total = sum, prop = mean))) %>% 
#       pivot_longer(everything(), names_to = c("variable", ".value"), names_pattern = "(.*)_(total|prop)") %$%
#       sum(total)
#     
#     expect_equal(out, 0)   # there should be no conflicting variables
#     
#   
#   
# })


# test_that("no contact tracing timings before infection timing", {
#   # (5) Why is contact_tracing stuff often added when negative? 
#   # i.e. the "first" time we see a contact tracing timing is when it's negative?
#   # - SOLVED - edited update_live_contact_tracing so that contact tracing timings are removed (set to NA) when they occur before infection time
# 
#     out <- variable_match %$% sum(contact_tracing_test_timing < infection_timing | contact_tracing_results_timing < infection_timing, na.rm = TRUE)
#     expect_equal(out, 0)
# 
# })





# 4. HH quarantining checks -----------------------------------------------

# (1) 
test_that("if in the same household, primary_quarantine should always be the same as secondary_quarnatined", {
  out <- secondary_cases_list %$% 
    sum(hh_id == secondary_hh_id & primary_quarantined != secondary_quarantined, na.rm = TRUE) # should be 0
  expect_equal(0, out)
})

test_that("All NAs for new_potential_case, susceptible, new_actual_case, pipped should come from NEW ENTRIES", {
  nas <- secondary_cases_list %>% select(new_potential_case, susceptible, new_actual_case, pipped) %>% count_nas(return_count = TRUE)
  new_entries <- secondary_cases_list %>% group_by(case_id, secondary_case_id) %>% filter(t == min(t))
  
  n_new_entries <- nrow(new_entries)
  # new_entries %>% select(new_potential_case, susceptible, new_actual_case, pipped) %>% count_nas(return_count = TRUE)
  
  expect_equal(sum(nas$missing == n_new_entries), 4)
})




primary_quarantine_check <- function(i) {
  
  # print(paste0("rep ", i))
  
  # Randomly select a case_id / secondary_id combination that's primary quarantined
  random_primary_quarantined <- secondary_cases_list %>% filter(primary_quarantined) %>% 
    mutate(secondary_case_t = t + secondary_case_timing, .after = secondary_case_id) %>% 
    select(case_id, secondary_case_id, secondary_case_t, primary_quarantined, hh_id) %>% 
    dups_drop(warn = FALSE) %>% 
    slice_sample(n = 1) #%>% 
  # {print(.$secondary_case_t); invisible(.)}
  
  # Look at the corresponding hh status data to check interval
  random_primary_hh <- hh_status_list %>% 
    semi_join(random_primary_quarantined, by = "hh_id") %>% 
    mutate(across(matches("isolate_interval"),
                  ~ t + .x)) %>% 
    select(matches("isolate_interval")) %>% 
    dups_drop(warn = FALSE) %>% 
    janitor::remove_empty(which = c("rows"))
  
  
  # Check whether the value's in the interval
  check_tf <- any(data.table::between(random_primary_quarantined$secondary_case_t, random_primary_hh$isolate_interval_1_l, random_primary_hh$isolate_interval_1_r, NAbounds = NA)) | 
    any(data.table::between(random_primary_quarantined$secondary_case_t, random_primary_hh$isolate_interval_2_l, random_primary_hh$isolate_interval_2_r, NAbounds = NA))
  
  if (is.na(check_tf)) browser()
  
  return(check_tf)
  
}

# Check 100 times
test_that("primary quarantine", {
  print("primary quarantine check (100 reps)")
  out <- tibble(i = 1:100) %>% 
    rowwise() %>% 
    mutate(primary_quarantine_check = primary_quarantine_check(i)) %$%
    sum(primary_quarantine_check == TRUE)
  
  expect_equal(out, 100)
})


# Secondary
# print("secondary quarantine check (100 reps)")
secondary_quarantine_check <- function(i) {
  
  # print(paste0("rep ", i))
  
  # Randomly select a case_id / secondary_id combination that's secondary quarantined
  random_secondary_quarantined <- secondary_cases_list %>% filter(secondary_quarantined) %>% 
    mutate(secondary_case_t = t + secondary_case_timing, .after = secondary_case_id) %>% 
    select(case_id, secondary_case_id, secondary_case_t, secondary_quarantined, secondary_hh_id) %>% 
    dups_drop(warn = FALSE) %>% 
    slice_sample(n = 1) #%>% 
  # {print(.$secondary_case_t); invisible(.)}
  
  # Look at the corresponding hh status data to check interval
  random_secondary_hh <- hh_status_list %>% 
    semi_join(random_secondary_quarantined, by = c("hh_id" = "secondary_hh_id")) %>% 
    mutate(across(matches("isolate_interval"),
                  ~ t + .x)) %>% 
    select(matches("isolate_interval")) %>% 
    dups_drop(warn = FALSE) %>% 
    janitor::remove_empty(which = c("rows"))
  
  
  # Check whether the value's in at least one of the intervals
  check_tf <- any(data.table::between(random_secondary_quarantined$secondary_case_t, random_secondary_hh$isolate_interval_1_l, random_secondary_hh$isolate_interval_1_r, NAbounds = NA)) | 
    any(data.table::between(random_secondary_quarantined$secondary_case_t, random_secondary_hh$isolate_interval_2_l, random_secondary_hh$isolate_interval_2_r, NAbounds = NA))
  
  
  
  return(check_tf)
  
}

# Check 100 times
test_that("primary quarantine", {
  print("secondary quarantine check (100 reps)")
  
  out <- tibble(i = 1:100) %>% 
    rowwise() %>% 
    mutate(secondary_quarantine_check = secondary_quarantine_check(i)) %$%
    sum(secondary_quarantine_check == TRUE)
  
  expect_equal(out, 100)
})









