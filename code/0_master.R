

    dir.create("data/processed")
    dir.create("data/temp")
    dir.create("figures")


# SET UP and IMPORT ------------------------------------------------------------------

    # 1 Run all functions
    source("code/individual_model_functions.R")
    

    # 2 Import all data
    source("code/import_data.R")

    # 3 Calibrate timing variables (serial interval, secondary timing etc.)
    source("code/calibrate_timings.R")

    # 3.5 Calculate initial R0 in the data
    source("code/r0_calc.R") # slope = 0.038; R0 used to be 1.242 [with gen interval = 6]; now R0 = 1.218

    # 4 choose initial value for outside of work factor, and match to R0 slope
    n_cores <- 25; n_sims <- 50; n_iterations <- 80
    source("code/r0_calibrate_systematic.R")

    # 5 Run the contact matrix with the selected outside_work_factor from (4) to save final contact_matrix
    # outside_work_factor <- 0.95    # Note: this scales outside of work contacts linearly
    source("code/calc_contact_matrix.R")
    
    
# OPTIONAL - TEST the model ----------------------------------------------------------
    
    test_model <- FALSE
    
    if (test_model) {
        
        # Clear the environment
        rm(list = ls())
        
        # 1 Import all functions and data from first step
        source("code/individual_model_functions.R")
        load("data/processed/data_save.RData")
        
        # 2 Run unit tests
        library("testthat")
        source("code/unit_tests_prep.R") # generates the data from a simulation to be used in unit_tests.R
        testthat::test_file("code/unit_tests.R")
        
        # 3 Ad-hoc tests can be run with this file:
        source("code/individual_model_testing.R")
        
    }
    


# RUN THE MODEL / generate results -----------------------------------------------------------

    # (1) Main models for results

        # Prep 
        rm(list = ls())
        source("code/individual_model_functions.R"); load("data/processed/data_save.RData")
        source("code/functions_group_diff.R"); source("code/params_baseline.R")

        # Models
        n_cores <- 25; n_sims <- 50; n_iterations <- 1000000
        # n_cores <- 4; n_sims <- 16; n_iterations <- 1000000
        # STARTED AT 13:00
        # source("code/run_model_counterfactuals.R")
        source("code/run_model_counterfactuals_2.R")

        # Graphs
        # source("code/graphs_counterfactuals.R")    # and output graphs from this
        source("code/graphs_counterfactuals_2.R")
    
    # (2) "Policy" style results
        
        # Prep 
        rm(list = ls())
        source("code/individual_model_functions.R"); load("data/processed/data_save.RData")
        source("code/functions_group_diff.R"); source("code/params_baseline.R")
        
        # Models
        n_cores <- 25; n_sims <- 50; n_iterations <- 1000000
        # n_cores <- 1; n_sims <- 1; n_iterations <- 1
        source("code/run_model_inequality_is_bad.R")
        
        rm(list = ls()); source("code/individual_model_functions.R"); load("data/processed/data_save.RData"); source("code/functions_group_diff.R"); source("code/params_baseline.R")
        n_cores <- 25; n_sims <- 50; n_iterations <- 1000000
        # n_cores <- 1; n_sims <- 1; n_iterations <- 1
        source("code/run_model_target.R")
        source("code/run_model_0_testing.R")
    
        # Graphs
        conf_level <- 0.95
        source("code/graphs_inequality_is_bad.R")
        source("code/graphs_target.R")
        source("code/graphs_policy_scenarios.R")
    
    
    # (3) Graphs for parameters
    
        # Prep 
        rm(list = ls())
        source("code/individual_model_functions.R"); load("data/processed/data_save.RData")
        source("code/functions_group_diff.R"); source("code/params_baseline.R")
    
        source("code/graphs_params.R")
    
    
    
    
    
    
    
# CALIBRATE MOBILITY ------------------------------------------------------
    
    
    # Prep 
    rm(list = ls())
    source("code/individual_model_functions.R"); load("data/processed/data_save.RData")
    source("code/functions_group_diff.R"); source("code/params_baseline.R")
    
    # 5 import case data (for mobility matching)
    source("code/import_case_data.R")
    
    # 6 Calibrate mobility
    n_cores <- 25; n_sims_per_round <- 100; n_loops <- 20; t_max <- 500
    # n_cores <- 1; n_sims_per_round <- 1; n_loops <- 1; t_max <- 500
    dir.create("figures/tweaking_process")
    source("code/calibrate_mobility_final.R")
    
    # 6.5 run a final sim, and save mobility to data_save
    # Choose final iteration
    # n_cores <- 1; n_sims_per_round <- 1; n_iterations <- 100000
    n_cores <- 25; n_sims_per_round <- 50; n_iterations <- 100000
    tweak_to_choose <- n_loops
    source("code/final_sim_check.R")
    
    # PLOT MATCH with baseline case (no mobility) and with mobility change
    # Still dependent on import_case_data.R
    conf_level <- 0.95 
    source("code/import_covida_data.R")
    source("code/graphs_model_match.R")
    
    
    # PLOT counterfactuals with mobility
    rm(list = ls()); source("code/individual_model_functions.R"); load("data/processed/data_save.RData"); source("code/functions_group_diff.R"); source("code/params_baseline.R")
    # n_cores <- 1; n_sims <- 1; n_iterations <- 1
    n_cores <- 25; n_sims_per_round <- 50; n_iterations <- 100000
    source("code/run_model_counterfactuals_mobility_change.R")
    source("code/graphs_counterfactuals_mobility_change.R")
    
    
    
    # OPTIONAL
    # 7 Check the new mobility (e.g. plots etc.)
    # source("code/calibrate_mobility_check.R")
    
    

