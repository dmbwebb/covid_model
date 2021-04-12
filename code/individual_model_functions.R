# Set up ------------------------------------------------------------------
  
# library(devtools)
# devtools::install_github("dmbwebb/dups", dep = FALSE)
# library("dups")
# ?view()

    library(tidyverse)
    library("svglite")
    library(magrittr) # for using %$% in testing
    library(zoo)
    library(devtools)
    library(mc2d)
    library(tictoc)
    library(sn)
    library(viridis) # for nice plots
    library("lubridate")
    library("parallel")
    library("pbmcapply")
    library("ggrepel")
    # source("~/Desktop/functions.R")
           
    # devtools::install_github("dmbwebb/dups", dep = FALSE)
    if (sum(str_detect(search(), "package:dups")) > 0)  detach("package:dups", unload=TRUE) # detach if already loaded
    library(dups)
    
    # devtools::install_github("dmbwebb/trackr", dep = FALSE)
    if (sum(str_detect(search(), "package:trackr")) > 0)  detach("package:trackr", unload=TRUE) # detach if already loaded
    library(trackr)
    


# Convenience functions ---------------------------------------------------

    # Theme for graphs
    theme_custom <- function(...) {
      theme_light() + 
        theme(plot.background = element_blank(),
              strip.background=element_rect(fill="#D1E0E6"), 
              strip.text = element_text(color = "black"),
              panel.grid.minor = element_blank(),
              ...)
    }
    
    # Quick histogram
    hist_basic <- function(dat, x, binwidth = NULL, boundary = 0) {
      plot <- ggplot(dat, aes(x = {{x}})) + 
        geom_histogram(boundary = boundary, binwidth = binwidth, colour = "grey", fill = "lightblue")
      
      print(plot)
      
      invisible(dat)
    }
    
    # Converts all NAs to false
    na_to_false <- function(x) {
      if (!is_logical(x)) stop("Not a logical vector")
      if_else(is.na(x), FALSE, x)
    }
    
    # Tests equality of two vectors, but deals with NAs as if they were values
    harsh_equal <- function(x, y) {
      st_equal <- x == y
      both_nas <- is.na(x) & is.na(y)
      dplyr::if_else(is.na(st_equal), both_nas, st_equal)
    }
    
    # Vectorised Bernoulli sample
    bern_vec <- function(n, p, NA_error = TRUE) {
      suppressWarnings(out <- as.logical(mc2d::rbern(n, p = p)))
      if (sum(is.na(out)) > 0 && NA_error) {
        print(n)
        print(p)
        stop("Bernoulli calc is causing NAs")
      }
      out
    }
    
    value <- function(x) {
      !is.na(x)
    }
    
    # Function to quickly see all columns in a df
    print_names <- function(x) {
      x %>% names %>% enframe %>% print_all
      invisible(x)
    }
    

    # Calculate delays
    interval_days <- function(x, y) {interval(x, y) %/% days(1)}
    
    
    recode_i_group <- function(x) {
      fct_recode(factor(x), "1&2" = "1", "3" = "2", "4" = "3", "5&6" = "4")
    }

# Distributions for individual cases --------------------------------------
    
    
    # Symptom severity randomly drawn from between 1 and 5 [TODO: add probability of each severity]
    symptom_severity <- function(i_group) {
      n <- length(i_group)
      sample(1:2, size = n, replace = TRUE, prob = c(0.2, 0.8)) # 
      # sample(1:5, size = n, replace = TRUE)
    }
    
    
    # tibble(
    #   symptom_severity = symptom_severity(sample(c(1,2,3), 100, replace = TRUE))
    # ) %>%
    #   ggplot(aes(x = symptom_severity)) + geom_bar()
    
    
    # Draws symptom timing "index" from a skewed normal 
    # uses skewed normal so it can be drawn in tandem with the secondary case timing as a multivariate skewed normal
    symptom_timing_index <- function(n, marg_dist_primary) {
      sn::rsn(n = n, dp = marg_dist_primary@dp)
    }
    
    # Converts symptom timing "index" into the actual symptom timing
    symptom_timing <- function(symptom_timing_index, marg_dist_primary, params_symptom_timing) {
      primary_q <- sn::psn(x = - symptom_timing_index, dp = marg_dist_primary@dp) # NOTE THE MINUS SIGN! (Used in secondary_case_timing to generate negative correlation between serial and primary symptoms)
      primary_converted <- qlnorm(primary_q, meanlog = params_symptom_timing$meanlog, sdlog = params_symptom_timing$sdlog)
      return(primary_converted)
    }
    

    # Calculates the parameter of the conditional SN distributions, 
    # when given the first parameter (VECTORISED)
    # Based on the function code for sn::conditionalSECdistr but vectorised so it's a LOT faster
    conditional_sn_vec <- function(object, fixed.comp, fixed.values) {
      
      # FOR TESTING:
      # object <- joint_dist
      # fixed.comp <- 1
      # fixed.values <- c(-2, -1, 0, 1, 2)
      
      # Extract parameters:
      dp <- slot(object, "dp")
      xi <- dp$xi
      Omega <- dp$Omega
      alpha <- dp$alpha
      d <- length(alpha)
      fix <- fixed.comp
      h <- length(fix)
      
      # Calculate stuff that's common to all distributions, ready to get tau_vec
      omega <- sqrt(diag(Omega))
      omega1 <- omega[fix]
      omega2 <- omega[-fix]
      R <- cov2cor(Omega)
      R11 <- R[fix, fix, drop = FALSE]
      R12 <- R[fix, -fix, drop = FALSE]
      R21 <- R[-fix, fix, drop = FALSE]
      R22 <- R[-fix, -fix, drop = FALSE]
      alpha1 <- matrix(alpha[fix], ncol = 1)
      alpha2 <- matrix(alpha[-fix], ncol = 1)
      iR11 <- mnormt::pd.solve(R11)
      R22.1 <- R22 - R21 %*% iR11 %*% R12
      a.sum <- as.vector(t(alpha2) %*% R22.1 %*% alpha2)
      alpha1_2 <- as.vector(alpha1 + iR11 %*% R12 %*% alpha2)/sqrt(1 + a.sum)
      
      # Calculate vector of tau's
      tau_new_vec <- alpha1_2 * (fixed.values - xi[fix])/omega1
      
      
      # Calculate other common stuff, ready to get xi_new_matrix
      O11 <- Omega[fix, fix, drop = FALSE]
      O12 <- Omega[fix, -fix, drop = FALSE]
      O21 <- Omega[-fix, fix, drop = FALSE]
      O22 <- Omega[-fix, -fix, drop = FALSE]
      iO11 <- (1/omega1) * iR11 * rep(1/omega1, each = h)
      reg <- O21 %*% iO11
      O22.1 <- O22 - reg %*% O12
      omega22.1 <- sqrt(diag(O22.1))
      alpha2.1 <- as.vector((omega22.1/omega2) * alpha2)
      
      # Calculate matrix of xi's 
      xi_new_matrix <- t(xi[-fix] + reg %*% (fixed.values - xi[fix]))
      
      # Output the parameters
      dp_out <- list(
        xi = xi_new_matrix,
        Omega = O22.1,
        alpha = alpha2.1,
        tau = tau_new_vec
      )
      
      return(dp_out)
      
    }
    
    
    # Calculate timings for secondary cases, including when they are infected and when they get symptoms
    # This has to be drawn jointly so that serial interval is realistic
    # The code in the while loop is based on a vectorised version of sn::rmsn
    secondary_case_timing <- function(case_df, timings_dist, params_serial, params_symptom_timing) {
      
      # FOR TESTING:
      # timings_dist <- outbreak_setup_test$params$timings_dist
      # case_df <- tibble(case_id = 1:1000, n_secondary = sample(1:10, 1000, TRUE)) %>%
      #   mutate(primary_symptoms_index = sn::rsn(n = 1000, dp = timings_dist$marg_dist_primary@dp))
      # params_serial <- outbreak_setup_test$params$params_serial
      # params_symptom_timing <- outbreak_setup_test$params$params_symptom_timing
      
      # Check case_df has correct format, extract the variables as vectors
      if (!setequal(names(case_df), c("case_id", "primary_symptoms_index", "n_secondary"))) stop("case_df should contain case_id, primary_symptoms_index, n_secondary")
      case_id <- case_df$case_id; primary_symptoms_index <- case_df$primary_symptoms_index; n_secondary <- case_df$n_secondary
      
      # Repeat the symptoms vector n_secondary times each
      primary_index_rep <- rep(primary_symptoms_index, times = n_secondary)
      case_id_rep <- rep(case_id, times = n_secondary)
      row_id <- seq_along(primary_index_rep)   # for recording position of each item
      
      # How many secondary cases should we sample (to be used in while loop)
      n_times <- sum(n_secondary)
      
      # Return empty tibble if there's nothing to use
      if (n_times == 0) {return(tibble(case_id = character(0), 
                                       primary_symptoms_timing = numeric(0),
                                       serial_interval = numeric(0),
                                       secondary_symptoms_delay = numeric(0),
                                       secondary_case_timing = numeric(0),
                                       secondary_symptoms_timing = numeric(0)))}

      # Set up the while loop
      n_remaining <- n_times   # number of times remaining
      all_data <- NULL         # Create empty dataset to fill with new case data
      
      # Create remaining primary symptoms index
      primary_index_remaining <- primary_index_rep
      row_id_remaining <- row_id
      
      while (n_remaining > 0) {   # while loop ensures there are no secondary timings below cutoff point
        
        # print(n_remaining)
        
        # Calculate TAU and XI for the conditional distributions based on those primary symptom times
        conditional_params <- conditional_sn_vec(object = timings_dist$joint_dist, 
                                                 fixed.comp = 1, # fix first value (primary symptoms)
                                                 fixed.values = primary_index_remaining)
        
        tau_vec <- conditional_params$tau   # only tau and xi change when conditioning the distribution, the rest is constant
        xi_matrix <- conditional_params$xi
        
        # Calculate auxiliary parameters for conditional distribution
        cond_dist <- sn::conditionalSECdistr(
          object = timings_dist$joint_dist,
          fixed.comp = 1,   # fixes primary symptoms 
          fixed.values = -5 # fixes at what value (this doesn't matter here, can be any value)
        )
        lot <- sn:::dp2cpMv(dp = cond_dist@dp, family = "SN", aux = TRUE)
        
        # Randomly generate 2 normal variables to be psi-ed (correlated) and then tau-ed
        norm_vals <- rnorm(n_remaining * (2))
        y <- matrix(norm_vals, n_remaining, 2, byrow = FALSE) %*% chol(lot$aux$Psi)   # Psi is invariant across conditional distributions
        
        truncN <- qnorm(runif(n = n_remaining, min = pnorm(-tau_vec), max = 1))
        truncN <- matrix(rep(truncN, 2), ncol = 2)
        delta <- lot$aux$delta      #  delta is invariant across conditional distributions
        z <- delta * t(truncN) + sqrt(1 - delta^2) * t(y)
        y <- t(t(xi_matrix) + lot$aux$omega * z)
        
        # Data to be output from the while loop:
        data_out <- tibble(
          row_id = row_id_remaining,
          primary_index = primary_index_remaining,
          secondary_symp = y[ , 1],
          serial = y[ , 2]
        ) %>% 
          
          # CONVERT DISTRIBUTIONS FROM "INDEXES" TO QUANTILES
          mutate(
            primary_symptoms_timing_q =  sn::psn(x = - primary_index, dp = timings_dist$marg_dist_primary@dp), # note the minus sign  - to generate negative correlation between primary symptoms and serial interval!!!
            secondary_symptoms_delay_q = sn::psn(x = secondary_symp, dp = timings_dist$marg_dist_secondary@dp),
            serial_interval_q =          sn::psn(x = serial, dp = timings_dist$marg_dist_serial@dp)
          ) %>% 
          
          # CONVERT TO the desired distributions
          mutate(
            primary_symptoms_timing = qlnorm(primary_symptoms_timing_q, meanlog = params_symptom_timing$meanlog, sdlog = params_symptom_timing$sdlog),
            secondary_symptoms_delay = qlnorm(secondary_symptoms_delay_q, meanlog = params_symptom_timing$meanlog, sdlog = params_symptom_timing$sdlog),
            serial_interval = qgamma(p = serial_interval_q, shape = params_serial$shape, rate = params_serial$rate) + params_serial$shift
          ) %>% 
          
          # Get the index for secondary symptoms delay 
          mutate(
            secondary_symptoms_index = - sn::qsn(p = secondary_symptoms_delay_q, dp = timings_dist$marg_dist_primary@dp) # the minus sign matches with the calculation of primary_symptoms_timing_q above
          ) %>% 
          select(-secondary_symp, -serial, -ends_with("_q")) %>% 
          mutate(secondary_case_timing = primary_symptoms_timing + serial_interval - secondary_symptoms_delay,
                 secondary_symptoms_timing = secondary_case_timing + secondary_symptoms_delay)
        
        # Update the counters
        n_remaining <- sum(data_out$secondary_case_timing < timings_dist$min_infection)
        primary_index_remaining <- data_out$primary_index[data_out$secondary_case_timing < timings_dist$min_infection]
        row_id_remaining <- data_out$row_id[data_out$secondary_case_timing < timings_dist$min_infection]
        
        # Export the data from the while loop
        all_data <- bind_rows(
          all_data,
          data_out %>% filter(secondary_case_timing >= timings_dist$min_infection)
        )
        
      }
      
      # Sort the new data and export
      all_data_sorted <- all_data %>% 
        arrange(row_id) %>% 
        mutate(case_id = case_id_rep, .before = "row_id") %>% 
        select(-row_id, -primary_index)
        
      return(all_data_sorted)
      
    }
  
    # TESTING 
    # secondary_timing_test <- tibble(
    #   symptom_timing_index = symptom_timing_index(n = 100, marg_dist_primary = marg_dist_primary),
    #   n_secondary = sample(1:20, size = 100, replace = TRUE)
    # ) %>% 
    #   mutate(row_id = row_number()) %>% 
    #   rowwise(row_id) %>% 
    #   summarise(secondary_case_timing(symptom_timing_index, n_secondary = n_secondary, timings_dist = timings_dist, params_serial = serial_parameters))
    
    
    
    # mean_symptom_timing <- tibble(
    #   x = symptom_timing(rep(1, 100000), marg_dist_primary = marg_dist_primary)
    # ) %>% 
    #   {mean(.$x)} %>% 
    #   print()
    

    
    
    # Function that takes in the probability df, and the symptom severity, and outputs a vector of Bernoulli draws
    # Used for testing decision, isolation decisions etc.
    symptom_severity_to_yn <- function(symptom_severity, i_group, probs_df) {
      
      # FOR DEBUGGING:
      # symptom_severity <- symptom_severity(n = 100)
      # i_group <- as.integer(cut(runif(n = 100), breaks = c(0, 0.1, 0.4, 1), labels = 1:3))
      # probs_default <- c(0, 0.04, 0.2, 0.3, 0.4)
      # probs_df <- tibble(
      #   i_group = c(rep(1, 5), rep(2, 5), rep(3, 5)),
      #   symptom_severity = as.factor(rep(1:5, 3)),
      #   prob = c(probs_default * 0.1, probs_default * 1, probs_default * 2.5)
      # )
      
      n_cats <- nlevels(symptom_severity)
      n_cats_df <- nlevels(probs_df$symptom_severity)
      if (!is.integer(probs_df$symptom_severity)) stop("symptom_severity in probs_df is not an integer")
      if (n_cats_df != n_cats) stop("probs_df doesn't have same number of levels as symptom_severity")
      
      require(mc2d)
      
      tibble(
        symptom_severity = symptom_severity,
        i_group = i_group
      ) %>% 
        left_join(probs_df, by = c("symptom_severity", "i_group")) %>% 
        mutate(yn = as.logical(mc2d::rbern(nrow(.), p = prob))) %>%       # bernoulli with vectorised p
        .$yn
      
    }
    
    
    # FUNCTION that takes samples from another dataframe according to i_group
    # Used to sample from the hh_data while making sure you sample from only correct group
    sample_by_i_group <- function(dat, i_group, sample_dat, sample_vars) {
      
      # FOR DEBUGGING:
      # dat <- tibble(
      #   i_group = sample(1:3, 100, replace = TRUE),
      #   x = 1:100
      # )
      # 
      # sample_dat <- tibble(
      #   i_group = c(rep(1, 30), rep(2, 40), rep(3, 30)),
      #   test_results_delay = c(runif(30, 0, 1), runif(40, 5, 7), runif(30, 10, 15)),
      #   y = 1:100
      # )
      # 
      # sample_vars <- c("test_results_delay", "y")
      
      # Check whether sample var is in the original data, and throw an error if it is
      named <- !is.null(names(sample_vars))
      if (named) error <- names(sample_vars) %in% names(dat)
      else if (!named) error <- sample_vars %in% names(dat)
      if (sum(error) > 0) {
        print(names(dat))
        if (named) print(names(sample_vars))
        else if (!named) print(sample_vars)
        stop("The variable to be sampled is already in the baseline data")
      }
      
      dat %>% 
        mutate(row_id = row_number()) %>% 
        mutate(i_group_ = {{i_group}}) %>%  # workaround for the fact that i_group can't be used in group_modify
        group_by({{i_group}}) %>% 
        group_modify(
          ~ {
            sample_dat_i <- sample_dat[sample_dat$i_group == .x$i_group_[[1]], ] # get the sample_dat with only relevant i_group
            
            new_sampled <- slice_sample(sample_dat_i, n = nrow(.x), replace = TRUE) # sample from this sample_dat
            
            new_sampled_vars <- new_sampled %>% select(all_of(sample_vars)) # extract only the vars you want
            
            bind_cols(.x, new_sampled_vars) # bind on new sampled vars to the original dataframe
          }
        ) %>% 
        ungroup %>% 
        arrange(row_id) %>% 
        select(-i_group_, -row_id)
      
    }
    
    
    # TEST the sample_by_i_group function
    # sample_by_i_group(
    #   dat = tibble(i_group = sample(1:2, size = 100, replace = TRUE)),
    #   i_group = i_group,
    #   sample_dat = tibble(
    #     i_group = c(rep(1, 100), rep(2, 100)),
    #     gdp = 1:200,
    #     bob = 1:200
    #   ),
    #   sample_vars = c("gdp", "bob")
    # )
    
    
    
    # Recovery timing based on when you exhibited symptoms
    recovery_timing <- function(symptom_timings, recov_val) {
      # Assume same by group
      symptom_timings + recov_val     # assume you recover [no longer infectious] 10 days after exhibiting symptoms
    }
    
  
    # Number of contacts individual makes outside of home
    n_contacts_out <- function(i_group, contact_means, contact_dispersion) {
      # i_group <- sample(1:2, size = 1000, replace = TRUE)
      # contact_means <- c(7, 3)
    
      contact_means <- tibble(
        i_group = i_group
      ) %>% 
        left_join(
          tibble(i_group = 1:length(contact_means), contact_means = contact_means),
          by = "i_group"
        ) %>% 
        .$contact_means
      
      rnbinom(n = length(i_group), size = contact_dispersion, mu = contact_means)
      
    }
    
    # CHECK - see what the shape looks like
    # tibble(
    #   i_group = sample(1:2, size = 10000, replace = TRUE)
    # ) %>% 
    #   mutate(n_contacts_out = n_contacts_out(i_group,
    #                                          contact_means = c(7, 3),
    #                                          contact_dispersion = 0.5)) %>% 
    #   ggplot(aes(x = n_contacts_out, fill = factor(i_group))) + 
    #   geom_histogram(colour = "white") + 
    #   facet_wrap(~ factor(i_group))
    
    
    
    # Calculate number of CASES from number of CONTACTS
    secondary_cases_n <- function(i_group, n_contacts, sar) {
      
      # DEBUGGING:
      # i_group <- as.integer(cut(runif(n = 100), breaks = c(0, 0.1, 0.4, 1), labels = 1:3))
      # n_contacts <- sample(3:10, length(i_group), replace = TRUE)
      # n_contacts <- c(rep(1, 50), rep(100, 50))
      # sar <- 0.2
      
      n <- length(i_group)
      rbinom(n = n, size = n_contacts, prob = sar)
    }
    
    
    
    # CHECK - see what the shape looks like
    # tibble(
    #   i_group = sample(1:2, size = 10000, replace = TRUE)
    # ) %>%
    #   mutate(n_contacts_out = n_contacts_out(i_group,
    #                                          contact_means = c(30, 20),
    #                                          contact_dispersion = 0.5)) %>%
    #   mutate(n_secondary_cases_out = secondary_cases_n(i_group, n_contacts = n_contacts_out, sar = 0.1)) %>% 
    #   ggplot(aes(x = n_secondary_cases_out, fill = factor(i_group))) +
    #   geom_histogram(colour = "white") +
    #   facet_wrap(~ factor(i_group))

  
    
    
    # Vectorised version to calculate IDs of the new secondary infections based on HH data
    secondary_ids <- function(case_id, i_group, hh_ind_id, hh_id, sar_home,
                                  # n_secondary_cases_home, 
                                  n_secondary_cases_out, 
                                  hh_data, contact_weights) {
      
      # DEBUGGING
      # n <- 1000
      # case_id <- 1:n
      # hh_data <- gen_hh_data(n_pop = 100000,
      #                        group_props = c(0.3, 0.4, 0.3),
      #                        hh_size_data = data_save$hh_data_bogota %>% filter(i_group %in% 1:3),
      #                        p_hh_quarantine = c(0.5, 0.5, 0.8))
      # 
      # hh_data_sampled <- slice_sample(hh_data, n = n, replace = FALSE)
      # hh_ind_id <- hh_data_sampled$hh_ind_id
      # hh_id <- hh_data_sampled$hh_id
      # # n_secondary_cases_home <- hh_data_sampled %>%
      # #   mutate(n_secondary_cases_home = if_else(hh_size > 3, hh_size - 2, hh_size - 1)) %>%
      # #   .$n_secondary_cases_home
      # i_group <- hh_data_sampled$i_group
      # sar_home <- rep(c(0.3, 0.2, 0.1, 0.25), n/4)
      # n_secondary_cases_out <- rep(c(10, 5, 20, 0, 1), n/5)
      # contact_weights <- crossing(
      #   to = 1:3, from = 1:3
      # ) %>%
      #   mutate(contact_weight = c(10, 5, 2, 5, 7, 1, 2, 1, 4))
      
      
      # Keep only bits of hh_data that are needed
      hh_data_trimmed <- hh_data %>% select(hh_id, hh_size, secondary_hh_ind_id = hh_ind_id)
      
      # Put all the IDs into a df
      live_ids <- tibble(
        case_id = case_id,
        hh_ind_id = hh_ind_id,
        hh_id = hh_id,
        i_group = i_group,
        sar_home = sar_home,
        n_secondary_cases_out = n_secondary_cases_out,
      )
      

      
 
      # (1) HOME INFECTIONS
    

      # Calculate who's infected
      secondary_ids_home <- live_ids %>% 
        left_join(hh_data_trimmed, by = "hh_id") %>% 
        
        # Individuals can't infect themselves
        filter(hh_ind_id != secondary_hh_ind_id ) %>% 
        
        # Work out who actually gets infected using Bernoulli simulation
        mutate(home_potential_infected = bern_vec(nrow(.), p = sar_home)) %>% 
        filter(home_potential_infected) %>% 
        select(-home_potential_infected, -hh_size, -sar_home)
      
      # Merge back on to get appropriate missing values
      secondary_ids_home_complete <- live_ids %>% 
        select(-sar_home) %>% 
        left_join(
          secondary_ids_home, by = c("case_id", "hh_id", "hh_ind_id", "i_group", "n_secondary_cases_out")
        ) %>% 
        mutate(
          infection_type = "home"
        )
    

      # (2) OUTSIDE HOME INFECTIONS
      
      # Create 1 version for each i_group of the hh_data with the appropriate probability weights (taken from beta matrix)
      hh_data_weighted <- map(
        .x = unique(contact_weights$from), 
        function(x) {
          left_join(hh_data, contact_weights[contact_weights$from == x, ], by = c("i_group" = "to"))
        }
      ) %>% 
        map(select, hh_id, hh_ind_id, contact_weight)

      # Sample the ID of new out infections 
      secondary_ids_out <- live_ids %>% 
        mutate(i_group_ = i_group,
               row_id = row_number()) %>% 
        group_by(i_group) %>% 
        group_modify(.f = ~ {
          hh_data_w <- hh_data_weighted[[ .x$i_group_[[1]] ]]
          probs <- hh_data_w$contact_weight
          # probs[hh_data_w$hh_id == .$hh_id] <- 0 # REMOVED setting probability to 0 when hh_id is the same [too big computation cost, little benefit when pop is large]
          sampled_hh_ind_id <- sample(x = hh_data_w$hh_ind_id, size = sum(.x$n_secondary_cases_out), prob = probs, replace = TRUE)
          
          tibble(
            case_id = rep(.x$case_id, times = .x$n_secondary_cases_out),
            hh_ind_id = rep(.x$hh_ind_id, times = .x$n_secondary_cases_out),
            hh_id = rep(.x$hh_id, times = .x$n_secondary_cases_out),
            secondary_hh_ind_id = sampled_hh_ind_id
          )
        }) %>% 
        # If same secondary case is selected twice for a given ind_case, remove it [only small effect]
        distinct(case_id, secondary_hh_ind_id, .keep_all = TRUE)
      
      # Join back onto live_ids
      secondary_ids_out_complete <- live_ids %>% 
        left_join(secondary_ids_out, by = c("case_id", "hh_ind_id", "hh_id", "i_group")) %>% 
        mutate(infection_type = "out")
      
      # (3) COMBINE TWO TYPES OF INFECTION
      # Bind the two types of infections and clean up
      secondary_infections <- bind_rows(
        secondary_ids_out_complete,
        secondary_ids_home_complete
      ) %>% 
        select(case_id, hh_ind_id, hh_id, secondary_hh_ind_id, infection_type) %>% 
        left_join(
          tibble(hh_ind_id = !!hh_ind_id) %>% mutate(row_id = row_number()),
          by = "hh_ind_id"
        ) %>% 
        arrange(row_id) %>% 
        left_join(
          hh_data %>% select(secondary_hh_id = hh_id, secondary_i_group = i_group, secondary_hh_ind_id = hh_ind_id, secondary_hh_size = hh_size, secondary_would_quarantine = would_quarantine),
          by = "secondary_hh_ind_id"
        ) %>% 
        
        # Remove cases where an "out" infection is from inside the household [causes bug with duplicates]
        filter(!(infection_type == "out" & secondary_hh_id == hh_id) | is.na(secondary_hh_id)) %>% 
        
        select(-row_id) %>% 
        # group_by(hh_ind_id) %>% 
        # Change the way secondary_case_id works (speed up calculation by not having to group)
        # but it's now (case_id)_(row_number_in_entire_dataset)
        mutate(secondary_case_id = if_else(!is.na(secondary_hh_ind_id), paste0(case_id, "_", row_number()), NA_character_)) %>% 
        ungroup
      

    
    
      return(secondary_infections)
      
    } 
    
    
  
    # CALCULATE ISOLATION TIMINGS dependent on whether it's pre or post symptomatic
    # Returns test timing when tested before symptoms, then symptom timing
    # Or test timing if tested after symptoms
    # Basically making use of the fact that probability of isolating when PRESYMPTOMATIC is same as probability of isolating when ASYMPTOMATIC
    # and the fact that people may be tested positive but not isolate while presymptomatic, but then isolate once they are symptomatic
    calc_isolation_timing <- function(test_timing, symptom_timing, isolate_presymp, isolate_symp) {
      
      # DEBUGGING
      # test_timing <- live_cases_w_random$random_test_timing + 5
      # symptom_timing <- live_cases_w_random$symptom_timing
      # isolate_presymp <- live_cases_w_random$isolate_after_test_presymp
      # isolate_symp <- live_cases_w_random$isolate_after_test_symp
      
      test_before <- test_timing < symptom_timing
      
      # CASE 1 - test is before symptoms
      
        # (i) Behaviour presymptomatic
        isolate_before_presymp <- if_else(
          test_before & isolate_presymp,
          test_timing, # isolate when receive test
          NA_real_
        )
      
        # (ii) Behaviour when symptomatic
        isolate_before_symp <- if_else(
          test_before & is.na(isolate_before_presymp) & isolate_symp,
          symptom_timing, # isolate when you get symptoms
          NA_real_
        )
      
      # CASE 2 - test is after symptoms
      isolate_after <- if_else(
        !test_before & isolate_symp,
        test_timing,
        NA_real_
      )
      
      values_by_row <- as.numeric(!is.na(isolate_before_presymp)) + as.numeric(!is.na(isolate_before_symp)) + as.numeric(!is.na(isolate_after))
      
      if (sum(values_by_row > 1) > 0) stop("non mutually exclusive cases in calc_isolation_timing")
      
      isolate_timing <- pmin(isolate_before_presymp, isolate_before_symp, isolate_after, na.rm = TRUE)
      
      return(isolate_timing)      
      
    }
    
    
    # Function to CALCULATE WHO IS QUARANTINED [primary and secondary]
    calc_quarantine <- function(secondary_cases_df, hh_data) {
      
      # secondary_cases_df <- new_secondary_merged
      # hh_data
      
      # Calculate list of primary quarantined at individual level
      primary_quarantined <- left_join(
        secondary_cases_df %>% select(case_id, secondary_hh_ind_id, hh_id, secondary_case_timing, would_quarantine) %>% filter(!is.na(secondary_hh_ind_id)),
        hh_data %>% select(hh_id, hh_ind_id_for_q  = hh_ind_id, matches("isolate_interval")),
        by = c("hh_id" = "hh_id"),
        # .merge = TRUE
      ) %>% 
        mutate(
          quarantined_1 = data.table::between(secondary_case_timing, lower = isolate_interval_1_l, upper = isolate_interval_1_r, NAbounds = NA),
          quarantined_2 = data.table::between(secondary_case_timing, lower = isolate_interval_2_l, upper = isolate_interval_2_r, NAbounds = NA)
        ) %>% 
        mutate(
          primary_quarantined = would_quarantine & (quarantined_1 | quarantined_2),
        ) %>% 
        filter(primary_quarantined == TRUE) %>% 
        select(case_id, secondary_hh_ind_id, primary_quarantined) %>% 
        distinct()
      
      # Match on details of all household members of the secondary case
      secondary_quarantined <- left_join(
        secondary_cases_df %>% select(case_id, secondary_hh_ind_id, secondary_hh_id, secondary_case_timing, secondary_would_quarantine) %>% filter(!is.na(secondary_hh_ind_id)),
        hh_data %>% select(hh_id, hh_ind_id_for_q  = hh_ind_id, matches("isolate_interval")),
        by = c("secondary_hh_id" = "hh_id"),
        # .merge = TRUE
      ) %>% 
        mutate(
          quarantined_1 = data.table::between(secondary_case_timing, lower = isolate_interval_1_l, upper = isolate_interval_1_r, NAbounds = NA),
          quarantined_2 = data.table::between(secondary_case_timing, lower = isolate_interval_2_l, upper = isolate_interval_2_r, NAbounds = NA)
        ) %>% 
        mutate(
          secondary_quarantined = secondary_would_quarantine & (quarantined_1 | quarantined_2)
        ) %>% 
        filter(secondary_quarantined == TRUE) %>% 
        select(case_id, secondary_hh_ind_id, secondary_quarantined) %>% 
        distinct()
      
      new_secondary_with_quarantine <- secondary_cases_df %>%
        left_join(primary_quarantined, by = c("case_id", "secondary_hh_ind_id"), .merge = TRUE) %>% 
        left_join(secondary_quarantined, by = c("case_id", "secondary_hh_ind_id")) %>% 
        mutate(across(c(primary_quarantined, secondary_quarantined), na_to_false))
    
      return(new_secondary_with_quarantine)
      
    }
    
    
    
# Functions for drawing multiple distributions at different stages ----------------------------
    
    # Function to update all timing variables by reducing t
    # updates timing on all variables that end in timing, and the interval_1_r-style variables as well
    update_timing <- function(df, dt) {
      
      if (is.null(df)) return(NULL) 
      else {
        df %>% 
          mutate(across(ends_with("_timing") | ends_with("_timing_2") | matches("_interval_(1|2)_(l|r)"), ~ .x - dt))
      }
      
    }
    
    # Function that calculates population size of each group, from group_props and n_pop
    calc_pops <- function(group_props, n_pop) {
      n_groups <- length(group_props)
      
      pops_but_one <- as.integer(round(n_pop * group_props[1:(n_groups - 1)], 0))
      pops <- c(
        pops_but_one,
        n_pop - sum(pops_but_one)
      )
      
      if (sum(pops) != n_pop) {
        print(n_pop)
        print(sum(pops))
        stop("Something going wrong with pops in gen_hh_data")
      }
      
      return(pops)
    }
    
    
    # Draw the characteristics of households
    gen_hh_data <- function(n_pop, group_props, hh_size_data, p_hh_quarantine) {
    
      # FOR DEBUGGING:
      # n_pop <- 80000; group_props <- data_save$group_props
      # hh_size_data <- data_save$hh_data_bogota
      # hh_size_data = tibble(
      #   i_group = c(rep(1, 10000), rep(2, 10000), rep(3, 10000)),
      #   hh_size = c(
      #     rbinom(10000, size = 10, prob = 0.3) + 1,
      #     rbinom(10000, size = 5, prob = 0.5) + 1,
      #     rbinom(10000, size = 5, prob = 0.5) + 1
      #   )
      # )
      # p_hh_quarantine = c(0.5, 0.5, 0.8, 0.7, 0.6, 0.6)
      
      pops <- calc_pops(group_props, n_pop)
      n_groups <- length(group_props)
      
      # GENERATE HH SIZE FOR EACH HH
      hh_sizes <- list()
      
      for (i in 1:n_groups) {
        
        data <- hh_size_data %>% filter(i_group == i) %>% .$hh_size
        sampled_data <- sample(data, size = pops[[i]], replace = TRUE) # if you sample at least pops[[i]] households, there will always definitely be enough to cover all individuals
        sampled_trimmed <- sampled_data[cumsum(sampled_data) < pops[[i]] ]
        
        # Make sure last item adds up to make a round number
        if (pops[[i]] - sum(sampled_trimmed) > 0) {
          sampled_trimmed[[length(sampled_trimmed) + 1]] <- pops[[i]] - sum(sampled_trimmed)
        }
        
        if (sum(sampled_trimmed) - pops[[i]] > 10e-3) stop("Error in hh size generation")
        
        hh_sizes[[i]] <- tibble(
          i_group = i,
          hh_size = sampled_trimmed
        )
      }
      
      hh_sizes_bind <- hh_sizes %>% 
        bind_rows() %>% 
        mutate(hh_id = as.character(row_number())) %>% 
        relocate(hh_id)
      
      # Add in quarantine variables
      hh_sizes_quarantine <- hh_sizes_bind %>% 
        left_join(
          tibble(i_group = 1:n_groups,
                 p_hh_quarantine = p_hh_quarantine),
          by = "i_group"
        ) %>% 
        mutate(
          would_quarantine = bern_vec(n = nrow(.), p = p_hh_quarantine)
        )
      
      hh_sizes_long <- hh_sizes_quarantine %>% 
        rowwise() %>% 
        mutate(hh_ind_id = list(1:hh_size)) %>% 
        unnest(hh_ind_id) %>% 
        mutate(hh_ind_id = paste0(hh_id, "_", hh_ind_id))
      
      return(hh_sizes_long)
    }
    
    # TEST THE FUNCTION
    # gen_hh_data(
    #   n_pop = 10000,
    #   group_props = c(0.2, 0.8),
    #   hh_size_data = list(
    #     `1` = rbinom(10000, size = 10, prob = 0.3) + 1,
    #     `2` = rbinom(10000, size = 5, prob = 0.5) + 1
    #   )
    # ) %>% 
    #   group_by(i_group) %>% 
    #   summarise(total = sum(hh_size))
    

    # Initial draw of symptoms, recovery time, isolation behaviour
    draw_symptoms_recovery <- function(ids,
                                       i_group,
                                       hh_id,
                                       hh_ind_id,
                                       hh_size,
                                       would_quarantine,
                                       symptom_timing_index,
                                       contact_tracing_test_timing = NULL,    
                                       contact_tracing_results_timing = NULL,
                                       marg_dist_primary,
                                       test_sensitivity_data,
                                       infection_timings,              # inputs from the other draws
                                       test_delay_data,
                                       # test_choice_delay_data,
                                       # test_results_delay_data, 
                                       params_symptom_timing,
                                       recov_val,
                                       probs_self_test_df, probs_isolate_symptoms_df, probs_isolate_test_df, probs_isolate_ct_df,
                                       infectiousness_by_symptoms) {
      
      # Needs to be run once for new cases only
      # DO NOT INCLUDE things that will need to be updated periodically
      # e.g. contact tracing, random testing, isolation details, 
      
      # DEBUGGING:
      # ids <- 1:1000
      # i_group <- as.integer(cut(runif(n = 1000), breaks = c(0, 0.33333333, 0.6666666, 1), labels = 1:3))
      # infection_timings <- 0
      # contact_tracing_test_timing = NULL
      # contact_tracing_results_timing = NULL
      # random_test_yn = NULL
      # random_test_timing = NULL
      # isolate_random_test_timing = NULL
      # test_results_delay_data <- 1.5; recov_val <- 5
      # probs_default <- c(0.3, 0.6, 0.8, 0.9, 1)
      # probs_df_basic <- tibble(
      #   i_group = c(rep(1, 5), rep(2, 5), rep(3, 5)),
      #   symptom_severity = as.factor(rep(1:5, 3))
      # )
      # test_sensitivity_data <- test_sensitivity_data
      # probs_self_test_df <- probs_df_basic %>% mutate(prob = c(probs_default * 0.1, probs_default * 0.5, probs_default * 1))
      # probs_isolate_symptoms_df <- probs_df_basic %>% mutate(prob = c(probs_default * 0.1, probs_default * 0.3, probs_default * 0.5))
      # probs_isolate_test_df <- probs_df_basic %>% mutate(prob = c(probs_default * 0.6, probs_default * 0.8, probs_default * 1))
      # probs_isolate_ct_df <- probs_df_basic %>% mutate(prob = c(probs_default * 0.3, probs_default * 0.4, probs_default * 0.5))
      # infectiousness_by_symptoms = c(0.7, 1, 1, 1, 1)
      # contact_tracing_test_timing <- sample(
      #   c(runif(50, min = -2, max = 20), NA, NA, NA, NA, NA), 
      #   size = 1000, 
      #   replace = TRUE
      # )
      # contact_tracing_results_timing <- contact_tracing_test_timing + 2
      
      # CHECK PARAMETERS 
      n <- length(ids)
      if (length(i_group) != length(ids)) stop ("length(i_group) != length(ids)")
      if (length(unique(ids)) != n) warning("Duplicate IDs when calculating symptoms")

      # Fill in the null values
      if (is.null(hh_id)) hh_id <- NA_real_
      if (is.null(hh_size)) hh_size <- NA_real_
      if (is.null(hh_ind_id)) hh_ind_id <- NA_character_
      if (is.null(would_quarantine)) would_quarantine <- NA
      if (is.null(contact_tracing_test_timing)) contact_tracing_test_timing <- NA_real_
      if (is.null(contact_tracing_results_timing)) contact_tracing_results_timing <- NA_real_
      if (is.null(symptom_timing_index)) symptom_timing_index <- NA_real_
      if (length(infectiousness_by_symptoms) != 2) stop ("Infectiousness by symptoms should be length 2 (same as number of symptom categories)")
      
      # Match infectiousness to symptom severity
      symptom_to_infectiousness <- tibble(
        symptom_severity = 1L:2L,
        infectiousness = infectiousness_by_symptoms
      )
      
      rm(infectiousness_by_symptoms) # remove the vector, so mutate knows to pick it up from the dataframe, not from the external object
      
      # CREATE SYMPTOMS dataframe
      tibble(
        case_id = as.character(ids),
        hh_id = hh_id,
        hh_ind_id = hh_ind_id,
        i_group = i_group,
        hh_size = hh_size,
        would_quarantine = would_quarantine,
        infection_timing = infection_timings,
        symptom_timing_index = symptom_timing_index,
        contact_tracing_test_timing = contact_tracing_test_timing,
        contact_tracing_results_timing = contact_tracing_results_timing,
        random_test_yn = NA,
        random_test_yn_ever = FALSE,
        random_test_timing = NA_real_,
        random_test_false_negative = NA,
        isolate_random_test_timing = NA_real_
      ) %>% 
        
        # CONTACT TRACING - test results
        # Contact tracing doesn't work [test negative] if they test people before they're infectious
        # or if the test is a false negative
        mutate(
          ct_test_negative = contact_tracing_test_timing < infection_timing,
          ct_test_days_exposure = round(contact_tracing_test_timing - infection_timing, digits = 0)
        ) %>% 
        left_join(test_sensitivity_data, by = c("ct_test_days_exposure" = "days_exposure")) %>% 
        rename(sensitivity_ct = sensitivity) %>% 
        mutate(ct_false_negative = bern_vec(n = nrow(.), p = 1 - sensitivity_ct, NA_error = FALSE)) %>% 
        
        # Symptoms & recovery
        mutate(
          symptom_severity = symptom_severity(i_group),
          symptom_timing_index = if_else(is.na(symptom_timing_index), # only fill in if it's empty (can be filled in from the secondary_ calculation)
                                         symptom_timing_index(n = nrow(.), marg_dist_primary = marg_dist_primary),
                                         symptom_timing_index),
          symptom_timing = infection_timing + symptom_timing(symptom_timing_index, marg_dist_primary = marg_dist_primary, params_symptom_timing = params_symptom_timing),
          recovery_timing = infection_timing + recovery_timing(symptom_timing, recov_val = recov_val)
        ) %>% 
        left_join(symptom_to_infectiousness, by = "symptom_severity") %>% 
        
        # SELF TESTING
        # When do you take the test?
        sample_by_i_group(
          i_group = i_group,
          sample_dat = test_delay_data,
          sample_vars = c("test_choice_delay", "test_results_delay")
        ) %>% 
        # sample_by_i_group(
        #   i_group = i_group,
        #   sample_dat = test_results_delay_data,
        #   sample_vars = "test_results_delay"
        # ) %>% 
        # How sensitive is the test
        mutate(
          self_test_days_exposure = round(symptom_timing + test_choice_delay - infection_timing, digits = 0)
        ) %>% 
        left_join(test_sensitivity_data, by = c("self_test_days_exposure" = "days_exposure")) %>% 
        rename(sensitivity_self_test = sensitivity) %>%
        mutate(
          self_test_false_negative = bern_vec(n = nrow(.), p = 1 - sensitivity_self_test, NA_error = FALSE),
          self_test_yn = symptom_severity_to_yn(symptom_severity = symptom_severity, i_group = i_group, probs_df = probs_self_test_df)
        ) %>% 
        # When do you get result?
        mutate(
          test_result_timing = if_else(self_test_yn, symptom_timing + test_choice_delay + test_results_delay, NA_real_),
        ) %>% 
        
        # Isolation behaviour:
        mutate(
          
          # Isolation behaviour when presymptomatic
          isolate_after_ct_presymp = symptom_severity_to_yn(symptom_severity = 1, 
                                                            i_group, probs_df = probs_isolate_ct_df), # isolation if you hear that you are contact traced
          isolate_after_test_presymp = symptom_severity_to_yn(symptom_severity = 1, 
                                                              i_group, probs_df = probs_isolate_test_df), # used for whether people isolate after contact tracaing test, self test, or random test
          
          # Isolation behaviour when symptomatic
          isolate_after_ct_symp = symptom_severity_to_yn(symptom_severity, i_group, probs_df = probs_isolate_ct_df), # isolation if you hear that you are contact traced
          isolate_after_test_symp = symptom_severity_to_yn(symptom_severity, i_group, probs_df = probs_isolate_test_df), # used for whether people isolate after contact tracaing test, self test, or random test
          
          # Isolation just on symptoms
          isolate_after_symptoms = symptom_severity_to_yn(symptom_severity, i_group, probs_df = probs_isolate_symptoms_df),
          
          # Isolation timings from self testing [always symptomatic]
          isolate_symptoms_timing = if_else(isolate_after_symptoms, symptom_timing, NA_real_),
          isolate_self_test_timing = if_else(self_test_yn & isolate_after_test_symp & !self_test_false_negative, test_result_timing, NA_real_)
          
        )
      
    }
    
    
    # FUNCTION that draws whether you are randomly tested at time t
    draw_random_testing <- function(live_cases, test_sensitivity_data, alpha, test_delay_data) {
      
      # Needs to be run on ALL live cases every period
      
      # DEBUGGING:
      # live_cases <- draw_symptoms_debug; test_sensitivity <- 0.85
      # alpha <- c(0.05, 0.1, 0.2);
      
      n_cases <- nrow(live_cases)
      
      alpha_df <- tibble(
        i_group = 1:length(alpha),
        alpha = alpha
      )
      rm(alpha)
      
      live_cases_w_random <- live_cases %>% 
        
        # Match to get individual alphas
        select(-any_of("alpha")) %>%  # remove existing alpha column so it can overwrite
        left_join(alpha_df, by = "i_group") %>% 
        
        # Only do random test if you're not currently being tested, and you haven't previously tested positive
        # Because we're only looking at infected people here, if you are ever randomly tested, you will never be tested again (because you were positive by construction)
        mutate(
          currently_testing = case_when(
            self_test_yn & symptom_timing <= 0 & test_result_timing >= 0    ~       TRUE, # testing period for self test
            contact_tracing_test_timing <= 0 & contact_tracing_results_timing >= 0 ~   TRUE,       # testing period for contact tracing
            random_test_yn & random_test_timing >= 0 ~ TRUE,                              # testing period for random test [no delay because test starts at 0]
            TRUE ~ FALSE
          ),
          
          # New assumption [used to be previously_tested_positive]: people are only randomly tested if they've never been tested before
          previously_tested = case_when(
            self_test_yn & test_result_timing < 0 ~ TRUE,
            contact_tracing_results_timing < 0 ~ TRUE,
            random_test_timing < 0 ~ TRUE,
            TRUE ~ FALSE
          )
        ) %>% 
        mutate(
          update_random_testing = (!currently_testing & !previously_tested) | is.na(random_test_yn)
        ) %>% 
        mutate(
          random_test_yn_lag1 = random_test_yn,
          random_test_yn = if_else(random_test_yn & random_test_timing < 0, FALSE, random_test_yn), # if random test has finished, set random_test_yn to false
          random_test_yn = if_else(update_random_testing,
                                   bern_vec(n_cases, p = alpha),
                                   random_test_yn)
        ) %>% 
        mutate(
          new_random_test = random_test_yn > random_test_yn_lag1 | (is.na(random_test_yn_lag1) & random_test_yn),
          random_test_days_exposure = round(0 - infection_timing, digits = 0)
        ) %>% 
        left_join(test_sensitivity_data, by = c("random_test_days_exposure" = "days_exposure")) %>% 
        rename(sensitivity_random_test = sensitivity) %>% 
        mutate(
          random_test_false_negative = if_else(new_random_test, bern_vec(nrow(.), p = 1 - sensitivity_random_test, NA_error = FALSE), random_test_false_negative),
          random_test_yn_ever = if_else(random_test_yn == TRUE & random_test_yn_ever == FALSE, TRUE, random_test_yn_ever)
        ) %>% 
        sample_by_i_group(
          i_group = i_group,
          sample_dat = test_delay_data,
          sample_vars = c("random_test_delay" = "test_results_delay")
        ) %>% 
          
        mutate(
          # Tested at t = now, receive results at random_test_timing
          random_test_timing = if_else(random_test_yn & new_random_test,
                                       random_test_delay,
                                       random_test_timing),
          
          # ISOLATION BEHAVIOUR
          # random_test_presymptomatic = random_test_timing < symptom_timing,
          new_random_test_positive =  random_test_yn & new_random_test & !random_test_false_negative,
          
          # Calculate isolation timing
          isolate_random_test_timing = if_else(new_random_test_positive,
                                               calc_isolation_timing(
                                                 test_timing = random_test_timing,
                                                 symptom_timing = symptom_timing,
                                                 isolate_presymp = isolate_after_test_presymp,
                                                 isolate_symp = isolate_after_test_symp
                                               ),
                                               isolate_random_test_timing) # otherwise return original value
        ) %>% 
        # Update currently testing
        mutate(
          currently_testing = case_when(
            self_test_yn & symptom_timing <= 0 & test_result_timing >= 0    ~       TRUE, # testing period for self test
            contact_tracing_test_timing <= 0 & contact_tracing_results_timing >= 0 ~   TRUE,       # testing period for contact tracing
            random_test_yn & random_test_timing >= 0 ~ TRUE,                              # testing period for random test [no delay because test starts at 0]
            TRUE ~ FALSE
          )
        ) %>% 
        select(-sensitivity_random_test, -random_test_delay)
      
      live_cases_w_random
      
    }
    

    # Function for removing cases that (i) occur after recovery, or (ii) have 0 potential cases
    filter_non_cases <- function(secondary_cases) {
      filtered <- secondary_cases %>% 
        filter(!is.na(secondary_case_id)) %>% 
        mutate(recovery_before_secondary = recovery_timing < secondary_case_timing) %>% 
        group_by(case_id) %>% 
        mutate(n_recovery_before_secondary_home = sum(recovery_before_secondary[infection_type == "home"], na.rm = TRUE)) %>% 
        mutate(n_recovery_before_secondary_out = sum(recovery_before_secondary[infection_type == "out"], na.rm = TRUE)) %>% 
        ungroup %>% 
        mutate(
          # n_secondary_cases_home = n_secondary_cases_home - n_recovery_before_secondary_home,
          n_secondary_cases_out = n_secondary_cases_out - n_recovery_before_secondary_out
        ) %>% 
        filter(!recovery_before_secondary | is.na(secondary_case_timing)) %>% 
        select(-recovery_before_secondary, -n_recovery_before_secondary_home, -n_recovery_before_secondary_out)
      
      if (nrow(filtered) == 0) return(NULL)
      else return(filtered)
    }
    
    
    
    # Calculates all the isolation timings 
    gen_isolation_all <- function(live_cases) {
      
      isolation_timings <- live_cases %>% 
        
        # Calculate isolation timings
        # Isolation timing based on CT [this has to be here so that it is updated as a new column in update_isolation]
        mutate(
          isolate_ct_test_timing = if_else(!is.na(contact_tracing_test_timing),
                                           calc_isolation_timing(
                                             test_timing = contact_tracing_test_timing,
                                             symptom_timing = symptom_timing,
                                             isolate_presymp = isolate_after_ct_presymp,
                                             isolate_symp = isolate_after_ct_symp
                                           ),
                                           NA_real_),
          
          isolate_ct_results_timing = if_else(!is.na(contact_tracing_test_timing) & !ct_false_negative & !ct_test_negative,
                                              calc_isolation_timing(
                                                test_timing = contact_tracing_results_timing,
                                                symptom_timing = symptom_timing,
                                                isolate_presymp = isolate_after_test_presymp,
                                                isolate_symp = isolate_after_test_symp
                                              ),
                                              NA_real_)
        ) %>% 
        
        # "Deisolation" timing - if people are isolating after ct, but then get a negative test result - they deisolate
        mutate(
          deisolate_timing = if_else(!is.na(isolate_ct_test_timing) & (ct_test_negative | ct_false_negative) & contact_tracing_results_timing > isolate_ct_test_timing, 
                                     contact_tracing_results_timing, 
                                     NA_real_),
          deisolate_timing = case_when(
            isolate_symptoms_timing < deisolate_timing ~ NA_real_,       # if people are already isolating due to symptoms, they won't deisolate
            isolate_self_test_timing < deisolate_timing ~ NA_real_,      # if people are already isolating due to another test, assume they won't deisolate
            isolate_random_test_timing < deisolate_timing ~ NA_real_, 
            TRUE ~ deisolate_timing
          )
        ) %>% 
        
        # Calculate isolation
        mutate(
          isolation_timing = pmin(isolate_random_test_timing, isolate_symptoms_timing, isolate_ct_test_timing, isolate_ct_results_timing, isolate_self_test_timing, na.rm = TRUE),
          isolation_timing = if_else(isolation_timing > recovery_timing, NA_real_, isolation_timing)
        ) %>% 
        
        # Second isolation after deisolation is possible
        mutate(
          across(c(isolate_random_test_timing, isolate_symptoms_timing, isolate_ct_test_timing, isolate_ct_results_timing, isolate_self_test_timing),
                 list(`2` = ~ if_else(. > deisolate_timing, ., NA_real_)))
        ) %>% 
        mutate(
          isolation_timing_2 = pmin(isolate_random_test_timing_2, isolate_symptoms_timing_2, 
                                    isolate_ct_test_timing_2, isolate_ct_results_timing_2, isolate_self_test_timing_2, na.rm = TRUE)
        ) %>% 
        
        # Only contact traced if you're detected, and if the row is a new potential case (so secondary case timing exists, and contact still occurs despite isolation)
        mutate(
          detected = self_test_yn | !is.na(contact_tracing_results_timing) | random_test_yn
        ) %>% 
        
        # Add hh-quarantining-relevant isolations (only isolate_ct_results_timing; isolate_self_test_timing)
        mutate(
          isolation_timing_hh = pmin(isolate_ct_results_timing, isolate_self_test_timing, na.rm = TRUE),
          isolation_timing_hh_2 = pmin(isolate_ct_results_timing_2, isolate_self_test_timing_2, na.rm = TRUE)
        ) %>% 
        
        # Calculate isolate intervals [used for the household-level quarantining]
        mutate(
          isolate_interval_1_l = if_else(!is.na(isolation_timing_hh), isolation_timing_hh, NA_real_),
          isolate_interval_1_r = case_when(
            !is.na(isolation_timing_hh) & is.na(deisolate_timing) ~ recovery_timing, # when they don't deisolate
            !is.na(isolation_timing_hh) & !is.na(deisolate_timing) ~ deisolate_timing, # when they do deisolate
            is.na(isolation_timing) ~ NA_real_
          ),
          isolate_interval_2_l = if_else(!is.na(deisolate_timing) & !is.na(isolation_timing_hh_2), isolation_timing_hh_2, NA_real_),
          isolate_interval_2_r = if_else(!is.na(deisolate_timing) & !is.na(isolation_timing_hh_2), recovery_timing, NA_real_)
        )
      

      return(isolation_timings)

    }
    
    
    # CONVENIENCE FUNCTION that calculates which columns are added when running gen_isolation_all
    # so that we know which ones to remove when recalculating isolation timings
    # CREATE DUMMY DATASETS to calculate which columns to remove
    infer_isolation_cols <- function(params) {
      probs_df_dummy <- crossing(i_group = 1:2, symptom_severity = as.integer(1:2)) %>% 
        mutate(prob = 0.5)
      
      marg_dist_primary_dummy <- sn::makeSECdistr(
        dp = list(xi = c(0, 0, 0), 
                  Omega = diag(c(1, 2, 3)), 
                  alpha = c(0, 0, 0)),
        family = "SN",
        compNames = c("primary_symp", "secondary_symp", "serial")
      ) %>% 
        marginalSECdistr(comp = 1)
      
      dummy_dataset <- draw_symptoms_recovery(ids = 1:100,
                                              i_group = rep(1:2, 50) ,
                                              hh_id = NULL,
                                              hh_ind_id = NULL,
                                              hh_size = NULL,
                                              would_quarantine = NULL,
                                              symptom_timing_index = NULL,
                                              contact_tracing_test_timing = NULL,    
                                              contact_tracing_results_timing = NULL,
                                              marg_dist_primary = marg_dist_primary_dummy,
                                              params_symptom_timing = list(meanlog = 1.5, sdlog = 0.5),
                                              test_sensitivity_data = tibble(days_exposure = as.integer(1:40), sensitivity = 0.8),
                                              infection_timings = 0,              # inputs from the other draws
                                              test_delay_data = tibble(i_group = 1:2, test_choice_delay = as.double(1:2), test_results_delay = as.double(3:4)),
                                              # test_choice_delay_data = tibble(i_group = 1:2, test_choice_delay = as.double(1:2)),
                                              # test_results_delay_data = tibble(i_group = 1:2, test_results_delay = as.double(3:4)), 
                                              recov_val = c(10, 10),
                                              probs_self_test_df = probs_df_dummy, probs_isolate_symptoms_df = probs_df_dummy, 
                                              probs_isolate_test_df = probs_df_dummy, probs_isolate_ct_df = probs_df_dummy,
                                              infectiousness_by_symptoms = c(0.7, 1)) %>% 
        draw_random_testing(test_sensitivity_data = tibble(days_exposure = as.integer(1:40), sensitivity = 0.8), alpha = c(0, 0, 0), test_delay_data = tibble(i_group = 1:2, test_choice_delay = as.double(1:2), test_results_delay = as.double(3:4)))
      
      dummy_dataset_w_isolation <- dummy_dataset %>% 
        gen_isolation_all()
      
      cols_no_isolation <- dummy_dataset %>% names
      cols_w_isolation <- dummy_dataset_w_isolation %>% names
      isolation_cols <- cols_w_isolation[! (cols_w_isolation %in% cols_no_isolation)]
      
      isolation_cols
      
    }
    
    
 
    
    

   
    
    

    
    # FUNCTION that takes live cases and draws all information on their potential
    # secondary infections
    draw_secondary_cases <- function(new_cases_w_isolation, hh_data,
                                     p_contact_if_isolated_home, p_contact_traced, 
                                     test_delay_data, ct_delay_data,
                                     contact_weights, group_props,
                                     contact_means,
                                     contact_dispersion,
                                     sar_home, sar_out,
                                     timings_dist,
                                     params_serial,
                                     params_symptom_timing
                                     # r0, dispersion
                                     ) {
      
      # FOR DEBUGGING:
      # new_cases_w_isolation <- live_with_hh
      # contact_means =  c(30, 20, 20)
      # sar_home = c(0.4, 0.3, 0.3); sar_out = c(0.2, 0.2, 0.2)
      # p_contact_if_isolated_home <- c(0.3, 0.2, 0.1)
      # p_contact_traced <- c(0.1, 0.2, 0.2)
      # contact_weights <- crossing(to = 1:3, from = 1:3) %>% mutate(contact_weight = c(10, 5, 2, 5, 7, 1, 2, 1, 4))
      # group_props <-  c(0.3, 0.2, 0.5)
      # # hh_data <-     timing_test$outbreak_t[[1]]$hh_status
      # hh_data <- hh_status %>% filter(i_group %in% 1:3)

      # Check lengths are the same
      if (length(p_contact_if_isolated_home) != length(p_contact_traced) | length(p_contact_traced) != length(contact_means)) {
        print(p_contact_if_isolated_home)
        print(p_contact_traced)
        print(contact_means)
        stop("Lengths of p_contact_if_isolated_home, p_contact_traced, and contact_means are not the same")
      }
      
      n_groups <- length(p_contact_if_isolated_home)

      # Match i groups to their parameters
      params_df <- tibble(
        i_group = 1:n_groups,  
        # r0 = r0,
        p_contact_if_isolated_home = p_contact_if_isolated_home,
        p_contact_traced = p_contact_traced,
        sar_home = sar_home,
        sar_out = sar_out
      )
      
      # Remove the values so that they're not accidentally taken instead of the ones in the dataframe
      rm(p_contact_if_isolated_home, p_contact_traced, sar_home, sar_out)
      
      n <- nrow(new_cases_w_isolation)
      
      # Generate number of potential cases
      # tic("n_secondary")
      n_secondary <- new_cases_w_isolation %>% 
        left_join(params_df, by = "i_group") %>% 
        mutate(sar_home = sar_home * infectiousness,
               sar_out = sar_out * infectiousness) %>% 
        mutate(
          n_contacts_out = n_contacts_out(i_group = i_group, contact_means = contact_means, contact_dispersion = contact_dispersion),
          # n_secondary_cases_home = secondary_cases_n(i_group = i_group, n_contacts = hh_size - 1, sar = sar_home),
          n_secondary_cases_out = secondary_cases_n(i_group = i_group, n_contacts = n_contacts_out, sar = sar_out)
        ) %>% 
        relocate(case_id, i_group, n_contacts_out, n_secondary_cases_out, symptom_timing) 
      # toc()
      
      # Generate IDs for new potential cases
      # tic("secondary ids")
      new_ids <- secondary_ids(
        case_id = n_secondary$case_id, i_group = n_secondary$i_group,
        hh_ind_id = n_secondary$hh_ind_id, hh_id = as.character(n_secondary$hh_id), 
        sar_home = n_secondary$sar_home,
        # n_secondary_cases_home = n_secondary$n_secondary_cases_home,
        n_secondary_cases_out = n_secondary$n_secondary_cases_out, hh_data = hh_data, contact_weights = contact_weights
      ) %>% 
        filter(!is.na(secondary_hh_ind_id))
      # toc()
      
      # Vector of how many new cases for each id
      new_cases_vec <- new_ids %>% group_by(case_id) %>% summarise(n_secondary = n(), .groups = "drop")
      
      # Calculate timings
      new_cases_timings <- n_secondary %>% 
        select(case_id, primary_symptoms_index = symptom_timing_index) %>% 
        inner_join(new_cases_vec, by = "case_id") %>% 
        secondary_case_timing(
          case_df = .,
          timings_dist = timings_dist,
          params_serial = params_serial,
          params_symptom_timing = params_symptom_timing
        )

      # MERGE the secondary cases with the IDs and the timings
      new_secondary_merged <- inner_join(n_secondary, new_ids, by = c("case_id", "hh_ind_id", "hh_id")) %>% 
        bind_cols(new_cases_timings %>% select(-case_id)) %>% 
        # Add infection time to new timing variables
        mutate(
          across(
            c(primary_symptoms_timing, secondary_case_timing, secondary_symptoms_timing),
            ~ infection_timing + .x
          )
        )
      

      # CALCULATE WHO IS QUARANTINED
      # tic("quarantined")
      if (nrow(new_secondary_merged %>% filter(!is.na(secondary_hh_ind_id))) == 0 ) {
        
        # If there are no actual new cases, return NA-filled dfs for primary and secondary quarantined
        new_secondary_transmissions <- new_secondary_merged %>% 
          mutate(primary_quarantined = NA, secondary_quarantined = NA)
  
        return(new_secondary_transmissions)
        
      } else {
        new_secondary_with_secondary_quarantine <- new_secondary_merged %>% 
          calc_quarantine(hh_data = hh_data)
      }
      # toc()
      
      # CALCULATE WHO IS ISOLATED
      # Work out whether transmission is isolated, and if contact happens despite that 
      # tic("new_secondary_transmissions")
      new_secondary_transmissions <- new_secondary_with_secondary_quarantine %>% 
        mutate(
          transmission_isolated = (isolation_timing < secondary_case_timing & (secondary_case_timing <= deisolate_timing | is.na(deisolate_timing))) | 
            isolation_timing_2 < secondary_case_timing | 
            primary_quarantined | secondary_quarantined,
          transmission_isolated = na_to_false(transmission_isolated),
          contact_if_isolated = case_when(
            transmission_isolated & hh_id == secondary_hh_id ~ bern_vec(n = nrow(.), p = p_contact_if_isolated_home),    # at home, some chance of transmission
            transmission_isolated & hh_id != secondary_hh_id ~ FALSE,                                                    # outside infections, no chance of transmission
            TRUE ~ NA
          )
        ) %>% 
        
        # CALCULATE WHO IS CONTACT TRACED
        # Only contact traced if you're detected, and if the row is a new potential case (so secondary case timing exists, and contact still occurs despite isolation)
        mutate(
          secondary_contact_traced = if_else(
            detected & !is.na(secondary_case_timing) & (is.na(contact_if_isolated) | contact_if_isolated == TRUE), # only when there is a real potential case, and the person got tested
            as.logical(mc2d::rbern(nrow(.), p = p_contact_traced)),
            NA
          )
        ) %>% 
        # Calculate CT timings
        mutate(
          secondary_contact_tracing_start_timing = if_else(secondary_contact_traced, pmin(test_result_timing, contact_tracing_results_timing, random_test_timing, na.rm = TRUE), NA_real_),
        ) %>% 
        sample_by_i_group(
          i_group = secondary_i_group,
          sample_dat = ct_delay_data,
          sample_vars = c("secondary_ct_delay" = "ct_delay")
        ) %>% 
        sample_by_i_group(
          i_group = secondary_i_group,
          sample_dat = test_delay_data,
          sample_vars = c("secondary_contact_tracing_results_delay" = "test_results_delay")
        ) %>% 
        mutate(
          secondary_contact_tracing_test_timing = if_else(secondary_contact_traced, secondary_contact_tracing_start_timing + secondary_ct_delay, NA_real_),
          secondary_contact_tracing_results_delay = if_else(secondary_contact_traced, secondary_contact_tracing_results_delay, NA_real_),
          secondary_contact_tracing_results_timing = secondary_contact_tracing_test_timing + secondary_contact_tracing_results_delay
        ) %>% 
        select(-secondary_ct_delay) %>% 
        ungroup
      # toc()
    
      return(new_secondary_transmissions)
      
    }
    
    
    
    # SPEED TEST OF new secondary
    # tic()
    # invisible(draw_secondary_cases(
    #   new_cases_w_isolation = timing_test$outbreak_t[[1]]$live_cases,
    #   sar_home = c(0.2, 0.1, 0.1), 
    #   sar_out =  c(0.05, 0.05, 0.05),
    #   p_contact_if_isolated_home = c(0.3, 0.2, 0.1),
    #   p_contact_traced = c(0.1, 0.2, 0.2),
    #   contact_means = k_matrix_trial %>% rowSums(),
    #   contact_weights = crossing(to = 1:3, from = 1:3) %>% mutate(contact_weight = c(10, 5, 2, 1, 5, 7, 4, 2, 2)),
    #   contact_dispersion = 0.5,
    #   timings_dist = timing_test$outbreak_t[[1]]$params$timings_dist,
    #   params_symptom_timing = timing_test$outbreak_t[[1]]$params$params_symptom_timing,
    #   params_serial = timing_test$outbreak_t[[1]]$params$params_serial,
    #   ct_delay_data = ct_delay_data_real,
    #   test_results_delay_data = data_save$test_results_delay_data %>% filter(i_group %in% 1:3),
    #   group_props = c(0.3, 0.2, 0.2),
    #   hh_data = timing_test$outbreak_t[[1]]$hh_status
    # ))
    # toc()  # starting point 3.6 sec
    # 2.4 sec is due to new_ids
    
  
    
    
   
    

# Functions used when random testing is live ------------------------------

    
    # outbreak_setup_test$secondary_cases
    
    update_isolation_secondary <- function(secondary_cases, live_cases_rerandom, isolation_cols, hh_data, test_delay_data, ct_delay_data) {
      
      # DEBUGGING:
      # secondary_cases <- secondary_cases_dt
      # hh_data <- hh_status_rerandom
      
      # After we've redrawn random testing on the live_cases df, we need to merge the new random testing data
      # back into the secondary cases data and recalculate isolation times.
      # This involves no drawing of new variables, just recalculations
      
      if (is.null(secondary_cases)) return(NULL)
      if (nrow(secondary_cases) == 0) return(secondary_cases)
      
      
      # Update the isolation timings in secondary cases to match the live_cases with new random testing data
      updated_secondary <- right_join(
        live_cases_rerandom %>% select(case_id, all_of(isolation_cols)),
        secondary_cases %>% 
          select(-all_of(isolation_cols), 
                 -any_of(
                   c(
                     "primary_quarantined",
                     "secondary_quarantined",
                     "transmission_isolated",
                     "contact_if_isolated"
                     # "secondary_contact_traced",
                     # "secondary_contact_tracing_start_timing",
                     # "secondary_contact_tracing_test_timing",
                     # "secondary_contact_tracing_results_delay",
                     # "secondary_contact_tracing_results_timing"
                   )
                 )),
        by = "case_id"
      )
      
      # Roughly performs the same actions as draw_secondary_cases but only recalculating rather than redrawing completely new cases
      # If there are no actual new cases, return NA-filled dfs for primary and secondary quarantined
      if (nrow(updated_secondary %>% filter(!is.na(secondary_hh_ind_id))) == 0 ) {
        
        new_secondary_with_secondary_quarantine <- updated_secondary %>%
          mutate(
            primary_quarantined = NA,
            secondary_quarantined = NA,
            # transmission_isolated = NA,
            # contact_if_isolated = NA
          )
        
      } else {
        
        new_secondary_with_secondary_quarantine <- updated_secondary %>% 
          calc_quarantine(hh_data = hh_data)
        
      }
      
      
      # Work out whether transmission is isolated, and if contact happens despite that 
      new_secondary_transmissions <- new_secondary_with_secondary_quarantine %>% 
        mutate(
          transmission_isolated = (isolation_timing < secondary_case_timing & (secondary_case_timing <= deisolate_timing | is.na(deisolate_timing))) | 
            isolation_timing_2 < secondary_case_timing | 
            primary_quarantined | secondary_quarantined,
          transmission_isolated = na_to_false(transmission_isolated),
          contact_if_isolated = case_when(
            transmission_isolated & hh_id == secondary_hh_id ~ bern_vec(n = nrow(.), p = p_contact_if_isolated_home),    # at home, some chance of transmission
            transmission_isolated & hh_id != secondary_hh_id ~ FALSE,                                                    # outside infections, no chance of transmission
            TRUE ~ NA
          )
        )
      
      # # Only contact traced if you're detected, and if the row is a new potential case (so secondary case timing exists, and contact still occurs despite isolation)
        # mutate(
        #   secondary_contact_traced = if_else(
        #     detected & !is.na(secondary_case_timing) & (is.na(contact_if_isolated) | contact_if_isolated == TRUE), # only when there is a real potential case, and the person got tested
        #     as.logical(mc2d::rbern(nrow(.), p = p_contact_traced)),
        #     NA
        #   )
        # ) %>% 
        # 
        # # Calculate CT timings
        # mutate(
        #   secondary_contact_tracing_start_timing = if_else(secondary_contact_traced, pmin(test_result_timing, contact_tracing_results_timing, random_test_timing, na.rm = TRUE), NA_real_),
        # ) %>% 
        # sample_by_i_group(
        #   i_group = secondary_i_group,
        #   sample_dat = ct_delay_data,
        #   sample_vars = c("secondary_ct_delay" = "ct_delay")
        # ) %>% 
        # sample_by_i_group(
        #   i_group = secondary_i_group,
        #   sample_dat = test_results_delay_data,
        #   sample_vars = c("secondary_contact_tracing_results_delay" = "test_results_delay")
        # ) %>% 
        # mutate(
        #   secondary_contact_tracing_test_timing = if_else(secondary_contact_traced, secondary_contact_tracing_start_timing + secondary_ct_delay, NA_real_),
        #   secondary_contact_tracing_results_delay = if_else(secondary_contact_traced, secondary_contact_tracing_results_delay, NA_real_),
        #   secondary_contact_tracing_results_timing = secondary_contact_tracing_test_timing + secondary_contact_tracing_results_delay
        # ) %>% 
        # select(-secondary_ct_delay)
      
      new_secondary_transmissions
      
    }
    
    
    
    
    
    redraw_contact_tracing <- function(secondary_cases_w_isolation, test_delay_data, ct_delay_data) {
      
      if (is.null(secondary_cases_w_isolation)) return(NULL)
      if (nrow(secondary_cases_w_isolation) == 0) return(secondary_cases_w_isolation)
      
      else {
        
        updated_secondary_contact_tracing <- secondary_cases_w_isolation %>% 
          
          # Are cases newly detected, or already detected?
          mutate(
            detected_lag1 = detected,
            detected = (self_test_yn & !self_test_false_negative) | 
              (!is.na(contact_tracing_results_timing) & !ct_test_negative & !ct_false_negative) | 
              (random_test_yn & !random_test_false_negative),
            newly_detected = detected > detected_lag1
          ) %>% 
          
          # If newly_detected, then calculate secondary contact tracing details anew
          mutate(
            secondary_contact_traced = if_else(
              newly_detected & !is.na(secondary_case_timing) & ( is.na(contact_if_isolated) | contact_if_isolated == TRUE ), # only when there is a real potential case, and the person got tested
              as.logical(mc2d::rbern(nrow(.), p = p_contact_traced)),
              secondary_contact_traced
            )
          ) %>% 
          mutate(
            secondary_contact_tracing_start_timing = if_else(
              newly_detected & secondary_contact_traced, 
              pmin(test_result_timing, contact_tracing_results_timing, random_test_timing, na.rm = TRUE), 
              secondary_contact_tracing_start_timing
            )
          ) %>% 
          sample_by_i_group(
            i_group = secondary_i_group,
            sample_dat = ct_delay_data,
            sample_vars = c("secondary_ct_delay" = "ct_delay")
          ) %>% 
          sample_by_i_group(
            i_group = secondary_i_group,
            sample_dat = test_delay_data,
            sample_vars = c("secondary_contact_tracing_results_delay_upd" = "test_results_delay")
          ) %>% 
          mutate(
            secondary_contact_tracing_test_timing = if_else(newly_detected & secondary_contact_traced, secondary_contact_tracing_start_timing + secondary_ct_delay, secondary_contact_tracing_test_timing),
            secondary_contact_tracing_results_delay = if_else(newly_detected & secondary_contact_traced, secondary_contact_tracing_results_delay_upd, secondary_contact_tracing_results_delay),
            secondary_contact_tracing_results_timing = secondary_contact_tracing_test_timing + secondary_contact_tracing_results_delay
          ) %>% 
          select(-secondary_contact_tracing_results_delay_upd, -secondary_ct_delay) %>% 
          
          
          # If was already detected, then check whether new detection time is earlier, and impute if it is
          mutate(
            secondary_contact_tracing_start_timing_lag1 = secondary_contact_tracing_start_timing,
            # Possible new start timing
            secondary_contact_tracing_start_timing_candidate = pmin(test_result_timing, contact_tracing_results_timing, random_test_timing, na.rm = TRUE),
            
            # Impute if the new start timing is earlier, and if newly detected and traced
            new_secondary_contact_tracing_start_timing = 
              secondary_contact_tracing_start_timing_candidate < secondary_contact_tracing_start_timing_lag1 & 
              !newly_detected & secondary_contact_traced,
            
            secondary_contact_tracing_start_timing = if_else(
              new_secondary_contact_tracing_start_timing, 
              secondary_contact_tracing_start_timing_candidate, 
              secondary_contact_tracing_start_timing_lag1
            )
          ) %>% 
          sample_by_i_group(
            i_group = secondary_i_group,
            sample_dat = ct_delay_data,
            sample_vars = c("secondary_ct_delay" = "ct_delay")
          ) %>% 
          mutate(
            # Update the other timing values
            secondary_contact_tracing_test_timing = if_else(new_secondary_contact_tracing_start_timing,
                                                            secondary_contact_tracing_start_timing + secondary_ct_delay,
                                                            secondary_contact_tracing_test_timing)
          ) %>% 
          sample_by_i_group(
            i_group = secondary_i_group,
            sample_dat = test_delay_data,
            sample_vars = c("secondary_contact_tracing_results_delay_upd" = "test_results_delay")
          ) %>% 
          mutate(
            secondary_contact_tracing_results_delay = if_else(new_secondary_contact_tracing_start_timing,
                                                              secondary_contact_tracing_results_delay_upd,
                                                              secondary_contact_tracing_results_delay),
            secondary_contact_tracing_results_timing = secondary_contact_tracing_test_timing + secondary_contact_tracing_results_delay
          ) %>% 
          select(-secondary_contact_tracing_results_delay_upd, -secondary_ct_delay)
        
        updated_secondary_contact_tracing
        
      }
      
    }
    
    
    
    
    # Takes the new contact tracing times from the updated secondary case df and feeds back into the live cases where applicable
    update_live_contact_tracing <- function(live_cases, updated_secondary_cases) {
      
      if (is.null(updated_secondary_cases)) return(live_cases)
      if (nrow(updated_secondary_cases)) return(live_cases)
      
      else {
        
        left_join(
          live_cases,
          updated_secondary_cases %>% select(case_id = secondary_case_id, 
                                             new_secondary_contact_tracing_start_timing,
                                             secondary_contact_tracing_start_timing,
                                             secondary_contact_tracing_test_timing,
                                             # secondary_contact_tracing_results_delay,
                                             secondary_contact_tracing_results_timing),
          by = "case_id"
        ) %>% 

          mutate(
            contact_tracing_test_timing = coalesce(secondary_contact_tracing_test_timing, contact_tracing_test_timing),
            contact_tracing_results_timing = coalesce(secondary_contact_tracing_results_timing, contact_tracing_results_timing)
          ) %>% 
          # If test occurs before infection, the test will return negative
          mutate(
            ct_test_negative = contact_tracing_test_timing < infection_timing
          ) %>% 
          select(-secondary_contact_tracing_start_timing, - secondary_contact_tracing_test_timing, -secondary_contact_tracing_results_timing, -new_secondary_contact_tracing_start_timing)
        
      }
      
    }
    
    

    
# Functions for running an outbreak simulation ----------------------------
    
    # Convenience function to convert the number of initial cases c(N1, N2, N3) into vector of the i_groups at individual level
    initial_cases_to_i_group <- function(n_initial_cases) {
      
      i_group <- c()
      
      for (i in 1:length(n_initial_cases)) {
        i_group <- c(i_group, rep(i, n_initial_cases[[i]]))
      }
      
      i_group
      
    }
    
    
    # OUTBREAK SETUP FUNCTION
    # Starts with a data set of initial cases
    outbreak_setup <- function(
      n_initial_cases, group_props, n_pop_total, dt_approx,
      hh_size_data, 
      sample_hh_data = TRUE,
      print_detail = TRUE,
      test_delay_data,
      # test_choice_delay_data,
      # test_results_delay_data, 
      ct_delay_data,
      recov_val, test_sensitivity_data,   # PARAMETERS HERE
      params_timing,
      probs_self_test_df, probs_isolate_symptoms_df, probs_isolate_test_df, probs_isolate_ct_df,
      params_serial,
      params_symptom_timing,
      p_contact_if_isolated_home, p_contact_traced, p_hh_quarantine,
      k_matrix,
      contact_dispersion,
      alpha, 
      p_immune = c(0, 0, 0, 0),
      sar_home, sar_out, infectiousness_by_symptoms
    ) {
      
      # DEBUGGING:
      # dt_approx <- 1
      # n_initial_cases <- c(5, 5, 5, 5)
      # group_props = data_save$group_props
      # ct_delay_data             = data_save$ct_delay_data
      # probs_self_test_df        = data_save$probs_self_test
      # probs_isolate_symptoms_df = data_save$probs_isolate_symptoms
      # probs_isolate_test_df     = data_save$probs_isolate_test
      # probs_isolate_ct_df       = data_save$probs_isolate_ct
      # p_contact_if_isolated_home= rep(1, 4)
      # p_contact_traced          = data_save$p_contact_traced
      # p_hh_quarantine           = data_save$params_data$p_hh_quarantine
      # hh_size_data              = data_save$hh_data_bogota
      # test_delay_data           = data_save$test_delay_data
      # sar_out                   = data_save$params_data$sar_out_input
      # sar_home                  = data_save$params_data$sar_home_input
      # alpha                     = rep(0, 4)
      # n_pop_total <- 100000
      # k_matrix                  = data_save$k_matrix * (n_pop_total / data_save$k_matrix_pop)
      # params_timing = data_save$params_timing
      # params_symptom_timing = data_save$params_symptom_timing
      # params_serial = data_save$params_serial
      # infectiousness_by_symptoms =  data_save$infectiouness_by_symptoms
      # test_sensitivity_data = data_save$test_sensitivity_data
      # recov_val = 10
      # sample_hh_data = TRUE
      # print_detail = TRUE
      
      n_groups <- length(group_props)
      
      # CHECK PARAMETERS - make sure they are same length as number of groups
      lst(n_initial_cases, p_contact_if_isolated_home, p_contact_traced, 
          unique(test_delay_data$i_group), unique(hh_size_data$i_group), 
          p_immune, p_hh_quarantine, alpha, sar_home, sar_out) %>% 
        imap(~ {
          if (length(.x) != n_groups) {
            print(paste0("n_groups = ", n_groups))
            print(paste0("length of ", .y, " = ", length(.x)))
            stop(paste0("Mismatch between number of groups in ", .y))
          }
        })
      if (sum(group_props) != 1) stop("group_props should sum to 1")
      
      
      if (dt_approx < 1) {
        if ((1/dt_approx) != as.integer(1/dt_approx)) stop("dt_approx needs to divide with no remainders into 1")
      } else if (dt_approx > 1) {
        if (dt_approx != as.integer(dt_approx)) stop("dt_approx needs to divide with no remainders into 1")
      }
      
      # Check that no negative values in testing delays
      if (sum(test_delay_data$test_choice_delay <= 0) > 0) stop("neg values in test_choice_delay")
      if (sum(test_delay_data$test_results_delay <= 0) > 0) stop("neg values in test_results_delay")
      if (sum(ct_delay_data$ct_delay <= 0) > 0) stop("neg values in ct_delay_data")

      # Check that k_matrix is symmetric
      if (!is.matrix(k_matrix) || !isSymmetric(k_matrix)) {
        print(k_matrix)
        stop("k_matrix is not symmetric")
      }
      if (dim(k_matrix)[[1]] != n_groups | dim(k_matrix)[[2]] != n_groups) stop("k_matrix doesn't have correct number of dimensions")
      
      # Calculate average next generation [R0] matrix
      pops <- calc_pops(group_props = group_props, n_pop = n_pop_total)
      contact_matrix_ind <- k_matrix / pops
      
      # Calculate contact_means
      contact_means <- rowSums(contact_matrix_ind)
      
      # Calculate contact_weights from k_matrix
      contact_weights <- k_matrix %>%
        {. / rowSums(.)} %>%     # get the probability of encountering each group (q)
        as.vector() %>% 
        {mutate(crossing(to = 1:n_groups, from = 1:n_groups), q = .)} %>% 
        arrange(from) %>% select(from, to, everything()) %>% 
        left_join(
          tibble(to = 1:n_groups, group_props = group_props), 
          by = "to"
        ) %>% 
        mutate(contact_weight = round(q / group_props, 2))
      
      # Calculate implied beta_matrix
      pop_matrix <- c(pops) %*% t(c(pops))
      beta_matrix <- k_matrix / pop_matrix

      total_initial_cases <- sum(n_initial_cases)
      
      # Calculate R0 value
      
          # Avg. HH size [at individual level]
          mean_hh_size <- hh_size_data %>% 
            group_by(i_group) %>% 
            summarise(mean_hh_size = mean(hh_size), .groups = "drop") %>% 
            pull(mean_hh_size)
          
          # Within home contact matrix
          contact_matrix_ind_home <- diag(mean_hh_size - 1) 
      
      next_gen_inf_out <- contact_matrix_ind * sar_out
      next_gen_inf_home <- contact_matrix_ind_home * sar_home
      next_gen_matrix <- next_gen_inf_out + next_gen_inf_home
      
      r0_groups <- rowSums(next_gen_matrix)
      r0_total <- eigen(next_gen_matrix)$values[[1]]
      
      print(paste0("r0 is approx ", round(r0_total, 2)))
      
      
      # PREPARE THE TIMINGS DISTIRBUTION PARAMETERS
      params_timing <- params_timing[1:7] %>% as.list() %>% set_names(
        c("var_symp", "var_serial", "cov_serial_primary", "cov_serial_secondary",
          "alpha_primary", "alpha_secondary", "alpha_serial")
      )
      
      # Calculate the Omega for the timings distributions
      params_timing$Omega <- matrix(
        c(params_timing$var_symp, 0, params_timing$cov_serial_primary, 
          0,  params_timing$var_symp, params_timing$cov_serial_secondary, 
          params_timing$cov_serial_primary, params_timing$cov_serial_secondary, params_timing$var_serial),
        nrow = 3
      )

      # Calculate distributions (used in secondary timings function, in draw secondary cases)
      joint_dist <- sn::makeSECdistr(
        dp = list(xi = c(0, 0, 0), 
                  Omega = params_timing$Omega, 
                  alpha = c(params_timing$alpha_primary, params_timing$alpha_secondary, params_timing$alpha_serial)),
        family = "SN",
        compNames = c("primary_symp", "secondary_symp", "serial")
      )
      
      marg_dist_primary <- marginalSECdistr(joint_dist, comp = 1)
      marg_dist_secondary <- marginalSECdistr(joint_dist, comp = 2)
      marg_dist_serial <- marginalSECdistr(joint_dist, comp = 3)
      
      timings_dist <- list(
        joint_dist = joint_dist,
        marg_dist_primary = marg_dist_primary,
        marg_dist_secondary = marg_dist_secondary,
        marg_dist_serial = marg_dist_serial,
        min_infection = 1
      )
      

      # Generate live cases (without HH data)
      live_cases_0 <- draw_symptoms_recovery(ids = 1:total_initial_cases,
                                             i_group = initial_cases_to_i_group(n_initial_cases),
                                             hh_id = NULL,
                                             hh_size = NULL,
                                             hh_ind_id = NULL,
                                             would_quarantine = NULL,
                                             infection_timings = 0             ,
                                             symptom_timing_index = NULL,
                                             contact_tracing_test_timing = NULL,
                                             contact_tracing_results_timing = NULL,
                                             # random_test_yn = NULL,
                                             # random_test_timing = NULL,
                                             # random_test_false_negative = NULL,
                                             # isolate_random_test_timing = NULL,
                                             marg_dist_primary = marg_dist_primary,
                                             params_symptom_timing = params_symptom_timing,
                                             test_sensitivity_data = test_sensitivity_data,
                                             test_delay_data = test_delay_data,
                                             # test_choice_delay_data = test_choice_delay_data,
                                             # test_results_delay_data = test_results_delay_data,
                                             recov_val = recov_val,
                                             probs_self_test_df = probs_self_test_df,
                                             probs_isolate_symptoms_df = probs_isolate_symptoms_df,
                                             probs_isolate_test_df = probs_isolate_test_df, 
                                             probs_isolate_ct_df = probs_isolate_ct_df,
                                             infectiousness_by_symptoms = infectiousness_by_symptoms) %>% 
        draw_random_testing(alpha = alpha, test_sensitivity_data = test_sensitivity_data, test_delay_data = test_delay_data) %>% 
        select(-hh_id, -hh_size, -would_quarantine) %>% 
        gen_isolation_all()
      
      # Generate hh data
      if (sample_hh_data) {
        if (print_detail) print("Generating sampled hh data")
        hh_details <- gen_hh_data(n_pop = n_pop_total, group_props = group_props, hh_size_data = hh_size_data, p_hh_quarantine = p_hh_quarantine)
      } else if (!sample_hh_data) {
        if (print_detail) print("Using full hh data")
        hh_details <- hh_size_data %>% 
          # Add in quarantine variables
          left_join(
            tibble(i_group = 1:n_groups,
                   p_hh_quarantine = p_hh_quarantine),
            by = "i_group"
          ) %>% 
          mutate(
            would_quarantine = bern_vec(n = nrow(.), p = p_hh_quarantine)
          )
      }
      
      # Split by group
      hh_details_list <- hh_details %>% 
        group_by(i_group) %>% 
        group_split()
      
      # Split live cases by group
      live_by_i_group <- live_cases_0 %>% 
        group_by(i_group) %>% 
        group_split()
      
      # Randomly draw the HH IND IDs of the starting live cases
      live_with_hh <- map2(.x = live_by_i_group,
           .y = hh_details_list,
           .f = ~ {
             .x %>% mutate(hh_ind_id = sample(x = .y$hh_ind_id, size = nrow(.), replace = FALSE))
           }) %>% 
        bind_rows() %>% 
        left_join(hh_details, by = c("hh_ind_id", "i_group")) %>% 
        relocate(starts_with("hh_"))
    
      # If immunity is not NULL, update the susceptibility status (for modellling vaccines)
      p_immune_df <- tibble(
        i_group = 1:n_groups,
        p_immune = p_immune
      )
      
      hh_details_immune <- hh_details %>% 
        left_join(p_immune_df, by = "i_group") %>% 
        mutate(susceptible = bern_vec(nrow(.), p = (1 - p_immune)))
      
      # Update who is susceptible and isolation timings in the hh status tracker
      hh_status <- hh_details_immune %>% 
        mutate(susceptible = if_else(hh_ind_id %in% live_with_hh$hh_ind_id, FALSE, susceptible)) %>%
        group_by(hh_id) %>% 
        mutate(n_susceptible = sum(susceptible)) %>% 
        ungroup %>% 
        left_join(select(live_with_hh, hh_ind_id, hh_id, matches("isolate_interval")), 
                  by = c("hh_id", "hh_ind_id"))

      # Generate secondary cases
      secondary_cases_0 <- draw_secondary_cases(live_with_hh, 
                                                # old_case = NULL,
                                                p_contact_if_isolated_home = p_contact_if_isolated_home, 
                                                p_contact_traced = p_contact_traced, 
                                                # r0 = r0, dispersion = dispersion,
                                                contact_means = contact_means,
                                                contact_dispersion = contact_dispersion,
                                                sar_home = sar_home, sar_out = sar_out,
                                                test_delay_data = test_delay_data,
                                                timings_dist = timings_dist,
                                                ct_delay_data = ct_delay_data,
                                                contact_weights = contact_weights, 
                                                group_props = group_props,
                                                params_serial = params_serial,
                                                params_symptom_timing = params_symptom_timing,
                                                hh_data = hh_status) %>% 
        filter_non_cases()

      secondary_cases_old <- NULL
      
      # Generate counters
      counters <- tibble(
        t = 0,
        i_group = 1:n_groups,
        n_pop = group_props * n_pop_total,
        n_cases_live = n_initial_cases,
        n_susceptible = n_pop - n_initial_cases,
        n_isolators = 0,
        n_recovered = 0,
        n_cases_cum = n_initial_cases,
        n_detected = 0,
        n_in_testing = 0,
        n_new_cases = 0,
        n_confirmed_cases = 0,
        n_undetected = n_cases_live,
        prop_susceptible = n_susceptible / n_pop
      )
      
      # Generate vector of isolation cols which is used in update_isolation_secondary
      isolation_cols <- infer_isolation_cols()
      
      return(
        list(live_cases = live_with_hh,
             secondary_cases = secondary_cases_0,
             secondary_cases_old = secondary_cases_old,
             recovered_cases = NULL,
             hh_status = hh_status,
             t = 0,
             counters = counters,
             ind_case_counter = NULL,
             n_pop_total = n_pop_total,
             model_end = FALSE,
             params = list(
               n_groups = n_groups,
               group_props = group_props,
               pops = pops,
               dt_approx = dt_approx,
               test_delay_data = test_delay_data,
               # test_choice_delay_data = test_choice_delay_data,
               # test_results_delay_data = test_results_delay_data,
               recov_val = recov_val,
               test_sensitivity_data = test_sensitivity_data,
               ct_delay_data = ct_delay_data,
               timings_dist = timings_dist,
               probs_self_test_df = probs_self_test_df, probs_isolate_symptoms_df = probs_isolate_symptoms_df, 
               probs_isolate_test_df = probs_isolate_test_df, probs_isolate_ct_df = probs_isolate_ct_df,
               p_contact_if_isolated_home = p_contact_if_isolated_home, p_contact_traced = p_contact_traced, 
               p_hh_quarantine = p_hh_quarantine,
               k_matrix = k_matrix,
               contact_matrix_ind = contact_matrix_ind,
               next_gen_inf_out = next_gen_inf_out,
               next_gen_inf_home = next_gen_inf_home,
               next_gen_matrix <- next_gen_inf_out + next_gen_inf_home,
               r0_groups = r0_groups,
               r0_total = r0_total,
               contact_weights = contact_weights, 
               beta_matrix = beta_matrix,
               hh_size_data = hh_size_data, 
               alpha = alpha, 
               contact_means = contact_means, 
               contact_dispersion = contact_dispersion,
               sar_home = sar_home, sar_out = sar_out,
               isolation_cols = isolation_cols,
               infectiousness_by_symptoms = infectiousness_by_symptoms,
               params_timing = params_timing,
               params_serial = params_serial,
               p_immune = p_immune,
               params_symptom_timing = params_symptom_timing
               # PARAMETERS HERE
             ))
      )
      
    }
    

    
    # OUTBREAK STEP function
    outbreak_step <- function(outbreak_t, record_times = FALSE, print_detail = TRUE) {
      
      # TESTING: 
      # outbreak_t <- outbreak_setup_test
      # 2nd loop:  outbreak_t <- out
      # record_times <- FALSE
      
      # Extract all aspects of the outbreak_t object
      live_cases <- outbreak_t$live_cases
      secondary_cases <- outbreak_t$secondary_cases
      secondary_cases_old <- outbreak_t$secondary_cases_old
      recovered_cases <- outbreak_t$recovered_cases
      hh_status <- outbreak_t$hh_status
      ind_case_counter <- outbreak_t$ind_case_counter
      t <- outbreak_t$t
      counters <- outbreak_t$counters
      n_pop_total <- outbreak_t$n_pop_total
      params <- outbreak_t$params
      
      # Recalculate contact_means (only relevant for policy_change)
      
      
      # Calculate average next generation [R0] matrix
      pops <- calc_pops(group_props = params$group_props, n_pop = n_pop_total)
      contact_matrix_ind <- params$k_matrix / pops
      
      # Calculate contact_means
      params$contact_means <- rowSums(contact_matrix_ind)
      # print(round(params$contact_means, 2))
      
      # FOR TESTING:
      # params$dt_approx <- 7
      

      # (0) UPDATE TIMING:
      t <- t + params$dt_approx
      
      live_cases_dt <- update_timing(live_cases, dt = params$dt_approx)
      secondary_cases_dt <- update_timing(secondary_cases, dt = params$dt_approx) %>% 
        filter_non_cases
      
      secondary_cases_old_dt <- update_timing(secondary_cases_old, dt = params$dt_approx)
      recovered_cases_dt <- update_timing(recovered_cases, dt = params$dt_approx)
      hh_status_dt <- update_timing(hh_status, dt = params$dt_approx)
      
      
      # END the model if no more secondary cases
      if (is.null(secondary_cases_dt)) {
        print("model end - no more cases")
        model_end <- TRUE
        out <- list(
          live_cases = NULL,
          secondary_cases = NULL,
          secondary_cases_old = NULL,
          recovered_cases = NULL,
          hh_status = NULL,
          counters = NULL,
          ind_case_counter = ind_case_counter,
          t = t,
          n_pop_total = n_pop_total,
          model_end = model_end,
          params = params
        )
        return(out)
      } else {
        model_end <- FALSE
      }
      
      # (1) RANDOM TESTING - and redraw the new isolation and contact tracing times
      if (record_times) tic(paste0("randomtesting_", t)) 
      
      # Only run if random testing is not all set to 0 AND if t is a whole number (once every WHOLE time step)
      if (sum(params$alpha) == 0 | round(t, 5) != round(t, 0))  {
        live_cases_rerandom <- live_cases_dt
        hh_status_rerandom <- hh_status_dt
        secondary_cases_old_rerandom <- secondary_cases_old_dt
        live_cases_ct_upd <- live_cases_rerandom
        secondary_cases_rerandom <- secondary_cases_dt
        
      } else if (sum(params$alpha) > 0 & round(t, 5) == round(t, 0)) {
        
        # Redraw random testing for old cases
        live_cases_rerandom <- live_cases_dt %>% 
          draw_random_testing(alpha = params$alpha, test_sensitivity_data = params$test_sensitivity_data,
                              test_delay_data = params$test_delay_data) %>%  # PARAMETERS HERE
          gen_isolation_all()
        
        # Update isolation timings in the hh status tracker
        hh_status_rerandom <- hh_status_dt %>% 
          select(-matches("isolate_interval")) %>% 
          left_join(select(live_cases_rerandom, hh_ind_id, hh_id, matches("isolate_interval")), 
                    by = c("hh_id", "hh_ind_id"))
      
        # if (t > 1)  browser()
        
        # And for old cases as well so that we can recalculate their contact tracing times
        # Do this first so contact tracing times are updated in live before updating the secondary_cases df
        secondary_cases_old_rerandom <- secondary_cases_old_dt %>% 
          update_isolation_secondary(live_cases_rerandom = live_cases_rerandom,
                                     isolation_cols = params$isolation_cols,
                                     hh_data = hh_status_rerandom,
                                     test_delay_data = params$test_delay_data,
                                     ct_delay_data = params$ct_delay_data) %>% 
          redraw_contact_tracing(test_delay_data = params$test_delay_data,
                                 ct_delay_data = params$ct_delay_data)
          
        # Input this new contact tracing time back into the live cases [where applicable]
        live_cases_ct_upd <- live_cases_rerandom %>% 
          update_live_contact_tracing(secondary_cases_old_rerandom) %>% 
          gen_isolation_all()
      
        
        # Recalculate isolation times and contact tracing for the secondary cases
        secondary_cases_rerandom <- secondary_cases_dt %>%
          update_isolation_secondary(live_cases_rerandom = live_cases_ct_upd,
                                     isolation_cols = params$isolation_cols,
                                     hh_data = hh_status_rerandom,
                                     test_delay_data = params$test_delay_data,
                                     ct_delay_data = params$ct_delay_data) %>% 
          redraw_contact_tracing(test_delay_data = params$test_delay_data,
                                 ct_delay_data = params$ct_delay_data)
        
        
      } else {
        stop("Something going wrong with alpha condition in the random testing section")
      }
      
      if (record_times) toc(log = TRUE, quiet = TRUE) 
      
      
      # (2) RECOVERIES
      
      # Calculate who recovers
      live_cases_no_recovered <- live_cases_ct_upd %>% 
        mutate(
          # symptomatic = symptom_timing <= 0,
          recovered = recovery_timing <= 0
        )
      
      # Count recoveries
      new_recoveries <- live_cases_no_recovered %>% group_by(i_group) %>% summarise(n_recovered = sum(recovered, na.rm = TRUE), .groups = "drop_last") %>% 
        full_join(tibble(i_group = 1:params$n_groups), by = "i_group") %>%
        mutate(n_recovered = if_else(is.na(n_recovered), 0L, n_recovered)) %>%
        arrange(i_group)
      n_new_recoveries <- sum(new_recoveries$n_recovered, na.rm = TRUE) #%>% print()
      
      if (n_new_recoveries > 0) {
        recovered_cases_upd <- bind_rows(
          recovered_cases_dt,
          live_cases_no_recovered %>% filter(recovered)
        )
        live_cases_no_recovered <- live_cases_no_recovered %>% filter(!recovered)

      } else if (n_new_recoveries == 0) {
        recovered_cases_upd <- recovered_cases_dt
        # Update the counter to 0 so that it adds to the counters_upd below
        # new_recoveries <- tibble(i_group = counters$i_group,
        #                          n_recovered = 0)
      }
      
      
      # (3) INFECTIONS: (S -> I)
      
      if (record_times) tic(paste0("infections_", t))
      
      if (record_times) tic(paste0("findpotentialcases_", t))
  
      # Calculate the secondary cases that could occur in this period
      cases_secondary_all <- secondary_cases_rerandom %>% 
        mutate(new_potential_case = secondary_case_timing <= 0 & secondary_case_timing < recovery_timing) %>% 
        # count_prop(new_potential_case, value(secondary_case_timing), value(recovery_timing)) %>% 
        select(-any_of("susceptible")) %>%
        left_join(hh_status_rerandom %>% select(secondary_hh_id = hh_id,
                                                secondary_i_group = i_group, 
                                                secondary_hh_ind_id = hh_ind_id,
                                                secondary_hh_size = hh_size,
                                                secondary_would_quarantine = would_quarantine,
                                                susceptible),
                  by = c("secondary_hh_ind_id", "secondary_hh_id", "secondary_i_group", "secondary_hh_size", "secondary_would_quarantine")) %>%         
        mutate(
          new_actual_case = if_else(
            new_potential_case & susceptible & (contact_if_isolated == TRUE | is.na(contact_if_isolated)), # get rid of people who don't contact because they are isolated [removes FALSE, keeps TRUE and NA]
            TRUE,
            FALSE
          )
        )
      
      # Calculate infections that are "pipped at the post" - they infected just after someone else already did [within the same dt_approx]
      # Creates a list of the true piped cases
      pipped_cases <- cases_secondary_all %>% 
        filter(new_actual_case) %>% 
        select(case_id, secondary_hh_ind_id, secondary_case_timing) %>% 
        arrange(secondary_hh_ind_id, secondary_case_timing) %>% 
        mutate(pipped = secondary_hh_ind_id == lag(secondary_hh_ind_id)) %>% 
        filter(pipped)
      
      cases_secondary_all <- cases_secondary_all %>% 
        select(-any_of("pipped")) %>% # remove old pipped variable when it exists
        left_join(pipped_cases, by = c("case_id", "secondary_hh_ind_id", "secondary_case_timing")) %>% 
        mutate(pipped = na_to_false(pipped)) %>% 
        # count_prop(pipped) %>% 
        mutate(
          new_actual_case = if_else(pipped, FALSE, new_actual_case)
        )
      
      if (record_times) toc(log = TRUE, quiet = TRUE) 

      # Are there new cases?
      n_new_cases <- sum(cases_secondary_all$new_actual_case, na.rm = TRUE) #%>% print()
      
      if (record_times) tic(paste0("newcases_", t))

      # If there are new cases, draw their details (symptoms, secondary cases)
      if (n_new_cases == 0) {
        
        new_cases <- NULL
        new_secondary_cases <- NULL
        live_cases_w_new <- live_cases_no_recovered
        new_cases_counter <- tibble(i_group = 1:params$n_groups, n_new_cases = 0)
        hh_status_upd <- hh_status_rerandom
        
      } else if (n_new_cases > 0) {
        
        # Secondary cases - new actual cases only
        secondary_cases_for_new <- cases_secondary_all %>% filter(new_actual_case) %>% 
          relocate(case_id, secondary_case_id) #%>% 

        # Optional: print number of infections that are traced in contact tracing
        # print(paste0("New transmissions that are traced: ",
        #              sum(!is.na(secondary_cases_for_new$secondary_contact_tracing_results_timing), na.rm = TRUE),
        #              " / ", nrow(secondary_cases_for_new)))
        
        if (record_times) tic(paste0("newcasesDrawNewCases_", t))
        
        # Draw symptoms and random testing for new cases
        new_cases <- draw_symptoms_recovery(ids = secondary_cases_for_new$secondary_case_id,
                                            hh_id = secondary_cases_for_new$secondary_hh_id,
                                            hh_ind_id = secondary_cases_for_new$secondary_hh_ind_id,
                                            hh_size = secondary_cases_for_new$secondary_hh_size,
                                            i_group = secondary_cases_for_new$secondary_i_group,
                                            would_quarantine = secondary_cases_for_new$secondary_would_quarantine,
                                            infection_timings = secondary_cases_for_new$secondary_case_timing,
                                            symptom_timing_index = secondary_cases_for_new$secondary_symptoms_index,
                                            contact_tracing_test_timing = secondary_cases_for_new$secondary_contact_tracing_test_timing,
                                            contact_tracing_results_timing = secondary_cases_for_new$secondary_contact_tracing_results_timing,
                                            test_delay_data = params$test_delay_data,
                                            # test_choice_delay_data = params$test_choice_delay_data,
                                            # test_results_delay_data = params$test_results_delay_data,
                                            marg_dist_primary = params$timings_dist$marg_dist_primary,
                                            params_symptom_timing = params$params_symptom_timing,
                                            recov_val = params$recov_val,
                                            test_sensitivity_data = params$test_sensitivity_data,
                                            probs_self_test_df = params$probs_self_test_df,
                                            probs_isolate_symptoms_df = params$probs_isolate_symptoms_df,
                                            probs_isolate_test_df = params$probs_isolate_test_df,
                                            probs_isolate_ct_df = params$probs_isolate_ct_df,
                                            infectiousness_by_symptoms = params$infectiousness_by_symptoms) %>%   # PARAMETERS HERE
          draw_random_testing(alpha = params$alpha, test_sensitivity_data = params$test_sensitivity_data, test_delay_data = params$test_delay_data) %>% 
          gen_isolation_all()

        # BIND TOGETHER THE OLD AND NEW CASES
        live_cases_w_new <- bind_rows(
          live_cases_no_recovered, 
          new_cases
        )
      
        # Update hh_status
        hh_status_upd <- hh_status_rerandom %>% 
          mutate(susceptible = if_else(
            hh_ind_id %in% new_cases$hh_ind_id, FALSE, susceptible # update new infections to be not susceptible
          )) %>% 
          # group_by(hh_id) %>%
          # mutate(n_susceptible = sum(susceptible)) %>% 
          # ungroup %>%
          select(-matches("isolate_interval")) %>%                 # Update new isolation intervals
          left_join(select(live_cases_w_new, hh_ind_id, hh_id, matches("isolate_interval")), 
                    by = c("hh_id", "hh_ind_id"))
        
        if (record_times) toc(log = TRUE, quiet = TRUE) 
        
        if (record_times) tic(paste0("newcasesNewCasesCounter_", t))
        
        # Record which group each new infection comes from (to use in the counter)
        new_cases_counter <- new_cases %>% count(i_group) %>% 
          rename(n_new_cases = n) %>% 
          # Make sure all the groups are represented in the counter
          full_join(tibble(i_group = 1:params$n_groups), by = "i_group") %>%
          mutate(n_new_cases = if_else(is.na(n_new_cases), 0L, n_new_cases)) %>%
          arrange(i_group)
        
        if (record_times) toc(log = TRUE, quiet = TRUE) 
        
        if (record_times) tic(paste0("newcasesNewSecondary_", t))
        
        # Draw secondary cases for the new cases
        new_secondary_cases <- new_cases %>% 
          draw_secondary_cases(
            hh_data = hh_status_upd,
            p_contact_if_isolated_home = params$p_contact_if_isolated_home,
            p_contact_traced = params$p_contact_traced,
            timings_dist = params$timings_dist,
            contact_weights = params$contact_weights,
            contact_means = params$contact_means,
            contact_dispersion = params$contact_dispersion,
            group_props = params$group_props,
            test_delay_data = params$test_delay_data,
            ct_delay_data = params$ct_delay_data,
            sar_home = params$sar_home,
            sar_out = params$sar_out,
            params_serial = params$params_serial,
            params_symptom_timing = params$params_symptom_timing
          )  # PARAMETERS HERE

        if (record_times) toc(log = TRUE, quiet = TRUE) 
        
      }
      
      if (record_times) toc(log = TRUE, quiet = TRUE) 
      
      
      # (4) BIND ALL THE NEW CASES TOGETHER WITH OLD
      if (record_times) tic(paste0("bindcases_", t))
      
      # Updated secondary cases 
      secondary_cases_w_new <- cases_secondary_all %>% 
        filter(!new_potential_case) %>%      # remove the cases that would have materialised this period
        bind_rows(new_secondary_cases)
      
      # Updated ("OLD") secondary cases i.e. ones that have already been realised, but may be required to 
      # change the contact_tracing time in the future
      secondary_cases_old_w_new <- bind_rows(
        secondary_cases_old_dt,
        cases_secondary_all %>% filter(new_potential_case)
      )
      
      # OPTIONAL:
      # Make sure no records in both secondary_cases_w_new and secondary_cases_old_w_new
      if (sum(secondary_cases_old_w_new$secondary_case_id %in% secondary_cases_w_new$secondary_case_id) > 0) {
        print("Error - some cases in secondary_cases_old_w_new AND secondary_cases_w_new:")
        print(sum(secondary_cases_old_w_new$secondary_case_id %in% secondary_cases_w_new$secondary_case_id))
        stop()
      }
      
      if (record_times) toc(log = TRUE, quiet = TRUE) # for bindcases
      
      if (record_times) toc(log = TRUE, quiet = TRUE) # for infections
      
      
      
      # (5) UPDATE COUNTERS
      
      if (record_times) tic(paste0("counters_", t)) 
      
      # OPTIONAL: Count the cases generated by each individual
      # ind_case_counter_new <- cases_secondary_all %>% 
      #   filter(new_potential_case) %>%
      #   select(case_id, i_group, secondary_case_id, secondary_i_group, new_potential_case, new_actual_case, infection_timing, transmission_isolated, susceptible, contact_if_isolated) %>% 
      #   # count_prop(new_potential_case, immune, contact_if_isolated, new_actual_case)
      #   mutate(
      #     isolated_no_contact = if_else(!immune, transmission_isolated & !contact_if_isolated, FALSE),
      #     # isolated_but_contact = if_else(!immune, na_to_false(transmission_isolated & contact_if_isolated), FALSE),
      #   ) %>% 
      #   # 3 mutually exclusive and exhaustive categories - new_actual_case, immune, and isolated_no_contact
      #   mutate(
      #     case_type = case_when(new_actual_case ~ "actual_case",
      #                           immune ~ "immune",
      #                           isolated_no_contact ~ "isolated_no_contact")
      #   ) %>% 
      #   select(case_id, i_group, secondary_case_id, secondary_i_group, infection_timing, case_type)
      
      # Count number of people detected at each time
      live_cases_w_detected <- live_cases_w_new %>% 
        mutate(
          detected_t_self = test_result_timing <= 0 & !self_test_false_negative,
          detected_t_ct = contact_tracing_results_timing <= 0 & !ct_test_negative & !ct_false_negative,
          detected_t_random = random_test_yn_ever & random_test_timing <= 0 & !random_test_false_negative,
          detected_at_t = detected_t_self | detected_t_ct | detected_t_random
        ) %>% 
        mutate(
          false_neg_t_self = test_result_timing <= 0 & self_test_false_negative
          # false_neg_t_ct = contact_tracing_results_timing <= 0 & (ct_test_negative | ct_false_negative),
          # false_neg_t_random = random_test_yn_ever & random_test_timing <= 0 & random_test_false_negative,
          # false_neg_at_t = false_neg_t_self | false_neg_t_ct | false_neg_t_random
        ) %>% 
        mutate(
          in_testing_t_self = self_test_yn & symptom_timing <= 0 & test_result_timing > 0,
          in_testing_t_ct = contact_tracing_test_timing <= 0 & contact_tracing_results_timing > 0,
          in_testing_t_random = random_test_yn & random_test_timing > 0,
          in_testing_at_t = if_else(!detected_at_t | is.na(detected_at_t), 
                                    in_testing_t_self | in_testing_t_ct | in_testing_t_random,
                                    FALSE)
        ) %>% 
        mutate(
          undetected_t_self = !self_test_yn | symptom_timing > 0 | (self_test_yn & test_result_timing <= 0 & self_test_false_negative),   # not (yet) detected through self testing
          undetected_t_ct = is.na(contact_tracing_test_timing) | contact_tracing_test_timing > 0 | (contact_tracing_results_timing <= 0 & (ct_test_negative | ct_false_negative)),  # not detected through contact tracing
          undetected_t_random = !random_test_yn_ever | is.na(random_test_yn_ever) | (random_test_timing <= 0 & random_test_false_negative),   # not detected through random testing
          undetected_at_t = undetected_t_self & undetected_t_ct & undetected_t_random
        ) %>% 
        mutate(across(c(detected_at_t, false_neg_t_self, undetected_at_t, in_testing_at_t), na_to_false))
      
      # Sum up all the detected people
      detection_counters <- live_cases_w_detected %>% 
        group_by(i_group) %>% 
        summarise(across(c(detected_at_t, false_neg_t_self, in_testing_at_t, undetected_at_t), sum),
                  n_cases_live = n(),
                  .groups = "drop_last") %>% 
        rename(n_detected = detected_at_t, n_in_testing = in_testing_at_t, n_undetected = undetected_at_t, n_false_neg_self = false_neg_t_self,
        ) %>% 
        # Make sure all the groups are represented in the counter
        full_join(tibble(i_group = 1:params$n_groups), by = "i_group") %>%
        mutate(across(c(n_detected, n_in_testing, n_undetected, n_cases_live, n_false_neg_self), ~ if_else(is.na(.x), 0L, .x))) %>%
        arrange(i_group)
      
      # Check - make sure the number of people in each detection-category sums up to total of live cases
      if (nrow(live_cases_w_detected) != sum(detection_counters$n_detected) + sum(detection_counters$n_in_testing) + sum(detection_counters$n_undetected)) {
        warning(paste0("detected counter is going wrong : n_cases = ", nrow(live_cases_w_detected), " and cases counted in detected is ", 
                       sum(detection_counters$n_detected) + sum(detection_counters$n_in_testing) + sum(detection_counters$n_undetected)))
        live_cases_w_detected %>% count_prop(
          detected_at_t,
          in_testing_at_t,
          undetected_at_t
        )
        # browser()
      }
      
      
      # Count number of confirmed cases at each time
      # (different to above - includes people who have already recovered by the time they receive results)
      confirmed_cases <- bind_rows(live_cases_w_new, recovered_cases_upd) %>% 
        mutate(
          confirmed_case_self = test_result_timing <= 0 & !self_test_false_negative,
          confirmed_case_ct = contact_tracing_results_timing <= 0 & !ct_test_negative & !ct_false_negative,
          confirmed_case_random = random_test_yn_ever & random_test_timing <= 0 & !random_test_false_negative,
          confirmed_case_t = na_to_false(confirmed_case_self | confirmed_case_ct | confirmed_case_random)
        ) %>% 
        group_by(i_group) %>% 
        summarise(confirmed_cases_self = as.integer(sum_na(confirmed_case_self)),
                  confirmed_cases_ct = as.integer(sum_na(confirmed_case_ct)),
                  confirmed_cases = as.integer(sum_na(confirmed_case_t)), 
                  all_self_test_yn = as.integer(sum_na(self_test_yn)),
                  eventually_confirmed_self = as.integer(sum_na(self_test_yn & !self_test_false_negative)),
                  eventually_false_neg_self = as.integer(sum_na(self_test_yn & self_test_false_negative)),
                  .groups = "drop") %>% 
        # Make sure all the groups are represented in the counter
        full_join(tibble(i_group = 1:params$n_groups), by = "i_group") %>%
        mutate(across(c(confirmed_cases_self, confirmed_cases_ct, confirmed_cases), ~ if_else(is.na(.x), 0L, .x))) %>%
        arrange(i_group)
      
      # Count number of ISOLATING PEOPLE
      hh_isolation <- hh_status_upd %>% 
        mutate(currently_isolating_1 = na_to_false(isolate_interval_1_l < 0 & isolate_interval_1_r > 0),
               currently_isolating_2 = na_to_false(isolate_interval_2_l < 0 & isolate_interval_2_r > 0),
               currently_isolating = currently_isolating_1 | currently_isolating_2) %>% 
        select(hh_id, i_group, hh_size, currently_isolating, would_quarantine)
      
      # Types of isolators, either single isolators, hh_quarantinors
      n_single_isolators <- hh_isolation %>% filter(!would_quarantine) %>% group_by(i_group) %>% summarise(n_single_isolators = sum(currently_isolating), .groups = "drop")
      n_hh_quarantinors <- hh_isolation %>% filter(would_quarantine) %>% filter(currently_isolating) %>% dups::dups_drop(hh_id, warn = FALSE) %>% group_by(i_group) %>% summarise(n_hh_quarantinors = sum(hh_size), .groups = "drop")
      n_isolators <- full_join(n_single_isolators, n_hh_quarantinors, by = "i_group") %>% 
        mutate(n_isolators = n_single_isolators + n_hh_quarantinors) %>% 
        full_join(tibble(i_group = 1:params$n_groups), by = "i_group") %>% # complete any missing groups
        mutate(n_isolators = if_else(is.na(n_isolators), 0, n_isolators)) %>%
        arrange(i_group) %>% select(i_group, n_isolators)
      
      # UPDATE ALL COUNTERS
      counters_upd <- tibble(
        t = counters$t + params$dt_approx,
        i_group = counters$i_group,
        n_new_cases = new_cases_counter$n_new_cases,
        n_susceptible = counters$n_susceptible - new_cases_counter$n_new_cases,
        n_recovered = counters$n_recovered + new_recoveries$n_recovered,
        n_cases_cum = counters$n_cases_cum + new_cases_counter$n_new_cases,
        n_isolators = n_isolators$n_isolators,
        n_confirmed_cases = confirmed_cases$confirmed_cases,
        all_self_test_yn = confirmed_cases$all_self_test_yn,
        eventually_confirmed_self = confirmed_cases$eventually_confirmed_self,
        eventually_false_neg_self = confirmed_cases$eventually_false_neg_self
      ) %>%
        left_join(detection_counters, by = "i_group") %>% 
        mutate(n_pop = n_cases_live + n_susceptible + n_recovered,
               prop_susceptible = n_susceptible / n_pop) %>% 
        relocate(t, i_group, n_pop, n_cases_live, n_susceptible, n_recovered, n_cases_cum, n_detected, n_in_testing, n_undetected, prop_susceptible)
      
      # Check that population size hasn't changed which would indicate an error
      if (sum(counters_upd$n_pop != counters$n_pop) > 0) {
        print(counters)
        print(counters_upd)
        stop("Error: population sizes have changed")
      }
      
      if (record_times) toc(log = TRUE, quiet = TRUE) 
      
      # Print progress to console
      if (print_detail) print(paste0("t = ", t, ", Live cases = ", nrow(live_cases_w_detected)))
    
      # Output the cases and counters
      out <- list(
        live_cases = live_cases_w_new,
        secondary_cases = secondary_cases_w_new,
        secondary_cases_old = secondary_cases_old_w_new,
        recovered_cases = recovered_cases_upd,
        hh_status = hh_status_upd,
        counters = counters_upd,
        confirmed_cases = confirmed_cases,
        t = t,
        n_pop_total = n_pop_total,
        model_end = model_end,
        params = params
      )
      
      return(out)
      
    }
    
    
    
    # Run a single simulation 
    outbreak_sim <- function(n_iterations = 100, keep_all_data = FALSE, record_times = FALSE, print_detail = TRUE, ...) {
      
      # TESTING: n_initial_cases = 100; n_pop = 10000; dt_approx = 1; n_iterations = 100
      
      # Set up using outbreak_setup
      outbreak_t <- outbreak_setup(...)
      
      # Record time series, and possibly all data
      cases_time_series <- list(
        outbreak_t$counters
      )
      if (keep_all_data) outbreak_t_record <- list(outbreak_t)
      
      # LOOP ROUND n_iterations
      for (i in 1:n_iterations) {
        
        outbreak_t <- outbreak_step(outbreak_t, record_times = record_times, print_detail = print_detail) # iterate the model
        cases_time_series[[i + 1]] <- outbreak_t$counters   # record time series
        if (keep_all_data) outbreak_t_record[[i + 1]] <- outbreak_t # record all data if applicable
        
        if (outbreak_t$model_end == TRUE) break # stop if model_end is TRUE
        
      }
      
      # Output data
      out <- list(
        outbreak_t = outbreak_t,
        time_series = cases_time_series %>% bind_rows()
      )
      if (keep_all_data) out$outbreak_t_record = outbreak_t_record
      
      return(out)
      
    }

    # RUN MULTIPLE SIMULATIONS
    outbreak_sims <- function(n_sims, parallel = TRUE, keep_all_data = FALSE, record_times = FALSE, print_detail = TRUE, ...) {
      loop_list <- list()
      
      for (i in 1:n_sims) {
        print(paste0("Running simulation ", i, " of ", n_sims))
        loop_list[[i]] <- outbreak_sim(keep_all_data = keep_all_data, record_times = record_times, print_detail = print_detail, ...)
      }
      
      outbreak_t <- loop_list %>% map("outbreak_t")
      
      time_series <- loop_list %>% map("time_series") %>%
        bind_rows(.id = "sim_id")
      
      # out <- list(
      #   outbreak_t = outbreak_t,
      #   params = outbreak_t[[1]]$params,
      #   time_series = time_series
      # )
      
      # class(out) <- c("outbreak_sims")
      
      # if (keep_all_data == TRUE) {
      #   out$outbreak_t_record <- loop_list %>% map("outbreak_t_record")
      # }
      
      out <- time_series
      
      return(out)
      
    }
  
    
    
    outbreak_sims_parallel <- function(n_sims, n_cores, keep_all_data = FALSE, record_times = FALSE, print_detail = TRUE, ...) {
      
      loop_list <- pbmclapply(1:n_sims, mc.cores = n_cores, function(i) {
        print(paste0("Running simulation ", i, " of ", n_sims))
        out <- outbreak_sim(keep_all_data = keep_all_data, record_times = record_times, print_detail = print_detail, ...)
        return(out$time_series)
      })
      
      # outbreak_t <- loop_list %>% map("outbreak_t")
      # browser()
      time_series <- loop_list %>% 
        set_names(1:n_sims) %>% 
        bind_rows(.id = "sim_id")
      
      out <- time_series
      
      # class(out) <- c("outbreak_sims")
      
      # if (keep_all_data == TRUE) {
      #   out$outbreak_t_record <- loop_list %>% map("outbreak_t_record")
      # }
      
      return(out)
      
    }
    

    
    # Function for testing the effect of different parameters
    outbreak_sims_by_params <- function(...) {
      crossing(...) %>% 
        mutate(epidemic_results = pmap(., outbreak_sims)) %>%
        "class<-"(c("outbreak_sims_by_params", "tbl_df", "tbl", "data.frame"))
    }
    
    
    
    
    
    
    
    
    
    # Run a simulation in which the parameters can change over time
    # time_varying_params becomes a (named) list for each set of parameters
    outbreak_sim_policy_change <- function(n_iterations, keep_all_data = FALSE, record_times = FALSE, print_detail = TRUE, 
                                           constant_params,
                                           time_varying_params,
                                           policy_start_times) {
      # START UP
      outbreak_t <- exec(outbreak_setup, !!!c(constant_params, time_varying_params[[1]]))
      cases_time_series <- list(
        outbreak_t$counters %>% mutate(policy = names(time_varying_params[1]), .after = t)
      )
      
      if (keep_all_data) outbreak_t_record <- list(outbreak_t)
      
      # LOOP ROUND n_iterations
      for (i in 1:n_iterations) {

        t <- outbreak_t$t + 1
        policy_set_index <- max(which(t >= policy_start_times)) # which policy parameter list to choose from?
        policy_set <- time_varying_params[[policy_set_index]]
        policy_set_name <- names(time_varying_params[policy_set_index])
        
        # Update parameters to reflect the policy
        outbreak_t$params <- outbreak_t$params %>% list_modify(!!!policy_set)
        
        outbreak_t <- outbreak_step(outbreak_t, record_times = record_times, print_detail = print_detail)
        
        # Add policy name to counter
        if (!is.null(outbreak_t$counters)) outbreak_t$counters <- outbreak_t$counters %>% mutate(policy = policy_set_name, .after = t)

        # Add counter to time series
        cases_time_series[[i + 1]] <- outbreak_t$counters
        
        if (keep_all_data) outbreak_t_record[[i + 1]] <- outbreak_t
        
        if (outbreak_t$model_end == TRUE) break
        
      }
      
      out <- list(
        outbreak_t = outbreak_t,
        time_series = cases_time_series %>% bind_rows()
      )
      
      if (keep_all_data) out$outbreak_t_record = outbreak_t_record
      
      return(out)
      
    }

    
    
    
    # RUN MULTIPLE SIMULATIONS
    outbreak_sims_policy_change <- function(n_sims, keep_all_data = FALSE, print_detail = TRUE, record_times = FALSE, ...) {
      loop_list <- list()
      
      for (i in 1:n_sims) {
        print(paste0("Running simulation ", i, " of ", n_sims))
        loop_list[[i]] <- outbreak_sim_policy_change(keep_all_data = keep_all_data, record_times = record_times, print_detail = print_detail, ...)
      }
      
      outbreak_t <- loop_list %>% map("outbreak_t")
      
      time_series <- loop_list %>% map("time_series") %>%
        bind_rows(.id = "sim_id")
      
      out <- list(
        outbreak_t = outbreak_t,
        # params = outbreak_t[[1]]$params,
        time_series = time_series
      )
      
      class(out) <- c("outbreak_sims")
      
      if (keep_all_data == TRUE) {
        out$outbreak_t_record <- loop_list %>% map("outbreak_t_record")
      }
      
      return(out)
      
    }
    
    # RUN MULTIPLE SIMULATIONS
    outbreak_sims_policy_change_parallel <- function(n_sims, n_cores, keep_all_data = FALSE, print_detail = TRUE, record_times = FALSE, ...) {
      loop_list <- list()
      
      loop_list <- pbmclapply(1:n_sims, mc.cores = n_cores, function(i) {
        print(paste0("Running simulation ", i, " of ", n_sims))
        out <- outbreak_sim_policy_change(keep_all_data = keep_all_data, record_times = record_times, print_detail = print_detail, ...)
        return(out$time_series)
      })
      
      # outbreak_t <- loop_list %>% map("outbreak_t")
      
      time_series <- loop_list %>% 
        set_names(1:n_sims) %>% 
        bind_rows(.id = "sim_id")
      
      out <- time_series
      
      # class(out) <- c("outbreak_sims")
      
      # if (keep_all_data == TRUE) {
      #   out$outbreak_t_record <- loop_list %>% map("outbreak_t_record")
      # }
      
      return(out)
      
    }
    
    
    
    
    # library("parallel")
 # n_cores <- detectCores()
    
    
    
    
    
    
    
# Functions for plotting --------------------------------------------------
    
    summarise_quantiles <- function(.data, ..., q = c(0.01, 0.025, 0.05, 0.25, 0.5, 0.75, 0.95, 0.975, 0.99)) {
      .data %>% 
        summarise(
          across(..., ~ quantile(.x, probs = q, na.rm = TRUE)),
          q = q
        ) %>% 
        relocate(q)
    }
    
    clean_up_group_names <- function(x) {
      if (is_double(x)) factor(x)
      else if (is_list(x)) map_chr(x, paste0, collapse = "_")
    }
    
    quantile_summarise <- function(data, ..., conf_level = 0.95) {
      if (is.null(conf_level)) stop("Please specify a conf_level")
      
      conf_lower <- (1 - conf_level) / 2
      conf_upper <- 1 - conf_lower
      
      data %>% 
        summarise(
          across(...,
                 list(median = ~ median(.x, na.rm = TRUE),
                      mean = ~ mean(.x, na.rm = TRUE),
                      upper = ~ quantile(.x, conf_upper, na.rm = TRUE),
                      lower = ~ quantile(.x, conf_lower, na.rm = TRUE))),
          .groups = "drop_last"
        )
    }
    
    
    plot.outbreak_sims <- function(data, indiv_sims = FALSE, conf_level = 0.95) {
      
      if (!indiv_sims) {
        data %>% 
          .$time_series %>% 
          group_by(sim_id, t) %>% 
          summarise(across(starts_with("n_"), sum), .groups = "drop_last") %>% 
          # mutate(new_cases = n_cases_cum - lag(n_cases_cum),
          #        new_cases_week = rollsum(new_cases, 7, align = "right", na.pad = TRUE)) %>%
          # group_by(sim_id) %>%
          group_by(t) %>% 
          select(t, sim_id, n_cases_live) %>% 
          filter(!is.na(n_cases_live)) %>% 
          # summarise(
          #   across(n_cases_live, 
          #          list(median = ~ median(.x),
          #               upper_95 = ~ quantile(.x, 0.975),
          #               lower_95 = ~ quantile(.x, 0.025)))
          # ) %>% 
          quantile_summarise(n_cases_live, conf_level = conf_level) %>% 
          ungroup %>% 
          ggplot(aes(x = t)) + 
          geom_line(aes(y = n_cases_live_median), size = 1.3) + 
          geom_ribbon(aes(ymax = n_cases_live_upper, ymin = n_cases_live_lower), alpha = 0.3) + 
          theme_custom() + theme(legend.position = "top") + 
          labs(y = "Live cases")
          
        
      } else if (indiv_sims) {
        data$time_series %>% 
          group_by(sim_id, t) %>% 
          summarise(across(starts_with("n_"), sum), .groups = "drop_last") %>% 
          group_by(sim_id) %>% 
          # filter(t < 100) %>%
          # mutate(new_cases = n_cases_cum - lag(n_cases_cum),
          # new_cases_week = rollsum(new_cases, 7, align = "right", na.pad = TRUE)) %>% 
          ggplot(aes(x = t, y = n_cases_live, colour = sim_id)) + 
          geom_line(alpha = 0.6, show.legend = FALSE) + 
          theme_custom() + 
          labs(y = "Live cases")
      }
      
      
    }
    
    is_named_list <- function(x) {
      is_list(x) & !is.null(names(x))
    }
    
    plot.outbreak_sims_by_params <- function(data, param_groups, conf_level = 0.95, indiv_sims = FALSE) {
      
      # if (missing(param_group)) stop("param_group is missing")
      param_groups_syms <- syms(param_groups)
      
      data_time_series <- data %>% 
        mutate(time_series = map(epidemic_results, "time_series")) %>% 
        select(-n_pop) %>% 
        unnest(time_series) %>% 
        # Extract named parameters where applicable
        mutate(across(where(is_named_list), 
                      ~ factor(names(.)))) %>% 
        mutate(param_group = paste(!!!param_groups_syms, sep = "_"))
      
      if (!indiv_sims) {
        
        data_time_series %>%
          # mutate("{{param_group}}" := factor({{param_group}})) %>%
          # mutate("{{param_group}}" := map_chr({{param_group}}, str_c, collapse = "_")) %>%
          group_by(param_group, sim_id, t) %>% 
          summarise(across(starts_with("n_") & !where(is.list), sum), .groups = "drop_last") %>% 
          group_by(param_group, sim_id) %>%
          # mutate(new_cases = n_cases_cum - lag(n_cases_cum),
          #        new_cases_week = rollsum(new_cases, 7, align = "right", na.pad = TRUE)) %>%
          group_by(param_group, t) %>%
          select(t, sim_id, n_cases_live) %>%
          filter(!is.na(n_cases_live)) %>%
          quantile_summarise(n_cases_live, conf_level = conf_level) %>% 
          # summarise(
          #   across(n_cases_live,
          #          list(median = ~ median(.x),
          #               upper = ~ quantile(.x, conf_upper),
          #               lower = ~ quantile(.x, conf_lower)))
          # ) %>%
          ungroup %>%
          ggplot(aes(x = t, group = param_group)) +
          geom_line(aes(y = n_cases_live_median,
                        colour = param_group), size = 1.3) +
          geom_ribbon(aes(ymax = n_cases_live_upper, ymin = n_cases_live_lower,
                          fill = param_group), alpha = 0.3) +
          theme_custom() + theme(legend.position = "top")
      } else if (indiv_sims) {
        
        data_time_series %>%
          group_by(param_group, sim_id, t) %>% 
          summarise(across(starts_with("n_") & !where(is.list), sum), .groups = "drop_last") %>% 
          mutate(sim_id_group = paste0(sim_id, param_group)) %>%
          ungroup %>%
          ggplot(aes(x = t, y = n_cases_live, colour = param_group, group = sim_id_group)) +
          geom_line(alpha = 0.6, show.legend = FALSE) + 
          theme_custom()
        
      }
    }
    
    # plot(data = speedy_test_params, param_groups = c("p_contact_traced"), indiv_sims = FALSE)
    
    
    
    
    
    plot_detected <- function(data, param_groups, prop = FALSE) {
      UseMethod("plot_detected")
    }
    
    plot_detected.outbreak_sims <- function(data, prop = FALSE) {
      # data <- speedy_test
      
      # BUG 
      
      print("Using method for outbreak_sims")
      
      data_detected <- data$time_series %>% 
        # mutate(
        #   across(c(n_detected, n_in_testing, n_undetected), list(prop = ~ . / n_cases_live))
        # ) %>%
        group_by(sim_id, t) %>% 
        summarise(across(starts_with("n_"), sum)) %>% 
        ungroup %>% 
        pivot_longer(c(n_undetected, n_in_testing, n_detected)) %>% 
        mutate(name = fct_rev(name)) %>% 
        group_by(t, name) %>% 
        quantile_summarise(value, conf_level = 0.95)
      
      if (prop) {
        data_detected <- data_detected %>% group_by(t) %>% 
          mutate(across(-c(name), ~ . / sum(.))) %>% 
          ungroup()
      }
      
      
      ggplot(data_detected, aes(x = t, y = value_median, fill = name)) + 
        geom_area()
      
    }
    
    plot_detected.outbreak_sims_by_params <- function(data, param_groups, prop = FALSE) {
      
      print("Using method for outbreak_sims_by_params")
      
      param_groups_syms <- syms(param_groups)
      
      data_time_series <- data %>% 
        mutate(time_series = map(epidemic_results, "time_series")) %>% 
        select(-n_pop) %>% 
        unnest(time_series) %>% 
        # Extract named parameters where applicable
        mutate(across(where(is_named_list), 
                      ~ factor(names(.)))) %>% 
        mutate(param_group = paste(!!!param_groups_syms, sep = "_")) %>% 
        group_by(sim_id, t, param_group, i_group) %>% 
        # print
        mutate(across(starts_with("n_") & where(is.numeric), sum)) %>% 
        ungroup
      
      data_detected <- data_time_series %>% 
        pivot_longer(c(n_undetected, n_in_testing, n_detected)) %>% 
        mutate(name = fct_rev(name)) %>% 
        group_by(t, param_group, name) %>% 
        quantile_summarise(value, conf_level = 0.95)
      
      if (prop) {
        data_detected <- data_detected %>% group_by(param_group, t) %>% 
          mutate(across(value_median, ~ . / sum(.))) %>% 
          ungroup()
      }
      
      ggplot(data_detected, aes(x = t, y = value_median, fill = name)) + 
        geom_area() + 
        facet_wrap(~ param_group)
      
    }
    
    
    
    plot_by_i_group <- function(data, stacked = FALSE, prop = FALSE, conf_level = 0.95) {
      
      data_summ <- data %>% .$time_series %>% 
        group_by(t, i_group) %>% 
        quantile_summarise(c(n_cases_live, n_detected), conf_level = conf_level) %>% 
        ungroup
      
      if (!prop && !stacked) {
        ggplot(data_summ, aes(x = t)) + 
          geom_line(aes(y = n_cases_live_median, colour = factor(i_group)), size = 1.3) + 
          geom_ribbon(aes(x = t, ymax = n_cases_live_upper, ymin = n_cases_live_lower, fill = factor(i_group)), alpha = 0.3) + 
          theme_custom() + theme(legend.position = "top") + labs(colour = "SES Group", fill = "SES Group", y = "Live cases")
      } else {
        
        if (prop) {
          data_summ <- data_summ %>% group_by(t) %>% 
            mutate(n_cases_live_median = n_cases_live_median / sum(n_cases_live_median)) %>% 
            ungroup()
        } 
        
        ggplot(data_summ, aes(x = t, y = n_cases_live_median, fill = factor(i_group))) + 
          geom_area()
        
      }
      
      
      
    }
    
    
    
    
    plot_incidence <- function(data) {
      data$time_series %>% 
        group_by(sim_id, t) %>% 
        summarise(weekly_new_cases = sum(n_cases_cum),
                  weekly_new_detected = sum(n_detected), .groups = "drop") %>% 
        mutate(across(starts_with("weekly_"), ~ .x - lag(.x, 7))) %>% 
        group_by(t) %>% 
        # drop_na_count() %>% 
        drop_na() %>% 
        quantile_summarise(c(weekly_new_cases, weekly_new_detected), conf_level = 0.8) %>% 
        
        ggplot(aes(x = t)) + 
        geom_line(aes(y = weekly_new_cases_median)) + 
        geom_line(aes(y = weekly_new_detected_median), linetype = "dashed")      
      
    }
    
    
    
    

      
    
    
    
    
    
    plot_by_i_group_policymaker_mode <- function(data, conf_level = 0.8) {
      data$time_series %>% 
        group_by(t, i_group) %>% 
        quantile_summarise(c(n_cases_live, n_detected), conf_level = conf_level) %>% 
        ungroup %>% 
        ggplot(aes(x = t)) + 
        # geom_line(aes(y = n_cases_live_median, colour = factor(i_group)), size = 1.3) +
        # geom_ribbon(aes(x = t, ymax = n_cases_live_upper, ymin = n_cases_live_lower, fill = factor(i_group)), alpha = 0.3) +
        geom_line(aes(y = n_detected_median, colour = factor(i_group)), linetype = "dashed", size = 1.3) + 
        geom_ribbon(aes(x = t, ymax = n_detected_upper, ymin = n_detected_lower, fill = factor(i_group)), alpha = 0.3) + 
        theme_custom() + theme(legend.position = "top") + labs(colour = "SES Group", fill = "SES Group", y = "Detected live cases")
    }
    
    
    plot_detected_by_i_group <- function(data, detected_only = FALSE, prop = TRUE) {
      
      data_detected <- data$time_series %>% 
        pivot_longer(c(n_undetected, n_in_testing, n_detected)) %>% 
        mutate(name = fct_rev(name)) %>% 
        group_by(t, i_group, name) %>% 
        quantile_summarise(value, conf_level = 0.95)
      
      if (prop) {
        data_detected <- data_detected %>% group_by(i_group, t) %>% 
          mutate(across(-c(name), ~ . / sum(.))) %>% 
          ungroup()
      }
      
      
      data_detected %>% 
        mutate(i_group = paste0("Group ", i_group)) %>% 
        ggplot(aes(x = t, y = value_median, fill = name, group = name)) + 
        geom_area() + 
        geom_line(position="stack", colour = "white") + 
        facet_wrap(~ i_group) + 
        theme_custom() + 
        labs(y = "Proportion of infected") + 
        scale_y_continuous(labels = scales::percent_format(accuracy = 1)) +  # for percentage
        scale_fill_discrete(name = element_blank(), labels = c("Undetected", "In testing", "Detected"))
      
    }
    
    
    # FUNCTION FOR REPORTING how many cases get confirmed
    report_confirmed_ratio <- function(time_series) {
      data_end <- time_series %>% 
        select(sim_id, t, i_group, n_cases_cum, n_confirmed_cases) %>% 
        filter(t == last(t))
      
      by_group <- data_end %>% 
        mutate(confirmed_ratio = n_confirmed_cases / n_cases_cum) %>% 
        group_by(i_group) %>% 
        quantile_summarise(confirmed_ratio) %>% 
        mutate(i_group = as.character(i_group))
      
      overall <- data_end %>% 
        group_by(sim_id) %>% 
        summarise(across(c(n_cases_cum, n_confirmed_cases), sum)) %>% 
        mutate(confirmed_ratio = n_confirmed_cases / n_cases_cum) %>% 
        quantile_summarise(confirmed_ratio) %>% 
        mutate(i_group = "overall")
      
      bind_rows(by_group, overall)
    }
    
    
    # FUNCTION FOR REPORTING HOW MANY NEGATIVE TESTS (to adjust with)
    report_p_negative_result <- function(data) {
      data %>% 
        group_by(sim_id, t) %>% 
        # sum up across i_groups
        summarise(across(-c(i_group, prop_susceptible, policy), sum)) %>% 
        mutate(prop_confirmed = eventually_confirmed_self / (eventually_false_neg_self + eventually_confirmed_self)) %>% 
        select(-c(n_isolators:n_confirmed_cases), n_cases_cum, n_confirmed_cases) %>% 
        mutate(prop_self_test_yn = all_self_test_yn / n_cases_cum) %>% 
        filter(t == last(t)) %>% 
        .$prop_confirmed %>% 
        {mean(1 - .)}
    }
    
    
    
    # data <- time_series_3
    # by_group <- TRUE
    
    plot_new_cases_pc <- function(data, by_group = TRUE) {
      
      data <- data %>% mutate(i_group = recode_i_group(i_group))
      
      # BY GROUP
      if (by_group) {
        data_by_group <- data %>% 
          select(sim_id, t, i_group, n_pop, n_cases_cum) %>% 
          arrange(sim_id, i_group, t) %>% 
          group_by(sim_id, i_group, n_pop) %>% 
          mutate(new_cases_week = n_cases_cum - lag(n_cases_cum, 7),
                 new_cases_week_pc = new_cases_week / n_pop) %>% 
          
          # Amalgamate across sims
          filter(!is.na(new_cases_week_pc)) %>% 
          group_by(t, i_group) %>% 
          quantile_summarise(new_cases_week_pc)
        
        
        p <- ggplot(data_by_group, aes(x = t)) + 
          geom_line(aes(y = new_cases_week_pc_median, colour = i_group), size = 1) + 
          # facet_wrap(~ name) + 
          geom_ribbon(aes(ymin = new_cases_week_pc_lower, 
                          ymax = new_cases_week_pc_upper,
                          fill = i_group),
                      alpha = 0.2) + 
          scale_colour_viridis(option = "viridis", begin = 0.3, discrete = TRUE) + 
          scale_fill_viridis(option = "viridis", begin = 0.3, discrete = TRUE) + 
          theme_custom() + 
          labs(colour = "SES Group")
      
      } else if (!by_group) {
        
        # OVERALL
        data_overall <- data %>% 
          select(sim_id, n_pop, t, i_group, n_cases_cum) %>% 
          # sum over groups
          group_by(sim_id, t) %>% 
          summarise(n_cases_cum = sum(n_cases_cum),
                    n_pop = sum(n_pop)) %>% 
          mutate(new_cases_week = n_cases_cum - lag(n_cases_cum, 7),
                 new_cases_week_pc = new_cases_week / n_pop) %>% 
          
          # Amalgamate across sims
          filter(!is.na(new_cases_week_pc)) %>% 
          group_by(t) %>% 
          quantile_summarise(new_cases_week_pc)
        
        
        p <- ggplot(data_overall, aes(x = t)) + 
          geom_line(aes(y = new_cases_week_pc_median), size = 1) + 
          # facet_wrap(~ name) + 
          geom_ribbon(aes(ymin = new_cases_week_pc_lower, 
                          ymax = new_cases_week_pc_upper),
                      alpha = 0.2) + 
          # scale_colour_viridis(option = "viridis", begin = 0.3, discrete = TRUE) + 
          # scale_fill_viridis(option = "viridis", begin = 0.3, discrete = TRUE) + 
          theme_custom() + 
          labs(colour = "SES Group")
        
        
      }
      
      return(p)
      
      
      
    }
    
    
    
    
    
    
    
    
    
    
    print.outbreak_sims <- function(data) {
      print(str(data, max.level = 2))
      print(data$time_series)
    }
    
    graph_remove_bottom_left <- function(plot) {
      require(grid)
      require(gridExtra)
      g <- ggplotGrob(plot)
      # get the grobs that must be removed
      # rm_grobs <- g$layout$name %in% c("panel-2-1", "panel-3-3", "strip-t-3-1", "strip-t-3-3")
      rm_grobs <- g$layout$name %in% c("panel-2-1", "axis-b-1")
      # remove grobs
      g$grobs[rm_grobs] <- NULL
      g$layout <- g$layout[!rm_grobs, ]
      ## move axis closer to panel
      
      g <- gtable::gtable_add_cols(g, unit(0.3, "null"), pos = 5)
      
      # g$layout[g$layout$name == "panel-1-1", c("t", "b")] <- c(, 8)
      grid.newpage()
      grid_plot <- grid.arrange(g)
      print(grid_plot)
      return(grid_plot)
    }

    
    
    save_plot <- function(file_name, width = 6, height = 4) {
      ggsave(str_glue("figures/{file_name}.pdf"), width = width, height = height, dpi = 75, device = cairo_pdf)
    }
    
    
    

# Functions for rescaling beta, contacts, etc. -------------------------------------------------

    
    
    
    
    # Function for combining the beta matrix and group proportions into one tibble 
    # Used in calculate_r0s and secondary_i_group
    # combine_beta_group_props <- function(contact_weights, group_props) {
    #   
    #   n_groups <- length(group_props)
    #   if (nrow(contact_weights) != n_groups^2) stop ("Mismatch between number of groups in beta matrix and group_props")
    #   
    #   group_props_df <- tibble(
    #     to = unique(contact_weights$to),
    #     group_pop_to = group_props
    #   )
    #   
    #   # Calculate the conditional probabilities of group A infecting group B
    #   combined <- contact_weights %>% 
    #     left_join(group_props_df, by = "to") %>% 
    #     mutate(beta_by_n = contact_weight * group_pop_to) %>% 
    #     group_by(from)
    #   
    #   combined
    # }
    
    
    k_matrix_trial <- Matrix::forceSymmetric(
      matrix(
        c(20, 1, 1, 1, 1, 1,
          10,  15, 1, 1, 1, 1,
          8,  9, 10, 1, 1, 1,
          5,  6, 7, 8, 1, 1,
          1,  1, 2, 2, 10, 1, 
          2,  2, 1, 1, 2, 6),
        nrow = 6
      )
    ) %>% as.matrix()
    

    # Calculates the contacts mean vector given the sar_out and contact_weights 
    # (ensuring consistency with contact_weights, sar_out, etc.)
    # scale_on and contact_means_scale chooses which group to set a scale parameter on, and what that group's contact mean is
    # calc_contact_means <- function(contact_weights, sar_out, group_props, contact_means_scale, scale_on = 1) {
    #   
    #   # contact_weights = crossing(to = 1:3, from = 1:3) %>% mutate(contact_weight = c(10, 5, 2, 5, 7, 1, 2, 1, 2))
    #   # group_props <- c(0.3, 0.4, 0.3)
    #   # contact_means_group1 <- 30
    #   # sar_out <- c(0.2, 0.1, 0.15)
    #   # scale_on <- 1
    #   
    #   
    #   beta_calcs <- combine_beta_group_props(contact_weights = contact_weights, group_props = group_props) %>% 
    #     summarise(total_beta_by_n = sum(beta_by_n), .groups = "drop") %>% 
    #     mutate(r0_ratio = total_beta_by_n / total_beta_by_n[from == scale_on])
    #   
    #   r0_group1 <- contact_means_scale * sar_out[[scale_on]]
    #   
    #   r0_all_groups <- beta_calcs$r0_ratio * r0_group1
    #   
    #   contact_means_all_groups <- r0_all_groups / sar_out
    #   
    #   return(contact_means_all_groups)
    #   
    # }
    
    
    # CONTACT_MEANS is now just rowSums(k_matrix)
    
  
    # Function that rescales number of each groups' contacts by a certain factor
    rescale_k <- function(k_matrix, rescale_factors) {
      
      # DEBUGGING:
      # k_matrix <- k_matrix_trial; rescale_factors <- c(1, 2, 3, 1, 1, 1)
      
      # OLD VERSION (converting to tibble)
      # k_tibble <- k_matrix %>% 
      #   as.vector() %>% 
      #   {mutate(crossing(to = 1:6, from = 1:6), k_val = .)} 
      # 
      # rescale_factors_df <- tibble(i_group = unique(k_tibble$to), scale_factor = rescale_factors)
      # 
      # k_tibble_rescaled <- k_tibble %>% 
      #   left_join(rescale_factors_df, by = c("to" = "i_group")) %>% 
      #   left_join(rescale_factors_df, by = c("from" = "i_group")) %>% 
      #   mutate(k_val_scaled = k_val * scale_factor.x * scale_factor.y) %>% 
      #   select(to, from, k_val = k_val_scaled)
      # 
      # k_matrix_rescaled <- matrix(k_tibble_rescaled$k_val, nrow = 6)
      
      if (nrow(k_matrix) != ncol(k_matrix)) stop("k_matrix must be square")
      
      # NEW VERSION: [matrix calc]
      n_groups <- length(rescale_factors)
      scale_x <- matrix(rep(rescale_factors, each = ncol(k_matrix)), nrow = n_groups, byrow = FALSE)
      scale_y <- matrix(rep(rescale_factors, each = nrow(k_matrix)), nrow = n_groups, byrow = TRUE)
      
      k_matrix * scale_x * scale_y
      
    }
    
 
    
    
    # Rescale the beta matrix given some rescale factors (e.g. double the contact rate of group 1)
    # rescale_beta <- function(contact_weights, rescale_factors) {
    #   # rescale_factors <- c(1.2, 1, 1)
    #   rescale_factors_df <- tibble(i_group = unique(contact_weights$to), scale_factor = rescale_factors)
    #   
    #   contact_weights %>% 
    #     left_join(rescale_factors_df, by = c("to" = "i_group")) %>% 
    #     left_join(rescale_factors_df, by = c("from" = "i_group")) %>% 
    #     mutate(contact_weight_scaled = contact_weight * scale_factor.x * scale_factor.y) %>% 
    #     select(to, from, contact_weight = contact_weight_scaled)
    # }
    

    
    
    
    # RESCALES THE CONTACT MEAN vector given some rescaling factors
    # rescale_contact_means <- function(contact_weights, rescale_factors, sar_out, group_props, contact_means_scale, scale_on = 1) {
    #   # r0_group_1 <- 2; group_props = c(0.5, 0.5)
    #   # group_props <- c(0.3, 0.4, 0.3)
    #   # contact_weights = crossing(to = 1:2, from = 1:2) %>% mutate(contact_weight = as.vector(matrix(c(10, 5, 5, 8), nrow = 2)))
    #   # contact_weights <- beta_default_3
    #   # scale_on <- 1; contact_means_scale <- 20
    #   # sar_out <- c(0.20, 0.10, 0.15)
    #   # sar_out <- c(0.2, 0.2, 0.2)
    #   # rescale_factors <- c(1.5, 1, 1)
    #   
    #   
    #   beta_by_n <- combine_beta_group_props(contact_weights = contact_weights,
    #                                         group_props = group_props) %>% 
    #     group_by(from) %>% 
    #     mutate(p_conditional = beta_by_n / sum(beta_by_n)) %>% 
    #     arrange(from) %>% 
    #     summarise(total_beta_by_n = sum(beta_by_n), .groups = "drop") %>% 
    #     .$total_beta_by_n
    #   
    #   contact_means_old <- calc_contact_means(
    #     contact_weights = contact_weights,
    #     sar_out = sar_out,
    #     group_props = group_props,
    #     contact_means_scale = contact_means_scale,
    #     scale_on = scale_on
    #   )
    #   
    #   
    #   beta_contacts_scale <- beta_by_n / (contact_means_old * sar_out) # scale factors to go from beta to contacts mean
    #   if (diff(range(beta_contacts_scale)) > 1e-5) stop("beta_contacts_scale factors should be the same across all groups")
    #   beta_contacts_scale <- beta_contacts_scale[[1]]
    #   
    #   rescaled_beta <- rescale_beta(
    #     contact_weights = contact_weights,
    #     rescale_factors = rescale_factors
    #   )
    #   
    #   beta_by_n_new <- combine_beta_group_props(contact_weights = rescaled_beta,
    #                                             group_props = group_props) %>% 
    #     group_by(from) %>% 
    #     mutate(p_conditional = beta_by_n / sum(beta_by_n)) %>% 
    #     arrange(from) %>% 
    #     summarise(total_beta_by_n = sum(beta_by_n), .groups = "drop") %>%
    #     .$total_beta_by_n
    #   
    #   contact_means_new <- beta_by_n_new  / (beta_contacts_scale * sar_out)
    #   
    #   return(contact_means_new)
    #   
    # }
    
    
    
    
    # THIS helps find the rescale factor that equalises the contacts across rich and poor group
    #(for the 1by1 parameter changes)
    # find_the_rescale_factor_OLD <- function(..., equalise = 1, equalise_to = 3) {
    #   
    #   rescale_testing <- tibble(
    #     rescale_1 = seq(0, 1, 0.01)
    #   ) %>% 
    #     rowwise() %>% 
    #     mutate(
    #       contact_means = list(rescale_contact_means(
    #         rescale_factors = c(rescale_1, 1, 1),
    #         ...
    #       ))
    #     )
    #   
    #   rescale_testing_w_condition <- rescale_testing %>%
    #     rowwise() %>% 
    #     mutate(i_group = list(c(1:length(contact_means)))) %>% 
    #     unnest(c(i_group, contact_means)) %>% 
    #     group_by(rescale_1) %>% 
    #     summarise(one_below_three = contact_means[i_group == equalise] <= contact_means[i_group == equalise_to]) 
    #   
    #   max(rescale_testing_w_condition$rescale_1[rescale_testing_w_condition$one_below_three])
    #   
    # }
    
    
    # find_the_rescale_factor(
    #   contact_weights = beta_default_3,
    #   # rescale_factors = c(1.5, 1, 1),
    #   sar_out = c(0.2, 0.1, 0.15),
    #   group_props = c(0.1, 0.4, 0.4),
    #   contact_means_scale = 8,
    #   scale_on = 3,
    #   equalise = 1,
    #   equalise_to = 3
    # )
    
    
    # FUnction to find the rescale factors that equalises a given groups' contact mean to the "equalise_to" group
    find_the_rescale_factor <- function(rescale_factors, k_matrix, pops, equalise, equalise_to) {
      
      # NOTE: only input the rescale_factors of the groups in "equalise" (not all the rescale_factors)
      if (length(rescale_factors) != length(equalise)) stop("length of rescale_factors and equalise should be the same")
      
      # TESTING
      # rescale_factors <- c(0.5); k_matrix <- data_save$k_matrix_5; equalise <- 1; equalise_to = 5; pops <- calc_pops(data_save$group_props, 10000)
      
      n_groups <- ncol(k_matrix)
      
      # CALCULATE rescale factors vector
      rescale_factors_all <- rep(1, n_groups) # fill with ones
      rescale_factors_all[equalise] <- rescale_factors # replace with tentative rescale factors
      
      # Calculate new k_matrix
      k_rescaled <- rescale_k(k_matrix, rescale_factors = rescale_factors_all)
      
      # Calculate implied contact means
      contact_matrix_ind <- k_rescaled / pops
      contact_means_rescaled <- rowSums(contact_matrix_ind) # remember, overall average of contacts_mean will change here!
      
      # Calculate difference (deviation from equalise_to group)
      diff <- sum(abs(contact_means_rescaled[equalise] - contact_means_rescaled[equalise_to]))   # this should be 0
      
      # Calculate penalty term (don't want to deviate more away from 1 in one group, spread evenly)
      penalty <- diff(range(rescale_factors))
      
      return(diff + penalty)
      
    }
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    # library("nloptr")
    # rescale_factor_opt <- nloptr(
    #   x0 = c(0.5, 0.7),
    #   eval_f = find_the_rescale_factor,
    #   ub = c(1, 1) + 1e-5,
    #   opts = list(
    #     "algorithm" = "NLOPT_LN_SBPLX",
    #     "xtol_abs"=1.0e-10,
    #     "maxeval" = 100,
    #     "print_level" = 3
    #   ),
    #   contact_weights = beta_default_3,
    #   sar_out = c(0.2, 0.2, 0.15),
    #   group_props = c(0.3, 0.4, 0.3),
    #   contact_means_scale = 8,
    #   scale_on = 3,
    #   equalise = c(1, 2),
    #   equalise_to = 3
    # )
    # 
    # rescale_factor_opt$solution
    
    # RESULTING RESCALED CONTACTS MEAN AND BETA MATRIX
    # rescale_contact_means(
    #   contact_weights = beta_default_3,
    #   rescale_factors = c(rescale_factor_opt$solution, 1),
    #   sar_out = c(0.2, 0.2, 0.15),
    #   group_props = c(0.3, 0.4, 0.3),
    #   contact_means_scale = 8,
    #   scale_on = 3
    # )
    
    # rescale_beta(
    #   contact_weights = beta_default_3,
    #   rescale_factors = c(rescale_factor_opt$solution, 1)
    # )
    
  
    
    
    # FUNCTION used for the "inequality is bad" argument
    
    
    # Function that takes parameters and tightens them towards the mean
    # (taking into account the non-uniform group sizes)
    # so that overall average stays the same
    rescale_tighten <- function(param, group_props, tighten_factor, tighten_to = NULL) {
      
      # when tighten_to is NULL, it tightens to the weighted.mean of population;
      # otherwise tightens to that group
      
      if (tighten_factor > 1) stop("scale factor can't be >1")
      if (tighten_factor < 0) stop("scale factor can't be <0")
      if (length(param) != length(group_props)) {
        print(param); print(group_props)
        stop("length of param needs to be same as length of group_props")
      }
      
      # Example:
      # param <- c(9, 8, 7, 6); group_props <- data_save$group_props; tighten_factor <- 0; tighten_to <- NULL
      if (is.null(tighten_to)) {
        weighted_mean <- weighted.mean(param, group_props)
        tighten_to_value <- weighted_mean
      } else if (tighten_to %in% 1:length(param)) {
        tighten_to_value <- param[tighten_to]
      } else {
        stop("Error: tighten_to is specified incorrectly")
      }
      
      deviation <- param - tighten_to_value
      scaled_deviation <- deviation * tighten_factor
      scaled_vals <- scaled_deviation + tighten_to_value
      
      # Check the new weighted mean is the same
      if (is.null(tighten_to)) {
        new_weighted_mean <- weighted.mean(scaled_vals, group_props)
        if (abs(new_weighted_mean - weighted_mean) > 1e-5) stop("weighted mean has changed")
      }
      
      return(scaled_vals)
      
    }
    
   
      
    find_the_rescale_factor_vec <- function(rescale_factors, # to try
                                            contact_means_aim, # to aim for
                                            which_groups = NULL,
                                            k_matrix,
                                            pops) {
      
      # TESTING:
      # k_matrix <- data_save$k_matrix_5; rescale_factors <- c(0.5, 0.7, 0.7, 0.7, 0.8); contact_means_aim <- c(30, 30, 20, 20, 20); pops <- calc_pops(group_props = data_save$group_props, n_pop = 10000)
      # which_groups <- 1:2; rescale_factors <- c(0.5, 0.7)
      
      # List of groups to manipulate
      n_groups <- nrow(k_matrix)
      all_groups <- 1:n_groups
      if (is.null(which_groups)) which_groups <- 1:n_groups
      if (length(which_groups) != length(rescale_factors)) {
        # print(which_groups)
        # print(rescale_factors)
        stop("length(which_groups) != length(rescale_factors)")
        
      }
      
      
      # Original contact means
      contact_means_original <- rowSums(k_matrix / pops)[which_groups]
      
      rescale_factors_all <- rep(1, n_groups)
      rescale_factors_all[which_groups] <- rescale_factors
      
      # New contact means after rescaling
      contact_means_rescaled <- rescale_k(k_matrix, rescale_factors = rescale_factors_all) %>% 
        {rowSums(. / pops)}
      
      # Difference between rescaled and the aim
      diff <- sum(abs(contact_means_rescaled[which_groups] - contact_means_aim))
      
      # Penalty term - make sure rescale factors are close to the ratio of the contact_means [avoids changing groups that shouldn't be changed]
      penalty <- sum((rescale_factors - (contact_means_aim / contact_means_original))^2)
      
      return(diff + penalty)
      
    }
    
    
    
    # FUNCTION THAT TAKES LIST OF PARAMETERS AS INPUTS, and TIGHTENS ALL OF THEM
    tighten_params <- function(params, tighten_factor, tighten_to = NULL, all_at_once) {
      
      # TESTING: 
      # tighten_factor <- 0; params <- params_baseline; tighten_to <- 4; all_at_once = FALSE
      # tighten_factor <- 0; params <- params_baseline; tighten_to <- NULL; all_at_once <- TRUE
      
      if (tighten_factor < 0 || tighten_factor > 1) stop("tighten factor must be 0 < t < 1")
      
      # Calculate contacts_mean from k_matrix
      pops <- calc_pops(group_props = params$group_props, n_pop = params$n_pop)
      params$contact_means <- rowSums(params$k_matrix / pops)
      
      # Extract group_props (shouldn't be tightened)
      group_props <- params$group_props
      params$group_props <- NULL; params$n_pop <- NULL
      
      # if (tighten_factor == 1) {
      #   return(params)
      # }   # tighten factor 1 means no tightening, just return the params
    
      n_groups <- length(group_props)
    
      # For the atomic ones
      atomic_tightened <- params %>% keep(~ is_atomic(.) & !is.matrix(.)) %>% 
        map(~ rescale_tighten(.x, group_props = group_props, tighten_factor = tighten_factor, tighten_to = tighten_to)) %>% print
      
      # For the probs_df ones
      probs_df_tightened <- params %>% keep(str_detect(names(.), "_df$")) %>% 
        bind_rows(.id = "param_type") %>% 
        group_by(param_type, symptom_severity) %>% 
        left_join(tibble(i_group = 1:n_groups, group_props = data_save$group_props), by = "i_group") %>% 
        mutate(prob_rescaled = rescale_tighten(prob, group_props = group_props, tighten_factor = tighten_factor, tighten_to = tighten_to)) %>% 
        select(param_type, i_group, symptom_severity, prob = prob_rescaled) %>% 
        split(.$param_type) %>% 
        map(ungroup) %>% 
        map(~ select(., -param_type))
      
      
      # For the isolation parameters
      probs_basic <-  crossing(i_group = 1L:n_groups, symptom_severity = 1L:2L)
      
          # baseline_isolation <- 1 - (params$days_work_ct / params$days_of_work)
          # tightened_isolation <- 1 - (atomic_tightened$days_work_ct / atomic_tightened$days_of_work)
      
          # In case of each channel treated separately, need to divide by OLD days_of_work (params$days_of_work)
          if (!all_at_once) {
            
            # (i) isolation CT
            p_isolate_ct <- 1 - (atomic_tightened$days_work_ct / params$days_of_work) # when changing isolation behaviour, keep days of work the same
            p_isolate_ct <- if_else(p_isolate_ct < 0, 0, p_isolate_ct)
            probs_isolate_ct <- probs_basic %>% mutate(prob = rep(p_isolate_ct, each = 2))   # independent of symptom severity
            
            # (ii) Isolate symptoms
            p_isolate_symp <- 1 - (atomic_tightened$days_work_symptoms / params$days_of_work)
            p_isolate_symp <- if_else(p_isolate_symp < 0, 0, p_isolate_symp)
            probs_isolate_symptoms <- probs_basic %>% 
              mutate(prob = rep(p_isolate_symp, each = 2),
                     prob = if_else(symptom_severity == 1, 0, prob)) 
            
            # (iii) P hh quarantine
            p_hh_quarantine <- 1 - (atomic_tightened$days_work_hh_quarantine / params$days_of_work)
            p_hh_quarantine <- if_else(p_hh_quarantine < 0, 0, p_hh_quarantine)
          
          } else if (all_at_once) {
            
            # In case of each channel treated all at once (for inequality_is_bad graph), need to divide by NEW days_of_work (atomic_tightened$days_of_work)
            
            # (i) isolation CT
            p_isolate_ct <- 1 - (atomic_tightened$days_work_ct / atomic_tightened$days_of_work) # when changing isolation behaviour, keep days of work the same
            p_isolate_ct <- if_else(p_isolate_ct < 0, 0, p_isolate_ct)
            probs_isolate_ct <- probs_basic %>% mutate(prob = rep(p_isolate_ct, each = 2))   # independent of symptom severity
            
            # (ii) Isolate symptoms
            p_isolate_symp <- 1 - (atomic_tightened$days_work_symptoms / atomic_tightened$days_of_work)
            p_isolate_symp <- if_else(p_isolate_symp < 0, 0, p_isolate_symp)
            probs_isolate_symptoms <- probs_basic %>% 
              mutate(prob = rep(p_isolate_symp, each = 2),
                     prob = if_else(symptom_severity == 1, 0, prob)) 
            
            # (iii) P hh quarantine
            p_hh_quarantine <- 1 - (atomic_tightened$days_work_hh_quarantine / atomic_tightened$days_of_work)
            p_hh_quarantine <- if_else(p_hh_quarantine < 0, 0, p_hh_quarantine)
          }
      
          
          
      
      
      # For the testing delays
      tighten_delays <- function(.data, var, tighten_factor, tighten_to, group_props) {
        means <- .data %>% 
          group_by(i_group) %>% 
          summarise(mean = mean({{var}}, na.rm = TRUE), .groups = "drop") %>% 
          .$mean
        
        means_tightened <- rescale_tighten(means, group_props = group_props, tighten_factor = tighten_factor, tighten_to = tighten_to)
        
        means_diff <- means_tightened - means
        
        data_with_means <- .data %>% 
          left_join(tibble(i_group = 1:length(group_props),
                           means_diff),
                    by = "i_group")
        
        data_upd <- data_with_means %>% 
          mutate("{{var}}" := if_else({{var}} + means_diff < 0, {{var}}, {{var}} + means_diff))
       # this ensures there's no new negative values (makes very little difference)
        
        # CHECK with a print
        # data_with_means %>% group_by(i_group) %>% summarise(across(everything(), mean)) %>% print
        
        data_upd %>% select(-means_diff)
        
      }
      
      
      test_delay_data_tightened <- params$test_delay_data %>% 
        tighten_delays(var = test_choice_delay, tighten_factor = tighten_factor, tighten_to = tighten_to, group_props = group_props) %>% 
        tighten_delays(var = test_results_delay, tighten_factor = tighten_factor, tighten_to = tighten_to, group_props = group_props)
      
      ct_delay_data_tightened <- params$ct_delay_data %>% 
        tighten_delays(var = ct_delay, tighten_factor = tighten_factor, tighten_to = tighten_to, group_props = group_props)
      
      # FOR HH SIZE DATA
      print("calculating new hh size data")
      
      hh_size_by_hh <- params$hh_size_data %>%
        select(hh_id, i_group, hh_size) %>%
        dups_drop(warn = FALSE)
      
      hh_size_means <- params$hh_size_data %>% 
        group_by(i_group) %>% 
        summarise(hh_size = mean(hh_size), .groups = "drop") %>% 
        .$hh_size
      
      hh_size_means_tightened <- rescale_tighten(hh_size_means, group_props = group_props, tighten_factor = tighten_factor, tighten_to = tighten_to)
      hh_size_means_diff <- hh_size_means_tightened - hh_size_means
      
      # i_group <- 2
      # direction <- -1
      # desired_mean <- hh_size_means_tightened[[2]]
      
      # iterate_hh_size(hh_size_by_hh, i_group = 2, direction = -1, desired_mean = hh_size_means_tightened[[2]])
      
      
      iterate_hh_size <- function(hh_size_by_hh, i_group, direction, desired_mean) {
        
        hh_i_group <- hh_size_by_hh %>% filter(i_group == !!i_group)
        
        if (direction == -1) {
          hh_i_group_filt <- hh_i_group %>% filter(hh_size >= 2)
        } else {
          hh_i_group_filt <- hh_i_group
        }
        
        # MEAN HH SIZE at the individual level is equal to 
        # sum (hh_sizes at hh level ^2 ) / sum(hh sizes at hh level)
        # e.g. because having a hh of 3 means that there are 3 people with that hh size
        
        # total pop in this group [original denominator on the mean]
        original_hh_sq <- sum(hh_i_group$hh_size ^ 2)
        original_N <- sum(hh_i_group$hh_size)
        original_mean <- original_hh_sq / original_N
          
        # Order randomly and keep only the first N before n_to_iterate
        
        
        hh_i_group_iterated <- hh_i_group_filt %>% 
          mutate(rand = runif(nrow(.))) %>% 
          arrange(rand) %>% 
          mutate(
            old_size = hh_size,
            new_size = hh_size + direction,
            added_to_numerator = new_size^2 - old_size^2,
            added_to_denominator = new_size - old_size,
            cum_num = cumsum(added_to_numerator),
            cum_den = cumsum(added_to_denominator),
            
            new_numerator = original_hh_sq + cum_num,
            new_denominator = original_N + cum_den,
            new_mean = new_numerator/ new_denominator
          )
        
        # Keep only the IDs required to get to the right mean
        if (direction == -1) {
          hh_i_group_required <- hh_i_group_iterated %>% filter(new_mean >= desired_mean)
        } else if (direction == 1) {
          hh_i_group_required <- hh_i_group_iterated %>% filter(new_mean <= desired_mean)
        } else {
          # warning("hh_i_group_required is NULL - no changes made")
          # print("no changes to hh_size made")
          return(hh_i_group)
        }
        
        hh_i_group_required <- hh_i_group_required %>% select(hh_id, i_group, new_size)
        
        
        hh_i_group_upd <- hh_i_group %>% 
          left_join(hh_i_group_required, by = c("hh_id", "i_group")) %>% 
          mutate(hh_size = if_else(!is.na(new_size), new_size, hh_size)) %>% 
          select(-new_size)
      
        return(hh_i_group_upd)
        
      }
      
      hh_iterated <- tibble(
        i_group = 1:n_groups,
        desired_mean = hh_size_means_tightened,
        direction = case_when(hh_size_means_diff < 0 ~ -1, hh_size_means_diff == 0 ~ 0, hh_size_means_diff > 0 ~ 1)
      ) %>% 
        rowwise() %>% 
        mutate(
          hh_new = list(iterate_hh_size(hh_size_by_hh = hh_size_by_hh,
                                        i_group = i_group,
                                        direction = direction,
                                        desired_mean = desired_mean))
        )
      
      
      # NEED TO MAKE SURE THERE ARE NO 0 HH SIZES
      
      hh_size_new <- hh_iterated %>% 
        select(hh_new) %>% 
        unnest(hh_new) %>% 
        arrange(hh_id)
      
      hh_size_new_ind <- tibble(
        hh_id =   rep(hh_size_new$hh_id, times = hh_size_new$hh_size),
        i_group = rep(hh_size_new$i_group, times = hh_size_new$hh_size),
        hh_size = rep(hh_size_new$hh_size, times = hh_size_new$hh_size),
      ) %>% 
        mutate(
          hh_ind_id = paste0(hh_id, "_", 1:nrow(.))
        )
      
      # CHECK WITH A GRAPH
      # hh_size_new_ind %>%
      #   group_by(i_group) %>%
      #   summarise(mean_se(hh_size)) %>%
      # 
      #   ggplot(aes(x = factor(i_group), y = y)) +
      #   geom_point() +
      #   geom_errorbar(aes(ymin = ymin, ymax = ymax), width = 0.2)
    
      
  
      # FOR CONTACT MATRIX
      # (contact_means has already been tightened above)
      # Calculate rescale_factors needed to align k_matrix with the contacts mean
      print("finding rescale factors for K matrix")
      rescale_factor_tighten <- nloptr::nloptr(x0 = rep(1, n_groups), 
                                               eval_f = find_the_rescale_factor_vec, 
                                               lb = rep(0, n_groups) + 1e-5,
                                               opts = list(
                                                 "algorithm" = "NLOPT_LN_SBPLX",
                                                 "xtol_abs"=1.0e-10,
                                                 "maxeval" = 1000,
                                                 "print_level" = 0
                                               ),
                                               contact_means_aim = atomic_tightened$contact_means,
                                               k_matrix = params$k_matrix, 
                                               pops = pops,
                                               which_groups = 1:n_groups)  
      
      rescale_tighten_solution <- rescale_factor_tighten$solution
      k_tightened <- rescale_k(params$k_matrix, rescale_tighten_solution)
      
      # Export all params tightened up
      params_tightened <- c(
        atomic_tightened,
        probs_df_tightened
      )
      
      params_tightened$test_delay_data <- test_delay_data_tightened
      params_tightened$ct_delay_data <- ct_delay_data_tightened
      params_tightened$k_matrix <- k_tightened
      params_tightened$hh_size_data <- hh_size_new_ind
      
      params_tightened$probs_isolate_symptoms_df <- probs_isolate_symptoms
      params_tightened$probs_isolate_ct_df <- probs_isolate_ct
      params_tightened$p_hh_quarantine <- p_hh_quarantine
      
      
      # if (any(map_lgl(params_tightened, is.null))) {
      #   browser()
      # }
      
      # Check same params in input and output
      # if (!setequal(names(params), names(params_tightened))) stop("different params in params and params_tightened")
      
      return(params_tightened)
      
    }
    
    
      
   
    
    # rescale_tighten(param = c(20, 30, 40, 100, 30, 40), 
    #                 group_props = data_save$group_props,
    #                 tighten_factor = 0.5)    
    
    # Function that repeats a tibble for all of the groups 1:5
    rep_tibble <- function(.data, n_groups) {
      map(1:n_groups, ~ .data %>% mutate(i_group = .x)) %>% bind_rows()
    }
    
    
    
    # FUNCTION THAT TAKES LIST OF PARAMETERS AS INPUTS, and EQUALISES THEM ALL TO LEVEL OF THE EQUALISE_TO group
    equalise_params_to <- function(params, equalise_to) {
      
      # TESTING:
      # params <- params_baseline; equalise_to <- 4
      
      n_groups <- length(params$group_props)
      
      # Calculate contacts_mean from k_matrix
      pops <- calc_pops(group_props = params$group_props, n_pop = params$n_pop)
      contact_means <- rowSums(params$k_matrix / pops)
      
      # Extract group_props (shouldn't be tightened)
      group_props <- params$group_props
      params$group_props <- NULL; params$n_pop <- NULL
  
      
      # DFs (probs and delays)
      param_dfs <- params %>% keep(is_tibble) %>% 
        purrr::list_modify("hh_size_data" = NULL) %>%  # don't do it for hh size data (needs to be treated separately)
        map(~ .x %>% filter(i_group == equalise_to) %>% rep_tibble(n_groups = n_groups))
      
      # Vectors
      param_vecs <- params %>% keep(~ is_atomic(.) & !is.matrix(.)) %>% 
        map(~ rep(.x[[equalise_to]], times = n_groups))
      

      
      
      # Calculate new isolation DFs (for CT, symptoms) [when just adjusting isolation, keep days of work the same]
      probs_basic <-  crossing(i_group = 1L:n_groups, symptom_severity = 1L:2L)
      p_isolate_ct <- 1 - (params$days_work_ct[[equalise_to]] / params$days_of_work)
      p_isolate_ct <- if_else(p_isolate_ct < 0, 0, p_isolate_ct)
      probs_isolate_ct <- probs_basic %>% mutate(prob = rep(p_isolate_ct, each = 2))
      p_isolate_symp <- 1 - (params$days_work_symptoms[[equalise_to]] / params$days_of_work)
      p_isolate_symp <- if_else(p_isolate_symp < 0, 0, p_isolate_symp)
      
      probs_isolate_symptoms <- probs_basic %>% 
        mutate(prob = rep(p_isolate_symp, each = 2),
               prob = if_else(symptom_severity == 1, 0, prob))
      
      # P hh quarantine
      p_hh_quarantine <- 1 - (params$days_work_hh_quarantine[[equalise_to]] / params$days_of_work)
      p_hh_quarantine <- if_else(p_hh_quarantine < 0, 0, p_hh_quarantine)
      
      
      # K matrix
      x0_start <- contact_means[[equalise_to]] / contact_means[-equalise_to]
      
      # Optimise to get rescale factors
      # library("nloptr")
      # print(find_the_rescale_factor)
      # print(x0_start)
      # print(pops)
      # print(equalise_to)
      
      rescale_factor_opt <- nloptr::nloptr(x0 = x0_start, 
                                           eval_f = find_the_rescale_factor, 
                                           # ub = c(1, 1) + 1e-5,
                                           opts = list(
                                             "algorithm" = "NLOPT_LN_SBPLX",
                                             "xtol_abs"=1.0e-15,
                                             "maxeval" = 1000,
                                             "print_level" = 0
                                           ),
                                           k_matrix = params$k_matrix,
                                           pops = pops,
                                           equalise = (1:n_groups)[-equalise_to],
                                           equalise_to = equalise_to)
      
      k_matrix_rescaled <- rescale_k(k_matrix = params$k_matrix,
                                     rescale_factors = c(rescale_factor_opt$solution, 1))
      
      
      # HH size data
      hh_size_vec <- params$hh_size_data %>% 
        select(hh_id, i_group, hh_size) %>% 
        dups_drop(warn = FALSE) %>% 
        filter(i_group == equalise_to) %>% 
        pull(hh_size)
      
      hh_size_data_equalised <- tibble(
        i_group = rep(1:n_groups, each = 100000),
        hh_size = sample(hh_size_vec, size = 100000 * n_groups, replace = TRUE)
      ) %>% 
        mutate(hh_id = row_number())
      
      hh_size_data_equalised_long <- tibble(
        i_group = rep(hh_size_data_equalised$i_group, hh_size_data_equalised$hh_size),
        hh_size = rep(hh_size_data_equalised$hh_size, hh_size_data_equalised$hh_size),
        hh_id = rep(hh_size_data_equalised$hh_id, hh_size_data_equalised$hh_size)
      ) %>% 
        mutate(hh_ind_id = paste0(hh_id, "_", row_number())) %>% 
        relocate(hh_id, i_group, hh_size, hh_ind_id)
      
      # Export all params tightened up
      params_tightened <- c(
        param_vecs,
        param_dfs
      )
      
      params_tightened$k_matrix <- k_matrix_rescaled
      params_tightened$hh_size_data <- hh_size_data_equalised_long
      
      
      params_tightened$probs_isolate_symptoms_df <- probs_isolate_symptoms
      params_tightened$probs_isolate_ct_df <- probs_isolate_ct
      params_tightened$p_hh_quarantine <- p_hh_quarantine
      
      # Check same params in input and output
      # if (!setequal(names(params), names(params_tightened))) stop("different params in params and params_tightened")
      
      return(params_tightened)
      
    }
    
    
    
    # close_params_gap <- function(params, equalise_to, gap_factor) {
    #   
    #   # TESTING:
    #   # params <- params_baseline; equalise_to <- 4; gap_factor <- 0.5
    #   
    #   if (gap_factor < 0 || gap_factor > 1) stop("gap factor must be 0 < g < 1")
    #   
    #   
    #   # Calculate contacts_mean from k_matrix
    #   pops <- calc_pops(group_props = params$group_props, n_pop = params$n_pop)
    #   params$contact_means <- rowSums(params$k_matrix / pops)
    #   
    #   # Extract group_props (shouldn't be tightened)
    #   group_props <- params$group_props
    #   params$group_props <- NULL; params$n_pop <- NULL
    #   
    #   if (gap_factor == 1) {
    #     return(params)
    #   }   # gap factor 1 means no tightening, just return the params
    #   
    #   n_groups <- length(group_props)
    #   
    #   
    #   # For the atomic ones
    #   atomic_tightened <- params %>% keep(~ is_atomic(.) & !is.matrix(.)) %>% 
    #     map(~ rescale_tighten(.x, group_props = group_props, tighten_factor = tighten_factor)) %>% print
    #   
    #   
    #   
    #   
    # }
    
    

    
