    # source("individual_model_functions.R")    # imports all functions
    # load("data/processed/data_save.RData")                         # loads all cleaned data
    



# Parameters --------------------------------------------------------------

    set.seed(12345)
    
    n_pop_for_k_calc <- 10000
    
    
    

# Calculate contacts and tracing probability ------------------------------

    contacts_outside_home <- data_save$params_data$outside_non_work + (data_save$params_data$days_of_work * outside_work_factor)
    
    ind_mean_hh_size <- data_save$hh_data_bogota %>% 
      group_by(i_group) %>% 
      summarise(ind_mean_hh_size = mean(hh_size, na.rm = TRUE)) %>% 
      .$ind_mean_hh_size
    
    total_contacts <- contacts_outside_home + (ind_mean_hh_size - 1)
    p_contact_traced <- data_save$params_data$contacts_traced / total_contacts
    
    # Save to data_save
    data_save$contacts_outside_home <- contacts_outside_home
    data_save$p_contact_traced <- p_contact_traced
    
  


# K with population size adjustment ----------------------------------


# FUNCTIONS ---------------------------------------------------------------



    # Takes k_vec and mu_vec as inputs, and outputs the full k_matrix
    k_vec_to_matrix <- function(k_vec, mu, pops) {
      
      # k_vec <- c(c(10, 5, 2),   # k_1
      #            c(50, 15),
      #            c(30)) * 1000
      # mu <- 10:7 - 1/3
      # group_props <- data_save$group_props
      # # group_props[[5]] <- group_props[[5]] + data_save$group_props[[6]]
      # n_pop <- 23000
      
      if (length(k_vec) != 3+2+1) {
        print(length(k_vec))
        stop("wrong length of k_vec")
      }
      
      # Total contacts in each group
      total_contacts <- as.integer(round(mu * pops))
      
      k_1_almost <- k_vec[1:3]
      k_2_almost <- k_vec[4:5]
      k_3_almost <- k_vec[6]
      # k_4 is pinned down by the mu
      
      k_1 <- append(k_1_almost, total_contacts[[1]] - sum(k_1_almost))
      k_2 <- c(k_1[[2]], k_2_almost) %>% append(total_contacts[[2]] - sum(.))
      k_3 <- c(k_1[[3]], k_2[[3]], k_3_almost) %>% append(total_contacts[[3]] - sum(.))
      k_4 <- c(k_1[[4]], k_2[[4]], k_3[[4]]) %>% append(total_contacts[[4]] - sum(.))
      # k_5 <- c(k_1[[5]], k_2[[5]], k_3[[5]], k_4[[5]]) %>% append(total_contacts[[5]] - sum(.))
      
      # Make the k matrix for clarity
      k_matrix <- matrix(
        c(
          k_1, k_2, k_3, k_4
        ),
        nrow = 4, 
        byrow = TRUE
      )
      
      return(k_matrix)
      
    }
    
    
    k_matrix_to_vec <- function(k_matrix) {
      
      i_row <- c(rep(1, 3), rep(2, 2), rep(3, 1))
      i_col <- c(1:3, 2:3, 3)
      
      k_vec <- list()
      for (i in seq_along(i_row)) {
        k_vec[[i]] <- as.vector(k_matrix[ i_row[[i]], i_col[[i]] ])
      }
      
      return(flatten_dbl(k_vec))
      
    }
    
    
 
    
    
    
    # Takes the parameter vector as an input, and outputs all the implied parameters
    # e.g. the k_matrix, the q_matrix, etc.
    read_params <- function(k_vec, mu, group_props, n_pop) {
      
      
      # k_vec <-  c(c(10, 5, 2, 2),   # k_1
      #             c(50, 15, 5),
      #             c(30, 5),
      #             c(5)) * 1000
      # mu <- rep(10, 5)
      # group_props <- data_save$group_props[1:5]
      # n_pop <- 23000
      
      # Calculate population size
      pops <- calc_pops(group_props, n_pop)
      
      
      k_matrix <- k_vec_to_matrix(k_vec = k_vec, mu = mu, pops = pops)
      
      # Check there's no negative values
      # if (any(k_matrix < 0, na.rm = TRUE) | any(is.na(k_matrix))) {
      #   # print("Negative or missing values in k")
      #   # return(Inf)
      #   # WILL NEED TO RETURN INFINITY HERE
      # }
    
      
      # mu_matrix <- matrix(
      #   rep(mu, each = 5), nrow = 5,
      #   byrow = FALSE
      # )
      
      
      q_matrix <- k_matrix / rowSums(k_matrix)     # rows sum to 1
    
      # Convert q_matrix into form that can be joined onto tibble
      q_tibble <- q_matrix %>% as_tibble() %>% 
        set_names(paste0("primary", 1:4)) %>% 
        mutate(secondary_i_group = row_number()) %>% 
        pivot_longer(-secondary_i_group, names_to = "i_group", names_prefix = "primary", values_to = "q_estimate") %>% 
        mutate(i_group = as.integer(i_group))
      
      
      # Calculate beta_matrix
      pop_matrix <- pops %*% t(pops)
      beta_matrix <- k_matrix / pop_matrix #* (n_pop^2)
      
      
      out <- list(
        mu = mu,
        k_vec = k_vec,
        k_matrix = k_matrix,
        # mu_matrix = mu_matrix,
        q_matrix = q_matrix,
        q_tibble = q_tibble,
        beta_matrix = beta_matrix,
        pops = pops
      )
      
      return(out)
      
    }
    
    log_na <- function(x) {
      ifelse(is.finite(log(x)), log(x), NA)
    }
  
    
    # Calculate the log likelihood of observing data given parameters
    contact_maxlik_k_only <- function(k_vec, mu, group_props, n_pop, data) {
      
      
      # INPUTS FOR TESTING
      # k_vec <-  c(c(10, 5, 2, 2),   # k_1
      #             c(50, 15, 5),
      #             c(30, 5),
      #             c(5)) * 1000
      # mu <- rep(10, 5)
      # group_props <- data_save$group_props[1:5]
      # n_pop <- 23000
      # data <- shitty_data_example
      
      params <- read_params(k_vec = k_vec, mu = mu, group_props = group_props, n_pop = n_pop)
      print(round(params$k_matrix, 2))
      
      if (any(is.na(params$k_matrix) | is.nan(params$k_matrix))) {
        return(1e11 + runif(0, 1000)) 
      }
      
      
      # Check there's no negative / missing values
      if (any(params$k_matrix < - 0.01, na.rm = TRUE)) {
        print("Negative or missing values in k")
        
        k_neg <- params$k_matrix %>% as.vector() %>% 
          .[. < -0.01] %>% 
          sum()
        
        return(1e7 + abs(k_neg)) # penalty to make sure it doesn't return negative values
      }
      
      # Check symmetric
      if (!isSymmetric(params$k_matrix)) {
        print(params$k_matrix)
        stop("matrix is not symmetric")
      }
      
      # Add on p_estimate for number of contacts
      # data_p <- data %>% dups_drop(id, warn = FALSE) %>% 
      #   left_join(tibble(i_group = 1:6, mu = params$mu), by = "i_group") %>% 
      #   mutate(
      #     p_estimate = dnbinom(x = contacts_total, size = params$size, mu = mu)
      #   )
      
      

      
      # And q_estimate on group of contact from q_matrix
      data_q <- data  %>% 
        left_join(params$q_tibble, by = c("i_group", "secondary_i_group")) %>% 
        mutate(q_estimate = replace_na(q_estimate, 1))
      
      # Negative log likelihood (to be minimised)
      neg_ll <- -sum(log(data_q$q_estimate))
      
      if (is.nan(neg_ll)) {
        print("neg_ll is nan")
        return(1e5 + -sum(log_na(data_q$q_estimate), na.rm = TRUE))
      }
      
      return(neg_ll)
      
    }
    
    
    # Constraint function to prevent negative values in K matrix
    k_constraint <- function(k_vec, mu, group_props, n_pop, data) {
      k_vals <- read_params(k_vec, mu, group_props, n_pop)$k_matrix %>% as.vector()
      return(-k_vals)
    }
    
    

# ............. -----------------------------------------------------------


# EXAMPLE / TEST [unused] -----------------------------------------------------------------


    
# Example data from matrix -----------------------------------------
    
#     q_data <- matrix(
#       c(0, 8, 5, 2, 0,
#         8, 73, 50, 4, 1,
#         6, 70, 165, 35, 8,
#         1, 3, 26, 22, 11,
#         0, 4, 3, 2, 5),
#       nrow = 5,
#       byrow = TRUE
#     ) %>% 
#       {. / rowSums(.)}
#     
#     q_tibble <- q_data %>% as_tibble() %>% 
#       set_names(paste0("primary", 1:5)) %>% 
#       mutate(secondary_i_group = row_number()) %>% 
#       pivot_longer(-secondary_i_group, names_to = "i_group", names_prefix = "primary", values_to = "q_estimate") %>% 
#       mutate(i_group = as.integer(i_group))
#     
#     # CREATE example data using those probabilities
#     data_example <- tibble(
#       i_group = sample(1:5, size = 300, replace = TRUE, prob = data_save$group_props[1:5])
#     ) %>% 
#       group_by(i_group) %>% 
#       group_split() %>% 
#       map2(.y = q_tibble %>% group_by(i_group) %>% group_split() %>% map("q_estimate"),
#            .f = ~ {
#              .x %>% mutate(secondary_i_group = sample(1:5, size = nrow(.), replace = TRUE, prob = .y))
#            }) %>% 
#       bind_rows()
#     
#     
#     
#     
# 
# # Test new functions -------------------------------------------------------
# 
#     mu_example <- 10:6
#     n_pop_total <- 10000
#     
#     # Calculate the best guess starting K matrix 
#     # assume mu is 8 for everyone, and they randomly mix
#     pops <- calc_pops(group_props = data_save$group_props, n_pop = n_pop_total)
#     pop_matrix <- pops %>% {. %*% t(.)} %>% print
#     beta_val <- pop_matrix %>% {10 * pops / rowSums(.)} 
#     k_start <- pop_matrix * beta_val[[1]]
#     
#     read_params(k_vec = k_matrix_to_vec(k_start),
#                 mu = mu_example,
#                 group_props = data_save$group_props,
#                 n_pop = n_pop_total)
#     
#     # TRY MAX LIKELIHOOD AGAIN
#     k_optimised <- nloptr::nloptr(
#       x0 =  k_matrix_to_vec(k_start),
#       eval_f = contact_maxlik_k_only, 
#       lb = rep(0, 10),  # entries to K can't go below 0 
#       eval_g_ineq = k_constraint,  # constraint ensures that no entailed elements of K are negative
#       opts = list(
#         "algorithm" = "NLOPT_LN_COBYLA",
#         "xtol_abs"=1.0e-6,
#         "maxeval" = 500,
#         "print_level" = 1
#       ),
#       n_pop = n_pop_total,
#       group_props = data_save$group_props,
#       mu = mu_example,
#       data = data_example
#     )
#     
# 
#     solution <- read_params(
#       k_vec = k_optimised$solution,
#       mu = mu_example,
#       group_props = data_save$group_props,
#       n_pop = 10000
#     )
#     
#     solution$q_matrix %>% round(2)
#     solution$k_matrix %>% round(0)
# 
# # Export example K matrix -------------------------------------------------
# 
# 
#     
#     data_save$k_matrix_5 <- round(solution$k_matrix)
    
    
    
    

    

# ...... ------------------------------------------------------------------


# REAL DATA ---------------------------------------------------------------
    
    # SIM
    n_groups <- 4
    
    q_data <- data_save$contact_matrix %>% 
      {. / rowSums(.)}
    
    q_tibble <- q_data %>% as_tibble() %>% 
      set_names(paste0("primary", 1:n_groups)) %>% 
      mutate(secondary_i_group = row_number()) %>% 
      pivot_longer(-secondary_i_group, names_to = "i_group", names_prefix = "primary", values_to = "q_estimate") %>% 
      mutate(i_group = as.integer(i_group))
    
    # CREATE example data using those probabilities
    # data_sim <- tibble(
    #   i_group = sample(1:n_groups, size = 10000, replace = TRUE, prob = data_save$group_props)
    # ) %>%
    #   group_by(i_group) %>%
    #   group_split() %>%
    #   map2(.y = q_tibble %>% group_by(i_group) %>% group_split() %>% map("q_estimate"),
    #        .f = ~ {
    #          .x %>% mutate(secondary_i_group = sample(1:n_groups, size = nrow(.), replace = TRUE, prob = .y))
    #        }) %>%
    #   bind_rows()
    

    # ACTUAL
    data_actual <- data_save$contact_matrix %>% as_tibble() %>% 
      set_names(1:4) %>% 
      mutate(from = row_number()) %>% 
      pivot_longer(-from, names_to = "to", values_to = "n") %>% 
      rowwise() %>% 
      mutate(n_list = list(1:n)) %>% 
      select(-n) %>% 
      unnest(n_list) %>% 
      select(i_group = from, secondary_i_group = to) %>% 
      mutate(i_group = as.integer(i_group),
             secondary_i_group = as.integer(secondary_i_group))
    
    
    
    
    
# Test new functions -------------------------------------------------------
    
    # mu_actual <- data_save$params_data$outside_home_contacts
    
    # TEMPORARY - comes from calibrate_mobility.R
    print("Contacts outside home:")
    print(contacts_outside_home)
    mu_actual <- contacts_outside_home
    
    
    # mu_example <- 10:6
    n_pop_total <- n_pop_for_k_calc
    
    # Calculate the best guess starting K matrix 
    # assume mu is 8 for everyone, and they randomly mix
    pops <- calc_pops(group_props = data_save$group_props, n_pop = n_pop_total)
    pop_matrix <- pops %>% {. %*% t(.)} %>% print
    beta_val <- pop_matrix %>% {mean(mu_actual) * pops / rowSums(.)} 
    k_start <- pop_matrix * beta_val[[1]]
    
    read_params(k_vec = k_matrix_to_vec(k_start),
                mu = mu_actual,
                group_props = data_save$group_props,
                n_pop = n_pop_total)
    
    # TRY MAX LIKELIHOOD AGAIN
    k_optimised <- nloptr::nloptr(
      x0 =  k_matrix_to_vec(k_start),
      eval_f = contact_maxlik_k_only, 
      lb = rep(0, length(k_matrix_to_vec(k_start))),  # entries to K can't go below 0 
      eval_g_ineq = k_constraint,  # constraint ensures that no entailed elements of K are negative
      opts = list(
        "algorithm" = "NLOPT_LN_COBYLA",
        "xtol_abs"= 0.1,
        "xtol_rel" = 1.0e-7,
        "maxeval" = Inf,
        "print_level" = 1
      ),
      n_pop = n_pop_total,
      group_props = data_save$group_props,
      mu = mu_actual,
      data = data_actual
    )
    
    
    k_optimised %>% print
    
    
    solution <- read_params(
      k_vec = k_optimised$solution,
      mu = mu_actual,
      group_props = data_save$group_props,
      n_pop = n_pop_for_k_calc
    )
    
    solution$q_matrix %>% round(3)
    solution$k_matrix %>% round(0)
    

# Export example K matrix -------------------------------------------------
    
    
    
    data_save$k_matrix <- solution$k_matrix
    
  
    round(data_save$k_matrix)
    # k_temp <- data_save$k_matrix
    # save(k_temp, file = "k_matrix_test_TEMP.RData")
    
    # data_save$k_matrix / k_temp
    # load("k_matrix_test_TEMP.RData")
    
    
    data_save$k_matrix_pop <- n_pop_for_k_calc
    
    k_scale_factor <- 1
    data_save$k_scale_factor <- k_scale_factor
    
    
    
    save(data_save, file = "data/processed/data_save.RData")
    