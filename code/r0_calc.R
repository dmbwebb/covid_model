
    library("haven")
    library("lubridate")
    
    load("data/processed/data_save.RData")
    load("data/processed/gen_data_clean_timings.RData")
    
# IMPORT ------------------------------------------------------------------
    
    sds_jobs <- read_dta(str_glue("data/casos_SDS_poblaciones_20Feb2021_slim.dta"))
    
    # Clean up
    sds_jobs_clean <- sds_jobs %>% 
      # sample_n(100) %>% 
      rename(
        any_of(c(
          "case_id" = "caso",
          "date_symptoms" = "fechainiciosintomas",
          "date_consultation" = "fecha_consulta",
          "date_results" = "fechadiagnostresultlaboratorio",
          "stratum" = "estratosocioeconomico",
          "stratum_2" = "estrato",
          "stratum_pop" = "poblacion_estrato"
        ))
      ) %>% 
      mutate(stratum = coalesce(as.numeric(stratum), stratum_2)) %>% 
      select(-stratum_2) %>% 
      count_prop(stratum) %>% 
      mutate_track(across(stratum, ~ if_else(as.numeric(.x) %in% 1:6, as.numeric(.x), NA_real_))) %>% 
      mutate(across(starts_with("date"), dmy)) %>% 
      count_prop(stratum, stratum_pop) %>% 
      count_prop(year(date_symptoms), year(date_results), year(date_consultation))
    
    
    # IMPUTE missing values using UPZ - check this
    stratum_upz <- sds_jobs_clean %>% 
      filter(str_detect(upzgeo, "UPZ")) %>% 
      mutate(stratum = case_when(stratum %in% c(1, 2) ~ 1,
                                 stratum %in% c(3, 4) ~ stratum - 1,
                                 stratum %in% c(5, 6) ~ 4)) %>% 
      group_by(upzgeo, stratum) %>% 
      summarise(n = n()) %>% 
      arrange(upzgeo) %>% 
      filter(!is.na(stratum)) %>% 
      mutate(p = n / sum(n)) %>% 
      filter(sum(n) > 0) %>%
      filter(p == max(p)) %>% 
      ungroup
    
    ggplot(stratum_upz, aes(x = p)) + geom_histogram(fill = "skyblue", colour = "darkblue", binwidth = 0.05, boundary = 0)
    
    stratum_impute_upz <- stratum_upz %>% 
      select(upzgeo, stratum_to_impute = stratum) %>% 
      dups_report(upzgeo, stratum_to_impute) # no dups
    
    
    sds_jobs_with_stratum_impute <- sds_jobs_clean %>% 
      left_join(stratum_impute_upz, by = "upzgeo") %>% 
      mutate(stratum = as.integer(stratum)) %>% 
      mutate(stratum = case_when(stratum %in% c(1, 2) ~ 1L,
                                 stratum %in% c(3, 4) ~ stratum - 1L,
                                 stratum %in% c(5, 6) ~ 4L)) %>% 
      mutate(stratum_imputed_yn = is.na(stratum) & !is.na(stratum_to_impute),
             stratum = if_else(is.na(stratum), as.integer(stratum_to_impute), stratum)) %>% 
      count_prop(stratum)
    
    
    sds_jobs_with_stratum_impute
    
# Exponential growth calc -------------------------------------------------

    
    
    overall_prev <- sds_jobs_with_stratum_impute %>% 
      select(case_id, date_results) %>% 
      group_by(date_results) %>% 
      summarise(new_cases = n()) %>% 
      full_join(tibble(date_results = ymd("2020-03-01") + days(0:300))) %>% 
      arrange(date_results) %>% 
      mutate(new_cases = if_else(is.na(new_cases), 0L, new_cases)) %>% print(n = 500) %>% 
      mutate(days_since_march_1 = interval_days("2020-03-01", date_results)) %>% 
      mutate(ln_new_cases = log(new_cases))
    
    
    # min_date <- ymd("2020-03-01") + days(30)
    min_date <- ymd("2020-04-01")
    max_date <- ymd("2020-04-01") + 60
      days(interval_days(data_save$start_date, data_save$lockdown_end_date[[1]]))
    
   
    
    
    # REGRESSION
    reg <- lm(ln_new_cases ~ days_since_march_1, data = overall_prev, 
              subset = date_results <= max_date & date_results >= min_date)
    
    # Coefficients
    r <- reg$coefficients[[2]]
    r_conf <- confint(reg)[2, ]
    
    (r_lab <- str_glue("r = {round(r, 3)}\n[{round(r_conf[[1]], 3)}, {round(r_conf[[2]], 3)}]"))
    
    # PLOT
    overall_prev %>% 
      filter(date_results >= ymd("2020-03-14")) %>% 
      filter(date_results <= ymd("2020-09-01")) %>% 
      filter(new_cases != 0) %>% 
      
      ggplot(aes(x = date_results, y = ln_new_cases))  + 
      geom_point(alpha = 0.2) + 
      geom_smooth(data = overall_prev %>% filter(date_results <= max_date & date_results >= min_date), method = "lm", fill = "skyblue", colour = "indianred", alpha = 0.4) + 
      # theme_custom(panel.grid = element_blank()) + 
      scale_x_date(breaks = "months", date_labels = "%b%est") + 
      labs(x = "Date (2020)", y = "ln(Daily Confirmed New Cases)") + 
      annotate(geom = "text", x = ymd("2020-04-01"), y = 7.5, label = r_lab, vjust = "inward", hjust = "inward") + 
      theme_custom(panel.grid = element_blank())
      # geom_text(label = r_lab)
    
    save_plot("r0_calibrate_data")
    
    
    
    
    
    
    

# Import generation interval dist -----------------------------------------

    
    
    ggplot(gen_data_clean, aes(x = secondary_timing)) + 
      geom_density(colour = "seagreen4", fill = "seagreen4", alpha = 0.5) + 
      scale_y_continuous(labels = scales::percent_format(accuracy = 1), minor_breaks = NULL) + 
      labs(y = "Density", x = "Generation Interval (days)") + 
      coord_cartesian(xlim = c(0, 20)) +
      theme_classic()
    
    

    dens <- density(gen_data_clean$secondary_timing)
    
    
    gen_int <- tibble(x = dens$x, y = dens$y) %>% 
      filter(x > 0) 
    
  
    ggplot(gen_int, aes(x = x, y = y)) + geom_line()
    

    
    
    

# Combine with exponential -----------------------------------------------

    
    r <- reg$coefficients[[2]]
    r_std_error <- summary(reg)$coefficients[2, 2]
    # r <- confint(reg)[2, ][[1]]
    
    integrand_df <- gen_int %>% 
      rename(t = x, g_x = y) %>% 
      mutate(exp_term = exp(-r * t),
             integrand = exp_term * g_x)
    
    
    
    
    # r0_stores$r0_est

# Numerically integrate ---------------------------------------------------

    # Create a trapezoid function that approximates the integrand
    integrand_trapezoid <- approxfun(integrand_df$t,
                                     integrand_df$integrand)
    
    
    r0_inv <- integrate(integrand_trapezoid, min(integrand_df$t), max(integrand_df$t))$value

    r0 <- 1 / r0_inv
    
    print(str_glue("R0 estimate from data is {round(r0, 3)}"))
    
    data_save$initial_r0 <- r0
    
    
    
    
# Quantify uncertainty ----------------------------------------------------
    
    # Function to go from r to R0 [has external dependencies]
    r_to_r0 <- function(r) {
      integrand_df <- gen_int %>% 
        rename(t = x, g_x = y) %>% 
        mutate(exp_term = exp(-r * t),
               integrand = exp_term * g_x)
      integrand_trapezoid <- approxfun(integrand_df$t,
                                       integrand_df$integrand)
      r0_inv <- integrate(integrand_trapezoid, min(integrand_df$t), max(integrand_df$t))$value
      r0 <- 1 / r0_inv
      return(r0)
    }
    
    
    # Assume that estimator of r is normally distributed with mean r and variance = std.error ^2 
    r_dist <- tibble(
      r = rnorm(n = 10000, mean = r, sd = r_std_error)
    ) %>% 
      hist_basic(x = r)
    
    # Takes about 2 mins
    r0_dist <- r_dist %>% 
      mutate(r0 = map2_dbl(r, row_number(), ~ {
        if (.y %% 100 == 0) print(.y)
        return(r_to_r0(.x))
      }))
    
    r0_dist %>% hist_basic(x = r0)
    
    # Calculate confidence interval
    round(quantile(r0_dist$r0, c(0.025, 0.5, 0.975)), 3)
    

    