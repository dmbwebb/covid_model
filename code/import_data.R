
# Setup -------------------------------------------------------------------

    # Packages
    library("haven")
    library("readxl")
    library("tidyverse")
    library("sn")
    library("nloptr")
    library("lubridate")


    # Functions
    hist_basic <- function(dat, x, binwidth = NULL, boundary = 0) {
      plot <- ggplot(dat, aes(x = {{x}})) + 
        geom_histogram(boundary = boundary, binwidth = binwidth, colour = "grey", fill = "lightblue")
      
      print(plot)
      
      invisible(dat)
    }
    
  

    
    mean_se_ci <- function(x, na.rm = TRUE) {
      
      
      reg <- lm(x ~ 1, data = data.frame(x = x))
      
      mean <- reg$coefficients[[1]]
      se <- summary(reg)$coefficients[, 2]
      conf_int <- stats::confint(reg, level = 0.95)
      
      return(tibble(
        stat = c("mean", "se", "conf_int_lower", "conf_int_upper"),
        value = c(mean, se, conf_int[[1]], conf_int[[2]])
      ))
      
    }
    
    f_p <- function(data, x) {
      f_fml <- as.formula(str_glue("{x} ~ i_group"))
      f_reg <- lm(f_fml, data = data)
      f_stats <- summary(f_reg)$fstatistic
      f_p <- pf(f_stats[1L], f_stats[2L], f_stats[3L], lower.tail = FALSE)
      return(f_p)
    }
    


    # List of data to be filled and saved to disk 
    data_save <- list()
    # load("data_save.RData")
    # data_save


# 1. HH size (Census data) -------------------------------------------------------------

    # Downloaded from
    # http://microdatos.dane.gov.co/index.php/catalog/643/data_dictionary
  
    # Directory for folder with Bogota census data
    dir_bogota_census <- "data/bogota census/11_Bogota_DTA"

    # Accommodation data (with stratum) - housing block
    # building <- read_dta(str_glue("{dir_bogota_census}/CNPV2018_1VIV_A2_11.DTA"))
    # 
    # building_cleaned <- building %>% 
    #   count_prop(l_tot_perl) %>% 
    #   count_prop(u_dpto) %>% 
    #   select(building_id = cod_encuestas,
    #          n_hh_in_housing = v_tot_hog,
    #          stratum = va1_estrato)
    
    # write_dta(building_cleaned, str_glue("{dir_bogota_census}/census_buildings.dta"))

    building_cleaned <- read_dta(str_glue("{dir_bogota_census}/census_buildings.dta"))
    
    
    # Household
    # households <- read_dta(str_glue("{dir_bogota_census}/CNPV2018_2HOG_A2_11.DTA"))
    # # 
    # households_cleaned <- households %>%
    #   select(building_id = cod_encuestas,
    #          hh_id = h_nrohog,
    #          n_rooms = h_nro_cuartos,
    #          n_bedrooms = h_nro_dormit,
    #          hh_size = ha_tot_per)
  
    # write_dta(households_cleaned %>% filter(row_number() <= 1e6), str_glue("{dir_bogota_census}/census_households_1.dta"))
    # write_dta(households_cleaned %>% filter(row_number() > 1e6), str_glue("{dir_bogota_census}/census_households_2.dta"))
    
    households_cleaned <- bind_rows(
      read_dta(str_glue("{dir_bogota_census}/census_households_1.dta")),
      read_dta(str_glue("{dir_bogota_census}/census_households_2.dta"))
    )
      
    # households_cleaned %>% dups_report(building_id)        # multiple hh per building
    # households_cleaned %>% dups_report(building_id, hh_id) # no dups
    
    
    # Individuals
    ind <- read_dta(str_glue("{dir_bogota_census}/CNPV2018_5PER_A2_11.DTA"))
    
    # Location (to restrict to metropolitan area)
    # geo <- read_dta(str_glue("{dir_bogota_census}/CNPV2018_MGN_A2_11.DTA"))
    
    ind_cleaned <- ind %>% 
      select(building_id = cod_encuestas,
             hh_id = p_nrohog,
             age_group = p_edadr,
             person_id = p_nro_per) %>% 
      mutate(age_lower = (age_group - 1) * 5)
      # dups_report(building_id, hh_id, person_id) this is the unique level
    
  
    
    
    # Merge building and household data
    hh_data_merged <- inner_join_track(building_cleaned, households_cleaned, by = "building_id")
    
    hh_data_cleaned <- hh_data_merged %>% 
      mutate(
        building_id = as.character(as.integer(building_id)),
        hh_id = as.character(as.integer(hh_id)),
        hh_id = paste0(building_id, hh_id)
      ) %>% 
      select(-building_id) %>% 
      relocate(hh_id) %>% 
      count_prop(stratum) %>% 
      filter_track(stratum %in% 1:6) %>% 
      mutate(stratum = factor(stratum))

    # Merge with ind
    ind_with_hh <- left_join_track(ind_cleaned, hh_data_merged, by = c("building_id", "hh_id"))
    
    ind_with_hh %>% 
      filter_track(stratum %in% 1:6) %>% 
      mutate(stratum = as.integer(stratum)) %>% 
      mutate(i_group = case_when(stratum %in% c(1, 2) ~ 1L,
                                 stratum %in% c(3, 4) ~ stratum - 1L,
                                 stratum %in% c(5, 6) ~ 4L)) %>% 
      ggplot(aes(x = age_lower + 0.1, fill = factor(i_group)), colour = "white") + 
      geom_histogram(aes(y=..density..), binwidth = 10) + 
      facet_wrap(~ factor(i_group))
    
    
    ages <- ind_with_hh %>% 
      filter_track(stratum %in% 1:6) %>% 
      mutate(stratum = as.integer(stratum)) %>% 
      mutate(i_group = case_when(stratum %in% c(1, 2) ~ 1L,
                                 stratum %in% c(3, 4) ~ stratum - 1L,
                                 stratum %in% c(5, 6) ~ 4L)) %>% 
      mutate(age_chunks = cut(age_lower, c(-1, 14, 29, 49, 69, 1000), labels = c("0-14", "15-29", "30-49", "50-69", "70+"))) %>% 
      # count_prop(age_lower, age_chunks)
      group_by(i_group, age_chunks) %>%
      summarise(n = n()) %>% 
      group_by(i_group) %>% 
      mutate(prop = n / sum(n))
    
    age_table <- ages %>% 
      mutate(i_group = factor(i_group)) %>% 
      mutate(i_group = recode_i_group(i_group)) %>% 
      select(-n) %>% 
      # mutate(prop = format(round(prop * 100, 1), nsmall = 1)) %>% 
      # mutate(prop = str_glue("{prop}%")) %>% 
      pivot_wider(names_from = i_group, values_from = prop) %>% 
      print %>% 
      rename(Age = age_chunks) %>% 
      mutate(Age = as.character(Age))
  
    total_perc_row <- c(
      "Total:",
      rep("100%", 4)
    )
    
    total_row <- c(
      "Observations:",
      ages %>% group_by(i_group) %>% summarise(n = sum(n)) %>% .$n
    )
    
    age_table_totals <- age_table %>% 
      rbind(total_perc_row) %>% 
      rbind(total_row)
    
    write_excel_csv(age_table_totals, path = "data/processed/age_table.csv")
    
    # library(xtable)
    # print(xtable(age_table_totals), include.rownames = FALSE, booktabs = TRUE)
    # ?print.xtable
    
    
    ages %>% 
      ggplot(aes(x = age_chunks, fill = factor(stratum))) + 
      geom_col(aes(y = prop), colour = "white") + 
      facet_wrap(~ factor(stratum))
    
    
    
    
    # PLOT THE CDFs of housing size
    hh_data_cleaned %>%
      ggplot(aes(x = hh_size, colour = factor(stratum))) +
      stat_ecdf(aes(colour = stratum), geom = "point") + 
      stat_ecdf(aes(colour = stratum), geom = "step", alpha = 0.4) + 
      facet_wrap(~ stratum) + 
      coord_cartesian(xlim = c(0, 6))
    
    # What are the mean hh_size, and size of each group
    hh_sizes <- hh_data_cleaned %>% group_by(stratum) %>% 
      summarise(mean_hh_size = mean(hh_size, na.rm = TRUE),
                group_n = sum(hh_size, na.rm = TRUE)) %>% 
      mutate(group_prop = group_n / sum(group_n))
    
    # EXPORT MEAN/SE/CI of hh size
    hh_data_export_mean_se <- hh_data_cleaned %>% 
      mutate(i_group = as.integer(stratum)) %>% 
      mutate(i_group = case_when(i_group %in% c(1, 2) ~ 1L,
                                 i_group %in% c(3, 4) ~ i_group - 1L,
                                 i_group %in% c(5, 6) ~ 4L)) %>% 
      mutate(i_group = recode_i_group(i_group))
    
    hh_data_export_mean_se %>% group_by(i_group) %>% summarise(mean_se_ci(hh_size)) %>% write_csv(path = "data/temp/hh_size.csv")
    hh_data_export_mean_se %>% summarise(mean_se_ci(hh_size)) %>% write_csv(path = "data/temp/hh_size_all.csv")
    
    hh_data_export_mean_se %>% f_p(x = "hh_size")
    
    group_props_bogota <- hh_sizes$group_prop
    n_pop_bogota <- sum(hh_sizes$group_n)
    
    # Measures of overcrowding
    hh_data_cleaned %>%
      mutate(people_per_room = hh_size / n_rooms,
             people_per_bedroom = hh_size / n_bedrooms) %>% 
      group_by(stratum) %>% 
      summarise(across(c(people_per_room, people_per_bedroom), mean))
    
    
    save(hh_data_cleaned, file = "data/processed/hh_data_cleaned.RData")
    
    
    
    # Create long version (at individual level)
    hh_data_long <- hh_data_cleaned %>% 
      group_by(hh_id, stratum, hh_size) %>% 
      summarise(
        hh_ind_id = 1:hh_size
      ) %>% 
      mutate(
        hh_ind_id = paste0(hh_id, "_", hh_ind_id)
      ) %>% 
      ungroup %>% 
      rename(i_group = stratum)
    
    hh_data_bogota <- hh_data_long %>% 
      mutate(i_group = as.integer(i_group), 
             hh_id = as.numeric(hh_id),
             i_group = case_when(i_group %in% c(1, 2) ~ 1L,
                                 i_group %in% c(3, 4) ~ i_group - 1L,
                                 i_group %in% c(5, 6) ~ 4L))
      
    # Add onto data_save
    data_save$hh_data_bogota <- hh_data_bogota
    
    

    # 
    # # OUTPUT A GRAPH of HH SIZE
    # hh_size_means <- data_save$hh_data_bogota %>% 
    #   select(-hh_ind_id) %>% 
    #   dups_drop(hh_id) %>% 
    #   group_by(i_group) %>% 
    #   summarise(mean_se(hh_size))
    # 
    # 
    # hh_size_means %>% 
    #   ggplot(aes(x = i_group, y = y, colour = factor(i_group))) + 
    #   geom_point() + 
    #   geom_errorbar(aes(ymin = ymin, ymax = ymax), width = 0.3) + 
    #   labs(x = "SES group", y = "Estimate") + 
    #   scale_fill_viridis(option = "viridis",
    #                      begin = 0.3, 
    #                      discrete = TRUE)
    
    
    
      

# 1b. Group size --------------------------------------------------------------

    group_props <- data_save$hh_data_bogota %>% 
      group_by(i_group) %>% summarise(n = n()) %>% 
      mutate(group_props = n / sum(n)) %>% 
      pull(group_props)
    
    data_save$group_props <- group_props
    
    
    
    
# 2. Import test sensitivity data --------------------------------------------
    
    test_sensitivity_data <- read_csv(file = "data/grassly_test_sensitivity.txt",
                                      col_types = cols(day_exposure = col_integer(),
                                                       sens = col_double())) %>% 
      rename(days_exposure = day_exposure, sensitivity = sens) %>% 
      bind_rows(
        tibble(
          days_exposure = c(0, 36:50),
          sensitivity = 0
        )
      ) %>% 
      arrange(days_exposure)
    
    
    # Add to data_save
    data_save$test_sensitivity_data <- test_sensitivity_data
    
    

    
    
    



    
    




    
    


    
    



    
    
    
    
    
    
    
# 3. Test delay data ---------------------------------------------------------
    
    
    # Import SDS delay data
    # dir_sds <- "/Users/user/Dropbox/covid-project/data/SDS"
    
    # WITH DATES
    # sds_old <- read_excel(str_glue("{dir_sds}/processed/casos_2020_08_06_cleaned_dates.xlsx"),
    #                       col_types = c(rep("guess", 7),
    #                                     "date",
    #                                     rep("guess", 25)))
    # sds <- read_excel(str_glue("{dir_sds}/originals/casos_2020_09_18.xlsx"))
    # n_max = 200
    # col_types = c("numeric", "text", "text", 
    #               "text", "text", "text", "text", "date", 
    #               "date", "date", "date", "text", "text", 
    #               "numeric", "text", "text", "text", 
    #               "text", "text", "text", "text", "text", 
    #               "text", "text", "text", "text", "text", 
    #               "text", "text", "text", "text", "text", 
    #               "numeric")
    # ?read_excel


    
    # sds_old %>% count_prop(`ESTRATO SOCIOECONOMICO`)
    
    # COUNT HOW MANY MISSING STRATA
    # sds %>% 
    #   rename_with(janitor::make_clean_names) %>% 
    #   select(starts_with("estrato")) %>% 
    #   mutate(across(starts_with("estrato"), ~ if_else(as.numeric(.) %in% 1:6, as.numeric(.), NA_real_))) %>% 
    #   count_nas()
    
    library("lubridate")
    library("data.table")
  
    sds_jobs <- read_dta(str_glue("data/casos_SDS_poblaciones_20Feb2021_slim.dta"))
    
    # sds_jobs %>%
      # select(caso, fechainiciosintomas, fecha_consulta, fechadiagnostresultlaboratorio, estratosocioeconomico, estrato, poblacion_estrato, recupsaluddatabog, fecharecuperadosaluddata, upzgeo) %>%
      # write_dta("data/casos_SDS_poblaciones_20Feb2021_slim.dta")

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
      
    
    
    
    # Calculate delays
    sds_delays <- sds_jobs_with_stratum_impute %>% 
      mutate(consultation_delay = interval_days(date_symptoms, date_consultation),
             results_delay = interval_days(date_consultation, date_results)) %>% 
      count_prop(consultation_delay, !is.na(date_symptoms)) %>% 
      count_prop(stratum) %>% 
      filter_track(stratum %in% 1:6) %>%     # 35% is unknown (0) !!
      # count_prop()
      count_prop(consultation_delay == 0) %>%   # lots of unkown
      mutate(consultation_delay = if_else(consultation_delay <= 0, NA_real_, consultation_delay),
             results_delay = if_else(results_delay < 0, NA_real_, results_delay)) %>% 
      mutate(total_delay = consultation_delay + results_delay)
    # filter_track(consultation_delay > 0)     # removes asympotmatic people
  
    
    
    
    # PLOT DISTRIBUTIONS
    # sds_delays %>% 
    #   filter(stratum %in% 1:6) %>% 
    #   # filter(total_delay %between% c(quantile(total_delay, 0.02, na.rm = TRUE), quantile(total_delay, 0.98, na.rm = TRUE))) %>% 
    #   ggplot(aes(x = total_delay)) + geom_density(aes(colour = factor(stratum)), adjust = 2) + 
    #   coord_cartesian(xlim = c(0, 25))
    # 
    # sds_delays %>% 
    #   filter(stratum %in% 1:6) %>% 
    #   filter(consultation_delay != 0) %>% 
    #   # filter(results_delay %between% c(quantile(total_delay, 0.02, na.rm = TRUE), quantile(total_delay, 0.98, na.rm = TRUE))) %>% 
    #   ggplot(aes(x = results_delay)) + geom_density(aes(colour = factor(stratum)), adjust = 2) + 
    #   coord_cartesian(xlim = c(0, 25))
    # 
    # sds_delays %>% 
    #   filter(stratum %in% 1:6) %>% 
    #   filter(consultation_delay != 0) %>% 
    #   ggplot(aes(x = total_delay)) + 
    #   geom_density(aes(colour = factor(stratum))) + 
    #   # geom_histogram(aes(y = ..density.., fill = factor(stratum)), binwidth = 0.5, boundary = -0.25, position = "dodge") + 
    #   coord_cartesian(xlim = c(0, 25))
    # 
    # sds_delays %>% 
    #   filter(stratum %in% 1:6) %>% 
    #   # filter(consultation_delay > 0) %>%
    #   ggplot(aes(x = consultation_delay)) + 
    #   geom_histogram(aes(y = ..density.., fill = factor(stratum)), colour = "white", size = 0.2, binwidth = 0.5, boundary = -0.25, position = "dodge")
    # # coord_cartesian(xlim = c(2, 10))
    
    
    
    # plot_density(sds_delays, x = consultation_delay, y = results_delay)
    
    # ggplot(sds_delays, aes(x = consultation_delay, y = results_delay)) + 
    #   geom_smooth()
    
    
    
    # EXPORT VALUES FOR THE TABLE OF PARAMETERS
    sds_delays_for_table <- sds_delays %>% 
      select(i_group = stratum, 
             test_choice_delay = consultation_delay,
             test_results_delay = results_delay) %>% 
      mutate(i_group = recode_i_group(i_group)) %>% 
      group_by(i_group)
    
    
    # se <- function(x, na.rm = FALSE) {
    #   if (na.rm) x <- na.omit(x)
    #   sqrt(var(x)/length(x))
    # }
    
  
  
 
    
    
    
    sds_delays_for_table %>% group_by(i_group) %>% summarise(mean_se_ci(test_choice_delay)) %>% write_csv(path = "data/temp/test_choice_delay.csv")
    sds_delays_for_table %>% ungroup %>% summarise(mean_se_ci(test_choice_delay)) %>% write_csv(path = "data/temp/test_choice_delay_all.csv")
    
    sds_delays_for_table %>% group_by(i_group) %>% summarise(mean_se_ci(test_results_delay)) %>% write_csv(path = "data/temp/test_results_delay.csv")
    sds_delays_for_table %>% ungroup %>% summarise(mean_se_ci(test_results_delay)) %>% write_csv(path = "data/temp/test_results_delay_all.csv")
    

    sds_delays_for_table %>% f_p(x = "test_choice_delay")
    sds_delays_for_table %>% f_p(x = "test_results_delay")
    
    # summary(lm(test_results_delay ~ i_group, data = sds_delays_for_table))
    
    # mean_se_ci_f(sds_delays_for_table, x = "test_choice_delay")
    
    
    # Export to data_save
    
    testing_delays_data <- sds_delays %>% 
      select(i_group = stratum, 
             test_choice_delay = consultation_delay,
             test_results_delay = results_delay) %>% 
      # mutate(i_group = case_when(i_group %in% 1:2 ~ 1,
      #                            i_group %in% 3:4 ~ i_group - 1,
      #                            i_group %in% 5:6 ~ 4)) %>% 
      mutate(test_choice_delay = if_else(test_choice_delay == 0, 0.05, test_choice_delay),
             test_results_delay = if_else(test_results_delay == 0, 0.05, test_results_delay)) # add a small delay to the 0s [avoids errors in the model]
  
    
    testing_delays_data %>% count_prop(value(test_choice_delay), value(test_results_delay))
    
    testing_delays_data %$% cor(test_choice_delay, test_results_delay, use = "pairwise.complete.obs") # basically no correlation
    
    data_save$test_delay_data <- testing_delays_data %>% drop_na()
    
   
    
    
# 4. Tracing delays (from testing delays) ---------------------------------------------------------------
    
    
    
    ct_delay_data <- tibble(
      i_group = c(rep(1, 1000), rep(2, 1000), rep(3, 1000), rep(4, 1000)),
      ct_delay = sample(testing_delays_data$test_choice_delay[!is.na(testing_delays_data$test_choice_delay)], 
                        4000, replace = TRUE)
    )
    
    
    data_save$ct_delay_data <- ct_delay_data
    
    
    # Contact tracing delay
    # ct_delay_data_real <- tibble(
    #   i_group = c(rep(1, 1000), rep(2, 1000), rep(3, 1000)),
    #   ct_delay = c(runif(1000, 5, 7), runif(1000, 3, 5), runif(1000, 1, 3))
    # )
    
    # ct_delay_data_real_6 <- tibble(
    #   i_group = c(rep(1, 1000), rep(2, 1000), rep(3, 1000), rep(4, 1000), rep(5, 1000), rep(6, 1000)),
    #   ct_delay = c(runif(1000, 5, 7), runif(1000, 3, 5), runif(1000, 1, 3), runif(1000, 1, 3), runif(1000, 1, 3), runif(1000, 1, 3))
    # )
    
    
    
# 5. Import calculated parameters (from Rachid) -----------------------------------------------------------
    
    library("readxl")
    
    param_file <- "data/params_final.xlsx"
    
    params_data <- read_excel(param_file, n_max = 5)
    
    contact_matrix <- read_excel(param_file, sheet = "contact_matrix") %>% 
      select(-1) %>% 
      as.matrix()
    
    
    # CALCULATE FULL PROBABILITY TIBBLES
    probs_basic_4 <- crossing(
      i_group = 1L:4L, 
      symptom_severity = 1L:2L
    )
    
  
    # Self test
    probs_self_test <- probs_basic_4 %>% 
      mutate(prob = rep(params_data$prob_self_test, each = 2),
             prob = if_else(symptom_severity == 1, 0, prob)) # asymptomatics (severity = 1) have 0 probability of isolating
    
    
    # Isolate Test
    probs_isolate_test <- probs_basic_4 %>% 
      mutate(prob = rep(params_data$p_isolate_test, each = 2)) # independent of symptom severity
    
    
  
    # Isolate CT
    p_isolate_ct <- 1 - (params_data$days_work_ct / params_data$days_of_work)
    p_isolate_ct <- if_else(p_isolate_ct < 0, 0, p_isolate_ct)
    
    
    probs_isolate_ct <- probs_basic_4 %>% 
      mutate(prob = rep(p_isolate_ct, each = 2))   # independent of symptom severity
    
    
    # Isolate symptoms
    p_isolate_symp <- 1 - (params_data$days_work_symptoms / params_data$days_of_work)
    
    probs_isolate_symptoms <- probs_basic_4 %>% 
      mutate(prob = rep(p_isolate_symp, each = 2),
             prob = if_else(symptom_severity == 1, 0, prob)) # asymptomatics (severity = 1) have 0 probability of isolating
    
    
    # P hh quarantine
    p_hh_quarantine <- 1 - (params_data$days_work_hh_quarantine / params_data$days_of_work)
    
    
    
    # Save SAR and contact matrix
    data_save$sar_out_input <- params_data$sar_out_input
    data_save$sar_home_input <- params_data$sar_home_input
    data_save$contact_matrix <- contact_matrix    
    
    # Save probs
    data_save$p_isolate_ct <- p_isolate_ct
    data_save$p_isolate_symp <- p_isolate_symp
    data_save$p_hh_quarantine <- p_hh_quarantine
    
    # Save probs DF
    data_save$probs_self_test <- probs_self_test
    data_save$probs_isolate_symptoms <- probs_isolate_symptoms
    data_save$probs_isolate_test <- probs_isolate_test
    data_save$probs_isolate_ct <- probs_isolate_ct
    
    # Save days of work (for recalculating isolation in counterfactuals)
    data_save$days_of_work <-       params_data$days_of_work
    data_save$days_work_ct <-       params_data$days_work_ct
    data_save$days_work_symptoms <- params_data$days_work_symptoms
    data_save$days_work_hh_quarantine <- params_data$days_work_hh_quarantine
    
    data_save$params_data <- params_data
    
    
# 6. Starting prevalence -----------------------------------------------------
    
    # sds_jobs_with_stratum_impute %>% 
    #   # print_names %>% 
    #   mutate(stratum = as.integer(stratum)) %>% 
    #   mutate(stratum = case_when(stratum %in% c(1, 2) ~ 1L,
    #                              stratum %in% c(3, 4) ~ stratum - 1L,
    #                              stratum %in% c(5, 6) ~ 4L)) %>% 
    #   arrange(date_results, stratum) %>% 
    #   group_by(date_results) %>% 
    #   # group_by(date_results, stratum) %>%
    #   summarise(n = n()) %>% 
    #   mutate(cum_n = cumsum(n)) %>% 
    #   print(n = 100)
    
    # Start date at April 1st - 489 confirmed cases
    # start_date <- ymd("2020-04-01")
    
    
    
    # best_guess_initial_cases <- sds_jobs_with_stratum_impute %>% 
    #   mutate(stratum = as.integer(stratum)) %>% 
    #   # mutate(stratum = case_when(stratum %in% c(1, 2) ~ 1L,
    #   #                            stratum %in% c(3, 4) ~ stratum - 1L,
    #   #                            stratum %in% c(5, 6) ~ 4L)) %>% 
    #   group_by(date_results, stratum, stratum_pop) %>%
    #   summarise(n = n()) %>% 
    #   ungroup %>% 
    #   filter(date_results <= start_date) %>% 
    #   group_by(stratum, stratum_pop) %>% 
    #   filter(!is.na(stratum)) %>% 
    #   summarise(n = sum(n)) %>% 
    #   group_by(stratum) %>% 
    #   summarise(stratum_pop = sum(stratum_pop),
    #             n = sum(n)) %>% 
    #   mutate(start_date = start_date)
    # 
    # 
    # data_save$best_guess_initial_cases <- best_guess_initial_cases
    # 
    # data_save$start_date <- start_date
    

# 7. Lockdown start / end, mobility -------------------------------------------------

    # Import
    mobility_data <- read_csv("data/bogota_mobility_report.csv") %>% 
      select(date:residential_percent_change_from_baseline) %>% 
      mutate(date = dmy(date)) %>% 
      pivot_longer(-date, names_to = "category", values_to = "perc_change") %>% 
      mutate(category = str_replace_all(category, "_percent_change_from_baseline", ""))
    
    
    # Quick plot
    mobility_data %>% 
      arrange(category, date) %>% 
      group_by(category) %>% 
      mutate(perc_change_roll = zoo::rollmean(perc_change, k = 7, na.pad = TRUE, align = "left")) %>% 
      ggplot(aes(x = date, y = perc_change_roll, colour = category)) + 
      scale_x_date(minor_breaks = "months") + 
      geom_line() + 
      facet_wrap(~ category)
    
    
    end_severe <- ymd("2020-06-01")
    end_weak <- ymd("2020-09-01")
    # lockdown_end_date <- ymd("2020-09-01")
    
    # Take the rolling average (only workplaces)
    mobility_summ <- mobility_data %>% 
      filter(category == "workplaces") %>% 
      mutate(perc_change_roll = zoo::rollmean(perc_change, k = 7, na.pad = TRUE, align = "left")) %>% 
      mutate(phase = case_when(
        date >= ymd("2020-04-01") & date < end_weak ~ "1_severe_lockdown",
        # date >= end_severe        & date < end_weak   ~ "2_lighter_lockdown",
        date >= end_weak                              ~ "3_loose"
      )) %>% 
      group_by(phase) %>% 
      mutate(phase_average = mean(perc_change))

    # Assume that in early april it's 
    
    
   
    
    ggplot(mobility_summ, aes(x = date, y = perc_change_roll, colour = category)) + 
      geom_line() + 
      geom_line(aes(y = phase_average), colour = "skyblue") + 
      # geom_vline(xintercept = start_date) +
      geom_vline(xintercept = end_severe) + 
      geom_vline(xintercept = end_weak) + 
      scale_x_date(breaks = "month") + 
      labs(y = "Mobility - Perc change from baseline in early 2020")
      # geom_vline(xintercept = ymd("2020-09-01")) # 1st sept
    

    mobility_increase_factor <- mobility_summ %>% select(phase, phase_average) %>% 
      ungroup %>% 
      dups_drop() %>% 
      filter(!is.na(phase)) %>% 
      mutate(pp = 1 + (phase_average / 100)) %>% 
      mutate(factor = pp /  pp[phase == "1_severe_lockdown"]) %>% 
      .$factor
    
    
    # data_save$lockdown_end_date <- c(end_severe, end_weak)
    # data_save$lockdown_end_mobility_increase <- mobility_increase_factor[2:3]
    
    data_save$lockdown_end_date <- end_weak
    data_save$lockdown_end_mobility_increase <- mobility_increase_factor[[2]]
    
    
    

# 8. Contact dispersion ---------------------------------------------------

    # From Bi et al 2020
    data_save$contact_dispersion <- 0.58
    
    

# 9. Infectiousness by symptoms -------------------------------------------

    # From https://www.ncbi.nlm.nih.gov/pmc/articles/PMC7508369/
    data_save$infectiouness_by_symptoms <- c(0.35, 1)
    


# ..... -------------------------------------------------------------------


    
# EXPORT TO DISK ----------------------------------------------------------
    
    save(data_save, file = "data/processed/data_save.RData")
    
    