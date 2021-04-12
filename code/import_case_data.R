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
# source("~/Desktop/functions.R")

# devtools::install_github("dmbwebb/dups", dep = FALSE)
if (sum(str_detect(search(), "package:dups")) > 0)  detach("package:dups", unload=TRUE) # detach if already loaded
library(dups)

# devtools::install_github("dmbwebb/trackr", dep = FALSE)
if (sum(str_detect(search(), "package:trackr")) > 0)  detach("package:trackr", unload=TRUE) # detach if already loaded
library(trackr)



    library("haven"); library("lubridate")

    max_date <- ymd("2021-02-20") - days(14) # the data stops at 23rd Jan, so case data is reliable about 2 weeks before
    
    
    

# IMPORT ------------------------------------------------------------------

    sds_jobs <- read_dta(str_glue("data/casos_SDS_poblaciones_20Feb2021_slim.dta"))
    # sds_jobs %>% print_names
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
          "stratum_pop" = "poblacion_estrato",
          "recovered" = "recupsaluddatabog",
          "date_death" = "fecharecuperadosaluddata"
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
    # sds_jobs_with_stratum_impute %>% filter(is.na(stratum)) %>% view_n()


# CLEAN CASE DATA --------------------------------------------------

    # Stratum populations
    stratum_pops <- sds_jobs_with_stratum_impute %>% 
      # count_prop(stratum, stratum_pop)
      group_by(stratum, stratum_pop) %>% 
      summarise() %>% 
      filter(!is.na(stratum), !is.na(stratum_pop)) %>% 
      summarise(stratum_pop = sum(stratum_pop))  
    
    
    
    
    # Confirmed cases by group
    confirmed_cases_by_group <- sds_jobs_with_stratum_impute %>% 
      filter_track(!is.na(stratum)) %>% 
      group_by(stratum, date_results) %>% 
      summarise(n_confirmed_new = n()) %>% 
      full_join(
        crossing(stratum = 1:4, date_results = ymd("2020-03-01") + days(1:400)), by = c("stratum", "date_results")
      ) %>% 
      mutate(n_confirmed_new = if_else(is.na(n_confirmed_new), 0L, n_confirmed_new)) %>% 
      left_join(
        stratum_pops, by = "stratum"
      ) %>% 
      arrange(stratum, date_results) %>% 
      group_by(stratum, stratum_pop) %>% 
      mutate(
        n_confirmed_new_pc = n_confirmed_new / stratum_pop,
        n_confirmed_cum = cumsum(n_confirmed_new),
        n_confirmed_cum_pc = n_confirmed_cum / stratum_pop
      ) %>% 
      select(i_group = stratum, date_results, starts_with("n_confirmed")) %>% 
      filter(date_results <= max_date) %>% 
      ungroup
    
    
    # Confirmed cases by total
    confirmed_cases_total <- confirmed_cases_by_group %>% 
      group_by(date_results) %>% 
      summarise(across(c(stratum_pop, n_confirmed_new, n_confirmed_cum), sum)) %>% 
      rename(total_pop = stratum_pop) %>% 
      mutate(across(c(n_confirmed_new, n_confirmed_cum), list(pc = ~ .x / total_pop))) %>% 
      filter(date_results <= max_date)
    
    
    # ggplot(confirmed_cases_by_group, aes(x = date_results,
    #                                      y = rollmean(n_confirmed_new_pc, na.pad = TRUE, align = "right", k = 14),
    #                                      colour = factor(i_group))) +
    #   geom_line()
    # 
    
    
    # Infection fatality ratio
    # ifr <- 0.0034
    ifr <- tibble(i_group = 1:4, ifr = c(0.0024, 0.0027, 0.0042, 0.032))
    
    # Deaths
    deaths_by_group <- sds_jobs_with_stratum_impute %>% 
      filter_track(!is.na(stratum)) %>% 
      # count_prop(recovered)
      group_by(stratum, date_death) %>%
      summarise(n_deaths = sum(recovered == "Fallecido")) %>% 
      full_join(
        crossing(stratum = 1:4, date_death = ymd("2020-03-01") + days(1:400)), by = c("stratum", "date_death")
      ) %>% 
      mutate(n_deaths = if_else(is.na(n_deaths), 0L, n_deaths)) %>% 
      left_join(
        stratum_pops, by = "stratum"
      ) %>%
      arrange(stratum, date_death) %>% 
      group_by(stratum, stratum_pop) %>% 
      mutate(
        n_deaths_pc = n_deaths / stratum_pop,
        n_deaths_cum = cumsum(n_deaths),
        n_deaths_cum_pc = n_deaths_cum / stratum_pop
      ) %>% 
      select(i_group = stratum, date_death, starts_with("n_deaths")) %>% 
      left_join(
        ifr, by = "i_group",
      ) %>% 
      mutate(across(starts_with("n_deaths"), list(inferred = ~ .x / ifr), .names = "{.fn}_{.col}")) %>% 
      rename_with(~ str_replace(.x, "inferred_n_deaths", "n_inferred")) %>% 
      filter(date_death <= max_date) %>% 
      ungroup %>% 
      print
    
    
    
    

