    library("ggpubr")
    
    load("data/processed/gen_data_clean_timings.RData")


    # test_run <- group_diff_mobility_change(n_sims = 2, n_cores = 1, n_iterations = 10, n_initial_cases = c(1, 1, 1, 1))
    # params_graph_df <- test_run$params
    # params_graph_df <- policy_change$outbreak_t[[1]]$params
    # l <- formals(group_diff_mobility_change)
    
    get_defaults <- function(.f) {
      .f <- group_diff
      list_formals <- formals(.f)
      
      list_vals <- list_formals %>% 
        discard(is.symbol) %>% 
        discard(is.null) %>% 
        map_if(.p = ~ is.language(.x), .f = eval)
      
      return(list_vals)
    }
    
    params_graph_df <- get_defaults(group_diff)
    
    

# Functions for plotting --------------------------------------------------


    
    param_bar_plot <- function(param_vec, y_lab) {
      tibble(i_group = 1:4, group_size = param_vec) %>% 
        mutate(i_group = recode_i_group(i_group)) %>% 
        ggplot(aes(x = i_group, y = group_size, fill = i_group)) + 
        geom_col() + 
        labs(x = "SES group", y = y_lab) + 
        scale_fill_viridis(option = "viridis",
                           begin = 0.3, 
                           discrete = TRUE, guide = FALSE) + 
        theme_custom()
    }
    
    
    param_box_plot <- function(data, y, y_lab = NULL) {
      
      # Calculate limits with no outliers
      ylim_no_outliers <- data %>% pull({{y}}) %>% boxplot.stats() %>% .$stats %>% .[c(1, 5)]
      
      data %>% 
        mutate(i_group = recode_i_group(i_group)) %>%
        ggplot(aes(x = factor(i_group), y = {{y}}, fill = factor(i_group))) + 
        geom_boxplot(outlier.shape = NA) + 
        scale_fill_viridis(option = "viridis",
                           begin = 0.3, 
                           discrete = TRUE,
                           guide = FALSE) + 
        coord_cartesian(ylim = ylim_no_outliers * 1.05) + 
        labs(x = "SES group", y = y_lab)
    }
    
    
    
    
    param_mean_point_plot <- function(data, y, y_lab = NULL) {
      data %>% 
        mutate(i_group = recode_i_group(i_group)) %>% 
        group_by(i_group) %>% 
        filter(!is.na({{y}})) %>% 
        summarise(mean_se({{y}})) %>%
        ggplot(aes(x = factor(i_group), y = y, colour = factor(i_group))) +
        geom_point(size = 2) +
        geom_errorbar(aes(ymin = ymin, ymax = ymax), width = 0.15) +
        labs(x = "SES group", y = y_lab) +
        scale_colour_viridis(option = "viridis",
                             begin = 0.3,
                             discrete = TRUE, guide = FALSE) + 
        theme_custom()
    }
    
    
    mean_sd <- function(data, x) {
      by_group <- data %>% 
        group_by(i_group) %>% 
        summarise(mean = mean_na({{x}}),
                  sd = sd({{x}}, na.rm = TRUE)) %>% 
        mutate(i_group = as.character(i_group))
      
      overall <- data %>% summarise(mean = mean_na({{x}}),
                                    sd = sd({{x}}, na.rm = TRUE)) %>% 
        mutate(i_group = "overall")
      
      bind_rows(by_group, overall)
    
    }
    



    
    
    
    # Group size
    param_bar_plot(params_graph_df$group_props, y_lab = "Group Share") + scale_y_continuous(labels = scales::percent_format(accuracy = 1), minor_breaks = NULL)
    
    

    
    # Test choice delay data
    tc_box <- param_box_plot(data_save$test_delay_data, y = test_choice_delay, y_lab = "Test consulation delay")
    tc_mean <- param_mean_point_plot(data_save$test_delay_data, y = test_choice_delay, y_lab = "Test consulation delay (mean)")
    data_save$test_delay_data %>% mean_sd(test_choice_delay)

    
    
    # Test results delay
    tr_box <- param_box_plot(data_save$test_delay_data, y = test_results_delay, y_lab = "Test results delay")
    tr_mean <- param_mean_point_plot(data_save$test_delay_data, y = test_results_delay, y_lab = "Test results delay (mean)")
    data_save$test_delay_data %>% mean_sd(test_results_delay)

    
    dist_ylim <- c(0, 14)
    mean_ylim <- c(4, 6.5)
    
    ggpubr::ggarrange(
      tc_box + rremove("xlab") + labs(title = "Consultation delay", y = "Days") + coord_cartesian(ylim = dist_ylim), 
      tr_box + rremove("xlab") + labs(title = "Results delay", y = element_blank()) + coord_cartesian(ylim = dist_ylim), 
      tc_mean + labs(y = "Mean Days (±1.96 SE)")+ coord_cartesian(ylim = mean_ylim), 
      tr_mean + rremove("ylab")+ coord_cartesian(ylim = mean_ylim),
      ncol = 2, nrow = 2, align = "hv"
    )
    
    
    
    ggsave("figures/test_delay_dist.pdf", device = cairo_pdf, width = 7, height = 6, scale = 0.7)
    
    # CT delay
    params_graph_df$ct_delay_data %>% 
      ggplot(aes(x = factor(i_group), y = ct_delay, fill = factor(i_group))) + 
      geom_boxplot() + 
      scale_fill_viridis(option = "viridis",
                         begin = 0.3, 
                         discrete = TRUE)
    
    
    # Probabilit of isolating
    tiled_grid <- function(data, colour_lab) {
      data %>% 
        ggplot(aes(x = factor(i_group), y = symptom_severity, fill = prob)) + 
        geom_tile(colour = "white") + 
        geom_text(aes(label = prob, colour = prob > 0.8), size = 4) + 
        scale_colour_manual(values = c("white", "black"), guide = FALSE) + 
        scale_fill_viridis(option = "plasma", limits = c(0, 1)) + 
        labs(y = "Symptom Severity Index", x = "SES Group", fill = colour_lab)
      
    }
    
    param_bar_plot(params_graph_df$probs_self_test_df %>% filter(symptom_severity == 2) %>% .$prob, y_lab = "P(Self Test)")
    param_bar_plot(params_graph_df$probs_isolate_symptoms_df %>% filter(symptom_severity == 2) %>% .$prob, y_lab = "P(Isolate at symptoms)")
    param_bar_plot(params_graph_df$probs_isolate_ct_df %>% filter(symptom_severity == 2) %>% .$prob, y_lab = "P(Isolate after CT)")
    
    # tiled_grid(params_graph_df$probs_self_test_df, "P(Self Test)")
    # tiled_grid(params_graph_df$probs_isolate_symptoms_df, "P(Isolate at Symptoms)")
    # tiled_grid(params_graph_df$probs_isolate_ct_df, "P(Isolate after CT)")
    
    
    
    
    # Contact if isolated
    # param_bar_plot(1 - params_graph_df$p_contact_if_isolated_home, y_lab = "P(Infection at home doesn't occur due to isolation)")
    
    # P contact traced
    param_bar_plot(params_graph_df$p_contact_traced, y_lab = "P(Contact traced)")
    
    # Household quarantine
    params_graph_df$p_hh_quarantine %>% param_bar_plot(y_lab = "P(HH Quarantine)")
    
    # CONtact matrix, contact means  - STILL TO DO 
    
    
    # (params_graph_df$k_matrix / (params_graph_df$pops %*%  t(params_graph_df$pops))) %>% 
    #   as.vector() %>%
    #   {mutate(crossing(to = 1:4, from = 1:4), k_val = .)} %>% 
    #   ggplot(aes(x = factor(from), y = fct_rev(factor(to)), fill = k_val)) + 
    #   geom_tile(colour = "white", size = 0) + 
    #   geom_text(
    #     aes(
    #       label = round(k_val, 5),
    #       colour = k_val > max(k_val) / 2
    #     ), 
    #     size = 4
    #   ) + 
    #   scale_x_discrete(position = "top") + 
    #   scale_colour_manual(values = c("white", "black"), guide = FALSE) + 
    #   scale_fill_viridis(option = "plasma", guide = FALSE) +
    #   labs(x = "Contacts From", y = "Contacts To", fill = "")
    
    
    load("data/processed/hh_data_cleaned.RData", verbose = TRUE)
    
    # HH SIZE
    hh_size_dedup <- hh_data_cleaned %>% 
      dups_drop(hh_id) %>% 
      mutate(i_group = as.integer(stratum), 
             hh_id = as.numeric(hh_id),
             i_group = case_when(i_group %in% c(1, 2) ~ 1L,
                                 i_group %in% c(3, 4) ~ i_group - 1L,
                                 i_group %in% c(5, 6) ~ 4L))
    
    
    # HISTOGRAM
    hh_size_dedup %>% 
      mutate(i_group = recode_i_group(i_group)) %>% 
      ggplot(aes(x = hh_size, fill = factor(i_group)), colour = "darkgrey") + 
      geom_histogram(aes(y=..density..), binwidth = 1, boundary = -0.499, colour = "white") + 
      facet_wrap(~ factor(i_group)) + 
      coord_cartesian(xlim = c(0, 8)) + 
      scale_fill_viridis(option = "viridis",
                         begin = 0.3, 
                         discrete = TRUE, 
                         guide = FALSE) + 
      scale_x_continuous(breaks = 1:8, minor_breaks = NULL) + 
      labs(x = "Household Size", y = "Density") + 
      theme_custom()
    
    ggsave("figures/hh_size_hist.pdf", device = cairo_pdf,  width = 6, height = 4)
    
    hh_size_dedup %>% 
      mutate(i_group = recode_i_group(i_group)) %>% 
      group_by(i_group) %>% 
      summarise(mean_se(hh_size)) %>%
      ggplot(aes(x = factor(i_group), y = y, colour = factor(i_group))) +
      geom_point(size = 3) +
      geom_errorbar(aes(ymin = ymin, ymax = ymax), width = 0) +
      labs(x = "SES group", y = "Mean HH Size") +
      scale_colour_viridis(option = "viridis",
                           begin = 0.3,
                           discrete = TRUE, guide = FALSE) + 
      theme_custom(panel.grid.major.x = element_blank())
    
    
    
    ggsave("figures/hh_size_means.pdf", device = cairo_pdf,  width = 4, height = 3)
    
    
    hh_size_dedup %>% 
      mutate(i_group = recode_i_group(i_group)) %>% 
      group_by(i_group) %>% 
      summarise(mean_se(n_rooms)) %>%
      ggplot(aes(x = factor(i_group), y = y, colour = factor(i_group))) +
      geom_point(size = 3) +
      geom_errorbar(aes(ymin = ymin, ymax = ymax), width = 0) +
      labs(x = "SES group", y = "Mean Number of Rooms") +
      scale_colour_viridis(option = "viridis",
                           begin = 0.3,
                           discrete = TRUE, guide = FALSE) + 
      theme_custom(panel.grid.major.x = element_blank())
    
    
    ggsave("figures/hh_n_rooms.pdf", device = cairo_pdf,  width = 4, height = 3)
    
    
    # data_save$hh_data_bogota %>% mean_sd(hh_size)
    # hh_size_dedup %>% mean_sd(hh_size)
    
    
    
    
    
    # geom_boxplot(outlier.shape = NA) + 
    #+ 
    #coord_cartesian(ylim = ylim_no_outliers * 1.05)
    
    # param_box_plot(y = hh_size)
    
    
    # SAR
    params_graph_df$sar_home %>% param_bar_plot(y_lab = "SAR (Within home)")
    params_graph_df$sar_out %>% param_bar_plot(y_lab = "SAR (Outside home)")
    
    
    
    
    
    

# Contact matrix outputs --------------------------------------------------
    
    
    # BROKEN - due to outbreak_ not outputting the params any more
    # params_graph_df$contact_matrix_ind %>% rowSums %>% View()
    # View(params_graph_df$contact_matrix_ind)
    # write.table(x, "clipboard", sep="\t", row.names=FALSE)
    
    
    # View(params_graph_df$k_matrix * data_save$)
    
    # params_graph_df$contact_weights
    
    # matrix(params_graph_df$contact_weights$contact_weight, byrow = TRUE, nrow = 4) %>% view()
    
    
# Timing distribution graphs ----------------------------------------------
    
    

    
    
    
    
    # DISTRIBUTION OF OUT-OF-HOME CONTACTS
    n_plot <- 100000
    
    
    out_of_home_contacts_plot_data <- tibble(
      i_group = rep(1:4, each = n_plot),
      mu = rep(data_save$contacts_outside_home, each = n_plot),
      contacts = rnbinom(n = n_plot * 4, size = data_save$contact_dispersion, mu = mu)
    ) %>% 
      group_by(i_group, contacts, mu) %>% 
      summarise(n = n()) %>% 
      group_by(i_group, mu) %>% 
      mutate(p = n / sum(n)) %>% 
      mutate(i_group = recode_i_group(i_group)) %>% 
      mutate(label = factor(str_glue("{as.character(i_group)}, mean = {round(mu, 1)}"))) %>% 
      print
      
    
    
    ggplot(out_of_home_contacts_plot_data, aes(x = contacts, y = p, colour = label)) + 
      # geom_density(alpha = 0.5) + 
      geom_line() + 
      geom_point(size = 1) + 
      # scale_x_continuous(trans = "log10") + 
      scale_y_continuous(labels = scales::percent, minor_breaks = NULL) + 
      labs(y = "Density", x = "Contacts outside of home", colour = "SES Group") + 
      coord_cartesian(xlim = c(0, 40)) +
      theme_classic() + 
      theme(legend.position = c(0.8, 0.8)) + 
      annotate("text", x = 35, y = c(0.20, 0.18), hjust = 1,
               label = c(
                 "Negative Binomial",
                 "k = 0.58 (Bi et al 2020)"
               )
      ) + 
      scale_colour_viridis(option = "viridis",
                           begin = 0.3, 
                           discrete = TRUE)
    
    
    save_plot("params_contacts_outside_home")
                                
    
    
    # TEST SENSITIVITY
    data_save$test_sensitivity_data %>% 
      ggplot(aes(x = days_exposure, y = sensitivity)) + 
      geom_line(colour = "indianred") +  
      geom_point(colour = "indianred") + 
      theme_classic() + 
      coord_cartesian(ylim = c(0, 1)) +
      scale_y_continuous(expand = c(0, 0), labels = scales::percent) +
      labs(x = "Days after exposure", y = "Test Sensitivity")
    
    save_plot("params_sensitivity")
    
    
    
    
    # params_timing <- data_save$params_timing
    # params_symptom_timing <- data
    
    
    # SYMPTOMS DISTRIBUTION
    symptom_timing_dist <- tibble(
      x = seq(0, 50, 0.1),
      y = dlnorm(x = x, 
                 meanlog = data_save$params_symptom_timing$meanlog,
                 sdlog = data_save$params_symptom_timing$sdlog)
    )
    
    max_x <- max(symptom_timing_dist$x[symptom_timing_dist$y > 0.001])
    
    ggplot(symptom_timing_dist, aes(x = x, y = y)) + 
      geom_area(colour = "indianred", fill = "indianred", alpha = 0.5) + 
      scale_y_continuous(labels = scales::percent_format(accuracy = 1), minor_breaks = NULL) + 
      labs(y = "Density", x = "Incubation Period (days)") + 
      annotate("text", x = 15, y = c(0.15, 0.13, 0.11), hjust = 1,
               label = c(
                 "Lognormal",
                 expression(paste(mu, " = 1.63, ", sigma, " = 0.50")),
                 "McAloon et al (2020)"
               )
      ) + 
      coord_cartesian(xlim = c(0, max_x)) + 
      theme_classic()
    
    
    save_plot("params_incubation_period")
    
    # SERIAL INTERVAL
    serial_interval_dist <- tibble(
      x = seq(-20, 50, 0.1),
      y = dgamma(x = x - data_save$params_serial$shift, 
                 shape = data_save$params_serial$shape,
                 rate = data_save$params_serial$rate) 
    )
    
    # pgamma(
    #   q = 0 - params_timing$params_serial$shift, 
    #   shape = params_timing$params_serial$shape,
    #   rate = params_timing$params_serial$rate
    # )
    
    min_x <- min(serial_interval_dist$x[serial_interval_dist$y > 0.001])
    max_x <- max(serial_interval_dist$x[serial_interval_dist$y > 0.001])
    

    ggplot(serial_interval_dist, aes(x = x, y = y)) + 
      geom_vline(xintercept = 0, linetype = "dotted") + 
      geom_area(colour = "dodgerblue", fill = "dodgerblue", alpha = 0.5) + 
      scale_y_continuous(labels = scales::percent_format(accuracy = 1), minor_breaks = NULL) + 
      labs(y = "Density", x = "Serial Interval (days)") + 
      annotate("text", x = 20, y = c(0.08, 0.07, 0.06), hjust = 1,
               label = c(
                 "Gamma",
                 expression(paste(alpha, " = 8.12, ", beta, " = 0.64, ", Delta, "x = -7.5")),
                 "He et al (2020)"
               )
      ) +
      coord_cartesian(xlim = c(min_x, max_x)) +
      theme_classic()
    
    save_plot("params_serial_interval")
    
    
    # IMPORT THE DATA FROM IMPORT_DATA - using the draw_timings function 
    # and the gen_data_clean stuff
    
    # Secondary timing
    ggplot(gen_data_clean, aes(x = secondary_timing)) + 
      geom_density(colour = "seagreen4", fill = "seagreen4", alpha = 0.5) + 
      scale_y_continuous(labels = scales::percent_format(accuracy = 1), minor_breaks = NULL) + 
      labs(y = "Density", x = "Generation Interval (days)") + 
      coord_cartesian(xlim = c(0, 20)) +
      theme_classic()
    
    
    save_plot("params_generation_interval")
    
    # Secondary timing - symptoms [infectiousness profile]
    prop_presymptomatic <- gen_data_clean %$% round(mean(secondary_timing < primary_symptoms) * 100)
    
    ggplot(gen_data_clean, aes(x = secondary_timing - primary_symptoms)) + 
      geom_vline(xintercept = 0, linetype = "dotted") + 
      geom_density(colour = "darkorange2", fill = "darkorange2", alpha = 0.5) + 
      scale_y_continuous(labels = scales::percent_format(accuracy = 1), minor_breaks = NULL) + 
      labs(y = "Density", x = "Day of secondary infection - day of symptom onset") + 
      annotate("text", x = 20, y = c(0.08), hjust = 1,
               label = c(
                 paste0("Prop presymptomatic = ", prop_presymptomatic, "%")
               )
      ) +
      coord_cartesian(xlim = c(-10, 20)) +
      theme_classic()
    
    save_plot("params_infectiousness_profile")
    
    
    

# Timing distribution calcs -----------------------------------------------

    # For generation interval
    gamma_shape <- data_save$params_timing[[8]]
    secondary_mean <- data_save$serial_mean
    pgamma(q = 1, shape = gamma_shape, scale = secondary_mean/gamma_shape) # prob of being truncated is 0.36%
    
    tibble(
      x = rgamma(n = 10000,
                 shape = gamma_shape,
                 scale = secondary_mean/gamma_shape)
    ) %>% 
      filter(x >= 1) %>% 
      {print(mean(.$x)); .} %>% 
      ggplot(aes(x = x)) + 
      geom_density()

    
    
    
    