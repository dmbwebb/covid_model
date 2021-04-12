
r0_store <- list()

# 4 choose initial value for outside of work factor, and match to R0 slope
# outside_work_factor <- 0.7
# source("code/calc_contact_matrix.R")
# source("code/r0_calibrate.R")
# 
# outside_work_factor <- 0.75
# source("code/calc_contact_matrix.R")
# source("code/r0_calibrate.R") 

outside_work_factor <- 0.8
source("code/calc_contact_matrix.R")
source("code/r0_calibrate.R") 

# outside_work_factor <- 0.85
# source("code/calc_contact_matrix.R")
# source("code/r0_calibrate.R")
# 
# outside_work_factor <- 0.9
# source("code/calc_contact_matrix.R")
# source("code/r0_calibrate.R") 


data_save$r0_store <- r0_store

r0_store %>% print


r0_stores <- tibble(
  outside_work_factor = names(r0_store),
  r0_est = map_dbl(r0_store, "r0")
) %>% 
  mutate(diff = abs(r0_est - data_save$initial_r0))

# r0_stores$r0_est
# data_save$initial_r0

print(r0_stores)

outside_work_factor <- as.numeric(r0_stores$outside_work_factor[r0_stores$diff == min(r0_stores$diff)])
print(str_glue("Chosen outside_work_factor is {outside_work_factor}"))
