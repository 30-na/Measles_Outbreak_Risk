library(ggplot2)
library(dplyr)
library(tidyr)
library(sf)


map_data <- readRDS("ProcessedData/map.rds")

########## Susceptible Proportion and  Local Outbreak Probabilit 

efficacy <- 0.97
R0 <- 18


map_data <- map_data %>%
  mutate(
    susceptible_prop = 1 - (efficacy * MMR), # Susceptible Proportion
    outbreak_prob = pmax(0, 1 - (1 / ((1 - efficacy * MMR) * R0))) #  Outbreak Probability
  )


######### Internal Infection Probability (average outbreak rate)

find_internal_infection_numeric <- function(Vj, efficacy = .97, R0 = 18, tol = 1e-4) {
  if (is.na(Vj)) return(NA_real_)
  
  a <- 1 - efficacy * Vj
  X_vals <- seq(0, 1, length.out = 10000)
  f_X <- X_vals - a * (1 - exp(-a * R0 * X_vals))
  
  close_to_zero <- abs(f_X) < tol
  
  return(max(X_vals[close_to_zero]))
}


map_data$internal_infection_prob <- sapply(map_data$MMR, find_internal_infection_numeric)
map_data$internal_infection_prob_log <- log(map_data$internal_infection_prob)
map_data$internal_infection_prob_log10 <- log10(map_data$internal_infection_prob)


######## Average outbreak size
map_data <- map_data %>%
  mutate(
    #pop_time_infection_under18 = internal_infection_prob * Under18_Population,
    pop_time_infection_total = internal_infection_prob * total_population,
    pop_time_infection_total_log10 = pmax(0, log10(pop_time_infection_total))
  )


# update map data
saveRDS(map_data, "ProcessedData/map_data.rds")








# plot_district_map(map_data, "internal_infection_prob_log", "Internal Infection Probability (Log)", file_name="internal_infection_prob_log")
# plot_district_map(map_data, "internal_infection_prob_log10", "Internal Infection Probability (Log10)", file_name="internal_infection_prob_log10")
# 
# plot_district_map(map_data, "internal_infection_prob", "Internal Infection Probability", file_name="internal_infection_prob")
# 
# plot_district_map(map_data, "pop_time_infection_under18", "Product of Internal Infection Probability and Under 18 Popualation (efficacy = .93)", file_name="infection_times_under18_popualtion_93")
# 
# plot_district_map(map_data, "pop_time_infection_total", "Product of Internal Infection Probability and Total Popualation (efficacy = .93)", file_name="infection_times_total_popualtion_93")
# plot_district_map(map_data, "pop_time_infection_total_log", "Product of Internal Infection Probability and Total Popualation (efficacy = .97), Log", file_name="infection_times_total_popualtion_97_log")
# plot_district_map(map_data, "pop_time_infection_total_log10", "Product of Internal Infection Probability and Total Popualation (efficacy = .97), Log10", file_name="infection_times_total_popualtion_97_log10")
