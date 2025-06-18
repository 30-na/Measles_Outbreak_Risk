library(ggplot2)
library(dplyr)
library(tidyr)
library(sf)

# Load county-level map data
map_data <- readRDS("ProcessedData/map_county.rds")

# Set vaccine efficacy and basic reproduction number
efficacy <- 0.97
R0 <- 18

# Compute susceptible proportion and outbreak probability
map_data <- map_data %>%
  mutate(
    susceptible_prop = 1 - (efficacy * MMR),
    outbreak_prob = pmax(0, 1 - (1 / ((1 - efficacy * MMR) * R0)))
  )

# Define function for internal infection probability
find_internal_infection_numeric <- function(Vj, efficacy_rate = efficacy, R0 = 18, tol = 1e-4) {
  if (is.na(Vj)) return(NA_real_)
  a <- 1 - efficacy * Vj
  X_vals <- seq(0, 1, length.out = 10000)
  f_X <- X_vals - a * (1 - exp(-a * R0 * X_vals))
  close_to_zero <- abs(f_X) < tol
  return(max(X_vals[close_to_zero]))
}

# Apply function to MMR values
map_data$internal_infection_prob <- sapply(map_data$MMR, find_internal_infection_numeric)
map_data$internal_infection_prob_log <- log(map_data$internal_infection_prob)
map_data$internal_infection_prob_log10 <- log10(map_data$internal_infection_prob)

# Compute average outbreak size
map_data <- map_data %>%
  mutate(
    pop_time_infection_total = internal_infection_prob * total,
    pop_time_infection_total_log10 = pmax(0, log10(pop_time_infection_total))
  )

# Save updated data
saveRDS(map_data, "ProcessedData/map_data_county.rds")
