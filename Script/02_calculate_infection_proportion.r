library(sf)
library(ggplot2)
library(reshape2)
library(lubridate)
library(dplyr)
library(tidyr)
library(readr)
library(tigris)

#------------------ Load datasets -------------------
mmr_county <- readRDS("ProcessedData/merged_map_county_df.rds")


# Set vaccine efficacy and basic reproduction number
efficacy <- 0.97
ratio <- 0.5
R0 <- 18

# Infection Proportion Functions ----
# Interface (Lin. )
find_internal_infection_lin <- function(v, R0_value = R0, tol = 1e-4, efficacy_rate=efficacy) {
  if (is.na(v)) return(NA_real_)
  a <- 1 - efficacy_rate  * v
  X_vals <- seq(0, 1, length.out = 10000)
  f_X <- X_vals - a * (1 - exp(-R0_value * X_vals))
  close_to_zero <- abs(f_X) < tol
  return(max(X_vals[close_to_zero]))
}


# Calculate Probability of infection for County ----

infection_county <- mmr_county %>%
  mutate(
    pop_above18 = pop_total - pop_under18,
  ) %>%
  select(
    County = county,
    MMR = mmr,
    pop_above18,
    pop_under18
  ) %>%
  mutate(
    MMR1 = ifelse(MMR < .90, .90, MMR),           # Strategy 1
    MMR2 = ifelse(MMR < .92, .92, MMR),           # Strategy 2
    MMR3 = ifelse(MMR < .95, .95, MMR),           # Strategy 3
    
    under18_infection_prob  = sapply(MMR, find_internal_infection_lin),
    under18_infection_prob1 = sapply(MMR1, find_internal_infection_lin),
    under18_infection_prob2 = sapply(MMR2, find_internal_infection_lin),
    under18_infection_prob3 = sapply(MMR3, find_internal_infection_lin),
    
    adult_infection_prob  = (pop_under18 / (ratio * (pop_above18))) * under18_infection_prob,
    adult_infection_prob1 = (pop_under18 / (ratio * (pop_above18))) * under18_infection_prob1,
    adult_infection_prob2 = (pop_under18 / (ratio * (pop_above18))) * under18_infection_prob2,
    adult_infection_prob3 = (pop_under18 / (ratio * (pop_above18))) * under18_infection_prob3,
    
    # P^S = (1−eV)
    susceptible_prop  = 1 - (efficacy * MMR), 
    susceptible_prop1 = 1 - (efficacy * MMR1),
    susceptible_prop2 = 1 - (efficacy * MMR2),
    susceptible_prop3 = 1 - (efficacy * MMR3),
    
    #P^(LO) = 1 − [1 / ((1−eV)R0)]
    outbreak_prob  = pmax(0, 1 - (1 / ((1 - efficacy * MMR) * R0))), 
    outbreak_prob1 = pmax(0, 1 - (1 / ((1 - efficacy * MMR1) * R0))),
    outbreak_prob2 = pmax(0, 1 - (1 / ((1 - efficacy * MMR2) * R0))),
    outbreak_prob3 = pmax(0, 1 - (1 / ((1 - efficacy * MMR3) * R0)))

  )


saveRDS(infection_county, "ProcessedData/map_county_infection_proportion_ratio0.5_efficacy0.97.rds")

# Visualization -----
# p_county_under18 <- ggplot(infection_county) +
#   geom_sf(
#     aes(fill = under18_infection_prob),
#     color = "gray30",
#     size = 0.1
#   ) +
#   scale_fill_gradientn(
#     colors = c("#1a9850", "#fee08b", "#d73027"),  
#     limits = c(0, .15),
#     na.value = "lightgray",
#     labels = scales::percent_format(accuracy = 1)
#   ) +
#   theme_minimal() +
#   labs(
#     title = "Texas counties infection probability (under 18 years old)",
#     fill = "Infection Probability"
#   )
# 
# #save to file
# ggsave("Figures/infection_probability_county_map_lin_under18.png",
#        plot = p_county_under18,
#        width = 10,
#        height = 8,
#        dpi = 300
# )
# 
# 
# p_county_adult <- ggplot(infection_county) +
#   geom_sf(
#     aes(fill = adult_infection_prob),
#     color = "gray30",
#     size = 0.1
#   ) +
#   scale_fill_gradientn(
#     colors = c("#1a9850", "#fee08b", "#d73027"),  
#     limits = c(0, .15),
#     na.value = "lightgray",
#     labels = scales::percent_format(accuracy = 1)
#   ) +
#   theme_minimal() +
#   labs(
#     title = "Texas counties infection probability (above 18 years old)",
#     fill = "Infection Probability"
#   )
# 
# #save to file
# ggsave("Figures/infection_probability_county_map_lin_adult.png",
#        plot = p_county_adult,
#        width = 10,
#        height = 8,
#        dpi = 300
# )
# 









