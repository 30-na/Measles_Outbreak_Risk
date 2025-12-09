library(sf)
library(ggplot2)
library(reshape2)
library(lubridate)
library(dplyr)
library(tidyr)
library(readr)
library(tigris)

#Load datasets ----
infection_county <- readRDS("ProcessedData/map_county_infection_proportion_ratio01_efficacy0.97.rds")

# Population Flow County
texas_flows <- read_csv("ProcessedData/texas_county_flows_2019.csv", 
                        col_types = cols(.default = "c"))



# Functions ----
contact_method7 <- function(flows_avg) {
  county_names <- sort(unique(c(flows_avg$county_o, flows_avg$county_d)))
  mat <- matrix(0, 
                length(county_names), 
                length(county_names),
                dimnames = list(county_names, county_names))
  
  for (k in seq_len(nrow(flows_avg))) {
    origin <- flows_avg$county_o[k]
    dest <- flows_avg$county_d[k]
    flow <- flows_avg$flow[k]
    
    if (origin %in% county_names && dest %in% county_names) {
      if (origin == dest) {
        mat[dest, origin] <- NA  # self-flows as NA
      } else {
        mat[dest, origin] <- flow
      }
    }
  }
  
  return(mat)
}

# Texas Population Flow County ----
# Method7
flows_avg <- texas_flows %>%
  mutate(
    pop_flows = as.numeric(pop_flows),
    week_start = mdy(str_sub(date_range, 1, 8))
  ) %>%
  # filter(month(week_start) >= 1 & month(week_start) <= 6) %>%
  group_by(
    geoid_o,
    geoid_d, 
    county_o, 
    county_d
  ) %>%
  summarize(
    flow = sum(pop_flows, na.rm = TRUE),
    .groups = "drop"
  )

C7 <- contact_method7(flows_avg)


#------------------ Transmission/ Outbreak Probability COUNTY ------------------

compute_transmission_matrix_county <- function(Cij, map_data, strategy = 0, q = .9) {
  suffix <- ifelse(strategy == 0, "", as.character(strategy))
  
  county_names <- map_data$County
  psi <- map_data[[paste0("susceptible_prop", suffix)]]
  pmo <- map_data[[paste0("outbreak_prob", suffix)]]
  Pj <- map_data[[paste0("adult_infection_prob", suffix)]]
  
  n <- length(county_names)
  transmission_mat <- matrix(NA, n, n, dimnames = list(county_names, county_names))
  
  for (i in seq_len(n)) {
    for (j in seq_len(n)) {
      name_i <- county_names[i]
      name_j <- county_names[j]
      
      base <- q * Pj[which(county_names == name_j)] *
        psi[which(county_names == name_i)] *
        pmo[which(county_names == name_i)]
      
      cij_val <- Cij[name_i, name_j]
      if (is.na(cij_val)) {
        transmission_mat[name_i, name_j] <- NA
      } else {
        transmission_mat[name_i, name_j] <- 1 - (1 - base)^cij_val
      }
      
    }
  }
  
  return(transmission_mat)
}


# Method 7 (Mobility Flows)
county_pij_M7_S0 <- compute_transmission_matrix_county(C7, infection_county)
county_pij_M7_S1 <- compute_transmission_matrix_county(C7, infection_county, strategy = 1)
county_pij_M7_S2 <- compute_transmission_matrix_county(C7, infection_county, strategy = 2)
county_pij_M7_S3 <- compute_transmission_matrix_county(C7, infection_county, strategy = 3)


saveRDS(county_pij_M7_S0, "ProcessedData/county_pij_M7_S0.rds")
saveRDS(county_pij_M7_S1, "ProcessedData/county_pij_M7_S1.rds")
saveRDS(county_pij_M7_S2, "ProcessedData/county_pij_M7_S2.rds")
saveRDS(county_pij_M7_S3, "ProcessedData/county_pij_M7_S3.rds")



