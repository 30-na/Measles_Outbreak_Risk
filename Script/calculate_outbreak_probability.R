library(sf)
library(ggplot2)
library(reshape2)
library(lubridate)
library(dplyr)
library(tidyr)
library(readr)

# Load county-level map data
# Load county-level data
map_data <- readRDS("ProcessedData/map_county.rds")
dist_matrix <- readRDS("ProcessedData/distance_matrix_haversine_county.rds")
texas_flows <- read_csv("ProcessedData/texas_county_flows_2019.csv", 
                        col_types = cols(.default = "c"))

# Set vaccine efficacy and basic reproduction number
efficacy <- 0.97
R0 <- 18


map_data <- map_data %>%
  mutate(
    MMR1 = ifelse(MMR < .90, .90, MMR),           # Strategy 1
    MMR2 = ifelse(MMR < .92, .92, MMR),           # Strategy 2
    MMR3 = pmin(MMR + .05, 1),                    # Strategy 3
    MMR4 = ifelse(County == "GAINES", 0.90, MMR),         # Strategy 4
    MMR5 = ifelse(County == "GAINES", 0.92, MMR),         # Strategy 5
    
    susceptible_prop = 1 - (efficacy * MMR),
    susceptible_prop1 = 1 - (efficacy * MMR1),
    susceptible_prop2 = 1 - (efficacy * MMR2),
    susceptible_prop3 = 1 - (efficacy * MMR3),
    susceptible_prop4 = 1 - (efficacy * MMR4),
    susceptible_prop5 = 1 - (efficacy * MMR5),
    
    susceptible_pop_size = susceptible_prop * total,
    susceptible_pop_size1 = susceptible_prop1 * total,
    susceptible_pop_size2 = susceptible_prop2 * total,
    susceptible_pop_size3 = susceptible_prop3 * total,
    susceptible_pop_size4 = susceptible_prop4 * total,
    susceptible_pop_size5 = susceptible_prop5 * total,
    
    outbreak_prob = pmax(0, 1 - (1 / ((1 - efficacy * MMR) * R0))),
    outbreak_prob1 = pmax(0, 1 - (1 / ((1 - efficacy * MMR1) * R0))),
    outbreak_prob2 = pmax(0, 1 - (1 / ((1 - efficacy * MMR2) * R0))),
    outbreak_prob3 = pmax(0, 1 - (1 / ((1 - efficacy * MMR3) * R0))),
    outbreak_prob4 = pmax(0, 1 - (1 / ((1 - efficacy * MMR4) * R0))),
    outbreak_prob5 = pmax(0, 1 - (1 / ((1 - efficacy * MMR5) * R0)))
    
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
map_data$internal_infection_prob1 <- sapply(map_data$MMR1, find_internal_infection_numeric)
map_data$internal_infection_prob2 <- sapply(map_data$MMR2, find_internal_infection_numeric)
map_data$internal_infection_prob3 <- sapply(map_data$MMR3, find_internal_infection_numeric)
map_data$internal_infection_prob4 <- sapply(map_data$MMR4, find_internal_infection_numeric)
map_data$internal_infection_prob5 <- sapply(map_data$MMR5, find_internal_infection_numeric)


# Compute average outbreak size
map_data <- map_data %>%
  mutate(
    avg_outbreak_size = internal_infection_prob * total,
    avg_outbreak_size1 = internal_infection_prob1 * total,
    avg_outbreak_size2 = internal_infection_prob2 * total,
    avg_outbreak_size3 = internal_infection_prob3 * total,
    avg_outbreak_size4 = internal_infection_prob4 * total,
    avg_outbreak_size5 = internal_infection_prob5 * total
  )

saveRDS(map_data, "ProcessedData/map_probability.rds")
############ Contact Matrix ###########

pop <- map_data$total  # county population
names(pop) <- map_data$County




# Method 1: Model A
contact_method1 <- function(pop, dist_matrix) {
  A <- 75.94
  B <- 278e-9
  C <- 1.85e4
  D <- 3.43e8
  alpha <- 1.80
  gamma <- 1.16
  
  county_names <- rownames(dist_matrix)
  mat <- matrix(NA, length(county_names), length(county_names),
                dimnames = list(county_names, county_names))
  
  for (i in county_names) {
    for (j in county_names) {
      dij <- dist_matrix[i, j]
      if (dij == 0) next  # skip self-distances or undefined cases
      mi <- pop[i]
      mj <- pop[j]
      Tij <- A * ((B * (mi * mj + C * mj + D)) / (dij^alpha) + 1)^gamma
      mat[i, j] <- 365 * Tij
    }
  }
  return(mat)
}


# Method 2: Model B
contact_method2 <- function(pop, dist_matrix) {
  A <- 4.10
  B <- 1240e-6
  C <- 61.2
  D <- 1.79e4
  beta <- 0.50
  xi <- 0.30
  
  county_names <- rownames(dist_matrix)
  mat <- matrix(NA, length(county_names), length(county_names),
                dimnames = list(county_names, county_names))
  
  for (i in county_names) {
    for (j in county_names) {
      dij <- dist_matrix[i, j]
      if (dij == 0) next
      mi <- pop[i]
      mj <- pop[j]
      term <- 1 + (B * ((mi + C) * (mj + D))^beta / dij)
      Tij <- exp(A * (term)^xi)
      mat[i, j] <- 365 * Tij
    }
  }
  return(mat)
}



# Method 3: Gravity Model
contact_method3 <- function(pop, dist_matrix) {
  county_names <- rownames(dist_matrix)
  mat <- matrix(NA, length(county_names), length(county_names),
                dimnames = list(county_names, county_names))
  
  for (i in county_names) {
    for (j in county_names) {
      dij <- dist_matrix[i, j]
      if (dij == 0 || is.na(dij)) next  # avoid divide by zero
      mi <- pop[i]
      mj <- pop[j]
      Tij <- (mi * mj) / (dij^2)
      mat[i, j] <- Tij * .33
    }
  }
  return(mat)
}



contact_method7 <- function(flows_avg) {
  county_names <- sort(unique(c(flows_avg$county_o, flows_avg$county_d)))
  mat <- matrix(0, length(county_names), length(county_names),
                dimnames = list(county_names, county_names))
  
  for (k in seq_len(nrow(flows_avg))) {
    origin <- flows_avg$county_o[k]
    dest <- flows_avg$county_d[k]
    flow <- flows_avg$flow[k]
    
    if (origin %in% county_names && dest %in% county_names) {
      if (origin == dest) {
        mat[origin, dest] <- NA  # self-flows as NA
      } else {
        mat[origin, dest] <- flow
      }
    }
  }
  
  return(mat)
}

  
  

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

C1 <- contact_method1(pop, dist_matrix)
C2 <- contact_method2(pop, dist_matrix)
C3 <- contact_method3(pop, dist_matrix)
C7 <- contact_method7(flows_avg)






############ Transmission/ Outbreak Probability ############


compute_transmission_matrix <- function(Cij, map_data, strategy = 0, q = 0.9) {
  suffix <- ifelse(strategy == 0, "", as.character(strategy))
  
  county_names <- map_data$County
  psi <- map_data[[paste0("susceptible_prop", suffix)]]
  pmo <- map_data[[paste0("outbreak_prob", suffix)]]
  Pj <- map_data[[paste0("internal_infection_prob", suffix)]]
  
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

# Method 1 (Model A)
pij_M1_S0 <- compute_transmission_matrix(C1, map_data)
pij_M1_S1 <- compute_transmission_matrix(C1, map_data, strategy = 1)
pij_M1_S2 <- compute_transmission_matrix(C1, map_data, strategy = 2)
pij_M1_S3 <- compute_transmission_matrix(C1, map_data, strategy = 3)
pij_M1_S4 <- compute_transmission_matrix(C1, map_data, strategy = 4)
pij_M1_S5 <- compute_transmission_matrix(C1, map_data, strategy = 5)

# Method 2 (Model B)
pij_M2_S0 <- compute_transmission_matrix(C2, map_data)
pij_M2_S1 <- compute_transmission_matrix(C2, map_data, strategy = 1)
pij_M2_S2 <- compute_transmission_matrix(C2, map_data, strategy = 2)
pij_M2_S3 <- compute_transmission_matrix(C2, map_data, strategy = 3)
pij_M2_S4 <- compute_transmission_matrix(C2, map_data, strategy = 4)
pij_M2_S5 <- compute_transmission_matrix(C2, map_data, strategy = 5)

# Compute pij for all methods and Strategies
# Method 3 (Gravity Model)
pij_M3_S0 <- compute_transmission_matrix(C3, map_data)
pij_M3_S1 <- compute_transmission_matrix(C3, map_data, strategy = 1)
pij_M3_S2 <- compute_transmission_matrix(C3, map_data, strategy = 2)
pij_M3_S3 <- compute_transmission_matrix(C3, map_data, strategy = 3)
pij_M3_S4 <- compute_transmission_matrix(C3, map_data, strategy = 4)
pij_M3_S5 <- compute_transmission_matrix(C3, map_data, strategy = 5)

# Method 7 (Mobility Flows)
pij_M7_S0 <- compute_transmission_matrix(C7, map_data)
pij_M7_S1 <- compute_transmission_matrix(C7, map_data, strategy = 1)
pij_M7_S2 <- compute_transmission_matrix(C7, map_data, strategy = 2)
pij_M7_S3 <- compute_transmission_matrix(C7, map_data, strategy = 3)
pij_M7_S4 <- compute_transmission_matrix(C7, map_data, strategy = 4)
pij_M7_S5 <- compute_transmission_matrix(C7, map_data, strategy = 5)




saveRDS(pij_M1_S0, "ProcessedData/pij_M1_S0.rds")
saveRDS(pij_M1_S1, "ProcessedData/pij_M1_S1.rds")
saveRDS(pij_M1_S2, "ProcessedData/pij_M1_S2.rds")
saveRDS(pij_M1_S3, "ProcessedData/pij_M1_S3.rds")
saveRDS(pij_M1_S4, "ProcessedData/pij_M1_S4.rds")
saveRDS(pij_M1_S5, "ProcessedData/pij_M1_S5.rds")

saveRDS(pij_M2_S0, "ProcessedData/pij_M2_S0.rds")
saveRDS(pij_M2_S1, "ProcessedData/pij_M2_S1.rds")
saveRDS(pij_M2_S2, "ProcessedData/pij_M2_S2.rds")
saveRDS(pij_M2_S3, "ProcessedData/pij_M2_S3.rds")
saveRDS(pij_M2_S4, "ProcessedData/pij_M2_S4.rds")
saveRDS(pij_M2_S5, "ProcessedData/pij_M2_S5.rds")


saveRDS(pij_M3_S0, "ProcessedData/pij_M3_S0.rds")
saveRDS(pij_M3_S1, "ProcessedData/pij_M3_S1.rds")
saveRDS(pij_M3_S2, "ProcessedData/pij_M3_S2.rds")
saveRDS(pij_M3_S3, "ProcessedData/pij_M3_S3.rds")
saveRDS(pij_M3_S4, "ProcessedData/pij_M3_S4.rds")
saveRDS(pij_M3_S5, "ProcessedData/pij_M3_S5.rds")


saveRDS(pij_M7_S0, "ProcessedData/pij_M7_S0.rds")
saveRDS(pij_M7_S1, "ProcessedData/pij_M7_S1.rds")
saveRDS(pij_M7_S2, "ProcessedData/pij_M7_S2.rds")
saveRDS(pij_M7_S3, "ProcessedData/pij_M7_S3.rds")
saveRDS(pij_M7_S4, "ProcessedData/pij_M7_S4.rds")
saveRDS(pij_M7_S5, "ProcessedData/pij_M7_S5.rds")








