library(dplyr)

# Load data and contact matrices
map_data <- readRDS("ProcessedData/map_data.rds")
C1 <- readRDS("ProcessedData/contact_matrix1.rds")
C2 <- readRDS("ProcessedData/contact_matrix2.rds")
C3 <- readRDS("ProcessedData/contact_matrix3.rds")

# Extract needed columns
district_name <- map_data$NAME
psi <- map_data$susceptible_prop         # 1 - efficacy * MMR
pmo <- map_data$outbreak_prob            # local outbreak probability
P_internal <- map_data$internal_infection_prob  # internal infection prob

# Transmission parameter
q <- 0.9

# Transmission matrix function
compute_transmission_matrix <- function(Cij, Pj, psi, pmo_i, q, names) {
  n <- length(Pj)
  mat <- matrix(0, n, n, dimnames = list(names, names))
  for (i in 1:n) {
    for (j in 1:n) {
      base <- q * Pj[j] * psi[i] * pmo_i[i]
      mat[i, j] <- 1 - (1 - base)^Cij[i, j]
    }
  }
  return(mat)
}

# Compute pij for all methods
pij_1 <- compute_transmission_matrix(C1, P_internal, psi, pmo, q, district_name)
pij_2 <- compute_transmission_matrix(C2, P_internal, psi, pmo, q, district_name)
pij_3 <- compute_transmission_matrix(C3, P_internal, psi, pmo, q, district_name)

# Save to disk
saveRDS(pij_1, "ProcessedData/transmission_matrix1.rds")
saveRDS(pij_2, "ProcessedData/transmission_matrix1.rds")
saveRDS(pij_3, "ProcessedData/transmission_matrix1.rds")
