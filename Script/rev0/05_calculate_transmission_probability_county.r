library(dplyr)

# Load data and contact matrices
map_data <- readRDS("ProcessedData/map_data_county.rds")
C1 <- readRDS("ProcessedData/contact_matrix1_county.rds")
C2 <- readRDS("ProcessedData/contact_matrix2_county.rds")
C3 <- readRDS("ProcessedData/contact_matrix3_county.rds")
C4 <- readRDS("ProcessedData/contact_matrix4_county.rds")
C5 <- readRDS("ProcessedData/contact_matrix5_county.rds")
C6 <- readRDS("ProcessedData/contact_matrix6_county.rds")
C7 <- readRDS("ProcessedData/contact_matrix7_county.rds")


compute_transmission_matrix <- function(Cij, map_data, q = 0.9) {
  county_names <- map_data$County
  psi <- map_data$susceptible_prop
  pmo <- map_data$outbreak_prob
  Pj <- map_data$internal_infection_prob
  
  n <- length(county_names)
  transmission_mat <- matrix(0, n, n, dimnames = list(county_names, county_names))
  
  for (i in seq_len(n)) {
    for (j in seq_len(n)) {
      name_i <- county_names[i]
      name_j <- county_names[j]
      
      base <- q * Pj[which(county_names == name_j)] *
        psi[which(county_names == name_i)] *
        pmo[which(county_names == name_i)]
      
      cij_val <- Cij[name_i, name_j]
      transmission_mat[name_i, name_j] <- 1 - (1 - base)^cij_val
    }
  }
  
  return(transmission_mat)
}


# Compute pij for all methods
pij_1 <- compute_transmission_matrix(C1, map_data)
pij_2 <- compute_transmission_matrix(C2, map_data)
pij_3 <- compute_transmission_matrix(C3, map_data)
pij_4 <- compute_transmission_matrix(C4, map_data)
pij_5 <- compute_transmission_matrix(C5, map_data)
pij_6 <- compute_transmission_matrix(C6, map_data)
pij_7 <- compute_transmission_matrix(C7, map_data)

# 
saveRDS(pij_1, "ProcessedData/transmission_matrix_county1.rds")
saveRDS(pij_2, "ProcessedData/transmission_matrix_county2.rds")
saveRDS(pij_3, "ProcessedData/transmission_matrix_county3.rds")
saveRDS(pij_4, "ProcessedData/transmission_matrix_county4.rds")
saveRDS(pij_5, "ProcessedData/transmission_matrix_county5.rds")
saveRDS(pij_6, "ProcessedData/transmission_matrix_county6.rds")
saveRDS(pij_7, "ProcessedData/transmission_matrix_county7.rds")



######### TEST SECTION
# pij_7["GAINES", "NEWTON"]
# C7["GAINES", "NEWTON"]
# 
# C1["GAINES", "NEWTON"]
# format(C7["GAINES", "NEWTON"], scientific = TRUE, digits = 20)
# C3["GAINES", "NEWTON"]

