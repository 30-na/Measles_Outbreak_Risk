library(tidygeocoder)
library(sf)
library(reshape2)
library(ggplot2)

map_data <- readRDS("ProcessedData/map.rds")
dist_matrix <- readRDS("ProcessedData/distance_matrix_haversine.rds")
pop <- map_data$total_population


# Method 1: Model A
contact_method1 <- function(pop, dist_matrix) {
  A <- 75.94
  B <- 278e-9
  C <- 1.85e4
  D <- 3.43e8
  alpha <- 1.80
  gamma <- 1.16
  
  n <- length(pop)
  mat <- matrix(NA, n, n)
  for (i in 1:n) {
    for (j in 1:n) {
      if (dist_matrix[i, j] == 0) next
      mi <- pop[i]
      mj <- pop[j]
      dij <- dist_matrix[i, j]
      Tij <- A * ((B * (mi * mj + C * mj + D)) / (dij^alpha) + 1)^gamma
      mat[i, j] <- Tij
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
  
  n <- length(pop)
  mat <- matrix(NA, n, n)
  for (i in 1:n) {
    for (j in 1:n) {
      if (dist_matrix[i, j] == 0) next
      mi <- pop[i]
      mj <- pop[j]
      dij <- dist_matrix[i, j]
      term <- 1 + (B * ((mi + C) * (mj + D))^beta / dij)
      Tij <- exp(A * (term)^xi)
      mat[i, j] <- Tij
    }
  }
  return(mat)
}

# Method 3: Gravity Model
contact_method3 <- function(pop, dist_matrix) {
  dist_matrix[dist_matrix == 0] <- NA  # avoid division by zero
  outer(pop, pop, FUN = function(m1, m2) m1 * m2) / (dist_matrix^2)
}


# Compute and save
contact_Haversine01 <- contact_method1(pop, dist_matrix)
saveRDS(contact_Haversine01, "ProcessedData/contact_matrix1.rds")

contact_Haversine02 <- contact_method2(pop, dist_matrix)
saveRDS(contact_Haversine02, "ProcessedData/contact_matrix2.rds")

contact_Haversine03 <- contact_method3(pop, dist_matrix)
saveRDS(contact_Haversine03, "ProcessedData/contact_matrix3.rds")


# Optional plots 
# Convert to long format
#contact_long <- melt(contact_Haversine03, varnames = c("District_i", "District_j"), value.name = "Contact")

# Heatmap
#g_h03 <- ggplot(contact_long, aes(x = District_j, y = District_i, fill = Contact)) +
#  geom_tile() +
#  scale_fill_viridis_c(trans = "log", option = "magma") +
#  theme_minimal() +
#  theme(axis.text.x = element_blank(), axis.text.y = element_blank()) +
#  labs(title = "Gravity-Based Contact Matrix (Method 3) using Haversine distance", fill = "Cij")


#ggsave("Figures/contact_matrix_haversine03.png", 
#       plot = g_h03, width = 6, height = 5, dpi = 300)




