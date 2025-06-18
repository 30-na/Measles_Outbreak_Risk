library(sf)
library(geosphere)
library(ggplot2)

# Load processed county-level spatial data
map_data <- readRDS("ProcessedData/map_county.rds")

# Ensure it's an sf object
map_data <- st_as_sf(map_data)

# Compute centroids
centroids <- st_centroid(map_data)

# Extract longitude/latitude
coords <- st_coordinates(centroids)
map_data$longitude <- coords[, 1]
map_data$latitude <- coords[, 2]

# Prepare coordinate input
locations <- map_data[, c("longitude", "latitude")] |> st_drop_geometry()

# Initialize distance matrix
n <- nrow(locations)
dist_matrix <- matrix(NA, n, n, dimnames = list(map_data$County, map_data$County))

# Compute Haversine distances in km
pb <- txtProgressBar(min = 0, max = n, style = 3)
for (i in 1:n) {
  for (j in 1:n) {
    dist_matrix[i, j] <- distHaversine(locations[i, ], locations[j, ]) / 1000
  }
  setTxtProgressBar(pb, i)
}
close(pb)

# Save the matrix
saveRDS(dist_matrix, "ProcessedData/distance_matrix_haversine_county.rds")



## Choose random county and check the distance (It was right)
set.seed(123)  # for reproducibility
random_indices <- sample(1:nrow(map_data), 2)
county1 <- map_data$County[random_indices[1]]
county2 <- map_data$County[random_indices[2]]

# Extract distance from the matrix
distance_km <- dist_matrix[toupper(county1), toupper(county2)]

# Print result
cat("Distance between", county1, "and", county2, "is", round(distance_km, 2), "km\n")



