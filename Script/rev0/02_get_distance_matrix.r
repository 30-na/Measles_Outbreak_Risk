

library(tidygeocoder)
library(sf)
library(httr)
library(jsonlite)
library(reshape2)
library(ggplot2)
library(geosphere)

map_data <- readRDS("ProcessedData/map_data.rds")

map_data <- st_as_sf(map_data)  

# Compute centroids of each multipolygon
centroids <- st_centroid(map_data)

# Extract coordinates (longitude = X, latitude = Y)
coords <- st_coordinates(centroids)
map_data$longitude <- coords[, 1]
map_data$latitude <- coords[, 2]


########################## Distance matrix (Haversine)

# Prepare coordinate input
locations <- map_data[, c("longitude", "latitude")] |> st_drop_geometry()

# Compute pairwise Haversine distances in kilometers
n <- nrow(locations)
dist_matrix <- matrix(NA, n, n, dimnames = list(map_data$NAME, map_data$NAME))

pb <- txtProgressBar(min = 0, max = n, style = 3)
for (i in 1:n) {
  for (j in 1:n) {
    dist_matrix[i, j] <- distHaversine(locations[i, ], locations[j, ]) / 1000  # in km
  }
  setTxtProgressBar(pb, i)
}
close(pb) 
# Preview

round(dist_matrix[1:5, 1:5], 2)


saveRDS(dist_matrix, "ProcessedData/distance_matrix_haversine.rds")



############################ (Driving Distance)
compute_distance_matrix <- function(locations, api_key, batch_size = 50, sleep_time = 2) {
  n <- nrow(locations)
  dist_matrix <- matrix(NA, n, n)
  n_batches <- ceiling(n / batch_size)
  total_batches <- n_batches^2
  pb <- txtProgressBar(min = 0, max = total_batches, style = 3)
  counter <- 0
  
  for (i in 1:n_batches) {
    for (j in 1:n_batches) {
      idx_i <- ((i - 1) * batch_size + 1):min(i * batch_size, n)
      idx_j <- ((j - 1) * batch_size + 1):min(j * batch_size, n)
      
      # Combine unique coordinates from idx_i and idx_j
      # Use only square batches to respect ORS limits (<= 50 locations)
      if (length(idx_i) != length(idx_j)) next  # skip non-square batches
      
      loc_batch <- data.frame(
        lon = as.numeric(locations[idx_i, 1]),
        lat = as.numeric(locations[idx_i, 2])
      )
      
      body <- list(
        locations = unname(as.matrix(loc_batch)),
        metrics = list("distance"),
        units = "km"
      )
      
      success <- FALSE
      attempt <- 1
      
      while (!success && attempt <= 3) {
        response <- tryCatch({
          resp <- POST(
            url = "https://api.openrouteservice.org/v2/matrix/driving-car",
            add_headers("Authorization" = api_key, "Content-Type" = "application/json"),
            body = toJSON(body, auto_unbox = TRUE),
            encode = "json"
          )
          cat("Status code:", resp$status_code, "\n")
          resp
        }, error = function(e) {
          message("POST failed: ", e$message)
          return(NULL)
        })
        
        if (!is.null(response) && response$status_code == 200) {
          result <- fromJSON(content(response, "text", encoding = "UTF-8"))
          dist_matrix[idx_i, idx_i] <- result$distances
          success <- TRUE
        } else {
          message(sprintf("Batch (%d, %d) failed. Retrying in 10s... (Attempt %d)", i, j, attempt))
          Sys.sleep(10)
          attempt <- attempt + 1
        }
      }
      
      counter <- counter + 1
      setTxtProgressBar(pb, counter)
      Sys.sleep(sleep_time)
    }
  }
  
  close(pb)
  return(dist_matrix)
}


# Assuming you've already computed centroids:
locations <- as.matrix(st_drop_geometry(map_data)[, c("longitude", "latitude")])
summary(locations)


# Your OpenRouteService API key
api_key <- '5b3ce3597851110001cf6248252b8924ec4a4a60b3dc78e2af683166'

# Compute full matrix
dist_matrix <- compute_distance_matrix(locations, api_key, batch_size = 5)







# Example for <= 5 districts
body <- list(
  locations = unname(as.data.frame(locations)),
  metrics = list("distance"),
  units = "km"
)

response <- POST(
  url = "https://api.openrouteservice.org/v2/matrix/driving-car",
  add_headers("Authorization" = openrout_service_api, "Content-Type" = "application/json"),
  body = toJSON(body, auto_unbox = TRUE)
)

dist_matrix <- fromJSON(content(response, "text", encoding = "UTF-8"))$distances


###############################
library(httr)
library(jsonlite)

# Use first 5 points only
test_coords <- data.frame(
  lon = locations[1:5, 1],
  lat = locations[1:5, 2]
)


body <- list(
  locations = unname(as.matrix(rbind(loc_i, loc_j))),
  sources = 0:(length(idx_i) - 1),
  destinations = length(idx_i):(length(idx_i) + length(idx_j) - 1),
  metrics = list("distance"),
  units = "km"
)


test_coords <- rbind(loc_i, loc_j)

test_body <- list(
  locations = unname(test_coords),
  metrics = list("distance"),
  units = "km"
)


test_resp <- tryCatch({
  httr::POST(
    url = "https://api.openrouteservice.org/v2/matrix/driving-car",
    add_headers("Authorization" = api_key, "Content-Type" = "application/json"),
    body = toJSON(body, auto_unbox = TRUE),
    encode = "json"
  )
}, error = function(e) {
  message("POST failed: ", e$message)
  return(NULL)
})



if (!is.null(test_resp)) {
  cat("Status code:", test_resp$status_code, "\n")
  cat(content(test_resp, "text"), "\n")
} else {
  cat("Request failed with no response.\n")
}




