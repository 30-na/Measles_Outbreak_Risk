
# Load necessary libraries
library(readxl)
library(dplyr)
library(tidyr)
library(stringr)
library(fuzzyjoin)
library(tigris)
library(sf)
library(stringdist)
library(purrr)
library(ggplot2)
library(readr)
library(geosphere)

# Set file path
mmr_path <- "RawData/2023-2024_School_Vaccination_Coverage_Levels_Kindergarten.xlsx"

# List available sheets
sheets <- excel_sheets(mmr_path)
print(sheets)

# Read the School County sheet
mmr_county <- read_excel(mmr_path, sheet=2, skip=0)


# Select relevant columns
mmr <- mmr_county %>%
  select(
    County,
    MMR
  ) %>%
  mutate(
    MMR = as.numeric(MMR),
    County = case_when(
      County == "Dewitt" ~ "DeWitt",
      County == "Mcculloch" ~ "McCulloch",
      County == "Mclennan" ~ "McLennan",
      County == "Mcmullen" ~ "McMullen",
      TRUE ~ County
    )
  )


tx_counties <- counties(state = "TX", year = 2024)


tx_counties_mmr <- left_join(
  tx_counties,
  mmr, 
  by = c("NAME" = "County")
) %>% 
  mutate(
    County = toupper(NAME)
  ) %>%
  select(
    County,
    MMR
  )



county_pop <- read_csv("rawData/Counties_Pop.txt") %>%
  mutate(
    FENAME = case_when(
      FENAME == "DE WITT" ~ "DEWITT",
      TRUE ~ FENAME)
  )


tx_counties_map <- tx_counties_mmr %>%
  left_join(
    county_pop,
    by = c("County" = "FENAME")
  ) %>%
  select(
    County,
    MMR,
    total
  )


saveRDS(tx_counties_map, "ProcessedData/map_county.rds")


###################### Population Flow
# Download code for the county2county weeky flow for 2019
# python download_weekly_data.py --start_year 2019 --start_month 1 --start_day 7 --end_year 2019 --end_month 12 --end_day 30 --output_folder weekly_flows  --county 



# Folder with weekly CSVs
file_list <- list.files("RawData/weekly_flows/county2county/", full.names = TRUE, pattern = "\\.csv$")

# Helper function to filter for Texas-only flows
is_texas_fips <- function(fips) {
  substr(fips, 1, 2) == "48"
}

# Read and filter all files
texas_flows <- lapply(file_list, function(file) {
  df <- read_csv(file, col_types = cols(.default = "c"))  # read as character to preserve leading 0s
  df <- df %>%
    filter(is_texas_fips(geoid_o) & is_texas_fips(geoid_d))
})

# Combine into one data frame
texas_flows_all <- bind_rows(texas_flows)

# Get county names from tigris
tx_counties <- counties(state = "TX", cb = TRUE, year = 2020) %>%
  st_drop_geometry() %>%
  select(GEOID, NAME)



# Add county names to origin and destination
texas_flows_named <- texas_flows_all %>%
  left_join(tx_counties, by = c("geoid_o" = "GEOID")) %>%
  rename(county_o = NAME) %>%
  left_join(tx_counties, by = c("geoid_d" = "GEOID")) %>%
  rename(county_d = NAME) %>%
  mutate(
    county_o = toupper(county_o),
    county_d = toupper(county_d)
  )

# Save final result with county names
write_csv(texas_flows_named, "ProcessedData/texas_county_flows_2019.csv")

################## Distance Matrix

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


