

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


# Set file path
mmr_path <- "RawData/2023-2024_School_Vaccination_Coverage_Levels_Kindergarten.xlsx"

# List available sheets
sheets <- excel_sheets(mmr_path)
print(sheets)

# Read the School District sheet
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



##################### Plot TEST


texas_flows_named <- read_csv("ProcessedData/texas_county_flows_2019.csv")

# Filter for rows where county_o is "GAINES"
gaines_flows <- texas_flows_named %>%
  filter(county_o == "GAINES") %>%
  group_by(county_d) %>%
  summarise(
    flow = sum(as.numeric(pop_flows)),
    .groups = "drop"
  )

# Load Texas counties directly from TIGRIS
tx_counties <- counties(state = "TX", cb = TRUE, year = 2020) %>% 
  st_transform(crs = 4326)  # ensure lat/lon projection


map_data_gaines <- tx_counties %>%
  mutate(
    NAME = toupper(NAME)
  ) %>%
  left_join(
    gaines_flows, 
    by = c("NAME" = "county_d"))

#  Plot
g <- ggplot(map_data_gaines) +
  geom_sf(fill = "white", color = "gray40", size = 0.1) +
  geom_sf_text(
    aes(label = ifelse(!is.na(flow) & flow > 0, round(flow), "")),
    size = 2.5,
    color = "black"
  ) +
  theme_minimal() +
  labs(
    title = "Population Flow from GAINES County (Text Labels)",
    subtitle = "2019 Flow values displayed as text",
    caption = "Source: texas_flows_named"
  ) +
  theme(
    plot.title = element_text(size = 16, face = "bold"),
    plot.subtitle = element_text(size = 12),
    axis.text = element_blank(),
    axis.ticks = element_blank()
  )

ggsave("Figures/Sum_weekly_flow.png", 
       g,
       width = 12,
       height = 8,
       dpi = 400)


####################################






names(texas_flows_all)
# Filter flows and remove intra-county
flow_subset <- texas_flows_all %>%
  mutate(pop_flows = as.numeric(pop_flows)) %>%
  filter(county_o == "GAINES") %>%
  filter(pop_flows > 0, geoid_o != geoid_d)


# Base plot
g <- ggplot() +
  # Black fill only inside Texas counties
  geom_sf(data = tx_counties, fill = "black", color = "gray30", size = 0.1) +
  
  # Flow lines
  geom_curve(data = flow_subset,
             aes(x = as.numeric(lng_o), y = as.numeric(lat_o),
                 xend = as.numeric(lng_d), yend = as.numeric(lat_d),
                 alpha = pop_flows),
             color = "yellow3", curvature = 0, linewidth = 0.4) +
  
  scale_alpha_continuous(range = c(0.01, 0.3), guide = "none") +
  coord_sf(xlim = c(-107, -93), ylim = c(25, 37)) +
  theme_void() +
  theme(
    panel.background = element_rect(fill = "black", color = NA),  # outside background
    plot.background = element_rect(fill = "white", color = NA)
  )

ggsave("Figures/flow_map_texas.png", 
       g,
       width = 12,
       height = 8,
       dpi = 400)


flow_subset_max <- texas_flows_all %>%
  mutate(pop_flows = as.numeric(pop_flows)) %>%
  filter(county_o == "GAINES", geoid_o != geoid_d) %>%
  group_by(geoid_d, county_d, lng_o, lat_o, lng_d, lat_d) %>%
  summarize(max_flow = max(pop_flows, na.rm = TRUE), .groups = "drop")


g <- ggplot() +
  geom_sf(data = tx_counties, fill = "black", color = "gray30", size = 0.1) +
  
  geom_curve(data = flow_subset_max,
             aes(x = as.numeric(lng_o), y = as.numeric(lat_o),
                 xend = as.numeric(lng_d), yend = as.numeric(lat_d),
                 alpha = max_flow),
             color = "yellow3", curvature = 0, linewidth = 1) +
  
  scale_alpha_continuous(range = c(0.1, 0.3), guide = "none") +
  coord_sf(xlim = c(-107, -93), ylim = c(25, 37)) +
  theme_void() +
  theme(
    panel.background = element_rect(fill = "black", color = NA),
    plot.background = element_rect(fill = "white", color = NA)
  )

ggsave("Figures/flow_map_texas_max.png", 
       g,
       width = 12,
       height = 8,
       dpi = 400)

