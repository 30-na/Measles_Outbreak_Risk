

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

# Set file path
mmr_path <- "RawData/2023-2024_School_Vaccination_Coverage_Levels_Kindergarten.xlsx"

# List available sheets
sheets <- excel_sheets(mmr_path)
print(sheets)

# Read the School District sheet
mmr_df <- read_excel(mmr_path, sheet=1, skip=2)


# Select relevant columns
mmr <- mmr_df %>%
  select(
    District  = "Facility Name",
    address = "Facility Address",
    County,
    MMR
    ) %>%
  mutate(
    District = str_to_upper(District),
    MMR = as.numeric(MMR) 
  )


# Set file path
pop_path <- "RawData/red635_schooldistrict_population_sy2324.xls"

# Read the file
pop_df <- read_excel(pop_path, skip=4)

total_pop <- pop_df %>%
  fill(
    `School District`,
    .direction = "down"
    ) %>%
  filter(
    ...4 == "Total:"
    ) %>%
  select(
    District = `School District`,
    total_population = Total
    )
  

u18_pop <- pop_df %>%
  fill(
    `School District`,
    .direction = "down"
  ) %>%
  filter(
    ...4 == "<18:"
    ) %>%
  select(
    District = `School District`, 
    Under18_Population = Total
    )

# NOTE:
# There is two population under 18 and total. Use both to check the difference


district_pop <- left_join(
  total_pop, 
  u18_pop,
  by = "District"
  ) %>%
  mutate(
    District = str_to_upper(District)
  )%>%
  mutate(
    District = case_when(
      District == "GOOSE CREEK CONSOLIDATED ISD" ~ "GOOSE CREEK CONS ISD",
      TRUE ~ District
      )
    )



## Merge all data 

# NOTE:
# Population data used as based and mmr district name modified to get match with based data (population)
# Even after modification there are 48 District that don't have data at mmr data 

# Check: rows in mmr not in population data
mmr_not_in_pop <- anti_join(
  mmr,
  district_pop,
  by = "District"
  )


# Do some modification to match both dataset
mmr_modified <- mmr %>%
  mutate(
    District = str_replace_all(District, "CISD", "CONS ISD")
    ) %>%
  mutate(
    District = case_when(
      District == "BIG SANDY ISD" ~ paste0(District, " (", toupper(County), ")"),
      District == "CENTERVILLE ISD" ~ paste0(District, " (", toupper(County), ")"),
      District == "CHAPEL HILL ISD" ~ paste0(District, " (", toupper(County), ")"),
      District == "DAWSON ISD" ~ paste0(District, " (", toupper(County), ")"),
      District == "EDGEWOOD ISD" ~ paste0(District, " (", toupper(County), ")"),
      District == "HIGHLAND PARK ISD" ~ paste0(District, " (", toupper(County), ")"),
      District == "HUBBARD ISD" ~ paste0(District, " (", toupper(County), ")"),
      District == "NORTHSIDE ISD" ~ paste0(District, " (", toupper(County), ")"),
      District == "MIDWAY ISD" ~ paste0(District, " (", toupper(County), ")"),
      District == "VALLEY VIEW ISD" ~ paste0(District, " (", toupper(County), ")"),
      District == "WYLIE ISD" ~ paste0(District, " (", toupper(County), ")"),
      
      District == "CROCKETT COUNTY CONSOLIDATED CSD" ~ "CROCKETT COUNTY CONS CSD",
      District == "EAGLE MT-SAGINAW ISD" ~ "EAGLE MOUNTAIN-SAGINAW ISD",
      District == "FT DAVIS ISD" ~ "FORT DAVIS ISD",
      District == "FT HANCOCK ISD" ~ "FORT HANCOCK ISD",
      District == "FT SAM HOUSTON ISD" ~ "FORT SAM HOUSTON ISD",
      District == "GOLD BURG ISD" ~ "GOLD-BURG ISD",
      District == "GOLDTHWAITE ISD" ~ "GOLDTHWAITE CONS ISD",
      District == "KENEDY COUNTY WIDE CSD" ~ "KENEDY COUNTY-WIDE CSD",
      District == "WEST RUSK COUNTY CONSOLIDATED ISD" ~ "WEST RUSK COUNTY CONS ISD",
      District == "SCHERTZ-CIBOLO-U CITY ISD" ~ "SCHERTZ-CIBOLO-UNIVERSAL CITY ISD",
      
      TRUE ~ District)
    )
  

# Check: rows in population data not in mmr
pop_not_in_mmr <- anti_join(
  district_pop,
  mmr_modified, 
  by = "District"
)

print(pop_not_in_mmr, n=50)

# Fuzzy match those unmatched
fuzzy_matches <- stringdist_join(
  pop_not_in_mmr,
  mmr_modified,
  by = "District",
  method = "jw",        # Jaro-Winkler
  max_dist = 0.2,       # Adjust threshold as needed
  distance_col = "dist"
) %>%
  group_by(District.x) %>%
  slice_min(
    order_by = dist,
    n = 1
  ) %>%
  ungroup() %>%
  select(District.x, District.y, dist)

# View top suggestions
print(fuzzy_matches, n = 50)



# Merge population and mmr datasets
merged_data <- left_join(
  district_pop,
  mmr_modified,
  by = "District"
  ) %>%
  filter(
    !is.na(MMR)
    )  # or !is.na(address)

# check the unmatch data()
sum(is.na(merged_data$address))


# Load Texas school district boundaries
# https://cran.r-project.org/web/packages/tigris/tigris.pdf?utm_source=chatgpt.com

options(tigris_use_cache = F)

tx_districts <- school_districts(state = "TX", year = 2024)


## https://catalog.data.gov/dataset/tiger-line-shapefile-2018-state-texas-current-unified-school-districts-shapefile-state-based
# tx_district has missing data so i got it from above link and made dummy data and add it to map (NOT THERE TOO) further investigation need


VYSEHRAD_districts <- st_read("RawData/tl_2018_48_unsd/tl_2018_48_unsd.shp")
vysehrad_geom <- VYSEHRAD_districts %>% 
  filter(
    str_detect(toupper(NAME), "VYSEHRAD")
    )


names(tx_districts)
tx_districts <- tx_districts %>%
  mutate(
    NAME = toupper(NAME),
    INTPTLAT = as.numeric(INTPTLAT),
    INTPTLON = as.numeric(INTPTLON)
    ) %>% 
  mutate(
    NAME = str_replace_all(NAME, "INDEPENDENT SCHOOL DISTRICT", "ISD"),
    NAME = str_replace_all(NAME, "CONSOLIDATED", "CONS"),
    NAME = str_replace_all(NAME, "COMMON SCHOOL DISTRICT", "CSD")
  ) %>%
  mutate(
    NAME = case_when(
      NAME == "LAPOYNER ISD" ~ "LAPOYNOR ISD",
      NAME == "STAFFORD MUNICIPAL SCHOOL DISTRICT" ~ "STAFFORD MSD",
      
      (NAME == "BIG SANDY ISD" & GEOID == 4810170) ~ "BIG SANDY ISD (POLK)",
      (NAME == "BIG SANDY ISD" & GEOID == 4810140) ~ "BIG SANDY ISD (UPSHUR)",
      
      (NAME == "CENTERVILLE ISD" & GEOID == 4813410) ~ "CENTERVILLE ISD (LEON)",
      (NAME == "CENTERVILLE ISD" & GEOID == 4813380) ~ "CENTERVILLE ISD (TRINITY)",
      
      (NAME == "CHAPEL HILL ISD" & GEOID == 4813680) ~ "CHAPEL HILL ISD (TITUS)",
      (NAME == "CHAPEL HILL ISD" & GEOID == 4813650) ~ "CHAPEL HILL ISD (SMITH)",
      
      (NAME == "MIDWAY ISD" & GEOID == 4830660) ~ "MIDWAY ISD (CLAY)",
      (NAME == "MIDWAY ISD" & GEOID == 4830640) ~ "MIDWAY ISD (MCLENNAN)",
      
      (NAME == "NORTHSIDE ISD" & GEOID == 4833090) ~ "NORTHSIDE ISD (WILBARGER)",
      (NAME == "NORTHSIDE ISD" & GEOID == 4833120) ~ "NORTHSIDE ISD (BEXAR)",
      
      (NAME == "VALLEY VIEW ISD" & GEOID == 4843800) ~ "VALLEY VIEW ISD (HIDALGO)",
      (NAME == "VALLEY VIEW ISD" & GEOID == 4843860) ~ "VALLEY VIEW ISD (COOKE)",
      
      (NAME == "DAWSON ISD" & GEOID == 4816350) ~ "DAWSON ISD (DAWSON)",
      (NAME == "DAWSON ISD" & GEOID == 4816380) ~ "DAWSON ISD (NAVARRO)",
      
      (NAME == "WYLIE ISD" & GEOID == 4846530) ~ "WYLIE ISD (COLLIN)",
      (NAME == "WYLIE ISD" & GEOID == 4846500) ~ "WYLIE ISD (TAYLOR)",
      
      (NAME == "EDGEWOOD ISD" & GEOID == 4818150) ~ "EDGEWOOD ISD (BEXAR)",
      (NAME == "EDGEWOOD ISD" & GEOID == 4818120) ~ "EDGEWOOD ISD (VAN ZANDT)",
      
      (NAME == "HUBBARD ISD" & GEOID == 4823730) ~ "HUBBARD ISD (HILL)",
      (NAME == "HUBBARD ISD" & GEOID == 4823700) ~ "HUBBARD ISD (BOWIE)",
      
      (NAME == "HIGHLAND PARK ISD" & GEOID == 4835560) ~ "HIGHLAND PARK ISD (POTTER)",
      (NAME == "HIGHLAND PARK ISD" & GEOID == 4823250) ~ "HIGHLAND PARK ISD (DALLAS)",
      
      TRUE ~ NAME
    )
  )
  

map_data <- left_join(
  tx_districts,
  merged_data, 
  by = c("NAME" = "District"))


saveRDS(map_data, "ProcessedData/map.rds")





##### Plots

g_map <- ggplot(map_data) +
  geom_sf(fill = "lightblue", color = "gray30", size = 0.1) +
  geom_sf_text(aes(label = NAME), size = .7, check_overlap = TRUE) +
  theme_minimal() +
  labs(
    title = "Texas Independent School Districts",
    subtitle = "Labeled by District Name",
    caption = "Source: TEA / Your Processed Map Data"
  )

ggsave("Figures/district_map.png", 
       plot = g_map, width = 20, height = 20, dpi = 900)


##################################################
ggplot(tx_districts) +
  geom_sf(fill = "lightblue", color = "gray30", size = 0.1) +
  theme_minimal() +
  labs(
    title = "Texas Independent School Districts (2024)",
    caption = "Source: US Census TIGER/Line"
  )
#########################################################
names(map_data)
# Get unmatched districts from map_data
unmatched_districts <- map_data %>%
  filter(is.na(MMR)) %>%
  select(NAME) %>%
  st_drop_geometry() %>%
  distinct()

# Find closest matches from merged_data
fuzzy_matches <- stringdist_join(
  unmatched_districts,
  merged_data,
  by = c("NAME" = "District"),
  method = "jw",          # Jaro-Winkler
  max_dist = 0.2,         # Adjust if needed
  distance_col = "dist"
) %>%
  group_by(
    NAME
    ) %>%
  slice_min(
    order_by = dist,
    n = 1
    ) %>%
  ungroup() %>%
  select(
    Unmatched = NAME,
    Closest_Match = District, 
    dist
    )



# Show results
print(fuzzy_matches, n = 70)

# Vector of unmatched district names
unmatched_names <- map_data %>%
  filter(is.na(MMR)) %>%
  pull(NAME)

# All district names
all_names <- merged_data$District

# Compute closest match for each unmatched name
matched_names <- map_chr(unmatched_names, function(name) {
  distances <- stringdist::stringdist(name, all_names, method = "jw")
  all_names[which.min(distances)]
})

# Combine into a tibble
result <- tibble(
  Unmatched = unmatched_names,
  Closest_Match = matched_names
)

print(result, n=70)



