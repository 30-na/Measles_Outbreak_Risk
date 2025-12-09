# ------------------ Libraries ------------------
library(readxl)
library(dplyr)
library(purrr)
library(ggplot2)
library(viridis)
library(tigris)
library(stringdist)
library(stringr)
library(sf)
library(tidyr)




# Functions ----

fuzzy_match <- function(col1, col2, threshold = 0.05, method = "jw") {
  vals1 <- col1 %>%
    as.character() %>%
    str_trim() %>%
    toupper() %>%
    unique()
  
  vals2 <- col2 %>%
    as.character() %>%
    str_trim() %>%
    toupper() %>%
    unique()
  
  unmatched1 <- setdiff(vals1, vals2)
  unmatched2 <- setdiff(vals2, vals1)
  
  expand.grid(column01 = unmatched1, column02 = unmatched2, stringsAsFactors = FALSE) %>%
    mutate(dist = stringdist(column01, column02, method = method)) %>%
    group_by(column01) %>%
    slice_min(dist, n = 1) %>%
    ungroup() %>%
    filter(dist < threshold) %>%
    arrange(dist)
}

keep_cols_district <- function(df) {
  df %>%
    dplyr::filter(
      grepl("public", `School Type`, ignore.case = TRUE)
    ) %>%
    select(
      id = "Facility Number",
      district = "Facility Name",
      county = matches("county", ignore.case = TRUE),
      MMR
    ) %>%
    mutate(
      MMR = as.numeric(MMR),
      county = toupper(county),
      district = toupper(district)
    ) %>%
    mutate(
      district = case_when(
        district == "ARANSAS COUNTY ISD" & county == "ARANSAS" ~ "ROCKPORT-FULTON ISD",
        district == "EDGEWOOD ISD (BEXAR)" & county == "BEXAR" ~ "EDGEWOOD ISD",
        district == "NORTHSIDE ISD - BEXAR" & county == "BEXAR" ~ "NORTHSIDE ISD",
        district == "HUBBARD ISD - DEKALB" & county == "BOWIE" ~ "HUBBARD ISD",
        district == "MIDWAY ISD - HENRIETTA" & county == "CLAY" ~ "MIDWAY ISD",
        district == "WYLIE ISD - COLLIN" & county == "COLLIN" ~ "WYLIE ISD",
        district == "DAWSON ISD - WELCH" & county == "DAWSON" ~ "DAWSON ISD",
        district == "LA GRANGEISD" & county == "FAYETTE" ~ "LA GRANGE ISD",
        district == "FLOYDADA ISD" & county == "FLOYD" ~ "FLOYDADA COLLEGIATE ISD",
        district == "PRINGLE-MORSE ISD" & county == "HANSFORD" ~ "PRINGLE-MORSE CISD",
        district == "VALLEY VIEW ISD-PHARR" & county == "HIDALGO" ~ "VALLEY VIEW ISD",
        district == "HAMLIN COLLEGIATE ISD" & county == "JONES" ~ "HAMLIN ISD",
        district == "KNOX CITY-OBRIEN CISD" & county == "KNOX" ~ "KNOX CITY-O'BRIEN CISD",
        district == "SUNRAY ISD" & county == "MOORE" ~ "SUNRAY COLLEGIATE ISD",
        district == "ROSCOE ISD" & county == "NOLAN" ~ "ROSCOE COLLEGIATE ISD",
        district == "IRAAN-SHEFFIELD ISD" & county == "PECOS" ~ "IRAAN-SHEFFIELD COLLEGIATE ISD",
        district == "HIGHLAND PARK ISD - AMARILLO" & county == "POTTER" ~ "HIGHLAND PARK ISD",
        district == "WEST RUSK ISD" & county == "RUSK" ~ "WEST RUSK COUNTY CONSOLIDATED ISD",
        district == "RIO GRANDE CITY CISD" & county == "STARR" ~ "RIO GRANDE CITY GRULLA ISD",
        district == "THROCKMORTON ISD" & county == "THROCKMORTON" ~ "THROCKMORTON COLLEGIATE ISD",
        district == "CENTERVILLE ISD - GROVETON" & county == "TRINITY" ~ "CENTERVILLE ISD",
        district == "EDGEWOOD ISD (VAN ZANDT)" & county == "VAN ZANDT" ~ "EDGEWOOD ISD",
        district == "CUMBY ISD" & county == "HOPKINS" ~ "CUMBY COLLEGIATE ISD",
        TRUE ~ district
      ),
      county = case_when(
        district == "LAGO VISTA ISD" & county == "WILLIAMSON" ~ "TRAVIS",
        TRUE ~ county
      )
    )
}

keep_cols_county <- function(df) {
  df %>%
    select(
      county = County,
      MMR
    ) %>%
    mutate(
      MMR = as.numeric(MMR),
      county = toupper(county)
    ) %>%
    filter(
      county != "TEXAS"
    )
}


# Load Texas Map ----

# Download 2024 Texas Independent School District boundaries
tx_school_district_map <- school_districts(
  state = "TX",
  year = 2024,
  type = "unified"
)

tx_counties <- counties(
  state = "TX",
  year = 2024
) %>%
  mutate(
    county = toupper(NAME)
  ) %>%
  select(
    county
  )

tx_school_district_map_clean <- tx_school_district_map %>%
  mutate(
    district = toupper(NAME),
    district = str_replace(district, " INDEPENDENT SCHOOL DISTRICT", " ISD")
  )%>%
  mutate(
    district = case_when(
      district == "EDGEWOOD ISD"      & GEOID == "4818150" ~ paste0(district, " (BEXAR)"),
      district == "EDGEWOOD ISD"      & GEOID == "4818120" ~ paste0(district, " (VAN ZANDT)"),
      district == "MIDWAY ISD"        & GEOID == "4830660" ~ paste0(district, " (CLAY)"),
      district == "MIDWAY ISD"        & GEOID == "4830640" ~ paste0(district, " (MCLENNAN)"),
      district == "CENTERVILLE ISD"   & GEOID == "4813380" ~ paste0(district, " (TRINITY)"),
      district == "CENTERVILLE ISD"   & GEOID == "4813410" ~ paste0(district, " (LEON)"),
      district == "CHAPEL HILL ISD"   & GEOID == "4813680" ~ paste0(district, " (TITUS)"),
      district == "CHAPEL HILL ISD"   & GEOID == "4813650" ~ paste0(district, " (SMITH)"),
      district == "BIG SANDY ISD"     & GEOID == "4810170" ~ paste0(district, " (POLK)"),
      district == "BIG SANDY ISD"     & GEOID == "4810140" ~ paste0(district, " (UPSHUR)"),
      district == "DAWSON ISD"        & GEOID == "4816350" ~ paste0(district, " (DAWSON)"),
      district == "DAWSON ISD"        & GEOID == "4816380" ~ paste0(district, " (NAVARRO)"),
      district == "NORTHSIDE ISD"     & GEOID == "4833090" ~ paste0(district, " (WILBARGER)"),
      district == "NORTHSIDE ISD"     & GEOID == "4833120" ~ paste0(district, " (BEXAR)"),
      district == "HIGHLAND PARK ISD" & GEOID == "4835560" ~ paste0(district, " (POTTER)"),
      district == "HIGHLAND PARK ISD" & GEOID == "4823250" ~ paste0(district, " (DALLAS)"),
      district == "WYLIE ISD"         & GEOID == "4846530" ~ paste0(district, " (TAYLOR)"),
      district == "WYLIE ISD"         & GEOID == "4846500" ~ paste0(district, " (COLLIN)"),
      district == "HUBBARD ISD"       & GEOID == "4823730" ~ paste0(district, " (HILL)"),
      district == "HUBBARD ISD"       & GEOID == "4823700" ~ paste0(district, " (BOWIE)"),
      district == "VALLEY VIEW ISD"   & GEOID == "4843800" ~ paste0(district, " (HIDALGO)"),
      district == "VALLEY VIEW ISD"   & GEOID == "4843860" ~ paste0(district, " (COOKE)"),
      TRUE ~ district
    )
  ) %>%
  select(
    district,
    geometry
  )



dup_rows <- tx_school_district_map_clean %>%
  group_by(district) %>%
  filter(n() > 1) %>%
  ungroup()



# Load Texas age pyramid ----

pyramid_age_data <- read.delim(
  "RawData/Texas.csv",
  fileEncoding = "UTF-16LE",
  check.names = FALSE
)


grade_weights <- pyramid_age_data %>%
  distinct(AGE, 
           Female, 
           Male, 
           .keep_all = TRUE
  ) %>%
  mutate(
    Female = as.numeric(gsub(",", "", Female)),
    Male   = as.numeric(gsub(",", "", Male)),
    age = as.numeric(AGE)
  ) %>%
  group_by(
    age,
    Year
  ) %>%
  summarise(
    Female = sum(Female, na.rm = TRUE),
    Male   = sum(Male,   na.rm = TRUE),
    .groups = "drop"
  ) %>%
  mutate(
    population = Female + Male
  ) %>%
  filter(
    age >= 5 & age <= 17
  ) %>%
  mutate(
    pop_proportion = population / sum(population)
  ) %>% 
  mutate(
    grade = case_when(
      age == 5  ~ "K",
      age == 6  ~ "G1",
      age == 7  ~ "G2",
      age == 8  ~ "G3",
      age == 9  ~ "G4",
      age == 10 ~ "G5",
      age == 11 ~ "G6",
      age == 12 ~ "G7",
      age == 13 ~ "G8",
      age == 14 ~ "G9",
      age == 15 ~ "G10",
      age == 16 ~ "G11",
      age == 17 ~ "G12"
    )
  ) %>%
  select(
    grade, 
    pop_proportion
  )%>%
  tibble::deframe()




# MMR coverage for public schools at District level ----

# File paths
Kindergarten_path_2025 <- "RawData/2024-2025_School_Vaccination_Coverage_Levels_Kindergarten.xlsx"
Kindergarten_path_2024 <- "RawData/2023-2024_School_Vaccination_Coverage_Levels_Kindergarten.xlsx"
Kindergarten_path_2023 <- "RawData/22-23-School-Vaccination-Coverage-by-District-and-County-K.xlsx"
Kindergarten_path_2022 <- "RawData/2021-2022-School-Vaccination-Coverage-by-District-and-County-Kindergarten.xlsx"
Kindergarten_path_2021 <- "RawData/2020-2021-School-Vaccination-Coverage-Levels-by-District-Private-School-and-County---Kindergarten.xlsx"
Kindergarten_path_2020 <- "RawData/2019-2020-School-Vaccination-Coverage-Levels---Kindergarten.xlsx"

Seventh_path_2025 <- "RawData/2024-2025_School_Vaccination_Coverage_Levels_Seventh_Grade.xlsx"
Seventh_path_2024 <- "RawData/2023-2024_School_Vaccination_Coverage_Levels_Seventh_Grade.xlsx"
Seventh_path_2023 <- "RawData/22-23-School-Vaccination-Coverage-by-District-and-County-Seventh-Grade.xlsx"
Seventh_path_2022 <- "RawData/2021-2022-School-Vaccination-Coverage-by-District-and-County-Seventh-Grade.xlsx"
Seventh_path_2021 <- "RawData/2020-2021-School-Vaccination-Coverage-Levels-by-District-Private-School-and-County---Seventh-Grade.xlsx"
Seventh_path_2020 <- "RawData/2019-2020-School-Vaccination-Coverage-Levels-by-District-Private-School-and-County---Seventh-Grade.xlsx"


#  Read and clean 
district_K_2025 <- read_excel(Kindergarten_path_2025, sheet = "Coverage by District", skip = 2)
district_K_2024 <- read_excel(Kindergarten_path_2024, sheet = "Coverage by District", skip = 2)
district_K_2023 <- read_excel(Kindergarten_path_2023, sheet = "Coverage by District", skip = 2)
district_K_2022 <- read_excel(Kindergarten_path_2022, sheet = "Coverage by District", skip = 2)
district_K_2021 <- read_excel(Kindergarten_path_2021, sheet = "Coverage by District", skip = 2)
district_K_2020 <- read_excel(Kindergarten_path_2020, sheet = "Coverage by District", skip = 2)

district_7th_2025 <- read_excel(Seventh_path_2025, sheet = "Coverage by District", skip = 2)
district_7th_2024 <- read_excel(Seventh_path_2024, sheet = "Coverage by District", skip = 2)
district_7th_2023 <- read_excel(Seventh_path_2023, sheet = "Coverage by District", skip = 2)
district_7th_2022 <- read_excel(Seventh_path_2022, sheet = "Coverage by District", skip = 2)
district_7th_2021 <- read_excel(Seventh_path_2021, sheet = "Coverage by District", skip = 2)
district_7th_2020 <- read_excel(Seventh_path_2020, sheet = "Coverage by District", skip = 2)


#  Keep only relevant columns 
district_K_2025_mmr <- keep_cols_district(district_K_2025)
district_K_2024_mmr <- keep_cols_district(district_K_2024)
district_K_2023_mmr <- keep_cols_district(district_K_2023)
district_K_2022_mmr <- keep_cols_district(district_K_2022)
district_K_2021_mmr <- keep_cols_district(district_K_2021)
district_K_2020_mmr <- keep_cols_district(district_K_2020)

district_7th_2025_mmr <- keep_cols_district(district_7th_2025)
district_7th_2024_mmr <- keep_cols_district(district_7th_2024)
district_7th_2023_mmr <- keep_cols_district(district_7th_2023)
district_7th_2022_mmr <- keep_cols_district(district_7th_2022)
district_7th_2021_mmr <- keep_cols_district(district_7th_2021)
district_7th_2020_mmr <- keep_cols_district(district_7th_2020)


#  Merge all years 
mmr_district_df <- district_K_2025_mmr %>% rename(K = MMR) %>%
  full_join(district_K_2024_mmr %>% rename(G1 = MMR), by = c("id", "district", "county")) %>%
  full_join(district_K_2023_mmr %>% rename(G2 = MMR), by = c("id", "district", "county")) %>%
  full_join(district_K_2022_mmr %>% rename(G3 = MMR), by = c("id", "district", "county")) %>%
  full_join(district_K_2021_mmr %>% rename(G4 = MMR), by = c("id", "district", "county")) %>%
  full_join(district_K_2020_mmr %>% rename(G5 = MMR), by = c("id", "district", "county")) %>%
  full_join(district_7th_2025_mmr %>% rename(G7 = MMR), by = c("id", "district", "county")) %>%
  full_join(district_7th_2024_mmr %>% rename(G8 = MMR), by = c("id", "district", "county")) %>%
  full_join(district_7th_2023_mmr %>% rename(G9 = MMR), by = c("id", "district", "county")) %>%
  full_join(district_7th_2022_mmr %>% rename(G10 = MMR), by = c("id", "district", "county")) %>%
  full_join(district_7th_2021_mmr %>% rename(G11 = MMR), by = c("id", "district", "county")) %>%
  full_join(district_7th_2020_mmr %>% rename(G12 = MMR), by = c("id", "district", "county")) %>%
  mutate(
    G6 = rowMeans(select(., G5, G7), na.rm = TRUE)
  ) %>%
  select(id, district, county, K, G1, G2, G3, G4, G5, G6, G7, G8, G9, G10, G11, G12)



mmr_district_df_clean <- mmr_district_df %>%
  mutate(
    district = str_replace_all(district, "\\bCISD\\b", "CONSOLIDATED ISD")
  ) %>%
  mutate(
    district = case_when(
      district == "EDGEWOOD ISD"      ~ paste0(district, " (", toupper(county), ")"),
      district == "MIDWAY ISD"        ~ paste0(district, " (", toupper(county), ")"),
      district == "CENTERVILLE ISD"   ~ paste0(district, " (", toupper(county), ")"),
      district == "CHAPEL HILL ISD"   ~ paste0(district, " (", toupper(county), ")"),
      district == "BIG SANDY ISD"     ~ paste0(district, " (", toupper(county), ")"),
      district == "DAWSON ISD"        ~ paste0(district, " (", toupper(county), ")"),
      district == "NORTHSIDE ISD"     ~ paste0(district, " (", toupper(county), ")"),
      district == "HIGHLAND PARK ISD" ~ paste0(district, " (", toupper(county), ")"),
      district == "WYLIE ISD"         ~ paste0(district, " (", toupper(county), ")"),
      district == "HUBBARD ISD"       ~ paste0(district, " (", toupper(county), ")"),
      district == "VALLEY VIEW ISD"   ~ paste0(district, " (", toupper(county), ")"),
      TRUE ~ district
    )
  )  %>% 
  mutate(
    district = case_when(
      district == "GOLD BURG ISD" ~ "GOLD-BURG ISD",
      district == "STAFFORD MSD" ~ "STAFFORD MUNICIPAL SCHOOL DISTRICT",
      district == "GOLDTHWAITE ISD" ~ "GOLDTHWAITE CONSOLIDATED ISD",
      district == "FT HANCOCK ISD" ~ "FORT HANCOCK ISD",
      district == "RAMIREZ CSD" ~ "RAMIREZ COMMON SCHOOL DISTRICT",
      district == "TERLINGUA CSD" ~ "TERLINGUA COMMON SCHOOL DISTRICT",
      district == "DOSS CONSOLIDATED CSD" ~ "DOSS CONSOLIDATED COMMON SCHOOL DISTRICT",
      district == "LAPOYNOR ISD" ~ "LAPOYNER ISD",
      district == "FT DAVIS ISD" ~ "FORT DAVIS ISD",
      district == "FT SAM HOUSTON ISD" ~ "FORT SAM HOUSTON ISD",
      district == "SCHERTZ-CIBOLO-U CITY ISD" ~ "SCHERTZ-CIBOLO-UNIVERSAL CITY ISD",
      district == "CROCKETT COUNTY CONSOLIDATED CSD" ~ "CROCKETT COUNTY CONSOLIDATED COMMON SCHOOL DISTRICT",
      district == "GUTHRIE CSD" ~ "GUTHRIE COMMON SCHOOL DISTRICT",
      district == "EAGLE MT-SAGINAW ISD" ~ "EAGLE MOUNTAIN-SAGINAW ISD",
      district == "KENEDY COUNTY WIDE CSD" ~ "KENEDY COUNTY-WIDE COMMON SCHOOL DISTRICT",
      TRUE ~ district
    )
  ) %>%
  rowwise() %>%
  mutate(
    mean_val = mean(c_across(K:G12), na.rm = TRUE),
    across(K:G12, ~ ifelse(is.na(.x), mean_val, .x)),
    mmr = rowSums(
      cbind(
        K  * grade_weights["K"],
        G1 * grade_weights["G1"],
        G2 * grade_weights["G2"],
        G3 * grade_weights["G3"],
        G4 * grade_weights["G4"],
        G5 * grade_weights["G5"],
        G6 * grade_weights["G6"],
        G7 * grade_weights["G7"],
        G8 * grade_weights["G8"],
        G9 * grade_weights["G9"],
        G10 * grade_weights["G10"],
        G11 * grade_weights["G11"],
        G12 * grade_weights["G12"]
      )
    )
  ) %>%
  select(
    district,
    county,
    mmr
  )

# duplicate_id <- mmr_district_df %>%
#   group_by(id) %>%
#   filter(n() > 1) %>%
#   arrange(id)




# Load District Population

district_pop_path <- "RawData/red635_schooldistrict_population_sy2324.xlsx"
pop_raw <- read_excel(district_pop_path, sheet = "Red635C", col_names = FALSE)

pop_district_df <- pop_raw %>%
  select(1, 4, 5) %>%                                 
  setNames(c("district", "marker", "value")) %>%      
  tidyr::fill(district) %>%                           
  filter(
    marker %in% c("Total:", "<18:")
  ) %>%         
  mutate(
    value = as.numeric(value),
    district = toupper(district),
    marker = ifelse(marker == "Total:", "pop_total", "pop_under18"),
    district = str_replace_all(district, "\\bCONS ISD\\b", "CISD")
  ) %>%
  pivot_wider(
    names_from = marker,
    values_from = value
  )

sum(duplicated(pop_district_df$district))

pop_district_df_clean <- pop_district_df %>%
  mutate(
    district = str_replace(district, "CISD", "CONSOLIDATED ISD"),
    district = str_replace(district, "CSD", "COMMON SCHOOL DISTRICT")
  ) %>%
  mutate(
    district = case_when(
      district == "STAFFORD MSD" ~ "STAFFORD MUNICIPAL SCHOOL DISTRICT",
      district == "DOSS CONS COMMON SCHOOL DISTRICT" ~ "DOSS CONSOLIDATED COMMON SCHOOL DISTRICT",
      district == "LAPOYNOR ISD" ~ "LAPOYNER ISD",
      district == "CROCKETT COUNTY CONS COMMON SCHOOL DISTRICT" ~ "CROCKETT COUNTY CONSOLIDATED COMMON SCHOOL DISTRICT",
      district == "HAMLIN COLLEGIATE ISD" ~ "HAMLIN ISD",
      TRUE ~ district
    )
  )



# Merge mmr and population District with map and plot

merged_map_district_df <- tx_school_district_map_clean %>%
  left_join(
    pop_district_df_clean,
    by = "district"
  ) %>%
  left_join(
    mmr_district_df_clean,
    by = "district"
  )


# saveRDS(merged_map_district_df, "ProcessedData/merged_map_district_df.rds")



p_district_mmr_weighted_map <- ggplot(merged_map_district_df) +
  geom_sf(
    aes(fill = mmr),
    color = "gray50",
    size = 0.1
  ) + 
  geom_sf(
    data = tx_counties,
    fill = NA,
    color = "gray20",
    size = 0.4
  ) + 
  scale_fill_gradientn(
    colors = c("#d73027", "#fee08b", "#1a9850"),  # soft red → light yellow → muted green
    limits = c(0.5, 1.0),
    na.value = "lightgray",
    labels = scales::percent_format(accuracy = 1)
  ) +
  theme_minimal() +
  labs(
    title = "Texas School Districts — MMR Coverage (2024)",
    fill = "MMR"
  )

# save to file
# ggsave("Figures/district_mmr_weighted_map.png",
#        plot = p_district_mmr_weighted_map,
#        width = 10,
#        height = 8, 
#        dpi = 300
# )
# 




# MMR coverage for public schools at County level ----

# Read and clean 
county_K_2025 <- read_excel(Kindergarten_path_2025, sheet = "Coverage by County", skip = 2)
county_K_2024 <- read_excel(Kindergarten_path_2024, sheet = "Coverage by County", skip = 0)
county_K_2023 <- read_excel(Kindergarten_path_2023, sheet = "Coverage by County", skip = 2)
county_K_2022 <- read_excel(Kindergarten_path_2022, sheet = "Coverage by County", skip = 2)
county_K_2021 <- read_excel(Kindergarten_path_2021, sheet = "Coverage by County", skip = 2)
county_K_2020 <- read_excel(Kindergarten_path_2020, sheet = "Coverage by County", skip = 2)

county_7th_2025 <- read_excel(Seventh_path_2025, sheet = "Coverage by County", skip = 2)
county_7th_2024 <- read_excel(Seventh_path_2024, sheet = "Coverage by County", skip = 2)
county_7th_2023 <- read_excel(Seventh_path_2023, sheet = "Coverage by County", skip = 2)
county_7th_2022 <- read_excel(Seventh_path_2022, sheet = "Coverage by County", skip = 2)
county_7th_2021 <- read_excel(Seventh_path_2021, sheet = "Coverage by County", skip = 2)
county_7th_2020 <- read_excel(Seventh_path_2020, sheet = "Coverage by County", skip = 2)


county_K_2025_mmr <- keep_cols_county(county_K_2025)
county_K_2024_mmr <- keep_cols_county(county_K_2024)
county_K_2023_mmr <- keep_cols_county(county_K_2023)
county_K_2022_mmr <- keep_cols_county(county_K_2022)
county_K_2021_mmr <- keep_cols_county(county_K_2021)
county_K_2020_mmr <- keep_cols_county(county_K_2020)

county_7th_2025_mmr <- keep_cols_county(county_7th_2025)
county_7th_2024_mmr <- keep_cols_county(county_7th_2024)
county_7th_2023_mmr <- keep_cols_county(county_7th_2023)
county_7th_2022_mmr <- keep_cols_county(county_7th_2022)
county_7th_2021_mmr <- keep_cols_county(county_7th_2021)
county_7th_2020_mmr <- keep_cols_county(county_7th_2020)


#  Merge all years 

mmr_county_df <- county_K_2025_mmr %>% rename(K  = MMR) %>%
  full_join(county_K_2024_mmr %>% rename(G1 = MMR), by = "county") %>%
  full_join(county_K_2023_mmr %>% rename(G2 = MMR), by = "county") %>%
  full_join(county_K_2022_mmr %>% rename(G3 = MMR), by = "county") %>%
  full_join(county_K_2021_mmr %>% rename(G4 = MMR), by = "county") %>%
  full_join(county_K_2020_mmr %>% rename(G5 = MMR), by = "county") %>%
  
  full_join(county_7th_2025_mmr %>% rename(G7  = MMR), by = "county") %>%
  full_join(county_7th_2024_mmr %>% rename(G8  = MMR), by = "county") %>%
  full_join(county_7th_2023_mmr %>% rename(G9  = MMR), by = "county") %>%
  full_join(county_7th_2022_mmr %>% rename(G10 = MMR), by = "county") %>%
  full_join(county_7th_2021_mmr %>% rename(G11 = MMR), by = "county") %>%
  full_join(county_7th_2020_mmr %>% rename(G12 = MMR), by = "county") %>%
  mutate(
    G6 = rowMeans(select(., G5, G7), na.rm = TRUE)
  ) %>%
  select(county, K, G1, G2, G3, G4, G5, G6, G7, G8, G9, G10, G11, G12)



mmr_county_df_clean <- mmr_county_df %>%
  rowwise() %>%
  mutate(
    mean_val = mean(c_across(K:G12), na.rm = TRUE),
    across(K:G12, ~ ifelse(is.na(.x), mean_val, .x)),
    mmr = rowSums(
      cbind(
        K  * grade_weights["K"],
        G1 * grade_weights["G1"],
        G2 * grade_weights["G2"],
        G3 * grade_weights["G3"],
        G4 * grade_weights["G4"],
        G5 * grade_weights["G5"],
        G6 * grade_weights["G6"],
        G7 * grade_weights["G7"],
        G8 * grade_weights["G8"],
        G9 * grade_weights["G9"],
        G10 * grade_weights["G10"],
        G11 * grade_weights["G11"],
        G12 * grade_weights["G12"]
      )
    )
  ) %>%
  ungroup() %>%
  select(
    county,
    mmr
  )



# county total and under18 population
county_pop <- merged_map_district_df %>%
  st_drop_geometry() %>%
  group_by(
    county
  ) %>% 
  summarize(
    pop_total = sum(pop_total),
    pop_under18 = sum(pop_under18)
  ) %>%
  na.omit() 



# Merge weighted mmr county with map and pop and plot
anti_join(tx_counties, county_pop, by = "county")


merged_map_county_df <- left_join(
  tx_counties,
  mmr_county_df_clean, 
  by = "county"
) %>%
  left_join(
    county_pop,
    by = "county"
  )


saveRDS(merged_map_county_df, "ProcessedData/merged_map_county_df.rds")



p_county_mmr_weighted_map <- ggplot(merged_map_county_df) +
  geom_sf(
    aes(fill = mmr),
    color = "gray30",
    size = 0.1
  ) +
  scale_fill_gradientn(
    colors = c("#d73027", "#fee08b", "#1a9850"),  # soft red → light yellow → muted green
    limits = c(0.8, 1.0),
    na.value = "lightgray",
    labels = scales::percent_format(accuracy = 1)
  ) +
  theme_minimal() +
  labs(
    title = "Texas Counties — MMR Coverage (2024)",
    fill = "MMR"
  )

# save to file
ggsave("Figures/county_mmr_weighted_map.png",
       plot = p_county_mmr_weighted_map,
       width = 10,
       height = 8, 
       dpi = 300
)






#Texas Population Flow-----

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












