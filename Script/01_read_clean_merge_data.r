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
library(tidyverse)



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
  )%>%
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
    ),
    
    # Scenario 01 vaccine uptake at kindergarten remains constant over the next 5 years
    mmr1 = rowSums(
      cbind(
        K  * grade_weights["K"],
        K * grade_weights["G1"],
        K * grade_weights["G2"],
        K * grade_weights["G3"],
        K * grade_weights["G4"],
        K * grade_weights["G5"],
        G1 * grade_weights["G6"],
        G2 * grade_weights["G7"],
        G3 * grade_weights["G8"],
        G4 * grade_weights["G9"],
        G5 * grade_weights["G10"],
        G6 * grade_weights["G11"],
        G7 * grade_weights["G12"]
      )
    ),
    # Scenario 02 vaccine uptake trend is reversed over the next 5 years
    k_slope = coef(lm(c(G4, G3, G2, G1, K) ~ c(1:5) ))[2],
    mmr2 = case_when(
      k_slope < 0 ~ rowSums(
        cbind(
          G4 * grade_weights["K"],
          G3 * grade_weights["G1"],
          G2 * grade_weights["G2"],
          G1 * grade_weights["G3"],
          K  * grade_weights["G4"],
          K  * grade_weights["G5"],
          G1 * grade_weights["G6"],
          G2 * grade_weights["G7"],
          G3 * grade_weights["G8"],
          G4 * grade_weights["G9"],
          G5 * grade_weights["G10"],
          G6 * grade_weights["G11"],
          G7 * grade_weights["G12"]
        )
      ),
      TRUE ~ mmr1
    ),
    
    # Scenario 3 kindergarten coverage to be 95% each year for the next 5 years. 
    mmr3 = rowSums(
      cbind(
        pmax(K, 0.95) * grade_weights["K"],
        pmax(K, 0.95) * grade_weights["G1"],
        pmax(K, 0.95) * grade_weights["G2"],
        pmax(K, 0.95) * grade_weights["G3"],
        pmax(K, 0.95) * grade_weights["G4"],
        K * grade_weights["G5"],
        G1 * grade_weights["G6"],
        G2 * grade_weights["G7"],
        G3 * grade_weights["G8"],
        G4 * grade_weights["G9"],
        G5 * grade_weights["G10"],
        G6 * grade_weights["G11"],
        G7 * grade_weights["G12"]
      )
    ),
    
    # scenario 04 all districts with total coverage below 95% are raised to 95% coverage
    mmr4 = ifelse(mmr < .95, .95, mmr)
  )%>%
  select(
    district,
    county,
    mmr,
    mmr1,
    mmr2,
    mmr3,
    mmr4
  )



# Load District Population

district_pop_path <- "RawData/StudPgmStateDistrict25state.csv"
pop_raw <- read_csv(district_pop_path, skip = 5)


pop_district_df <- pop_raw %>%
select(
  district = "DISTRICT NAME",
  enrollment = "ALL ENROLLMENT",
  district_number = "DISTRICT NUMBER"
) %>%
  mutate(
    enrollment = as.numeric(enrollment),
    enrollment = if_else(district == "SAN VICENTE ISD", 10, enrollment)
  ) %>%
  drop_na()

print(pop_raw[
  duplicated(pop_raw$`DISTRICT NAME`) |
    duplicated(pop_raw$`DISTRICT NAME`, fromLast = TRUE),
], n=30)



pop_district_df_clean <- pop_district_df %>%
  mutate(
    district = str_replace(district, "CISD", "CONSOLIDATED ISD"),
    district = str_replace(district, "CSD", "COMMON SCHOOL DISTRICT")
  ) %>%
  mutate(
    district = case_when(
      district == "GOLD BURG ISD" ~ "GOLD-BURG ISD",
      district == "STAFFORD MSD" ~ "STAFFORD MUNICIPAL SCHOOL DISTRICT",
      district == "LAPOYNOR ISD" ~ "LAPOYNER ISD",
      district == "HAMLIN COLLEGIATE ISD" ~ "HAMLIN ISD",
      
      district == "EAGLE MT-SAGINAW ISD" ~ "EAGLE MOUNTAIN-SAGINAW ISD",
      
      district == "FT DAVIS ISD" ~ "FORT DAVIS ISD",
      district == "FT HANCOCK ISD" ~ "FORT HANCOCK ISD",
      district == "FT SAM HOUSTON ISD" ~ "FORT SAM HOUSTON ISD",
      
      district == "GOLDTHWAITE ISD" ~ "GOLDTHWAITE CONSOLIDATED ISD",
      district == "KENEDY COUNTY WIDE COMMON SCHOOL DISTRICT" ~ "KENEDY COUNTY-WIDE COMMON SCHOOL DISTRICT",
      district == "SCHERTZ-CIBOLO-U CITY ISD" ~ "SCHERTZ-CIBOLO-UNIVERSAL CITY ISD",
      district == "HAMLIN COLLEGIATE ISD" ~ "HAMLIN ISD",
      
      TRUE ~ district
    )
  ) %>%
  mutate(
    district = case_when(
      # BIG SANDY ISD
      district == "BIG SANDY ISD" & district_number == "230901" ~ "BIG SANDY ISD (UPSHUR)",
      district == "BIG SANDY ISD" & district_number == "187901" ~ "BIG SANDY ISD (POLK)",
      
      # CENTERVILLE ISD
      district == "CENTERVILLE ISD" & district_number == "145902" ~ "CENTERVILLE ISD (LEON)",
      district == "CENTERVILLE ISD" & district_number == "228904" ~ "CENTERVILLE ISD (TRINITY)",
      
      # CHAPEL HILL ISD
      district == "CHAPEL HILL ISD" & district_number == "212909" ~ "CHAPEL HILL ISD (SMITH)",
      district == "CHAPEL HILL ISD" & district_number == "225906" ~ "CHAPEL HILL ISD (TITUS)",
      
      # DAWSON ISD
      district == "DAWSON ISD" & district_number == "175904" ~ "DAWSON ISD (NAVARRO)",
      district == "DAWSON ISD" & district_number == "058902" ~ "DAWSON ISD (DAWSON)",
      
      # EDGEWOOD ISD
      district == "EDGEWOOD ISD" & district_number == "234903" ~ "EDGEWOOD ISD (VAN ZANDT)",
      district == "EDGEWOOD ISD" & district_number == "015905" ~ "EDGEWOOD ISD (BEXAR)",
      
      # HIGHLAND PARK ISD
      district == "HIGHLAND PARK ISD" & district_number == "057911" ~ "HIGHLAND PARK ISD (DALLAS)",
      district == "HIGHLAND PARK ISD" & district_number == "188903" ~ "HIGHLAND PARK ISD (POTTER)",
      
      # HUBBARD ISD
      district == "HUBBARD ISD" & district_number == "019913" ~ "HUBBARD ISD (BOWIE)",
      district == "HUBBARD ISD" & district_number == "109905" ~ "HUBBARD ISD (HILL)",
      
      # MIDWAY ISD
      district == "MIDWAY ISD" & district_number == "039905" ~ "MIDWAY ISD (CLAY)",
      district == "MIDWAY ISD" & district_number == "161903" ~ "MIDWAY ISD (MCLENNAN)",
      
      # NORTHSIDE ISD
      district == "NORTHSIDE ISD" & district_number == "244905" ~ "NORTHSIDE ISD (WILBARGER)",
      district == "NORTHSIDE ISD" & district_number == "015915" ~ "NORTHSIDE ISD (BEXAR)",
      
      # VALLEY VIEW ISD
      district == "VALLEY VIEW ISD" & district_number == "108916" ~ "VALLEY VIEW ISD (HIDALGO)",
      district == "VALLEY VIEW ISD" & district_number == "049903" ~ "VALLEY VIEW ISD (COOKE)",
      
      # WYLIE ISD
      district == "WYLIE ISD" & district_number == "043914" ~ "WYLIE ISD (COLLIN)",
      district == "WYLIE ISD" & district_number == "221912" ~ "WYLIE ISD (TAYLOR)",
      
      TRUE ~ district
    )
  ) %>%
  select(
    district,
    enrollment
  )


missing_districts <- tx_school_district_map_clean %>%
  distinct(district) %>%
  anti_join(
    pop_district_df_clean %>% distinct(district),
    by = "district"
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


saveRDS(merged_map_district_df, "ProcessedData/merged_map_district_df.rds")



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
    limits = c(0.6, 1.0),
    na.value = "lightgray",
    labels = scales::percent_format(accuracy = 1)
  ) +
  labs(
    title = "School district-level MMR vaccination coverage",
    fill = "MMR"
  ) +
    theme(
      panel.background = element_rect(fill = "white", color = NA),
      plot.background  = element_rect(fill = "white", color = NA),
      panel.grid.major = element_blank(),
      panel.grid.minor = element_blank(),
      plot.title = element_text(size = 20, face = "bold", hjust = 0)
    ) +
  coord_sf(datum = NA)

#save to file
ggsave("Figures/district_mmr_weighted_map.png",
       plot = p_district_mmr_weighted_map,
       width = 10,
       height = 8,
       dpi = 300
)
