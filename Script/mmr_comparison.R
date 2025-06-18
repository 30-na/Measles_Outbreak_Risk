library(readxl)
library(gdata)
library(writexl)
library(gdata)

mmr_2022 <- read.xls(mmr_path_2022, sheet = 2, stringsAsFactors = FALSE)
# File paths
mmr_path_2024 <- "RawData/2023-2024_School_Vaccination_Coverage_Levels_Kindergarten.xlsx"
mmr_path_2023 <- "RawData/22-23-School-Vaccination-Coverage-by-District-and-County-K.xlsx"
#mmr_path_2022 <- "RawData/2021-2022-School-Vaccination-Coverage-by-District-and-County-Kindergarten.xls"
mmr_path_2021 <- "RawData/2020-2021-School-Vaccination-Coverage-Levels-by-District-Private-School-and-County---Kindergarten.xlsx"
mmr_path_2020 <- "RawData/2019-2020-School-Vaccination-Coverage-Levels---Kindergarten.xlsx"

# Read sheets
sheet_2024 <- excel_sheets(mmr_path_2024)[2]
sheet_2023 <- excel_sheets(mmr_path_2023)[2]
#mmr_2022 <- read_xls(mmr_path_2022, sheet = 2)
sheet_2021 <- excel_sheets(mmr_path_2021)[2]
sheet_2020 <- excel_sheets(mmr_path_2020)[2]

# Read each year and keep only County and MMR columns
mmr_2024 <- read_excel(mmr_path_2024, sheet = sheet_2024)
mmr_2023 <- read_excel(mmr_path_2023, sheet = sheet_2023, skip=2)
#mmr_2022 <- read_excel(mmr_path_2022, sheet = sheet_2022)
mmr_2021 <- read_excel(mmr_path_2021, sheet = sheet_2021, skip=2)
mmr_2020 <- read_excel(mmr_path_2020, sheet = sheet_2020, skip=2)

# Clean and select manually
mmr_2024 <- mmr_2024[, c("County", "MMR")]
mmr_2023 <- mmr_2023[, c("County", "MMR")]
#mmr_2022 <- mmr_2022[, c("County", "MMR")]
mmr_2021 <- mmr_2021[, c("County", "MMR")]
mmr_2020 <- mmr_2020[, c("County", "MMR")]

# Rename MMR columns
names(mmr_2024)[2] <- "MMR_2024"
names(mmr_2023)[2] <- "MMR_2023"
#names(mmr_2022)[2] <- "MMR_2022"
names(mmr_2021)[2] <- "MMR_2021"
names(mmr_2020)[2] <- "MMR_2020"

# Merge all years by County
df <- merge(mmr_2024, mmr_2023, by = "County", all = TRUE)

#df <- merge(df, mmr_2022, by = "County", all = TRUE)
df <- merge(df, mmr_2021, by = "County", all = TRUE)
df <- merge(df, mmr_2020, by = "County", all = TRUE)

names(df)

df <- df %>%
  mutate(
    MMR_2024 = as.numeric(MMR_2024),
    MMR_2023 = as.numeric(MMR_2023),
    MMR_2021 = as.numeric(MMR_2021),
    MMR_2020 = as.numeric(MMR_2020)
  )


# Compute 5-year average (2020–2023)
df$Average_4yr <- rowMeans(df[, c("MMR_2020", "MMR_2021", "MMR_2023", "MMR_2024")], na.rm = TRUE)

# Compute proportional difference
df$Proportional_Diff <- (df$MMR_2024 - df$Average_4yr) / df$Average_4yr

# Show result
print(df[, c("County", "MMR_2024", "Average_4yr", "Proportional_Diff")])



write_xlsx(
  df[, c("County", "MMR_2024", "Average_4yr", "Proportional_Diff")],
  path = "processedData/MMR_Comparison_by_County_4years_average.xlsx"
)
