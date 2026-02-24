# Measles Outbreak Risk Modeling in Texas School Districts

This repository contains a framework for modeling the risk of measles outbreaks across **K–12 public school districts in Texas**. The approach uses district-level school vaccination coverage and enrollment data to estimate outbreak probability, attack rates, and expected outbreak sizes, and evaluates the impact of several vaccination improvement strategies.

## Objectives

- Estimate the risk of measles outbreaks in Texas K–12 public school districts.
- Evaluate how school-based vaccination improvements reduce outbreak risk.
- Visualize outbreak risk using district-level geographic plots.

## Data Sources

#### Vaccination Coverage
- **Source**: Texas Department of State Health Services  
- **Link**: https://www.dshs.texas.gov/immunizations/data/school/coverage  
- **Files used**:
  - Kindergarten MMR vaccination coverage by district (2019–2025)
  - Seventh-grade MMR vaccination coverage by district (2019–2025)

#### School Enrollment
- **Source**: Texas Education Agency  
- **Link**: https://tea.texas.gov/texas-schools/accountability/academic-accountability/performance-reporting/texas-academic-performance-reports  
- **File used**: `StudPgmStateDistrict25state.csv`  
- Contains total K–12 public school enrollment by district.

#### Population Age Structure
- **Source**: Texas Population Pyramids (2020)
- **Link**: https://demographics.texas.gov/Interactive/2021/CBEstimates
- **File used**: `Texas.csv`
- Used to derive grade-level population weights (ages 5–17) for estimating district-level K–12 MMR coverage.

#### Geometries Data
- **School District Geometries**: TIGER/Line shapefiles via the `tigris` R package
- **County Geometries**: TIGER/Line county shapefiles (used for map overlays)

## Vaccination Strategies

Four vaccination scenarios are analyzed:

- **Scenario 0 (Baseline)**: Uses the observed district-level K–12 MMR vaccination coverage.
- **Scenario 1 (Continued Decline)**: Recent kindergarten coverage trend is extrapolated forward 5 years if declining.
- **Scenario 2 (Constant Uptake)**: Kindergarten MMR vaccination coverage remains constant at the 2024–2025 level over the next five years.
- **Scenario 3 (Trend Reversal)**: Any decline in kindergarten MMR coverage between 2020 and 2024 is reversed over the next five years.
- **Scenario 4 (Kindergarten ≥95%**: Kindergarten MMR vaccination coverage is increased to at least 95% over the next five years.
- **Scenario 5 (District ≥95%)**: All districts achieve and maintain at least 95% MMR vaccination coverage across K–12.

These strategies are implemented in `02_calculate_infection_proportion_and_plots.r`.

## Modeling Process

- **Susceptible Population**: Calculated from district level K–12 MMR coverage and vaccine efficacy (97%).
- **Outbreak Probability**: Probability that a measles outbreak takes off in each district.
- **Internal Infection Probability**: Solves an implicit equation for the attack rate.
- **Expected Outbreak Size**: Expected number of infected students in each district.
- **Transmission Assumption**: Independent district level outbreaks (no inter district transmission).

Analyses are performed for R0 = 11.7 (baseline) and R0 = 13 (sensitivity analysis).

## Main Scripts

- `01_read_clean_merge_data.r` :
  Reads and cleans multi year vaccination data, reconstructs grade specific K–12 coverage, applies population based weighting, merges enrollment and district geometry data, and saves processed datasets.

- `02_calculate_infection_proportion_and_plots.r` :
  Computes outbreak probability, attack rate, and expected outbreak size for all vaccination scenarios and generates district-level maps and figures.

## Outputs

- **ProcessedData/**:
  - `merged_map_district_df.rds` : District-level vaccination, enrollment, and geometry data
  - `map_district_infection_proportion.rds` : Outbreak probability, attack rate, and expected outbreak size

- **Figures/**:
  - District-level maps of MMR vaccination coverage
  - Expected attack rate maps
  - Expected outbreak size maps
  - Maps showing changes in attack rate under vaccination scenarios
 
## Visual Abstract
[https://tidbitapp.io/tidbits/risk-of-measles-outbreak-in-texas-school-districts-a-modeling-study?utm_campaign=tidbitlinkshare&utm_source=IO](PASTE_YOUR_LINK_HERE)

<img width="5060" height="2560" alt="Risk of Measles Outbreak in Texas School Districts_ A Modeling Study (1)" src="https://github.com/user-attachments/assets/a762c276-cf97-4027-b29a-4797a88d582e" />

