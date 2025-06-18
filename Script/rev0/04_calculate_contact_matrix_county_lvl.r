library(sf)
library(ggplot2)
library(reshape2)
library(lubridate)

# Load county-level data
map_data <- readRDS("ProcessedData/map_data_county.rds")
dist_matrix <- readRDS("ProcessedData/distance_matrix_haversine_county.rds")
texas_flows <- read_csv("ProcessedData/texas_county_flows_2019.csv", 
                            col_types = cols(.default = "c"))

pop <- map_data$total  # county population
names(pop) <- map_data$County


# Method 1: Model A
contact_method1 <- function(pop, dist_matrix) {
  A <- 75.94
  B <- 278e-9
  C <- 1.85e4
  D <- 3.43e8
  alpha <- 1.80
  gamma <- 1.16
  
  county_names <- rownames(dist_matrix)
  mat <- matrix(NA, length(county_names), length(county_names),
                dimnames = list(county_names, county_names))
  
  for (i in county_names) {
    for (j in county_names) {
      dij <- dist_matrix[i, j]
      if (dij == 0) next  # skip self-distances or undefined cases
      mi <- pop[i]
      mj <- pop[j]
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
  
  county_names <- rownames(dist_matrix)
  mat <- matrix(NA, length(county_names), length(county_names),
                dimnames = list(county_names, county_names))
  
  for (i in county_names) {
    for (j in county_names) {
      dij <- dist_matrix[i, j]
      if (dij == 0) next
      mi <- pop[i]
      mj <- pop[j]
      term <- 1 + (B * ((mi + C) * (mj + D))^beta / dij)
      Tij <- exp(A * (term)^xi)
      mat[i, j] <- Tij
    }
  }
  return(mat)
}


# Method 3: Gravity Model
contact_method3 <- function(pop, dist_matrix) {
  county_names <- rownames(dist_matrix)
  mat <- matrix(NA, length(county_names), length(county_names),
                dimnames = list(county_names, county_names))
  
  for (i in county_names) {
    for (j in county_names) {
      dij <- dist_matrix[i, j]
      if (dij == 0 || is.na(dij)) next  # avoid divide by zero
      mi <- pop[i]
      mj <- pop[j]
      Tij <- (mi * mj) / (dij^2)
      mat[i, j] <- Tij
    }
  }
  return(mat)
}



# Method 4: Distance-Sensitive Interaction Model
contact_method4 <- function(pop, dist_matrix) {
  county_names <- rownames(dist_matrix)
  mat <- matrix(NA, length(county_names), length(county_names),
                dimnames = list(county_names, county_names))
  
  for (i in county_names) {
    for (j in county_names) {
      dij <- dist_matrix[i, j]
      if (dij == 0 || is.na(dij)) next
      mi <- pop[i]
      mj <- pop[j]
      if (dij < 300) {
        Tij <- (mi^0.46) * (mj^0.64) * exp(-dij / 82)
      } else {
        Tij <- (mi^0.35) * (mj^0.37) * exp(-dij)
      }
      mat[i, j] <- Tij
    }
  }
  return(mat)
}




# Method 5: Simplified Distance-Sensitive Model
contact_method5 <- function(pop, dist_matrix) {
  county_names <- rownames(dist_matrix)
  mat <- matrix(NA, length(county_names), length(county_names),
                dimnames = list(county_names, county_names))
  
  for (i in county_names) {
    for (j in county_names) {
      dij <- dist_matrix[i, j]
      if (dij == 0 || is.na(dij)) next
      mi <- pop[i]
      mj <- pop[j]
      if (dij < 300) {
        Tij <- mi * mj * exp(-dij / 82)
      } else {
        Tij <- mi * mj * exp(-dij)
      }
      mat[i, j] <- Tij
    }
  }
  return(mat)
}


# Method 6: Piecewise Power-Law Decay Model
contact_method6 <- function(pop, dist_matrix) {
  county_names <- rownames(dist_matrix)
  mat <- matrix(NA, length(county_names), length(county_names),
                dimnames = list(county_names, county_names))
  
  for (i in county_names) {
    for (j in county_names) {
      dij <- dist_matrix[i, j]
      if (dij == 0 || is.na(dij)) next
      mi <- pop[i]
      mj <- pop[j]
      if (dij < 390) {
        Tij <- (mi^0.48) * (mj^0.58) / (dij^4.35)
      } else {
        Tij <- (mi^0.40) * (mj^0.33) / (dij^0.49)
      }
      mat[i, j] <- Tij
    }
  }
  return(mat)
}



# Method7
flows_avg <- texas_flows %>%
  mutate(
    pop_flows = as.numeric(pop_flows),
    week_start = mdy(str_sub(date_range, 1, 8))
    ) %>%
  # filter(month(week_start) >= 1 & month(week_start) <= 6) %>%
  group_by(
    geoid_o,
    geoid_d, 
    county_o, 
    county_d
    ) %>%
  summarize(
    flow = sum(pop_flows, na.rm = TRUE),
    .groups = "drop"
    )




contact_method7 <- function(flows_avg) {
  county_names <- sort(unique(c(flows_avg$county_o, flows_avg$county_d)))
  mat <- matrix(0, length(county_names), length(county_names),
                dimnames = list(county_names, county_names))
  
  for (k in seq_len(nrow(flows_avg))) {
    origin <- flows_avg$county_o[k]
    dest <- flows_avg$county_d[k]
    flow <- flows_avg$flow[k]
    
    if (origin %in% county_names && dest %in% county_names) {
      if (origin == dest) {
        mat[origin, dest] <- NA  # self-flows as NA
      } else {
        mat[origin, dest] <- flow
      }
    }
  }
  
  return(mat)
}




# Compute and save
contact_Haversine01 <- contact_method1(pop, dist_matrix)
saveRDS(contact_Haversine01, "ProcessedData/contact_matrix1_county.rds")

contact_Haversine02 <- contact_method2(pop, dist_matrix)
saveRDS(contact_Haversine02, "ProcessedData/contact_matrix2_county.rds")

contact_Haversine03 <- contact_method3(pop, dist_matrix)
saveRDS(contact_Haversine03, "ProcessedData/contact_matrix3_county.rds")

contact_Haversine04 <- contact_method4(pop, dist_matrix)
saveRDS(contact_Haversine04, "ProcessedData/contact_matrix4_county.rds")

contact_Haversine05 <- contact_method5(pop, dist_matrix)
saveRDS(contact_Haversine05, "ProcessedData/contact_matrix5_county.rds")

contact_Haversine06 <- contact_method6(pop, dist_matrix)
saveRDS(contact_Haversine06, "ProcessedData/contact_matrix6_county.rds")

contact_Haversine07 <- contact_method7(flows_avg)
saveRDS(contact_Haversine07, "ProcessedData/contact_matrix7_county.rds")





############# TEST 
# Check that Gaines is in the row names
# which(rownames(contact_Haversine03) == "GAINES")
# 
# # Extract the row for Gaines (flows from Gaines to all others)
# gaines_row <- contact_Haversine03["GAINES", ]
# 
# # Remove self-flow (usually NA or 0 on diagonal)
# gaines_row <- gaines_row[names(gaines_row) != "GAINES"]
# 
# # Plot histogram
# hist(gaines_row,
#      main = "Gravity from Gaines County (Method 3)",
#      xlab = "Gravity",
#      col = "steelblue",
#      border = "gray40",
#      breaks=50)
# 
# 
# # Check that Gaines is in the row names
# which(rownames(contact_Haversine03) == "GAINES")
# 
# # Extract the row for Gaines (flows from Gaines to all others)
# gaines_row <- contact_Haversine03["GAINES", ]
# 
# # Remove self-flow (usually NA or 0 on diagonal)
# gaines_row <- gaines_row[names(gaines_row) != "GAINES"]
# 
# # Plot histogram
# hist(gaines_row,
#      main = "Outgoing Flows from Gaines County (Method 7)",
#      xlab = "Estimated Flow (max weekly pop)",
#      col = "steelblue",
#      border = "gray40",
#      breaks=50)
# 
# 
# 
# 
# # SUM of flow for gain
# flows_gaines <- texas_flows %>%
#   mutate(
#     pop_flows = as.numeric(pop_flows),
#     start_date = mdy(sub(" -.*", "", date_range))  # extract and parse start date
#   ) %>%
#   #filter(start_date >= as.Date("2019-01-01") & start_date <= as.Date("2019-12-30")) %>%
#   group_by(
#     geoid_o, 
#     geoid_d, 
#     county_o, 
#     county_d
#   ) %>%
#   summarize(
#     flow = sum(pop_flows, na.rm = TRUE),
#     .groups = "drop"
#   ) %>%
#   filter(county_o == "GAINES")
# 
# 
# 
# ######### TEST 
# # Load county boundaries
# tx_counties <- counties(state = "TX", cb = TRUE, year = 2020) %>%
#   st_transform(4326) %>%
#   mutate(
#     NAME = toupper(NAME)
#   )
# 
# 
# # Merge with flow data
# gaines_flows_map <- tx_counties %>%
#   left_join(
#     gaines_flows, 
#     by = c("NAME" = "county_d"))
# mutate(
#   flow_label = ifelse(!is.na(avg_flow), as.character(round(avg_flow)), NA),
#   color_fill = case_when(
#     county_d == "GAINES" ~ "blue",
#     !is.na(avg_flow) ~ "white",
#     TRUE ~ "gray80"
#   )
# )
# 
# map_counties <- unique(tx_map$county_d)
# 
# # Assuming flow_gaines is your filtered flow data from Gaines County
# flow_counties <- unique(avg_flows$county_d)
# flow_tx <- unique(texas_flows$county_d)
# # Counties in map but not in flow data
# missing_in_flow <- setdiff(map_counties, flow_counties)
# 
# # Counties in flow data but not in map
# missing_in_map <- setdiff(flow_counties, map_counties)
# 
# 
# # Compute centroids for label placement
# centroids <- st_centroid(gaines_flows$geometry)
# centroid_coords <- st_coordinates(centroids)
# gaines_flows$x <- centroid_coords[, 1]
# gaines_flows$y <- centroid_coords[, 2]
# gaines_flows <- st_as_sf(gaines_flows)
# # Plot
# g <- ggplot(gaines_flows) +
#   geom_sf(aes(fill = color_fill), color = "gray30", size = 0.2) +
#   geom_text(aes(x = x, y = y, label = flow_label), size = 1.5, color = "black") +
#   scale_fill_identity() +
#   theme_void() +
#   labs(title = "Flows from Gaines County to Texas Counties (2019)")
# 
# ggsave("Figures/flow_county_gaines_sum1-12.png", 
#        g,
#        width = 12,
#        height = 8,
#        dpi = 400)
# 
# 
