library(sf)
library(ggplot2)
library(reshape2)
library(lubridate)
library(dplyr)
library(tidyr)
library(readr)
library(tigris)

#Load datasets ----
infection_county <- readRDS("ProcessedData/map_county_infection_proportion.rds")
infection_district <- readRDS("ProcessedData/map_district_infection_proportion.rds")

# Population Flow County
texas_flows <- read_csv("ProcessedData/texas_county_flows_2019.csv", 
                        col_types = cols(.default = "c"))



# Functions ----
contact_method7 <- function(flows_avg) {
  county_names <- sort(unique(c(flows_avg$county_o, flows_avg$county_d)))
  mat <- matrix(0, 
                length(county_names), 
                length(county_names),
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

# Texas Population Flow County ----
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

C7 <- contact_method7(flows_avg)


#------------------ Transmission/ Outbreak Probability COUNTY ------------------

compute_transmission_matrix_county <- function(Cij, map_data, strategy = 0, q = .9) {
  suffix <- ifelse(strategy == 0, "", as.character(strategy))
  
  county_names <- map_data$County
  psi <- map_data[[paste0("susceptible_prop", suffix)]]
  pmo <- map_data[[paste0("outbreak_prob", suffix)]]
  Pj <- map_data[[paste0("adult_infection_prob", suffix)]]
  
  n <- length(county_names)
  transmission_mat <- matrix(NA, n, n, dimnames = list(county_names, county_names))
  
  for (i in seq_len(n)) {
    for (j in seq_len(n)) {
      name_i <- county_names[i]
      name_j <- county_names[j]
      
      base <- q * Pj[which(county_names == name_j)] *
        psi[which(county_names == name_i)] *
        pmo[which(county_names == name_i)]
      
      cij_val <- Cij[name_i, name_j]
      if (is.na(cij_val)) {
        transmission_mat[name_i, name_j] <- NA
      } else {
        transmission_mat[name_i, name_j] <- 1 - (1 - base)^cij_val
      }
      
    }
  }
  
  return(transmission_mat)
}


# Method 7 (Mobility Flows)
county_pij_M7_S0 <- compute_transmission_matrix_county(C7, infection_county)
county_pij_M7_S1 <- compute_transmission_matrix_county(C7, infection_county, strategy = 1)
county_pij_M7_S2 <- compute_transmission_matrix_county(C7, infection_county, strategy = 2)
county_pij_M7_S3 <- compute_transmission_matrix_county(C7, infection_county, strategy = 3)


saveRDS(county_pij_M7_S0, "ProcessedData/county_pij_M7_S0.rds")
saveRDS(county_pij_M7_S1, "ProcessedData/county_pij_M7_S1.rds")
saveRDS(county_pij_M7_S2, "ProcessedData/county_pij_M7_S2.rds")
saveRDS(county_pij_M7_S3, "ProcessedData/county_pij_M7_S3.rds")



#-------------------Transmission/ Outbreak Probability District------------------


compute_transmission_matrix_district <- function(Cij, map_data, origin_scaler=NULL, dest_scaler=NULL, strategy = 0, q = 0.9) {
  
  district_pop_ratio <- map_data %>%
    st_drop_geometry() %>%                
    group_by(
      County
    ) %>%
    mutate(
      county_pop_total = sum(pop_total, na.rm = TRUE),
      scaler = pop_total / county_pop_total
    ) %>%
    ungroup() %>%
    filter(
      !is.na(County)
    )
  
  suffix <- ifelse(strategy == 0, "", as.character(strategy))
  
  district_names <- district_pop_ratio$district
  county_of_district <- district_pop_ratio$County
  
  psi <- district_pop_ratio[[paste0("susceptible_prop", suffix)]]
  pmo <- district_pop_ratio[[paste0("outbreak_prob", suffix)]]
  Pj <- district_pop_ratio[[paste0("adult_infection_prob", suffix)]]
  
  n <- length(district_names)
  scaler <- district_pop_ratio$scaler
  transmission_mat <- matrix(NA, n, n, dimnames = list(district_names, district_names))
  
  for (district_i in district_names) {
    for (district_j in district_names) {
      county_i <- county_of_district[district_names == district_i]
      county_j <- county_of_district[district_names == district_j]
      
      
      base <- q *
        Pj[district_names == district_j] *
        psi[district_names == district_i] *
        pmo[district_names == district_i]
      
      cij_base <- Cij[county_i, county_j]
      
      if (is.null(dest_scaler)) {
        scaler_i <- scaler[district_names == district_i]
      } else {
        scaler_i <- dest_scaler
      }
      
      if (is.null(origin_scaler)) {
        scaler_j <- scaler[district_names == district_j]
      } else {
        scaler_j <- origin_scaler
      }
      
      
      cij_scaled <- cij_base * scaler_i * scaler_j
      
      if (is.na(cij_base)) {
        transmission_mat[district_i, district_j] <- NA
      } else {
        transmission_mat[district_i, district_j] <- 1 - (1 - base)^cij_scaled
      }
      
    }
  }
  
  return(transmission_mat)
}


district_pij_M7_S0_CountyToDistrict <- compute_transmission_matrix_district(C7, infection_district, origin_scaler=1, dest_scaler=NULL, strategy = 0)
district_pij_M7_S1_CountyToDistrict <- compute_transmission_matrix_district(C7, infection_district, origin_scaler=1, dest_scaler=NULL, strategy = 1)
district_pij_M7_S2_CountyToDistrict <- compute_transmission_matrix_district(C7, infection_district, origin_scaler=1, dest_scaler=NULL, strategy = 2)
district_pij_M7_S3_CountyToDistrict <- compute_transmission_matrix_district(C7, infection_district, origin_scaler=1, dest_scaler=NULL, strategy = 3)

saveRDS(district_pij_M7_S0_CountyToDistrict, "ProcessedData/district_pij_M7_S0_CountyToDistrict.rds")
saveRDS(district_pij_M7_S1_CountyToDistrict, "ProcessedData/district_pij_M7_S1_CountyToDistrict.rds")
saveRDS(district_pij_M7_S2_CountyToDistrict, "ProcessedData/district_pij_M7_S2_CountyToDistrict.rds")
saveRDS(district_pij_M7_S3_CountyToDistrict, "ProcessedData/district_pij_M7_S3_CountyToDistrict.rds")



district_pij_M7_S0_CountyToCounty <- compute_transmission_matrix_district(C7, infection_district, origin_scaler=1, dest_scaler=1, strategy = 0)
district_pij_M7_S1_CountyToCounty <- compute_transmission_matrix_district(C7, infection_district, origin_scaler=1, dest_scaler=1, strategy = 1)
district_pij_M7_S2_CountyToCounty <- compute_transmission_matrix_district(C7, infection_district, origin_scaler=1, dest_scaler=1, strategy = 2)
district_pij_M7_S3_CountyToCounty <- compute_transmission_matrix_district(C7, infection_district, origin_scaler=1, dest_scaler=1, strategy = 3)

saveRDS(district_pij_M7_S0_CountyToCounty, "ProcessedData/district_pij_M7_S0_CountyToCounty.rds")
saveRDS(district_pij_M7_S1_CountyToCounty, "ProcessedData/district_pij_M7_S1_CountyToCounty.rds")
saveRDS(district_pij_M7_S2_CountyToCounty, "ProcessedData/district_pij_M7_S2_CountyToCounty.rds")
saveRDS(district_pij_M7_S3_CountyToCounty, "ProcessedData/district_pij_M7_S3_CountyToCounty.rds")




#---------------------- Get the indirect outbreak risk -----------------------
# 
# compute_indirect_risk_county <- function(trans_mat, county_name = "GAINES", threshold = .5) {
#   county_upper <- toupper(county_name)
#   county_names <- rownames(trans_mat)
#   
#   if (!(county_upper %in% county_names)) stop("County not found in matrix")
#   
#   direct_vec <- trans_mat[county_upper, ]
#   ik_vec <- trans_mat[county_upper, ]
#   k_indices <- which(ik_vec >= threshold & !is.na(ik_vec))
#   
#   combined_vec <- rep(NA_real_, length(direct_vec))
#   names(combined_vec) <- names(direct_vec)
#   
#   for (j in names(direct_vec)) {
#     pj <- direct_vec[j]
#     sum_indirect <- 0
#     valid_ks <- names(k_indices)
#     
#     if (length(valid_ks) > 0) {
#       for (k in valid_ks) {
#         pik <- ik_vec[k]
#         pkj <- trans_mat[k, j]
#         if (!is.na(pkj)) sum_indirect <- sum_indirect + (pik * pkj)
#       }
#       sum_indirect <- sum_indirect / length(valid_ks)
#       #print(length(valid_ks))
#     } else {
#       sum_indirect <- 0
#     }
#     
#     combined_vec[j] <- min(1, pj + sum_indirect)
#   }
#   
#   combined_vec[county_upper] <- NA  # Hide self-transmission
#   return(combined_vec)
# }
# 
# 
# get_indirect_risk_list_county <- function(strategy=0, county_name = "GAINES", method = 7, threshold = 0.5) {
#   # Load transmission matrix
#   file_path <- paste0("ProcessedData/county_pij_M", method, "_S", strategy, ".rds")
#   trans_mat <- readRDS(file_path)
#   
#   # Compute indirect risk
#   risk_vec <- compute_indirect_risk_county(trans_mat, county_name = county_name, threshold = threshold)
#   
#   # Return as sorted data.frame
#   out_df <- data.frame(
#     county = names(risk_vec),
#     indirect_risk = round(risk_vec, 5)
#   ) %>%
#     arrange(desc(indirect_risk))
#   
#   return(out_df)
# }
# 
# 
# 
# figure03 <- function(method = 7, counties = c("Gaines"), strategies = c(1, 2, 3),
#                      map_data = map_probability, threshold = 0.5, out_dir = "Figures/") {
#   
#   dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)
#   
#   # Discrete bins
#   breaks <- c(0, 0.2, 0.4, 0.6, 0.8, 1)
#   bin_labels <- c("0–0.2", "0.2–0.4", "0.4–0.6", "0.6–0.8", "0.8–1.0")
#   color_palette <- rev(c("#d73027", "#fc8d59", "#fee08b", "#91cf60", "#1a9850"))
#   names(color_palette) <- bin_labels
#   
#   full_labels <- c(
#     "Strategy 0: Baseline",
#     "Strategy 1: MMR < 90% → 90%", 
#     "Strategy 2: MMR < 92% → 92%", 
#     "Strategy 3: MMR +5%",
#     "Strategy 4: Gaines → 90%",
#     "Strategy 5: Gaines → 92%"
#   )
#   
#   strategy_labels <- full_labels[strategies + 1]
#   label_prefix <- LETTERS[seq_along(strategies)]
#   
#   county_names <- toupper(map_data$County)
#   mats <- lapply(strategies, function(s) readRDS(paste0("ProcessedData/county_pij_M", method, "_S", s, ".rds")))
#   
#   for (county in counties) {
#     county_upper <- toupper(county)
#     highlighted_geom <- map_data %>% filter(toupper(County) == county_upper)
#     
#     maps <- list()
#     
#     for (i in seq_along(strategies)) {
#       s <- strategies[i]
#       mat <- mats[[i]]
#       indirect_vec <- compute_indirect_risk_county(mat, county_name = county_upper, threshold = threshold)
#       
#       bin_vec <- cut(indirect_vec, breaks = breaks, labels = bin_labels, include.lowest = TRUE, right = FALSE)
#       bin_fact <- factor(bin_vec[match(toupper(map_data$County), names(indirect_vec))], levels = bin_labels)
#       map_data$bin <- bin_fact
#       
#       p <- ggplot(map_data) +
#         geom_sf(aes(fill = bin), color = "gray40", size = 0.1) +
#         geom_sf(data = highlighted_geom, fill = "blue", color = "black", size = 0.3) +
#         scale_fill_manual(
#           values = color_palette,
#           name = "Outbreak Probability",
#           drop = FALSE
#         ) +
#         labs(
#           #title = paste0(label_prefix[i], ": ", strategy_labels[i])
#           title = "Second generation outbreak probability"
#         ) +
#         theme_void() +
#         theme(
#           plot.title = element_text(hjust = 0.5, size = 13),
#           legend.position = "none"
#         )
#       
#       maps[[i]] <- p
#     }
#     
#     # legend <- get_legend(
#     #   maps[[2]] + theme(legend.position = "right", legend.title = element_text(size = 10))
#     # )
#     legend <- get_legend(
#       maps[[length(maps)]] + theme(legend.position = "right", legend.title = element_text(size = 10))
#     )
#     
#     row_plot <- plot_grid(plotlist = lapply(maps, \(p) p + theme(legend.position = "none")), nrow = 1)
#     final_plot <- plot_grid(row_plot, legend, nrow = 1, rel_widths = c(length(maps), 0.4))
#     
#     file_name <- paste0(out_dir, "figure03_Indirect_", tolower(county), "_method", method, "_S", paste(strategies, collapse = "_"), ".png")
#     ggsave(file_name, final_plot, width = 6 * length(strategies), height = 6, dpi = 400)
#   }
# }
# 
# 
# figure03_district <- function(method = 7, districts = c("LOOP ISD"), strategies = c(1, 2, 3),
#                               map_data = map_probability, threshold = 0.5, out_dir = "Figures/") {
#   
#   dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)
#   
#   # Discrete bins
#   breaks <- c(0, 0.2, 0.4, 0.6, 0.8, 1)
#   bin_labels <- c("0–0.2", "0.2–0.4", "0.4–0.6", "0.6–0.8", "0.8–1.0")
#   color_palette <- rev(c("#d73027", "#fc8d59", "#fee08b", "#91cf60", "#1a9850"))
#   names(color_palette) <- bin_labels
#   
#   full_labels <- c(
#     "Strategy 0: Baseline",
#     "Strategy 1: MMR < 90% → 90%", 
#     "Strategy 2: MMR < 92% → 92%", 
#     "Strategy 2: MMR < 95% → 95%",
#     "Strategy 4: Gaines → 90%",
#     "Strategy 5: Gaines → 92%"
#   )
#   
#   strategy_labels <- full_labels[strategies + 1]
#   label_prefix <- LETTERS[seq_along(strategies)]
#   
#   district_names <- toupper(map_data$district)
#   mats <- lapply(strategies, function(s) readRDS(paste0("ProcessedData/district_pij_M", method, "_S", s, ".rds")))
#   
#   for (district in districts) {
#     district_upper <- toupper(district)
#     highlighted_geom <- map_data %>% filter(toupper(district) == district_upper)
#     
#     maps <- list()
#     
#     for (i in seq_along(strategies)) {
#       s <- strategies[i]
#       mat <- mats[[i]]
#       indirect_vec <- compute_indirect_risk_county(mat, county_name = district_upper, threshold = threshold)
#       
#       bin_vec <- cut(indirect_vec, breaks = breaks, labels = bin_labels, include.lowest = TRUE, right = FALSE)
#       bin_fact <- factor(bin_vec[match(toupper(map_data$district), names(indirect_vec))], levels = bin_labels)
#       map_data$bin <- bin_fact
#       
#       p <- ggplot(map_data) +
#         geom_sf(
#           aes(fill = bin), 
#           color = "gray40",
#           size = 0.1
#         ) +
#         geom_sf(
#           data = infection_county,
#           fill = NA,
#           color = "gray20",
#           size = 0.4
#         )+
#         geom_sf(data = highlighted_geom, fill = "blue", color = "black", size = 0.3) +
#         scale_fill_manual(
#           values = color_palette,
#           name = "Outbreak Probability",
#           drop = FALSE
#         ) +
#         labs(
#           title = paste0(label_prefix[i], ": ", strategy_labels[i])
#           # title = "Second generation outbreak probability"
#         ) +
#         theme_void() +
#         theme(
#           plot.title = element_text(hjust = 0.5, size = 13),
#           legend.position = "none"
#         )
#       
#       maps[[i]] <- p
#     }
#     
#     # legend <- get_legend(
#     #   maps[[2]] + theme(legend.position = "right", legend.title = element_text(size = 10))
#     # )
#     legend <- get_legend(
#       maps[[1]] + theme(legend.position = "right", legend.title = element_text(size = 10))
#     )
#     
#     row_plot <- plot_grid(plotlist = lapply(maps, \(p) p + theme(legend.position = "none")), nrow = 1)
#     final_plot <- plot_grid(row_plot, legend, nrow = 1, rel_widths = c(length(maps), 0.4))
#     
#     file_name <- paste0(out_dir, "figure03_Indirect_", tolower(district), "_method", method, "_S", paste(strategies, collapse = "_"), ".png")
#     ggsave(file_name, final_plot, width = 6 * length(strategies), height = 6, dpi = 400)
#   }
# }
# 
# 
# figure03_district(method = 7, infection_district, districts =  c("LOOP ISD"), strategies = c(1,2,3),threshold = .5)
# 
# 
