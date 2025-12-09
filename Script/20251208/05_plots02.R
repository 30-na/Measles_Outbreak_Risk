library(cowplot)
library(ggplot2)
library(dplyr)
library(grid)
library(sf)

# Load datasets ----
infection_county   <- readRDS("ProcessedData/map_county_infection_proportion.rds")
infection_district <- readRDS("ProcessedData/map_district_infection_proportion.rds")

tx_counties <- counties(
  state = "TX",
  year  = 2024
) %>%
  mutate(county = toupper(NAME)) %>%
  select(county)


# Functions -----

compute_indirect_risk_county_new <- function(trans_mat, county_name = "GAINES", threshold = 0.5) {
  county_upper <- toupper(county_name)
  county_names <- rownames(trans_mat)
  
  if (!(county_upper %in% county_names))
    stop("County not found in matrix")
  
  # Direct transmission vector P_j→i
  Pji_vec <- trans_mat[county_upper, ]
  
  # Identify candidate first-generation counties k
  Pjk_vec  <- trans_mat[county_upper, ]          # Same row: j → k
  k_indices <- which(Pjk_vec >= threshold & !is.na(Pjk_vec))
  K         <- names(k_indices)
  
  # Output vector
  SG_vec <- rep(NA_real_, length(Pji_vec))
  names(SG_vec) <- names(Pji_vec)
  
  for (i in names(Pji_vec)) {
    
    Pji <- Pji_vec[i]  # direct probability j → i
    
    # If no k, then SG probability = Pji
    if (length(K) == 0) {
      SG_vec[i] <- Pji
      next
    }
    
    # Part 1: (1 - Pji)
    term1 <- (1 - Pji)
    
    # Part 2: product over k of (1 - Pjk * Pki)
    prod_term <- 1
    for (k in K) {
      Pjk <- Pjk_vec[k]  # j → k
      Pki <- trans_mat[k, i]  # k → i
      if (is.na(Pki)) next
      prod_term <- prod_term * (1 - Pjk * Pki)
    }
    
    # Final formula:
    # P_SG_ij = 1 - term1 * prod_term
    SG_vec[i] <- 1 - term1 * prod_term
  }
  
  # Remove self-transmission
  SG_vec[county_upper] <- NA
  
  return(SG_vec)
}


figure_MMR <- function() {
  
  # Plot A: County MMR
  p1 <- ggplot(infection_county) +
    geom_sf(aes(fill = MMR), color = "gray30", size = 0.1) +
    scale_fill_gradientn(
      colors = c("#d73027", "#fee08b", "#1a9850"),
      limits = c(0.6, 1.0),
      na.value = "lightgray",
      labels = scales::percent_format(accuracy = 1)
    ) +
    theme_minimal() +
    labs(
      title = "A: Weighted MMR Coverage by Texas County",
      fill  = "MMR"
    ) +
    theme_void() +
    theme(
      plot.title      = element_text(hjust = 0.5, size = 13),
      legend.position = "none"
    )
  
  # Plot B: District MMR
  p2 <- ggplot(infection_district) +
    geom_sf(aes(fill = MMR), color = "gray50", size = 0.1) +
    geom_sf(
      data  = tx_counties,
      fill  = NA,
      color = "gray20",
      size  = 0.4
    ) +
    scale_fill_gradientn(
      colors = c("#d73027", "#fee08b", "#1a9850"),
      limits = c(0.6, 1.0),
      na.value = "lightgray",
      labels   = scales::percent_format(accuracy = 1)
    ) +
    theme_minimal() +
    labs(
      title = "B: Weighted MMR Coverage by Texas School Districts",
      fill  = "MMR"
    ) +
    theme_void() +
    theme(
      plot.title = element_text(hjust = 0.5, size = 13)
    )
  
  legend   <- get_legend(p2)
  p2_clean <- p2 + theme(legend.position = "none")
  
  row_plot   <- plot_grid(p1, p2_clean, nrow = 1, rel_widths = c(1, 1))
  final_plot <- plot_grid(row_plot, legend, nrow = 1, rel_widths = c(2, 0.4))
  
  file_name <- "Figures/figure_B01.png"
  ggsave(file_name, final_plot, width = 12, height = 6, dpi = 400)
}

figure_MMR()


figure_outbreaks_1st_2nd_county <- function(
    method    = 7,
    strategy  = 0,
    county    = "Gaines",
    threshold = 0.5,
    map_data  = infection_county,
    out_dir   = "Figures/"
) {
  
  dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)
  
  # Define breaks and labels
  breaks     <- c(0, 0.2, 0.4, 0.6, 0.8, 1)
  bin_labels <- c("0–0.2", "0.2–0.4", "0.4–0.6", "0.6–0.8", "0.8–1.0")
  color_palette <- rev(c("#d73027", "#fc8d59", "#fee08b", "#91cf60", "#1a9850"))
  
  mat <- readRDS(paste0("ProcessedData/county_pij_M", method, "_S", strategy, ".rds"))
  
  county_upper <- toupper(county)
  county_names <- toupper(map_data$County)
  
  direct_vec <- mat[county_upper, match(county_names, colnames(mat))]
  
  indirect_vec <- compute_indirect_risk_county_new(
    trans_mat   = mat,
    county_name = county_upper,
    threshold   = threshold
  )
  
  highlighted_geom <- map_data %>% filter(toupper(County) == county_upper)
  
  # Cut values into bins
  map_data <- map_data %>%
    mutate(
      direct_bin = cut(
        direct_vec[match(toupper(County), names(direct_vec))],
        breaks = breaks,
        labels = bin_labels,
        include.lowest = TRUE,
        right = FALSE
      ),
      indirect_bin = cut(
        indirect_vec[match(toupper(County), names(indirect_vec))],
        breaks = breaks,
        labels = bin_labels,
        include.lowest = TRUE,
        right = FALSE
      )
    )
  
  # Plot A: Direct
  p1 <- ggplot(map_data) +
    geom_sf(aes(fill = direct_bin), color = "gray40", size = 0.1) +
    geom_sf(data = highlighted_geom, fill = "blue", color = "black", size = 0.3) +
    scale_fill_manual(values = color_palette, drop = FALSE, name = "Outbreak Probability") +
    labs(title = "A: Probability of First-generation Outbreak") +
    theme_void() +
    theme(
      plot.title      = element_text(hjust = 0.5, size = 13),
      legend.position = "none"
    )
  
  # Plot B: Indirect
  p2 <- ggplot(map_data) +
    geom_sf(aes(fill = indirect_bin), color = "gray40", size = 0.1) +
    geom_sf(data = highlighted_geom, fill = "blue", color = "black", size = 0.3) +
    scale_fill_manual(values = color_palette, drop = FALSE, name = "Outbreak Probability") +
    labs(title = "B: Probability of First & Second-generation Outbreak") +
    theme_void() +
    theme(
      panel.grid  = element_blank(),
      axis.text   = element_blank(),
      axis.ticks  = element_blank(),
      axis.title  = element_blank()
    )
  
  legend   <- get_legend(p2)
  p2_clean <- p2 + theme(legend.position = "none")
  
  row_plot   <- plot_grid(p1, p2_clean, nrow = 1, rel_widths = c(1, 1))
  final_plot <- plot_grid(row_plot, legend, nrow = 1, rel_widths = c(2, 0.4))
  
  file_name <- paste0(out_dir, "figure_B02.png")
  ggsave(file_name, final_plot, width = 12, height = 6, dpi = 400)
}

figure_outbreaks_1st_2nd_county(method = 7, strategy = 0, county = "Gaines")
