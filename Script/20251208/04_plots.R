library(cowplot)
library(ggplot2)
library(dplyr)
library(grid)
library(ggplot2)
library(cowplot)
library(dplyr)
library(grid)
library(sf)


#Load datasets ----
infection_county <- readRDS("ProcessedData/map_county_infection_proportion.rds")
infection_district <- readRDS("ProcessedData/map_district_infection_proportion.rds")

# Functions -----

compute_indirect_risk_county <- function(trans_mat, county_name = "GAINES", threshold = .5) {
  county_upper <- toupper(county_name)
  county_names <- rownames(trans_mat)
  
  if (!(county_upper %in% county_names)) stop("County not found in matrix")
  
  direct_vec <- trans_mat[county_upper, ]
  ik_vec <- trans_mat[county_upper, ]
  k_indices <- which(ik_vec >= threshold & !is.na(ik_vec))
  
  combined_vec <- rep(NA_real_, length(direct_vec))
  names(combined_vec) <- names(direct_vec)
  
  for (j in names(direct_vec)) {
    pj <- direct_vec[j]
    sum_indirect <- 0
    valid_ks <- names(k_indices)
    
    if (length(valid_ks) > 0) {
      for (k in valid_ks) {
        pik <- ik_vec[k]
        pkj <- trans_mat[k, j]
        if (!is.na(pkj)) sum_indirect <- sum_indirect + (pik * pkj)
      }
      sum_indirect <- sum_indirect / length(valid_ks)
      #print(length(valid_ks))
    } else {
      sum_indirect <- 0
    }
    
    combined_vec[j] <- min(1, pj + sum_indirect)
  }
  
  combined_vec[county_upper] <- NA  # Hide self-transmission
  return(combined_vec)
}


# Figure 01 County -----

figure01_county <- function(method = 7, strategy = 0, county = "Gaines",
                     threshold = 0.5, map_data = infection_county,
                     out_dir = "Figures/") {
  
  
  dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)
  
  # Define breaks and labels
  breaks <- c(0, 0.2, 0.4, 0.6, 0.8, 1)
  bin_labels <- c("0–0.2", "0.2–0.4", "0.4–0.6", "0.6–0.8", "0.8–1.0")
  color_palette <- rev(c("#d73027", "#fc8d59", "#fee08b", "#91cf60", "#1a9850"))  # red → yellow → green
  
  mat <- readRDS(paste0("ProcessedData/county_pij_M", method, "_S", strategy, ".rds"))
  county_upper <- toupper(county)
  county_names <- toupper(map_data$County)
  
  direct_vec <- mat[county_upper, match(county_names, colnames(mat))]
  indirect_vec <- compute_indirect_risk_county(mat, county_name = county_upper, threshold = threshold)
  highlighted_geom <- map_data %>% filter(toupper(County) == county_upper)
  
  # Cut values into bins
  map_data <- map_data %>%
    mutate(
      direct_bin = cut(direct_vec[match(toupper(County), names(direct_vec))], breaks = breaks, labels = bin_labels, include.lowest = TRUE, right = FALSE),
      indirect_bin = cut(indirect_vec[match(toupper(County), names(indirect_vec))], breaks = breaks, labels = bin_labels, include.lowest = TRUE, right = FALSE)
    )
  
  # Plot A: Direct
  p1 <- ggplot(map_data) +
    geom_sf(aes(fill = direct_bin), color = "gray40", size = 0.1) +
    geom_sf(data = highlighted_geom, fill = "blue", color = "black", size = 0.3) +
    scale_fill_manual(values = color_palette, drop = FALSE, name = "Outbreak Probability") +
    labs(title = "A: Probability of First-generation Outbreak") +
    theme_void() +
    theme(
      plot.title = element_text(hjust = 0.5, size = 13),
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
      panel.grid = element_blank(),
      axis.text = element_blank(),
      axis.ticks = element_blank(),
      axis.title = element_blank()
    )
  
  legend <- get_legend(p2)
  p2_clean <- p2 + theme(legend.position = "none")
  row_plot <- plot_grid(p1, p2_clean, nrow = 1, rel_widths = c(1, 1))
  final_plot <- plot_grid(row_plot, legend, nrow = 1, rel_widths = c(2, 0.4))
  
  file_name <- paste0(out_dir, "figure01.png")
  ggsave(file_name, final_plot, width = 12, height = 6, dpi = 400)
  
}


figure01_county()


# Figure 01 District -----

figure01_district <- function(method = 7, strategy = 0, districts = c("LOOP ISD"),
                            threshold = 0.5, map_data = infection_district,
                            out_dir = "Figures/") {
  
  
  dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)
  
  # Define breaks and labels
  breaks <- c(0, 0.2, 0.4, 0.6, 0.8, 1)
  bin_labels <- c("0–0.2", "0.2–0.4", "0.4–0.6", "0.6–0.8", "0.8–1.0")
  color_palette <- rev(c("#d73027", "#fc8d59", "#fee08b", "#91cf60", "#1a9850"))  # red → yellow → green
  
  mat <- readRDS(paste0("ProcessedData/district_pij_M", method, "_S", strategy, ".rds"))
  county_upper <- toupper(county)
  district_names <- toupper(map_data$County)
  
  direct_vec <- mat[county_upper, match(county_names, colnames(mat))]
  indirect_vec <- compute_indirect_risk_county(mat, county_name = county_upper, threshold = threshold)
  highlighted_geom <- map_data %>% filter(toupper(County) == county_upper)
  
  # Cut values into bins
  map_data <- map_data %>%
    mutate(
      direct_bin = cut(direct_vec[match(toupper(County), names(direct_vec))], breaks = breaks, labels = bin_labels, include.lowest = TRUE, right = FALSE),
      indirect_bin = cut(indirect_vec[match(toupper(County), names(indirect_vec))], breaks = breaks, labels = bin_labels, include.lowest = TRUE, right = FALSE)
    )
  
  # Plot A: Direct
  p1 <- ggplot(map_data) +
    geom_sf(aes(fill = direct_bin), color = "gray40", size = 0.1) +
    geom_sf(data = highlighted_geom, fill = "blue", color = "black", size = 0.3) +
    scale_fill_manual(values = color_palette, drop = FALSE, name = "Outbreak Probability") +
    labs(title = "A: Probability of First-generation Outbreak") +
    theme_void() +
    theme(
      plot.title = element_text(hjust = 0.5, size = 13),
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
      panel.grid = element_blank(),
      axis.text = element_blank(),
      axis.ticks = element_blank(),
      axis.title = element_blank()
    )
  
  legend <- get_legend(p2)
  p2_clean <- p2 + theme(legend.position = "none")
  row_plot <- plot_grid(p1, p2_clean, nrow = 1, rel_widths = c(1, 1))
  final_plot <- plot_grid(row_plot, legend, nrow = 1, rel_widths = c(2, 0.4))
  
  file_name <- paste0(out_dir, "figure01.png")
  ggsave(file_name, final_plot, width = 12, height = 6, dpi = 400)
  
}


figure01_county()

# Figure 02-----

figure02 <- function(method = 7, counties = c("Gaines"), strategies = c(1, 2, 3),
                     map_data = infection_county, out_dir = "Figures/") {
  
  dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)
  
  color_palette <- c("white", "#1a9850")  # white to green
  county_names <- toupper(map_data$County)
  
  full_labels <- c(
    "Strategy 0: Baseline",
    "Strategy 1: MMR < 90% → 90%", 
    "Strategy 2: MMR < 92% → 92%", 
    "Strategy 3: MMR < 95% → 95%",
    "Strategy 4: Gaines → 90%",
    "Strategy 5: Gaines → 92%"
  )
  
  strategy_labels <- full_labels[strategies + 1]
  label_prefix <- LETTERS[seq_along(strategies)]
  
  baseline_mat <- readRDS(paste0("ProcessedData/county_pij_M", method, "_S0.rds"))
  
  for (county in counties) {
    county_upper <- toupper(county)
    highlighted_geom <- map_data %>% filter(toupper(County) == county_upper)
    baseline_vec <- baseline_mat[county_upper, match(county_names, colnames(baseline_mat))]
    
    maps <- list()
    
    # Compute global max reduction across all strategies for this county
    global_max <- 0
    for (s in strategies) {
      mat_s <- readRDS(paste0("ProcessedData/county_pij_M", method, "_S", s, ".rds"))
      vec_s <- mat_s[county_upper, match(county_names, colnames(mat_s))]
      red_s <- (baseline_vec - vec_s)*100
      global_max <- max(global_max, max(red_s, na.rm = TRUE))
    }
    
    
    for (i in seq_along(strategies)) {
      s <- strategies[i]
      mat <- readRDS(paste0("ProcessedData/county_pij_M", method, "_S", s, ".rds"))
      vec <- mat[county_upper, match(county_names, colnames(mat))]
      reduction <- (baseline_vec - vec)*100
      #reduction[reduction < 0] <- 0
      
      map_data$reduction <- reduction[match(toupper(map_data$County), names(reduction))]
      
      p <- ggplot(map_data) +
        geom_sf(aes(fill = reduction), color = "gray40", size = 0.1) +
        geom_sf(data = highlighted_geom, fill = "blue", color = "black", size = 0.3) +
        scale_fill_gradient(
          low = "white",
          high = "#1a9850", 
          #limits = c(0, max(reduction, na.rm = TRUE)),
          limits = c(0, global_max),
          name = "Reduction in\nOutbreak Probability (%)", na.value = "gray80"
        ) +
        labs(title = paste0(label_prefix[i], ": ", strategy_labels[i])) +
        theme_void() +
        theme(
          plot.title = element_text(hjust = 0.5, size = 13),
          legend.position = "none"
        )
      
      maps[[i]] <- p
    }
    
    legend <- get_legend(
      maps[[2]] + theme(legend.position = "right", legend.title = element_text(size = 10))
    )
    
    row_plot <- plot_grid(plotlist = lapply(maps, \(p) p + theme(legend.position = "none")), nrow = 1)
    final_plot <- plot_grid(row_plot, legend, nrow = 1, rel_widths = c(length(maps), 0.4))
    
    file_name <- paste0(out_dir, "figure02_", tolower(county), "_method", method, "_S", paste(strategies, collapse = "_"), ".png")
    ggsave(file_name, final_plot, width = 6 * length(maps), height = 6, dpi = 400)
  }
}


figure02(method = 7, counties = c("Gaines"), strategies = c(1, 2, 3), map_data = infection_county)
















# Figure 03-----
# Plot low vaccine rate counties in a roew indirect risk

plot_indirect_transmission_row_multiple <- function(method, map_data, counties,
                                                    strategy = 0, threshold = 0.5,
                                                    out_dir = "Figures/") {
  
  base_colors <- c("#1a9850", "#91cf60", "#d9ef8b", "#fee08b", "#fc8d59", "#d73027")
  county_names <- toupper(map_data$County)
  
  mat <- readRDS(paste0("ProcessedData/county_pij_M", method, "_S", strategy, ".rds"))
  maps <- list()
  
  for (county in counties) {
    county_upper <- toupper(county)
    highlighted_geom <- map_data %>% filter(toupper(County) == county_upper)
    
    indirect_vec <- compute_indirect_risk_county(mat, county_name = county_upper, threshold = threshold)
    map_data$indirect_risk <- indirect_vec[match(county_names, names(indirect_vec))]
    
    p <- ggplot(map_data) +
      geom_sf(aes(fill = indirect_risk), color = "gray40", size = 0.1) +
      geom_sf(data = highlighted_geom, fill = "blue", color = "black", size = 0.3) +
      scale_fill_gradientn(colors = base_colors, na.value = "gray80", limits = c(0, 1)) +
      labs(
        title = str_to_title(county),
        fill = "Outbreak Probability"
      ) +
      theme_minimal() +
      theme(
        plot.title = element_text(hjust = 0.5, size = 13),
        legend.title = element_text(size = 10),
        legend.text = element_text(size = 9)
      )
    
    maps[[length(maps) + 1]] <- p
  }
  
  combined_plot <- plot_grid(plotlist = maps, nrow = 1)
  
  file_name <- paste0(out_dir, "indirect_risk_multiple_method", method,
                      "_S", strategy, "_", paste(tolower(counties), collapse = "_"), ".png")
  ggsave(file_name, combined_plot, width = 6 * length(counties), height = 6, dpi = 400)
  
  print(combined_plot)
}



plot_indirect_transmission_row_multiple(
  method = 7,
  map_data = infection_county,
  #counties = c("POLK", "MONTAGUE", "LIMESTONE"),
  #counties = c("KING", "HALL", "THROCKMORTON"),
  counties = c("DALLAS", "TARRANT", "COLLIN"),
  strategy = 0,
  threshold = 0.5
)






# Figure 04 ---- 
# compare vaccine rate, infection risk and outbreak
figure04 <- function(method = 7, strategy = 0, county = "Gaines",
                            threshold = 0.5, map_data = infection_county,
                            out_dir = "Figures/") {
  
  
  dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)
  
  # Define breaks and labels
  breaks <- c(0, 0.2, 0.4, 0.6, 0.8, 1)
  bin_labels <- c("0–0.2", "0.2–0.4", "0.4–0.6", "0.6–0.8", "0.8–1.0")
  color_palette <- rev(c("#d73027", "#fc8d59", "#fee08b", "#91cf60", "#1a9850"))  # red → yellow → green
  
  mat <- readRDS(paste0("ProcessedData/county_pij_M", method, "_S", strategy, ".rds"))
  county_upper <- toupper(county)
  county_names <- toupper(map_data$County)
  
  direct_vec <- mat[county_upper, match(county_names, colnames(mat))]
  indirect_vec <- compute_indirect_risk_county(mat, county_name = county_upper, threshold = threshold)
  highlighted_geom <- map_data %>% filter(toupper(County) == county_upper)
  
  # Cut values into bins
  map_data <- map_data %>%
    mutate(
      direct_bin = cut(direct_vec[match(toupper(County), names(direct_vec))], breaks = breaks, labels = bin_labels, include.lowest = TRUE, right = FALSE),
      indirect_bin = cut(indirect_vec[match(toupper(County), names(indirect_vec))], breaks = breaks, labels = bin_labels, include.lowest = TRUE, right = FALSE)
    )
  
  
  # Plot A: Vaccine Coverage
  p1 <- ggplot(map_data) +
    geom_sf(
      aes(fill = MMR),
      color = "gray30",
      size = 0.1
    ) +
    scale_fill_gradientn(
      colors = c("#d73027", "#fee08b", "#1a9850"),  # soft red → light yellow → muted green
      limits = c(0.8, 1.0),
      na.value = "lightgray",
      labels = scales::percent_format(accuracy = 1)
    ) +
    theme_void() +
    labs(
      title = "MMR Coverage (2024)",
      fill = "MMR                     "
    )
  
  p2 <- ggplot(map_data) +
    geom_sf(
      aes(fill = adult_infection_prob),
      color = "gray30",
      size = 0.1
    ) +
    scale_fill_gradientn(
      colors = c("#1a9850", "#fee08b", "#d73027"),  
      limits = c(0, .15),
      na.value = "lightgray",
      labels = scales::percent_format(accuracy = 1)
    ) +
    theme_void() +
    labs(
      title = "Infection probability (above 18 years old)",
      fill = "Infection Probability"
    )
  
  
  # Plot C: Direct
  p3 <- ggplot(map_data) +
    geom_sf(
      aes(fill = direct_bin), 
      color = "gray40", 
      size = 0.1
      ) +
    geom_sf(data = highlighted_geom, 
            fill = "blue", 
            color = "black", 
            size = 0.3
            ) +
    scale_fill_manual(
      values = color_palette, 
      drop = FALSE, 
      name = "Outbreak Probability"
      ) +
    labs(title = "Probability of First-generation Outbreak") +
    theme_void() +
    theme(
      plot.title = element_text(hjust = 0.5, size = 13)
      #legend.position = "none"  
    )
  
  
  
  #legend <- get_legend(p2)
  #p2_clean <- p2 + theme(legend.position = "none")
  row_plot <- plot_grid(p1, p2, p3, nrow = 1, rel_widths = c(1, 1, 1))
  final_plot <- plot_grid(
    row_plot, 
    #legend, 
    nrow = 1, 
    rel_heights = c(1, 0.12)
    )
  
  file_name <- paste0(out_dir, "figure04.png")
  ggsave(file_name, final_plot, width = 12, height = 6, dpi = 400)
  
}


figure04()


compute_indirect_risk_county_new <- function(trans_mat, county_name = "GAINES", threshold = .5) {
  county_upper <- toupper(county_name)
  county_names <- rownames(trans_mat)
  
  if (!(county_upper %in% county_names))
    stop("County not found in matrix")
  
  # Direct transmission vector j→i, county_upper is the origin county
  Pij_vec <- trans_mat[ ,county_upper]
  
  # Identify candidate first-generation counties k
  Pjk_vec  <- trans_mat[, county_upper]          # Same row: j → k
  k_indices <- which(Pjk_vec >= threshold & !is.na(Pjk_vec))
  K         <- names(k_indices)
  
  # Output vector
  SG_vec <- rep(NA_real_, length(Pij_vec))
  names(SG_vec) <- names(Pij_vec)
  
  for (i in names(Pij_vec)) {
    
    Pji <- Pij_vec[i]  # direct probability j → i
    
    
    compute_indirect_risk_county_new <- function(trans_mat, county_name = "GAINES", threshold = -1) {
      county_upper <- toupper(county_name)
      county_names <- rownames(trans_mat)
      
      if (!(county_upper %in% county_names))
        stop("County not found in matrix")
      
      # Direct transmission vector P_j→i
      Pij_vec <- trans_mat[county_upper, ]
      
      # Identify candidate first-generation counties k
      Pkj_vec  <- trans_mat[county_upper, ]          # Same row: j → k
      k_indices <- which(Pkj_vec >= threshold & !is.na(Pkj_vec))
      K         <- names(k_indices)
      
      # Output vector
      SG_vec <- rep(NA_real_, length(Pij_vec))
      names(SG_vec) <- names(Pij_vec)
      
      for (i in names(Pij_vec)) {
        
        pij <- Pij_vec[i]  # direct probability j → i
        
        all_k <- names(Pkj_vec)[!is.na(Pkj_vec)]
        K <- setdiff(all_k, c(county_upper, i))
        
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
          Pjk <- Pkj_vec[k]  # j → k
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
      Pjk <- Pkj_vec[k]  # j → k
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

