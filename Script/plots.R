library(ggplot2)
library(sf)
library(cowplot)
library(stringr)
library(dplyr)


map_data <- readRDS("ProcessedData/map_county.rds")

## Plot transmission risk


plot_transmission_row <- function(method, map_data, counties, out_dir = "Figures/") {
  dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)
  
  strategy_labels <- c(
    "Strategy 0: Baseline",
    "Strategy 1: MMR < 90% → 90%", 
    "Strategy 2: MMR < 92% → 92%", 
    "Strategy 3: MMR +5%"
  )
  
  method_label <- ifelse(method == 3, "Gravity Model", "52-Week Mobility Flow")
  
  # Load matrices
  mat_s0 <- readRDS(paste0("ProcessedData/pij_M", method, "_S0.rds"))
  mat_s1 <- readRDS(paste0("ProcessedData/pij_M", method, "_S1.rds"))
  mat_s2 <- readRDS(paste0("ProcessedData/pij_M", method, "_S2.rds"))
  mat_s3 <- readRDS(paste0("ProcessedData/pij_M", method, "_S3.rds"))
  
  county_names <- toupper(map_data$County)
  base_colors <- c("#1a9850", "#91cf60", "#d9ef8b", "#fee08b", "#fc8d59", "#d73027")
  
  for (county in counties) {
    county_upper <- toupper(county)
    county_title <- str_to_title(tolower(county))
    
    maps <- list()
    
    # Strategy 0 (baseline)
    strategy_col0 <- paste0("method", method, "_s0")
    map_data[[strategy_col0]] <- mat_s0[county_upper, match(county_names, colnames(mat_s0))]
    highlighted_geom <- map_data %>% filter(toupper(County) == county_upper)
    
    baseline_plot <- ggplot(map_data) +
      geom_sf(aes_string(fill = strategy_col0), color = "gray40", size = 0.1) +
      scale_fill_gradientn(colors = base_colors, na.value = "gray80", limits = c(0, 1)) +
      labs(
        title = strategy_labels[1],
        fill = "Outbreak Probability"
      ) +
      geom_sf(data = highlighted_geom, fill = "blue", color = "black", size = 0.3) +
      theme_minimal() +
      theme(
        plot.title = element_text(hjust = 0.5, size = 13),
        legend.title = element_text(size = 10),
        legend.text = element_text(size = 9)
      )
    
    # Save the baseline plot separately
    ggsave(
      paste0(out_dir, "transmission_", tolower(county), "_method", method, "_baseline.png"),
      baseline_plot, width = 6, height = 6, dpi = 400
    )
    
    # Strategies 1–3
    for (i in 1:3) {
      mat <- list(mat_s1, mat_s2, mat_s3)[[i]]
      strategy_col <- paste0("method", method, "_s", i)
      
      map_data[[strategy_col]] <- mat[county_upper, match(county_names, colnames(mat))]
      
      p <- ggplot(map_data) +
        geom_sf(aes_string(fill = strategy_col), color = "gray40", size = 0.1) +
        scale_fill_gradientn(colors = base_colors, na.value = "gray80", limits = c(0, 1)) +
        labs(
          title = strategy_labels[i + 1],
          fill = NULL
        ) +
        geom_sf(data = highlighted_geom, fill = "blue", color = "black", size = 0.3) +
        theme_minimal() +
        theme(
          plot.title = element_text(hjust = 0.5, size = 11),
          legend.position = "none"
        )
      
      maps[[i]] <- p
    }
    
    # Combine the three maps into one row
    row_plot <- plot_grid(plotlist = maps, nrow = 1)
    
    # Save and show
    file_name <- paste0(out_dir, "row_transmission_", tolower(county), "_method", method, ".png")
    ggsave(file_name, row_plot, width = 16, height = 6, dpi = 400)
    
    print(row_plot)
    print(baseline_plot)
  }
}



plot_transmission_row(method = 3, map_data, counties = c("Gaines"))
plot_transmission_row(method = 7, map_data, counties = c("Gaines"))



############################# Histograms ############
plot_transmission_histograms <- function(method, map_data, counties, out_dir = "Figures/") {
  stopifnot(method %in% c(3, 7))
  dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)
  
  strategy_labels <- c(
    "Strategy 0: Baseline",
    "Strategy 1: MMR < 90% → 90%", 
    "Strategy 2: MMR < 92% → 92%", 
    "Strategy 3: MMR +5% (max 100%)"
  )
  
  method_label <- ifelse(method == 3, "Gravity Model", "52-Week Mobility Flow")
  base_colors <- "#fc8d59"
  
  # Load matrices
  mats <- lapply(0:3, function(s) readRDS(paste0("ProcessedData/pij_M", method, "_S", s, ".rds")))
  
  county_names <- toupper(map_data$County)
  
  for (county in counties) {
    county_upper <- toupper(county)
    county_title <- str_to_title(tolower(county))
    
    hist_plots <- list()
    
    for (i in 1:4) {
      prob_vec <- mats[[i]][county_upper, ]
      prob_vec <- prob_vec[!is.na(prob_vec) & prob_vec >= 0.05]
      
      df <- data.frame(probability = prob_vec)
      
      p <- ggplot(df, aes(x = probability)) +
        geom_histogram(binwidth = 0.02, fill = base_colors, color = "black", boundary = 0) +
        xlim(0, 1) +
        labs(
          title = strategy_labels[i],
          x = "Transmission Probability",
          y = "Number of Counties"
        ) +
        theme_minimal() +
        theme(
          plot.title = element_text(size = 12, hjust = 0.5),
          axis.text = element_text(size = 10),
          axis.title = element_text(size = 11)
        )
      
      hist_plots[[i]] <- p
    }
    
    combined_plot <- plot_grid(plotlist = hist_plots, ncol = 1)
    
    file_name <- paste0(out_dir, "histogram_transmission_", tolower(county), "_method", method, ".png")
    ggsave(file_name, combined_plot, width = 6, height = 10, dpi = 400)
    
    print(combined_plot)
  }
}


plot_transmission_histograms(method = 3, map_data, counties = c("Gaines"))
plot_transmission_histograms(method = 7, map_data, counties = c("Gaines"))



########## Histogram suceptible

plot_transmission_population_risk_histograms <- function(method, map_probability, counties,
                                                         breaks = seq(0, 1, by = 0.1),
                                                         out_dir = "Figures/") {
  stopifnot(method %in% c(3, 7))
  dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)
  
  strategy_labels <- c(
    "Strategy 0: Baseline",
    "Strategy 1: MMR < 90% → 90%", 
    "Strategy 2: MMR < 92% → 92%", 
    "Strategy 3: MMR +5%"
  )
  
  method_label <- ifelse(method == 3, "Gravity Model", "52-Week Mobility Flow")
  base_colors <- "#fc8d59"
  
  # Load transmission matrices for all 4 strategies
  mats <- lapply(0:3, function(s) readRDS(paste0("ProcessedData/pij_M", method, "_S", s, ".rds")))
  county_names <- toupper(map_probability$County)
  
  full_range_levels <- levels(cut(breaks[-length(breaks)], breaks = breaks, include.lowest = TRUE, right = FALSE))
  
  for (county in counties) {
    county_upper <- toupper(county)
    county_title <- str_to_title(tolower(county))
    hist_plots <- list()
    
    for (i in 1:4) {
      pij <- mats[[i]]
      prob_vec <- pij[county_upper, ]
      prob_vec <- prob_vec[!is.na(prob_vec)]
      target_names <- names(prob_vec)
      
      susceptible_col <- paste0("susceptible_pop_size", ifelse(i == 1, "", i - 1))
      
      df <- map_probability %>%
        mutate(county_upper = toupper(County)) %>%
        filter(county_upper %in% target_names) %>%
        mutate(prob = prob_vec[match(county_upper, target_names)]) %>%
        filter(!is.na(prob), !is.na(!!sym(susceptible_col))) %>%
        mutate(range = cut(prob, breaks = breaks, include.lowest = TRUE, right = FALSE)) %>%
        group_by(range) %>%
        summarise(total_at_risk = sum(!!sym(susceptible_col)), .groups = "drop") %>%
        complete(range = factor(full_range_levels, levels = full_range_levels), fill = list(total_at_risk = 0))
      
      p <- ggplot(df, aes(x = range, y = total_at_risk)) +
        geom_col(fill = base_colors, color = "black") +
        labs(
          title = strategy_labels[i],
          x = "Transmission Probability Range",
          y = "Total Susceptible Population"
        ) +
        theme_minimal() +
        theme(
          plot.title = element_text(size = 12, hjust = 0.5),
          axis.text.x = element_text(angle = 45, hjust = 1),
          axis.text = element_text(size = 10),
          axis.title = element_text(size = 11)
        )
      
      hist_plots[[i]] <- p
    }
    
    combined_plot <- plot_grid(plotlist = hist_plots, ncol = 1)
    
    file_name <- paste0(out_dir, "population_risk_histogram_", tolower(county), "_method", method, ".png")
    ggsave(file_name, combined_plot, width = 6, height = 12, dpi = 400)
    
    print(combined_plot)
  }
}



plot_transmission_population_risk_histograms(method = 3, map_probability, counties = c("Gaines"))
plot_transmission_population_risk_histograms(method = 7, map_probability, counties = c("Gaines"))



#### Plot MMR 


plot_map <- function(data, value_col, plot_title, file_name,
                              limits = NULL, legend_title = NULL, reverse_colors = FALSE) {
  
  base_colors <- c("#1a9850", "#91cf60", "#d9ef8b", "#fee08b", "#fc8d59", "#d73027")
  #base_colors <- c("#3288bd","#99d594", "#e6f598", "#fee08b", "#fc8d59", "#d53e4f")
  custom_colors <- if (reverse_colors) rev(base_colors) else base_colors
  
  if (is.null(legend_title)) legend_title <- value_col
  
  p <- ggplot(data) +
    geom_sf(aes_string(fill = value_col), color = "gray40", size = 0.1) +
    scale_fill_gradientn(colors = custom_colors, na.value = "gray80", limits = limits) +
    labs(title = plot_title, fill = legend_title) +
    theme_minimal() +
    theme(
      plot.title = element_text(hjust = 0.5, size = 14),
      legend.title = element_text(size = 10),
      legend.text = element_text(size = 9)
    )
  
  ggsave(paste0("Figures/", file_name, ".png"), plot = p, width = 10, height = 6, dpi = 300)
  
  return(p)
}



plot_map(map_data,
         "MMR",
         "Measles, Mumps, and Rubella Vaccine Rate", 
         file_name="mmr",
         legend_title = "Vaccine Rate",
         reverse_colors = TRUE)


plot_map(map_data,
         "outbreak_prob",
         "Local Outbreak Probability", 
         file_name="outbreak",
         legend_title = "Local Outbreak Probability",
         reverse_colors = FALSE)



