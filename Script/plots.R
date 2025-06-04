library(ggplot2)
library(sf)
library(cowplot)
library(stringr)
library(dplyr)


map_probability <- readRDS("ProcessedData/map_probability.rds")

## Plot transmission risk


plot_transmission_row <- function(method, map_data, counties, strategies = c(0,1,2,3), out_dir = "Figures/") {
  dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)
  
  method_label <- ifelse(method == 3, "Gravity Model", "52-Week Mobility Flow")
  base_colors <- c("#1a9850", "#91cf60", "#d9ef8b", "#fee08b", "#fc8d59", "#d73027")
  county_names <- toupper(map_data$County)
  
  full_labels <- c(
    "Strategy 0: Baseline",
    "Strategy 1: MMR < 90% → 90%", 
    "Strategy 2: MMR < 92% → 92%", 
    "Strategy 3: MMR +5%",
    "Strategy 4: Gaines → 90%",
    "Strategy 5: Gaines → 92%"
  )
  
  strategy_labels <- full_labels[strategies + 1]
  
  # Load matrices for selected strategies
  mats <- lapply(strategies, function(s) readRDS(paste0("ProcessedData/pij_M", method, "_S", s, ".rds")))
  
  for (county in counties) {
    county_upper <- toupper(county)
    county_title <- str_to_title(tolower(county))
    highlighted_geom <- map_data %>% filter(toupper(County) == county_upper)
    
    maps <- list()
    
    for (i in seq_along(strategies)) {
      s <- strategies[i]
      mat <- mats[[i]]
      strategy_col <- paste0("method", method, "_s", s)
      
      map_data[[strategy_col]] <- mat[county_upper, match(county_names, colnames(mat))]
      
      p <- ggplot(map_data) +
        geom_sf(aes_string(fill = strategy_col), color = "gray40", size = 0.1) +
        scale_fill_gradientn(colors = base_colors, na.value = "gray80", limits = c(0, 1)) +
        labs(
          title = strategy_labels[i],
          fill = "Outbreak Probability"
        ) +
        geom_sf(data = highlighted_geom, fill = "blue", color = "black", size = 0.3) +
        theme_minimal() +
        theme(
          plot.title = element_text(hjust = 0.5, size = 13),
          legend.title = element_text(size = 10),
          legend.text = element_text(size = 9)
        )
      
      maps[[i]] <- p
    }
    
    row_plot <- plot_grid(plotlist = maps, nrow = 1)
    
    file_name <- paste0(out_dir, "row_transmission_", tolower(county), "_method", method, "_S", paste(strategies, collapse = "_"), ".png")
    ggsave(file_name, row_plot, width = 6 * length(strategies), height = 6, dpi = 400)
    
    print(row_plot)
  }
}




plot_transmission_row(method = 3, map_probability, counties = c("Gaines"), strategies = c(0, 4, 5))
plot_transmission_row(method = 7, map_probability, counties = c("Gaines"), strategies = c(0, 4, 5))

plot_transmission_row(method = 3, map_probability, counties = c("Gaines"), strategies = c(1, 2, 3))
plot_transmission_row(method = 7, map_probability, counties = c("Gaines"), strategies = c(1, 2, 3))

plot_transmission_row(method = 3, map_probability, counties = c("Gaines"), strategies = c(0))
plot_transmission_row(method = 7, map_probability, counties = c("Gaines"), strategies = c(0))


plot_transmission_row(method = 1, map_probability, counties = c("Gaines"), strategies = c(0, 4, 5))
plot_transmission_row(method = 2, map_probability, counties = c("Gaines"), strategies = c(0, 4, 5))

plot_transmission_row(method = 1, map_probability, counties = c("Gaines"), strategies = c(1, 2, 3))
plot_transmission_row(method = 2, map_probability, counties = c("Gaines"), strategies = c(1, 2, 3))


plot_transmission_row(method = 1, map_probability, counties = c("Gaines"), strategies = c(0))
plot_transmission_row(method = 2, map_probability, counties = c("Gaines"), strategies = c(0))

############################# Histograms ############
plot_transmission_histograms <- function(method, map_data, counties, strategies = 0:3, out_dir = "Figures/") {
  dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)
  
  # Full strategy labels including 4 and 5
  full_labels <- c(
    "Strategy 0: Baseline",
    "Strategy 1: MMR < 90% → 90%", 
    "Strategy 2: MMR < 92% → 92%", 
    "Strategy 3: MMR +5% (max 100%)",
    "Strategy 4: Gaines → 90%",
    "Strategy 5: Gaines → 92%"
  )
  
  strategy_labels <- full_labels[strategies + 1]
  base_colors <- "#fc8d59"
  county_names <- toupper(map_data$County)
  
  # Load only selected matrices
  mats <- lapply(strategies, function(s) readRDS(paste0("ProcessedData/pij_M", method, "_S", s, ".rds")))
  
  for (county in counties) {
    county_upper <- toupper(county)
    hist_plots <- list()
    
    # Compute global max y across all strategies
    all_bin_counts <- c()
    for (k in seq_along(strategies)) {
      prob_vec <- mats[[k]][county_upper, ]
      prob_vec <- prob_vec[!is.na(prob_vec) & prob_vec >= 0.05]
      hist_data <- hist(prob_vec, breaks = seq(0, 1, by = 0.02), plot = FALSE)
      all_bin_counts <- c(all_bin_counts, hist_data$counts)
    }
    max_y <- max(all_bin_counts)
    
    # Generate plots
    for (i in seq_along(strategies)) {
      prob_vec <- mats[[i]][county_upper, ]
      prob_vec <- prob_vec[!is.na(prob_vec) & prob_vec >= 0.05]
      df <- data.frame(probability = prob_vec)
      
      p <- ggplot(df, aes(x = probability)) +
        geom_histogram(binwidth = 0.02, fill = base_colors, color = "black", boundary = 0) +
        xlim(0, 1) +
        ylim(0, max_y) +
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
    
    file_name <- paste0(out_dir, "histogram_transmission_", tolower(county), "_method", method, "_S", paste(strategies, collapse = "_"), ".png")
    ggsave(file_name, combined_plot, width = 6, height = 3 * length(strategies), dpi = 400)
    
    print(combined_plot)
  }
}



plot_transmission_histograms(method = 3, map_probability, counties = c("Gaines"), strategies = c(0, 4, 5))
plot_transmission_histograms(method = 7, map_probability, counties = c("Gaines"), strategies = c(0, 4, 5))


plot_transmission_histograms(method = 3, map_probability, counties = c("Gaines"), strategies = c(0, 1, 2, 3))
plot_transmission_histograms(method = 7, map_probability, counties = c("Gaines"), strategies = c(0, 1, 2, 3))


plot_transmission_histograms(method = 3, map_probability, counties = c("Gaines"), strategies = c(0))
plot_transmission_histograms(method = 7, map_probability, counties = c("Gaines"), strategies = c(0))


plot_transmission_histograms(method = 1, map_probability, counties = c("Gaines"), strategies = c(0, 4, 5))
plot_transmission_histograms(method = 2, map_probability, counties = c("Gaines"), strategies = c(0, 4, 5))


plot_transmission_histograms(method = 1, map_probability, counties = c("Gaines"), strategies = c(0, 1, 2, 3))
plot_transmission_histograms(method = 2, map_probability, counties = c("Gaines"), strategies = c(0, 1, 2, 3))


plot_transmission_histograms(method = 1, map_probability, counties = c("Gaines"), strategies = c(0))
plot_transmission_histograms(method = 2, map_probability, counties = c("Gaines"), strategies = c(0))


########## Histogram suceptible

plot_transmission_population_risk_histograms <- function(method, map_probability, counties,
                                                         strategies = 0:3,
                                                         breaks = seq(0, 1, by = 0.1),
                                                         out_dir = "Figures/") {

  dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)
  
  full_labels <- c(
    "Strategy 0: Baseline",
    "Strategy 1: MMR < 90% → 90%", 
    "Strategy 2: MMR < 92% → 92%", 
    "Strategy 3: MMR +5%",
    "Strategy 4: Gaines → 90%",
    "Strategy 5: Gaines → 92%"
  )
  strategy_labels <- full_labels[strategies + 1]
  base_colors <- "#fc8d59"
  
  # Load only selected transmission matrices
  mats <- lapply(strategies, function(s) readRDS(paste0("ProcessedData/pij_M", method, "_S", s, ".rds")))
  county_names <- toupper(map_probability$County)
  full_range_levels <- levels(cut(breaks[-length(breaks)], breaks = breaks, include.lowest = TRUE, right = FALSE))
  
  for (county in counties) {
    county_upper <- toupper(county)
    hist_plots <- list()
    all_y_vals <- c()
    
    # First pass to compute global max y
    y_vals_per_strategy <- list()
    
    for (i in seq_along(strategies)) {
      pij <- mats[[i]]
      prob_vec <- pij[county_upper, ]
      prob_vec <- prob_vec[!is.na(prob_vec)]
      target_names <- names(prob_vec)
      susceptible_col <- paste0("susceptible_pop_size", ifelse(strategies[i] == 0, "", strategies[i]))
      
      df <- map_probability %>%
        mutate(county_upper = toupper(County)) %>%
        filter(county_upper %in% target_names) %>%
        mutate(prob = prob_vec[match(county_upper, target_names)]) %>%
        filter(!is.na(prob), !is.na(!!sym(susceptible_col))) %>%
        mutate(range = cut(prob, breaks = breaks, include.lowest = TRUE, right = FALSE)) %>%
        group_by(range) %>%
        summarise(total_at_risk = sum(!!sym(susceptible_col)), .groups = "drop") %>%
        complete(range = factor(full_range_levels, levels = full_range_levels), fill = list(total_at_risk = 0))
      
      y_vals_per_strategy[[i]] <- df
      all_y_vals <- c(all_y_vals, df$total_at_risk)
    }
    
    max_y <- max(all_y_vals)
    
    # Plot with unified y-axis
    for (i in seq_along(strategies)) {
      df <- y_vals_per_strategy[[i]]
      
      p <- ggplot(df, aes(x = range, y = total_at_risk)) +
        geom_col(fill = base_colors, color = "black") +
        ylim(0, max_y) +
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
    
    file_name <- paste0(out_dir, "population_risk_histogram_", tolower(county), "_method", method, "_S", paste(strategies, collapse = "_"), ".png")
    ggsave(file_name, combined_plot, width = 6, height = 3 * length(strategies), dpi = 400)
    
    print(combined_plot)
  }
}




plot_transmission_population_risk_histograms(method = 3, map_probability, counties = c("Gaines"), strategies = c(0, 4, 5))
plot_transmission_population_risk_histograms(method = 7, map_probability, counties = c("Gaines"), strategies = c(0, 4, 5))

plot_transmission_population_risk_histograms(method = 3, map_probability, counties = c("Gaines"), strategies = c(0, 1, 2, 3))
plot_transmission_population_risk_histograms(method = 7, map_probability, counties = c("Gaines"), strategies = c(0, 1, 2, 3))

plot_transmission_population_risk_histograms(method = 3, map_probability, counties = c("Gaines"), strategies = c(0))
plot_transmission_population_risk_histograms(method = 7, map_probability, counties = c("Gaines"), strategies = c(0))



plot_transmission_population_risk_histograms(method = 1, map_probability, counties = c("Gaines"), strategies = c(0, 4, 5))
plot_transmission_population_risk_histograms(method = 2, map_probability, counties = c("Gaines"), strategies = c(0, 4, 5))

plot_transmission_population_risk_histograms(method = 1, map_probability, counties = c("Gaines"), strategies = c(0, 1, 2, 3))
plot_transmission_population_risk_histograms(method = 2, map_probability, counties = c("Gaines"), strategies = c(0, 1, 2, 3))

plot_transmission_population_risk_histograms(method = 1, map_probability, counties = c("Gaines"), strategies = c(0))
plot_transmission_population_risk_histograms(method = 2, map_probability, counties = c("Gaines"), strategies = c(0))

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



plot_map(map_probability,
         "MMR",
         "Measles, Mumps, and Rubella Vaccine Rate", 
         file_name="mmr",
         legend_title = "Vaccine Rate",
         reverse_colors = TRUE)


plot_map(map_probability,
         "outbreak_prob",
         "Local Outbreak Probability", 
         file_name="outbreak",
         legend_title = "Local Outbreak Probability",
         reverse_colors = FALSE)



