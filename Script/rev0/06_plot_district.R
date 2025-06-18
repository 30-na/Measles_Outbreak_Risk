library(gridExtra)
library(ggplot2)
library(cowplot)
library(sf)
library(stringr)


# Load data
#map_data <- readRDS("ProcessedData/map_data.rds")
map_data_county <- readRDS("ProcessedData/map_data_county.rds")

plot_district_map <- function(data, value_col, plot_title, file_name,
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
  
  #ggsave(paste0("Figures/", file_name, ".png"), plot = p, width = 10, height = 6, dpi = 300)
  
  return(p)
}





g1 <- plot_district_map(map_data,
                  "MMR",
                  "Measles, Mumps, and Rubella Vaccine Rate", 
                  file_name="mmr",
                  legend_title = "Vaccine Rate",
                  reverse_colors = TRUE)


g2 <- plot_district_map(map_data,
                        "outbreak_prob",
                        "Local Outbreak Probability", 
                        file_name="outbreak",
                        legend_title = "Outbreak Probability")


# g3 <- plot_district_map(map_data,
#                        "internal_infection_prob_log10",
#                        "Log10 of Average Attack Rate", 
#                        file_name="internal_infection_probability_97_log10",
#                        legend_title = "Average Attack Rate")


g3 <- plot_district_map(map_data,
                        "pop_time_infection_total_log10",
                        "Log10 of Average Outbreak Size", 
                        file_name="infection_times_total_popualtion_97_log10",
                        legend_title = "Average Outbreak Size")





# Align and combine plots
combined_plot1 <- plot_grid(g1, g2, g3,
                            ncol = 3, 
                            align = "hv",
                            axis = "tblr",
                            rel_widths = c(1,1,1))

ggsave("Figures/combined_map1.png", 
       combined_plot1,
       width = 16,
       height = 6,
       dpi = 400)


####  COUTNY LVL

g1 <- plot_district_map(map_data_county,
                        "MMR",
                        "Measles, Mumps, and Rubella Vaccine Rate", 
                        file_name="mmr",
                        legend_title = "Vaccine Rate",
                        reverse_colors = TRUE)


g2 <- plot_district_map(map_data_county,
                        "outbreak_prob",
                        "Local Outbreak Probability", 
                        file_name="outbreak",
                        legend_title = "Outbreak Probability")


g3 <- plot_district_map(map_data,
                       "internal_infection_prob_log10",
                       "Log10 of Average Attack Rate",
                       file_name="internal_infection_probability_97_log10",
                       legend_title = "Average Attack Rate")


g4 <- plot_district_map(map_data,
                        "pop_time_infection_total_log10",
                        "Log10 of Average Outbreak Size", 
                        file_name="infection_times_total_popualtion_97_log10",
                        legend_title = "Average Outbreak Size")





# Align and combine plots
combined_plot1 <- plot_grid(g1, g2, g3, g4,
                            ncol = 2, 
                            align = "hv",
                            axis = "tblr",
                            rel_widths = c(1,1,1))

ggsave("Figures/combined_county_map1_0.95.png", 
       combined_plot1,
       width = 16,
       height = 6,
       dpi = 400)



# Load map and matrices

p1 <- readRDS("processedData/transmission_matrix_county1.rds")
p2 <- readRDS("processedData/transmission_matrix_county2.rds")
p3 <- readRDS("processedData/transmission_matrix_county3.rds")
p4 <- readRDS("processedData/transmission_matrix_county4.rds")
p5 <- readRDS("processedData/transmission_matrix_county5.rds")
p6 <- readRDS("processedData/transmission_matrix_county6.rds")
p7 <- readRDS("processedData/transmission_matrix_county7.rds")


# Standardize county names
county_names <- toupper(map_data_county$County)

# Add transmission values
map_data_county$method1 <- p1["GAINES", match(county_names, colnames(p1))]
map_data_county$method2 <- p2["GAINES", match(county_names, colnames(p2))]
map_data_county$method3 <- p3["GAINES", match(county_names, colnames(p3))]
map_data_county$method4 <- p4["GAINES", match(county_names, colnames(p4))]
map_data_county$method5 <- p5["GAINES", match(county_names, colnames(p5))]
map_data_county$method6 <- p6["GAINES", match(county_names, colnames(p6))]
map_data_county$method7 <- p7["GAINES", match(county_names, colnames(p7))]


gaines_geom <- map_data_county %>% filter(toupper(County) == "GAINES")

# Create individual plots

g <- plot_district_map(map_data_county, 
                  "method7", 
                  "Probability of Large Outbreak Initiation from Gains County",
                  file_name="",
                  legend_title = "Outbreak Probability",
                  reverse_colors = F) +
  geom_sf(data = gaines_geom, fill = "blue", color = "black", size = 0.3)



ggsave("Figures/transmission_probability_county_gaines_method7_sum1-12.png", 
       g,
       width = 12,
       height = 8,
       dpi = 400)







plot_transmission_maps <- function(map_data, method_number, counties, out_dir = "Figures/county_transmission/") {
  base_colors <- c("#1a9850", "#91cf60", "#d9ef8b", "#fee08b", "#fc8d59", "#d73027")
  method <- paste0("method", method_number)
  county_names <- toupper(map_data$County)
  
  plots <- list()
  
  for (county in counties) {
    county_upper <- toupper(county)
    county_title <- str_to_title(tolower(county))
    
    # Read matrix inside the loop based on method_number
    file_path <- paste0("processedData/transmission_matrix_county", method_number, ".rds")
    p_matrix <- readRDS(file_path)
    
    # Add dynamic column to map_data
    map_data[[method]] <- p_matrix[county_upper, match(county_names, colnames(p_matrix))]
    
    # Highlight the source county
    highlighted_geom <- map_data %>% filter(toupper(County) == county_upper)
    
    # Plot
    p <- ggplot(map_data) +
      geom_sf(aes_string(fill = method), color = "gray40", size = 0.1) +
      scale_fill_gradientn(colors = base_colors, na.value = "gray80", limits = c(0, 1)) +
      labs(
        title = paste("Probability of Large Outbreak Initiation from", county_title, "County"),
        fill = "Outbreak Probability"
      ) +
      geom_sf(data = highlighted_geom, fill = "blue", color = "black", size = 0.3) +
      theme_minimal() +
      theme(
        plot.title = element_text(hjust = 0.5, size = 14),
        legend.title = element_text(size = 10),
        legend.text = element_text(size = 9)
      )
    
    # Save
    file_name <- paste0(out_dir, "transmission_probability_", tolower(county), "_method", method_number, ".png")
    ggsave(file_name, p, width = 12, height = 8, dpi = 400)
    
    plots[[county]] <- p
  }
  
  return(plots)
}

low_mmr_counties <- map_data_county %>%
  arrange(MMR) %>%
  slice_head(n = 5) %>%
  pull(County) %>%
  toupper()

plot_transmission_maps(map_data_county, method_number = 7, counties = low_mmr_counties)







