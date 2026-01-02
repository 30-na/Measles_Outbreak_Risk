library(sf)
library(ggplot2)
library(reshape2)
library(lubridate)
library(dplyr)
library(tidyr)
library(readr)
library(tigris)
library(patchwork)
library(cowplot)
library(grid)
library(tidyverse)
library(grid)

#------------------ Load datasets -------------------
mmr_district <- readRDS("ProcessedData/merged_map_district_df.rds")

# Set vaccine efficacy and basic reproduction number
efficacy <- 0.97
R0 <- 15

# Infection Proportion Functions ----
# Interface (Lin. )
find_internal_infection_lin <- function(v, R0_value = R0, tol = 1e-4, efficacy_rate=efficacy) {
  if (is.na(v)) return(NA_real_)
  a <- 1 - efficacy_rate  * v
  X_vals <- seq(0, 1, length.out = 10000)
  f_X <- X_vals - a * (1 - exp(-R0_value * X_vals))
  close_to_zero <- abs(f_X) < tol
  return(max(X_vals[close_to_zero]))
}


# Calculate Probability of infection for District ----

infection_district <- mmr_district %>%
  rename(
    County = county,
    MMR = mmr,
    MMR1 = mmr1,
    MMR2 = mmr2,
    MMR3 = mmr3,
    MMR4 = mmr4
  ) %>%
  mutate(

    #  p^I attack rate 
    # Fraction of students infected if an outbreak occurs at 5-17
    under18_infection_prob = sapply(MMR, find_internal_infection_lin),
    under18_infection_prob1 = sapply(MMR1, find_internal_infection_lin),
    under18_infection_prob2 = sapply(MMR2, find_internal_infection_lin),
    under18_infection_prob3 = sapply(MMR3, find_internal_infection_lin),
    under18_infection_prob4 = sapply(MMR4, find_internal_infection_lin),
    
    
    #P^(LO) = 1 − [1 / ((1−eV)R0)] outbreak probability
    # Probability an outbreak takes off
    outbreak_prob = pmax(0, 1 - (1 / ((1 - efficacy * MMR) * R0))), 
    outbreak_prob1 = pmax(0, 1 - (1 / ((1 - efficacy * MMR1) * R0))),
    outbreak_prob2 = pmax(0, 1 - (1 / ((1 - efficacy * MMR2) * R0))),
    outbreak_prob3 = pmax(0, 1 - (1 / ((1 - efficacy * MMR3) * R0))),
    outbreak_prob4 = pmax(0, 1 - (1 / ((1 - efficacy * MMR4) * R0))),
    
    # Outbreak size
    #Expected number of infected students
    expected_outbreak_size = log10(outbreak_prob * (under18_infection_prob * enrollment)),
    expected_outbreak_size1 = log10(outbreak_prob1 * (under18_infection_prob1 * enrollment)),
    expected_outbreak_size2 = log10(outbreak_prob2 * (under18_infection_prob2 * enrollment)),
    expected_outbreak_size3 = log10(outbreak_prob3 * (under18_infection_prob3 * enrollment)),
    expected_outbreak_size4 = log10(outbreak_prob4 * (under18_infection_prob4 * enrollment)),
    
    # Expected attach rate p^I * p^LO 
    # expected attack rate = attack rate * outbreak probability
    expected_attack_rate = outbreak_prob * under18_infection_prob,
    expected_attack_rate1 = outbreak_prob1 * under18_infection_prob1,
    expected_attack_rate2 = outbreak_prob2 * under18_infection_prob2,
    expected_attack_rate3 = outbreak_prob3 * under18_infection_prob3,
    expected_attack_rate4 = outbreak_prob4 * under18_infection_prob4
  )


saveRDS(infection_district, "ProcessedData/map_district_infection_proportion_v02.rds")

# Visualization ------
# Load County boundries
tx_counties <- counties(
  state = "TX",
  year = 2024
) %>%
  mutate(
    county = toupper(NAME)
  )




#########################


p_expected_size <- ggplot(infection_district) +
  geom_sf(
    aes(fill = expected_outbreak_size),
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
    colors = c("#41ae76", "#fee08b", "#d73027"),
    na.value = "lightgray",
    limits = c(-3, 4)
  ) +
  theme_minimal() +
  labs(
    title = "Expected size of an unabated outbreak",
    fill = expression(log[10]*"(persons)")
  )

p_expected_attackrate <- ggplot(infection_district) +
  geom_sf(
    aes(fill = expected_attack_rate),
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
    colors = c("#41ae76", "#fee08b", "#d73027"),
    na.value = "lightgray",
    limits = c(0, 0.4)
  ) +
  theme_minimal() +
  labs(
    title = "Expected attack rate",
    fill = ""
  )

p_expected_size <- p_expected_size +
  labs(title = "B) Expected size of an unabated outbreak")+
  theme(
    plot.title = element_text(size = 20, face = "bold", hjust = 0)
  )

p_expected_attackrate <- p_expected_attackrate +
  labs(title = "A) Expected attack rate of an unabated outbreak")+
  theme(
    plot.title = element_text(size = 20, face = "bold", hjust = 0)
  )

p_combined <- (p_expected_attackrate + p_expected_size) +
  plot_annotation(
    title = "",
    theme = theme(
      plot.title = element_text(size = 22, face = "bold", hjust = 0.5)
    )
  )

ggsave(
  "Figures/expected_outbreak_combined_R15_senario00.png",
  plot = p_combined,
  width = 16,
  height = 8,
  dpi = 300
)



############################ change in attack rate

infection_district_scenario <- infection_district %>%
  mutate(
    perc_diff_attack_rate1 = expected_attack_rate1 - expected_attack_rate,
    perc_diff_attack_rate2 = expected_attack_rate2 - expected_attack_rate,
    perc_diff_attack_rate3 = expected_attack_rate3 - expected_attack_rate,
    perc_diff_attack_rate4 = expected_attack_rate4 - expected_attack_rate
  )

vals <- c(
  infection_district_scenario$perc_diff_attack_rate1,
  infection_district_scenario$perc_diff_attack_rate2,
  infection_district_scenario$perc_diff_attack_rate3,
  infection_district_scenario$perc_diff_attack_rate4
)

global_min <- min(vals, na.rm = TRUE)
global_max <- max(vals, na.rm = TRUE)

limits_global <- c(global_min, global_max)

max_abs <- max(abs(limits_global))
limits_sym <- c(-max_abs, max_abs)
fill_scale <- scale_fill_gradientn(
  colours = c('#006837','#1a9850',"#66bd63", '#ffffbf', "#f46d43",'#d73027','#a50026'),
  # values = scales::rescale(c(-max_abs, -max_abs/2, 0, max_abs/2, max_abs), from = limits_sym),
  limits = limits_sym,
  name   = "Change in attack rate",
  guide = guide_colorbar(
    direction  = "horizontal",
    barwidth   = unit(6, "cm"),
    barheight  = unit(0.4, "cm"),
    title.position = "top"
  )
)

# -------------------------------
# Common map theme
# -------------------------------
map_theme <- theme(
  axis.text        = element_blank(),
  axis.ticks       = element_blank(),
  axis.title       = element_blank(),
  panel.grid       = element_blank(),
  panel.background = element_rect(fill = "white", color = NA),
  plot.background  = element_rect(fill = "white", color = NA),
  legend.position  = c(0.5, 0.05),   # inside bottom-center
  legend.justification = c(0.5, 0),
  legend.background = element_rect(fill = "white", color = NA)
)

# -------------------------------
# Plots
# -------------------------------
p1 <- ggplot(infection_district_scenario) +
  geom_sf(aes(fill = perc_diff_attack_rate1), color = "gray40") +
  geom_sf(data = tx_counties, fill = NA, color = "gray20", size = 0.4) +
  fill_scale +
  labs(title = "A) Scenario 1") +
  coord_sf(datum = NA) +
  map_theme +
  theme(legend.position = "none")

p2 <- ggplot(infection_district_scenario) +
  geom_sf(aes(fill = perc_diff_attack_rate2), color = "gray40") +
  geom_sf(data = tx_counties, fill = NA, color = "gray20", size = 0.4) +
  fill_scale +
  labs(title = "B) Scenario 2") +
  coord_sf(datum = NA) +
  map_theme +
  theme(legend.position = "none")

p3 <- ggplot(infection_district_scenario) +
  geom_sf(aes(fill = perc_diff_attack_rate3), color = "gray40") +
  geom_sf(data = tx_counties, fill = NA, color = "gray20", size = 0.4) +
  fill_scale +
  labs(title = "C) Scenario 3") +
  coord_sf(datum = NA) +
  map_theme +
  theme(legend.position = "none")

p4 <- ggplot(infection_district_scenario) +
  geom_sf(aes(fill = perc_diff_attack_rate4), color = "gray40") +
  geom_sf(data = tx_counties, fill = NA, color = "gray20", size = 0.4) +
  fill_scale +
  labs(title = "D) Scenario 4") +
  coord_sf(datum = NA) +
  map_theme +
  theme(legend.position = "none")




# -------------------------------
# Combine and save
# -------------------------------
legend_plot <- ggplot(infection_district_scenario) +
  geom_sf(aes(fill = perc_diff_attack_rate4)) +
  fill_scale +
  coord_sf(datum = NA) +
  
  theme(
    legend.position = "bottom",
    legend.background = element_blank()
  ) +
  guides(
    fill = guide_colorbar(
      direction = "horizontal",
      barwidth  = unit(10, "cm"),
      barheight = unit(0.6, "cm"),
      title.position = "top"
    )
  )




g <- ggplotGrob(legend_plot)

legend <- g$grobs[[which(sapply(g$grobs, function(x) x$name) == "guide-box")]]



final_plot <- plot_grid(
  plot_grid(p1, p2, p3, p4, nrow = 2),
  legend,
  ncol = 1,
  rel_heights = c(1, 0.12)
)

ggsave(
  "Figures/change_attackrate_R15.png",
  final_plot,
  width = 10,
  height = 10,
  dpi = 400,
  bg = "white"
)

######################## Attach Rate ####


vals <- c(
  infection_district_scenario$expected_attack_rate1,
  infection_district_scenario$expected_attack_rate2,
  infection_district_scenario$expected_attack_rate3,
  infection_district_scenario$expected_attack_rate4
)

global_min <- min(vals, na.rm = TRUE)
global_max <- max(vals, na.rm = TRUE)

limits_global <- c(global_min, global_max)

fill_scale <- scale_fill_gradientn(
  colours = c("#66bd63", '#ffffbf', "#f46d43",'#d73027','#a50026'),
  limits = limits_global,
  name   = "Attack rate",
  guide = guide_colorbar(
    direction  = "horizontal",
    barwidth   = unit(6, "cm"),
    barheight  = unit(0.4, "cm"),
    title.position = "top"
  )
)

# -------------------------------
# Common map theme
# -------------------------------
map_theme <- theme(
  axis.text        = element_blank(),
  axis.ticks       = element_blank(),
  axis.title       = element_blank(),
  panel.grid       = element_blank(),
  panel.background = element_rect(fill = "white", color = NA),
  plot.background  = element_rect(fill = "white", color = NA),
  legend.position  = c(0.5, 0.05),   # inside bottom-center
  legend.justification = c(0.5, 0),
  legend.background = element_rect(fill = "white", color = NA)
)

# -------------------------------
# Plots
# -------------------------------
p1 <- ggplot(infection_district_scenario) +
  geom_sf(aes(fill = expected_attack_rate1), color = "gray40") +
  geom_sf(data = tx_counties, fill = NA, color = "gray20", size = 0.4) +
  fill_scale +
  labs(title = "A) Scenario 1") +
  coord_sf(datum = NA) +
  map_theme +
  theme(legend.position = "none")

p2 <- ggplot(infection_district_scenario) +
  geom_sf(aes(fill = expected_attack_rate2), color = "gray40") +
  geom_sf(data = tx_counties, fill = NA, color = "gray20", size = 0.4) +
  fill_scale +
  labs(title = "B) Scenario 2") +
  coord_sf(datum = NA) +
  map_theme +
  theme(legend.position = "none")

p3 <- ggplot(infection_district_scenario) +
  geom_sf(aes(fill = expected_attack_rate3), color = "gray40") +
  geom_sf(data = tx_counties, fill = NA, color = "gray20", size = 0.4) +
  fill_scale +
  labs(title = "C) Scenario 3") +
  coord_sf(datum = NA) +
  map_theme +
  theme(legend.position = "none")

p4 <- ggplot(infection_district_scenario) +
  geom_sf(aes(fill = expected_attack_rate4), color = "gray40") +
  geom_sf(data = tx_counties, fill = NA, color = "gray20", size = 0.4) +
  fill_scale +
  labs(title = "D) Scenario 4") +
  coord_sf(datum = NA) +
  map_theme +
  theme(legend.position = "none")




# -------------------------------
# Combine and save
# -------------------------------
legend_plot <- ggplot(infection_district_scenario) +
  geom_sf(aes(fill = expected_attack_rate4)) +
  fill_scale +
  coord_sf(datum = NA) +
  
  theme(
    legend.position = "bottom",
    legend.background = element_blank()
  ) +
  guides(
    fill = guide_colorbar(
      direction = "horizontal",
      barwidth  = unit(10, "cm"),
      barheight = unit(0.6, "cm"),
      title.position = "top"
    )
  )




g <- ggplotGrob(legend_plot)

legend <- g$grobs[[which(sapply(g$grobs, function(x) x$name) == "guide-box")]]



final_plot <- plot_grid(
  plot_grid(p1, p2, p3, p4, nrow = 2),
  legend,
  ncol = 1,
  rel_heights = c(1, 0.12)
)

ggsave(
  "Figures/attackrate_R15.png",
  final_plot,
  width = 10,
  height = 10,
  dpi = 400,
  bg = "white"
)



######################## Expected Outbreak size ####

# For visualization
infection_district_scenario <- infection_district_scenario %>%
  mutate(
    expected_outbreak_size = ifelse(is.infinite(expected_outbreak_size),
                                     -4, expected_outbreak_size),
    expected_outbreak_size1 = ifelse(is.infinite(expected_outbreak_size1),
                                     -4, expected_outbreak_size1),
    expected_outbreak_size2 = ifelse(is.infinite(expected_outbreak_size2),
                                     -4, expected_outbreak_size2),
    expected_outbreak_size3 = ifelse(is.infinite(expected_outbreak_size3),
                                     -4, expected_outbreak_size3),
    expected_outbreak_size4 = ifelse(is.infinite(expected_outbreak_size4),
                                     -4, expected_outbreak_size4)
  )

vals <- c(
  infection_district_scenario$expected_outbreak_size1,
  infection_district_scenario$expected_outbreak_size2,
  infection_district_scenario$expected_outbreak_size3,
  infection_district_scenario$expected_outbreak_size4
)

global_min <- min(vals, na.rm = TRUE)
global_max <- max(vals, na.rm = TRUE)

limits_global <- c(global_min, global_max)

fill_scale <- scale_fill_gradientn(
  colours = c("#66bd63", '#ffffbf', "#f46d43",'#a50026'),
  limits = limits_global,
  name   = expression(log[10]*"(persons)"),
  guide = guide_colorbar(
    direction  = "horizontal",
    barwidth   = unit(6, "cm"),
    barheight  = unit(0.4, "cm"),
    title.position = "top"
  )
)

# -------------------------------
# Common map theme
# -------------------------------
map_theme <- theme(
  axis.text        = element_blank(),
  axis.ticks       = element_blank(),
  axis.title       = element_blank(),
  panel.grid       = element_blank(),
  panel.background = element_rect(fill = "white", color = NA),
  plot.background  = element_rect(fill = "white", color = NA),
  legend.position  = c(0.5, 0.05),   # inside bottom-center
  legend.justification = c(0.5, 0),
  legend.background = element_rect(fill = "white", color = NA)
)

# -------------------------------
# Plots
# -------------------------------
p1 <- ggplot(infection_district_scenario) +
  geom_sf(aes(fill = expected_outbreak_size1), color = "gray40") +
  geom_sf(data = tx_counties, fill = NA, color = "gray20", size = 0.4) +
  fill_scale +
  labs(title = "A) Scenario 1") +
  coord_sf(datum = NA) +
  map_theme +
  theme(legend.position = "none")

p2 <- ggplot(infection_district_scenario) +
  geom_sf(aes(fill = expected_outbreak_size2), color = "gray40") +
  geom_sf(data = tx_counties, fill = NA, color = "gray20", size = 0.4) +
  fill_scale +
  labs(title = "B) Scenario 2") +
  coord_sf(datum = NA) +
  map_theme +
  theme(legend.position = "none")

p3 <- ggplot(infection_district_scenario) +
  geom_sf(aes(fill = expected_outbreak_size3), color = "gray40") +
  geom_sf(data = tx_counties, fill = NA, color = "gray20", size = 0.4) +
  fill_scale +
  labs(title = "C) Scenario 3") +
  coord_sf(datum = NA) +
  map_theme +
  theme(legend.position = "none")

p4 <- ggplot(infection_district_scenario) +
  geom_sf(aes(fill = expected_outbreak_size4), color = "gray40") +
  geom_sf(data = tx_counties, fill = NA, color = "gray20", size = 0.4) +
  fill_scale +
  labs(title = "D) Scenario 4") +
  coord_sf(datum = NA) +
  map_theme +
  theme(legend.position = "none")




# -------------------------------
# Combine and save
# -------------------------------
legend_plot <- ggplot(infection_district_scenario) +
  geom_sf(aes(fill = expected_outbreak_size4)) +
  fill_scale +
  coord_sf(datum = NA) +
  
  theme(
    legend.position = "bottom",
    legend.background = element_blank()
  ) +
  guides(
    fill = guide_colorbar(
      direction = "horizontal",
      barwidth  = unit(10, "cm"),
      barheight = unit(0.6, "cm"),
      title.position = "top"
    )
  )




g <- ggplotGrob(legend_plot)

legend <- g$grobs[[which(sapply(g$grobs, function(x) x$name) == "guide-box")]]



final_plot <- plot_grid(
  plot_grid(p1, p2, p3, p4, nrow = 2),
  legend,
  ncol = 1,
  rel_heights = c(1, 0.12)
)

ggsave(
  "Figures/expected_outbreak_size_R15.png",
  final_plot,
  width = 10,
  height = 10,
  dpi = 400,
  bg = "white"
)
