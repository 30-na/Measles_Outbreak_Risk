library(readr)
library(ggplot2)
library(dplyr)
library(tidyr)
library(tigris)
library(ggpubr)
library(readxl)
library(caret)    
library(broom)  
library(pROC)
library(pROC)
library(ggplot2)

map_probability <- readRDS("ProcessedData/map_probability.rds")

west_texas_jun17 <- data.frame(
  # Report_Date = rep("2025-06-17", 35),
  county = c(
    "Andrews", "Atascosa", "Bailey", "Borden", "Brewster", "Brown", "Carson",
    "Cochran", "Collin", "Dallam", "Dawson", "Eastland", "Ector", "El Paso",
    "Erath", "Gaines", "Garza", "Hale", "Hardeman", "Hockley", "Lamar", "Lamb",
    "Lubbock", "Lynn", "Martin", "McLennan", "Midland", "Parmer", "Potter",
    "Randall", "Reeves", "Rockwall", "Terry", "Upshur", "Yoakum"
  ),
  Confirmed = c(
    3, 1, 2, 1, 1, 1, 1,
    14, 1, 7, 26, 2, 12, 58,
    1, 413, 2, 5, 1, 6, 23, 1,
    53, 2, 3, 9, 6, 5, 1,
    1, 2, 1, 60, 5, 20
  ),
  Percent = c(
    0.0040, 0.0013, 0.0027, 0.0013, 0.0013, 0.0013, 0.0013,
    0.0187, 0.0013, 0.0093, 0.0347, 0.0027, 0.0160, 0.0773,
    0.0013, 0.5507, 0.0027, 0.0067, 0.0013, 0.0080, 0.0307, 0.0013,
    0.0707, 0.0027, 0.0040, 0.0120, 0.0080, 0.0067, 0.0013,
    0.0013, 0.0027, 0.0013, 0.0800, 0.0067, 0.0267
  )
) %>%
  mutate(
    county = toupper(county)
  )

efficacy <- 0.97

# Get the population
# county_pop <- read_csv("rawData/Counties_Pop.txt") %>%
#   mutate(
#     FENAME = case_when(
#       FENAME == "DE WITT" ~ "DEWITT",
#       TRUE ~ FENAME)
#   )


# Get the map and merge with population
# tx_counties_map <- counties(state = "TX", year = 2024) %>%
#   mutate(
#     NAME = toupper(NAME)
#   ) %>%
#   select(
#     "county" = NAME
#   ) %>%
#   left_join(
#     county_pop,
#     by = c("county" = "FENAME")
#   ) %>%
#   select(
#     county,
#     "population" = total
#   )


# Merge  with cases
# Merge cases

merged_df <- map_probability %>%
  left_join(
    west_texas_jun17,
    c("County" = "county")
    ) %>%
  select(
    "county" = County,
    MMR,
    "population" = total,
    Confirmed,
    Percent
  )


########## get the indirect outbreak risk



compute_indirect_risk <- function(trans_mat, county_name = "GAINES", threshold = .5) {
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


get_indirect_risk_list <- function(strategy=0, county_name = "GAINES", method = 7, threshold = 0.5) {
  # Load transmission matrix
  file_path <- paste0("ProcessedData/pij_M", method, "_S", strategy, ".rds")
  trans_mat <- readRDS(file_path)
  
  # Compute indirect risk
  risk_vec <- compute_indirect_risk(trans_mat, county_name = county_name, threshold = threshold)
  
  # Return as sorted data.frame
  out_df <- data.frame(
    county = names(risk_vec),
    indirect_risk = round(risk_vec, 5)
  ) %>%
    arrange(desc(indirect_risk))
  
  return(out_df)
}



indirect_df <- get_indirect_risk_list(strategy = 0, county_name = "Gaines")
 

# Merge with the NA-filled data

merged_with_indirect_na <- merged_df %>%
  left_join(
    indirect_df,
    by = "county"
    ) %>%
  filter(
    !is.na(indirect_risk)
    ) %>%
  mutate(
    Confirmed = replace_na(Confirmed, 0),
    Percent = replace_na(Percent, 0),
    cases_per_capita = Confirmed / population,
    effective_pop = population * (1 - efficacy * MMR),
    cases_per_effpop = Confirmed / effective_pop
  )

write.csv(sf::st_drop_geometry(merged_with_indirect_na), "ProcessedData/merged_with_indirect_na.csv", row.names = FALSE)
################# Corelation aanalysis 

# 1. All counties (even those with Confirmed = 0)
df_all <- merged_with_indirect_na

# 2. Only counties with Confirmed > 0
df_nonzero <- df_all %>% filter(Confirmed > 0)

# 3. Only counties with Confirmed ≥ 3
df_gt3 <- df_all %>% filter(Confirmed >= 3)

##############PLOTS

save_correlation_plot <- function(data, x_col, y_col, title = NULL, filename = "correlation_plot.png") {
  # Create the ggplot
  p <- ggplot(data, aes_string(x = x_col, y = y_col)) +
    geom_point(
      alpha = 0.6
      ) +
    geom_smooth(
      method = "lm", 
      se = TRUE, 
      color = "blue"
      ) +
    stat_cor(method = "pearson",
             label.x = 0, 
             label.y = 0.9) +
    labs(
      title = ifelse(is.null(title), paste0(x_col, " vs ", y_col), title),
      x = x_col,
      y = y_col
    ) +
    theme_minimal()
  
  # Save the plot
  ggsave(filename, plot = p, width = 6, height = 5, dpi = 300)
}


# Plot for Confirmed cases vs indirect risk
save_correlation_plot(df_all,     "cases_per_effpop", "indirect_risk", "All Counties",     "Figures/corr_cases_effpop_all.png")
save_correlation_plot(df_nonzero, "cases_per_effpop", "indirect_risk", "Confirmed > 0",    "Figures/corr_cases_effpop_nonzero.png")
save_correlation_plot(df_gt3,     "cases_per_effpop", "indirect_risk", "Confirmed ≥ 3",    "Figures/corr_cases_effpop_gt3.png")



########## Logistic Regresion
# 1. Define binary outcome: 1 if Confirmed ≥ 6, else 0
df_bin <- merged_with_indirect_na %>%
  mutate(
    high_case = ifelse(Confirmed >= 6, 1, 0)
  )

# 2. Fit logistic regression model
# You can include more predictors if desired
logit_model <- glm(high_case ~ cases_per_effpop + indirect_risk, 
                   data = df_bin,
                   family = binomial)

# 3. Predict probabilities
df_bin <- df_bin %>%
  mutate(
    prob = predict(logit_model, type = "response"),
    pred = ifelse(prob >= 0.1, 1, 0)  # threshold can be adjusted
  )

# 4. Evaluate using confusion matrix
confusion <- confusionMatrix(
  factor(df_bin$pred),
  factor(df_bin$high_case),
  positive = "1"
)

# 5. Print results
print(confusion)
summary(logit_model)


# Compute ROC
roc_obj <- roc(df_bin$high_case, df_bin$prob)

# Get full ROC coordinates with thresholds
# Get all coordinates as a data frame
roc_coords <- coords(
  roc_obj,
  x = "all",
  ret = c("threshold", "specificity", "sensitivity"),
  transpose = FALSE
) %>%
  as.data.frame() %>%
  mutate(
    fpr = 1 - specificity
  )

# Choose ~10 evenly spaced thresholds to label
label_idx <- round(seq(1, nrow(roc_coords), length.out = 20))

# Plot
p <- ggplot(roc_coords, aes(x = fpr, y = sensitivity)) +
  geom_line(color = "blue", size = 1.2) +
  geom_abline(linetype = "dashed", color = "gray") +
  geom_text(data = roc_coords[label_idx, ],
            aes(label = round(threshold, 2)),
            vjust = -0.7, hjust = 1.1, size = 3, color = "red") +
  annotate("text", x = 0.6, y = 0.1,
           label = paste0("AUC = ", round(auc(roc_obj), 3)),
           color = "blue", size = 5) +
  labs(
    title = "ROC Curve with Thresholds Annotated",
    x = "False Positive Rate (1 - Specificity)",
    y = "True Positive Rate (Sensitivity)"
  ) +
  xlim(0, 1) +
  ylim(0, 1) +
  theme_minimal(base_size = 14)

# Save it
ggsave("Figures/roc_curve_threshold_labels.png", plot = p, width = 6, height = 5, dpi = 300)




















compute_correlation_metrics <- function(data, x_col = "Confirmed", y_col = "indirect_risk", label = "df") {
  results <- list()
  
  # Pearson
  pear <- cor.test(data[[x_col]], data[[y_col]], method = "pearson")
  results[["Pearson"]] <- c(r = round(pear$estimate, 3), p = signif(pear$p.value, 3))
  
  # Spearman
  spear <- cor.test(data[[x_col]], data[[y_col]], method = "spearman")
  results[["Spearman"]] <- c(r = round(spear$estimate, 3), p = signif(spear$p.value, 3))
  
  # Kendall
  kend <- cor.test(data[[x_col]], data[[y_col]], method = "kendall")
  results[["Kendall"]] <- c(r = round(kend$estimate, 3), p = signif(kend$p.value, 3))
  
  # Logistic (modeling outbreak probability as a function of cases)
  # Convert indirect_risk to binary: high vs low (threshold = 0.5)
  log_df <- data %>%
    mutate(binary_risk = ifelse(indirect_risk >= 0.5, 1, 0))
  
  log_fit <- glm(binary_risk ~ get(x_col), data = log_df, family = binomial())
  log_p <- summary(log_fit)$coefficients[2, 4]
  log_odds <- exp(coef(log_fit)[2])  # odds ratio
  
  results[["Logistic"]] <- c("odds_ratio" = round(log_odds, 3), "p" = signif(log_p, 3))
  
  # Combine
  out <- do.call(rbind, results)
  out <- data.frame(type = rownames(out), out, row.names = NULL)
  out$dataset <- label
  return(out)
}

# Run for all three
res_all <- compute_correlation_metrics(df_all, "Confirmed", "indirect_risk", label = "All")
res_nonzero <- compute_correlation_metrics(df_nonzero, "Confirmed", "indirect_risk", label = "Confirmed > 0")
res_gt3 <- compute_correlation_metrics(df_gt3, "Confirmed", "indirect_risk", label = "Confirmed ≥ 3")

# Combine and view
correlation_summary <- bind_rows(res_all, res_nonzero, res_gt3)
print(correlation_summary)
