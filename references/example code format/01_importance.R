library(tidyverse)
library(stringr)
library(Matrix)     
library(sandwich)   
library(lmtest)     
library(glmnet)   
library(fixest)
library(ggplot2)



source("/Users/liuqianyi/Dropbox/swb/Epinion/code/01_importance_functions.R")



# ============================================
# PART 1: Read data and estimate MUs
# ============================================

# Read the data for MU estimation
df_mu <- read.csv("/Users/liuqianyi/Dropbox/swb/Epinion/raw/SWB wave 2 (Nov 4)_November 10, 2025.csv")

# Process the data
processed_data <- purrr::map_df(1:nrow(df_mu), function(i) {
  if (is.na(df_mu$allComparisons[i]) || df_mu$allComparisons[i] == "") return(NULL)
  ac <- parse_comparison(df_mu$allComparisons[i], df_mu$ResponseId[i])
  if (is.null(ac) || nrow(ac) == 0) return(NULL)
  ac %>%
    group_by(respondent_id, aspect_j, aspect_k, delta_j, delta_k) %>%
    slice_tail(n = 1) %>%
    ungroup()
}) %>% filter(!is.na(choice))

# Filter to only include core aspects
processed_data_core <- processed_data %>%
  filter(aspect_j %in% core_aspects & aspect_k %in% core_aspects)

# Estimate the model for core aspects
estimate <- estimate_model_glm(processed_data_core)

# Extract MU estimates for core aspects
baseline_MU <- estimate$MU[core_aspects]
baseline_SE <- estimate$SE[core_aspects]

# ============================================
# PART 2: Extract Importance Ratings
# ============================================

# Function to extract importance data
extract_importance_data <- function(df) {
  # Get the shuffled statements to match with importance ratings
  shuffled_cols <- paste0("shuffled_statement_", 1:16)
  importance_cols <- paste0("core_importance_", 1:16)
  
  # Initialize results list
  importance_data <- list()
  
  for (i in 1:nrow(df)) {
    respondent_id <- df$ResponseId[i]
    
    # Extract shuffled statements for this respondent
    shuffled_statements <- df[i, shuffled_cols]
    
    # Extract importance ratings for this respondent
    importance_ratings <- df[i, importance_cols]
    
    # Create mapping between aspects and importance ratings
    for (j in 1:16) {
      aspect_text <- as.character(shuffled_statements[1, j])
      importance_value <- suppressWarnings(as.numeric(importance_ratings[1, j]))
      
      # Only include if both aspect and rating are valid
      if (!is.na(aspect_text) && aspect_text != "" && 
          !is.na(importance_value) && aspect_text %in% core_aspects) {
        
        importance_data[[length(importance_data) + 1]] <- tibble(
          respondent_id = respondent_id,
          aspect = aspect_text,
          importance_rating = importance_value
        )
      }
    }
  }
  
  # Combine all data
  if (length(importance_data) > 0) {
    bind_rows(importance_data)
  } else {
    NULL
  }
}

# Read the data for importance ratings
df_importance <- read.csv("/Users/liuqianyi/Dropbox/swb/Epinion/raw/SWB wave 1 (Nov 4)_November 10, 2025.csv")

# Extract importance data
importance_df <- extract_importance_data(df_importance)

write_csv(importance_df, "/Users/liuqianyi/Dropbox/swb/Epinion/raw/importance.csv")

# Calculate average importance by aspect
avg_importance <- importance_df %>%
  group_by(aspect) %>%
  summarise(
    mean_importance = mean(importance_rating, na.rm = TRUE),
    sd_importance = sd(importance_rating, na.rm = TRUE),
    n_ratings = n(),
    se_importance = sd_importance / sqrt(n_ratings)
  ) %>%
  arrange(desc(mean_importance))

# ============================================
# PART 3: Merge and Analyze Relationship
# ============================================

# Merge importance ratings with MU estimates
comparison_df <- avg_importance %>%
  mutate(
    log_MU = baseline_MU[aspect],
    MU = exp(log_MU),
    SE_MU = baseline_SE[aspect]
  ) %>%
  filter(!is.na(log_MU))

# Calculate correlation between importance and log MU
cor_importance_logMU <- cor(comparison_df$mean_importance, comparison_df$log_MU, use = "complete.obs")

# Fit linear model
lm_model <- lm(log_MU ~ mean_importance, data = comparison_df)
cat("\nLinear regression: log_MU ~ importance_rating\n")
print(summary(lm_model))

# Standardize importance ratings
comparison_df <- comparison_df %>%
  mutate(
    importance_z = (mean_importance - mean(mean_importance)) / sd(mean_importance)
  )

# ============================================
# PART 4: Create Visualization
# ============================================

# Create scatterplot
p_importance <- ggplot(comparison_df, aes(x = importance_z, y = log_MU)) +
  geom_point(size = 3, alpha = 0.7) +
  geom_smooth(method = "lm", se = TRUE, alpha = 0.2) +
  geom_text(aes(label = substr(aspect, 1, 30)), 
            vjust = -0.5, hjust = 0, size = 2.5, 
            check_overlap = TRUE) +
  labs(
    title = "Importance Ratings vs. Log Marginal Utility",
    subtitle = paste0("Correlation = ", round(cor_importance_logMU, 3)),
    x = "Standardized Importance Rating",
    y = "Log Marginal Utility (from tradeoffs)"
  ) +
  theme_minimal(base_size = 12) +
  theme(
    plot.title = element_text(size = 14, face = "bold"),
    plot.subtitle = element_text(size = 11),
    axis.title = element_text(size = 11)
  )

print(p_importance)

# Save the plot
ggsave("/Users/liuqianyi/Dropbox/swb/Epinion/graph/importance_vs_logMU.png", p_importance, width = 10, height = 8, dpi = 300)



results <- importance_predictive_power_fast("/Users/liuqianyi/Dropbox/swb/Epinion/raw/choice_level_data_20251111.csv")

lollipop_above_below <- results$results %>%
  mutate(aspect_label = paste0(substr(aspect, 1, 25), "\n(", round(predictive_ratio, 3), ")")) %>%
  ggplot(aes(x = reorder(aspect_label, predictive_ratio))) +
  geom_segment(aes(xend = aspect_label, y = log_MU_low, yend = log_MU_high), 
               arrow = arrow(length = unit(0.1, "cm"))) +
  geom_point(aes(y = log_MU_low), color = "red", size = 2) +
  geom_point(aes(y = log_MU_high), color = "blue", size = 2) +
  coord_flip() +
  labs(y = "Log MU", x = "Aspect (predictive ratio)", 
       title = "Log MU: Below-median (red) to above-median (Blue) importance")

ggsave("/Users/liuqianyi/Dropbox/swb/Epinion/graph/importance_vs_logMU_above_below.png", lollipop_above_below, width = 10, height = 8, dpi = 300)
