library(tidyverse)
library(sandwich)
library(lmtest)
library(ggplot2)
library(cowplot)   # for arranging 4x4 pages

# avoid scientific notation globally for printed numbers
options(scipen = 999)

# helper: compute clustered SE of a linear combination of model coefficients
# usage: lincomb_se(model, c("term1","term2"), weights = c(1,1))
lincomb_se <- function(model, terms, weights = NULL) {
  vc <- vcovCL(model, cluster = ~ respondent_id)
  coefn <- names(coef(model))
  L <- rep(0, length(coefn)); names(L) <- coefn
  if (is.null(weights)) weights <- rep(1, length(terms))
  for (i in seq_along(terms)) {
    if (terms[i] %in% coefn) L[terms[i]] <- weights[i]
  }
  sqrt(as.numeric(t(L) %*% vc %*% L))
}

# Notation mapping for plots/text:
# j -> 1 (top), k -> 2 (bottom). delta_j -> Δ1, delta_k -> Δ2.
# choice -> D, Pr(D = 1) = Pr(choose top)

# Create output folder if it doesn't exist
output_dir <- "/Users/liuqianyi/Dropbox/swb/Epinion/graph"
if (!dir.exists(output_dir)) dir.create(output_dir, recursive = TRUE)

# Load data - UPDATE THIS PATH to use the merged file
data <- read.csv("/Users/liuqianyi/Dropbox/swb/Epinion/raw/choice_level_data_20251111.csv", stringsAsFactors = FALSE)

# Create above-median importance indicators
data$am_imp_j <- as.numeric(data$importance_j > median(data$importance_j, na.rm = TRUE))
data$am_imp_k <- as.numeric(data$importance_k > median(data$importance_k, na.rm = TRUE))

# Create above-median *level* indicators (use level_j/level_k if present, otherwise fallback to importance)
if ("level_j" %in% names(data) && "level_k" %in% names(data)) {
  data$am_level_j <- as.numeric(data$level_j > median(data$level_j, na.rm = TRUE))
  data$am_level_k <- as.numeric(data$level_k > median(data$level_k, na.rm = TRUE))
} else {
  cat("Note: 'level_j'/'level_k' not found — using 'importance_j'/'importance_k' as a fallback for above-median level splits.\n")
  data$am_level_j <- as.numeric(data$importance_j > median(data$importance_j, na.rm = TRUE))
  data$am_level_k <- as.numeric(data$importance_k > median(data$importance_k, na.rm = TRUE))
}

# Create aspect groups (factor encoding)
data$group_aj <- as.factor(data$aspect_j)
data$group_ak <- as.factor(data$aspect_k)

cat("Sample size:", nrow(data), "\n")
cat("Median importance_j:", median(data$importance_j, na.rm = TRUE), "\n")
cat("Median importance_k:", median(data$importance_k, na.rm = TRUE), "\n")

# =============================================================================
# MODELS
# =============================================================================

# Model 1: No controls
model1 <- lm(choice ~ delta_j + delta_k, data = data)
coeftest(model1, vcov = vcovCL(model1, cluster = ~ respondent_id))

# Model 2: Choice by aspect (pooled with aspect interactions)
model2 <- lm(choice ~ delta_j + delta_j:group_aj + delta_k + delta_k:group_ak, data = data)

# Model 3: Choice by importance
model3 <- lm(choice ~ delta_j + delta_j:am_imp_j + delta_k, data = data)
coeftest(model3, vcov = vcovCL(model3, cluster = ~ respondent_id))

# Model 4: Choice by importance (both j and k interactions)
model4 <- lm(choice ~ delta_j + delta_j:am_imp_j + delta_k + delta_k:am_imp_k, data = data)

# Model 5: Aspect-specific importance
model5 <- lm(choice ~ delta_j + delta_j:group_aj:am_imp_j + delta_k + delta_k:group_ak, data = data)

# =============================================================================
# EXTRACT RESULTS WITH CLUSTERED SEs (helper)
# =============================================================================

get_clustered_results <- function(model, pattern = NULL) {
  vcov_clustered <- vcovCL(model, cluster = ~ respondent_id)
  coefs <- coef(model)
  ses <- sqrt(diag(vcov_clustered))
  if (!is.null(pattern)) {
    matches <- grepl(pattern, names(coefs))
    coefs <- coefs[matches]
    ses <- ses[matches]
  }
  tibble(
    term = names(coefs),
    coefficient = as.numeric(coefs),
    se = as.numeric(ses),
    ci_lower = coefficient - 1.96 * se,
    ci_upper = coefficient + 1.96 * se,
    significant = abs(coefficient / se) > 1.96
  )
}

# Model 4 effects (label with Δ1/Δ2 notation)
model4_j <- get_clustered_results(model4, "delta_j:am_imp_j") %>%
  mutate(effect_type = "Δ1 × 1(above median importance 1)")
model4_k <- get_clustered_results(model4, "delta_k:am_imp_k") %>%
  mutate(effect_type = "Δ2 × 1(above median importance 2)")
model4_effects <- bind_rows(model4_j, model4_k)

# Model 2 interactions (for main effects)
model2_int_j <- get_clustered_results(model2, "delta_j:group_aj") %>%
  mutate(aspect = str_extract(term, "group_aj(.+)$") %>% str_remove("^group_aj"))

# Model 2 main effects (compute linear combos' point estimates; SE via lincomb_se helper)
model2_coeftest <- coeftest(model2, vcov = vcovCL(model2, cluster = ~ respondent_id))
baseline_j <- model2_coeftest["delta_j", c("Estimate", "Std. Error")]
model2_main_j <- model2_int_j %>%
  mutate(
    main_effect = baseline_j[1] + coefficient,
    se_main = map2_dbl(term, coefficient, function(term_name, inter_coef) {
      # compute SE of delta_j + interaction_term (uses vcovCL under the hood)
      lincomb_se(model2, c("delta_j", term_name))
    }),
    ci_lower_main = main_effect - 1.96 * se_main,
    ci_upper_main = main_effect + 1.96 * se_main,
    significant_main = abs(main_effect / se_main) > 1.96
  )

# Model 5 effects
model5_effects <- get_clustered_results(model5, "delta_j:group_aj.*am_imp_j") %>%
  mutate(aspect = str_extract(term, "group_aj(.+):am_imp_j") %>% str_remove("^group_aj") %>% str_remove(":am_imp_j$"))

# =============================================================================
# Binscatter data (global)
# =============================================================================

binscatter_data_j <- data %>%
  mutate(delta_j_bin = ntile(delta_j, 20)) %>%
  group_by(delta_j_bin) %>%
  summarise(mean_delta_j = mean(delta_j), mean_choice = mean(choice), n_obs = n(), .groups = "drop")

binscatter_data_k <- data %>%
  mutate(delta_k_bin = ntile(delta_k, 20)) %>%
  group_by(delta_k_bin) %>%
  summarise(mean_delta_k = mean(delta_k), mean_choice = mean(choice), n_obs = n(), .groups = "drop")

binscatter_data_j_by_imp <- data %>%
  filter(!is.na(am_imp_j)) %>%
  group_by(am_imp_j) %>%
  mutate(delta_j_bin = ntile(delta_j, 20)) %>%
  group_by(am_imp_j, delta_j_bin) %>%
  summarise(mean_delta_j = mean(delta_j), mean_choice = mean(choice), n_obs = n(), .groups = "drop") %>%
  mutate(importance_group = ifelse(am_imp_j == 1, "Above median", "Below median"))

# Binscatter by above-median *level* (Model 3b requested)
binscatter_data_j_by_level <- data %>%
  filter(!is.na(am_level_j)) %>%
  group_by(am_level_j) %>%
  mutate(delta_j_bin = ntile(delta_j, 20)) %>%
  group_by(am_level_j, delta_j_bin) %>%
  summarise(mean_delta_j = mean(delta_j), mean_choice = mean(choice), n_obs = n(), .groups = "drop") %>%
  mutate(level_group = ifelse(am_level_j == 1, "Above median level", "Below median level"))

# =============================================================================
# PLOTTING HELPERS
# =============================================================================

# helper to format numbers: no scientific notation, fixed decimal places
fmt <- function(x, digits = 4) {
  if (is.null(x) || length(x) == 0 || all(is.na(x))) return(NA_character_)
  formatC(as.numeric(x), format = "f", digits = digits)
}

# MODEL 4 importance effects plot (with Δ1/Δ2)
plot_model4 <- function(effects) {
  effects %>%
    mutate(effect_type = fct_reorder(effect_type, coefficient)) %>%
    ggplot(aes(x = coefficient, y = effect_type)) +
    geom_vline(xintercept = 0, linetype = "dashed", color = "gray50") +
    geom_pointrange(aes(xmin = ci_lower, xmax = ci_upper, color = significant),
                    size = 0.8, fatten = 3, show.legend = FALSE) +
    scale_color_manual(values = c("FALSE" = "gray60", "TRUE" = "red")) +
    labs(title = "Model 4: Importance Effects (Δ1 / Δ2)", x = "Coefficient (Pr(D = 1))", y = "") +
    theme_minimal() +
    theme(panel.background = element_rect(fill = "white", color = NA),
          plot.background = element_rect(fill = "white", color = NA),
          panel.grid = element_blank())
}

plot_interactions <- function(effects, title) {
  effects %>%
    mutate(aspect_clean = str_trunc(aspect, 35),
           aspect_clean = fct_reorder(aspect_clean, coefficient)) %>%
    ggplot(aes(x = coefficient, y = aspect_clean)) +
    geom_vline(xintercept = 0, linetype = "dashed", color = "gray50") +
    geom_pointrange(aes(xmin = ci_lower, xmax = ci_upper, color = significant),
                    size = 0.5, fatten = 2, show.legend = FALSE) +
    scale_color_manual(values = c("FALSE" = "gray60", "TRUE" = "red")) +
    labs(title = title, x = "Coefficient (Pr(D = 1))", y = "Aspect") +
    theme_minimal(base_size = 10) +
    theme(axis.text.y = element_text(size = 8),
          panel.background = element_rect(fill = "white", color = NA),
          plot.background = element_rect(fill = "white", color = NA),
          panel.grid = element_blank())
}

plot_main_effects <- function(effects, title) {
  effects %>%
    mutate(aspect_clean = str_trunc(aspect, 35),
           aspect_clean = fct_reorder(aspect_clean, main_effect)) %>%
    ggplot(aes(x = main_effect, y = aspect_clean)) +
    geom_vline(xintercept = 0, linetype = "dashed", color = "gray50") +
    geom_pointrange(aes(xmin = ci_lower_main, xmax = ci_upper_main, color = significant_main),
                    size = 0.5, fatten = 2, show.legend = FALSE) +
    scale_color_manual(values = c("FALSE" = "gray60", "TRUE" = "red")) +
    labs(title = title, x = "Main Effect (Pr(D = 1))", y = "Aspect") +
    theme_minimal(base_size = 10) +
    theme(axis.text.y = element_text(size = 8),
          panel.background = element_rect(fill = "white", color = NA),
          plot.background = element_rect(fill = "white", color = NA),
          panel.grid = element_blank())
}

# Binscatter functions (Δ1 / Δ2 notation in labels). SE only (no 2xSE).
plot_binscatter_j <- function(binscatter_data) {
  model1_coef <- coef(model1)["delta_j"]
  model1_se <- sqrt(diag(vcovCL(model1, cluster = ~ respondent_id)))["delta_j"]
  subtitle_txt <- paste0("OLS slope = ", fmt(model1_coef, 4), " (SE = ", fmt(model1_se, 4), ")")
  ggplot() +
    geom_point(data = binscatter_data, aes(x = mean_delta_j, y = mean_choice),
               size = 3, alpha = 0.7, color = "steelblue") +
    geom_smooth(data = data, aes(x = delta_j, y = choice), method = "lm", se = TRUE, color = "red", linetype = "dashed") +
    labs(title = "Model 1: Binscatter Plot - Δ1",
         subtitle = subtitle_txt,
         x = expression(Delta[1] ~ "(within bin)"),
         y = expression(Pr(D == 1))) +
    theme_minimal() +
    theme(panel.background = element_rect(fill = "white", color = NA),
          plot.background = element_rect(fill = "white", color = NA),
          panel.grid = element_blank())
}

plot_binscatter_k <- function(binscatter_data) {
  model1_coef_k <- coef(model1)["delta_k"]
  model1_se_k <- sqrt(diag(vcovCL(model1, cluster = ~ respondent_id)))["delta_k"]
  subtitle_txt <- paste0("OLS slope = ", fmt(model1_coef_k, 4), " (SE = ", fmt(model1_se_k, 4), ")")
  ggplot() +
    geom_point(data = binscatter_data, aes(x = mean_delta_k, y = mean_choice),
               size = 3, alpha = 0.7, color = "steelblue") +
    geom_smooth(data = data, aes(x = delta_k, y = choice), method = "lm", se = TRUE, color = "red", linetype = "dashed") +
    labs(title = "Model 1: Binscatter Plot - Δ2",
         subtitle = subtitle_txt,
         x = expression(Delta[2] ~ "(within bin)"),
         y = expression(Pr(D == 1))) +
    theme_minimal() +
    theme(panel.background = element_rect(fill = "white", color = NA),
          plot.background = element_rect(fill = "white", color = NA),
          panel.grid = element_blank())
}

plot_binscatter_by_importance <- function(binscatter_data) {
  model3_baseline <- coef(model3)["delta_j"]
  vcov_m3 <- vcovCL(model3, cluster = ~ respondent_id)
  model3_se_baseline <- sqrt(diag(vcov_m3))["delta_j"]
  subtitle_txt <- paste0("Below median slope = ", fmt(model3_baseline, 4),
                         "; Above median slope = ", fmt(model3_baseline + coef(model3)["delta_j:am_imp_j"], 4),
                         " (SE baseline = ", fmt(model3_se_baseline, 4), ")")
  ggplot() +
    geom_point(data = binscatter_data, aes(x = mean_delta_j, y = mean_choice, color = importance_group),
               size = 3, alpha = 0.7) +
    # removed CI ribbons as requested for by-aspect importance plots
    geom_smooth(data = data %>% filter(!is.na(am_imp_j)) %>% mutate(importance_group = ifelse(am_imp_j == 1, "Above median", "Below median")),
                aes(x = delta_j, y = choice, color = importance_group),
                method = "lm", se = FALSE, linetype = "dashed") +
    scale_color_manual(values = c("Above median" = "red", "Below median" = "blue")) +
    labs(title = "Model 3: Binscatter Plot - Δ1 by Importance",
         subtitle = subtitle_txt,
         x = expression(Delta[1] ~ "(within bin)"),
         y = expression(Pr(D == 1)),
         color = "Importance") +
    theme_minimal() +
    theme(panel.background = element_rect(fill = "white", color = NA),
          plot.background = element_rect(fill = "white", color = NA),
          panel.grid = element_blank(), legend.position = "bottom")
}

# Combined binscatter (Δ1 & Δ2)
plot_binscatter_combined <- function(binscatter_data_j, binscatter_data_k) {
  model1_coef_j <- coef(model1)["delta_j"]
  model1_coef_k <- coef(model1)["delta_k"]
  model1_ses <- sqrt(diag(vcovCL(model1, cluster = ~ respondent_id)))
  subtitle_txt <- paste0("Δ1 slope = ", fmt(model1_coef_j,4), " (SE=", fmt(model1_ses["delta_j"],4), "); ",
                         "Δ2 slope = ", fmt(model1_coef_k,4), " (SE=", fmt(model1_ses["delta_k"],4), ")")
  combined_bins <- bind_rows(
    binscatter_data_j %>% select(mean_delta = mean_delta_j, mean_choice) %>% mutate(variable = "Δ1"),
    binscatter_data_k %>% select(mean_delta = mean_delta_k, mean_choice) %>% mutate(variable = "Δ2")
  )
  combined_raw <- bind_rows(
    data %>% select(delta = delta_j, choice) %>% mutate(variable = "Δ1"),
    data %>% select(delta = delta_k, choice) %>% mutate(variable = "Δ2")
  )
  ggplot() +
    geom_point(data = combined_bins, aes(x = mean_delta, y = mean_choice, color = variable), size = 3, alpha = 0.7) +
    geom_smooth(data = combined_raw, aes(x = delta, y = choice, color = variable), method = "lm", se = TRUE, linetype = "dashed") +
    scale_color_manual(values = c("Δ1" = "red", "Δ2" = "blue"),
                       labels = c(expression(Delta[1]), expression(Delta[2]))) +
    labs(title = "Model 1: Combined Binscatter Plot - Δ1 and Δ2",
         subtitle = subtitle_txt,
         x = expression(Delta ~ "value (within bin)"),
         y = expression(Pr(D == 1)), color = "Variable") +
    theme_minimal() +
    theme(panel.background = element_rect(fill = "white", color = NA),
          plot.background = element_rect(fill = "white", color = NA),
          panel.grid = element_blank(), legend.position = "bottom")
}

# -----------------------------------------------------------------------------
# Aspect-specific helpers: focal-as-j subset for both series (own vs other),
# binned means only. SE for pooled linear combos via lincomb_se helper.
# -----------------------------------------------------------------------------

create_aspect_binscatter_data <- function(aspect_name, position = "j", n_bins = 10) {
  aspect_name <- as.character(aspect_name)
  df <- data %>% filter(as.character(aspect_j) == aspect_name)
  if (nrow(df) == 0) {
    return(tibble(delta_bin = integer(), mean_delta = numeric(), mean_choice = numeric(),
                  n_obs = integer(), position = character(), aspect = character()))
  }
  if (position == "j") {
    aspect_data <- df %>%
      mutate(delta_bin = ntile(delta_j, n_bins)) %>%
      group_by(delta_bin) %>%
      summarise(mean_delta = mean(delta_j), mean_choice = mean(choice), n_obs = n(), .groups = "drop")
  } else {
    aspect_data <- df %>%
      mutate(delta_bin = ntile(delta_k, n_bins)) %>%
      group_by(delta_bin) %>%
      summarise(mean_delta = mean(delta_k), mean_choice = mean(choice), n_obs = n(), .groups = "drop")
  }
  aspect_data$position <- position
  aspect_data$aspect <- aspect_name
  return(aspect_data)
}

plot_aspect_binscatter <- function(aspect_name, n_bins = 10, digits = 4) {
  aspect_name <- as.character(aspect_name)
  # pooled (model2) linear combos -> use lincomb_se helper
  coefs2 <- coef(model2)
  j_term <- "delta_j"
  k_term <- "delta_k"
  j_inter_term <- paste0("delta_j:group_aj", aspect_name)
  k_inter_term <- paste0("delta_k:group_ak", aspect_name)
  j_total <- coefs2[j_term] + ifelse(j_inter_term %in% names(coefs2), coefs2[j_inter_term], 0)
  k_total <- coefs2[k_term] + ifelse(k_inter_term %in% names(coefs2), coefs2[k_inter_term], 0)
  se_j_total <- lincomb_se(model2, c(j_term, j_inter_term))
  se_k_total <- lincomb_se(model2, c(k_term, k_inter_term))
  
  # binned series computed from focal-J subset
  j_bins <- create_aspect_binscatter_data(aspect_name, "j", n_bins = n_bins)
  k_bins <- create_aspect_binscatter_data(aspect_name, "k", n_bins = n_bins)
  combined_bins <- bind_rows(j_bins, k_bins) %>%
    mutate(series = factor(position, levels = c("j", "k"), labels = c("Own (Δ1)", "Other (Δ2)")))
  
  # binned-fit slopes (regressions on the binned means)
  slope_j_binned <- tryCatch(coef(lm(mean_choice ~ mean_delta, data = combined_bins %>% filter(series == "Own (Δ1)")))[2],
                             error = function(e) NA_real_)
  slope_k_binned <- tryCatch(coef(lm(mean_choice ~ mean_delta, data = combined_bins %>% filter(series == "Other (Δ2)")))[2],
                             error = function(e) NA_real_)
  
  subtitle_txt <- paste0(
    "coef (resid.): Δ1 ", fmt(j_total, digits), " (SE=", fmt(se_j_total, digits), "); ",
    "Δ2 ", fmt(k_total, digits), " (SE=", fmt(se_k_total, digits), "); ",
    "Binned slopes — own: ", fmt(slope_j_binned, digits), ", other: ", fmt(slope_k_binned, digits)
  )
  
  p <- ggplot(combined_bins, aes(x = mean_delta, y = mean_choice, color = series)) +
    geom_point(size = 2.5, alpha = 0.9) +
    scale_color_manual(values = c("Own (Δ1)" = "red", "Other (Δ2)" = "blue"),
                       labels = c(expression(Delta[1]), expression(Delta[2]))) +
    labs(title = paste(" ", str_trunc(aspect_name, 40)),
         subtitle = subtitle_txt,
         x = expression(Delta ~ "(within bin)"),
         y = expression(Pr(D == 1)),
         color = "") +
    theme_minimal() +
    theme(
      legend.position = "bottom",
      panel.background = element_rect(fill = "white", color = NA),
      plot.background  = element_rect(fill = "white", color = NA),
      panel.grid = element_blank()
    )
  
  # add binned-fit lines (from regressions on binned means)
  fit_j <- tryCatch(coef(lm(mean_choice ~ mean_delta, data = combined_bins %>% filter(series == "Own (Δ1)"))),
                    error = function(e) NULL)
  fit_k <- tryCatch(coef(lm(mean_choice ~ mean_delta, data = combined_bins %>% filter(series == "Other (Δ2)"))),
                    error = function(e) NULL)
  if (!is.null(fit_j)) p <- p + geom_abline(intercept = fit_j[1], slope = fit_j[2], color = "red", linetype = "dashed")
  if (!is.null(fit_k)) p <- p + geom_abline(intercept = fit_k[1], slope = fit_k[2], color = "blue", linetype = "dashed")
  
  return(p)
}

plot_aspect_importance_binscatter <- function(aspect_name, position = "j", n_bins = 10, digits = 4) {
  aspect_name <- as.character(aspect_name)
  if (position == "j") {
    aspect_data <- data %>%
      filter(as.character(aspect_j) == aspect_name, !is.na(am_imp_j)) %>%
      group_by(am_imp_j) %>%
      mutate(delta_bin = ntile(delta_j, n_bins)) %>%
      group_by(am_imp_j, delta_bin) %>%
      summarise(mean_delta = mean(delta_j), mean_choice = mean(choice), n_obs = n(), .groups = "drop") %>%
      mutate(importance_group = ifelse(am_imp_j == 1, "Above median", "Below median"))
    raw_data <- data %>%
      filter(as.character(aspect_j) == aspect_name, !is.na(am_imp_j)) %>%
      select(delta = delta_j, choice, am_imp_j) %>%
      mutate(importance_group = ifelse(am_imp_j == 1, "Above median", "Below median"))
  } else {
    aspect_data <- data %>%
      filter(as.character(aspect_j) == aspect_name, !is.na(am_imp_j)) %>%
      group_by(am_imp_j) %>%
      mutate(delta_bin = ntile(delta_k, n_bins)) %>%
      group_by(am_imp_j, delta_bin) %>%
      summarise(mean_delta = mean(delta_k), mean_choice = mean(choice), n_obs = n(), .groups = "drop") %>%
      mutate(importance_group = ifelse(am_imp_j == 1, "Above median", "Below median"))
    raw_data <- data %>%
      filter(as.character(aspect_j) == aspect_name, !is.na(am_imp_j)) %>%
      select(delta = delta_k, choice, am_imp_j) %>%
      mutate(importance_group = ifelse(am_imp_j == 1, "Above median", "Below median"))
  }
  
  vcov_m3 <- vcovCL(model3, cluster = ~ respondent_id)
  model3_baseline <- coef(model3)["delta_j"]
  model3_se_baseline <- if (!is.null(vcov_m3) && "delta_j" %in% rownames(vcov_m3)) sqrt(vcov_m3["delta_j","delta_j"]) else NA_real_
  
  slope_above <- NA_real_; slope_below <- NA_real_
  if ("Above median" %in% aspect_data$importance_group) {
    slope_above <- tryCatch(coef(lm(mean_choice ~ mean_delta, data = aspect_data %>% filter(importance_group == "Above median")))[2],
                            error = function(e) NA_real_)
  }
  if ("Below median" %in% aspect_data$importance_group) {
    slope_below <- tryCatch(coef(lm(mean_choice ~ mean_delta, data = aspect_data %>% filter(importance_group == "Below median")))[2],
                            error = function(e) NA_real_)
  }
  
  subtitle_txt <- paste0("Baseline (model3) = ", fmt(model3_baseline, digits),
                         " (SE=", fmt(model3_se_baseline, digits), "); ",
                         "Binned slopes — Above: ", fmt(slope_above, digits),
                         ", Below: ", fmt(slope_below, digits))
  
  ggplot() +
    geom_point(data = aspect_data, aes(x = mean_delta, y = mean_choice, color = importance_group), size = 2, alpha = 0.7) +
    geom_smooth(data = raw_data, aes(x = delta, y = choice, color = importance_group), method = "lm", se = FALSE, linetype = "dashed") +
    scale_color_manual(values = c("Above median" = "red", "Below median" = "blue")) +
    labs(title = paste("Importance:", str_trunc(aspect_name, 30), "- Pos", toupper(position)),
         subtitle = subtitle_txt,
         x = expression(Delta[1] ~ "(within bin)"),
         y = expression(Pr(D == 1)),
         color = "Importance") +
    theme_minimal() +
    theme(panel.background = element_rect(fill = "white", color = NA),
          plot.background = element_rect(fill = "white", color = NA),
          panel.grid = element_blank(), legend.position = "bottom")
}

# -----------------------------------------------------------------------------
# New: Model3-like plot but split by above-median *level* (not importance)
# -----------------------------------------------------------------------------
plot_binscatter_by_level <- function(binscatter_data_level) {
  model3_baseline <- coef(model3)["delta_j"]
  vcov_m3 <- vcovCL(model3, cluster = ~ respondent_id)
  model3_se_baseline <- if (!is.null(vcov_m3) && "delta_j" %in% rownames(vcov_m3)) sqrt(vcov_m3["delta_j","delta_j"]) else NA_real_
  subtitle_txt <- paste0("Baseline (model3) = ", fmt(model3_baseline, 4), " (SE=", fmt(model3_se_baseline, 4), ")")
  ggplot() +
    geom_point(data = binscatter_data_level, aes(x = mean_delta_j, y = mean_choice, color = level_group), size = 3, alpha = 0.7) +
    geom_smooth(data = data %>% filter(!is.na(am_level_j)) %>% mutate(level_group = ifelse(am_level_j == 1, "Above median level", "Below median level")),
                aes(x = delta_j, y = choice, color = level_group), method = "lm", se = FALSE, linetype = "dashed") +
    scale_color_manual(values = c("Above median level" = "red", "Below median level" = "blue")) +
    labs(title = "Model 3b: Binscatter Plot - Δ1 by Above-median Level",
         subtitle = subtitle_txt,
         x = expression(Delta[1] ~ "(within bin)"),
         y = expression(Pr(D == 1)),
         color = "Level") +
    theme_minimal() +
    theme(panel.background = element_rect(fill = "white", color = NA),
          plot.background = element_rect(fill = "white", color = NA),
          panel.grid = element_blank(), legend.position = "bottom")
}

# =============================================================================
# CREATE AND SAVE PLOTS - NOW SAVING TO output_dir
# =============================================================================

p1 <- plot_model4(model4_effects)
print(p1)
ggsave(file.path(output_dir, "model4_importance_effects.png"), p1, width = 10, height = 6, dpi = 300)

p2 <- plot_interactions(model2_int_j, "Model 2: Δ1 Interactions by Aspect")
print(p2)
ggsave(file.path(output_dir, "model2_interactions_delta_j.png"), p2, width = 12, height = 8, dpi = 300)

p3 <- plot_main_effects(model2_main_j, "Main Effects: Δ1")
print(p3)
ggsave(file.path(output_dir, "model2_main_effects_delta_j.png"), p3, width = 12, height = 8, dpi = 300)

p4 <- plot_interactions(model5_effects, "Model 5: Aspect-Specific Importance Effects")
print(p4)
ggsave(file.path(output_dir, "model5_aspect_importance.png"), p4, width = 12, height = 8, dpi = 300)

p5 <- plot_binscatter_j(binscatter_data_j)
print(p5)
ggsave(file.path(output_dir, "model1_binscatter_delta_j.png"), p5, width = 10, height = 6, dpi = 300)

p6 <- plot_binscatter_k(binscatter_data_k)
print(p6)
ggsave(file.path(output_dir, "model1_binscatter_delta_k.png"), p6, width = 10, height = 6, dpi = 300)

p7 <- plot_binscatter_by_importance(binscatter_data_j_by_imp)
print(p7)
ggsave(file.path(output_dir, "model3_binscatter_by_importance.png"), p7, width = 10, height = 6, dpi = 300)

# new Model 3b by above-median level
p3b <- plot_binscatter_by_level(binscatter_data_j_by_level)
print(p3b)
ggsave(file.path(output_dir, "model3b_binscatter_by_level.png"), p3b, width = 10, height = 6, dpi = 300)

p8 <- plot_binscatter_combined(binscatter_data_j, binscatter_data_k)
print(p8)
ggsave(file.path(output_dir, "model1_binscatter_combined.png"), p8, width = 10, height = 6, dpi = 300)

# =============================================================================
# ASPECT-BY-ASPECT BINSCATTERS (print + save individual + collect for 4x4 pages)
# =============================================================================

plots_positions <- list()
plots_imp_j     <- list()
plots_imp_k     <- list()

unique_aspects <- sort(unique(c(data$aspect_j, data$aspect_k)))
cat("\nGenerating aspect-by-aspect binscatters for", length(unique_aspects), "aspects...\n")

for (i in seq_along(unique_aspects)) {
  aspect <- unique_aspects[i]
  safe_name <- gsub("[^A-Za-z0-9]", "_", aspect)
  safe_name <- substr(safe_name, 1, 50)
  
  p_aspect <- plot_aspect_binscatter(aspect)
  print(p_aspect)
  ggsave(file.path(output_dir, paste0("aspect_", sprintf("%02d", i), "_", safe_name, "_positions.png")),
         p_aspect, width = 10, height = 6, dpi = 300)
  
  p_importance_j <- plot_aspect_importance_binscatter(aspect, "j")
  print(p_importance_j)
  ggsave(file.path(output_dir, paste0("aspect_", sprintf("%02d", i), "_", safe_name, "_importance_j.png")),
         p_importance_j, width = 10, height = 6, dpi = 300)
  
  p_importance_k <- plot_aspect_importance_binscatter(aspect, "k")
  print(p_importance_k)
  ggsave(file.path(output_dir, paste0("aspect_", sprintf("%02d", i), "_", safe_name, "_importance_k.png")),
         p_importance_k, width = 10, height = 6, dpi = 300)
  
  plots_positions <- c(plots_positions, list(p_aspect))
  plots_imp_j     <- c(plots_imp_j,     list(p_importance_j))
  plots_imp_k     <- c(plots_imp_k,     list(p_importance_k))
  
  if (i %% 5 == 0) cat("Completed", i, "of", length(unique_aspects), "aspects\n")
}

cat("Generated", length(unique_aspects) * 3, "aspect-specific binscatter plots\n")

save_grid_pages <- function(plot_list, base_fname, ncol = 4, nrow = 4, width = 16, height = 12, dpi = 300) {
  chunk_size <- ncol * nrow
  n_pages <- ceiling(length(plot_list) / chunk_size)
  for (page in seq_len(n_pages)) {
    start <- (page - 1) * chunk_size + 1
    end <- min(start + chunk_size - 1, length(plot_list))
    chunk <- plot_list[start:end]
    if (length(chunk) < chunk_size) {
      n_missing <- chunk_size - length(chunk)
      chunk <- c(chunk, replicate(n_missing, ggplot() + theme_void(), simplify = FALSE))
    }
    grid <- plot_grid(plotlist = chunk, ncol = ncol, nrow = nrow)
    outname <- file.path(output_dir, sprintf("%s_page%02d.png", base_fname, page))
    ggsave(outname, grid, width = width, height = height, dpi = dpi)
    cat("Saved", outname, "\n")
  }
}

save_grid_pages(plots_positions, "aspects_positions_4x4")
save_grid_pages(plots_imp_j,     "aspects_importance_j_4x4")
save_grid_pages(plots_imp_k,     "aspects_importance_k_4x4")

cat("Saved 4x4 combined pages for positions and importance plots.\n")

# Summaries
cat("\n=== MODEL 4 SUMMARY ===\n")
for (i in 1:nrow(model4_effects)) {
  cat(model4_effects$effect_type[i], ":", round(model4_effects$coefficient[i], 4),
      "(SE:", round(model4_effects$se[i], 4), ")\n")
}

cat("\n=== MODEL 5 SUMMARY ===\n")
cat("Aspect-specific importance effects:", nrow(model5_effects), "coefficients\n")
cat("Positive coefficients:", sum(model5_effects$coefficient > 0), "\n")
cat("Significant:", sum(model5_effects$significant), "\n")

# Model 4 variant: above-median *level* interactions (Δ1 / Δ2)
model4_level <- lm(choice ~ delta_j + delta_j:am_level_j + delta_k + delta_k:am_level_k, data = data)
coeftest(model4_level, vcov = vcovCL(model4_level, cluster = ~ respondent_id))

# Extract clustered results and label with Δ notation
model4_level_j <- get_clustered_results(model4_level, "delta_j:am_level_j") %>%
  mutate(effect_type = "Δ1 × 1(above median level 1)")
model4_level_k <- get_clustered_results(model4_level, "delta_k:am_level_k") %>%
  mutate(effect_type = "Δ2 × 1(above median level 2)")

model4_level_effects <- bind_rows(model4_level_j, model4_level_k)

# quick print / plot
print(model4_level_effects)
p <- plot_model4(model4_level_effects)
print(p)
ggsave(file.path(output_dir, "model4_level_effects.png"), p, width = 8, height = 4, dpi = 300)

##### raw coefs 
model_beta <- lm(choice ~ delta_j + delta_k, data = data)
vc_beta <- vcovCL(model_beta, cluster = ~ respondent_id)
est <- coef(model_beta)[c("delta_j","delta_k")]
se  <- sqrt(diag(vc_beta))[c("delta_j","delta_k")]

df_beta <- tibble(
  id = factor(c("Delta1","Delta2"), levels = c("Delta1","Delta2")),
  estimate = as.numeric(est),
  se = as.numeric(se),
  ci_lower = estimate - 1.96 * se,
  ci_upper = estimate + 1.96 * se
)

# reverse axis so Delta1 appears on top
p_beta <- ggplot(df_beta, aes(y = id, x = estimate)) +
  geom_pointrange(aes(xmin = ci_lower, xmax = ci_upper), size = 0.8, fatten = 2) +
  geom_vline(xintercept = 0, linetype = "dashed", color = "gray50") +
  scale_y_discrete(limits = rev(levels(df_beta$id)),
                   labels = rev(c(expression(Delta[1]), expression(Delta[2])))) +
  labs(title = expression("Estimates:" ~ beta[1] ~ "and" ~ beta[2]),
       x = expression("Estimate (Pr(D == 1))"), y = "") +
  theme_minimal() +
  theme(panel.background = element_rect(fill = "white", color = NA),
        plot.background = element_rect(fill = "white", color = NA),
        panel.grid = element_blank())

print(p_beta)
ggsave(file.path(output_dir, "beta1_beta2_ci.png"), p_beta, width = 6, height = 3, dpi = 300)

### importance coefficients 
model4 <- lm(choice ~ delta_j + delta_j:am_imp_j + delta_k + delta_k:am_imp_k, data = data)
coefs4 <- coef(model4)

desired_order_imp <- c("Delta1_imp_above","Delta1_imp_below","Delta2_imp_above","Delta2_imp_below")
combos_imp <- tibble(
  label = desired_order_imp,
  terms = list(
    c("delta_j","delta_j:am_imp_j"), # Δ1 imp above
    c("delta_j"),                     # Δ1 imp below
    c("delta_k","delta_k:am_imp_k"), # Δ2 imp above
    c("delta_k")                      # Δ2 imp below
  ),
  weights = list(c(1,1),c(1),c(1,1),c(1))
)

combos_imp <- combos_imp %>%
  mutate(
    estimate = map2_dbl(terms, weights, ~ sum(coefs4[.x] * .y, na.rm = TRUE)),
    se = map2_dbl(terms, weights, ~ lincomb_se(model4, .x, .y)),
    ci_lower = estimate - 1.96 * se,
    ci_upper = estimate + 1.96 * se,
    label_f = factor(label, levels = label)
  )

y_labels_imp <- c(expression(Delta[1] ~ "(imp=above)"),
                  expression(Delta[1] ~ "(imp=below)"),
                  expression(Delta[2] ~ "(imp=above)"),
                  expression(Delta[2] ~ "(imp=below)"))

p4_imp <- ggplot(combos_imp, aes(y = label_f, x = estimate)) +
  geom_pointrange(aes(xmin = ci_lower, xmax = ci_upper), size = 0.8, fatten = 2) +
  geom_vline(xintercept = 0, linetype = "dashed", color = "gray50") +
  scale_y_discrete(limits = rev(levels(combos_imp$label_f)), labels = rev(y_labels_imp)) +
  labs(title = "Model4 (importance): Δ coefficients by importance",
       x = expression("Estimate (Pr(D == 1))"), y = "") +
  theme_minimal() +
  theme(panel.background = element_rect(fill = "white", color = NA),
        plot.background = element_rect(fill = "white", color = NA),
        panel.grid = element_blank())

print(p4_imp)
ggsave(file.path(output_dir, "model4_importance_4coefs.png"), p4_imp, width = 6.5, height = 4, dpi = 300)

#### BINSCATTERS 
# ddiff = Δ2 - Δ1 (user requested)
data <- data %>% mutate(ddiff = delta_k - delta_j)

# ensure importance indicators exist (am_imp_j / am_imp_k)
if (!"am_imp_j" %in% names(data)) data <- data %>% mutate(am_imp_j = as.numeric(importance_j > median(importance_j, na.rm = TRUE)))
if (!"am_imp_k" %in% names(data)) data <- data %>% mutate(am_imp_k = as.numeric(importance_k > median(importance_k, na.rm = TRUE)))

# helper to get clustered SE and slope for ddiff ~ choice (raw)
cluster_slope_ddiff_on_choice <- function(df) {
  if (nrow(df) < 5) return(list(slope = NA_real_, se = NA_real_))
  m <- lm(ddiff ~ choice, data = df)
  vc <- vcovCL(m, cluster = ~ respondent_id)
  coefm <- coef(m)
  slope <- if ("choice" %in% names(coefm)) coefm["choice"] else NA_real_
  se    <- if (!is.null(vc) && "choice" %in% rownames(vc)) sqrt(diag(vc))["choice"] else NA_real_
  list(slope = slope, se = se)
}

# Use fixed global bins so comparisons are consistent
n_bins <- 20
breaks <- quantile(data$ddiff, probs = seq(0,1,length.out = n_bins+1), na.rm = TRUE)
# make sure breaks are unique (if ties, jitter slightly)
if (any(duplicated(breaks))) {
  breaks <- jitter(breaks, amount = 1e-10)
  breaks <- sort(breaks)
}
bin_label_vec <- seq_len(n_bins)

# ---------- 1) Overall binscatter (x = Pr(D=1), y = Δ2 - Δ1) ----------
bins_overall <- data %>%
  mutate(bin = cut(ddiff, breaks = breaks, include.lowest = TRUE, labels = bin_label_vec)) %>%
  filter(!is.na(bin)) %>%
  group_by(bin) %>%
  summarise(mean_ddiff = mean(ddiff), mean_choice = mean(choice), n = n(), .groups = "drop")

# overall raw slope (ddiff ~ choice) with clustered SE
overall_res <- cluster_slope_ddiff_on_choice(data)

fit_bins_overall <- tryCatch(coef(lm(mean_ddiff ~ mean_choice, data = bins_overall)), error = function(e) NULL)

p_overall <- ggplot(bins_overall, aes(x = mean_choice, y = mean_ddiff)) +
  geom_point(size = 3, alpha = 0.85, color = "steelblue") +
  { if (!is.null(fit_bins_overall)) geom_abline(intercept = fit_bins_overall[1], slope = fit_bins_overall[2], linetype = "dashed", color = "red") } +
  labs(
    title = expression("Binscatter: " ~ Delta[2] - Delta[1] ~ " vs Pr(D == 1) (overall)"),
    subtitle = paste0("Raw slope (ddiff ~ choice) = ", fmt(overall_res$slope,4), " (SE=", fmt(overall_res$se,4), ")"),
    x = expression(Pr(D == 1)),
    y = expression(Delta[2] - Delta[1] ~ "(binned means)")
  ) +
  theme_minimal() +
  theme(panel.background = element_rect(fill = "white", color = NA),
        plot.background  = element_rect(fill = "white", color = NA),
        panel.grid = element_blank())

print(p_overall)
ggsave(file.path(output_dir, "binscatter_ddiff_overall_y_vs_x_amimp.png"), p_overall, width = 8, height = 4, dpi = 300)

# ---------- 2) By above/below median IMPORTANCE of aspect 1 (am_imp_j) ----------
bins_by_imp <- data %>%
  filter(!is.na(am_imp_j)) %>%
  mutate(bin = cut(ddiff, breaks = breaks, include.lowest = TRUE, labels = bin_label_vec)) %>%
  filter(!is.na(bin)) %>%
  group_by(am_imp_j, bin) %>%
  summarise(mean_ddiff = mean(ddiff), mean_choice = mean(choice), n = n(), .groups = "drop") %>%
  mutate(imp_label = ifelse(am_imp_j == 1, "Above median importance", "Below median importance"))

# subgroup raw slopes
res_imp_above <- cluster_slope_ddiff_on_choice(filter(data, am_imp_j == 1))
res_imp_below <- cluster_slope_ddiff_on_choice(filter(data, am_imp_j == 0))

fits_by_imp <- bins_by_imp %>%
  group_by(imp_label) %>%
  do({
    fit <- tryCatch(coef(lm(mean_ddiff ~ mean_choice, data = .)), error = function(e) NULL)
    tibble(intercept = ifelse(is.null(fit), NA_real_, fit[1]), slope = ifelse(is.null(fit), NA_real_, fit[2]))
  }) %>% ungroup()

p_by_imp <- ggplot(bins_by_imp, aes(x = mean_choice, y = mean_ddiff, color = imp_label, shape = imp_label)) +
  geom_point(size = 2.8, alpha = 0.9) +
  geom_abline(data = fits_by_imp %>% filter(!is.na(slope)), aes(intercept = intercept, slope = slope, color = imp_label), linetype = "dashed") +
  scale_color_manual(values = c("Above median importance" = "red", "Below median importance" = "blue")) +
  scale_shape_manual(values = c(16, 17)) +
  labs(
    title = expression("Binscatter: " ~ Delta[2] - Delta[1] ~ " vs Pr(D == 1) — Above vs Below median importance 1"),
    subtitle = paste0("Above: slope=", fmt(res_imp_above$slope,4), " (SE=", fmt(res_imp_above$se,4), "); ",
                      "Below: slope=", fmt(res_imp_below$slope,4), " (SE=", fmt(res_imp_below$se,4), ")"),
    x = expression(Pr(D == 1)),
    y = expression(Delta[2] - Delta[1] ~ "(binned means)"),
    color = "", shape = ""
  ) +
  theme_minimal() +
  theme(legend.position = "bottom",
        panel.background = element_rect(fill = "white", color = NA),
        plot.background = element_rect(fill = "white", color = NA),
        panel.grid = element_blank())

print(p_by_imp)
ggsave(file.path(output_dir, "binscatter_ddiff_by_imp_y_vs_x.png"), p_by_imp, width = 8, height = 4, dpi = 300)

# ---------- 3) Contrast binscatter (importance-based groups) ----------
# Group A: am_imp_j == 1 & am_imp_k == 0
# Group B: am_imp_j == 0 & am_imp_k == 1
data <- data %>%
  mutate(contrast_group = case_when(
    am_imp_j == 1 & am_imp_k == 0 ~ "A: imp1=above & imp2=below",
    am_imp_j == 0 & am_imp_k == 1 ~ "B: imp1=below & imp2=above",
    TRUE ~ NA_character_
  ))

bins_contrast <- data %>%
  filter(!is.na(contrast_group)) %>%
  mutate(bin = cut(ddiff, breaks = breaks, include.lowest = TRUE, labels = bin_label_vec)) %>%
  filter(!is.na(bin)) %>%
  group_by(contrast_group, bin) %>%
  summarise(mean_ddiff = mean(ddiff), mean_choice = mean(choice), n = n(), .groups = "drop")

# subgroup raw slopes (ddiff ~ choice)
res_A <- cluster_slope_ddiff_on_choice(filter(data, contrast_group == "A: imp1=above & imp2=below"))
res_B <- cluster_slope_ddiff_on_choice(filter(data, contrast_group == "B: imp1=below & imp2=above"))

fits_contrast <- bins_contrast %>%
  group_by(contrast_group) %>%
  do({
    fit <- tryCatch(coef(lm(mean_ddiff ~ mean_choice, data = .)), error = function(e) NULL)
    tibble(intercept = ifelse(is.null(fit), NA_real_, fit[1]), slope = ifelse(is.null(fit), NA_real_, fit[2]))
  }) %>% ungroup()

p_contrast <- ggplot(bins_contrast, aes(x = mean_choice, y = mean_ddiff, color = contrast_group, shape = contrast_group)) +
  geom_point(size = 2.8, alpha = 0.95) +
  geom_abline(data = fits_contrast %>% filter(!is.na(slope)), aes(intercept = intercept, slope = slope, color = contrast_group), linetype = "dashed") +
  scale_color_manual(values = c("A: imp1=above & imp2=below" = "purple", "B: imp1=below & imp2=above" = "darkgreen")) +
  scale_shape_manual(values = c(16, 17)) +
  labs(
    title = expression("Binscatter: " ~ Delta[2] - Delta[1] ~ " vs Pr(D == 1) — Contrast (importance)"),
    subtitle = paste0("A slope=", fmt(res_A$slope,4), " (SE=", fmt(res_A$se,4), "); ",
                      "B slope=", fmt(res_B$slope,4), " (SE=", fmt(res_B$se,4), ")"),
    x = expression(Pr(D == 1)),
    y = expression(Delta[2] - Delta[1] ~ "(binned means)"),
    color = "", shape = ""
  ) +
  theme_minimal() +
  theme(legend.position = "bottom",
        panel.background = element_rect(fill = "white", color = NA),
        plot.background = element_rect(fill = "white", color = NA),
        panel.grid = element_blank())

print(p_contrast)
ggsave(file.path(output_dir, "binscatter_ddiff_contrast_imp_y_vs_x.png"), p_contrast, width = 8, height = 4, dpi = 300)

cat("\n=== ALL PLOTS SAVED TO:", output_dir, "===\n")