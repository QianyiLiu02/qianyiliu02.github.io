# =============================================================================
# Temporal Discounting: Vertical vs. Horizontal Design Comparison
# March 5, 2026 pilot data
# =============================================================================

setwd("/Users/qianyi/Library/CloudStorage/Dropbox/swb/temporal discounting/temporal_design_3_5")

options(stringsAsFactors = FALSE)

library(dplyr)
library(ggplot2)

# -----------------------------------------------------------------------------
# PART 1: LOAD AND FILTER TO MARCH 5
# -----------------------------------------------------------------------------

raw_df <- read.csv(
  "/Users/qianyi/Library/CloudStorage/Dropbox/swb/temporal discounting/temporal_design_3_5/temporal_3_5.csv",
  stringsAsFactors = FALSE
)

# Remove header rows (row 1 = labels, row 2 = import IDs)
df <- raw_df[-(1:2), ]

# Keep only March 5, 2026
df$StartDate <- as.POSIXct(df$StartDate, format = "%Y-%m-%d %H:%M:%S")
df <- df[format(df$StartDate, "%Y-%m-%d") == "2026-03-05", ]

cat("Respondents from March 5:", nrow(df), "\n")

# Keep those with temporal discounting data
df <- df[df$allTemporalComparisons != "" & !is.na(df$allTemporalComparisons), ]

# Recode design variable
df$is_horizontal <- as.integer(df$isHorizontalDesign)

cat("Design breakdown:\n")
print(table(df$is_horizontal, useNA = "ifany"))

# -----------------------------------------------------------------------------
# PART 2: PARSE TEMPORAL COMPARISONS
# -----------------------------------------------------------------------------

parse_trial <- function(trial_str, respondent_id, is_horizontal) {
  parts <- strsplit(trial_str, ":")[[1]]
  if (length(parts) < 9) return(NULL)
  
  data.frame(
    respondent_id    = respondent_id,
    is_horizontal    = is_horizontal,
    aspect           = parts[1],
    tradeoff_type    = parts[2],
    baseline         = as.numeric(parts[3]),
    delta_w1         = as.numeric(parts[4]),
    delta_w2         = as.numeric(parts[5]),
    level1           = as.numeric(parts[6]),
    level2           = as.numeric(parts[7]),
    choice           = parts[8],
    response_time_ms = as.numeric(parts[9])
  )
}

all_trials <- list()
for (i in seq_len(nrow(df))) {
  resp_id       <- df$ResponseId[i]
  is_horiz      <- df$is_horizontal[i]
  trials_str    <- df$allTemporalComparisons[i]
  trials        <- strsplit(trials_str, "\\|")[[1]]
  for (trial in trials) {
    parsed <- parse_trial(trial, resp_id, is_horiz)
    if (!is.null(parsed)) all_trials[[length(all_trials) + 1]] <- parsed
  }
}
trials_df <- do.call(rbind, all_trials)

cat("\nTotal parsed trials:", nrow(trials_df), "\n")

# -----------------------------------------------------------------------------
# PART 3: CLEAN AND PREPARE
# -----------------------------------------------------------------------------

trials_df$chose_earlier <- ifelse(trials_df$choice == "Option1", 1, 0)

valid_trials <- trials_df[
  trials_df$choice %in% c("Option1", "Option2") &
    !is.na(trials_df$delta_w1) &
    !is.na(trials_df$delta_w2) &
    trials_df$delta_w1 != 0 &
    trials_df$delta_w2 != 0,
]

# Same-sign filter (main analysis)
valid_trials$same_sign <- (valid_trials$delta_w1 > 0 & valid_trials$delta_w2 > 0) |
  (valid_trials$delta_w1 < 0 & valid_trials$delta_w2 < 0)

same_sign_trials <- valid_trials[valid_trials$same_sign, ]

cat("\nSame-sign trials by design:\n")
print(table(same_sign_trials$is_horizontal))

# Subsets by design
trials_V <- same_sign_trials[same_sign_trials$is_horizontal == 0, ]
trials_H <- same_sign_trials[same_sign_trials$is_horizontal == 1, ]

# -----------------------------------------------------------------------------
# PART 4: PROBIT ESTIMATION BY DESIGN
# -----------------------------------------------------------------------------

run_probit <- function(data, label) {
  probit <- glm(chose_earlier ~ delta_w2 + delta_w1 - 1,
                data = data,
                family = binomial(link = "probit"))
  
  b2 <- coef(probit)["delta_w2"]
  b1 <- coef(probit)["delta_w1"]
  
  vc   <- vcov(probit)
  grad <- c(b2 / b1^2, -1 / b1)
  se_delta <- sqrt(t(grad) %*%
                     vc[c("delta_w1", "delta_w2"),
                        c("delta_w1", "delta_w2")] %*%
                     grad)
  
  delta_est <- -b2 / b1
  
  cat("\n--- Design:", label, "---\n")
  cat("N trials:", nrow(data), "\n")
  cat("gamma_1 (delta_w1):", round(b1, 5),
      " SE:", round(sqrt(vc["delta_w1", "delta_w1"]), 5), "\n")
  cat("gamma_2 (delta_w2):", round(b2, 5),
      " SE:", round(sqrt(vc["delta_w2", "delta_w2"]), 5), "\n")
  cat("delta estimate:    ", round(delta_est, 4),
      " SE:", round(se_delta, 4), "\n")
  cat("95% CI: [", round(delta_est - 1.96 * se_delta, 4),
      ",", round(delta_est + 1.96 * se_delta, 4), "]\n")
  
  list(
    model      = probit,
    label      = label,
    n          = nrow(data),
    gamma1     = b1,
    gamma2     = b2,
    se_gamma1  = sqrt(vc["delta_w1", "delta_w1"]),
    se_gamma2  = sqrt(vc["delta_w2", "delta_w2"]),
    delta      = delta_est,
    se_delta   = se_delta
  )
}

res_V <- run_probit(trials_V, "Vertical (0)")
res_H <- run_probit(trials_H, "Horizontal (1)")

# -----------------------------------------------------------------------------
# PART 5: LR TEST — POOLED UNRESTRICTED VS RESTRICTED
# -----------------------------------------------------------------------------

# Restricted: same coefficients for both designs
probit_restricted <- glm(
  chose_earlier ~ delta_w2 + delta_w1 - 1,
  data   = same_sign_trials,
  family = binomial(link = "probit")
)

# Unrestricted: coefficients allowed to differ by design
same_sign_trials$dw1_H <- same_sign_trials$delta_w1 * same_sign_trials$is_horizontal
same_sign_trials$dw2_H <- same_sign_trials$delta_w2 * same_sign_trials$is_horizontal

probit_unrestricted <- glm(
  chose_earlier ~ delta_w2 + delta_w1 + dw2_H + dw1_H - 1,
  data   = same_sign_trials,
  family = binomial(link = "probit")
)

lr_stat <- 2 * (logLik(probit_unrestricted) - logLik(probit_restricted))
lr_pval <- pchisq(as.numeric(lr_stat), df = 2, lower.tail = FALSE)

cat("\n--- LR Test: H0: gamma_1^V = gamma_1^H AND gamma_2^V = gamma_2^H ---\n")
cat("LR statistic:", round(lr_stat, 4), "\n")
cat("Degrees of freedom: 2\n")
cat("p-value:", round(lr_pval, 4), "\n")

# -----------------------------------------------------------------------------
# PART 6: COEFFICIENT PLOT — gamma_1 and gamma_2 by design
# -----------------------------------------------------------------------------

coef_df <- data.frame(
  design    = rep(c("Vertical", "Horizontal"), each = 2),
  param     = rep(c("gamma[1]~(Delta*w[1])", "gamma[2]~(Delta*w[2])"), 2),
  estimate  = c(res_V$gamma1, res_V$gamma2,
                res_H$gamma1, res_H$gamma2),
  se        = c(res_V$se_gamma1, res_V$se_gamma2,
                res_H$se_gamma1, res_H$se_gamma2)
)

coef_df$ci_lo <- coef_df$estimate - 1.96 * coef_df$se
coef_df$ci_hi <- coef_df$estimate + 1.96 * coef_df$se
coef_df$design <- factor(coef_df$design, levels = c("Vertical", "Horizontal"))

p_coef <- ggplot(coef_df, aes(x = design, y = estimate, color = design)) +
  geom_point(size = 3.5, position = position_dodge(0.4)) +
  geom_errorbar(aes(ymin = ci_lo, ymax = ci_hi),
                width = 0.15, linewidth = 0.8,
                position = position_dodge(0.4)) +
  geom_hline(yintercept = 0, linetype = "dashed", color = "gray50") +
  facet_wrap(~ param, scales = "free_y", labeller = label_parsed) +
  scale_color_manual(values = c("Vertical" = "#2166ac", "Horizontal" = "#d6604d")) +
  labs(
    title    = "Probit Coefficient Estimates by Design",
    subtitle = paste0("LR test p = ", round(lr_pval, 3),
                      "  |  V: n = ", res_V$n, ", H: n = ", res_H$n),
    x        = NULL,
    y        = "Estimate (95% CI)",
    color    = "Design"
  ) +
  theme_bw(base_size = 13) +
  theme(
    legend.position  = "none",
    strip.background = element_rect(fill = "gray95"),
    strip.text       = element_text(size = 12),
    plot.title       = element_text(face = "bold")
  )

ggsave("coef_plot_design_compare.pdf", p_coef, width = 7, height = 4.5)
ggsave("coef_plot_design_compare.png", p_coef, width = 7, height = 4.5, dpi = 300)

# -----------------------------------------------------------------------------
# PART 7: BINSCATTER — exact style matching temporal.R, one plot per design
# -----------------------------------------------------------------------------

make_binscatter_design <- function(data, design_label, n_bins = 10) {
  # Keep only positive delta_w1 and delta_w2 (same filter as temporal.R)
  pos <- data[data$delta_w1 > 0 & data$delta_w2 > 0, ]
  pos$rho <- pos$delta_w1 / pos$delta_w2
  pos <- pos[is.finite(pos$rho) & !is.na(pos$rho), ]
  
  if (nrow(pos) < 5) {
    message("Too few rows for ", design_label, " - skipping")
    return(invisible(NULL))
  }
  
  rho_breaks <- quantile(pos$rho,
                         probs = seq(0, 1, length.out = n_bins + 1),
                         na.rm = TRUE)
  rho_breaks <- unique(rho_breaks)
  if (length(rho_breaks) <= 2L)
    rho_breaks <- pretty(pos$rho, n = n_bins + 1)
  
  pos$rho_bin <- cut(pos$rho,
                     breaks         = rho_breaks,
                     include.lowest = TRUE,
                     labels         = FALSE)
  
  bin_df <- pos %>%
    dplyr::filter(!is.na(rho_bin)) %>%
    dplyr::group_by(rho_bin) %>%
    dplyr::summarise(
      rho_mean = mean(rho),
      p_early  = mean(chose_earlier),
      n_trials = dplyr::n(),
      .groups  = "drop"
    ) %>%
    dplyr::arrange(rho_mean)
  
  # Interpolated delta (crossing of 0.5 line)
  diff_sign <- sign(bin_df$p_early - 0.5)
  cross_idx <- which(diff(diff_sign) != 0)
  if (length(cross_idx) == 0L) {
    delta_bins <- NA
  } else {
    i1 <- cross_idx[1]; i2 <- i1 + 1
    delta_bins <- approx(
      x    = bin_df$p_early[c(i1, i2)],
      y    = bin_df$rho_mean[c(i1, i2)],
      xout = 0.5
    )$y
  }
  cat("Binscatter delta (", design_label, "):", round(delta_bins, 4), "\n")
  
  max_rho <- ceiling(max(bin_df$rho_mean, na.rm = TRUE))
  
  title_expr <- bquote(
    "Binscatter: Pr(chose earlier) vs " * rho *
      " (" * Delta * w[1] * ", " * Delta * w[2] *
      " > 0, all tradeoffs, " * .(design_label) * " design)"
  )
  
  p <- ggplot(bin_df, aes(x = rho_mean, y = p_early)) +
    geom_point() +
    geom_line() +
    geom_hline(yintercept = 0.5, linetype = "dashed") +
    scale_x_continuous(
      breaks = 0:max_rho,
      limits = c(0, max_rho)
    ) +
    scale_y_continuous(
      limits = c(0, 1),
      breaks = seq(0, 1, by = 0.2)
    ) +
    labs(
      x     = expression(rho == Delta * w[1] / Delta * w[2]),
      y     = "Pr(chose earlier)",
      title = title_expr
    )
  
  fname <- paste0("binscatter_all_", tolower(design_label), ".pdf")
  ggsave(filename = fname, plot = p, device = "pdf", width = 7, height = 4)
  cat("Saved:", fname, "\n")
  
  invisible(list(bin_df = bin_df, delta_bins = delta_bins, plot = p))
}

# One binscatter per design, pooled across all tradeoff types
make_binscatter_design(trials_V, "Vertical")
make_binscatter_design(trials_H, "Horizontal")

# -----------------------------------------------------------------------------
# PART 8: DELTA COMPARISON PLOT
# -----------------------------------------------------------------------------

delta_df <- data.frame(
  design   = c("Vertical", "Horizontal"),
  delta    = c(res_V$delta, res_H$delta),
  se_delta = c(res_V$se_delta, res_H$se_delta)
)
delta_df$ci_lo <- delta_df$delta - 1.96 * delta_df$se_delta
delta_df$ci_hi <- delta_df$delta + 1.96 * delta_df$se_delta
delta_df$design <- factor(delta_df$design, levels = c("Vertical", "Horizontal"))

p_delta <- ggplot(delta_df, aes(x = design, y = delta, color = design)) +
  geom_point(size = 4) +
  geom_errorbar(aes(ymin = ci_lo, ymax = ci_hi), width = 0.1, linewidth = 1) +
  scale_color_manual(values = c("Vertical" = "#2166ac", "Horizontal" = "#d6604d")) +
  labs(
    title    = expression(hat(delta) ~ "Estimates by Design (Delta Method SEs)"),
    subtitle = paste0(" V: n = ", res_V$n, ", H: n = ", res_H$n),
    x        = NULL,
    y        = expression(hat(delta)),
    color    = "Design"
  ) +
  theme_bw(base_size = 13) +
  theme(legend.position = "none",
        plot.title = element_text(face = "bold"))

ggsave("delta_compare.pdf", p_delta, width = 5, height = 4.5)
ggsave("delta_compare.png", p_delta, width = 5, height = 4.5, dpi = 300)

# -----------------------------------------------------------------------------
# PART 9: SUMMARY TABLE (printed to console)
# -----------------------------------------------------------------------------

cat("\n\n=== SUMMARY: Vertical vs. Horizontal Design ===\n")
cat(sprintf("%-25s %10s %10s\n", "Parameter", "Vertical", "Horizontal"))
cat(strrep("-", 47), "\n")
cat(sprintf("%-25s %10d %10d\n",   "N trials",
            res_V$n, res_H$n))
cat(sprintf("%-25s %10.5f %10.5f\n", "gamma_1",
            res_V$gamma1, res_H$gamma1))
cat(sprintf("%-25s %10.5f %10.5f\n", "  SE(gamma_1)",
            res_V$se_gamma1, res_H$se_gamma1))
cat(sprintf("%-25s %10.5f %10.5f\n", "gamma_2",
            res_V$gamma2, res_H$gamma2))
cat(sprintf("%-25s %10.5f %10.5f\n", "  SE(gamma_2)",
            res_V$se_gamma2, res_H$se_gamma2))
cat(sprintf("%-25s %10.4f %10.4f\n", "delta_hat",
            res_V$delta, res_H$delta))
cat(sprintf("%-25s %10.4f %10.4f\n", "  SE(delta_hat)",
            res_V$se_delta, res_H$se_delta))
cat(strrep("-", 47), "\n")
cat(sprintf("LR statistic: %.4f  (chi2(2), p = %.4f)\n", lr_stat, lr_pval))

cat("\nPlots saved:\n")
cat("  coef_plot_design_compare.pdf/.png\n")
cat("  binscatter_design_compare.pdf/.png\n")
cat("  delta_compare.pdf/.png\n")




# -----------------------------------------------------------------------------
# PART 4b: OLS WITH CLUSTERED SEs (individual level) BY DESIGN
# -----------------------------------------------------------------------------

library(fixest)

# Estimate by design arm
ols_V <- feols(chose_earlier ~ delta_w2 + delta_w1 | 0,
               data    = trials_V,
               cluster = ~respondent_id)

ols_H <- feols(chose_earlier ~ delta_w2 + delta_w1 | 0,
               data    = trials_H,
               cluster = ~respondent_id)

# Coefficient comparison plot — gamma_1 and gamma_2 side by side
ols_coef_df <- data.frame(
  design   = rep(c("Vertical", "Horizontal"), each = 2),
  param    = rep(c("gamma[1]~(Delta*w[1])", "gamma[2]~(Delta*w[2])"), 2),
  estimate = c(coef(ols_V)["delta_w1"], coef(ols_V)["delta_w2"],
               coef(ols_H)["delta_w1"], coef(ols_H)["delta_w2"]),
  se       = c(se(ols_V)["delta_w1"],   se(ols_V)["delta_w2"],
               se(ols_H)["delta_w1"],   se(ols_H)["delta_w2"])
)

ols_coef_df$ci_lo  <- ols_coef_df$estimate - 1.96 * ols_coef_df$se
ols_coef_df$ci_hi  <- ols_coef_df$estimate + 1.96 * ols_coef_df$se
ols_coef_df$design <- factor(ols_coef_df$design, levels = c("Vertical", "Horizontal"))

p_ols_coef <- ggplot(ols_coef_df, aes(x = design, y = estimate, color = design)) +
  geom_point(size = 3.5) +
  geom_errorbar(aes(ymin = ci_lo, ymax = ci_hi), width = 0.15, linewidth = 0.8) +
  geom_hline(yintercept = 0, linetype = "dashed", color = "gray50") +
  facet_wrap(~ param, scales = "free_y", labeller = label_parsed) +
  scale_color_manual(values = c("Vertical" = "#2166ac", "Horizontal" = "#d6604d")) +
  labs(
    title    = "OLS Coefficient Estimates by Design",
    subtitle = paste0("SEs clustered at individual-level | V: n = ", res_V$n, ", H: n = ", res_H$n),
    x        = NULL,
    y        = "Estimate (95% CI)",
    color    = "Design"
  ) +
  theme_bw(base_size = 13) +
  theme(
    legend.position  = "none",
    strip.background = element_rect(fill = "gray95"),
    strip.text       = element_text(size = 12),
    plot.title       = element_text(face = "bold")
  )

ggsave("ols_coef_plot_design_compare.pdf", p_ols_coef, width = 7, height = 4.5)

# Print side-by-side summary
etable(ols_V, ols_H,
       title      = "OLS: Chose Earlier ~ delta_w1 + delta_w2, by Design",
       headers    = list("Vertical", "Horizontal"),
       se.below   = TRUE,
       digits     = 5)




# -----------------------------------------------------------------------------
# MONOTONICITY VIOLATION ANALYSIS (Option 3)
# -----------------------------------------------------------------------------

mono_data <- same_sign_trials %>%
  mutate(
    rho      = delta_w1 / delta_w2,
    sign_grp = ifelse(delta_w1 > 0, "positive", "negative")
  ) %>%
  filter(is.finite(rho), !is.na(rho)) %>%
  arrange(respondent_id, rho)

# Count violations per person
violations <- mono_data %>%
  group_by(respondent_id, is_horizontal) %>%
  summarise(
    n_trials     = n(),
    n_adjacent   = n() - 1,
    # For each adjacent pair, only count if both trials have the same sign
    # positive pairs: violation is 1->0; negative pairs: violation is 0->1
    n_violations = {
      ce   <- chose_earlier
      sg   <- sign_grp
      pairs_same_sign <- sg[-n()] == sg[-1]
      is_pos_pair     <- sg[-n()] == "positive"
      viol_pos <- diff(ce) < 0  & is_pos_pair  & pairs_same_sign
      viol_neg <- diff(ce) > 0  & !is_pos_pair & pairs_same_sign
      sum(viol_pos | viol_neg)
    },
    n_comparable = {
      sg <- sign_grp
      sum(sg[-n()] == sg[-1])   # number of same-sign adjacent pairs
    },
    viol_rate = ifelse(n_comparable > 0,
                       n_violations / n_comparable,
                       NA_real_),
    .groups = "drop"
  ) %>%
  filter(n_trials >= 3, n_comparable >= 2)

violations$design <- factor(
  ifelse(violations$is_horizontal == 0, "Vertical", "Horizontal"),
  levels = c("Vertical", "Horizontal")
)

# Summary by design
viol_summary <- violations %>%
  group_by(design) %>%
  summarise(
    n_respondents  = n(),
    mean_viol_rate = mean(viol_rate, na.rm = TRUE),
    sd_viol_rate   = sd(viol_rate,   na.rm = TRUE),
    se_viol_rate   = sd(viol_rate,   na.rm = TRUE) / sqrt(n()),
    .groups        = "drop"
  )

print(viol_summary)

# T-test
viol_V <- violations$viol_rate[violations$is_horizontal == 0]
viol_H <- violations$viol_rate[violations$is_horizontal == 1]
ttest  <- t.test(viol_H, viol_V)

cat("\nT-test: Horizontal vs Vertical violation rates\n")
cat("Difference:", round(mean(viol_H, na.rm = TRUE) - mean(viol_V, na.rm = TRUE), 4), "\n")
cat("p-value:",    round(ttest$p.value, 4), "\n")

# -----------------------------------------------------------------------------
# PLOT 1: Histogram (relative frequency)
# -----------------------------------------------------------------------------

p_viol <- ggplot(violations, aes(x = viol_rate, fill = design)) +
  geom_histogram(aes(y = after_stat(count / sum(count))),
                 binwidth = 0.1, boundary = 0,
                 color = "white", alpha = 0.8) +
  geom_vline(data = viol_summary,
             aes(xintercept = mean_viol_rate, color = design),
             linetype = "dashed", linewidth = 1) +
  facet_wrap(~ design, ncol = 1) +
  scale_fill_manual(values  = c("Vertical"   = "#2166ac",
                                "Horizontal" = "#d6604d")) +
  scale_color_manual(values = c("Vertical"   = "#2166ac",
                                "Horizontal" = "#d6604d")) +
  scale_x_continuous(limits = c(0, 1), breaks = seq(0, 1, 0.2)) +
  scale_y_continuous(labels = scales::percent_format(accuracy = 1)) +
  labs(
    title    = "Distribution of Monotonicity Violation Rates by Design",
    subtitle = paste0("p-value (t-test): ", round(ttest$p.value, 3)),
    x        = "Violation rate (per person)",
    y        = "Proportion of respondents"
  ) +
  theme_bw(base_size = 13) +
  theme(
    legend.position  = "none",
    strip.background = element_rect(fill = "gray95"),
    strip.text       = element_text(face = "bold"),
    plot.title       = element_text(face = "bold")
  )

ggsave("monotonicity_violations.pdf", p_viol, width = 7, height = 5)

# -----------------------------------------------------------------------------
# PLOT 2: Overlapping empirical CDF (FOSD check)
# -----------------------------------------------------------------------------

p_cdf <- ggplot(violations, aes(x = viol_rate, color = design)) +
  stat_ecdf(linewidth = 1) +
  scale_color_manual(values = c("Vertical"   = "#2166ac",
                                "Horizontal" = "#d6604d")) +
  scale_x_continuous(limits = c(0, 1), breaks = seq(0, 1, 0.2)) +
  scale_y_continuous(limits = c(0, 1), breaks = seq(0, 1, 0.2)) +
  labs(
    title    = "Empirical CDF of Monotonicity Violation Rates by Design",
    subtitle = "FOSD: the lower-noise design should have its CDF shifted left everywhere",
    x        = "Violation rate (per person)",
    y        = "Cumulative proportion",
    color    = "Design"
  ) +
  theme_bw(base_size = 13) +
  theme(
    legend.position = c(0.8, 0.3),
    plot.title      = element_text(face = "bold")
  )

ggsave("monotonicity_violations_cdf.pdf", p_cdf, width = 7, height = 4.5)