library(dplyr)
library(tidyr)

# Calculate overall correlations
overall_swb_corr <- cor(aspect_data$SWB_level_W1, aspect_data$SWB_level_W2, 
                        use = "complete.obs")
overall_imp_corr <- cor(aspect_data$Importance_W1, aspect_data$Importance_W2, 
                        use = "complete.obs")
overall_cali_corr <- cor(calibration_data$Response_W1, calibration_data$Response_W2, 
                         use = "complete.obs")

# Calculate rankings within each individual for core aspects
aspect_data_with_ranks <- aspect_data %>%
  filter(Core_aspect == TRUE) %>%
  group_by(PROLIFIC_PID) %>%
  mutate(
    Imp_rank_W1 = rank(-Importance_W1, ties.method = "average"),
    Imp_rank_W2 = rank(-Importance_W2, ties.method = "average")
  ) %>%
  ungroup()

# Calculate by-aspect correlations for core aspects
aspect_correlations <- aspect_data_with_ranks %>%
  group_by(Aspect_Code) %>%
  summarise(
    Aspect_name = first(Aspect),
    SWB_corr = cor(SWB_level_W1, SWB_level_W2, use = "complete.obs"),
    Imp_rank_corr = cor(Imp_rank_W1, Imp_rank_W2, use = "complete.obs"),
    Importance_corr = cor(Importance_W1, Importance_W2, use = "complete.obs"),
    .groups = "drop"
  ) %>%
  arrange(Aspect_Code) %>%
  select(Aspect_name, SWB_corr, Imp_rank_corr, Importance_corr)

# Manual LaTeX generation
latex_table <- bind_rows(
  tibble(
    Aspect_name = "Calibration Questions",
    SWB_corr = overall_cali_corr,
    Imp_rank_corr = NA,
    Importance_corr = NA
  ),
  tibble(
    Aspect_name = "Pooled",
    SWB_corr = overall_swb_corr,
    Imp_rank_corr = NA,
    Importance_corr = overall_imp_corr
  ),
  aspect_correlations
) %>%
  mutate(
    SWB_corr = sprintf("%.2f", SWB_corr),
    Imp_rank_corr = ifelse(is.na(Imp_rank_corr), "", sprintf("%.2f", Imp_rank_corr)),
    Importance_corr = ifelse(is.na(Importance_corr), "", sprintf("%.2f", Importance_corr))
  )

# Create the LaTeX output as a string
latex_output <- c(
  "\\begin{table}[htbp]",
  "\\centering",
  "\\caption{Test-Retest Correlations}",
  "\\begin{tabular}{lccc}",
  "\\hline\\hline",
  " & SWB Level & \\multicolumn{2}{c}{Importance} \\\\",
  "\\cmidrule(lr){3-4}",
  "Aspect & Level & Ranking & Level \\\\",
  "\\hline"
)

for(i in 1:nrow(latex_table)) {
  latex_output <- c(latex_output,
                    sprintf("%s & %s & %s & %s \\\\", 
                            latex_table$Aspect_name[i],
                            latex_table$SWB_corr[i],
                            latex_table$Imp_rank_corr[i],
                            latex_table$Importance_corr[i]))
}

latex_output <- c(latex_output,
                  "\\hline\\hline",
                  "\\end{tabular}",
                  "\\end{table}")

# Save to file
writeLines(latex_output, "/Users/liuqianyi/Dropbox/swb/Epinion/graph/test_retest_correlations.tex")

cat("LaTeX table saved to /Users/liuqianyi/Dropbox/swb/Epinion/graph/test_retest_correlations.tex\n")