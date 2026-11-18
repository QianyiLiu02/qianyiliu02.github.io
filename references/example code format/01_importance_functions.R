# -----------------------------
# Define core aspects
# -----------------------------
core_aspects <- c(
  "How happy you feel",
  "How much you enjoy your life",
  "How satisfied you are with your life",
  "The absence of sadness in your life",
  "The extent to which you feel the things you do in your life are worthwhile",
  "The happiness of your family",
  "You being a good person",
  "You being able to support your family financially",
  "You feeling that you have enough time for the things that are most important to you",
  "You not being lonely",
  "You not feeling anxious",
  "Your mental health",
  "Your physical health",
  "Your rating of your life on a ladder where the lowest rung is worst possible life for you and the highest rung is the best possible life for you",
  "Your sense of control over your life",
  "Your sense of purpose"
)

# -----------------------------
# Parsing function
# -----------------------------
aspects_with_hyphen <- c(
  "You feeling like you have work-life balance",
  "Your children growing up to follow the golden rule – treating others as they would like to be treated",
  "You feeling well-rested" 
)

parse_comparison <- function(comp_string, respondent_id) {
  if (is.na(comp_string) || comp_string == "") return(NULL)
  
  comparisons <- str_split(comp_string, "\\|")[[1]]
  comparisons <- comparisons[comparisons != ""]
  parsed_list <- list()
  
  for (i in seq_along(comparisons)) {
    comp <- comparisons[i]
    if (!str_detect(comp, "Option")) next
    
    parts <- str_split(comp, ":")[[1]]
    if (length(parts) < 6) next
    
    aspect_pair <- parts[1]
    aspect_j <- NULL
    aspect_k <- NULL
    
    found_problematic <- FALSE
    for (prob_aspect in aspects_with_hyphen) {
      if (str_detect(aspect_pair, fixed(prob_aspect))) {
        if (str_starts(aspect_pair, fixed(prob_aspect))) {
          aspect_j <- prob_aspect
          remainder <- str_remove(aspect_pair, fixed(paste0(prob_aspect, "-")))
          aspect_k <- remainder
        } else if (str_ends(aspect_pair, fixed(prob_aspect))) {
          aspect_k <- prob_aspect
          remainder <- str_remove(aspect_pair, fixed(paste0("-", prob_aspect)))
          aspect_j <- remainder
        }
        found_problematic <- TRUE
        break
      }
    }
    
    if (!found_problematic) {
      last_hyphen_pos <- str_locate_all(aspect_pair, "-")[[1]]
      if (nrow(last_hyphen_pos) > 0) {
        last_pos <- last_hyphen_pos[nrow(last_hyphen_pos), "start"]
        aspect_j <- str_sub(aspect_pair, 1, last_pos - 1)
        aspect_k <- str_sub(aspect_pair, last_pos + 1)
      }
    }
    
    if (!is.null(aspect_j) && !is.null(aspect_k)) {
      delta_j <- suppressWarnings(as.numeric(parts[2]))
      delta_k <- suppressWarnings(as.numeric(parts[3]))
      
      if (is.na(delta_j) || is.na(delta_k) || delta_j == 0 || delta_k == 0) next
      
      choice_part <- parts[length(parts)]
      choice <- ifelse(str_detect(choice_part, "Option1"), 1,
                       ifelse(str_detect(choice_part, "Option2"), 0, NA))
      
      if (!is.na(choice)) {
        aspect_j <- str_trim(aspect_j)
        aspect_k <- str_trim(aspect_k)
        
        parsed_list[[length(parsed_list) + 1]] <- tibble(
          respondent_id = respondent_id,
          comparison_num = i,
          aspect_j = aspect_j,
          aspect_k = aspect_k,
          delta_j = delta_j,
          delta_k = delta_k,
          choice = choice
        )
      }
    }
  }
  
  if (length(parsed_list) > 0) bind_rows(parsed_list) else NULL
}

# -----------------------------
# Estimation function
# -----------------------------
estimate_model_glm <- function(data_subset) {
  
  stopifnot(all(c("respondent_id","aspect_j","aspect_k","delta_j","delta_k","choice") %in% names(data_subset)))
  
  # Stable ordering so columns map reproducibly to aspects
  aspects_subset <- sort(unique(c(data_subset$aspect_j, data_subset$aspect_k)))
  n_asp <- length(aspects_subset)
  n_obs <- nrow(data_subset)
  
  if (n_asp <= 1L) stop("Need at least two aspects to estimate the model.")
  if (!all(data_subset$choice %in% c(0,1))) stop("`choice` must be 0/1.")
  
  # one column per aspect except the anchor (first in sorted order)
  X_sparse <- Matrix(0, n_obs, n_asp - 1L, sparse = TRUE)
  
  j_idx <- match(data_subset$aspect_j, aspects_subset)  # 1..n_asp
  k_idx <- match(data_subset$aspect_k, aspects_subset)
  
  rows_j <- which(j_idx > 1L)  # not anchor
  rows_k <- which(k_idx > 1L)
  
  if (length(rows_j)) X_sparse[cbind(rows_j, j_idx[rows_j] - 1L)] <- data_subset$delta_j[rows_j]
  if (length(rows_k)) X_sparse[cbind(rows_k, k_idx[rows_k] - 1L)] <- -data_subset$delta_k[rows_k]
  
  # Convert to dense only if manageable.
  X_dense <- as.matrix(X_sparse)
  colnames(X_dense) <- aspects_subset[-1L]  # columns labeled by original aspect names (minus anchor)
  
  # Build model frame
  dat <- cbind(
    data.frame(
      choice        = data_subset$choice,
      respondent_id = data_subset$respondent_id,
      stringsAsFactors = FALSE
    ),
    as.data.frame(X_dense, check.names = TRUE)
  )
  
  # Ensure all names are syntactic & unique (also affects choice/respondent_id)
  names(dat) <- make.names(names(dat), unique = TRUE)
  
  # Build formula from the names actually present in the data frame
  rhs_vars <- setdiff(names(dat), c("choice", "respondent_id"))
  if (length(rhs_vars) == 0L) stop("No regressors (only anchor present).")
  fml <- reformulate(rhs_vars, response = "choice")
  fml <- update(fml, . ~ . - 1)  # no intercept
  model_glm <- glm(fml, data = dat, family = binomial(link = "probit"))
  
  # Clustered VCOV at respondent level; using model frame ensures perfect alignment
  vc <- sandwich::vcovCL(model_glm, cluster = ~ respondent_id, data = dat, type = "HC0")
  se_cluster <- sqrt(diag(vc))
  
  # Map coefficients back to original aspect labels (including anchor)
  syntactic_names <- rhs_vars
  original_names  <- aspects_subset[-1L]
  name_map <- setNames(original_names, syntactic_names)
  
  # Collect estimates in original-name order
  beta_hat <- coef(model_glm)
  beta_se  <- se_cluster
  
  # Reorder to match rhs_vars (already matches), then relabel
  beta_hat <- beta_hat[syntactic_names]
  beta_se  <-  beta_se[syntactic_names]
  names(beta_hat) <- names(beta_se) <- name_map[names(beta_hat)]
  
  # Add anchor (beta_1 = 0; SE = NA)
  MU <- c(setNames(0, aspects_subset[1L]), beta_hat)
  SE <- c(setNames(NA_real_, aspects_subset[1L]), beta_se)
  
  # Ensure names are exactly the original aspect labels for all entries
  MU <- MU[aspects_subset]
  SE <- SE[aspects_subset]
  
  list(
    MU     = MU,
    SE     = SE,
    model  = model_glm,
    vcov   = vc,
    anchor = aspects_subset[1L],
    method = "glm_probit_cluster_by_respondent"
  )
}

# Put this *below* your existing helpers (core_aspects, parse_comparison, estimate_model_glm)

stratified_importance_vs_mu <- function(csv_path, var_name, aspects = core_aspects,
                                        save_plot = TRUE, width = 10, height = 8, dpi = 300) {
  library(tidyverse)
  library(Matrix)
  library(sandwich)
  library(lmtest)
  library(ggplot2)

  # ---------- load ----------
  df <- read.csv(csv_path, check.names = FALSE, stringsAsFactors = FALSE)

  # get college_plus variable 
  df$college_plus <- ifelse(df$education %in% c("Associate's degree (AA, AS)", 
                                             "Bachelor's degree (BA, BS)", 
                                             "Doctoral degree (PhD, EdD)", 
                                             "Master's degree (MA, MS, MBA)", 
                                             "Professional degree (JD, MD)", 
                                             "Some college, no degree"), 1, 0)

  df$above_median_duration <- as.numeric(!is.na(as.numeric(df$`Duration (in seconds)`)) & as.numeric(df$`Duration (in seconds)`) > median(as.numeric(df$`Duration (in seconds)`), na.rm = TRUE))
  df$female <- as.numeric(df$gender == "Female" & !is.na(df$gender))
  df$less_than_40 <- ifelse(as.numeric(df$age) < 40 & !is.na(as.numeric(df$age)), 1, 0)
  df$income_less_than_70k <- ifelse(as.numeric(df$income) < 70000 & !is.na(as.numeric(df$income)), 1, 0)

  if (!all(c("ResponseId", "allComparisons", var_name) %in% names(df))) {
    stop("csv must include ResponseId, allComparisons, and the binary variable: ", var_name)
  }

  # generate education variable 
  
  # coerce group to numeric 0/1
  grp <- suppressWarnings(as.numeric(df[[var_name]]))
  if (!all(grp %in% c(0, 1), na.rm = TRUE)) stop(var_name, " must be binary 0/1.")
  df$..group.. <- grp

  # ---------- build pairwise-comparison rows w/ group ----------
  processed_data <- purrr::map_df(seq_len(nrow(df)), function(i) {
    if (is.na(df$allComparisons[i]) || df$allComparisons[i] == "") return(NULL)
    ac <- parse_comparison(df$allComparisons[i], df$ResponseId[i])
    if (is.null(ac) || nrow(ac) == 0) return(NULL)
    ac$group <- df$..group..[i]
    ac %>%
      group_by(respondent_id, aspect_j, aspect_k, delta_j, delta_k, group) %>%
      slice_tail(n = 1) %>%
      ungroup()
  }) %>% 
    filter(!is.na(choice), !is.na(group))

  processed_data_core <- processed_data %>%
    filter(aspect_j %in% aspects & aspect_k %in% aspects)

  # ---------- estimate MU by group (0,1) ----------
  groups <- sort(intersect(c(0, 1), unique(processed_data_core$group)))

  safe_est <- function(d) tryCatch(estimate_model_glm(d), error = function(e) NULL)

  mu_by_group <- purrr::map_dfr(groups, function(g) {
    est <- safe_est(processed_data_core %>% filter(group == g))
    if (is.null(est)) return(tibble())
    tibble(
      aspect = names(est$MU),
      log_MU = as.numeric(est$MU),
      SE_MU  = as.numeric(est$SE),
      group  = g
    )
  }) %>% filter(aspect %in% aspects, !is.na(log_MU))

  if (nrow(mu_by_group) == 0) stop("No MU estimates produced. Check data and filters.")

  # ---------- extract importance w/ group ----------
  extract_importance_by_group <- function(df_full) {
    shuffled_cols   <- paste0("shuffled_statement_", 1:16)
    importance_cols <- paste0("core_importance_", 1:16)
    stopifnot(all(shuffled_cols %in% names(df_full)),
              all(importance_cols %in% names(df_full)))

    res <- vector("list", nrow(df_full))
    for (i in seq_len(nrow(df_full))) {
      st  <- df_full[i, shuffled_cols, drop = FALSE]
      imp <- df_full[i, importance_cols, drop = FALSE]
      rr  <- vector("list", 16)
      for (j in 1:16) {
        aspect_text <- as.character(st[[j]])
        importance_value <- suppressWarnings(as.numeric(imp[[j]]))
        if (!is.na(aspect_text) && nzchar(aspect_text) &&
            !is.na(importance_value) && aspect_text %in% aspects) {
          rr[[j]] <- tibble(
            respondent_id    = df_full$ResponseId[i],
            aspect           = aspect_text,
            importance_rating = importance_value,
            group            = df_full$..group..[i]
          )
        }
      }
      res[[i]] <- bind_rows(rr)
    }
    bind_rows(res) %>% filter(!is.na(group))
  }

  importance_df <- extract_importance_by_group(df)

  avg_importance <- importance_df %>%
    group_by(group, aspect) %>%
    summarise(
      mean_importance = mean(importance_rating, na.rm = TRUE),
      sd_importance   = sd(importance_rating, na.rm = TRUE),
      n_ratings       = dplyr::n(),
      se_importance   = sd_importance / sqrt(n_ratings),
      .groups = "drop"
    )

  # ---------- merge ----------
  comparison_df <- avg_importance %>%
    inner_join(mu_by_group, by = c("group", "aspect")) %>%
    mutate(
      MU = exp(log_MU)  # keep for reference
    ) %>%
    group_by(group) %>%
    mutate(
      importance_z = (mean_importance - mean(mean_importance, na.rm = TRUE)) /
                     sd(mean_importance,  na.rm = TRUE)
    ) %>%
    ungroup()

  # correlations and labels
  cors <- comparison_df %>%
    group_by(group) %>%
    summarise(r = cor(mean_importance, log_MU, use = "complete.obs"), .groups = "drop")

  lbls <- setNames(
    sprintf("%s = %d (rho = %.3f)", var_name, cors$group, cors$r),
    cors$group
  )

  comparison_df <- comparison_df %>%
    mutate(group_label = factor(group, levels = cors$group, labels = unname(lbls)))


  # ---------- plot ----------
p <- ggplot(comparison_df, aes(x = importance_z, y = log_MU)) +
   geom_point(size = 3, alpha = 0.7) +
   geom_smooth(method = "lm", se = TRUE, alpha = 0.2) +
   geom_text(aes(label = substr(aspect, 1, 30)),
             vjust = -0.5, hjust = 0, size = 2.5, check_overlap = TRUE) +
   facet_wrap(~ group_label) +
   labs(
     title = sprintf("Importance vs. Log Marginal Utility by %s", var_name),
     x = "Standardized Importance (within group)",
     y = "Log Marginal Utility (from tradeoffs)"
   ) +
   theme_minimal(base_size = 12) +
   theme(
     plot.title   = element_text(size = 14, face = "bold"),
     axis.title   = element_text(size = 11),
     strip.text   = element_text(face = "bold"),
     panel.background = element_rect(fill = "white", color = NA),
     plot.background = element_rect(fill = "white", color = NA)
   )
   
  file_out <- paste0("importance_vs_logMU_by_", var_name, ".png")
  if (save_plot) ggsave(file_out, p, width = width, height = height, dpi = dpi)

  # also return per-group OLS if you want the slopes
  lm_by_group <- purrr::map(setNames(groups, groups), function(g) {
    lm(log_MU ~ mean_importance, data = comparison_df %>% filter(group == g))
  })

  list(
    comparison = comparison_df,
    correlations = cors,
    lm_models = lm_by_group,
    plot = p,
    file = if (save_plot) file_out else NULL
  )
}

extract_importance_data <- function(df) {
  # Get the shuffled statements to match with importance ratings
  shuffled_cols <- paste0("shuffled_statement_", 1:16)
  importance_cols <- paste0("core_importance_", 1:16)
  
  # Initialize results list
  importance_data <- list()
  
  for (i in 1:nrow(df)) {
    prolific_id <- df$prolific_ID[i]
    
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
          prolific_id = prolific_id,
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

importance_predictive_power_fast <- function(csv_path, aspects = core_aspects) {
  library(tidyverse)
  library(Matrix)
  
  processed_data <- read.csv(csv_path, check.names = FALSE, stringsAsFactors = FALSE)
  
  # Estimate log(MU) for all aspects using full data (no importance splits)
  full_data_mu <- tryCatch({
    est <- estimate_model_glm(processed_data)
    tibble(
      aspect = names(est$MU),
      log_MU = as.numeric(est$MU)
    ) %>%
      filter(aspect %in% aspects)
  }, error = function(e) tibble())
  
  # Calculate pooled standard deviation from full data estimates
  pooled_sd_log_mu <- sd(full_data_mu$log_MU, na.rm = TRUE)
  
  # Second pass: calculate final results with pooled SD
  results <- map_dfr(aspects, function(asp) {
    
    # Get data where this aspect appears
    asp_data <- processed_data %>%
      filter(aspect_j == asp | aspect_k == asp) %>%
      filter(!is.na(importance_j) | !is.na(importance_k))
    
    if(nrow(asp_data) == 0) return(tibble())
    
    # Create importance variable for this aspect
    asp_data <- asp_data %>%
      mutate(
        focal_importance = case_when(
          aspect_j == asp ~ importance_j,
          aspect_k == asp ~ importance_k,
          TRUE ~ NA_real_
        )
      ) %>%
      filter(!is.na(focal_importance))
    
    # Get median importance for this aspect
    median_imp <- median(asp_data$focal_importance, na.rm = TRUE)
    
    # Split into above/below median
    above_med <- asp_data %>%
      filter(focal_importance > median_imp)
    
    below_med <- asp_data %>%
      filter(focal_importance <= median_imp)
    
    # Estimate MU for each group
    mu_results <- map_dfr(list(below = below_med, above = above_med), function(data) {
      if(nrow(data) < 10) return(tibble())  # Need minimum observations
      
      tryCatch({
        est <- estimate_model_glm(data)
        if(asp %in% names(est$MU)) {
          tibble(
            log_MU = as.numeric(est$MU[asp]),
            SE_MU = as.numeric(est$SE[asp]),
            n_obs = nrow(data)
          )
        } else tibble()
      }, error = function(e) tibble())
    }, .id = "group")
    
    if(nrow(mu_results) == 2) {
      mu_high <- mu_results$log_MU[mu_results$group == "above"]
      mu_low <- mu_results$log_MU[mu_results$group == "below"]
      mu_diff <- mu_high - mu_low
      
      # Calculate importance difference between groups
      imp_high <- mean(above_med$focal_importance, na.rm = TRUE)
      imp_low <- mean(below_med$focal_importance, na.rm = TRUE)
      imp_diff <- imp_high - imp_low
      
      tibble(
        aspect = asp,
        log_MU_high = mu_high,
        log_MU_low = mu_low,
        log_MU_diff = mu_diff,
        importance_high = imp_high,
        importance_low = imp_low,
        importance_diff = imp_diff,
        predictive_ratio = mu_diff / pooled_sd_log_mu,
        n_high = mu_results$n_obs[mu_results$group == "above"],
        n_low = mu_results$n_obs[mu_results$group == "below"]
      )
    } else {
      tibble()
    }
  })
  
  # Summary
  summary_stats <- results %>%
    filter(!is.na(predictive_ratio)) %>%
    summarise(
      n_aspects = n(),
      pooled_sd_log_mu = pooled_sd_log_mu,
      mean_predictive_ratio = mean(predictive_ratio, na.rm = TRUE),
      median_predictive_ratio = median(predictive_ratio, na.rm = TRUE),
      positive_ratios = sum(predictive_ratio > 0, na.rm = TRUE),
      pct_positive = positive_ratios / n_aspects * 100
    )
  
  list(
    results = results,
    summary = summary_stats
  )
}

simple_demand_curve <- function(csv_path, focal_aspect, split_by = "importance", n_bins = 10) {
  library(tidyverse)
  library(ggplot2)
  
  # Load data
  df <- read.csv(csv_path, check.names = FALSE, stringsAsFactors = FALSE)
  
  if (split_by == "importance" || split_by == "importance_quartile") {
    # Extract importance for focal aspect only - simple approach
    shuffled_cols   <- paste0("shuffled_statement_", 1:16)
    importance_cols <- paste0("core_importance_", 1:16)
    
    df$focal_importance <- NA_real_
    
    for (i in seq_len(nrow(df))) {
      for (j in 1:16) {
        aspect_text <- as.character(df[i, shuffled_cols[j]])
        importance_value <- suppressWarnings(as.numeric(df[i, importance_cols[j]]))
        if (!is.na(aspect_text) && nzchar(aspect_text) && 
            aspect_text == focal_aspect && !is.na(importance_value)) {
          df$focal_importance[i] <- importance_value
          break
        }
      }
    }
    
    # Show distribution
    cat("Importance distribution:\n")
    print(table(df$focal_importance, useNA = "always"))
    cat("Mean:", mean(df$focal_importance, na.rm = TRUE), "\n")
    cat("Median:", median(df$focal_importance, na.rm = TRUE), "\n")
    
    # Create splits
    median_importance <- median(df$focal_importance, na.rm = TRUE)
    q25 <- quantile(df$focal_importance, 0.25, na.rm = TRUE)
    q75 <- quantile(df$focal_importance, 0.75, na.rm = TRUE)
    
    if (split_by == "importance_quartile") {
      df$split_var <- ifelse(!is.na(df$focal_importance) & df$focal_importance >= q75, 1,
                            ifelse(!is.na(df$focal_importance) & df$focal_importance <= q25, 0, NA))
      split_label <- "Top vs Bottom Quartile Importance"
      cat("Using quartile split: <=", q25, "vs >=", q75, "\n")
    } else {
      df$split_var <- ifelse(!is.na(df$focal_importance) & df$focal_importance > median_importance, 1, 0)
      split_label <- "High Importance"
      cat("Using median split: >", median_importance, "\n")
    }
    
  } else {
    # Handle characteristic variables
    if (split_by == "income_low") {
      income_numeric <- suppressWarnings(as.numeric(as.character(df$income)))
      df$split_var <- ifelse(!is.na(income_numeric) & income_numeric < 70000, 1, 0)
      split_label <- "Income < 70k"
    } else if (split_by == "college_plus") {
      df$split_var <- ifelse(df$education %in% c("Associate's degree (AA, AS)", 
                                                "Bachelor's degree (BA, BS)", 
                                                "Doctoral degree (PhD, EdD)", 
                                                "Master's degree (MA, MS, MBA)", 
                                                "Professional degree (JD, MD)", 
                                                "Some college, no degree"), 1, 0)
      split_label <- "College+"
    } else if (split_by == "female") {
      df$split_var <- as.numeric(df$gender == "Female" & !is.na(df$gender))
      split_label <- "Female"
    } else {
      stop("Unknown split_by option: ", split_by)
    }
  }
  
  # Parse all comparisons
  processed_data <- map_df(seq_len(nrow(df)), function(i) {
    if (is.na(df$allComparisons[i]) || df$allComparisons[i] == "") return(NULL)
    ac <- parse_comparison(df$allComparisons[i], df$ResponseId[i])
    if (is.null(ac) || nrow(ac) == 0) return(NULL)
    ac$split_var <- df$split_var[i]
    ac
  }) %>% 
    filter(!is.na(choice), !is.na(split_var))
  
  cat("Total comparisons:", nrow(processed_data), "\n")
  
  # Filter to focal aspect
  focal_data <- processed_data %>%
    filter(aspect_j == focal_aspect | aspect_k == focal_aspect) %>%
    mutate(
      qi = case_when(
        aspect_j == focal_aspect & choice == 1 ~ 1,
        aspect_k == focal_aspect & choice == 0 ~ 1,
        TRUE ~ 0
      ),
      pi = case_when(
        aspect_j == focal_aspect ~ delta_k - delta_j,
        aspect_k == focal_aspect ~ delta_j - delta_k,
        TRUE ~ NA_real_
      ),
      focal_position = ifelse(aspect_j == focal_aspect, "j", "k")
    ) %>%
    filter(!is.na(pi))
  
  cat("Focal aspect comparisons:", nrow(focal_data), "\n")
  cat("Split var = 0:", sum(focal_data$split_var == 0), "\n")
  cat("Split var = 1:", sum(focal_data$split_var == 1), "\n")
  
  if (split_by == "importance" || split_by == "importance_quartile") {
    # Check average choice rates by group
    choice_by_group <- focal_data %>%
      group_by(split_var) %>%
      summarise(
        mean_choice = mean(qi),
        mean_price = mean(pi),
        .groups = "drop"
      )
    cat("Choice rates by group:\n")
    print(choice_by_group)
  }
  
  if (nrow(focal_data) == 0) {
    stop("No data for focal aspect")
  }
  
  # Simple binning (n_bins per group)
  binned_data <- focal_data %>%
    group_by(split_var) %>%
    mutate(price_bin = ntile(pi, n_bins)) %>%
    group_by(split_var, price_bin) %>%
    summarise(
      mean_qi = mean(qi),
      mean_pi = mean(pi),
      n_obs = n(),
      .groups = "drop"
    )
  
  cat("Binned data rows:", nrow(binned_data), "\n")
  cat("Average obs per bin:", mean(binned_data$n_obs), "\n")
  
  # Plot
  p <- ggplot(binned_data, aes(x = mean_qi, y = mean_pi, color = factor(split_var))) +
    geom_point(size = 2) +
    geom_smooth(method = "lm", se = FALSE) +
    labs(
      title = paste("Demand Curve for:", substr(focal_aspect, 1, 40), "..."),
      x = "Proportion choosing focal aspect",
      y = "Price",
      color = split_label
    ) +
    theme_minimal()
  
  return(list(plot = p, binned_data = binned_data, focal_data = focal_data))
}

simple_regression_test <- function(csv_path, focal_aspect) {
  library(tidyverse)
  library(broom)
  
  # Load data
  df <- read.csv(csv_path, check.names = FALSE, stringsAsFactors = FALSE)
  
  # Extract importance for focal aspect
  shuffled_cols   <- paste0("shuffled_statement_", 1:16)
  importance_cols <- paste0("core_importance_", 1:16)
  
  df$focal_importance <- NA_real_
  
  for (i in seq_len(nrow(df))) {
    for (j in 1:16) {
      aspect_text <- as.character(df[i, shuffled_cols[j]])
      importance_value <- suppressWarnings(as.numeric(df[i, importance_cols[j]]))
      if (!is.na(aspect_text) && nzchar(aspect_text) && 
          aspect_text == focal_aspect && !is.na(importance_value)) {
        df$focal_importance[i] <- importance_value
        break
      }
    }
  }
  
  # Parse comparisons
  processed_data <- map_df(seq_len(nrow(df)), function(i) {
    if (is.na(df$allComparisons[i]) || df$allComparisons[i] == "") return(NULL)
    ac <- parse_comparison(df$allComparisons[i], df$ResponseId[i])
    if (is.null(ac) || nrow(ac) == 0) return(NULL)
    ac$focal_importance <- df$focal_importance[i]
    ac
  }) %>% 
    filter(!is.na(choice))
  
  # Create focal aspect data
  focal_data <- processed_data %>%
    filter(aspect_j == focal_aspect | aspect_k == focal_aspect) %>%
    mutate(
      qi = case_when(
        aspect_j == focal_aspect & choice == 1 ~ 1,
        aspect_k == focal_aspect & choice == 0 ~ 1,
        TRUE ~ 0
      ),
      pi = case_when(
        aspect_j == focal_aspect ~ delta_k - delta_j,
        aspect_k == focal_aspect ~ delta_j - delta_k,
        TRUE ~ NA_real_
      )
    ) %>%
    filter(!is.na(pi), !is.na(focal_importance))
  
  cat("N observations:", nrow(focal_data), "\n")
  cat("Mean choice rate:", mean(focal_data$qi), "\n")
  cat("Mean price:", mean(focal_data$pi), "\n")
  cat("Mean importance:", mean(focal_data$focal_importance), "\n")
  
  # Run regression: qi ~ importance + price
  reg1 <- lm(qi ~ focal_importance + pi, data = focal_data)
  
  cat("\n=== Regression: qi ~ importance + price ===\n")
  print(summary(reg1))
  
  # Also run without price control for comparison
  reg2 <- lm(qi ~ focal_importance, data = focal_data)
  
  cat("\n=== Regression: qi ~ importance (no price) ===\n")  
  print(summary(reg2))
  
  # Correlation for reference
  cor_qi_imp <- cor(focal_data$qi, focal_data$focal_importance)
  cor_qi_price <- cor(focal_data$qi, focal_data$pi)
  cor_imp_price <- cor(focal_data$focal_importance, focal_data$pi)
  
  cat("\n=== Correlations ===\n")
  cat("qi ~ importance:", round(cor_qi_imp, 4), "\n")
  cat("qi ~ price:", round(cor_qi_price, 4), "\n") 
  cat("importance ~ price:", round(cor_imp_price, 4), "\n")
  
  return(list(
    data = focal_data,
    reg_with_price = reg1,
    reg_no_price = reg2,
    correlations = c(qi_imp = cor_qi_imp, qi_price = cor_qi_price, imp_price = cor_imp_price)
  ))
}

# debugging of importance 
debug_importance_detailed <- function(csv_path, focal_aspect) {
  library(tidyverse)
  
  df <- read.csv(csv_path, check.names = FALSE, stringsAsFactors = FALSE)
  
  # Method 1: Your working approach from importance_predictive_power_fast
  extract_importance_working <- function(df_full) {
    shuffled_cols   <- paste0("shuffled_statement_", 1:16)
    importance_cols <- paste0("core_importance_", 1:16)
    
    res <- vector("list", nrow(df_full))
    for (i in seq_len(nrow(df_full))) {
      st  <- df_full[i, shuffled_cols, drop = FALSE]
      imp <- df_full[i, importance_cols, drop = FALSE]
      rr  <- vector("list", 16)
      for (j in 1:16) {
        aspect_text <- as.character(st[[j]])
        importance_value <- suppressWarnings(as.numeric(imp[[j]]))
        if (!is.na(aspect_text) && nzchar(aspect_text) &&
            !is.na(importance_value) && aspect_text %in% core_aspects) {
          rr[[j]] <- tibble(
            respondent_id = df_full$ResponseId[i],
            aspect = aspect_text,
            importance = importance_value
          )
        }
      }
      res[[i]] <- bind_rows(rr)
    }
    bind_rows(res)
  }
  
  # Method 2: My simple approach
  extract_importance_simple <- function(df_full, target_aspect) {
    shuffled_cols   <- paste0("shuffled_statement_", 1:16)
    importance_cols <- paste0("core_importance_", 1:16)
    
    df_full$focal_importance <- NA_real_
    
    for (i in seq_len(nrow(df_full))) {
      for (j in 1:16) {
        aspect_text <- as.character(df_full[i, shuffled_cols[j]])
        importance_value <- suppressWarnings(as.numeric(df_full[i, importance_cols[j]]))
        if (!is.na(aspect_text) && nzchar(aspect_text) && 
            aspect_text == target_aspect && !is.na(importance_value)) {
          df_full$focal_importance[i] <- importance_value
          break
        }
      }
    }
    return(df_full$focal_importance)
  }
  
  # Extract using both methods
  cat("=== METHOD 1: Working approach ===\n")
  importance_df_working <- extract_importance_working(df)
  focal_working <- importance_df_working %>% filter(aspect == focal_aspect)
  cat("Found importance for", nrow(focal_working), "people\n")
  cat("Range:", range(focal_working$importance, na.rm = TRUE), "\n")
  cat("Mean:", mean(focal_working$importance, na.rm = TRUE), "\n")
  
  cat("\n=== METHOD 2: Simple approach ===\n")
  importance_simple <- extract_importance_simple(df, focal_aspect)
  cat("Found importance for", sum(!is.na(importance_simple)), "people\n")
  cat("Range:", range(importance_simple, na.rm = TRUE), "\n")
  cat("Mean:", mean(importance_simple, na.rm = TRUE), "\n")
  
  # Compare the two
  df$imp_working <- NA_real_
  df$imp_working[match(focal_working$respondent_id, df$ResponseId)] <- focal_working$importance
  df$imp_simple <- importance_simple
  
  comparison <- df %>%
    select(ResponseId, imp_working, imp_simple) %>%
    filter(!is.na(imp_working) | !is.na(imp_simple))
  
  cat("\n=== COMPARISON ===\n")
  cat("Both methods agree:", sum(!is.na(comparison$imp_working) & !is.na(comparison$imp_simple) & 
                                  comparison$imp_working == comparison$imp_simple, na.rm = TRUE), "\n")
  cat("Only working method:", sum(!is.na(comparison$imp_working) & is.na(comparison$imp_simple)), "\n")
  cat("Only simple method:", sum(is.na(comparison$imp_working) & !is.na(comparison$imp_simple)), "\n")
  
  # Show discrepancies
  discrepancies <- comparison %>%
    filter((!is.na(imp_working) & is.na(imp_simple)) | 
           (is.na(imp_working) & !is.na(imp_simple)) |
           (!is.na(imp_working) & !is.na(imp_simple) & imp_working != imp_simple))
  
  if (nrow(discrepancies) > 0) {
    cat("\n=== DISCREPANCIES (first 10) ===\n")
    print(head(discrepancies, 10))
  }
  
  return(list(
    working_method = focal_working,
    simple_method = data.frame(ResponseId = df$ResponseId, importance = importance_simple),
    comparison = comparison
  ))
}

export_choice_level_data <- function(csv_path, output_file = "/Users/liuqianyi/Dropbox/swb/Epinion/raw/choice_level_data_20250820.csv") {
  library(tidyverse)
  
  df <- read.csv(csv_path, check.names = FALSE, stringsAsFactors = FALSE)
  
  # Create demographic variables
  df$college_plus <- ifelse(df$education %in% c("Associate's degree (AA, AS)", 
                                               "Bachelor's degree (BA, BS)", 
                                               "Doctoral degree (PhD, EdD)", 
                                               "Master's degree (MA, MS, MBA)", 
                                               "Professional degree (JD, MD)", 
                                               "Some college, no degree"), 1, 0)
  
  duration_numeric <- suppressWarnings(as.numeric(df$`Duration (in seconds)`))
  df$above_median_duration <- as.numeric(!is.na(duration_numeric) & 
                                        duration_numeric > median(duration_numeric, na.rm = TRUE))
  
  df$female <- as.numeric(df$gender == "Female" & !is.na(df$gender))
  
  age_numeric <- suppressWarnings(as.numeric(df$age))
  df$less_than_40 <- ifelse(!is.na(age_numeric) & age_numeric < 40, 1, 0)
  
  income_numeric <- suppressWarnings(as.numeric(df$income))
  df$income_less_than_70k <- ifelse(!is.na(income_numeric) & income_numeric < 70000, 1, 0)
  
  # Extract importance ratings for all respondents
  extract_importance <- function(df_full) {
    shuffled_cols   <- paste0("shuffled_statement_", 1:16)
    importance_cols <- paste0("core_importance_", 1:16)
    
    res <- vector("list", nrow(df_full))
    for (i in seq_len(nrow(df_full))) {
      st  <- df_full[i, shuffled_cols, drop = FALSE]
      imp <- df_full[i, importance_cols, drop = FALSE]
      rr  <- vector("list", 16)
      for (j in 1:16) {
        aspect_text <- as.character(st[[j]])
        importance_value <- suppressWarnings(as.numeric(imp[[j]]))
        if (!is.na(aspect_text) && nzchar(aspect_text) && !is.na(importance_value)) {
          rr[[j]] <- data.frame(
            respondent_id = df_full$ResponseId[i],
            aspect = aspect_text,
            importance = importance_value,
            stringsAsFactors = FALSE
          )
        }
      }
      res[[i]] <- do.call(rbind, rr)
    }
    do.call(rbind, res)
  }
  
  # Read the data for importance ratings
  df_importance <- read.csv("/Users/liuqianyi/Dropbox/swb/Epinion/raw/SWB wave 1 (Nov 4)_November 10, 2025.csv")
  
  # Extract importance data
  importance_df <- extract_importance_data(df_importance)
  
  df1 <- read.csv("/Users/liuqianyi/Dropbox/swb/Epinion/raw/SWB wave 2 (Nov 4)_November 10, 2025.csv", check.names = FALSE, stringsAsFactors = FALSE)
  
  # Parse all comparisons into choice-level data
  choice_data <- map_df(seq_len(nrow(df1)), function(i) {
    if (is.na(df1$allComparisons[i]) || df1$allComparisons[i] == "") return(NULL)
    
    ac <- parse_comparison(df1$allComparisons[i], df1$ResponseId[i])
    if (is.null(ac) || nrow(ac) == 0) return(NULL)
    
    # Add demographic variables
    ac$college_plus <- df$college_plus[i]
    ac$above_median_duration <- df$above_median_duration[i]
    ac$female <- df$female[i]
    ac$less_than_40 <- df$less_than_40[i]
    ac$income_less_than_70k <- df$income_less_than_70k[i]
    ac$education <- df$education[i]
    ac$gender <- df$gender[i]
    ac$age <- df$age[i]
    ac$income <- df$income[i]
    
    return(ac)
  }) %>%
    filter(!is.na(choice)) %>%
    # Filter to only core aspects
    filter(aspect_j %in% core_aspects & aspect_k %in% core_aspects) %>%
  
  
  # Add importance ratings - first for aspect_j

  choice_data$importance_j <- 0
  choice_data$importance_k <- 0
  
  
  for (i in 1:nrow(choice_data)) {
    rid <- choice_data$prolific_ID[i]
    asp_j <- choice_data$aspect_j[i]
    asp_k <- choice_data$aspect_k[i]
    
    # Find importance for aspect_j
    imp_j_match <- which(importance_df$prolific_id == rid & importance_df$aspect == asp_j)
    choice_data$importance_j[i] <- importance_df$importance_rating[imp_j_match[1]]
    
    # Find importance for aspect_k
    imp_k_match <- which(importance_df$prolific_id == rid & importance_df$aspect == asp_k)
    choice_data$importance_k[i] <- importance_df$importance_rating[imp_k_match[1]]
  }
  
  # Write to CSV

  write.csv(choice_data, output_file, row.names = FALSE)
  
  return(choice_data)
}