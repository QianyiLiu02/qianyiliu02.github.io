library(tidyverse)
library(stringr)

# Define core aspects
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

# Create aspect code mapping
aspect_code_mapping <- data.frame(
  Aspect = core_aspects,
  Aspect_Code = 1:16,
  stringsAsFactors = FALSE
)

# ========================================
# FUNCTION 1: Process respondent-aspect level data
# ========================================
process_respondent_aspect_data <- function(wave1_file, wave2_file) {
  
  # Read Wave 1 data
  cat("Reading Wave 1 data from:", wave1_file, "\n")
  df_w1 <- read.csv(wave1_file, stringsAsFactors = FALSE)
  df_w1 <- df_w1[-1, ]  # Drop first row (question text)
  cat("Wave 1 loaded:", nrow(df_w1), "rows\n\n")
  
  # Read Wave 2 data  
  cat("Reading Wave 2 data from:", wave2_file, "\n")
  df_w2 <- read.csv(wave2_file, stringsAsFactors = FALSE)
  df_w2 <- df_w2[-1, ]  # Drop first row (question text)
  cat("Wave 2 loaded:", nrow(df_w2), "rows\n\n")
  
  # ========================================
  # WAVE 1 PROCESSING
  # ========================================
  
  
  # Initialize list for Wave 1 aspect-level data
  w1_aspect_data <- list()
  
  for (i in 1:nrow(df_w1)) {
    row_data <- df_w1[i, ]
    response_id <- row_data$ResponseId
    prolific_id <- row_data$PROLIFIC_PID
    
    # Skip if ResponseId is empty or invalid
    if (is.na(response_id) || response_id == "" || 
        grepl("ImportId|_recordId|Response ID", response_id)) {
      next
    }
    
    # Extract SWB levels and importance for 16 core aspects
    for (j in 1:16) {
      # Get the aspect name from shuffled_statement
      aspect_col <- paste0("shuffled_statement_", j)
      aspect_name <- row_data[[aspect_col]]
      
      # Only process core aspects
      if (is.na(aspect_name) || !(aspect_name %in% core_aspects)) {
        next
      }
      
      # Get SWB level
      swb_col <- paste0("swb_item_", j)
      swb_value <- as.numeric(row_data[[swb_col]])
      
      # Get importance rating
      imp_col <- paste0("core_importance_", j)
      imp_value <- as.numeric(row_data[[imp_col]])
      
      # Only include if BOTH swb and importance are non-missing
      if (!is.na(swb_value) && !is.na(imp_value)) {
        aspect_row <- data.frame(
          ResponseId = response_id,
          PROLIFIC_PID = prolific_id,
          Wave = 1,
          Aspect = aspect_name,
          SWB_level = swb_value,
          Importance = imp_value,
          stringsAsFactors = FALSE
        )
        w1_aspect_data[[length(w1_aspect_data) + 1]] <- aspect_row
      }
    }
    
    # Also process supplementary aspects (16 more)
    for (j in 1:16) {
      # Get the aspect name from supplementary_statement
      aspect_col <- paste0("supplementary_statement_", j)
      aspect_name <- row_data[[aspect_col]]
      
      # Skip if empty
      if (is.na(aspect_name) || aspect_name == "") {
        next
      }
      
      # Get SWB level
      swb_col <- paste0("supplementary_swb_item_", j)
      swb_value <- as.numeric(row_data[[swb_col]])
      
      # Get importance rating
      imp_col <- paste0("supplementary_importance_", j)
      imp_value <- as.numeric(row_data[[imp_col]])
      
      # Only include if BOTH swb and importance are non-missing
      if (!is.na(swb_value) && !is.na(imp_value)) {
        aspect_row <- data.frame(
          ResponseId = response_id,
          PROLIFIC_PID = prolific_id,
          Wave = 1,
          Aspect = aspect_name,
          SWB_level = swb_value,
          Importance = imp_value,
          stringsAsFactors = FALSE
        )
        w1_aspect_data[[length(w1_aspect_data) + 1]] <- aspect_row
      }
    }
  }
  
  w1_panel <- bind_rows(w1_aspect_data)
  cat("Wave 1 aspect-level data:", nrow(w1_panel), "rows\n")
  cat("  - Unique respondents:", length(unique(w1_panel$ResponseId)), "\n\n")
  
  # ========================================
  # WAVE 2 PROCESSING
  # ========================================
  
  
  # Initialize list for Wave 2 aspect-level data
  w2_aspect_data <- list()
  
  for (i in 1:nrow(df_w2)) {
    row_data <- df_w2[i, ]
    response_id <- row_data$ResponseId
    prolific_id <- row_data$PROLIFIC_PID
    
    # Skip if ResponseId is empty or invalid
    if (is.na(response_id) || response_id == "" || 
        grepl("ImportId|_recordId|Response ID", response_id)) {
      next
    }
    
    # Extract SWB levels and importance for 16 core aspects
    for (j in 1:16) {
      # Get the aspect name from shuffled_statement
      aspect_col <- paste0("shuffled_statement_", j)
      aspect_name <- row_data[[aspect_col]]
      
      # Only process core aspects
      if (is.na(aspect_name) || !(aspect_name %in% core_aspects)) {
        next
      }
      
      # Get SWB level
      swb_col <- paste0("swb_item_", j)
      swb_value <- as.numeric(row_data[[swb_col]])
      
      # Get importance rating
      imp_col <- paste0("core_importance_", j)
      imp_value <- as.numeric(row_data[[imp_col]])
      
      # Only include if BOTH swb and importance are non-missing
      if (!is.na(swb_value) && !is.na(imp_value)) {
        aspect_row <- data.frame(
          ResponseId = response_id,
          PROLIFIC_PID = prolific_id,
          Wave = 2,
          Aspect = aspect_name,
          SWB_level = swb_value,
          Importance = imp_value,
          stringsAsFactors = FALSE
        )
        w2_aspect_data[[length(w2_aspect_data) + 1]] <- aspect_row
      }
    }
    
    # Also process supplementary aspects (16 more)
    for (j in 1:16) {
      # Get the aspect name from supplementary_statement
      aspect_col <- paste0("supplementary_statement_", j)
      aspect_name <- row_data[[aspect_col]]
      
      # Skip if empty
      if (is.na(aspect_name) || aspect_name == "") {
        next
      }
      
      # Get SWB level
      swb_col <- paste0("supplementary_swb_item_", j)
      swb_value <- as.numeric(row_data[[swb_col]])
      
      # Get importance rating
      imp_col <- paste0("supplementary_importance_", j)
      imp_value <- as.numeric(row_data[[imp_col]])
      
      # Only include if BOTH swb and importance are non-missing
      if (!is.na(swb_value) && !is.na(imp_value)) {
        aspect_row <- data.frame(
          ResponseId = response_id,
          PROLIFIC_PID = prolific_id,
          Wave = 2,
          Aspect = aspect_name,
          SWB_level = swb_value,
          Importance = imp_value,
          stringsAsFactors = FALSE
        )
        w2_aspect_data[[length(w2_aspect_data) + 1]] <- aspect_row
      }
    }
  }
  
  w2_panel <- bind_rows(w2_aspect_data)
  cat("Wave 2 aspect-level data:", nrow(w2_panel), "rows\n")
  cat("  - Unique respondents:", length(unique(w2_panel$ResponseId)), "\n\n")
  
  # ========================================
  # MERGE WAVE 1 AND WAVE 2
  # ========================================
  
  cat("Merging Wave 1 and Wave 2 data...\n")
  
  # Merge by PROLIFIC_PID and Aspect (not ResponseId, since that changes across waves)
  merged_data <- w1_panel %>%
    inner_join(
      w2_panel, 
      by = c("PROLIFIC_PID", "Aspect"),
      suffix = c("_W1", "_W2")
    ) %>%
    select(
      ResponseId_W1,
      ResponseId_W2,
      PROLIFIC_PID,
      Aspect,
      SWB_level_W1,
      SWB_level_W2,
      Importance_W1,
      Importance_W2
    )
  
  # Add aspect codes
  merged_data <- merged_data %>%
    left_join(aspect_code_mapping, by = "Aspect")
  
  # Calculate differences
  merged_data <- merged_data %>%
    mutate(
      SWB_diff = SWB_level_W2 - SWB_level_W1,
      Importance_diff = Importance_W2 - Importance_W1,
      Core_aspect = Aspect %in% core_aspects  # Flag for core vs supplementary
    )
  
  cat("Merged data:\n")
  cat("  - Total rows:", nrow(merged_data), "\n")
  cat("  - Unique respondents:", length(unique(merged_data$PROLIFIC_PID)), "\n")
  cat("  - Unique aspects:", length(unique(merged_data$Aspect)), "\n")
  cat("  - Core aspects:", sum(merged_data$Core_aspect), "rows\n")
  cat("  - Supplementary aspects:", sum(!merged_data$Core_aspect), "rows\n\n")
  
  return(merged_data)
}

# ========================================
# FUNCTION 2: Process calibration questions
# ========================================
process_calibration_data <- function(wave1_file, wave2_file) {
  
  cat("========================================\n")
  cat("PROCESSING CALIBRATION QUESTIONS\n")
  cat("========================================\n\n")
  
  # Read Wave 1 data
  cat("Reading Wave 1 data...\n")
  df_w1 <- read.csv(wave1_file, stringsAsFactors = FALSE)
  df_w1 <- df_w1[-1, ]  # Drop first row
  
  # Read Wave 2 data  
  cat("Reading Wave 2 data...\n")
  df_w2 <- read.csv(wave2_file, stringsAsFactors = FALSE)
  df_w2 <- df_w2[-1, ]  # Drop first row
  
  # Calibration questions are:
  # - line1_rating_1, line2_rating_1, line3_rating_1
  # - circle1_rating_1, circle2_rating_1, circle3_rating_1
  # - area1_rating_1, area2_rating_1, area3_rating_1
  
  # Wave 1 calibration
  calib_cols <- c("line1_rating_1", "line2_rating_1", "line3_rating_1",
                  "circle1_rating_1", "circle2_rating_1", "circle3_rating_1",
                  "area1_rating_1", "area2_rating_1", "area3_rating_1")
  
  # New naming in the data (without _1 suffix)
  calib_cols_new <- c("line_1", "line_2", "line_3",
                      "circle_1", "circle_2", "circle_3",
                      "region_1", "region_2", "region_3")
  
  w1_calib_data <- list()
  
  for (i in 1:nrow(df_w1)) {
    row_data <- df_w1[i, ]
    response_id <- row_data$ResponseId
    prolific_id <- row_data$PROLIFIC_PID
    
    # Skip if ResponseId is empty or invalid
    if (is.na(response_id) || response_id == "" || 
        grepl("ImportId|_recordId|Response ID", response_id)) {
      next
    }
    
    # Try both old and new naming conventions
    for (q_idx in 1:9) {
      # Try old naming first
      calib_val <- as.numeric(row_data[[calib_cols[q_idx]]])
      
      # If NA, try new naming
      if (is.na(calib_val) && q_idx <= length(calib_cols_new)) {
        calib_val <- as.numeric(row_data[[calib_cols_new[q_idx]]])
      }
      
      # Only include if non-missing
      if (!is.na(calib_val)) {
        question_name <- ifelse(q_idx <= length(calib_cols_new), 
                                calib_cols_new[q_idx], 
                                calib_cols[q_idx])
        
        calib_row <- data.frame(
          ResponseId = response_id,
          PROLIFIC_PID = prolific_id,
          Wave = 1,
          Calibration_Question = question_name,
          Response = calib_val,
          stringsAsFactors = FALSE
        )
        w1_calib_data[[length(w1_calib_data) + 1]] <- calib_row
      }
    }
  }
  
  w1_calib <- bind_rows(w1_calib_data)
  cat("Wave 1 calibration data:", nrow(w1_calib), "rows\n")
  cat("  - Unique respondents:", length(unique(w1_calib$ResponseId)), "\n\n")
  
  # Wave 2 calibration
  w2_calib_data <- list()
  
  for (i in 1:nrow(df_w2)) {
    row_data <- df_w2[i, ]
    response_id <- row_data$ResponseId
    prolific_id <- row_data$PROLIFIC_PID
    
    # Skip if ResponseId is empty or invalid
    if (is.na(response_id) || response_id == "" || 
        grepl("ImportId|_recordId|Response ID", response_id)) {
      next
    }
    
    # Try both old and new naming conventions
    for (q_idx in 1:9) {
      # Try old naming first
      calib_val <- as.numeric(row_data[[calib_cols[q_idx]]])
      
      # If NA, try new naming
      if (is.na(calib_val) && q_idx <= length(calib_cols_new)) {
        calib_val <- as.numeric(row_data[[calib_cols_new[q_idx]]])
      }
      
      # Only include if non-missing
      if (!is.na(calib_val)) {
        question_name <- ifelse(q_idx <= length(calib_cols_new), 
                                calib_cols_new[q_idx], 
                                calib_cols[q_idx])
        
        calib_row <- data.frame(
          ResponseId = response_id,
          PROLIFIC_PID = prolific_id,
          Wave = 2,
          Calibration_Question = question_name,
          Response = calib_val,
          stringsAsFactors = FALSE
        )
        w2_calib_data[[length(w2_calib_data) + 1]] <- calib_row
      }
    }
  }
  
  w2_calib <- bind_rows(w2_calib_data)
  cat("Wave 2 calibration data:", nrow(w2_calib), "rows\n")
  cat("  - Unique respondents:", length(unique(w2_calib$ResponseId)), "\n\n")
  
  # Merge Wave 1 and Wave 2 calibration by PROLIFIC_PID
  cat("Merging calibration data...\n")
  
  merged_calib <- w1_calib %>%
    inner_join(
      w2_calib,
      by = c("PROLIFIC_PID", "Calibration_Question"),
      suffix = c("_W1", "_W2")
    ) %>%
    select(
      ResponseId_W1,
      ResponseId_W2,
      PROLIFIC_PID,
      Calibration_Question,
      Response_W1,
      Response_W2
    ) %>%
    mutate(
      Response_diff = Response_W2 - Response_W1
    )
  
  cat("Merged calibration data:\n")
  cat("  - Total rows:", nrow(merged_calib), "\n")
  cat("  - Unique respondents:", length(unique(merged_calib$PROLIFIC_PID)), "\n")
  cat("  - Unique questions:", length(unique(merged_calib$Calibration_Question)), "\n\n")
  
  return(merged_calib)
}

# ========================================
# MAIN EXECUTION
# ========================================

# Set file paths
wave1_file <- "/Users/liuqianyi/Dropbox/swb/Epinion/raw/SWB wave 1 (Nov 10).csv"
wave2_file <- "/Users/liuqianyi/Dropbox/swb/Epinion/raw/SWB wave 2 (Nov 10).csv"

# Process respondent-aspect level data
aspect_data <- process_respondent_aspect_data(wave1_file, wave2_file)

# Save respondent-aspect data
output_file_1 <- "test_retest_aspect_level.csv"
write.csv(aspect_data, output_file_1, row.names = FALSE)
cat("Saved respondent-aspect data to:", output_file_1, "\n\n")

# Process calibration questions
calibration_data <- process_calibration_data(wave1_file, wave2_file)

# Save calibration data
output_file_2 <- "test_retest_calibration.csv"
write.csv(calibration_data, output_file_2, row.names = FALSE)
cat("Saved calibration data to:", output_file_2, "\n\n")
