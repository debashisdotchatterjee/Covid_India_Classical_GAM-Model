# R Code for GAM Analysis of COVID-19 Data (Total India - TT)

# --- 0. Install and Load Necessary Libraries ---
# Ensure these packages are installed: install.packages(c("tidyverse", "mgcv", "lubridate", "here", "gridExtra", "knitr"))

library(tidyverse) # For data manipulation (dplyr, tidyr, ggplot2) and reading CSV
library(mgcv)      # For Generalized Additive Models (GAMs)
library(lubridate) # For date manipulation
library(here)      # For path management (optional, makes path handling easier)
library(gridExtra) # For arranging multiple plots (if needed)
library(knitr)     # For creating tables (kable)

# --- 1. Setup Output Directory ---
output_dir <- "R_GAM_COVID_Analysis_Outputs_TT" # Changed directory name for TT
if (!dir.exists(output_dir)) {
  dir.create(output_dir, recursive = TRUE)
  cat("Created output directory:", file.path(getwd(), output_dir), "\n")
} else {
  cat("Output directory already exists:", file.path(getwd(), output_dir), "\n")
}

# --- 2. Load Data ---
# Assuming the data file "India_state_wise_daily.csv" is in the current working directory or a known path
file_path <- "India_state_wise_daily.csv" # MODIFY THIS PATH IF THE FILE IS ELSEWHERE

if (!file.exists(file_path)) {
  stop("Data file not found at: ", file_path, ". Please check the path.")
}

cat("Loading data from:", file_path, "\n")
raw_data <- read_csv(file_path, show_col_types = FALSE)
cat("Data loaded successfully. Dimensions:", dim(raw_data)[1], "rows,", dim(raw_data)[2], "columns.\n")
# print(head(raw_data)) # Print head of raw data if needed

# --- 3. Preprocess Data for Total India (TT) ---
cat("\n--- Preprocessing Data for Total India (TT) ---\n")

# Filter for "Confirmed" status and select Date_YMD and TT columns
# Rename TT to Cases for consistency with the model formula
total_india_data <- raw_data %>%
  mutate(Date_YMD = ymd(Date_YMD)) %>%
  filter(Status == "Confirmed") %>% # Focus on Confirmed cases
  select(Date_YMD, Cases = TT) %>%   # Select TT column and rename to Cases
  arrange(Date_YMD)

if (nrow(total_india_data) == 0) {
  stop("No data found for 'Confirmed' status and 'TT' column. Please check the input CSV and filters.")
}

cat("Data for Total India (TT) - Confirmed cases extracted. Dimensions:", dim(total_india_data)[1], "rows,", dim(total_india_data)[2], "columns.\n")
print(head(total_india_data))

# Add time-related variables and lagged variables
total_india_data <- total_india_data %>%
  mutate(
    time_numeric = as.numeric(Date_YMD - min(Date_YMD)) + 1, # Days since start
    day_of_week = factor(wday(Date_YMD, label = TRUE, week_start = 1)), # Monday as start
    lagged_cases_1 = lag(Cases, 1),
    lagged_cases_7 = lag(Cases, 7),
    # Add a small constant for log transformation to avoid log(0)
    log_lagged_cases_1 = log(lagged_cases_1 + 1),
    log_lagged_cases_7 = log(lagged_cases_7 + 1)
  )

# Remove rows with NA in lagged predictors (important for model fitting)
data_model_ready_tt <- total_india_data %>%
  filter(!is.na(log_lagged_cases_1) & !is.na(log_lagged_cases_7))

cat("Data prepared for modeling Total India (TT). Removed rows with NA in lagged variables.\n")
print(head(data_model_ready_tt))

# --- 4. Fit GAM and Generate Outputs for Total India (TT) ---
current_analysis_label <- "TT" # Label for file names and titles

cat("\n\n=========================================================\n")
cat("Processing Total India (", current_analysis_label, ")\n")
cat("=========================================================\n\n")

if (nrow(data_model_ready_tt) < 50) { # Basic check for sufficient data
  stop("Insufficient data points for Total India (", current_analysis_label, ") after preprocessing (", nrow(data_model_ready_tt), "rows).\nAnalysis halted.")
}

# Define k for splines.
k_val_time = 30 # Number of knots for the time spline for TT.
# For a national aggregate over a longer period, more knots might be appropriate.
# Adjust as needed, e.g., 10-15 knots per year of data.

cat("--- Fitting GAM for", current_analysis_label, "---\n")
cat("Number of observations:", nrow(data_model_ready_tt), "\n")
cat("Using k =", k_val_time, "for time spline.\n")

# GAM formula
gam_formula_tt <- Cases ~ s(time_numeric, bs = "ps", k = k_val_time) +
  log_lagged_cases_1 +
  log_lagged_cases_7 +
  day_of_week

# Fit the GAM model
gam_model_tt <- NULL
tryCatch({
  gam_model_tt <- gam(gam_formula_tt,
                      family = nb(), # Negative Binomial family
                      data = data_model_ready_tt,
                      method = "REML") # REML for smoothness selection
}, error = function(e) {
  cat("Error fitting GAM for", current_analysis_label, ":", e$message, "\n")
  stop("GAM fitting failed for ", current_analysis_label, ". Analysis halted.")
})

if (is.null(gam_model_tt)) {
  stop("GAM model is NULL after fitting attempt for ", current_analysis_label, ". Analysis halted.")
}

# --- 4.1. Model Summary ---
cat("\n--- GAM Summary for", current_analysis_label, "---\n")
model_summary_tt <- summary(gam_model_tt)
print(model_summary_tt)

summary_file_tt <- file.path(output_dir, paste0("gam_summary_", current_analysis_label, ".txt"))
sink(summary_file_tt)
print(model_summary_tt)
cat("\n\nCoefficients (parametric part):\n")
print(model_summary_tt$p.table)
cat("\n\nSmooth terms (sbtable):\n")
print(model_summary_tt$s.table)
sink()
cat("Model summary saved to:", summary_file_tt, "\n")

# --- 4.2. Coefficients Table ---
cat("\n--- Coefficients Table for", current_analysis_label, "---\n")
coeffs_parametric_tt <- as.data.frame(model_summary_tt$p.table)
colnames(coeffs_parametric_tt) <- c("Estimate", "Std.Error", "t.value", "Pr(>|t|)")
cat("Parametric Coefficients:\n")
print(kable(coeffs_parametric_tt, format = "pipe", digits = 4))

coeffs_smooth_tt <- as.data.frame(model_summary_tt$s.table)
colnames(coeffs_smooth_tt) <- c("EDF", "Ref.df", "F.value", "p.value")
cat("\nSmooth Terms:\n")
print(kable(coeffs_smooth_tt, format = "pipe", digits = 4))

coeffs_file_parametric_tt <- file.path(output_dir, paste0("gam_coeffs_parametric_", current_analysis_label, ".csv"))
write.csv(coeffs_parametric_tt, coeffs_file_parametric_tt)
cat("Parametric coefficients saved to:", coeffs_file_parametric_tt, "\n")

coeffs_file_smooth_tt <- file.path(output_dir, paste0("gam_coeffs_smooth_", current_analysis_label, ".csv"))
write.csv(coeffs_smooth_tt, coeffs_file_smooth_tt)
cat("Smooth term details saved to:", coeffs_file_smooth_tt, "\n")

# --- 4.3. Diagnostic Plots (gam.check) ---
cat("\n--- Generating gam.check() plots for", current_analysis_label, "---\n")
diag_plot_file_tt <- file.path(output_dir, paste0("gam_diag_plots_", current_analysis_label, ".png"))
png(diag_plot_file_tt, width = 800, height = 800)
par(mfrow = c(2, 2))
gam.check(gam_model_tt)
dev.off()
par(mfrow = c(1, 1)) # Reset par
cat("Diagnostic plots saved to:", diag_plot_file_tt, "\n")
# To show in console (if interactive):
# par(mfrow = c(2, 2)); gam.check(gam_model_tt); par(mfrow = c(1, 1))

# --- 4.4. Fitted vs. Observed Plot ---
cat("\n--- Generating Fitted vs. Observed Plot for", current_analysis_label, "---\n")
data_model_ready_tt_with_fits <- data_model_ready_tt %>%
  mutate(Fitted_Cases = predict(gam_model_tt, type = "response"))

fitted_plot_tt <- ggplot(data_model_ready_tt_with_fits, aes(x = Date_YMD)) +
  geom_line(aes(y = Cases, color = "Observed"), alpha = 0.7) +
  geom_line(aes(y = Fitted_Cases, color = "Fitted"), linewidth = 0.8, linetype = "dashed") +
  labs(title = paste("Observed vs. Fitted COVID-19 Cases:", current_analysis_label),
       x = "Date", y = "Number of Confirmed Cases",
       color = "Legend") +
  scale_color_manual(values = c("Observed" = "steelblue", "Fitted" = "darkred")) +
  theme_minimal() +
  theme(legend.position = "top")

print(fitted_plot_tt)
fitted_plot_file_tt <- file.path(output_dir, paste0("gam_fitted_vs_observed_", current_analysis_label, ".png"))
ggsave(fitted_plot_file_tt, plot = fitted_plot_tt, width = 10, height = 6)
cat("Fitted vs. Observed plot saved to:", fitted_plot_file_tt, "\n")

# --- 4.5. Smooth Term Plot (Time Trend) ---
cat("\n--- Generating Smooth Term Plot (Time Trend) for", current_analysis_label, "---\n")
smooth_plot_file_tt <- file.path(output_dir, paste0("gam_smooth_time_trend_", current_analysis_label, ".png"))
png(smooth_plot_file_tt, width = 800, height = 600)
plot(gam_model_tt, select = 1, residuals = TRUE, pch = 1, cex = 0.8, shade = TRUE,
     main = paste("Non-linear Time Trend (s(time_numeric)) for", current_analysis_label),
     xlab = "Time (Numeric Index)", ylab = "Effect on Log(Cases)")
dev.off()
cat("Smooth term plot saved to:", smooth_plot_file_tt, "\n")
# To show in console:
# plot(gam_model_tt, select = 1, residuals = TRUE, pch = 1, cex = 0.8, shade = TRUE,
#      main = paste("Non-linear Time Trend (s(time_numeric)) for", current_analysis_label),
#      xlab = "Time (Numeric Index)", ylab = "Effect on Log(Cases)")

# --- 4.6. Plotting Day of the Week Effect ---
cat("\n--- Generating Day of the Week Effect Plot for", current_analysis_label, "---\n")
day_effect_plot_file_tt <- file.path(output_dir, paste0("gam_day_of_week_effect_", current_analysis_label, ".png"))
png(day_effect_plot_file_tt, width = 800, height = 600)
termplot(gam_model_tt, terms = "day_of_week", se = TRUE, rug = FALSE,
         main = paste("Day of the Week Effect for", current_analysis_label),
         ylab = "Effect on Log(Cases)")
dev.off()
cat("Day of week effect plot saved to:", day_effect_plot_file_tt, "\n")

cat("\nFinished processing for Total India:", current_analysis_label, "\n")

# --- 5. Illustrative Time-Series Cross-Validation for Total India (TT) ---
cat("\n\n=========================================================\n")
cat("Illustrative Time-Series Cross-Validation for", current_analysis_label, "\n")
cat("=========================================================\n\n")

if (nrow(data_model_ready_tt) > 100) { # Need enough data for CV
  
  n_obs_tt <- nrow(data_model_ready_tt)
  initial_train_window_prop_tt <- 0.7
  initial_train_samples_tt <- floor(n_obs_tt * initial_train_window_prop_tt)
  forecast_horizon_h_tt <- 14 # Predict 14 days ahead
  step_m_tt <- 14 # Roll origin by 14 days
  
  cv_results_tt <- list()
  forecast_origins_tt <- seq(initial_train_samples_tt, n_obs_tt - forecast_horizon_h_tt, by = step_m_tt)
  
  cat("Number of observations for CV:", n_obs_tt, "\n")
  cat("Initial training samples:", initial_train_samples_tt, "\n")
  cat("Forecast horizon (h):", forecast_horizon_h_tt, "\n")
  cat("Rolling step (m):", step_m_tt, "\n")
  cat("Number of CV folds to run:", length(forecast_origins_tt), "\n")
  
  fold_count_tt <- 0
  for (train_end_idx_tt in forecast_origins_tt) {
    fold_count_tt <- fold_count_tt + 1
    cat("CV Fold:", fold_count_tt, "- Training up to index", train_end_idx_tt, "\n")
    
    train_data_tt <- data_model_ready_tt[1:train_end_idx_tt, ]
    test_data_tt <- data_model_ready_tt[(train_end_idx_tt + 1):(min(train_end_idx_tt + forecast_horizon_h_tt, n_obs_tt)), ]
    
    if(nrow(test_data_tt) == 0) next
    
    k_val_time_cv_tt = k_val_time # Use the same k as main model or adjust if needed for smaller training sets
    
    gam_formula_cv_tt <- Cases ~ s(time_numeric, bs = "ps", k = k_val_time_cv_tt) +
      log_lagged_cases_1 +
      log_lagged_cases_7 +
      day_of_week
    
    cv_gam_model_tt <- NULL
    tryCatch({
      cv_gam_model_tt <- gam(gam_formula_cv_tt,
                             family = nb(),
                             data = train_data_tt,
                             method = "REML")
    }, error = function(e) {
      cat("Error in CV GAM fit for fold", fold_count_tt, ":", e$message, "\n")
    })
    
    if (!is.null(cv_gam_model_tt)) {
      predictions_tt <- predict(cv_gam_model_tt, newdata = test_data_tt, type = "response")
      
      valid_preds_tt <- !is.na(predictions_tt)
      actuals_tt <- test_data_tt$Cases[valid_preds_tt]
      preds_tt <- predictions_tt[valid_preds_tt]
      
      if(length(actuals_tt) > 0 && length(preds_tt) > 0) {
        mae_tt <- mean(abs(actuals_tt - preds_tt))
        rmse_tt <- sqrt(mean((actuals_tt - preds_tt)^2))
        
        cv_results_tt[[fold_count_tt]] <- data.frame(
          Fold = fold_count_tt,
          TrainEndIndex = train_end_idx_tt,
          TestStartIndex = train_end_idx_tt + 1,
          TestEndIndex = min(train_end_idx_tt + forecast_horizon_h_tt, n_obs_tt),
          MAE = mae_tt,
          RMSE = rmse_tt
        )
        cat("Fold", fold_count_tt, "- MAE:", round(mae_tt,2), "RMSE:", round(rmse_tt,2), "\n")
      } else {
        cat("Fold", fold_count_tt, "- No valid predictions or actuals to calculate error.\n")
      }
    } else {
      cat("Fold", fold_count_tt, "- GAM fitting failed, skipping error calculation.\n")
    }
  }
  
  if (length(cv_results_tt) > 0) {
    cv_results_df_tt <- bind_rows(cv_results_tt)
    cat("\n--- CV Results Summary for", current_analysis_label, "---\n")
    print(kable(cv_results_df_tt, format="pipe", digits=3))
    
    cv_summary_file_tt <- file.path(output_dir, paste0("gam_cv_summary_", current_analysis_label, ".csv"))
    write.csv(cv_results_df_tt, cv_summary_file_tt, row.names = FALSE)
    cat("CV summary saved to:", cv_summary_file_tt, "\n")
    
    cat("\nAverage MAE across folds:", mean(cv_results_df_tt$MAE, na.rm=TRUE), "\n")
    cat("Average RMSE across folds:", mean(cv_results_df_tt$RMSE, na.rm=TRUE), "\n")
  } else {
    cat("No results from CV to summarize for", current_analysis_label, "\n")
  }
  
} else {
  cat("Skipping CV for", current_analysis_label, "due to insufficient data.\n")
}

cat("\n\n--- R Script Execution Finished ---\n")
cat("All outputs saved in directory:", file.path(getwd(), output_dir), "\n")

##########################################

# --- 6. Short-Term Forecast Performance (1 to 4 Days Ahead) ---

short_term_results <- list()
short_term_plot_dir <- file.path(output_dir, "short_term_forecasts")
dir.create(short_term_plot_dir, showWarnings = FALSE)

for (h in 1:4) {
  cat("\n\n--- Running Short-Term CV for h =", h, "day(s) ahead ---\n")
  
  cv_results_short <- list()
  fold_counter <- 0
  for (train_end_idx in forecast_origins_tt) {
    fold_counter <- fold_counter + 1
    train_data <- data_model_ready_tt[1:train_end_idx, ]
    test_data <- data_model_ready_tt[(train_end_idx + 1):(min(train_end_idx + h, n_obs_tt)), ]
    
    if (nrow(test_data) == 0) next
    
    gam_model_cv <- tryCatch({
      gam(Cases ~ s(time_numeric, bs = "ps", k = k_val_time) +
            log_lagged_cases_1 + log_lagged_cases_7 + day_of_week,
          family = nb(), data = train_data, method = "REML")
    }, error = function(e) NULL)
    
    if (!is.null(gam_model_cv)) {
      preds <- predict(gam_model_cv, newdata = test_data, type = "response")
      valid <- !is.na(preds)
      actual <- test_data$Cases[valid]
      predicted <- preds[valid]
      
      if (length(actual) > 0 && length(predicted) > 0) {
        mae <- mean(abs(actual - predicted))
        rmse <- sqrt(mean((actual - predicted)^2))
        
        cv_results_short[[fold_counter]] <- data.frame(
          Fold = fold_counter,
          Horizon = h,
          MAE = mae,
          RMSE = rmse
        )
        
        # Save forecast plot
        # Save forecast plot
        df_plot <- tibble(
          Date      = test_data$Date_YMD[valid],
          Actual    = actual,
          Predicted = predicted
        )
        
        p <- ggplot(df_plot, aes(x = Date, group = 1)) +        # ← add group = 1
          geom_line(aes(y = Actual,    color = "Actual"),    size = 1) +
          geom_line(aes(y = Predicted, color = "Predicted"),
                    linetype = "dashed", size = 1) +
          scale_color_manual(values = c("Actual" = "steelblue",
                                        "Predicted" = "darkred")) +
          labs(title = paste("h =", h, "day ahead – Fold", fold_counter),
               y = "Confirmed cases", color = NULL) +
          theme_minimal()
        
        ggsave(
          filename = file.path(short_term_plot_dir,
                               sprintf("forecast_h%d_fold%02d.png", h, fold_counter)),
          plot    = p, width = 8, height = 5
        )
      }
    }
  }
  
  # Bind and summarize short-term results
  df_short <- bind_rows(cv_results_short)
  write.csv(df_short, file.path(short_term_plot_dir, paste0("cv_short_term_h", h, ".csv")), row.names = FALSE)
  
  if (nrow(df_short) > 0) {
    cat("\nSummary for h =", h, "\n")
    print(summary(df_short[, c("MAE", "RMSE")]))
    
    short_term_results[[h]] <- data.frame(
      Horizon = h,
      Avg_MAE = mean(df_short$MAE, na.rm = TRUE),
      Avg_RMSE = mean(df_short$RMSE, na.rm = TRUE)
    )
  }
}

# Final summary table
df_short_summary <- bind_rows(short_term_results)
write.csv(df_short_summary, file.path(short_term_plot_dir, "cv_short_term_summary.csv"), row.names = FALSE)
cat("Short-term forecast summary saved to:", file.path(short_term_plot_dir, "cv_short_term_summary.csv"), "\n")
