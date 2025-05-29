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
k_val_time = 20 # Number of knots for the time spline for TT.
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


####################################################################



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

#######################################################################################


# --- 6. Create Grid Plots for Each Short-Term Forecast Horizon via gridExtra ---

library(png)      # to read PNG files
library(grid)     # to convert them into grobs
library(gridExtra)

cat("\n--- Assembling grid plots for short-term forecasts (h = 1 to 4) ---\n")

plot_dir <- file.path(output_dir, "short_term_forecasts")
if (!dir.exists(plot_dir)) stop("No such plot directory: ", plot_dir)

for (h in 1:7) {
  cat(" • Horizon h =", h, "\n")
  # find all PNGs for this horizon
  files <- list.files(plot_dir,
                      pattern = sprintf("forecast_h%d_fold.*\\.png$", h),
                      full.names = TRUE)
  if (length(files) == 0) {
    warning("  no files found for h =", h); next
  }
  
  # read each into a rasterGrob
  grobs <- lapply(files, function(f) {
    img <- readPNG(f)
    rasterGrob(img, interpolate = TRUE)
  })
  
  # assemble into a 5×ceiling(n/5) grid
  ncol <- 5
  nrow <- ceiling(length(grobs) / ncol)
  grid_plot <- arrangeGrob(grobs = grobs, ncol = ncol, nrow = nrow)
  
  # save to PNG
  out_file <- file.path(plot_dir, sprintf("grid_h%d.png", h))
  ggsave(out_file, grid_plot,
         width  = ncol * 3,   # 3 inches per panel
         height = nrow * 2.5, # 2.5 inches per panel
         dpi    = 150)
  cat("   → saved grid to", out_file, "\n")
}


# Load required packages
library(png)
library(grid)
library(gridExtra)
library(ggplot2)

cat("\n--- Creating Enhanced Grid Plots for Short-Term Forecasts (h = 1 to 7) ---\n")

# Define directory with individual forecast plots
plot_dir <- file.path(output_dir, "short_term_forecasts")
if (!dir.exists(plot_dir)) stop("No such plot directory: ", plot_dir)

# Define layout settings
ncol <- 4  # number of columns in grid
scale_factor <- 3  # adjust for better quality

# Loop over each forecast horizon
for (h in 1:7) {
  cat(" → Horizon h =", h, "\n")
  
  # Find forecast PNGs for current horizon
  files <- list.files(plot_dir,
                      pattern = sprintf("forecast_h%d_fold.*\\.png$", h),
                      full.names = TRUE)
  if (length(files) == 0) {
    warning("  No forecast plots found for h =", h)
    next
  }
  
  # Convert each PNG to a raster grob
  grobs <- lapply(files, function(f) {
    tryCatch({
      img <- readPNG(f)
      rasterGrob(img, interpolate = FALSE)
    }, error = function(e) {
      cat("   Error reading file:", f, "\n")
      NULL
    })
  })
  grobs <- Filter(Negate(is.null), grobs)  # remove failed loads
  
  # Set up grid dimensions
  nrow <- ceiling(length(grobs) / ncol)
  
  # Arrange the grid
  grid_plot <- arrangeGrob(grobs = grobs, ncol = ncol, nrow = nrow,
                           top = textGrob(paste("Short-Term Forecasts (h =", h, ")"),
                                          gp = gpar(fontsize = 18, fontface = "bold")))
  
  # Save to file
  grid_file <- file.path(plot_dir, sprintf("grid_h%d.png", h))
  ggsave(grid_file, grid_plot,
         width = ncol * scale_factor,
         height = nrow * scale_factor,
         dpi = 200)
  
  cat("    ✔ Saved:", grid_file, "\n")
}


##############################################################

## --- 6. Generate 1–7 day forecasts and plots ---
for (h in 1:7) {
  cat("\n\n--- Running Short-Term CV for h =", h, "days ahead ---\n")
  cv_results_short <- list()
  fold_counter <- 0
  for (train_end_idx in forecast_origins_tt) {
    fold_counter <- fold_counter + 1
    train_data <- data_model_ready_tt[1:train_end_idx, ]
    test_data  <- data_model_ready_tt[(train_end_idx + 1):(min(train_end_idx + h, n_obs_tt)), ]
    if (!nrow(test_data)) next
    #… fit the GAM, predict, calc MAE/RMSE …
    # save individual plots to:
     #   file.path(short_term_plot_dir, sprintf("forecast_h%d_fold%02d.png", h, fold_counter))
  }
   #write.csv(cv_results_short, file.path(short_term_plot_dir, paste0("cv_short_term", h, ".csv")))
}

if (is.null(gam_model_tt)) {
  stop("GAM model is NULL after fitting attempt for ", current_analysis_label, ". Analysis halted.")
}


###################

## --- Assembling grid plots for short-term forecasts (h = 1 to 7) ---
      for (h in 1:7) {
        cat(" • Horizon h =", h, "\n")
        files <- list.files(plot_dir,
                            pattern = sprintf("forecast_h%d_fold.*\\.png$", h),
                            full.names = TRUE)
        if (!length(files)) {
          warning("  no files found for h =", h)
          next
        }
        grobs <- lapply(files, function(f) {
          img <- png::readPNG(f)
          grid::rasterGrob(img, interpolate = FALSE)
        })
        ncol <- 5
        nrow <- ceiling(length(grobs) / ncol)
        grid_plot <- gridExtra::arrangeGrob(grobs = grobs, ncol = ncol, nrow = nrow,
                                            top = grid::textGrob(
                                              paste("Short-Term Forecasts (h =", h, ")"),
                                              gp = grid::gpar(fontsize = 16, fontface = "bold")))
        ggsave(file.path(plot_dir, sprintf("grid_h%d.png", h)),
                          grid_plot,
                          width = ncol * 2.5, height = nrow * 2.5, dpi = 150)
        cat("   → saved grid to", file.path(plot_dir, sprintf("grid_h%d.png", h)), "\n")
      }
    
#####################################################

################################################################################
## 6 · Short-term (h = 1–7) rolling-origin CV -- plots, csv, mosaics
################################################################################

library(ggplot2)
library(png)        # readPNG()
library(grid)       # rasterGrob()
library(gridExtra)  # arrangeGrob()

short_term_plot_dir <- file.path(output_dir, "short_term_forecasts")
dir.create(short_term_plot_dir, showWarnings = FALSE, recursive = TRUE)

short_term_results <- list()   # to collect horizon-wise summaries

## ── run CV for each horizon ────────────────────────────────────────────────
for (h in 1:7) {
  
  cat("\n\n──────────────────────────────────────────────────────────\n")
  cat("Short-term CV  |  horizon h =", h, "day(s)\n")
  cat("──────────────────────────────────────────────────────────\n")
  
  cv_res  <- list()   # fold-level results for this horizon
  fold_id <- 0
  
  ## rolling origins – reuse the ones you built earlier,
  ## but make sure there are at least ‘h’ observations left
  origins_h <- forecast_origins_tt[forecast_origins_tt <= n_obs_tt - h]
  
  for (train_end in origins_h) {
    
    fold_id <- fold_id + 1
    
    train_dat <- data_model_ready_tt[ 1:train_end , ]
    test_dat  <- data_model_ready_tt[ (train_end + 1):(train_end + h) , ]
    
    ## --- fit ---------------------------------------------------------------
    m_cv <- tryCatch(
      gam(gam_formula_tt, family = nb(),
          data = train_dat, method = "REML"),
      error = function(e) NULL)
    
    if (is.null(m_cv)) next   # skip fold if it failed
    
    ## --- predict & score ---------------------------------------------------
    pred <- as.numeric(predict(m_cv, newdata = test_dat, type = "response"))
    act  <- test_dat$Cases
    mae  <- mean(abs(act - pred))
    rmse <- sqrt(mean((act - pred)^2))
    
    cv_res[[fold_id]] <- data.frame(
      Fold = fold_id,
      TrainEnd = train_end,
      TestStart = train_end + 1,
      TestEnd = train_end + h,
      MAE  = mae,
      RMSE = rmse
    )
    
    ## --- save per-fold plot -------------------------------------------------
    p <- ggplot(test_dat, aes(Date_YMD)) +
      geom_line(aes(y = act, colour = "Actual"), linewidth = 0.9) +
      geom_line(aes(y = pred, colour = "Forecast"), linewidth = 0.9,
                linetype = "dashed") +
      scale_colour_manual(values = c(Actual = "#2C77B8", Forecast = "#D5402B")) +
      labs(title = paste0("h = ", h,
                          " - Fold ", fold_id,
                          "  (", format(min(test_dat$Date_YMD)), " → ",
                          format(max(test_dat$Date_YMD)), ")"),
           y = "Confirmed cases", x = NULL, colour = NULL) +
      theme_minimal(base_size = 10) +
      theme(legend.position = "bottom",
            plot.title = element_text(size = 10, face = "bold"))
    
    ggsave(file.path(short_term_plot_dir,
                     sprintf("forecast_h%d_fold%02d.png", h, fold_id)),
           p, width = 4.5, height = 3, dpi = 150)
  } # end folds
  
  ## ---------- write fold-level CSV & horizon summary -----------------------
  df_h <- bind_rows(cv_res)
  write.csv(df_h,
            file.path(short_term_plot_dir,
                      sprintf("cv_short_term_h%d.csv", h)),
            row.names = FALSE)
  
  short_term_results[[h]] <- df_h %>%
    summarise(Horizon = h,
              Min_MAE  = min(MAE , na.rm = TRUE),
              Min_RMSE = min(RMSE, na.rm = TRUE))
  
  ## ---------- build mosaic of fold plots ----------------------------------
  png_files <- list.files(short_term_plot_dir,
                          pattern = sprintf("^forecast_h%d_fold.*png$", h),
                          full.names = TRUE)
  
  if (length(png_files)) {
    grobs <- lapply(png_files, function(f) {
      rasterGrob(readPNG(f), interpolate = FALSE)
    })
    ncol <- 5
    nrow <- ceiling(length(grobs)/ncol)
    ggrid <- arrangeGrob(grobs = grobs, ncol = ncol, nrow = nrow,
                         top = textGrob(
                           paste("Short-term forecasts | horizon =", h, "day(s)"),
                           gp = gpar(fontsize = 16, fontface = "bold")))
    ggsave(file.path(short_term_plot_dir,
                     sprintf("grid_h%d.png", h)),
           ggrid,
           width = ncol*2.2, height = nrow*2.0, dpi = 160)
  }
} # end horizons

## ---------- overall summary table -----------------------------------------
summary_df <- bind_rows(short_term_results)
write.csv(summary_df,
          file.path(short_term_plot_dir, "cv_short_term_summary.csv"),
          row.names = FALSE)

cat("\n✔  Short-term forecasting (h = 1…7) finished.",
    "\n   Results & plots in:", short_term_plot_dir, "\n")


#########################################


# ----------------------------------------------------------------------
# 7 · Generate LaTeX Table for Each h (1–7) Based on CV Results
# ----------------------------------------------------------------------

library(xtable)
latex_table_dir <- file.path(short_term_plot_dir, "latex_tables")
dir.create(latex_table_dir, showWarnings = FALSE, recursive = TRUE)

for (h in 1:7) {
  cv_file <- file.path(short_term_plot_dir, sprintf("cv_short_term_h%d.csv", h))
  if (file.exists(cv_file)) {
    cv_data <- read.csv(cv_file)
    
    if (nrow(cv_data) > 0) {
      # Create xtable
      caption_text <- sprintf("Time-series cross-validation results for the GAM model (%d-day forecast horizon) for Total India (TT) confirmed COVID-19 cases.", h)
      label_text <- sprintf("tab:cv_results_tt_h%d", h)
      
      cv_xtable <- xtable(cv_data,
                          caption = caption_text,
                          label = label_text,
                          digits = c(0, 0, 0, 0, 0, 3, 3))
      
      # Save to .tex file
      tex_file <- file.path(latex_table_dir, sprintf("cv_results_h%d.tex", h))
      print(cv_xtable,
            include.rownames = FALSE,
            tabular.environment = "longtable",
            floating = TRUE,
            include.colnames = TRUE,
            hline.after = c(-1, 0, nrow(cv_data)),
            file = tex_file)
      
      cat(sprintf("✔ LaTeX table for h = %d saved to: %s\n", h, tex_file))
    } else {
      cat(sprintf("⚠ No data found in CV result file for h = %d\n", h))
    }
  } else {
    cat(sprintf("⚠ CV result file missing for h = %d: %s\n", h, cv_file))
  }
}


