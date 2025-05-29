###############################################################################
##   Robust state-wise GAM of COVID-19 deaths  –  *complete working script*
###############################################################################

## ── 0 · packages ────────────────────────────────────────────────────────────
pkgs <- c("tidyverse", "lubridate", "mgcv",
          "png", "grid", "gridExtra", "ggplot2", "xtable")
invisible(lapply(pkgs, \(p){
  if (!requireNamespace(p, quietly = TRUE)) install.packages(p)
  library(p, character.only = TRUE)}))

## ── 1 · user settings ───────────────────────────────────────────────────────
DATA_CSV  <- "India_state_wise_daily.csv"      # <- adapt if needed
OUT_ROOT  <- "R_GAM_COVID_Death_Outputs"
dir.create(OUT_ROOT, showWarnings = FALSE)

K_TIME    <- 20              # knots for s(time_numeric)
SHORT_H   <- 1:7             # 1–7-day horizons
H_LONG    <- 14              # 14-day horizon
META_COLS <- c("Date", "Status", "Date_YMD")

## · helpers ------------------------------------------------------------------
silent_gam <- function(...) {            # run GAM – suppress console noise,
  suppressMessages(suppressWarnings(gam(...)))  # but *return* the model
}
safe_gs <- function(file, plt, ...) {    # wrapper so ggsave never aborts loop
  tryCatch(ggsave(filename = file, plot = plt, ...),
           error = \(e) message("✘ ggsave(): ", e$message))
}

## ── 2 · load master CSV ─────────────────────────────────────────────────────
raw <- read_csv(DATA_CSV, show_col_types = FALSE) |>
  mutate(Date_YMD = ymd(Date_YMD))
state_cols <- setdiff(names(raw), META_COLS)

## ── 3 · loop over every state / UT ──────────────────────────────────────────
for (st in state_cols) {
  
  cat("\n=======================================================\n",
      "▶  STATE :", st, "(Deaths)\n",
      "=======================================================\n", sep = "")
  
  ## folders ------------------------------------------------------------------
  out_dir   <- file.path(OUT_ROOT, paste0("DEATH_", st))
  short_dir <- file.path(out_dir,  "short_term_forecasts")
  dir.create(short_dir, recursive = TRUE, showWarnings = FALSE)
  
  ## 3·1 · preprocess ---------------------------------------------------------
  dat <- raw |>
    filter(Status == "Deceased") |>
    select(Date_YMD, Deaths = !!sym(st)) |>
    mutate(Deaths = as.numeric(Deaths)) |>
    drop_na(Deaths) |>                      # <-- REMOVE missing/non-numeric rows
    arrange(Date_YMD) |>
    mutate(time_numeric = as.numeric(Date_YMD - min(Date_YMD)) + 1,
           day_of_week  = factor(wday(Date_YMD, label = TRUE, week_start = 1)),
           lag1         = lag(Deaths),
           lag7         = lag(Deaths, 7),
           log_lag1     = log(lag1 + 1),
           log_lag7     = log(lag7 + 1)) |>
    drop_na(log_lag1, log_lag7)
  
  ## sanity check -------------------------------------------------------------
  if (nrow(dat) < 50 || all(dat$Deaths == 0)) {
    cat("⚠  skipped – too few (or zero) deaths.\n");  next
  }
  
  ## 3·2 · main GAM fit -------------------------------------------------------
  fmla <- Deaths ~ s(time_numeric, bs = "ps", k = K_TIME) +
    log_lag1 + log_lag7 + day_of_week
  mod  <- tryCatch(silent_gam(fmla, family = nb(), data = dat, method = "REML"),
                   error = \(e) {message("✘ gam(): ", e$message); NULL})
  if (is.null(mod)) next                   # jump to next state on failure
  
  ## 3·3 · diagnostics & static plots ----------------------------------------
  png(file.path(out_dir, "gam_diag.png"), 800, 800)
  try(suppressMessages(gam.check(mod)), silent = TRUE)
  dev.off()
  
  p_fit <- ggplot(dat, aes(Date_YMD)) +
    geom_line(aes(y = Deaths,  colour = "Observed"), alpha = .7) +
    geom_line(aes(y = predict(mod, type = "response"),
                  colour = "Fitted"), linetype = "dashed") +
    scale_colour_manual(values = c(Observed = "#2C77B8",
                                   Fitted   = "#D5402B")) +
    labs(title = paste("Observed vs Fitted deaths –", st),
         y = "Deaths", x = NULL, colour = NULL) +
    theme_minimal() + theme(legend.position = "top")
  safe_gs(file.path(out_dir, "obs_fit.png"), p_fit, width = 8, height = 4)
  
  png(file.path(out_dir, "smooth_time.png"), 800, 600)
  try(plot(mod, select = 1, residuals = TRUE,
           main = paste("s(time_numeric) –", st)), silent = TRUE)
  dev.off()
  
  write_csv(as.data.frame(summary(mod)$p.table),
            file.path(out_dir, "coeff_param.csv"))
  write_csv(as.data.frame(summary(mod)$s.table),
            file.path(out_dir, "coeff_smooth.csv"))
  
  ## 3·4 · 14-day rolling-origin CV ------------------------------------------
  n        <- nrow(dat)
  initial  <- floor(0.70 * n)
  origins  <- seq(initial, n - H_LONG, by = H_LONG)
  
  cv_long <- purrr::imap_dfr(origins, \(o, idx){
    tr <- dat[1:o, ];  te <- dat[(o+1):(o+H_LONG), ]
    m  <- tryCatch(silent_gam(fmla, family = nb(), data = tr, method = "REML"),
                   error = \(e) NULL)
    if (is.null(m)) return(tibble())                  # skip fold on error
    pr <- predict(m, newdata = te, type = "response")
    tibble(Fold = idx, TrainEnd = o, TestStart = o+1, TestEnd = o+H_LONG,
           MAE  = mean(abs(te$Deaths - pr)),
           RMSE = sqrt(mean((te$Deaths - pr)^2)))
  })
  write_csv(cv_long, file.path(out_dir, "cv14_long.csv"))
  
  ## 3·5 · short-term CV (h = 1…7) -------------------------------------------
  short_summary <- list()
  
  for (h in SHORT_H) {
    cat("   • horizon h =", h, "\n")
    origins_h <- origins[origins <= n - h]
    
    cv_h <- purrr::imap_dfr(origins_h, \(o, fid){
      tr <- dat[1:o, ];  te <- dat[(o+1):(o+h), ]
      m  <- tryCatch(silent_gam(fmla, family = nb(), data = tr, method = "REML"),
                     error = \(e) NULL)
      if (is.null(m)) return(tibble())
      pr <- predict(m, newdata = te, type = "response")
      
      ## per-fold plot --------------------------------------------------------
      p_fold <- ggplot(te, aes(Date_YMD)) +
        geom_line(aes(y = Deaths, colour = "Actual")) +
        geom_line(aes(y = pr,     colour = "Forecast"),
                  linetype = "dashed") +
        scale_colour_manual(values = c(Actual   = "#2C77B8",
                                       Forecast = "#D5402B")) +
        labs(title = paste0(st, " | h = ", h, "  Fold ", fid),
             y = "Deaths", x = NULL, colour = NULL) +
        theme_minimal(base_size = 9) +
        theme(legend.position = "bottom")
      safe_gs(file.path(short_dir,
                        sprintf("%s_h%d_fold%02d.png", st, h, fid)),
              p_fold, width = 4, height = 3, dpi = 160)
      
      tibble(Fold = fid, TrainEnd = o, TestStart = o+1, TestEnd = o+h,
             MAE = mean(abs(te$Deaths - pr)),
             RMSE = sqrt(mean((te$Deaths - pr)^2)))
    })
    
    write_csv(cv_h, file.path(short_dir, sprintf("%s_cv_h%d.csv", st, h)))
    
    if (nrow(cv_h))
      short_summary[[h]] <- cv_h |>
      summarise(Horizon = h,
                Min_MAE  = min(MAE),
                Min_RMSE = min(RMSE))
    
    ## grid mosaic -----------------------------------------------------------
    pngs <- list.files(short_dir,
                       pattern  = sprintf("^%s_h%d_fold.*png$", st, h),
                       full.names = TRUE)
    if (length(pngs)) {
      gbs <- lapply(pngs, \(f) rasterGrob(readPNG(f)))
      nc  <- 5; nr <- ceiling(length(gbs)/nc)
      grd <- arrangeGrob(grobs = gbs, ncol = nc, nrow = nr,
                         top = textGrob(paste(st, "| h =", h, "days"),
                                        gp = gpar(fontsize = 14,
                                                  fontface = "bold")))
      safe_gs(file.path(short_dir,
                        sprintf("%s_grid_h%d.png", st, h)),
              grd, width = nc*2, height = nr*1.8, dpi = 160)
    }
  } # end h loop
  
  if (length(short_summary))
    write_csv(bind_rows(short_summary),
              file.path(short_dir, sprintf("%s_short_summary.csv", st)))
  
  ## 3·6 · LaTeX tables -------------------------------------------------------
  latex_dir <- file.path(short_dir, "latex_tables")
  dir.create(latex_dir, showWarnings = FALSE)
  
  for (h in SHORT_H) {
    csv <- file.path(short_dir, sprintf("%s_cv_h%d.csv", st, h))
    if (!file.exists(csv)) next
    tab <- read_csv(csv, show_col_types = FALSE)
    if (!nrow(tab)) next
    xt  <- xtable(tab,
                  caption = sprintf("Cross-validation results — %s, %d-day horizon.",
                                    st, h),
                  label   = sprintf("tab:%s_h%d", st, h),
                  digits  = c(0,0,0,0,0,3,3))
    print(xt, include.rownames = FALSE, tabular.environment = "longtable",
          floating = TRUE, hline.after = c(-1, 0, nrow(tab)),
          file = file.path(latex_dir,
                           sprintf("%s_cv_h%d.tex", st, h)))
  }
  
  cat("✔  Finished state ", st,
      " — outputs in ", out_dir, "\n", sep = "")
}
