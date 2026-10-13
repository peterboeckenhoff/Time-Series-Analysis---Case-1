## ============================================================
## Case 1: Forecasting Citi Bike Bicycle Trips
## ARIMA and SARIMA (s = 168, weekly cycle) for area 6 demand
## ============================================================

## ---- 0. Packages -------------------------------------------
## install.packages("forecast")   # run once if needed
library(forecast)  

## ---- 0b. Working directory -----------------------------------
## One working directory for everything (data, caches, ETS csv files, outputs)
setwd("C:/Users/janne/OneDrive/Documents/Time Series Analytics/Case")

## All outputs that are handed in are saved in this one folder
results_dir <- "results"
dir.create(results_dir, showWarnings = FALSE, recursive = TRUE)

## ---- 1. Settings ---------------------------------------------
s <- 168   # seasonal period: weekly cycle in hourly data (168 = 7 * 24)
h <- 24    # forecast horizon in hours; also the step between forecast origins

## ---- 2. Load data and split ----------------------------------
load("citibike.RData")

data_citibike$datetime <- as.POSIXct(
  paste(
    data_citibike$year,
    data_citibike$month,
    data_citibike$day,
    data_citibike$hour
  ),
  format = "%Y %m %d %H"
)

data_citibike <- data_citibike[order(data_citibike$datetime), ]

# Training sample: April 1 - July 31, 2023 (4 months)
train_df <- subset(data_citibike,
                   datetime >= as.POSIXct("2023-04-01 00:00:00") &
                     datetime <= as.POSIXct("2023-07-31 23:00:00"))

# Test sample: August 2023
test_df  <- subset(data_citibike,
                   datetime >= as.POSIXct("2023-08-01 00:00:00") &
                     datetime <= as.POSIXct("2023-08-31 23:00:00"))


y_train <- ts(train_df$`area 6`, frequency = s)

# Full series (train + test), used in Part B
y_all <- c(train_df$`area 6`, test_df$`area 6`)


## ============================================================
## PART A: IDENTIFICATION AND ESTIMATION ON THE TRAINING SAMPLE
## ============================================================

## ---- 3. ACF/PACF plots for identifying p, q, P, Q ------------

## 3a. Raw series (relevant for D = 0)
acf(y_train, lag.max = 48, main = "ACF - area 6 (training)")
pacf(y_train, lag.max = 48, main = "PACF - area 6 (training)")

acf(y_train, lag.max = s * 5, main = "ACF - seasonal lags (s = 168)")
abline(v = seq(s, s * 5, by = s), col = "red", lty = 2)

pacf(y_train, lag.max = s * 5, main = "PACF - seasonal lags (s = 168)")
abline(v = seq(s, s * 5, by = s), col = "red", lty = 2)

## 3b. Seasonally differenced series (relevant for D = 1)
y_train_sdiff <- diff(y_train, lag = s, differences = 1)

acf(y_train_sdiff, lag.max = 48, main = "ACF - seasonally differenced (short lags)")
pacf(y_train_sdiff, lag.max = 48, main = "PACF - seasonally differenced (short lags)")

acf(y_train_sdiff, lag.max = s * 5, main = "ACF - seasonally differenced, seasonal lags")
abline(v = seq(s, s * 5, by = s), col = "red", lty = 2)

pacf(y_train_sdiff, lag.max = s * 5, main = "PACF - seasonally differenced, seasonal lags")
abline(v = seq(s, s * 5, by = s), col = "red", lty = 2)

## ---- 4. ARIMA models -------------------------------------------

## 4a. Manually specified ARIMA (non-seasonal).
## p, q chosen by AIC grid search.
p_range <- 0:3
q_range <- 0:3

grid_results_pq <- data.frame()
for (p_try in p_range) {
  for (q_try in q_range) {
    fit_try <- tryCatch(
      Arima(y_train, order = c(p_try, 0, q_try)),
      error = function(e) NULL
    )
    if (!is.null(fit_try) && is.finite(AIC(fit_try))) {
      grid_results_pq <- rbind(grid_results_pq,
                               data.frame(p = p_try, q = q_try, AIC = AIC(fit_try)))
    }
  }
}
if (nrow(grid_results_pq) == 0) stop("No ARIMA(p, 0, q) fit succeeded in the grid.")

grid_results_pq <- grid_results_pq[order(grid_results_pq$AIC), ]
print(grid_results_pq)   # lowest AIC at the top

p_hat <- grid_results_pq$p[1]
q_hat <- grid_results_pq$q[1]
cat("Selected p =", p_hat, ", q =", q_hat, "(lowest AIC in grid)\n")

fit_arima_manual <- Arima(y_train, order = c(p_hat, 0, q_hat))
summary(fit_arima_manual)
checkresiduals(fit_arima_manual)

## 4b. Automatic search over the ARIMA class (Hyndman & Khandakar, 2008)
fit_arima_auto <- auto.arima(y_train, d = 0, seasonal = FALSE)
summary(fit_arima_auto)
checkresiduals(fit_arima_auto)

## ---- 5. SARIMA models --------------------------------------------

grid_method <- "CSS"
P_range <- 0:3
Q_range <- 0:3

grid_criterion <- function(fit, D, P) {
  if (grid_method == "CSS") {
    n_eff <- length(y_train) - D * s - (p_hat + P * s)
    n_eff * log(fit$sigma2) + 2 * length(fit$coef)
  } else {
    AIC(fit)
  }
}


sarima_grid <- function(D) {
  cache_file <- paste0("grid_cache_s", s, "_", grid_method,
                       "_p", p_hat, "q", q_hat, "_D", D, ".rds")
  cache <- if (file.exists(cache_file)) readRDS(cache_file) else list()
  
  out      <- data.frame()
  best_fit <- NULL
  best_crit <- Inf
  for (P_try in P_range) {
    for (Q_try in Q_range) {
      key <- paste0("P", P_try, "_Q", Q_try)
      
      if (!is.null(cache[[key]])) {
        fit_try <- if (isFALSE(cache[[key]])) NULL else cache[[key]]
        cat("  grid D =", D, "- P =", P_try, ", Q =", Q_try, "loaded from cache\n")
      } else {
        t0 <- Sys.time()
        fit_try <- tryCatch(
          Arima(y_train, order = c(p_hat, 0, q_hat),
                seasonal = list(order = c(P_try, D, Q_try), period = s),
                method = grid_method),
          error = function(e) {
            cat("  [FAILED] D =", D, "P =", P_try, "Q =", Q_try, ":", conditionMessage(e), "\n")
            NULL
          }
        )
        cache[[key]] <- if (is.null(fit_try)) FALSE else fit_try
        saveRDS(cache, cache_file)
        el <- Sys.time() - t0
        cat("  grid D =", D, "- P =", P_try, ", Q =", Q_try, "done in",
            round(as.numeric(el, units = "mins"), 1), "min\n")
      }
      
      if (!is.null(fit_try)) {
        crit <- grid_criterion(fit_try, D, P_try)
        if (is.finite(crit)) {
          out <- rbind(out, data.frame(P = P_try, Q = Q_try, crit = crit))
          if (crit < best_crit) {
            best_crit <- crit
            best_fit  <- fit_try
          }
        }
      }
    }
  }
  if (nrow(out) == 0) stop("No SARIMA fit succeeded in the grid for D = ", D, ".")
  list(table = out[order(out$crit), ], fit = best_fit)
}

## 5a. Manually specified SARIMA, D = 0
grid_D0 <- sarima_grid(D = 0)
cat("\n=== SARIMA grid, D = 0 ===\n")
print(grid_D0$table)

P_hat_D0 <- grid_D0$table$P[1]
Q_hat_D0 <- grid_D0$table$Q[1]
cat("Selected P =", P_hat_D0, ", Q =", Q_hat_D0, "(lowest criterion in grid, D = 0)\n")

fit_sarima_manual_D0 <- grid_D0$fit
summary(fit_sarima_manual_D0)
checkresiduals(fit_sarima_manual_D0)

## 5b. Manually specified SARIMA, D = 1
grid_D1 <- sarima_grid(D = 1)
cat("\n=== SARIMA grid, D = 1 ===\n")
print(grid_D1$table)

P_hat_D1 <- grid_D1$table$P[1]
Q_hat_D1 <- grid_D1$table$Q[1]
cat("Selected P =", P_hat_D1, ", Q =", Q_hat_D1, "(lowest criterion in grid, D = 1)\n")

fit_sarima_manual_D1 <- grid_D1$fit
summary(fit_sarima_manual_D1)
checkresiduals(fit_sarima_manual_D1)

## 5c. Automatic SARIMA search - lets auto.arima choose D on its own,
## independent of the D = 0 / D = 1 assumptions above.
fit_sarima_auto <- auto.arima(y_train, d = 0, seasonal = TRUE)
summary(fit_sarima_auto)
checkresiduals(fit_sarima_auto)


## ============================================================
## PART B: OUT-OF-SAMPLE FORECAST EXERCISE
## ============================================================

## ---- 6. Model Specifications (taken from Part A) ---------------------
## Turns an auto.arima fit into an order / seasonal specification.
## arimaorder() returns (p, d, q) for a non-seasonal model and
## (p, d, q, P, D, Q, period) for a seasonal one.
spec_from_fit <- function(fit) {
  ao <- arimaorder(fit)
  seasonal <- NULL
  if (length(ao) >= 6 && any(ao[4:6] > 0)) {
    seasonal <- list(order = unname(ao[4:6]), period = s)
  }
  list(order = unname(ao[1:3]), seasonal = seasonal)
}

specs <- list(
  arima_manual     = list(order = c(p_hat, 0, q_hat), seasonal = NULL),
  arima_auto       = spec_from_fit(fit_arima_auto),
  sarima_manual_D0 = list(order    = c(p_hat, 0, q_hat),
                          seasonal = list(order = c(P_hat_D0, 0, Q_hat_D0), period = s)),
  sarima_manual_D1 = list(order    = c(p_hat, 0, q_hat),
                          seasonal = list(order = c(P_hat_D1, 1, Q_hat_D1), period = s)),
  sarima_auto      = spec_from_fit(fit_sarima_auto)
)
models <- names(specs)

format_spec <- function(sp) {
  txt <- paste0("(", paste(sp$order, collapse = ","), ")")
  if (!is.null(sp$seasonal)) {
    txt <- paste0(txt, "(", paste(sp$seasonal$order, collapse = ","), ")[", sp$seasonal$period, "]")
  }
  txt
}

cat("\n=== Models used in Part B (from Part A) ===\n")
for (m in models) cat(sprintf("%-18s %s\n", m, format_spec(specs[[m]])))
cat("\n")

## ---- 7. Setup & Forecast Origins ------------------------------------
n_train    <- length(train_df$`area 6`)
n_test     <- length(test_df$`area 6`)
n_days_all <- n_test / h       # 31 forecast origins in August

stopifnot(n_test %% h == 0)

TEST_MODE <- FALSE             # set to TRUE for a quick test of 3 days
test_days <- 3
n_days    <- if (TEST_MODE) test_days else n_days_all

make_storage <- function() matrix(NA, nrow = n_days, ncol = h)

## ---- 7b. Safe Forecast Helper Function ------------------------------
safe_forecast <- function(window_data, order, seasonal = NULL, h, label, day) {
  fit <- tryCatch(
    if (is.null(seasonal)) {
      Arima(window_data, order = order, method = "CSS")
    } else {
      Arima(window_data, order = order, seasonal = seasonal, method = "CSS")
    },
    error = function(e) {
      cat("  [FAILED]", label, "- day", day, ":", conditionMessage(e), "\n")
      NULL
    }
  )
  if (is.null(fit)) return(rep(NA_real_, h))
  forecast(fit, h = h)$mean
}

## ---- 7c. Sanity Check (optional, slow at s = 168) --------------------
## Refits the D = 1 SARIMA on the initial training set with CSS-ML
## (exact likelihood) to compare against the CSS fits used in the loops.
RUN_SANITY_CHECK <- FALSE
if (RUN_SANITY_CHECK) {
  check_fit <- Arima(
    ts(y_all[1:n_train], frequency = s),
    order    = specs$sarima_manual_D1$order,
    seasonal = specs$sarima_manual_D1$seasonal,
    method   = "CSS-ML"
  )
  cat("=== Sanity check: CSS-ML fit on the initial training set ===\n")
  print(check_fit)
  cat("\n")
}

## ---- 8. Rolling Window Loop -----------------------------------------
forecasts <- setNames(lapply(models, function(m) make_storage()), models)
actuals   <- make_storage()

t_start <- Sys.time()

for (day in 1:n_days) {
  end_idx   <- n_train + (day - 1) * h
  start_idx <- end_idx - n_train + 1
  window_data <- ts(y_all[start_idx:end_idx], frequency = s)
  
  actual_idx <- (end_idx + 1):(end_idx + h)
  actuals[day, ] <- y_all[actual_idx]
  
  for (m in models) {
    forecasts[[m]][day, ] <- safe_forecast(window_data,
                                           order    = specs[[m]]$order,
                                           seasonal = specs[[m]]$seasonal,
                                           h = h, label = m, day = day)
  }
  
  cat("Rolling - day", day, "of", n_days, "done\n")
}

elapsed <- Sys.time() - t_start
cat("\nTime for rolling window:", elapsed, units(elapsed), "\n")

forecasts_rolling <- forecasts
actuals_rolling   <- actuals

## ---- 9. Expanding Window Loop ---------------------------------------
forecasts <- setNames(lapply(models, function(m) make_storage()), models)
actuals   <- make_storage()

for (day in 1:n_days) {
  end_idx     <- n_train + (day - 1) * h
  window_data <- ts(y_all[1:end_idx], frequency = s)
  
  actual_idx  <- (end_idx + 1):(end_idx + h)
  actuals[day, ] <- y_all[actual_idx]
  
  for (m in models) {
    forecasts[[m]][day, ] <- safe_forecast(window_data,
                                           order    = specs[[m]]$order,
                                           seasonal = specs[[m]]$seasonal,
                                           h = h, label = m, day = day)
  }
  
  cat("Expanding - day", day, "of", n_days, "done\n")
}

forecasts_expanding <- forecasts
actuals_expanding   <- actuals

## ---- 10. Accuracy Summary (RMSE, MAE, MAPE) -------------------------
summarize_accuracy <- function(forecasts, actuals) {
  out <- data.frame()
  for (m in names(forecasts)) {
    err  <- actuals - forecasts[[m]]
    rmse <- sqrt(mean(err^2, na.rm = TRUE))
    mae  <- mean(abs(err), na.rm = TRUE)
    mape <- mean(abs(err / actuals), na.rm = TRUE) * 100
    out  <- rbind(out, data.frame(model = m, RMSE = rmse, MAE = mae, MAPE = mape))
  }
  out[order(out$RMSE), ]
}

cat("\n=== Rolling window accuracy (pooled across all horizons) ===\n")
print(summarize_accuracy(forecasts_rolling, actuals_rolling))

cat("\n=== Expanding window accuracy (pooled across all horizons) ===\n")
print(summarize_accuracy(forecasts_expanding, actuals_expanding))

## ---- 11. Accuracy per Forecast Horizon (h = 1 ... 24) ---------------
accuracy_by_horizon <- function(forecasts, actuals) {
  out <- data.frame()
  for (m in names(forecasts)) {
    for (hh in 1:h) {
      err  <- actuals[, hh] - forecasts[[m]][, hh]
      rmse <- sqrt(mean(err^2, na.rm = TRUE))
      mae  <- mean(abs(err), na.rm = TRUE)
      mape <- mean(abs(err / actuals[, hh]), na.rm = TRUE) * 100
      out  <- rbind(out, data.frame(model = m, horizon = hh, RMSE = rmse, MAE = mae, MAPE = mape))
    }
  }
  out
}

by_h_rolling   <- accuracy_by_horizon(forecasts_rolling, actuals_rolling)
by_h_expanding <- accuracy_by_horizon(forecasts_expanding, actuals_expanding)

cat("\n=== Rolling Window accuracy per horizon (hours 1-3 vs 22-24) ===\n")
print(subset(by_h_rolling, horizon %in% c(1, 2, 3, 22, 23, 24)))

# Graph: RMSE across horizons h = 1 to 24
plot(
  1, type = "n",
  xlim = c(1, h),
  ylim = range(by_h_rolling$RMSE, na.rm = TRUE),
  xlab = "Forecast horizon (hours ahead)",
  ylab = "RMSE",
  main = "Rolling window: RMSE per forecast horizon (s = 168)"
)

colors <- rainbow(length(models))
for (i in seq_along(models)) {
  m <- models[i]
  lines(1:h, by_h_rolling$RMSE[by_h_rolling$model == m], col = colors[i], lwd = 2)
}
legend("topleft", legend = models, col = colors, lwd = 2, cex = 0.7)

## ============================================================
## PART C: ASSESSING THE SIGNIFICANCE OF FORECAST RESULTS
##         Diebold-Mariano (DM) test and Model Confidence Set (MCS)
## ============================================================


## ---- 12. Packages ----------------------------------------------
## dm.test() is in the forecast package (already loaded).
## MCS: package "MCS" (Bernardi and Catania), function MCSprocedure().
if (!requireNamespace("MCS", quietly = TRUE)) install.packages("MCS")
library(MCS)

## ---- 13. Prepare forecast errors -------------------------------
get_errors <- function(forecasts, actuals) {
  E <- lapply(forecasts, function(f) as.vector(t(actuals - f)))
  ok <- complete.cases(do.call(cbind, E))
  if (any(!ok)) cat("Dropping", sum(!ok), "hours with missing forecasts\n")
  lapply(E, function(e) e[ok])
}

## Duplicates are removed and reported.
drop_duplicates <- function(E) {
  keep <- rep(TRUE, length(E))
  for (i in seq_along(E)) {
    if (!keep[i]) next
    for (j in seq_along(E)) {
      if (j > i && keep[j] && isTRUE(all.equal(E[[i]], E[[j]]))) {
        cat("Model", names(E)[j], "gives the same forecasts as", names(E)[i],
            "- dropped from the tests\n")
        keep[j] <- FALSE
      }
    }
  }
  E[keep]
}

errors_list <- list(
  rolling   = drop_duplicates(get_errors(forecasts_rolling,   actuals_rolling)),
  expanding = drop_duplicates(get_errors(forecasts_expanding, actuals_expanding))
)

## ---- 14. Pairwise Diebold-Mariano tests -------------------------
## Settings:
##   power = 2 : squared-error loss (matches RMSE); power = 1: absolute loss (MAE)
##   h_dm  = 24: forecast horizon used for the HAC variance. The hourly
##               errors within one 24-hour forecast are correlated, and
##               errors from consecutive days can be as well.
##   two-sided test: H0 equal accuracy, H1 different accuracy.
h_dm <- h

run_dm <- function(E, power = 2, h_dm = 24) {
  nm  <- names(E)
  out <- data.frame()
  for (i in 1:(length(nm) - 1)) {
    for (j in (i + 1):length(nm)) {
      res <- tryCatch(
        dm.test(E[[i]], E[[j]], alternative = "two.sided",
                h = h_dm, power = power),
        error = function(e) NULL
      )
      loss_i <- mean(abs(E[[i]])^power)
      loss_j <- mean(abs(E[[j]])^power)
      out <- rbind(out, data.frame(
        model_1   = nm[i],
        model_2   = nm[j],
        loss_1    = loss_i,
        loss_2    = loss_j,
        DM_stat   = if (is.null(res)) NA else unname(res$statistic),
        p_value   = if (is.null(res)) NA else res$p.value,
        better    = if (loss_i < loss_j) nm[i] else nm[j]
      ))
    }
  }
  out$p_holm <- p.adjust(out$p_value, method = "holm")   # 10 pairs -> adjust
  out
}

dm_results <- list()
for (w in names(errors_list)) {
  for (pw in c(2, 1)) {
    key <- paste0(w, "_", if (pw == 2) "squared" else "absolute")
    dm_results[[key]] <- run_dm(errors_list[[w]], power = pw, h_dm = h_dm)
    cat("\n=== DM tests:", w, "window,",
        if (pw == 2) "squared-error" else "absolute-error", "loss ===\n")
    print(dm_results[[key]], digits = 4, row.names = FALSE)
  }
}

## p-value matrix (unadjusted), for a quick overview. Upper triangle only.
dm_matrix <- function(res) {
  nm <- unique(c(res$model_1, res$model_2))
  M  <- matrix(NA, length(nm), length(nm), dimnames = list(nm, nm))
  for (r in seq_len(nrow(res))) M[res$model_1[r], res$model_2[r]] <- res$p_value[r]
  round(M, 4)
}
cat("\n=== DM p-value matrix: rolling window, squared loss (unadjusted) ===\n")
print(dm_matrix(dm_results$rolling_squared))

## Robustness check: one loss value per forecast origin (daily mean loss).
## Each origin covers 24 different hours, so the 31 daily losses are
## less dependent than the hourly ones. 
run_dm_daily <- function(E, power = 2) {
  nm  <- names(E)
  D   <- sapply(E, function(e) rowMeans(matrix(abs(e)^power, ncol = h, byrow = TRUE)))
  out <- data.frame()
  for (i in 1:(length(nm) - 1)) {
    for (j in (i + 1):length(nm)) {
      res <- tryCatch(dm.test(D[, i], D[, j], h = 1, power = 1,
                              alternative = "two.sided"),
                      error = function(e) NULL)
      out <- rbind(out, data.frame(model_1 = nm[i], model_2 = nm[j],
                                   p_value = if (is.null(res)) NA else res$p.value))
    }
  }
  out$p_holm <- p.adjust(out$p_value, method = "holm")
  out
}
cat("\n=== DM tests on daily mean squared error (rolling window) ===\n")
print(run_dm_daily(errors_list$rolling, power = 2), digits = 4, row.names = FALSE)

## ---- 15. Model Confidence Set -----------------------------------
## Loss matrix: one row per hour, one column per model.
## alpha = 0.10 (a common choice in the MCS literature).
## statistic = "Tmax" (alternative: "TR").
## k = 24 : block length of the block bootstrap. It should be at least as
##          long as the dependence in the losses; 24 hours is a
##          judgement call tied to the daily forecast cycle. 
## B = 5000 bootstrap replications.
mcs_alpha <- 0.10
mcs_B     <- 5000
mcs_k     <- 24

set.seed(123)
mcs_results <- list()
for (w in names(errors_list)) {
  for (pw in c(2, 1)) {
    key  <- paste0(w, "_", if (pw == 2) "squared" else "absolute")
    Loss <- sapply(errors_list[[w]], function(e) abs(e)^pw)
    cat("\n=== MCS:", w, "window,",
        if (pw == 2) "squared-error" else "absolute-error", "loss ===\n")
    mcs_results[[key]] <- tryCatch(
      MCSprocedure(Loss = Loss, alpha = mcs_alpha, B = mcs_B,
                   statistic = "Tmax", k = mcs_k, verbose = FALSE),
      error = function(e) { cat("MCS failed:", conditionMessage(e), "\n"); NULL }
    )
    if (!is.null(mcs_results[[key]])) print(mcs_results[[key]])
  }
}

## Models that remain in the MCS. This uses the @show slot of the MCS
## object, whose row names are the surviving models. 
mcs_members <- function(x) {
  if (is.null(x)) return(NA_character_)
  tryCatch(rownames(x@show), error = function(e) NA_character_)
}

## ---- 16. Summary table: accuracy, DM and MCS -------------------
cat("\n=== Summary ===\n")
for (key in names(dm_results)) {
  res  <- dm_results[[key]]
  w    <- strsplit(key, "_")[[1]][1]
  pw   <- if (grepl("squared", key)) 2 else 1
  mean_loss <- sapply(errors_list[[w]], function(e) mean(abs(e)^pw))
  cat("\n--", key, "--\n")
  cat("Lowest mean loss:", names(which.min(mean_loss)), "\n")
  cat("Pairs significant at 5% (Holm-adjusted):",
      sum(res$p_holm < 0.05, na.rm = TRUE), "of", nrow(res), "\n")
  cat("Pairs significant at 5% (unadjusted):",
      sum(res$p_value < 0.05, na.rm = TRUE), "of", nrow(res), "\n")
  cat("Models in the", 100 * (1 - mcs_alpha), "% MCS:",
      paste(mcs_members(mcs_results[[key]]), collapse = ", "), "\n")
}



## ---- 17. Final conclusion (printed at the end of the run) ----------
## Collects accuracy measures, DM results and MCS membership in one place.
## Measures: MSE, RMSE, MAE, MAPE and MASE.

mase_scale <- mean(abs(diff(as.numeric(y_train), lag = s)))

full_accuracy <- function(forecasts, actuals) {
  out <- data.frame()
  for (m in names(forecasts)) {
    err <- actuals - forecasts[[m]]
    mse <- mean(err^2, na.rm = TRUE)
    mae <- mean(abs(err), na.rm = TRUE)
    out <- rbind(out, data.frame(
      model = m,
      spec  = format_spec(specs[[m]]),
      MSE   = mse,
      RMSE  = sqrt(mse),
      MAE   = mae,
      MAPE  = mean(abs(err / actuals), na.rm = TRUE) * 100,
      MASE  = mae / mase_scale
    ))
  }
  out[order(out$RMSE), ]
}

## Maps a model to the kept model with identical forecasts.
representative <- function(m, w) {
  kept <- names(errors_list[[w]])
  if (m %in% kept) return(m)
  fc <- if (w == "rolling") forecasts_rolling else forecasts_expanding
  for (k in kept) if (isTRUE(all.equal(fc[[m]], fc[[k]]))) return(k)
  m
}

conclude_window <- function(w) {
  fc  <- if (w == "rolling") forecasts_rolling else forecasts_expanding
  act <- if (w == "rolling") actuals_rolling   else actuals_expanding
  acc <- full_accuracy(fc, act)
  
  cat("\n", strrep("=", 70), "\n", sep = "")
  cat(toupper(w), "WINDOW\n")
  cat(strrep("=", 70), "\n", sep = "")
  
  cat("\nAccuracy (sorted by RMSE; pooled over all 24 horizons, August 2023):\n")
  print(acc, digits = 4, row.names = FALSE)
  
  best_rmse <- acc$model[which.min(acc$RMSE)]
  best_mae  <- acc$model[which.min(acc$MAE)]
  best_mape <- acc$model[which.min(acc$MAPE)]
  best_mase <- acc$model[which.min(acc$MASE)]
  cat("\nBest by RMSE:", best_rmse, "\n")
  cat("Best by MAE :", best_mae, "\n")
  cat("Best by MAPE:", best_mape, "\n")
  cat("Best by MASE:", best_mase, "\n")
  
  ## DM: best-by-RMSE model against every other model (squared loss)
  b   <- representative(best_rmse, w)
  res <- dm_results[[paste0(w, "_squared")]]
  sub <- res[res$model_1 == b | res$model_2 == b, ]
  sub$opponent <- ifelse(sub$model_1 == b, sub$model_2, sub$model_1)
  sub$best_has_lower_loss <- sub$better == b
  cat("\nDM tests (squared loss), ", best_rmse, " versus each other model:\n", sep = "")
  print(sub[, c("opponent", "p_value", "p_holm", "best_has_lower_loss")],
        digits = 4, row.names = FALSE)
  n_sig <- sum(sub$p_holm < 0.05 & sub$best_has_lower_loss, na.rm = TRUE)
  cat("Significantly better than", n_sig, "of", nrow(sub),
      "competitors (5% level, Holm-adjusted).\n")
  
  ## MCS membership under both losses
  for (pw in c("squared", "absolute")) {
    mem <- mcs_members(mcs_results[[paste0(w, "_", pw)]])
    cat("\nMCS (", pw, " loss, ", 100 * (1 - mcs_alpha), "%): ",
        paste(mem, collapse = ", "), "\n", sep = "")
    if (!any(is.na(mem)))
      cat("  ", best_rmse, if (b %in% mem) "IS" else "is NOT", "in the MCS.\n")
  }
  invisible(list(acc = acc, best = best_rmse))
}

cat("\n\n", strrep("#", 70), "\n", sep = "")
cat("FINAL CONCLUSION\n")
cat(strrep("#", 70), "\n", sep = "")

concl_roll <- conclude_window("rolling")
concl_exp  <- conclude_window("expanding")

cat("\n", strrep("=", 70), "\n", sep = "")
cat("OVERALL\n")
cat(strrep("=", 70), "\n", sep = "")
if (concl_roll$best == concl_exp$best) {
  cat("The same model has the lowest RMSE under both window schemes:",
      concl_roll$best, "\n")
} else {
  cat("The best model by RMSE differs by window scheme: rolling =",
      concl_roll$best, ", expanding =", concl_exp$best, "\n")
}
cat("Specification of the rolling-window winner:",
    format_spec(specs[[concl_roll$best]]), "\n")
cat("Read the DM and MCS output above before calling any model 'best':\n",
    "a lowest RMSE only counts as a clear win if the DM tests against the\n",
    "other models are significant (Holm-adjusted) and the MCS is small.\n", sep = "")


results_all <- list(
  forecasts_rolling = forecasts_rolling,
  forecasts_expanding = forecasts_expanding,
  actuals_rolling = actuals_rolling,
  actuals_expanding = actuals_expanding
)

saveRDS(results_all, file = file.path(results_dir, "all_forecast_results.rds"))


## ============================================================
## PART C (continued): COMPARISON AGAINST THE ETS BENCHMARK
##   Diebold-Mariano tests and MCS including ETS
## ============================================================
## Uses the forecasts that are already in memory from Part B
## (no need to reload all_forecast_results.rds).
## The ETS forecasts come from two csv files with the columns
## day, horizon, actual, ets_forecast.

fc_rolling        <- forecasts_rolling
fc_expanding      <- forecasts_expanding

## ---- 18. Load ETS forecasts --------------------------------------
ets_files <- c("rolling_ets_mcs.csv", "expanding_ets_mcs.csv")
ets_missing <- ets_files[!file.exists(ets_files)]
if (length(ets_missing) > 0) {
  stop("ETS file(s) not found in ", getwd(), ": ",
       paste(ets_missing, collapse = ", "))
}

read_ets <- function(file) {
  ets <- read.csv(file, check.names = FALSE)
  
  if (!all(c("day", "horizon", "actual", "ets_forecast") %in% names(ets))) {
    ets <- read.delim(file, check.names = FALSE)
  }
  
  names(ets) <- trimws(tolower(names(ets)))
  
  required <- c("day", "horizon", "actual", "ets_forecast")
  
  if (!all(required %in% names(ets))) {
    stop(paste("Check the ETS file:", file))
  }
  
  ets <- ets[order(ets$day, ets$horizon), ]
  ets
}

ets_rolling   <- read_ets("rolling_ets_mcs.csv")
ets_expanding <- read_ets("expanding_ets_mcs.csv")


## ---- 19. DM tests: each model against ETS ------------------------
## (named run_dm_ets so it does not overwrite run_dm from section 14)
run_dm_ets <- function(fc, actuals, ets, power) {
  
  my_actuals <- as.vector(t(as.matrix(actuals)))
  
  if (nrow(ets) != length(my_actuals)) {
    stop("The number of ETS forecasts does not match your saved forecasts.")
  }
  
  if (!isTRUE(all.equal(
    as.numeric(my_actuals),
    as.numeric(ets$actual),
    tolerance = 1e-8
  ))) {
    stop("Actual values do not match or are in a different order.")
  }
  
  e_ets <- as.numeric(ets$actual - ets$ets_forecast)
  
  results <- lapply(names(fc), function(m) {
    
    e_model <- as.vector(t(as.matrix(actuals - fc[[m]])))
    
    valid <- is.finite(e_model) & is.finite(e_ets)
    
    test <- dm.test(
      e1 = e_model[valid],
      e2 = e_ets[valid],
      alternative = "two.sided",
      h = 24,
      power = power
    )
    
    data.frame(
      Model = m,
      Loss_model = mean(abs(e_model[valid])^power),
      Loss_ETS = mean(abs(e_ets[valid])^power),
      DM_statistic = as.numeric(test$statistic),
      p_value = test$p.value
    )
  })
  
  results <- do.call(rbind, results)
  results$p_holm <- p.adjust(results$p_value, method = "holm")
  
  results[order(results$p_holm), ]
}

# ROLLING: Squared-error loss and absolute-error loss
DM_rolling_squared <- run_dm_ets(
  fc_rolling, actuals_rolling, ets_rolling, power = 2
)

DM_rolling_absolute <- run_dm_ets(
  fc_rolling, actuals_rolling, ets_rolling, power = 1
)

# EXPANDING: Squared-error loss and absolute-error loss
DM_expanding_squared <- run_dm_ets(
  fc_expanding, actuals_expanding, ets_expanding, power = 2
)

DM_expanding_absolute <- run_dm_ets(
  fc_expanding, actuals_expanding, ets_expanding, power = 1
)

# Display results
cat("\nDM TEST VS ETS: ROLLING - SQUARED-ERROR LOSS\n")
print(DM_rolling_squared, row.names = FALSE)

cat("\nDM TEST VS ETS: ROLLING - ABSOLUTE-ERROR LOSS\n")
print(DM_rolling_absolute, row.names = FALSE)

cat("\nDM TEST VS ETS: EXPANDING - SQUARED-ERROR LOSS\n")
print(DM_expanding_squared, row.names = FALSE)

cat("\nDM TEST VS ETS: EXPANDING - ABSOLUTE-ERROR LOSS\n")
print(DM_expanding_absolute, row.names = FALSE)

# Save results (results folder)
write.csv(DM_rolling_squared, file.path(results_dir, "DM_ETS_rolling_squared.csv"),
          row.names = FALSE)
write.csv(DM_rolling_absolute, file.path(results_dir, "DM_ETS_rolling_absolute.csv"),
          row.names = FALSE)

write.csv(DM_expanding_squared, file.path(results_dir, "DM_ETS_expanding_squared.csv"),
          row.names = FALSE)
write.csv(DM_expanding_absolute, file.path(results_dir, "DM_ETS_expanding_absolute.csv"),
          row.names = FALSE)


## ---- 20. MCS including the ETS benchmark -------------------------
## (named run_mcs_ets so it does not clash with the MCS objects from
## section 15; the settings below are those of the original MCS script)
run_mcs_ets <- function(fc, actuals, ets, window_name) {
  
  actuals <- as.matrix(actuals)
  
  # Check actual values match
  stopifnot(nrow(ets) == length(as.vector(t(actuals))))
  
  stopifnot(isTRUE(all.equal(
    as.numeric(as.vector(t(actuals))),
    as.numeric(ets$actual),
    tolerance = 1e-8
  )))
  
  # Create model forecast errors
  errors <- sapply(names(fc), function(m) {
    as.vector(t(actuals - as.matrix(fc[[m]])))
  })
  
  # Add ETS forecasts for the matching window
  errors <- cbind(
    errors,
    ETS = as.numeric(ets$actual - ets$ets_forecast)
  )
  
  # Create loss matrices
  loss_squared <- errors^2
  loss_absolute <- abs(errors)
  
  # MCS for squared-error loss
  set.seed(94991)
  mcs_squared <- MCSprocedure(
    Loss = loss_squared,
    alpha = 0.10,
    B = 5000,
    statistic = "Tmax"
  )
  
  # MCS for absolute-error loss
  set.seed(94991)
  mcs_absolute <- MCSprocedure(
    Loss = loss_absolute,
    alpha = 0.10,
    B = 5000,
    statistic = "Tmax"
  )
  
  # Print results
  cat("\n", toupper(window_name), "- SQUARED LOSS\n", sep = "")
  print(mcs_squared)
  
  cat("\n", toupper(window_name), "- ABSOLUTE LOSS\n", sep = "")
  print(mcs_absolute)
  
  # Save results (results folder)
  write.csv(
    mcs_squared@show,
    file.path(results_dir, paste0("MCS_ETS_", window_name, "_squared.csv")),
    row.names = FALSE
  )
  
  write.csv(
    mcs_absolute@show,
    file.path(results_dir, paste0("MCS_ETS_", window_name, "_absolute.csv")),
    row.names = FALSE
  )
  
  saveRDS(
    list(
      squared = mcs_squared,
      absolute = mcs_absolute
    ),
    file.path(results_dir, paste0("MCS_ETS_", window_name, "_both_losses.rds"))
  )
  
  return(list(
    squared = mcs_squared,
    absolute = mcs_absolute
  ))
}

# ROLLING MCS (with ETS)
mcs_rolling <- run_mcs_ets(
  fc = fc_rolling,
  actuals = actuals_rolling,
  ets = ets_rolling,
  window_name = "rolling"
)

# EXPANDING MCS (with ETS)
mcs_expanding <- run_mcs_ets(
  fc = fc_expanding,
  actuals = actuals_expanding,
  ets = ets_expanding,
  window_name = "expanding"
)


## ============================================================
## PART D: SAVE ALL RESULTS
## ============================================================

# Create a folder in the current working directory
results_dir <- "results"
dir.create(results_dir, showWarnings = FALSE, recursive = TRUE)

# ------------------------------------------------------------
# 1. Collect the main results into one R object
# ------------------------------------------------------------

all_results <- list(
  
  # Model specifications
  model_specs = specs,
  grid_results_arima = grid_results_pq,
  grid_results_sarima_D0 = grid_D0$table,
  grid_results_sarima_D1 = grid_D1$table,
  
  # Fitted models from Part A
  fit_arima_manual = fit_arima_manual,
  fit_arima_auto = fit_arima_auto,
  fit_sarima_manual_D0 = fit_sarima_manual_D0,
  fit_sarima_manual_D1 = fit_sarima_manual_D1,
  fit_sarima_auto = fit_sarima_auto,
  
  # Forecasts and actual values
  forecasts_rolling = forecasts_rolling,
  forecasts_expanding = forecasts_expanding,
  actuals_rolling = actuals_rolling,
  actuals_expanding = actuals_expanding,
  
  # Accuracy results
  accuracy_rolling = full_accuracy(
    forecasts_rolling, actuals_rolling
  ),
  accuracy_expanding = full_accuracy(
    forecasts_expanding, actuals_expanding
  ),
  accuracy_by_horizon_rolling = by_h_rolling,
  accuracy_by_horizon_expanding = by_h_expanding,
  
  # Statistical tests
  errors_list = errors_list,
  dm_results = dm_results,
  mcs_results = mcs_results,
  
  # Tests including the ETS benchmark
  dm_ets_results = list(
    rolling_squared    = DM_rolling_squared,
    rolling_absolute   = DM_rolling_absolute,
    expanding_squared  = DM_expanding_squared,
    expanding_absolute = DM_expanding_absolute
  ),
  mcs_ets_results = list(rolling = mcs_rolling, expanding = mcs_expanding),
  
  # Final conclusions
  conclusion_rolling = concl_roll,
  conclusion_expanding = concl_exp
)

# Save everything in one RDS file
saveRDS(
  all_results,
  file = file.path(results_dir, "all_results.rds")
)


# ------------------------------------------------------------
# 2. Export accuracy tables as CSV files
# ------------------------------------------------------------

write.csv(
  all_results$accuracy_rolling,
  file.path(results_dir, "accuracy_rolling.csv"),
  row.names = FALSE
)

write.csv(
  all_results$accuracy_expanding,
  file.path(results_dir, "accuracy_expanding.csv"),
  row.names = FALSE
)

write.csv(
  by_h_rolling,
  file.path(results_dir, "accuracy_by_horizon_rolling.csv"),
  row.names = FALSE
)

write.csv(
  by_h_expanding,
  file.path(results_dir, "accuracy_by_horizon_expanding.csv"),
  row.names = FALSE
)


# ------------------------------------------------------------
# 3. Export all DM test results as separate CSV files
# ------------------------------------------------------------

for (nm in names(dm_results)) {
  
  write.csv(
    dm_results[[nm]],
    file.path(results_dir, paste0("DM_", nm, ".csv")),
    row.names = FALSE
  )
}


# ------------------------------------------------------------
# 4. Export MCS membership in a readable CSV file
# ------------------------------------------------------------

mcs_summary <- data.frame(
  analysis = names(mcs_results),
  models_in_MCS = vapply(
    mcs_results,
    function(x) paste(mcs_members(x), collapse = ", "),
    character(1)
  )
)

write.csv(
  mcs_summary,
  file.path(results_dir, "MCS_membership.csv"),
  row.names = FALSE
)


# ------------------------------------------------------------
# 4b. Export MCS membership including the ETS benchmark
# ------------------------------------------------------------

mcs_ets_summary <- data.frame(
  analysis = c("rolling_squared", "rolling_absolute",
               "expanding_squared", "expanding_absolute"),
  models_in_MCS = c(
    paste(mcs_members(mcs_rolling$squared),    collapse = ", "),
    paste(mcs_members(mcs_rolling$absolute),   collapse = ", "),
    paste(mcs_members(mcs_expanding$squared),  collapse = ", "),
    paste(mcs_members(mcs_expanding$absolute), collapse = ", ")
  )
)

write.csv(
  mcs_ets_summary,
  file.path(results_dir, "MCS_ETS_membership.csv"),
  row.names = FALSE
)


# ------------------------------------------------------------
# 5. Export model specifications
# ------------------------------------------------------------

model_specifications <- data.frame(
  model = names(specs),
  specification = vapply(
    specs,
    format_spec,
    character(1)
  )
)

write.csv(
  model_specifications,
  file.path(results_dir, "model_specifications.csv"),
  row.names = FALSE
)


# ------------------------------------------------------------
# 6. Save a short text summary of the main conclusions
# ------------------------------------------------------------

summary_lines <- c(
  "FINAL FORECAST RESULTS",
  "======================",
  "",
  paste(
    "Rolling-window best model by RMSE:",
    concl_roll$best
  ),
  paste(
    "Expanding-window best model by RMSE:",
    concl_exp$best
  ),
  "",
  "Rolling-window accuracy:",
  capture.output(print(all_results$accuracy_rolling, row.names = FALSE)),
  "",
  "Expanding-window accuracy:",
  capture.output(print(all_results$accuracy_expanding, row.names = FALSE)),
  "",
  "MCS membership (ARIMA/SARIMA models only):",
  capture.output(print(mcs_summary, row.names = FALSE)),
  "",
  "MCS membership (including ETS):",
  capture.output(print(mcs_ets_summary, row.names = FALSE))
)

writeLines(
  summary_lines,
  con = file.path(results_dir, "results_summary.txt")
)


# ------------------------------------------------------------
# 7. Confirm where the results were saved
# ------------------------------------------------------------

cat(
  "\nAll results have been saved to:",
  normalizePath(results_dir),
  "\n"
)

print(list.files(results_dir))