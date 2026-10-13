library(forecast)
library(FoReco)
library(FoRecoML)
library(ggplot2)

# ---- 0. Settings ----
# Set the working directory to the folder containing the data files first:
#setwd("path/to/folder")
data_file  <- "citibike.RData"
recon_file <- "ReconData.RData"
run_sarima <- FALSE   # TRUE: refit all SARIMA models (slow); FALSE: load saved results

load(data_file)

set.seed(123)
s <- 168    # weekly seasonality in hourly data
h <- 168    # forecast horizon: one week
t_total <- Sys.time()

# ---- 1. Data prep ----
data_citibike$datetime <- as.POSIXct(
  paste(data_citibike$year, data_citibike$month, data_citibike$day, data_citibike$hour),
  format = "%Y %m %d %H")
data_citibike <- data_citibike[order(data_citibike$datetime), ]

train_df <- subset(data_citibike, datetime >= as.POSIXct("2023-04-01 00:00:00") &
                     datetime <= as.POSIXct("2023-07-31 23:00:00"))
test_df  <- subset(data_citibike, datetime >= as.POSIXct("2023-08-01 00:00:00") &
                     datetime <= as.POSIXct("2023-08-07 23:00:00"))  # test week

# ---- 2. Hierarchy ----
# y9 = total, y7 = areas 1-4, y8 = areas 5-6
agg_mat <- rbind(
  c(1,1,1,1,1,1),   # y9
  c(1,1,1,1,0,0),   # y7
  c(0,0,0,0,1,1)    # y8
)
colnames(agg_mat) <- paste0("y", 1:6)
rownames(agg_mat) <- c("y9", "y7", "y8")

bottom    <- colnames(agg_mat)
all_names <- c(rownames(agg_mat), bottom)   # FoReco order: upper first
ord       <- paste0("y", 1:9)               # display order: y1 ... y9

make_all <- function(df) {
  bts <- as.matrix(df[, paste("area", 1:6)])
  colnames(bts) <- bottom
  out <- cbind(bts %*% t(agg_mat), bts)     # upper series = sum of bottom series
  colnames(out) <- all_names
  out
}
train_all <- make_all(train_df)
test_all  <- make_all(test_df)

# ---- 3./4. Base forecasts, residuals and validation forecasts ----
# Slow part. By default the saved results are loaded instead.
if (run_sarima) {
  
  # ---- 3. One SARIMA (2,0,3)(2,0,2)[168] per series, CSS ----
  t3 <- Sys.time()
  y_order    <- c(2, 0, 3)    # (p, d, q)
  y_seasonal <- c(2, 0, 2)    # (P, D, Q), period = 168
  
  fit_one <- function(x, nm, h) {
    x <- ts(x, frequency = s)
    fit <- Arima(x,
                 order    = y_order,
                 seasonal = list(order = y_seasonal, period = s),
                 method   = "CSS")
    cat("  ", nm, "done:", arimaorder(fit)[1:3], "\n")
    list(fit = fit,
         fc  = as.numeric(forecast(fit, h = h)$mean),
         res = as.numeric(residuals(fit, type = "response")))
  }
  
  fits <- setNames(lapply(all_names, function(nm) fit_one(train_all[, nm], nm, h)), all_names)
  yhat <- sapply(fits, `[[`, "fc");  colnames(yhat) <- all_names   # 168 x 9
  res  <- sapply(fits, `[[`, "res"); colnames(res)  <- all_names   # T x 9
  res  <- res[-(1:s), ]    # drop first 168 residuals (burn-in from seasonal differencing)
  cat("Section 3 (fit 9 SARIMAs):",
      round(as.numeric(Sys.time() - t3, units = "secs"), 1), "s\n")
  
  # SARIMA(p,d,q)(P,D,Q)[s] configuration of each series
  sarima_label <- function(fit) {
    a <- fit$arma   # c(p, q, P, Q, period, d, D)
    sprintf("SARIMA(%d,%d,%d)(%d,%d,%d)[%d]", a[1], a[6], a[2], a[3], a[7], a[4], a[5])
  }
  sarima_cfg <- data.frame(series = ord,
                           model  = sapply(fits[ord], function(f) sarima_label(f$fit)),
                           row.names = NULL)
  
  # ---- 4. Validation forecasts (only needed for the ML method) ----
  t4 <- Sys.time()
  n_val <- 1
  T_tr  <- nrow(train_all)
  val_hat_list <- val_obs_list <- vector("list", n_val)
  
  total_steps <- n_val * length(all_names)
  step <- 0
  pb <- txtProgressBar(min = 0, max = total_steps, style = 3)
  
  for (d in 1:n_val) {
    end_idx <- T_tr - h * (n_val - d + 1)
    idx     <- (end_idx + 1):(end_idx + h)
    
    hat <- matrix(NA_real_, nrow = h, ncol = length(all_names),
                  dimnames = list(NULL, all_names))
    for (nm in all_names) {
      x <- ts(train_all[1:end_idx, nm], frequency = s)
      hat[, nm] <- as.numeric(forecast(Arima(x,
                                             model = fits[[nm]]$fit,
                                             fixed = fits[[nm]]$fit$coef),
                                       h = h)$mean)
      gc(verbose = FALSE)   # free the large state-space objects right away
      
      step <- step + 1
      setTxtProgressBar(pb, step)
    }
    val_hat_list[[d]] <- hat
    val_obs_list[[d]] <- train_all[idx, bottom]
  }
  close(pb)
  
  val_hat <- do.call(rbind, val_hat_list); colnames(val_hat) <- all_names
  val_obs <- do.call(rbind, val_obs_list)
  cat("\nSection 4 (validation loop):",
      round(as.numeric(Sys.time() - t4, units = "secs"), 1), "s\n")
  
  save(yhat, res, val_hat, val_obs, sarima_cfg, file = recon_file)
  
} else {
  load(recon_file)
}

cat("\n=== SARIMA configurations ===\n")
print(sarima_cfg, row.names = FALSE)

# ---- 5. Reconciliation ----
t5 <- Sys.time()
comb_ls <- "wls"

reco_bu <- csbu(base = yhat[, bottom], agg_mat = agg_mat)                        # bottom-up
reco_ls <- csrec(base = yhat, agg_mat = agg_mat, res = res, comb = comb_ls)      # regression-based
fit_ml  <- csrml_fit(hat = val_hat, obs = val_obs, agg_mat = agg_mat,
                     approach = "randomForest")                                  # machine learning
reco_ml <- csrml(base = yhat, fit = fit_ml, agg_mat = agg_mat)

to_mat <- function(x) { m <- as.matrix(x); colnames(m) <- all_names; m }
bu_mat <- to_mat(reco_bu); ls_mat <- to_mat(reco_ls); ml_mat <- to_mat(reco_ml)
cat("Section 5 (reconciliation):",
    round(as.numeric(Sys.time() - t5, units = "secs"), 1), "s\n")

# ---- 6. Weight matrices G (6 x 9, columns y1 ... y9) ----
S <- rbind(agg_mat, diag(6))
colnames(S) <- bottom; rownames(S) <- all_names

# Bottom-up
G_bu <- cbind(matrix(0, 6, 3), diag(6))
colnames(G_bu) <- all_names; rownames(G_bu) <- bottom

# WLS: reconcile the 9 unit vectors, the bottom part gives the columns of G
I9 <- diag(9); colnames(I9) <- all_names
G_ls <- t(to_mat(csrec(base = I9, agg_mat = agg_mat, res = res, comb = comb_ls))[, bottom])
colnames(G_ls) <- all_names; rownames(G_ls) <- bottom

# ML: no fixed G exists, so this is a linear approximation of the forest
b_ml <- to_mat(csrml(base = val_hat, fit = fit_ml, agg_mat = agg_mat))[, bottom]
G_ml <- t(lm.fit(val_hat, b_ml)$coefficients)
colnames(G_ml) <- all_names; rownames(G_ml) <- bottom
r2_ml <- 1 - colSums((b_ml - val_hat %*% t(G_ml))^2) /
  colSums(scale(b_ml, scale = FALSE)^2)

# Reorder columns to y1 ... y9
#G_bu <- G_bu[, ord]; G_ls <- G_ls[, ord]; G_ml <- G_ml[, ord]

cat("\n=== G bottom-up ===\n");                   print(round(G_bu, 3))
cat("\n=== G WLS ===\n");                         print(round(G_ls, 3))
cat("\n=== G ML (linear approximation) ===\n");   print(round(G_ml, 3))
cat("\nR^2 of the ML approximation:\n");           print(round(r2_ml, 3))

# Check: S G yhat reproduces the package output (both should be TRUE)
print(all.equal(as.numeric(S %*% G_bu %*% yhat[1, ord]), as.numeric(bu_mat[1, ])))
print(all.equal(as.numeric(S %*% G_ls %*% yhat[1, ord]), as.numeric(ls_mat[1, ])))

# ---- 7. First time point (1 Aug 2023, 00:00): base vs reconciled vs actual ----
first_tbl <- rbind(Base               = yhat[1, ord],
                   `Bottom-up`        = bu_mat[1, ord],
                   WLS                = ls_mat[1, ord],
                   `Machine learning` = ml_mat[1, ord],
                   Actual             = test_all[1, ord])
cat("\n=== First time point ===\n")
print(round(first_tbl, 1))

# Coherence gap at the first time point (0 = coherent)
gap <- function(m) c(`y9-(y7+y8)`  = m[1, "y9"] - m[1, "y7"] - m[1, "y8"],
                     `y7-(y1..y4)` = m[1, "y7"] - sum(m[1, paste0("y", 1:4)]),
                     `y8-(y5+y6)`  = m[1, "y8"] - sum(m[1, paste0("y", 5:6)]))
cat("\n=== Coherence gap ===\n")
print(round(rbind(Base = gap(yhat), BU = gap(bu_mat), WLS = gap(ls_mat), ML = gap(ml_mat)), 2))

# ---- 8. RMSE over the test week (ratio to base, <1 = better) ----
rmse <- function(f, a) sqrt(colMeans((f - a)^2))
rmse_tbl <- rbind(Base = rmse(yhat,   test_all),
                  BU   = rmse(bu_mat, test_all),
                  WLS  = rmse(ls_mat, test_all),
                  ML   = rmse(ml_mat, test_all))[, ord]
ratio <- sweep(rmse_tbl[-1, ], 2, rmse_tbl["Base", ], "/")
cat("\n=== RMSE ===\n");         print(round(rmse_tbl, 1))
cat("\n=== RMSE ratio ===\n");   print(round(ratio, 3))

# ---- 9. Plots ----
panel_labels <- c(paste0("y", 1:6), "y7 (areas 1-4)", "y8 (areas 5-6)", "y9 (all areas)")
methods      <- c("Base", "Bottom-up", "WLS", "Machine learning")
pal          <- setNames(scales::hue_pal()(4), methods)   # one colour per method

# 9a. First time point: bars = methods, red line = actual
plot_df <- data.frame(series = rep(ord, 4), method = rep(methods, each = 9),
                      value  = c(yhat[1, ord], bu_mat[1, ord], ls_mat[1, ord], ml_mat[1, ord]))
actual_df <- data.frame(series = ord, actual = test_all[1, ord])
plot_df$series   <- factor(plot_df$series,   levels = ord, labels = panel_labels)
actual_df$series <- factor(actual_df$series, levels = ord, labels = panel_labels)
plot_df$method   <- factor(plot_df$method, levels = rev(methods))

print(
  ggplot(plot_df, aes(y = method, x = value, fill = method)) +
    geom_col(width = 0.7) +
    geom_vline(data = actual_df, aes(xintercept = actual, colour = "actual"),
               linewidth = 0.9) +
    scale_colour_manual(values = c(actual = "red")) +
    guides(fill  = guide_legend(reverse = TRUE),
           colour = guide_legend(order = 2)) +
    facet_wrap(~ series, ncol = 3, scales = "free_x") +
    theme(legend.position = "bottom",
          legend.title = element_blank(),
          axis.text.y = element_blank(),
          axis.ticks.y = element_blank())
)

# 9b. RMSE ratio per series (dashed line = base)
rdf <- data.frame(series = factor(rep(ord, each = 3), levels = ord, labels = panel_labels),
                  method = factor(rep(methods[-1], 9), levels = methods[-1]),
                  ratio  = as.numeric(ratio))
cols_ratio <- c(pal[methods[-1]], Base = "red")

print(
  ggplot(rdf, aes(x = method, y = ratio, fill = method)) +
    geom_col(width = 0.7) +
    geom_hline(aes(yintercept = 1, colour = "Base"),
               linetype = "dashed", linewidth = 0.9) +
    scale_fill_manual(values = cols_ratio, breaks = names(cols_ratio)) +
    scale_colour_manual(values = cols_ratio, breaks = names(cols_ratio)) +
    facet_wrap(~ series, ncol = 3) +
    labs(x = NULL, y = "RMSE ratio", fill = NULL, colour = NULL) +
    theme_minimal() +
    theme(legend.position = "bottom",
          legend.title = element_blank(),
          axis.text.x = element_blank(),
          axis.ticks.x = element_blank())
)

# 9c. Absolute RMSE per series (reconciled methods only)
rmse_df <- data.frame(series = factor(rep(ord, each = 3), levels = ord, labels = panel_labels),
                      method = factor(rep(methods[-1], 9), levels = methods[-1]),
                      rmse   = as.numeric(rmse_tbl[-1, ]))

print(
  ggplot(rmse_df, aes(x = method, y = rmse, fill = method)) +
    geom_col(width = 0.7) +
    scale_fill_manual(values = pal[methods[-1]], breaks = methods[-1]) +
    facet_wrap(~ series, ncol = 3, scales = "free_y") +
    labs(x = NULL, y = "RMSE", fill = NULL) +
    theme_minimal() +
    theme(legend.position = "bottom",
          legend.title = element_blank(),
          axis.text.x = element_blank(),
          axis.ticks.x = element_blank())
)

# 9d. 1 August 2023, hour by hour: actual vs base and reconciled forecasts
stopifnot(nrow(test_all) >= 24,
          format(test_df$datetime[1], "%Y-%m-%d %H") == "2023-08-01 00")

mk_day <- function(m, label) {
  data.frame(hour   = rep(0:23, times = 9),
             series = rep(ord, each = 24),
             method = label,
             value  = as.numeric(m[1:24, ord]))
}
day_df <- rbind(mk_day(test_all, "Actual"),
                mk_day(yhat,     "Base"),
                mk_day(bu_mat,   "Bottom-up"),
                mk_day(ls_mat,   "WLS"),
                mk_day(ml_mat,   "Machine learning"))
day_df$series <- factor(day_df$series, levels = ord, labels = panel_labels)
day_df$method <- factor(day_df$method, levels = c(methods, "Actual"))

print(
  ggplot(day_df, aes(x = hour, y = value, colour = method)) +
    geom_line(aes(linewidth = method == "Actual")) +
    facet_wrap(~ series, ncol = 3, scales = "free_y") +
    scale_linewidth_manual(values = c(`TRUE` = 1.2, `FALSE` = 0.7), guide = "none") +
    scale_colour_manual(values = c(pal, Actual = "black")) +
    scale_x_continuous(breaks = seq(0, 24, 6)) +
    labs(x = "Hour of day", y = "Count", colour = NULL) +
    guides(colour = guide_legend(reverse = TRUE)) +
    theme_minimal() +
    theme(legend.position = "bottom")
)

# 9e. Test week (August 2023) for three selected series
aug <- which(format(test_df$datetime, "%Y-%m") == "2023-08")
stopifnot(length(aug) > 0)

mk_week <- function(m, label) {
  data.frame(datetime = test_df$datetime[aug],
             series   = rep(ord, each = length(aug)),
             method   = label,
             value    = as.numeric(m[aug, ord]))
}
week_df <- rbind(mk_week(test_all, "Actual"),
                 mk_week(yhat,     "Base"),
                 mk_week(bu_mat,   "Bottom-up"),
                 mk_week(ls_mat,   "WLS"),
                 mk_week(ml_mat,   "Machine learning"))
week_df$series <- factor(week_df$series, levels = ord, labels = panel_labels)
week_df$method <- factor(week_df$method, levels = c(methods, "Actual"))

selected_series <- c("y1", "y7 (areas 1-4)", "y9 (all areas)")
week_sub <- subset(week_df, series %in% selected_series)
week_sub$series <- droplevels(week_sub$series)

print(
  ggplot(week_sub, aes(x = datetime, y = value, colour = method, linetype = method)) +
    geom_line(aes(linewidth = method == "Actual")) +
    facet_wrap(~ series, ncol = 3, scales = "free_y") +
    scale_linewidth_manual(values = c(`TRUE` = 0.4, `FALSE` = 0.8), guide = "none") +
    scale_colour_manual(values = c(pal, Actual = "red")) +
    scale_linetype_manual(values = c(Base = "solid", `Bottom-up` = "solid",
                                     WLS = "solid", `Machine learning` = "solid",
                                     Actual = "dashed")) +
    scale_x_datetime(date_breaks = "3 days", date_labels = "%d %b") +
    labs(x = NULL, y = "Count", colour = NULL, linetype = NULL) +
    theme_minimal() +
    theme(legend.position = "bottom")
)

                           
# Mean Absolute Error
mae <- function(f, a) colMeans(abs(f - a))
# MAE over the test week
mae_tbl <- rbind(
  Base = mae(yhat, test_all),
  BU   = mae(bu_mat, test_all),
  WLS  = mae(ls_mat, test_all),
  ML   = mae(ml_mat, test_all)
)[, ord]

# Ratio relative to the Base model
mae_ratio <- sweep(
  mae_tbl[-1, ],
  2,
  mae_tbl["Base", ],
  "/"
)

cat("\n=== MAE ===\n")
print(round(mae_tbl, 1))

cat("\n=== MAE ratio ===\n")
print(round(mae_ratio, 3))

                           
# ---- Total runtime ----
cat("\nTOTAL runtime:",
    round(as.numeric(Sys.time() - t_total, units = "secs"), 1), "s\n")
