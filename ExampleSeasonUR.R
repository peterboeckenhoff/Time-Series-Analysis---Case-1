set.seed(555)

n          <- 24 * 2  # hourly data over four days: 96 observations
s          <- 12      # seasonal period
shock_t    <- 20      # time of the one-off shock
shock_size <- 5       # size of the one-off shock

# ---- DGP: seasonal AR(1), y_t = rho * y_{t-4} + eps_t ----
simulate_seasonal <- function(rho, n, s, shock_t, shock_size, sd = 0.3) {
  y <- numeric(n)
  y[1:s] <- c(2, -1, 1, -2)
  for (t in (s + 1):n) {
    eps <- rnorm(1, 0, sd)
    if (t == shock_t) eps <- eps + shock_size
    y[t] <- rho * y[t - s] + eps
  }
  y
}

y1 <- simulate_seasonal(rho = 1.0, n, s, shock_t, shock_size)  # seasonal unit root
y2 <- simulate_seasonal(rho = 0.5, n, s, shock_t, shock_size)  # stationary seasonality

# Shared y-range so the two plots are directly comparable
ylim <- range(c(y1, y2))

# ---- Plot helper: raw series + seasonal moving average + shock marker ----
# The moving average is what makes "does it revert?" visible despite the
# seasonal zig-zag; the dotted line is the pre-shock average level.
plot_seasonal <- function(y, main, shock_t, s, ylim) {
  plot(y, type = "o", pch = 16,
       xlab = "Time", ylab = "y[t]",
       main = main, ylim = ylim)
  
  # Seasonal moving average spanning the full series:
  # trailing MA for t >= s, running mean for the warm-up t = 1, ..., s - 1
  ma <- stats::filter(y, rep(1 / s, s), sides = 1)
  ma[1:(s - 1)] <- cumsum(y[1:(s - 1)]) / (1:(s - 1))
  lines(ma, col = "red", lwd = 2)
  
  abline(h = mean(y[1:(shock_t - 1)]), lty = 3, col = "gray40")
  abline(v = shock_t, lty = 2, col = "steelblue")
  text(x = shock_t, y = ylim[2], "Shock", pos = 4-1, col = "steelblue", xpd = NA)
}
# ---- Plot ----
plot_seasonal(y1, "Seasonal Unit Root (rho = 1, s = 12)", shock_t, s, ylim)
plot_seasonal(y2, "Stationary Seasonality (rho = 0.5, s = 12)", shock_t, s, ylim)

library(bootUR)

# ---- Pantula sequence: test d = dmax, dmax-1, ..., 0, stopping logic applied after ----
pantula_adf <- function(y, dmax = 2, level_deterministics = "intercept") {
  results <- data.frame(d = integer(0), deterministics = character(0),
                        statistic = numeric(0), p_value = numeric(0))
  
  for (d in dmax:0) {
    if (d > 0) {
      series   <- diff(y, differences = d)
      det_spec <- "none"
    } else {
      series   <- y
      det_spec <- level_deterministics
    }
    test <- boot_adf(series, deterministics = det_spec)
    results <- rbind(results, data.frame(
      d = d, deterministics = det_spec,
      statistic = test$statistic, p_value = round(test$p.value, 3)
    ))
  }
  results
}

# ---- Apply the Pantula principle: reject moving down from d = dmax; stop at
#      the first FAILURE to reject -- that tells you how many differences you need ----
decide_d <- function(results, alpha = 0.05) {
  ord <- results[order(-results$d), ]
  for (i in seq_len(nrow(ord))) {
    if (ord$p_value[i] > alpha) {                 # fail to reject: unit root still there
      if (ord$d[i] == max(ord$d)) {
        return(sprintf("d >= %d (raise d_max, inconclusive)", ord$d[i]))
      }
      return(ord$d[i] + 1)                        # one more difference than this step needed
    }
  }
  return(0)                                        # rejected at every step, incl. the level
}

# ---- Run on both series ----
results_y1 <- pantula_adf(y1)
results_y2 <- pantula_adf(y2)

cat("Pantula sequence -- y1 (seasonal unit root, rho = 1)\n")
print(results_y1)
cat("Pantula-selected d for y1:", decide_d(results_y1), "\n\n")

cat("Pantula sequence -- y2 (stationary seasonality, rho = 0.5)\n")
print(results_y2)
cat("Pantula-selected d for y2:", decide_d(results_y2), "\n")

# We are using the Boot-ADF (Augmented Dickey-Fuller) Unit Root Test introduced in the lecture

# ============================================================
# Pantula sequence for y1 (seasonal unit root, rho = 1)
# ============================================================
cat("Pantula sequence -- y1 (seasonal unit root, rho = 1)\n")

# d_max = 2
diff2_y1 <- diff(y1, differences = 2)
adf_d2_y1 <- boot_adf(diff2_y1, deterministics = "none")
print(adf_d2_y1)

# d = 1
diff1_y1 <- diff(y1, differences = 1)
adf_d1_y1 <- boot_adf(diff1_y1, deterministics = "none")
print(adf_d1_y1)

# d = 0 with intercept
adf_d0_y1 <- boot_adf(y1, deterministics = "intercept")
print(adf_d0_y1)

# Create a summary table with d, deterministic components, test statistic, and p-value
results <- data.frame(
  d = c(2, 1, 0),
  deterministics = c("none", "constant", "trend"),
  statistic = c(adf_d2_y1$statistic, adf_d1_y1$statistic, adf_d0_y1$statistic),
  p_value = c(adf_d2_y1$p.value, adf_d1_y1$p.value, adf_d0_y1$p.value)
)
# Round the p_value to 3 decimal places
results$p_value <- round(results$p_value, 3)

print(results)

# ============================================================
# Pantula sequence for y2 (stationary seasonality, rho = 0.5)
# ============================================================
cat("Pantula sequence -- y2 (stationary seasonality, rho = 0.5)\n")

# d_max = 2
diff2_y2 <- diff(y2, differences = 2)
adf_d2_y2 <- boot_adf(diff2_y2, deterministics = "none")
print(adf_d2_y2)

# d = 1
diff1_y2 <- diff(y2, differences = 1)
adf_d1_y2 <- boot_adf(diff1_y2, deterministics = "none")
print(adf_d1_y2)

# d = 0 with intercept
adf_d0_y2 <- boot_adf(y2, deterministics = "intercept")
print(adf_d0_y2)

# Create a summary table with d, deterministic components, test statistic, and p-value
results <- data.frame(
  d = c(2, 1, 0),
  deterministics = c("none", "constant", "trend"),
  statistic = c(adf_d2_y2$statistic, adf_d1_y2$statistic, adf_d0_y2$statistic),
  p_value = c(adf_d2_y2$p.value, adf_d1_y2$p.value, adf_d0_y2$p.value)
)
# Round the p_value to 3 decimal places
results$p_value <- round(results$p_value, 3)

print(results)