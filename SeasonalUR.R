# Case 1 - Forecasting Citi Bike Bicycle Trips: unit root tests (area 6)

library(bootUR)
library(uroot)
library(forecast)
library(ggplot2)

# ---- 0. Settings ----
# Set the working directory to the folder containing the data file first:
# setwd("path/to/folder")
data_file <- "/Users/peterbockenhoff/Downloads/citibike.RData"

nb_daily  <- 999   # bootstrap replicates, daily HEGY
nb_weekly <- 100   # bootstrap replicates, weekly HEGY (kept low because of runtime)

set.seed(123)

# ---- 1. Data ----
load(data_file)

data_citibike$datetime <- as.POSIXct(
  paste(data_citibike$year, data_citibike$month, data_citibike$day, data_citibike$hour),
  format = "%Y %m %d %H"
)
data_citibike <- data_citibike[order(data_citibike$datetime), ]

area6_trips <- data_citibike[["area 6"]]

# ---- 2. Regular unit root test (ADF, Pantula principle) ----
# Test d = 2, 1, 0 with the deterministic terms none / intercept / trend
diff_d2 <- diff(area6_trips, differences = 2)
diff_d1 <- diff(area6_trips, differences = 1)

adf_d2 <- boot_adf(diff_d2,     deterministics = "none")        # d = 2
adf_d1 <- boot_adf(diff_d1,     deterministics = "intercept")   # d = 1
adf_d0 <- boot_adf(area6_trips, deterministics = "trend")       # d = 0

results_adf <- data.frame(
  d              = c(2, 1, 0),
  deterministics = c("None", "Intercept", "Trend"),
  statistic      = round(c(adf_d2$statistic, adf_d1$statistic, adf_d0$statistic), 3),
  p_value        = round(c(adf_d2$p.value,   adf_d1$p.value,   adf_d0$p.value),   3)
)
cat("\n=== ADF (Pantula principle) ===\n")
print(results_adf)

# ---- 3. Seasonal unit root test (HEGY) ----
trips_ts_daily  <- ts(area6_trips, frequency = 24)    # daily cycle
trips_ts_weekly <- ts(area6_trips, frequency = 168)   # weekly cycle

hegy_24_trend  <- hegy.test(trips_ts_daily,  deterministic = c(1, 1, 0))
hegy_24_const  <- hegy.test(trips_ts_daily,  deterministic = c(1, 0, 0))
hegy_168_trend <- hegy.test(trips_ts_weekly, deterministic = c(1, 1, 0))
hegy_168_const <- hegy.test(trips_ts_weekly, deterministic = c(1, 0, 0))

cat("\n=== HEGY daily (constant + trend) ===\n");    print(hegy_24_trend)
cat("\n=== HEGY daily (constant only) ===\n");       print(hegy_24_const)
cat("\n=== HEGY weekly (constant + trend) ===\n");   print(hegy_168_trend)
cat("\n=== HEGY weekly (constant only) ===\n");      print(hegy_168_const)

# Canova-Hansen test (complementary check, not used)
# ch_24  <- ch.test(trips_ts_daily)
# ch_168 <- ch.test(trips_ts_weekly)

# ---- 4. Bootstrapped HEGY p-values (constant + trend) ----
# Daily cycle (s = 24)
hegy_24_boot <- hegy.boot.pval(
  x             = trips_ts_daily,
  model0        = hegy_24_trend$fit,
  stats0        = hegy_24_trend$stat,
  deterministic = c(1, 1, 0),
  lag.method    = "fixed",
  maxlag        = 0,
  byseason      = FALSE,
  nb            = nb_daily
)
cat("\n=== Bootstrapped HEGY daily ===\n")
print(hegy_24_boot)

# Weekly cycle (s = 168)
hegy_168_boot <- hegy.boot.pval(
  x             = trips_ts_weekly,
  model0        = hegy_168_trend$fit,
  stats0        = hegy_168_trend$stat,
  deterministic = c(1, 1, 0),
  lag.method    = "fixed",
  maxlag        = 0,
  byseason      = FALSE,
  nb            = nb_weekly
)
cat("\n=== Bootstrapped HEGY weekly ===\n")
print(hegy_168_boot)

# ---- 5. MSTL decomposition (daily + weekly seasonality) ----
trips_msts <- msts(area6_trips, seasonal.periods = c(24, 168))
mstl_fit   <- mstl(trips_msts)

print(
  autoplot(mstl_fit) +
    labs(title = "MSTL Decomposition of Area 6 Trips")
)