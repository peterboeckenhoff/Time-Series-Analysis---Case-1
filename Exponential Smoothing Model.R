# CASE 1: FORECASTING CITI BIKE BICYCLE TRIPS

library(tidyr) # for reshaping the data - transform wide to long dataset
library(ggplot2) # for customizable graphs
library(dplyr)
library(forecast)


# DATA
# Import citibike dataset
load("C:/Users/thomas/OneDrive/E&OR Msc/Time Series Analysis/Case 1 Citi Bike/Case 1/citibike.RData")

# Number of observations
n_train <- nrow(train_df)
n_test  <- nrow(test_df)

# Combine training and test data chronologically
y_all <- c(
  train_df$`area 6`,
  test_df$`area 6`
)

# Create datetime
data_citibike$datetime <- as.POSIXct(
  paste(
    data_citibike$year,
    data_citibike$month,
    data_citibike$day,
    data_citibike$hour
  ),
  format = "%Y %m %d %H"
)

# Sort observations chronologically
data_citibike <- data_citibike[
  order(data_citibike$datetime),
]

# Training period: April 1 - July 31, 2023
train_df <- subset(
  data_citibike,
  datetime >= as.POSIXct("2023-04-01 00:00:00") &
    datetime < as.POSIXct("2023-08-01 00:00:00")
)

# Test period: observations after July 2023
test_df <- subset(
  data_citibike,
  datetime >= as.POSIXct("2023-08-01 00:00:00")
)

# Create time series for area 6
y_train <- ts(
  train_df$`area 6`,
  frequency = 24
)


# DATA EXPLORATION

# Create the dataset
data_long <- pivot_longer(
  data_citibike,
  cols = `area 1`:`area 6`,
  names_to = "area",
  values_to = "trips"
)



# Look at Area 6 training data and determine which ETS specification fits these data
# best according to the automatic selection procedure."
ets_auto <- ets(y_train)
summary(ets_auto)
# Automatic ETS model selected: ETS(M,A,M)

# residual diagnostics. We want to check whether the ETS model has captured 
# the systematic patterns in the training data.
checkresiduals(ets_auto)



# ROLLING WINDOW
# Training data: April-July

# Rolling window: 168 hours (7 days)
window_size <- 168
h <- 24

n_days <- n_test %/% h

# Store one 24-hour forecast per day
max_h <- 24
forecasts_ets <- matrix(
  NA_real_,
  nrow = n_days,
  ncol = h
)

colnames(forecasts_ets) <- paste0("h", 1:h)

for (day in 1:n_days) {
  
  # Move forward by 24 hours each day
  end_window <- n_train + (day - 1) * h
  
  start_window <- end_window - window_size + 1
  
  y_window <- ts(
    y_all[start_window:end_window],
    frequency = 24
  )
  
  fit <- ets(y_window, model = "MAM")
  
  fc <- forecast(fit, h = h)
  
  forecasts_ets[day, ] <- as.numeric(fc$mean)
}

# Actual observations, arranged into days
actuals_matrix <- matrix(
  test_df$`area 6`,
  ncol = h,
  byrow = TRUE
)

# Calculate accuracy separately for each horizon
ets_results <- data.frame(
  horizon = 1:h,
  RMSE = NA_real_,
  MAE = NA_real_,
  MAPE = NA_real_
)

for (hh in 1:h) {
  
  actual <- actuals_matrix[, hh]
  predicted <- forecasts_ets[, hh]
  
  errors <- actual - predicted
  
  ets_results$RMSE[hh] <- sqrt(mean(errors^2))
  ets_results$MAE[hh] <- mean(abs(errors))
  
  nonzero <- actual != 0
  ets_results$MAPE[hh] <- mean(
    abs(errors[nonzero] / actual[nonzero])
  ) * 100
}

ets_results



# EXPANDING WINDOW
# Expanding window ETS forecasts

# Expanding window: use all observations available so far
forecasts_ets_expanding <- matrix(
  NA_real_,
  nrow = n_days,
  ncol = h
)

colnames(forecasts_ets_expanding) <- paste0("h", 1:h)

for (day in 1:n_days) {
  
  # Move forward by 24 hours each day
  end_window <- n_train + (day - 1) * h
  
  y_window <- ts(
    y_all[1:end_window],
    frequency = 24
  )
  
  fit <- ets(y_window, model = "MAM")
  
  fc <- forecast(fit, h = h)
  
  forecasts_ets_expanding[day, ] <- as.numeric(fc$mean)
}

# Calculate accuracy separately for each horizon
ets_expanding_results <- data.frame(
  horizon = 1:h,
  RMSE = NA_real_,
  MAE = NA_real_,
  MAPE = NA_real_
)

for (hh in 1:h) {
  
  actual <- actuals_matrix[, hh]
  predicted <- forecasts_ets_expanding[, hh]
  
  errors <- actual - predicted
  
  ets_expanding_results$RMSE[hh] <- sqrt(mean(errors^2))
  ets_expanding_results$MAE[hh] <- mean(abs(errors))
  
  nonzero <- actual != 0
  ets_expanding_results$MAPE[hh] <- mean(
    abs(errors[nonzero] / actual[nonzero])
  ) * 100
}

ets_expanding_results

# Function to calculate pooled accuracy metrics
calculate_pooled_metrics <- function(forecasts, actuals) {
  
  errors <- actuals - forecasts
  
  errors_vec <- as.vector(errors)
  actuals_vec <- as.vector(actuals)
  
  # Remove missing or non-finite values
  valid <- is.finite(errors_vec) & is.finite(actuals_vec)
  
  errors_vec <- errors_vec[valid]
  actuals_vec <- actuals_vec[valid]
  
  # MAPE: exclude observations where actual demand is zero
  nonzero <- actuals_vec != 0
  
  c(
    RMSE = sqrt(mean(errors_vec^2)),
    MAE  = mean(abs(errors_vec)),
    MAPE = mean(
      abs(errors_vec[nonzero] / actuals_vec[nonzero])
    ) * 100
  )
}

# Pooled metrics: Rolling ETS
rolling_metrics <- calculate_pooled_metrics(
  forecasts_ets,
  actuals_matrix
)

# Pooled metrics: Expanding ETS
expanding_metrics <- calculate_pooled_metrics(
  forecasts_ets_expanding,
  actuals_matrix
)

# Display results in one table
pooled_ets_results <- data.frame(
  Window = c("Rolling ETS", "Expanding ETS"),
  RMSE = c(rolling_metrics["RMSE"], expanding_metrics["RMSE"]),
  MAE  = c(rolling_metrics["MAE"], expanding_metrics["MAE"]),
  MAPE = c(rolling_metrics["MAPE"], expanding_metrics["MAPE"])
)

pooled_ets_results




# Rolling-window ETS forecasts
rolling_ets_mcs <- data.frame(
  day = rep(1:n_days, each = h),
  horizon = rep(1:h, times = n_days),
  actual = as.vector(t(actuals_matrix)),
  ets_forecast = as.vector(t(forecasts_ets))
)

# Expanding-window ETS forecasts
expanding_ets_mcs <- data.frame(
  day = rep(1:n_days, each = h),
  horizon = rep(1:h, times = n_days),
  actual = as.vector(t(actuals_matrix)),
  ets_forecast = as.vector(t(forecasts_ets_expanding))
)


# ASSESSING SIGNIFICANCE OF FORECAST RESULTS
# Calculate forecast errors:
# Error = actual demand - forecast demand
# Positive error means underprediction; negative means overprediction.

# Both forecast matrices and actuals_matrix have:
# Rows    = daily forecast origins (31 days)
# Columns = forecast horizons (h = 1, ..., 24)

forecast_errors_ets <- actuals_matrix - forecasts_ets

forecast_errors_ets_expanding <-
  actuals_matrix - forecasts_ets_expanding

# Give the columns meaningful names
colnames(forecast_errors_ets) <- paste0("h", 1:max_h)
colnames(forecast_errors_ets_expanding) <- paste0("h", 1:max_h)

# The DM test compares the predictive accuracy of two models.
# H0: The models have equal predictive accuracy.
# H1: Their predictive accuracy differs.
#
# We test rolling ETS against expanding ETS separately
# for each forecast horizon.

dm_results <- data.frame(
  horizon = 1:max_h,
  DM = NA_real_,
  p_value = NA_real_
)

for (hh in 1:max_h) {
  
  # Select errors for the same horizon from both models
  errors_rolling <- forecast_errors_ets[, hh]
  errors_expanding <- forecast_errors_ets_expanding[, hh]
  
  # Retain only forecast origins with errors from both models
  valid <- complete.cases(errors_rolling, errors_expanding)
  
  errors_rolling <- errors_rolling[valid]
  errors_expanding <- errors_expanding[valid]
  
  # Run the two-sided DM test using squared forecast errors
  dm <- dm.test(
    errors_rolling,
    errors_expanding,
    alternative = "two.sided",
    h = 1,
    power = 2,
    varestimator = "bartlett"
  )
  
  # Store the test statistic and p-value
  dm_results$DM[hh] <- as.numeric(dm$statistic)
  dm_results$p_value[hh] <- dm$p.value
}

# Determine significance at the 5% level
dm_results$significant <- dm_results$p_value < 0.05


# Adjust the 24 DM-test p-values
dm_results$p_adjusted <- p.adjust(
  dm_results$p_value,
  method = "holm"
)

# Determine significance after correction
dm_results$significant_adjusted <-
  dm_results$p_adjusted < 0.05

dm_results

