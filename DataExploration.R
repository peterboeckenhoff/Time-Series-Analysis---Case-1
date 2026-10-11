# Case 1 - Forecasting Citi Bike Bicycle Trips: data exploration

library(tidyr)     # reshape wide -> long
library(ggplot2)   # plots

# ---- 0. Settings ----
# Set the working directory to the folder containing the data file first:
# setwd("path/to/folder")
data_file <- "citibike.RData"

# ---- 1. Data ----
load(data_file)

data_citibike$datetime <- as.POSIXct(
  paste(data_citibike$year, data_citibike$month, data_citibike$day, data_citibike$hour),
  format = "%Y %m %d %H"
)
data_citibike <- data_citibike[order(data_citibike$datetime), ]
head(data_citibike)

# Long format: one row per hour and area
data_long <- pivot_longer(
  data_citibike,
  cols = `area 1`:`area 6`,
  names_to = "area",
  values_to = "trips"
)
data_long$date         <- as.Date(data_long$datetime)
data_long$weekday      <- as.numeric(format(data_long$datetime, "%u"))  # 1 = Monday, 7 = Sunday
data_long$hour         <- as.numeric(format(data_long$datetime, "%H"))
data_long$hour_of_week <- (data_long$weekday - 1) * 24 + data_long$hour

# Total trips per hour across all areas
data_total <- aggregate(trips ~ datetime, data = data_long, FUN = sum)
data_total$date         <- as.Date(data_total$datetime)
data_total$weekday      <- as.numeric(format(data_total$datetime, "%u"))
data_total$hour         <- as.numeric(format(data_total$datetime, "%H"))
data_total$hour_of_week <- (data_total$weekday - 1) * 24 + data_total$hour
data_total$month        <- as.numeric(format(data_total$datetime, "%m"))
data_total$day_type     <- ifelse(data_total$weekday >= 6, "Weekend", "Weekday")

week_labels <- c("Mon", "Tue", "Wed", "Thu", "Fri", "Sat", "Sun")

# ---- 2. Total demand over time (grey) with 7-day rolling mean (blue) ----
data_total$roll_mean <- as.numeric(stats::filter(data_total$trips, rep(1 / (24 * 7), 24 * 7), sides = 2))

print(
  ggplot(data_total, aes(x = datetime)) +
    geom_line(aes(y = trips), colour = "grey70") +
    geom_line(aes(y = roll_mean), colour = "blue", na.rm = TRUE) +
    labs(x = "Time", y = "Trips")
)

# ---- 3. Average daily demand by day of the week (all areas) ----
daily <- aggregate(trips ~ date, data = data_total, FUN = sum)
daily$weekday <- factor(weekdays(daily$date),
                        levels = c("Monday", "Tuesday", "Wednesday", "Thursday",
                                   "Friday", "Saturday", "Sunday"))
weekly <- aggregate(trips ~ weekday, data = daily, FUN = mean)

print(
  ggplot(weekly, aes(x = weekday, y = trips, group = 1)) +
    geom_line() +
    geom_point() +
    labs(title = "Average daily demand by day of week", x = "", y = "Average trips per day")
)

# ---- 4. Average hourly demand over a day (all areas) ----
hourly <- aggregate(trips ~ hour, data = data_total, FUN = mean)

print(
  ggplot(hourly, aes(x = hour, y = trips)) +
    geom_line() +
    geom_point() +
    labs(title = "Average hourly demand over a day", x = "Hour of day", y = "Average trips per hour")
)

# ---- 5. Average demand per hour of the week, for each of the 6 areas ----
week_hour <- aggregate(trips ~ hour_of_week + area, data = data_long, FUN = mean)

print(
  ggplot(week_hour, aes(x = hour_of_week, y = trips)) +
    geom_line() +
    facet_wrap(~ area, ncol = 2, scales = "free_y") +
    scale_x_continuous(breaks = seq(0, 144, 24), labels = week_labels) +
    labs(x = "", y = "Average trips per hour")
)

# ---- 6. Weekday vs weekend hourly demand (all areas) ----
hourly_daytype <- aggregate(trips ~ hour + day_type, data = data_total, FUN = mean)

print(
  ggplot(hourly_daytype, aes(x = hour, y = trips)) +
    geom_line() +
    geom_point() +
    facet_wrap(~ day_type) +
    scale_x_continuous(breaks = seq(0, 23, 3)) +
    labs(x = "Hour", y = "Average trips per hour")
)

# ---- 7. Average demand per month (all areas) ----
monthly <- aggregate(trips ~ month, data = data_total, FUN = mean)

print(
  ggplot(monthly, aes(x = month, y = trips)) +
    geom_line() +
    geom_point() +
    scale_x_continuous(breaks = 1:12, labels = month.abb) +
    labs(x = "", y = "Average trips per hour")
)

# ---- 8. Average demand per hour of the week (all areas) ----
week_hour_total <- aggregate(trips ~ hour_of_week, data = data_total, FUN = mean)

print(
  ggplot(week_hour_total, aes(x = hour_of_week, y = trips)) +
    geom_line() +
    scale_x_continuous(breaks = seq(0, 144, 24), labels = week_labels) +
    labs(x = "", y = "Average trips per hour")
)

# ---- 9. Average demand per day of the week, for each of the 6 areas ----
daily_area <- aggregate(trips ~ date + area, data = data_long, FUN = sum)
daily_area$weekday <- as.numeric(format(daily_area$date, "%u"))

week_day <- aggregate(trips ~ weekday + area, data = daily_area, FUN = mean)

print(
  ggplot(week_day, aes(x = weekday, y = trips)) +
    geom_line() +
    geom_point() +
    facet_wrap(~ area, ncol = 2, scales = "free_y") +
    scale_x_continuous(breaks = 1:7, labels = week_labels) +
    labs(x = "", y = "Average trips per day")
)

# ---- 10. Average demand per day of the week for all 9 series ----
# y1-y6 = areas, y7 = areas 1-4, y8 = areas 5-6, y9 = all areas
daily_bottom <- transform(daily_area, series = sub("area ", "y", area))[, c("date", "series", "trips")]

agg_daily <- function(areas, label) {
  d <- aggregate(trips ~ date, data = subset(daily_area, area %in% paste("area", areas)), FUN = sum)
  d$series <- label
  d[, c("date", "series", "trips")]
}

labels_all <- c(paste0("y", 1:6), "y7 (areas 1-4)", "y8 (areas 5-6)", "y9 (all areas)")

daily_all <- rbind(
  daily_bottom,
  agg_daily(1:4, "y7 (areas 1-4)"),
  agg_daily(5:6, "y8 (areas 5-6)"),
  agg_daily(1:6, "y9 (all areas)")
)
daily_all$series  <- factor(daily_all$series, levels = labels_all)   # keeps y1 ... y9 order
daily_all$weekday <- as.numeric(format(daily_all$date, "%u"))

week_day_all <- aggregate(trips ~ weekday + series, data = daily_all, FUN = mean)

print(
  ggplot(week_day_all, aes(x = weekday, y = trips)) +
    geom_line() +
    geom_point() +
    facet_wrap(~ series, ncol = 3, scales = "free_y") +
    scale_x_continuous(breaks = 1:7, labels = week_labels) +
    labs(x = "", y = "Average trips per day")
)