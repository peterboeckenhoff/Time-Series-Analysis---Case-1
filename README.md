# Citi Bike – ETS Forecasting

`<ETS_Forecasting>.R` performs ETS forecasting for Citi Bike Area 6. It compares rolling and expanding window approaches using 24-hour forecasts, evaluates forecast accuracy, and tests whether the differences between the models are statistically significant.

1. Put `<ETS_Forecasting>.R` and `citibike.RData` in the appropriate folders. Make sure the file path in the script points to the dataset.

2. Install the required packages:

```r
install.packages(c("tidyr", "ggplot2", "dplyr", "forecast"))
```

3. In RStudio, open the script and ensure that the dataset is loaded and the training and test datasets are defined. Run the script from top to bottom.

4. The script fits an ETS(M,A,M) model and performs residual diagnostics. It then generates 24-hour forecasts using:
   - **Rolling window:** Uses the most recent 168 hours of observations.
   - **Expanding window:** Uses all observations available up to each forecast origin.

5. The script calculates RMSE, MAE and MAPE for each forecast horizon and across all forecasts. It also performs horizon-specific Diebold–Mariano tests and applies the Holm correction for multiple testing.

The resulting tables contain the forecast accuracy metrics, DM statistics, p-values and Holm-adjusted significance results.
