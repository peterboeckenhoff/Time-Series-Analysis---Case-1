# README: Citi Bike forecasting case (ARIMA / SARIMA, DM tests, MCS)

This folder contains one R script, `Complete_code_combined.R`. It forecasts hourly Citi Bike trips for area 6 with ARIMA and SARIMA models (seasonal period s = 168), evaluates them out of sample, and tests the forecasts against each other and against an ETS benchmark with Diebold-Mariano (DM) tests and the Model Confidence Set (MCS).

## 1. What you need

**Software:** R with the packages `forecast` and `MCS`. The script installs `MCS` automatically if it is missing. `forecast` has to be installed beforehand (`install.packages("forecast")`).

**Files in the working directory:**

| File | Used for |
|---|---|
| `citibike.RData` | Data object `data_citibike` with columns `year`, `month`, `day`, `hour` and `area 6` |
| `rolling_ets_mcs.csv` | ETS forecasts, rolling window |
| `expanding_ets_mcs.csv` | ETS forecasts, expanding window |

Both ETS files need the columns `day`, `horizon`, `actual` and `ets_forecast`. The ETS forecasts are not computed in this script. They are read from these files, and the `actual` column must match the actual values used in the script (the script stops with an error if it does not).

## 2. How to run

1. Open `Complete_code_combined.R`.
2. Change the path in the `setwd(...)` line near the top (section 0b) to the folder that contains the three input files.
3. Run the whole script from top to bottom (`source("Complete_code_combined.R")` or Ctrl+Shift+Enter in RStudio).

The SARIMA models with s = 168 are slow to estimate. The grid search results are cached in `grid_cache_*.rds` files in the working directory, so a rerun reuses them. Delete these files to force a fresh grid search.

For a quick test, set `TEST_MODE <- TRUE` in section 7. This only forecasts 3 days. The ETS part will then stop with a length-mismatch error, because the ETS files cover all 31 days. Set it back to `FALSE` for the full run.

## 3. What the script does

| Part | Content |
|---|---|
| Settings and data | Training sample: 1 April to 31 July 2023. Test sample: August 2023. Forecast horizon h = 24 hours. |
| A | Identification (ACF/PACF plots) and estimation on the training sample. ARIMA(p,0,q) with p, q in 0:3 chosen by AIC. SARIMA with D = 0 and D = 1, P, Q in 0:3 chosen by a CSS-based criterion. Automatic `auto.arima` versions of both. |
| B | Out-of-sample exercise: 31 daily forecast origins in August, 24-hour-ahead forecasts. Done with a rolling window (fixed length) and an expanding window. Accuracy is reported as MSE, RMSE, MAE, MAPE and MASE, pooled and per horizon. |
| C | Pairwise DM tests between the five models (squared and absolute loss, Holm-adjusted p-values) and the MCS (alpha = 0.10, Tmax statistic, 5000 bootstrap replications). A final conclusion is printed for each window scheme. |
| C (continued) | The same DM tests and MCS, now comparing the models against the ETS benchmark. |
| D | Saves all results (see below). |

The five models are `arima_manual`, `arima_auto`, `sarima_manual_D0`, `sarima_manual_D1` and `sarima_auto`. The specifications that were selected are printed in the console and saved in `model_specifications.csv`.

## 4. Output

Everything that is saved ends up in the folder `results/`, inside the working directory.

| File | Content |
|---|---|
| `all_results.rds` | All main objects in one R list (fitted models, forecasts, accuracy, DM, MCS, ETS tests) |
| `all_forecast_results.rds` | Forecasts and actual values for both window schemes |
| `accuracy_rolling.csv`, `accuracy_expanding.csv` | Accuracy per model, pooled over all horizons |
| `accuracy_by_horizon_rolling.csv`, `accuracy_by_horizon_expanding.csv` | Accuracy per model and forecast horizon |
| `DM_<window>_<loss>.csv` | Pairwise DM tests between the ARIMA/SARIMA models (e.g. `DM_rolling_squared.csv`) |
| `MCS_membership.csv` | Models in the MCS, ARIMA/SARIMA models only |
| `DM_ETS_<window>_<loss>.csv` | DM tests of each model against ETS |
| `MCS_ETS_<window>_<loss>.csv`, `MCS_ETS_<window>_both_losses.rds` | MCS including ETS |
| `MCS_ETS_membership.csv` | Models in the MCS, including ETS |
| `model_specifications.csv` | Selected order of each model |
| `results_summary.txt` | Short summary of the main conclusions |

Plots (ACF/PACF, residual checks, RMSE per horizon) are shown in R and are not saved to files.

## 5. Things to be aware of

- **Two MCS runs use different settings.** The MCS on the ARIMA/SARIMA models alone uses block length k = 24 and seed 123 and removes models with identical forecasts first. The MCS including ETS uses the package default block length, seed 94991, and does not remove duplicates. Results of the two are therefore not strictly like-for-like.
- **Estimation method in the forecast loops.** The rolling and expanding window models are refitted with `method = "CSS"` for speed. A slower CSS-ML check can be switched on with `RUN_SANITY_CHECK <- TRUE`.
- **Models with identical forecasts.** If two models turn out to give identical forecasts, the pairwise DM tests and the first MCS drop the duplicate, and the console reports this.
- **Reproducibility.** The MCS bootstrap uses fixed seeds. The grid search caches depend on the selected p and q, so changing the training sample or the ranges requires deleting the cache files.
