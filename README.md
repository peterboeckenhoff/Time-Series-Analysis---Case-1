# Time Series Analysis – Case 1

This repository is organized around a single case study in time series analysis, with each major topic kept in its own branch. The goal is to keep the workflow readable, modular, and easy to navigate without mixing different analytical methods in the same branch.

## Repository structure

The project is structured as follows:

- `main` — landing page and repository overview; this is the branch you are reading now.
- `Data-Exploration` — initial data inspection, cleaning, visual exploration, and understanding of the time series behavior.
- `Exponential-Smoothing-Model` — implementation and evaluation of exponential smoothing approaches.
- `ARIMA-and-SARIMA` — ARIMA and Seasonal ARIMA modeling, diagnostics, and interpretation.
- `Forecast-Reconciliation` — forecasting reconciliation methods and aggregation of forecast outputs.
- `Seasonal-Unit-Root-test` — tests for seasonal stationarity and unit-root-related diagnostics.

## Why this structure?

Each branch represents one analytical topic of the case study. This approach helps us:

- keep methods isolated and easy to compare,
- avoid mixing exploratory, modeling, and validation work in one place,
- maintain a clear progression from data understanding to forecasting,
- make it easier for collaborators to review and discuss one topic at a time.

## Recommended workflow

1. Start in `main` to get the overall context of the project.
2. Switch to the topic branch relevant to the analysis you want to review or continue.
3. Work on that topic independently before moving to the next modeling step.
4. Use the repository overview in `main` as the map for the full case study.

## Typical case-study flow

A natural progression is:

- `Data-Exploration` → understand the data and identify patterns
- `Seasonal-Unit-Root-test` → assess stationarity and seasonal dynamics
- `ARIMA-and-SARIMA` → model autocorrelation and seasonality
- `Exponential-Smoothing-Model` → test baseline forecasting approaches
- `Forecast-Reconciliation` → combine or reconcile forecast outputs

## Notes

This branch-based layout is intentional: the repository is not organized as a single monolithic project, but as a collection of focused analytical tracks. Each branch documents a distinct part of the same case study while preserving a clean and structured project history.

If you are new to the project, begin with `main` and then navigate to the relevant topic branch for the section you want to study.
