# Citi Bike Forecast Reconciliation

1. Put `citibike_reconciliation.R`, `citibike.RData` and `ReconData.RData` in one folder.
2. Install the packages:

```r
install.packages(c("forecast", "FoReco", "FoRecoML", "ggplot2",
                   "randomForest", "scales"))
```

3. In RStudio, open the script, set the working directory to that folder
   (*Session > Set Working Directory > To Source File Location*) and run everything
   (*Ctrl/Cmd + Alt + R*).

By default the saved SARIMA results in `ReconData.RData` are used, so the slow model
fitting is skipped. To refit the models yourself, set `run_sarima <- TRUE` at the top
of the script (takes 5 minutes on my laptop).

Note: To reduce complexity and running time, the same SARIMA(2,0,3)(2,0,2)[168] is
used for all 9 series.