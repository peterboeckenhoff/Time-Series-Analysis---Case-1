# Citi Bike – Unit Root Tests (Area 6)

`SeasonalUR.R` runs the regular (ADF) and seasonal (HEGY) unit root tests for area 6,
including bootstrapped p-values, and plots an MSTL decomposition.

1. Put `SeasonalUR.R` and `citibike.RData` in one folder.
2. Install the packages:

```r
install.packages(c("bootUR", "uroot", "forecast", "ggplot2"))
```

3. In RStudio, open the script, set the working directory to that folder
   (*Session > Set Working Directory > To Source File Location*) and run everything
   (*Ctrl/Cmd + Alt + R*) or update the lines
   ```r
   "data_file  <- "citibike.RData"
   "recon_file <- "ReconData.RData"
   ```
   manually.

The bootstrapped HEGY tests are slow. To speed them up, reduce
`nb_daily` and `nb_weekly` at the top of the script.
