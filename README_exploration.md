# Citi Bike – Data Exploration

`<DataExploration>.R` creates the exploratory plots used in the presentation
(demand over time, by hour, by weekday, by month and by area). It does not fit any models.

1. Put `<DataExploration>.R` and `citibike.RData` in one folder.
2. Install the packages:

```r
install.packages(c("tidyr", "ggplot2"))
```

3. In RStudio, open the script, set the working directory to that folder
   (*Session > Set Working Directory > To Source File Location*) and run everything
   (*Ctrl/Cmd + Alt + R*) or update the lines
   ```r
   "data_file  <- "citibike.RData"
   "recon_file <- "ReconData.RData"
   ```
   manually.

The plots appear one after another in the *Plots* pane.

`<ExampleSeasonalUR>.R` runs the UR tests for an example time series and generates the plots used in the presentation.
