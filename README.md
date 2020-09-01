# tshacks

This package creates Time Series Hacks to enhance time series problems, such as forecasting, scenario generation among others. List of methods:

- [Bagged.Cluster.ETS](Bagged.md)
- TBA



First, I haven't uploaded a stable version to CRAN but you can download from GitHub using devtools, type:

```r
library(devtools)
install_github("tiagomendesdantas/tshacks")
```

### Examples:

#### Generate forecast using Bagged.Cluster.ETS
```r
library(tshacks)
method <- baggedClusterETS(gas)
f_method <- forecast(method,h=12)
plot(f_method)
```

<center><img src="figures/forecast.png" alt="Forecast using BaggedClusterETS" style="width: 750px;"/></center>


#### Details:

To access the bootstrap versions before the cluster phase:
```r
bootstrap_before_selection <- method$bootstrapped_series_ori
```


To access the number o clusters selected:
```r
n_clusters <- method$k
```

### To access the bootstrap versions selected after the cluster phase
```r
bootstrap_after_selection <- method$bootstrapped_series
```

### To access the selected model for each selected bootstrapped version  
```r
selected_models <- method$models

```

### To access the fitted series and the residuals
```r
fitted_series <- method$fitted
residuals <- method$residuals
```
