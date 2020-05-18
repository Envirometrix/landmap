# landmap package for R

----

- TOC
{:toc}

----

Package provides methodology for automated mapping i.e. spatial interpolation and/or 
prediction using **Ensemble Machine Learning** (extends functionality of the [mlr package](https://mlr.mlr-org.com/)). Key functionality includes:

* `train.spLearner` --- train a spatial prediction and/or interpolation model using Ensemble Machine Learning (works with numeric, binomial and factor-type variables),
* `buffer.dist` --- derive buffer (geographical) distances that can be used as covariates in spLearner, 
* `spc` --- derive Principal Components using stack of spatial layers,
* `tile` --- tile spatial layers so they can be used to run processing in parallel,
* `spsample.prob` --- determine inclusion probability / representation of a given point sample based on feature space analysis (maxlike function) and kernel density analysis,
* `download.landgis` --- access and download LandGIS layers from www.openlandmap.org,

Warning: most of functions are optimized to run in parallel by default. This might result in high RAM and CPU usage.

Spatial prediction using [Ensemble Machine Learning](https://koalaverse.github.io/machine-learning-in-R/stacking.html#stacking-software-in-r) with geographical distances 
is explained in detail in:

- Hengl, T., MacMillan, R.A., (2019). 
   [Predictive Soil Mapping with R](https://soilmapper.org/soilmapping-using-mla.html). 
   OpenGeoHub foundation, Wageningen, the Netherlands, 370 pages, www.soilmapper.org, 
   ISBN: 978-0-359-30635-0.
- Hengl, T., Nussbaum, M., Wright, M. N., Heuvelink, G. B., and Gräler, B. (2018). 
   [Random Forest as a generic framework for predictive modeling of spatial and spatio-temporal variables](https://doi.org/10.7717/peerj.5518). PeerJ 6:e5518.

## Installing

Install development versions from github:

```r
library(devtools)
install_github("envirometrix/landmap")
```

Under construction. Use for testing purposes only.

## Functionality

### Automated mapping using Ensemble Machine Learning

The following examples demostrates spatial prediction using the meuse data set:

```r
library(rgdal)
library(geoR)
library(plotKML)
library(raster)
library(glmnet)
library(xgboost)
library(kernlab)
library(deepnet)
library(Cubist)
demo(meuse, echo=FALSE)
m <- train.spLearner(meuse["lead"], covariates=meuse.grid[,c("dist","ffreq")], lambda = 1)
```

this runs several steps:

```
Converting ffreq to indicators...
Converting covariates to principal components...
Deriving oblique coordinates...TRUE
Fitting a variogram using 'linkfit' and trend model...TRUE
Estimating block size ID for spatial Cross Validation...TRUE
Starting parallelization in mode=socket with cpus=8.
Using learners: regr.ranger, regr.ksvm, regr.glmnet, regr.cubist...TRUE
Fitting a spatial learner using 'mlr::makeRegrTask'...TRUE
Exporting objects to slaves for mode socket: .mlr.slave.options
Mapping in parallel: mode = socket; cpus = 8; elements = 5.
Exporting objects to slaves for mode socket: .mlr.slave.options
Mapping in parallel: mode = socket; cpus = 8; elements = 5.
Exporting objects to slaves for mode socket: .mlr.slave.options
Mapping in parallel: mode = socket; cpus = 8; elements = 5.
Exporting objects to slaves for mode socket: .mlr.slave.options
Mapping in parallel: mode = socket; cpus = 8; elements = 5.
Stopped parallelization. All cleaned up.
```

The variogram model is only fitted to estimate effective range of spatial dependence.
Spatial Prediction models are based only on fitting the [Ensemble Machine Learning](https://koalaverse.github.io/machine-learning-in-R/stacking.html#stacking-software-in-r) 
(by default landmap uses `c("regr.ranger", "regr.ksvm", "regr.glmnet", "regr.cubist")`; see [a complete list of learners available via mlr](https://mlr.mlr-org.com/articles/tutorial/integrated_learners.html)) 
with oblique coordinates (rotated coordinates) as described in [Moller et al. (2019) 
"Oblique Coordinates as Covariates for Digital Soil Mapping"](https://www.soil-discuss.net/soil-2019-83/) to account for spatial autocorrelation in 
values. Geographical distances to ALL points can be added 
by specifying `buffer.dist=TRUE`; this is however not recommended for large point data sets.
The meta-learning i.e. the SuperLearner model shows which individual learners are most important:

```r
summary(m@spModel$learner.model$super.model$learner.model)
```
```
Call:
stats::lm(formula = f, data = d)

Residuals:
     Min       1Q   Median       3Q      Max 
-173.960  -35.549   -6.493   16.507  308.414 

Coefficients:
             Estimate Std. Error t value Pr(>|t|)    
(Intercept) -15.65153   13.62986  -1.148    0.253    
regr.ranger   0.91523    0.21575   4.242 3.86e-05 ***
regr.ksvm     0.33415    0.24691   1.353    0.178    
regr.glmnet   0.01482    0.15172   0.098    0.922    
regr.cubist  -0.13882    0.13639  -1.018    0.310    
---
Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1

Residual standard error: 72.28 on 150 degrees of freedom
Multiple R-squared:  0.5894,	Adjusted R-squared:  0.5784 
F-statistic: 53.82 on 4 and 150 DF,  p-value: < 2.2e-16
```

in this case `regr.ranger` seems to be most important for predicting lead concentration, 
while `regr.glmnet` is the least important. Overall this ensemble model explains ca 58% of variance (based on repeated cross-validation):

Next we can generate predictions using:

```r
meuse.lead <- predict(m)
```

Note that, based on the current set-up with `method = "stack.cv"`, every time you re-run the model training you 
might get somewhat different models / different betas. On the other hand, the final ensemble predictions (map) should visually not differ too much.

<center><img src="https://github.com/thengl/GeoMLA/blob/master/RF_vs_kriging/results/meuse/Fig_meuse_EML.png" width="650"></center> 
<center>_Figure: Predicted lead content for the Meuse data set. Model error is derived as weighted standard deviation from multiple model predictions._</center>


Animated predictions by 9 models (3x independently fitted random forest, SVM and Xgboost) looks like this 
(the coefficients are beta coefficients from the metalearner fit: the higher the coefficient, more important the model for the ensemble merge):

<center><img src="https://github.com/thengl/GeoMLA/blob/master/RF_vs_kriging/results/meuse/meuse_lead_ensemble.gif" width="400" /></center>
<center>_Figure: examples of independently generated predictions._</center>


The predictions shown in the image above incorporate spatial correlation between values, 
and hence can be used as a possible replacement for kriging methods ([Hengl et al. 2018](https://doi.org/10.7717/peerj.5518)). Automation comes, however, at the high computing and RAM usage costs.

In the following example we use somewhat larger data set from the SIC1997 exercise.

```r
data("sic1997")
X <- sic1997$swiss1km[c("CHELSA_rainfall","DEM")]
mR <- train.spLearner(sic1997$daily.rainfall, covariates=X, lambda=1)
rainfall1km <- predict(mR)
```

The processing is now much more computational because the data set consists from 467 points (hence 467 buffer distance maps need to be produced).
This will make the regression matrix becoming extensive, and also 5x3 models need to be fitted.
At the moment, using `train.spLearner` for point data set with >>1000 points should be done with caution.

The final results also shows quite similar results to universal kriging in [geoR](http://leg.ufpr.br/~paulojus/geoR/). The model error map above, however, shows more spatial contrast and helps detect areas of especially high errors.

<center><img src="https://github.com/thengl/GeoMLA/blob/master/RF_vs_kriging/results/rainfall/Fig_SIC1997_EML.png" width="900"></center> 
<center>_Figure: Predicted daily rainfall for the SIC1997 data set._</center>


The same function can also be used to interpolate factor-type variables:

```r
library(plotKML)
data(eberg_grid)
gridded(eberg_grid) <- ~x+y
proj4string(eberg_grid) <- CRS("+init=epsg:31467")
data(eberg)
coordinates(eberg) <- ~X+Y
proj4string(eberg) <- CRS("+init=epsg:31467")
X <- eberg_grid[c("PRMGEO6","DEMSRT6","TWISRT6","TIRAST6")]
mF <- train.spLearner(eberg["TAXGRSC"], covariates=X, buffer.dist=FALSE)
TAXGRSC <- predict(mF)
```

<center><img src="https://github.com/thengl/GeoMLA/blob/master/RF_vs_kriging/results/eberg/predicted_classes_eberg.png" width="900"></center> 
<center>_Figure: Predicted Ebergotzen soil types (probabilities)._</center>

Note that in the case of factor variables, prediction are based on ensemble stacking
based on the following three classification algorithms `c("classif.ranger", "classif.multinom", "classif.svm")`. See mlr documentation on how to add additional [learners](https://mlr.mlr-org.com/articles/tutorial/integrated_learners.html).

In summary: package mlr provides a comprehensive environment for Machine Learning:

- Ensemble predictions are based on the `mlr::makeStackedLearner` function,
- Additional [learners](https://mlr.mlr-org.com/articles/tutorial/integrated_learners.html) can be added,
- Processing can be parallelized using the [parallelMap package](https://mlr.mlr-org.com/articles/tutorial/parallelization.html),

Ensemble Machine Learning is also available via the [subsemble](https://github.com/ledell/subsemble) and the [SuperLearner](https://github.com/ecpolley/SuperLearner) packages (not used here).

### Accessing LandGIS layers

Landmap package also provides functionality to access and download LandGIS layers
from www.openlandmap.org. Recommend process is to first search the coverage ID 
names and file names e.g.:

```r
search.landgis(pattern=c("clay", "10..10cm"))
```

This shows that a clay map at 10 cm depth of the world is available via:

```
[[1]]
                                                   predicted250m.file 
"sol_clay.wfraction_usda.3a1a1a_m_250m_b10..10cm_1950..2017_v0.2.tif" 

[[2]]
[1] "https://www.zenodo.org/api/files/d95e82f3-203d-4ae5-86dc-b1ddd65ff8b2/sol_clay.wfraction_usda.3a1a1a_m_250m_b10..10cm_1950..2017_v0.2.tif" 
[2] "https://www.zenodo.org/api/files/d95e82f3-203d-4ae5-86dc-b1ddd65ff8b2/sol_clay.wfraction_usda.3a1a1a_md_250m_b10..10cm_1950..2017_v0.2.tif"
```

Web Coverage Service functionality and zenodo.org API are explained in detail [here](https://github.com/Envirometrix/LandGISmaps#accessing-data).
Next we can download only clay map for Switzerland using the Web Coverage Service 
functionality of LandGIS:

```r
coverageId = "predicted250m:sol_clay.wfraction_usda.3a1a1a_m_250m_b10..10cm_1950..2017_v0.2"
swiss1km.ll <- raster::projectRaster(raster(sic1997$swiss1km[2]), crs = "+init=epsg:4326", res=c(1/120, 1/120))
s1 = paste0("Lat(", swiss1km.ll@extent@ymin, ",", swiss1km.ll@extent@ymax,")")
s2 = paste0("Long(", swiss1km.ll@extent@xmin, ",", swiss1km.ll@extent@xmax,")")
sf = (1/480)/(1/120)
download.landgis(coverageId, filename = "clay_ch1km.tif", subset = c(s1,s2), scalefactor = sf)
swiss1km.ll1km <- as(swiss1km.ll, "SpatialGridDataFrame")
swiss1km.ll1km$clay_10..10cm <- readGDAL("clay_ch1km.tif")$band1
swiss1km.ll1km$clay_10..10cm <- ifelse(is.na(swiss1km.ll1km$DEM), NA, swiss1km.ll1km$clay_10..10cm)
mapview(swiss1km.ll1km["clay_10..10cm"])
```

<center><img src="https://github.com/thengl/GeoMLA/blob/master/RF_vs_kriging/results/rainfall/Fig_download_LandGIS_swiss1km.jpg" width="650"></center> 
<center>_Figure: Clay content map for Switzerland._</center>


This takes few steps because you have to determine:

* bounding box,
* scaling factor,
* mask out pixels of interest,

For smaller areas (<500Mb in size) download of data using WCS is fast and efficient.
For accessing and using global layers larger than 1GB we recommend directly downloading data from [zenodo.org](https://zenodo.org/search?page=1&size=20&q=LandGIS).

## Contributions

* Contributions to landmap are welcome. Issues and pull requests are the preferred ways of sharing them.
* We are interested in your results and experiences with using the `train.spLearner` function 
  for generating spatial predictions with your own data. Share your data sets, 
  code and results either using github issues and/or R-sig-geo mailing list.
