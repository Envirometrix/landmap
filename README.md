# landmap package for R

Provides methodology for automated mapping i.e. spatial interpolation and/or 
prediction using Ensemble Machine Learning (extends functionality of the 
[subsemble](https://github.com/ledell/subsemble) and the [SuperLearner](https://github.com/ecpolley/SuperLearner) packages). Key functionality includes:

* `train.spLearner` --- train a spatial prediction and/or interpolation model using Ensemble Machine Learning,
* `buffer.dist` --- derive buffer (geographical) distances that can be used as covariates in spLearner, 
* `spc` --- derive Principal Components using stack of spatial layers,
* `tile` --- tile spatial layers so they can be used to run processing in parallel,
* `download.landgis` --- access and download LandGIS layers from www.openlandmap.org,

Spatial prediction using Ensemble Machine Learning methodology is explained in 
detail in:

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

## Functionality

### Automated mapping using Ensemble Machine Learning

The following examples demostrates spatial prediction using the meuse data set:

```r
library(rgdal)
library(geoR)
library(plotKML)
library(raster)
library(SuperLearner)
library(subsemble)
demo(meuse, echo=FALSE)
m <- train.spLearner(meuse["lead"], covariates=meuse.grid[,c("dist","ffreq")], lambda = 1)
```

this runs several steps:

```
Converting ffreq to indicators...
Converting covariates to principal components...
Deriving buffer distances to points...TRUE
Fitting a variogram using 'linkfit' and trend model...TRUE
Using block size ID for spatial Cross Validation...TRUE
Fitting a spatial learner using 'subsemble'...TRUE
Loading required package: parallel
```

Note that the variogram model is only fitted to estimate effective range of spatial dependence.
Spatial Prediction models are based only on fitting the Ensemble Machine Learning 
(by default uses `c("SL.xgboost", "SL.ranger", "SL.ksvm")`) with geographical distances 
as additional covariates. To check modelling success we can look at the summary model properties:

```r
print(m)
```

which shows that the Cross Validation R-square is about 55%. Next we can generate predictions using:

```r
meuse.lead <- predict(m)
```

![figure](https://github.com/thengl/GeoMLA/blob/master/RF_vs_kriging/results/meuse/Fig_meuse_EML.png) *Figure: Predicted lead content for the Meuse data set.*

Notice that the predictions reflect spatial correlation between values, and hence can be used
as a possible replacement for kriging methods (read more in [Hengl et al. 2018](https://doi.org/10.7717/peerj.5518)). Automation comes, however, at the high computing and RAM usage costs.

In the following example we use somewhat larger data set from the SIC1997 exercise.

```r
data("sic1997")
mR <- train.spLearner(sic1997$daily.rainfall, covariates=sic1997$swiss1km[c("CHELSA_rainfall","DEM")], lambda=1)
rainfall1km <- predict(mR)
```

The processing is much more computational because the data set consists from 467 points.
This will make the regression matrix becoming extensive, and also 5x3 models need to be fitted.
At the moment we do not recommend using `train.spLearner` for point data set with >1000 points.

The final results also shows quite similar results to universal kriging in geoR.
The model error map, however, shows more spatial contrast and helps detect areas of 
especially high errors.

![figure](https://github.com/thengl/GeoMLA/blob/master/RF_vs_kriging/results/rainfall/Fig_SIC1997_EML.png) *Figure: Predicted daily rainfall for the SIC1997 data set.*

The same function can also be used to interpolate facto-type variables:

```r
library(plotKML)
data(eberg_grid)
gridded(eberg_grid) <- ~x+y
proj4string(eberg_grid) <- CRS("+init=epsg:31467")
data(eberg)
coordinates(eberg) <- ~X+Y
proj4string(eberg) <- CRS("+init=epsg:31467")
mF <- train.spLearner(eberg["TAXGRSC"], covariates=eberg_grid[c("PRMGEO6","DEMSRT6","TWISRT6","TIRAST6")])
TAXGRSC <- predict(mF)
```

![figure](https://github.com/Envirometrix/PredictiveSoilMapping/blob/master/figures/predicted_classes_eberg.png) *Figure: Predicted Ebergotzen soil types (probabilities).*

Note that in the case of factor variables, prediction are based on simple average from
`ranger::range`, `e1071::svm` and `nnet::multinom`, which might be suboptimal.

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

![figure](https://github.com/thengl/GeoMLA/blob/master/RF_vs_kriging/results/rainfall/Fig_download_LandGIS_swiss1km.jpg) *Figure: Clay content map for Switzerland.*

This takes few steps because we have to determine:

* bounding box,
* scaling factor,
* mask out pixels of interest,

For smaller areas (<500Mb in size) download of data using WCS is fast and efficient.
For accessing and using global layers larger than 1GB we recommend directly downloading data from zenodo.org.

## Contributions

* Contributions to landmap are most welcome. Issues and pull requests are the preferred ways of sharing them.
* We are interested in results and experiences of using the train.spLearner function for generating spatial predictions with your own data sets.
