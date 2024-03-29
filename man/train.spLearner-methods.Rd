% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/train.spLearner.R
\docType{methods}
\name{train.spLearner,SpatialPointsDataFrame,ANY,SpatialPixelsDataFrame-method}
\alias{train.spLearner,SpatialPointsDataFrame,ANY,SpatialPixelsDataFrame-method}
\alias{train.spLearner}
\alias{train.spLearner,data.frame,formula,SpatialPixelsDataFrame-method}
\title{Train a spatial prediction and/or interpolation model using Ensemble Machine Learning}
\usage{
\S4method{train.spLearner}{SpatialPointsDataFrame,ANY,SpatialPixelsDataFrame}(
  observations,
  formulaString,
  covariates,
  SL.library,
  family = stats::gaussian(),
  method = "stack.cv",
  predict.type,
  super.learner = "regr.lm",
  subsets = 5,
  lambda = 0.5,
  cov.model = "exponential",
  subsample = 10000,
  parallel = "multicore",
  oblique.coords = TRUE,
  nearest = FALSE,
  buffer.dist = FALSE,
  theta.list = seq(0, 180, length.out = 14) * pi/180,
  spc = TRUE,
  id = NULL,
  weights = NULL,
  n.obs = 10,
  ...
)
}
\arguments{
\item{observations}{SpatialPointsDataFrame.}

\item{formulaString}{ANY.}

\item{covariates}{SpatialPixelsDataFrame.}

\item{SL.library}{List of learners,}

\item{family}{Family e.g. gaussian(),}

\item{method}{Ensemble stacking method (see makeStackedLearner) usually \code{stack.cv},}

\item{predict.type}{Prediction type 'prob' or 'response',}

\item{super.learner}{Ensemble stacking model usually \code{regr.lm},}

\item{subsets}{Number of subsets for repeated CV,}

\item{lambda}{Target variable transformation (0.5 or 1),}

\item{cov.model}{Covariance model for variogram fitting,}

\item{subsample}{For large datasets consider random subsetting training data,}

\item{parallel}{logical, Initiate parallel processing,}

\item{oblique.coords}{Specify whether to use oblique coordinates as covariates,}

\item{nearest}{Specify whether to use nearest values and distances i.e. the method of Sekulic et al. (2020),}

\item{buffer.dist}{Specify whether to use buffer distances to points as covariates,}

\item{theta.list}{List of angles (in radians) used to derive oblique coordinates,}

\item{spc}{specifies whether to apply principal components transformation.}

\item{id}{Id column name to control clusters of data,}

\item{weights}{Optional weights (per row) that learners will use to account for variable data quality,}

\item{n.obs}{Number of nearest observations to be found in \code{meteo::near.obs} (by default 10),}

\item{...}{other arguments that can be passed on to \code{mlr::makeStackedLearner},}
}
\value{
object of class \code{spLearner}, which contains fitted model, variogram model and spatial grid
used for Cross-validation.
}
\description{
Automated spatial predictions and/or interpolation using Ensemble Machine Learning. Extends functionality of the \href{https://github.com/mlr-org/mlr}{mlr} package. Suitable for predicting numeric, binomial and factor-type variables.
}
\note{
By default uses oblique coordinates (rotated coordinates) as described in Moller et al. (2020; \doi{10.5194/soil-6-269-2020}) to account for geographical distribution of values.
By setting the \code{nearest = TRUE}, distances to nearest observations and values of nearest neighbors will be used (see: Sekulic et al, 2020; \doi{10.3390/rs12101687}). This method closely resembles geostatistical interpolators such as kriging.
Buffer geographical distances can be added by setting \code{buffer.dist=TRUE}.
Using either oblique coordinates and/or buffer distances is not recommended for point data set with distinct spatial clustering.
Effects of adding geographical distances into modeling are explained in detail in Hengl et al. (2018; \doi{10.7717/peerj.5518}) and Sekulic et al. (2020; \doi{10.3390/rs12101687}).
Default learners used for regression are: \code{c("regr.ranger", "regr.ksvm", "regr.nnet", "regr.cvglmnet")}.
Default learners used for classification / binomial variables are: \code{c("classif.ranger", "classif.svm", "classif.multinom")}, with \code{predict.type="prob"}.
When using \code{method = "stack.cv"} each training and prediction round could produce somewhat different results due to randomization of CV.
Prediction errors are derived by default using the \code{forestError} package method described in Lu & Hardin (2021).
Optionally, the quantreg (Quantile Regression) option from the ranger package (\href{https://jmlr.org/papers/v7/meinshausen06a.html}{Meinshausen, 2006}) can also be used.
}
\examples{
library(rgdal)
library(mlr)
library(rpart)
library(nnet)
demo(meuse, echo=FALSE)
## Regression:
sl = c("regr.rpart", "regr.nnet", "regr.glm")
system.time( m <- train.spLearner(meuse["lead"],
      covariates=meuse.grid[,c("dist","ffreq")],
      oblique.coords = FALSE, lambda=0,
      parallel=FALSE, SL.library=sl) )
summary(m@spModel$learner.model$super.model$learner.model)
\dontrun{
library(plotKML)
## regression-matrix:
str(m@vgmModel$observations@data)
meuse.y <- predict(m, error.type="weighted.sd")
plot(raster::raster(meuse.y$pred["response"]), col=plotKML::R_pal[["rainbow_75"]][4:20],
   main="Predictions spLearner", axes=FALSE, box=FALSE)

## Regression with default settings:
m <- train.spLearner(meuse["zinc"], covariates=meuse.grid[,c("dist","ffreq")],
        parallel=FALSE, lambda = 0)
## Ensemble model (meta-learner):
summary(m@spModel$learner.model$super.model$learner.model)
meuse.y <- predict(m)
## Plot of predictions and prediction error (RMSPE)
op <- par(mfrow=c(1,2), oma=c(0,0,0,1), mar=c(0,0,4,3))
plot(raster::raster(meuse.y$pred["response"]), col=plotKML::R_pal[["rainbow_75"]][4:20],
   main="Predictions spLearner", axes=FALSE, box=FALSE)
points(meuse, pch="+")
plot(raster::raster(meuse.y$pred["model.error"]), col=rev(bpy.colors()),
   main="Prediction errors", axes=FALSE, box=FALSE)
points(meuse, pch="+")
par(op)
while (!is.null(dev.list())) dev.off()
## Plot of prediction intervals:
pts = list("sp.points", meuse, pch = "+", col="black")
spplot(meuse.y$pred[,c("q.lwr","q.upr")], col.regions=plotKML::R_pal[["rainbow_75"]][4:20],
   sp.layout = list(pts),
   main="Prediction intervals (alpha = 0.318)")
while (!is.null(dev.list())) dev.off()

## Method from https://doi.org/10.3390/rs12101687
#library(meteo)
mN <- train.spLearner(meuse["zinc"], covariates=meuse.grid[,c("dist","ffreq")],
        parallel=FALSE, lambda=0, nearest=TRUE)
meuse.N <- predict(mN)
## Plot of predictions and prediction error (RMSPE)
op <- par(mfrow=c(1,2), oma=c(0,0,0,1), mar=c(0,0,4,3))
plot(raster::raster(meuse.N$pred["response"]), col=plotKML::R_pal[["rainbow_75"]][4:20],
   main="Predictions spLearner meteo::near.obs", axes=FALSE, box=FALSE)
points(meuse, pch="+")
plot(raster::raster(meuse.N$pred["model.error"]), col=rev(bpy.colors()),
   main="Prediction errors", axes=FALSE, box=FALSE)
points(meuse, pch="+")
par(op)
while (!is.null(dev.list())) dev.off()

## Classification:
SL.library <- c("classif.ranger", "classif.xgboost", "classif.nnTrain")
mC <- train.spLearner(meuse["soil"], covariates=meuse.grid[,c("dist","ffreq")],
   SL.library = SL.library, super.learner = "classif.glmnet", parallel=FALSE)
meuse.soil <- predict(mC)
spplot(meuse.soil$pred[grep("prob.", names(meuse.soil$pred))],
        col.regions=plotKML::SAGA_pal[["SG_COLORS_YELLOW_RED"]], zlim=c(0,1))
spplot(meuse.soil$pred[grep("error.", names(meuse.soil$pred))],
         col.regions=rev(bpy.colors()))

## SIC1997
data("sic1997")
X <- sic1997$swiss1km[c("CHELSA_rainfall","DEM")]
mR <- train.spLearner(sic1997$daily.rainfall, covariates=X, lambda=1,
         nearest = TRUE, parallel=FALSE)
summary(mR@spModel$learner.model$super.model$learner.model)
rainfall1km <- predict(mR, what="mspe")
op <- par(mfrow=c(1,2), oma=c(0,0,0,1), mar=c(0,0,4,3))
plot(raster::raster(rainfall1km$pred["response"]), col=plotKML::R_pal[["rainbow_75"]][4:20],
    main="Predictions spLearner", axes=FALSE, box=FALSE)
points(sic1997$daily.rainfall, pch="+")
plot(raster::raster(rainfall1km$pred["model.error"]), col=rev(bpy.colors()),
    main="Prediction errors", axes=FALSE, box=FALSE)
points(sic1997$daily.rainfall, pch="+")
par(op)
while (!is.null(dev.list())) dev.off()

## Ebergotzen data set
data(eberg_grid)
gridded(eberg_grid) <- ~x+y
proj4string(eberg_grid) <- CRS("+init=epsg:31467")
data(eberg)
eb.s <- sample.int(nrow(eberg), 1400)
eberg <- eberg[eb.s,]
coordinates(eberg) <- ~X+Y
proj4string(eberg) <- CRS("+init=epsg:31467")
## Binomial variable
summary(eberg$TAXGRSC)
eberg$Parabraunerde <- ifelse(eberg$TAXGRSC=="Parabraunerde", 1, 0)
X <- eberg_grid[c("PRMGEO6","DEMSRT6","TWISRT6","TIRAST6")]
mB <- train.spLearner(eberg["Parabraunerde"], covariates=X,
   family=binomial(), cov.model = "nugget", parallel=FALSE)
eberg.Parabraunerde <- predict(mB)
plot(raster::raster(eberg.Parabraunerde$pred["prob.1"]),
   col=plotKML::SAGA_pal[["SG_COLORS_YELLOW_RED"]], zlim=c(0,1))
points(eberg["Parabraunerde"], pch="+")

## Factor variable:
data(eberg)
coordinates(eberg) <- ~X+Y
proj4string(eberg) <- CRS("+init=epsg:31467")
X <- eberg_grid[c("PRMGEO6","DEMSRT6","TWISRT6","TIRAST6")]
mF <- train.spLearner(eberg["TAXGRSC"], covariates=X, parallel=FALSE)
TAXGRSC <- predict(mF)
plot(raster::stack(TAXGRSC$pred[grep("prob.", names(TAXGRSC$pred))]),
    col=plotKML::SAGA_pal[["SG_COLORS_YELLOW_RED"]], zlim=c(0,1))
plot(raster::stack(TAXGRSC$pred[grep("error.", names(TAXGRSC$pred))]),
    col=plotKML::SAGA_pal[["SG_COLORS_YELLOW_BLUE"]], zlim=c(0,0.45))
while (!is.null(dev.list())) dev.off()
}
}
\references{
\itemize{
  \item Moller, A. B., Beucher, A. M., Pouladi, N., and Greve, M. H. (2020). Oblique geographic coordinates as covariates for digital soil mapping. SOIL, 6, 269–289. \doi{10.5194/soil-6-269-2020}
  \item Hengl, T., Nussbaum, M., Wright, M. N., Heuvelink, G. B., and Graler, B. (2018) Random Forest as a generic framework for predictive modeling of spatial and spatio-temporal variables. PeerJ 6:e5518. \doi{10.7717/peerj.5518}
  \item Lu, B., & Hardin, J. (2021). A unified framework for random forest prediction error estimation. Journal of Machine Learning Research, 22(8), 1–41. \url{https://jmlr.org/papers/v22/18-558.html}
  \item Meinshausen, N. (2006). \href{https://jmlr.org/papers/v7/meinshausen06a.html}{Quantile regression forests}. Journal of Machine Learning Research, 7(Jun), 983–999. \url{https://jmlr.org/papers/v7/meinshausen06a.html}
  \item Sekulic, A., Kilibarda, M., Heuvelink, G. B., Nikolic, M. & Bajat, B. (2020). Random Forest Spatial Interpolation. Remote. Sens. 12, 1687, \doi{10.3390/rs12101687}
}
}
\author{
\href{https://opengeohub.org/people/tom-hengl/}{Tom Hengl}
}
