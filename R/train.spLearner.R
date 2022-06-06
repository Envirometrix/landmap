#' Train a spatial prediction and/or interpolation model using Ensemble Machine Learning
#' from a regression/classification matrix
#'
#' @param observations Data frame regression matrix,
#' @param formulaString Model formula,
#' @param covariates SpatialPixelsDataFrame object,
#' @param SL.library List of learners,
#' @param family Family e.g. \code{gaussian()},
#' @param method Ensemble stacking method (see makeStackedLearner),
#' @param predict.type Prediction type 'prob' or 'response',
#' @param super.learner Ensemble stacking model usually \code{regr.lm},
#' @param subsets Number of subsets for repeated CV,
#' @param lambda Target variable transformation lambda (0.5 or 1),
#' @param cov.model Covariance model for variogram fitting,
#' @param subsample For large datasets consider random subsetting training data,
#' @param parallel Initiate parellel processing,
#' @param cell.size Block size for spatial Cross-validation,
#' @param id Id column name to control clusters of data,
#' @param weights Optional weights (per row) that learners will use to account for variable data quality,
#' @param quantreg Fit additional ranger model as meta-learner to allow for derivation of prediction intervals,
#' @param ... other arguments that can be passed on to \code{mlr::makeStackedLearner},
#'
#' @return Object of class \code{spLearner}
#' @export
#'
#' @author \href{https://opengeohub.org/people/tom-hengl}{Tom Hengl}
train.spLearner.matrix <- function(observations, formulaString, covariates, SL.library, family = stats::gaussian(), method = "stack.cv", predict.type, super.learner, subsets = 5, lambda = 0.5, cov.model = "exponential", subsample = 10000, parallel = "multicore", cell.size, id = NULL, weights = NULL, quantreg = TRUE, ...){
  #if(!.Platform$OS.type=="unix") { parallel <- "seq" }
  tv <- all.vars(formulaString)[1]
  if(family$family == "binomial"){
    observations[,tv] <- as.factor(observations[,tv])
    if(!length(levels(observations[,tv]))==2){
      stop("Binary variable expected with two classes")
    }
  }
  Y <- observations[,tv]
  xyn <- attr(covariates@bbox, "dimnames")[[1]]
  if(missing(SL.library)){
    if(is.numeric(Y) & family$family == "gaussian"){
      if(is.null(weights)){
        SL.library <- c("regr.ranger", "regr.xgboost", "regr.nnet", "regr.ksvm", "regr.cvglmnet")
      } else {
        SL.library <- c("regr.ranger", "regr.xgboost", "regr.nnet")
      }
      if(missing(predict.type)){ predict.type <- "response" }
    }
    if(is.factor(Y) | family$family == "binomial"){
      SL.library <- c("classif.ranger", "classif.xgboost", "classif.nnTrain")
      if(missing(predict.type)){ predict.type <- "prob" }
    }
  }
  if(is.numeric(Y)){
    if(missing(super.learner)){ super.learner <- "regr.lm" }
  } else {
    if(missing(super.learner)){ super.learner <- "classif.glmnet" }
  }
  points <- observations
  sp::coordinates(points) <- stats::as.formula(paste("~", paste(xyn, collapse = "+"), sep=""))
  sp::proj4string(points) <- covariates@proj4string
  rvgm <- list(vgm=list(practicalRange=NA, cov.model="nugget", lambda=NA), observations=points)
  if(missing(cell.size)){
    ## automatically determine cell.size using fitted range
    cell.size <- abs(diff(covariates@bbox[1,]))/40
  }
  ## spatial ID for CV:
  if(is.null(id)){
    grd <- sp::GridTopology(cellcentre.offset=covariates@bbox[,1], cellsize=rep(cell.size,2), cells.dim=c(ceiling(abs(diff(covariates@bbox[1,])/cell.size)), ncols=ceiling(abs(diff(covariates@bbox[2,])/cell.size))))
    r.sp <- sp::SpatialGridDataFrame(grd, proj4string = covariates@proj4string, data=data.frame(gid=1:(grd@cells.dim[1] * grd@cells.dim[2])))
    id <- sp::over(sp::SpatialPoints(observations[,attr(covariates@bbox, "dimnames")[[1]]], proj4string = covariates@proj4string), r.sp)$gid
    id[is.na(id)] = "NULL overlay"
    id <- as.factor(make.names(id, unique=TRUE))
    message("Estimating block size ID for spatial Cross Validation...", immediate. = TRUE)
  }
  r.sel <- stats::complete.cases(observations[, all.vars(formulaString)])
  observations.s = observations[which(r.sel),all.vars(formulaString)]
  ## fit the mlr model:
  mlr::configureMlr()
  if(parallel=="multicore"){
    parallelMap::parallelStartSocket(parallel::detectCores())
  }
  message(paste0("Using learners: ", paste(SL.library, collapse = ", "), "..."), immediate. = TRUE)
  lrns <- lapply(SL.library, mlr::makeLearner)
  if(is.factor(Y) | family$family == "binomial"){
    message("Fitting a spatial learner using 'mlr::makeClassifTask'...", immediate. = TRUE)
    if(missing(predict.type)){ predict.type = "prob" }
    if(is.null(weights)){
      tsk <- mlr::makeClassifTask(data = observations.s, target = tv, coordinates = observations[which(r.sel),xyn])
    } else {
      tsk <- mlr::makeClassifTask(data = observations.s, target = tv, weights = weights[which(r.sel)], coordinates = observations[which(r.sel),xyn])
    }
    lrns <- lapply(lrns, mlr::setPredictType, "prob")
    init.m <- mlr::makeStackedLearner(base.learners = lrns, predict.type = predict.type, method = method, super.learner = super.learner, resampling=mlr::makeResampleDesc(method = "CV", blocking.cv=TRUE), ...)
  } else {
    message("Fitting a spatial learner using 'mlr::makeRegrTask'...", immediate. = TRUE)
    if(missing(predict.type)){ predict.type = "response" }
    if(any(SL.library %in% "regr.xgboost")){
      ## https://github.com/dmlc/xgboost/issues/4599
      lrns[[which(SL.library %in% "regr.xgboost")]] <- mlr::makeLearner("regr.xgboost", par.vals = list(objective ='reg:squarederror'))
    }
    if(is.null(weights)){
      tsk <- mlr::makeRegrTask(data = observations.s, target = tv, coordinates = observations[which(r.sel),xyn], blocking = id[which(r.sel)])
    } else {
      tsk <- mlr::makeRegrTask(data = observations.s, target = tv, weights = weights[which(r.sel)], coordinates = observations[which(r.sel),xyn], blocking = id[which(r.sel)])
    }
    init.m <- mlr::makeStackedLearner(base.learners = lrns, predict.type = predict.type, method = method, super.learner = super.learner, resampling=mlr::makeResampleDesc(method = "CV", blocking.cv=TRUE), ...)
  }
  suppressWarnings( m <- mlr::train(init.m, tsk) )
  if(quantreg ==TRUE & !is.factor(Y) & m$learner.model$super.model$learner$id=="regr.lm"){
    message("Fitting a quantreg model using 'ranger::ranger'...", immediate. = TRUE)
    quantregModel <- ranger::ranger(m$learner.model$super.model$learner.model$terms, m$learner.model$super.model$learner.model$model, num.trees=85, importance="impurity", quantreg=TRUE, keep.inbag = TRUE)
  } else {
    quantregModel = NULL
  }
  if(parallel=="multicore"){
    parallelMap::parallelStop()
  }
  out <- methods::new("spLearner", spModel = m, vgmModel = rvgm, covariates = covariates, spID = r.sp, quantregModel = quantregModel)
  return(out)
}

setMethod("train.spLearner", signature(observations = "data.frame", formulaString = "formula", covariates = "SpatialPixelsDataFrame"), train.spLearner.matrix)

#' Train a spatial prediction and/or interpolation model using Ensemble Machine Learning
#'
#' @rdname train.spLearner-methods
#' @aliases train.spLearner train.spLearner,data.frame,formula,SpatialPixelsDataFrame-method
#'
#' @description Automated spatial predictions and/or interpolation using Ensemble Machine Learning. Extends functionality of the \href{https://github.com/mlr-org/mlr}{mlr} package. Suitable for predicting numeric, binomial and factor-type variables.
#'
#' @param observations SpatialPointsDataFrame.
#' @param formulaString ANY.
#' @param covariates SpatialPixelsDataFrame.
#' @param SL.library List of learners,
#' @param family Family e.g. gaussian(),
#' @param method Ensemble stacking method (see makeStackedLearner) usually \code{stack.cv},
#' @param predict.type Prediction type 'prob' or 'response',
#' @param super.learner Ensemble stacking model usually \code{regr.lm},
#' @param subsets Number of subsets for repeated CV,
#' @param lambda Target variable transformation (0.5 or 1),
#' @param cov.model Covariance model for variogram fitting,
#' @param subsample For large datasets consider random subsetting training data,
#' @param parallel logical, Initiate parallel processing,
#' @param oblique.coords Specify whether to use oblique coordinates as covariates,
#' @param nearest Specify whether to use nearest values and distances i.e. the method of Sekulic et al. (2020),
#' @param buffer.dist Specify whether to use buffer distances to points as covariates,
#' @param theta.list List of angles (in radians) used to derive oblique coordinates,
#' @param spc specifies whether to apply principal components transformation.
#' @param id Id column name to control clusters of data,
#' @param weights Optional weights (per row) that learners will use to account for variable data quality,
#' @param n.obs Number of nearest observations to be found in \code{meteo::near.obs} (by default 10),
#' @param ... other arguments that can be passed on to \code{mlr::makeStackedLearner},
#'
#' @importClassesFrom sp SpatialPixelsDataFrame SpatialPointsDataFrame
#'
#' @return object of class \code{spLearner}, which contains fitted model, variogram model and spatial grid
#' used for Cross-validation.
#'
#' @author \href{https://opengeohub.org/people/tom-hengl}{Tom Hengl}
#'
#' @note By default uses oblique coordinates (rotated coordinates) as described in Moller et al. (2020; \doi{10.5194/soil-6-269-2020}) to account for geographical distribution of values.
#' By setting the \code{nearest = TRUE}, distances to nearest observations and values of nearest neighbors will be used (see: Sekulic et al, 2020; \doi{10.3390/rs12101687}). This method closely resembles geostatistical interpolators such as kriging.
#' Buffer geographical distances can be added by setting \code{buffer.dist=TRUE}.
#' Using either oblique coordinates and/or buffer distances is not recommended for point data set with distinct spatial clustering.
#' Effects of adding geographical distances into modeling are explained in detail in Hengl et al. (2018; \doi{10.7717/peerj.5518}) and Sekulic et al. (2020; \doi{10.3390/rs12101687}).
#' Default learners used for regression are: \code{c("regr.ranger", "regr.ksvm", "regr.nnet", "regr.cvglmnet")}.
#' Default learners used for classification / binomial variables are: \code{c("classif.ranger", "classif.svm", "classif.multinom")}, with \code{predict.type="prob"}.
#' When using \code{method = "stack.cv"} each training and prediction round could produce somewhat different results due to randomization of CV.
#' Prediction errors are derived by default using the \code{forestError} package method described in Lu & Hardin (2021).
#' Optionally, the quantreg (Quantile Regression) option from the ranger package (\href{https://jmlr.org/papers/v7/meinshausen06a.html}{Meinshausen, 2006}) can also be used.
#'
#' @references
#' \itemize{
#'   \item Moller, A. B., Beucher, A. M., Pouladi, N., and Greve, M. H. (2020). Oblique geographic coordinates as covariates for digital soil mapping. SOIL, 6, 269–289. \doi{10.5194/soil-6-269-2020}
#'   \item Hengl, T., Nussbaum, M., Wright, M. N., Heuvelink, G. B., and Graler, B. (2018) Random Forest as a generic framework for predictive modeling of spatial and spatio-temporal variables. PeerJ 6:e5518. \doi{10.7717/peerj.5518}
#'   \item Lu, B., & Hardin, J. (2021). A unified framework for random forest prediction error estimation. Journal of Machine Learning Research, 22(8), 1–41. \url{https://jmlr.org/papers/v22/18-558.html}
#'   \item Meinshausen, N. (2006). \href{https://jmlr.org/papers/v7/meinshausen06a.html}{Quantile regression forests}. Journal of Machine Learning Research, 7(Jun), 983–999. \url{https://jmlr.org/papers/v7/meinshausen06a.html}
#'   \item Sekulic, A., Kilibarda, M., Heuvelink, G. B., Nikolic, M. & Bajat, B. (2020). Random Forest Spatial Interpolation. Remote. Sens. 12, 1687, \doi{10.3390/rs12101687}
#' }
#'
#' @examples
#' library(rgdal)
#' library(mlr)
#' library(rpart)
#' library(nnet)
#' demo(meuse, echo=FALSE)
#' ## Regression:
#' sl = c("regr.rpart", "regr.nnet", "regr.glm")
#' system.time( m <- train.spLearner(meuse["lead"],
#'       covariates=meuse.grid[,c("dist","ffreq")],
#'       oblique.coords = FALSE, lambda=0,
#'       parallel=FALSE, SL.library=sl) )
#' summary(m@spModel$learner.model$super.model$learner.model)
#' \dontrun{
#' library(plotKML)
#' ## regression-matrix:
#' str(m@vgmModel$observations@data)
#' meuse.y <- predict(m, error.type="weighted.sd")
#' plot(raster::raster(meuse.y$pred["response"]), col=plotKML::R_pal[["rainbow_75"]][4:20],
#'    main="Predictions spLearner", axes=FALSE, box=FALSE)
#'
#' ## Regression with default settings:
#' m <- train.spLearner(meuse["zinc"], covariates=meuse.grid[,c("dist","ffreq")],
#'         parallel=FALSE, lambda = 0)
#' ## Ensemble model (meta-learner):
#' summary(m@spModel$learner.model$super.model$learner.model)
#' meuse.y <- predict(m)
#' ## Plot of predictions and prediction error (RMSPE)
#' op <- par(mfrow=c(1,2), oma=c(0,0,0,1), mar=c(0,0,4,3))
#' plot(raster::raster(meuse.y$pred["response"]), col=plotKML::R_pal[["rainbow_75"]][4:20],
#'    main="Predictions spLearner", axes=FALSE, box=FALSE)
#' points(meuse, pch="+")
#' plot(raster::raster(meuse.y$pred["model.error"]), col=rev(bpy.colors()),
#'    main="Prediction errors", axes=FALSE, box=FALSE)
#' points(meuse, pch="+")
#' par(op)
#' while (!is.null(dev.list())) dev.off()
#' ## Plot of prediction intervals:
#' pts = list("sp.points", meuse, pch = "+", col="black")
#' spplot(meuse.y$pred[,c("q.lwr","q.upr")], col.regions=plotKML::R_pal[["rainbow_75"]][4:20],
#'    sp.layout = list(pts),
#'    main="Prediction intervals (alpha = 0.318)")
#' while (!is.null(dev.list())) dev.off()
#'
#' ## Method from https://doi.org/10.3390/rs12101687
#' #library(meteo)
#' mN <- train.spLearner(meuse["zinc"], covariates=meuse.grid[,c("dist","ffreq")],
#'         parallel=FALSE, lambda=0, nearest=TRUE)
#' meuse.N <- predict(mN)
#' ## Plot of predictions and prediction error (RMSPE)
#' op <- par(mfrow=c(1,2), oma=c(0,0,0,1), mar=c(0,0,4,3))
#' plot(raster::raster(meuse.N$pred["response"]), col=plotKML::R_pal[["rainbow_75"]][4:20],
#'    main="Predictions spLearner meteo::near.obs", axes=FALSE, box=FALSE)
#' points(meuse, pch="+")
#' plot(raster::raster(meuse.N$pred["model.error"]), col=rev(bpy.colors()),
#'    main="Prediction errors", axes=FALSE, box=FALSE)
#' points(meuse, pch="+")
#' par(op)
#' while (!is.null(dev.list())) dev.off()
#'
#' ## Classification:
#' SL.library <- c("classif.ranger", "classif.xgboost", "classif.nnTrain")
#' mC <- train.spLearner(meuse["soil"], covariates=meuse.grid[,c("dist","ffreq")],
#'    SL.library = SL.library, super.learner = "classif.glmnet", parallel=FALSE)
#' meuse.soil <- predict(mC)
#' spplot(meuse.soil$pred[grep("prob.", names(meuse.soil$pred))],
#'         col.regions=plotKML::SAGA_pal[["SG_COLORS_YELLOW_RED"]], zlim=c(0,1))
#' spplot(meuse.soil$pred[grep("error.", names(meuse.soil$pred))],
#'          col.regions=rev(bpy.colors()))
#'
#' ## SIC1997
#' data("sic1997")
#' X <- sic1997$swiss1km[c("CHELSA_rainfall","DEM")]
#' mR <- train.spLearner(sic1997$daily.rainfall, covariates=X, lambda=1,
#'          nearest = TRUE, parallel=FALSE)
#' summary(mR@spModel$learner.model$super.model$learner.model)
#' rainfall1km <- predict(mR, what="mspe")
#' op <- par(mfrow=c(1,2), oma=c(0,0,0,1), mar=c(0,0,4,3))
#' plot(raster::raster(rainfall1km$pred["response"]), col=plotKML::R_pal[["rainbow_75"]][4:20],
#'     main="Predictions spLearner", axes=FALSE, box=FALSE)
#' points(sic1997$daily.rainfall, pch="+")
#' plot(raster::raster(rainfall1km$pred["model.error"]), col=rev(bpy.colors()),
#'     main="Prediction errors", axes=FALSE, box=FALSE)
#' points(sic1997$daily.rainfall, pch="+")
#' par(op)
#' while (!is.null(dev.list())) dev.off()
#'
#' ## Ebergotzen data set
#' data(eberg_grid)
#' gridded(eberg_grid) <- ~x+y
#' proj4string(eberg_grid) <- CRS("+init=epsg:31467")
#' data(eberg)
#' eb.s <- sample.int(nrow(eberg), 1400)
#' eberg <- eberg[eb.s,]
#' coordinates(eberg) <- ~X+Y
#' proj4string(eberg) <- CRS("+init=epsg:31467")
#' ## Binomial variable
#' summary(eberg$TAXGRSC)
#' eberg$Parabraunerde <- ifelse(eberg$TAXGRSC=="Parabraunerde", 1, 0)
#' X <- eberg_grid[c("PRMGEO6","DEMSRT6","TWISRT6","TIRAST6")]
#' mB <- train.spLearner(eberg["Parabraunerde"], covariates=X,
#'    family=binomial(), cov.model = "nugget", parallel=FALSE)
#' eberg.Parabraunerde <- predict(mB)
#' plot(raster::raster(eberg.Parabraunerde$pred["prob.1"]),
#'    col=plotKML::SAGA_pal[["SG_COLORS_YELLOW_RED"]], zlim=c(0,1))
#' points(eberg["Parabraunerde"], pch="+")
#'
#' ## Factor variable:
#' data(eberg)
#' coordinates(eberg) <- ~X+Y
#' proj4string(eberg) <- CRS("+init=epsg:31467")
#' X <- eberg_grid[c("PRMGEO6","DEMSRT6","TWISRT6","TIRAST6")]
#' mF <- train.spLearner(eberg["TAXGRSC"], covariates=X, parallel=FALSE)
#' TAXGRSC <- predict(mF)
#' plot(raster::stack(TAXGRSC$pred[grep("prob.", names(TAXGRSC$pred))]),
#'     col=plotKML::SAGA_pal[["SG_COLORS_YELLOW_RED"]], zlim=c(0,1))
#' plot(raster::stack(TAXGRSC$pred[grep("error.", names(TAXGRSC$pred))]),
#'     col=plotKML::SAGA_pal[["SG_COLORS_YELLOW_BLUE"]], zlim=c(0,0.45))
#' while (!is.null(dev.list())) dev.off()
#' }
#' @export
#' @docType methods
setMethod("train.spLearner", signature(observations = "SpatialPointsDataFrame", formulaString = "ANY", covariates = "SpatialPixelsDataFrame"), function(observations, formulaString, covariates, SL.library, family = stats::gaussian(), method = "stack.cv", predict.type, super.learner = "regr.lm", subsets = 5, lambda = 0.5, cov.model = "exponential", subsample = 10000, parallel = "multicore", oblique.coords = TRUE, nearest = FALSE, buffer.dist = FALSE, theta.list=seq(0, 180, length.out = 14)*pi/180, spc = TRUE, id = NULL, weights = NULL, n.obs = 10, ...){
  if(missing(formulaString)){
    tv <- names(observations)[1]
    observations <- observations[!is.na(observations@data[,tv]),]
    if(is.factor(observations@data[,1])){
      ## Drop levels with <5 observations
      xg <- summary(observations@data[,1], maxsum=length(levels(observations@data[,1])))
      selg.levs <- attr(xg, "names")[xg > 5]
      if(length(selg.levs)<length(levels(observations@data[,1]))){
        flev <- as.factor(observations@data[,1])
        flev[which(!observations@data[,1] %in% selg.levs)] <- NA
        observations@data[,1] <- droplevels(flev)
      }
      if(missing(predict.type)){ predict.type  <- "prob" }
    }
    if(spc==TRUE&ncol(covariates)>1){
      covariates <- spc(covariates)@predicted
      if(ncol(covariates)>4){
        ## remove last PC as this usually only carries noise?
        covariates <- covariates[-ncol(covariates)]
      }
    }
    if(nearest==TRUE){
      message("Deriving nearest observations using meteo::near.obs...", immediate. = TRUE)
      nearest_obs <- .near.obs(locations = covariates, observations = observations[tv], zcol = tv, n.obs = n.obs, rm.dupl = TRUE)
      nearest_obs.dev <- .near.obs(locations = observations[tv], observations = observations[tv], zcol = tv, n.obs = n.obs, rm.dupl = TRUE)
      covariates <- sp::cbind.Spatial(covariates, sp::SpatialPixelsDataFrame(covariates@coords, nearest_obs, proj4string = sp::CRS(sp::proj4string(covariates))))
      oblique.coords <- FALSE; buffer.dist <- FALSE
    }
    if(buffer.dist==TRUE){
      message("Deriving buffer distances to points...", immediate. = TRUE)
      classes <- as.factor(1:nrow(observations))
      covariates <- sp::cbind.Spatial(covariates, buffer.dist(observations[tv], covariates[1], classes))
    }
    if(oblique.coords==TRUE){
      message("Deriving oblique coordinates...", immediate. = TRUE)
      if(parallel=="multicore" & .Platform$OS.type=="unix"){
        #oblique.xy <- parallel::mclapply(theta.list, function(i){ spdep::Rotation(covariates[1]@coords, angle = i) }, mc.cores = parallel::detectCores())
        ## hangs
        oblique.xy <- lapply(theta.list, function(i){ spdep::Rotation(covariates[1]@coords, angle = i) })
      } else {
        oblique.xy <- list(NULL)
        for(i in 1:length(theta.list)){
          oblique.xy[[i]] <- spdep::Rotation(covariates[1]@coords, angle = theta.list[i])
        }
      }
      oblique.xy <- sp::SpatialPixelsDataFrame(covariates@coords, data=as.data.frame(do.call(cbind, oblique.xy)), proj4string = covariates@proj4string)
      R <- expand.grid(c("rX","rY"), round(theta.list, 1))
      names(oblique.xy) <- paste(R$Var1, R$Var2, sep="_")
      covariates <- sp::cbind.Spatial(covariates, oblique.xy)
    }
    formulaString <- stats::as.formula(paste(tv, " ~ ", paste(names(covariates), collapse="+")))
    if(length(all.vars(formulaString))>1000){
      warning("Number of covariarates exceeds 1000. Consider using 'classes' argument for buffer distances.", immediate. = TRUE)
    }
  } else {
    observations <- observations[!is.na(observations@data[,all.vars(formulaString)[1]]),]
    if(any(!all.vars(formulaString)[-1] %in% names(covariates))){
      stop("Covariates from 'formulaString' missing in the 'covariates' object")
    }
  }
  ov <- model.data(observations, formulaString, covariates)
  if(nearest==TRUE){
    ## replace values
    ov[,names(nearest_obs.dev)] <- nearest_obs.dev
  }
  r.sel <- stats::complete.cases(ov[,all.vars(formulaString)])
  n.r.sel <- sum(r.sel)
  if(!n.r.sel==nrow(observations)){
    message(paste0("Subsetting observations to ", round(n.r.sel/nrow(observations), 2)*100, "% complete cases..."), immediate. = TRUE)
  }
  if(nearest==TRUE){
    ## skip variogram modeling
    m <- train.spLearner.matrix(observations = ov, formulaString = formulaString, covariates = covariates, SL.library = SL.library, family = family, subsets = subsets, lambda = lambda, cov.model = "nugget", subsample = subsample, parallel = parallel, id=id, weights=weights, predict.type = predict.type, ...)
  } else {
    m <- train.spLearner.matrix(observations = ov, formulaString = formulaString, covariates = covariates, SL.library = SL.library, family = family, subsets = subsets, lambda = lambda, cov.model = cov.model, subsample = subsample, parallel = parallel, id=id, weights=weights, predict.type = predict.type, ...)
  }
  return(m)
})

#' Overlay points and grids and prepare regression matrix
#'
#' @param observations SpatialPointsDataFrame
#' @param formulaString Model definition
#' @param covariates List of covariates column names
#' @param dimensions 2D, 3D models
#'
#' @return Regression matrix data.frame
#' @export
model.data <- function(observations, formulaString, covariates, dimensions=c("2D","3D","2D+T","3D+T")){
  if(missing(dimensions)){ dimensions <- dimensions[1] }
  if(dimensions=="2D"){
    ## spatial overlay:
    ov <- sp::over(observations, covariates)
    ## subset to variables of interest:
    tv <- all.vars(formulaString)[1]
    seln <- names(covariates) %in% all.vars(formulaString)[-1]
    ## check if all covariates are available:
    if(sum(!is.na(seln))==0){
      stop("Covariates in the 'formulaString' do not match names in the 'covariates' object")
    }
    ov <- cbind(data.frame(observations[,tv]), ov)
    ## copy coordinate column names for consistency:
    xyn <- which(names(ov) %in% attr(observations@bbox, "dimnames")[[1]])
    if(is.null(attr(covariates@bbox, "dimnames")[[1]])) {
      dim.names = attr(observations@bbox, "dimnames")[[1]]
    } else {
      dim.names = attr(covariates@bbox, "dimnames")[[1]]
    }
    if(length(xyn)==2) {
      names(ov)[xyn] <- dim.names[1:2]
    } else {
      names(ov)[xyn] <- dim.names
    }
    ## check the size of output:
    if(nrow(ov)==0|is.null(ov[,tv])) {
      stop("The 'over' operations resulted in an empty set.")
    }
  }
  return(ov)
}


#' Print object of type 'spLearner'
#'
#' @param x of type \code{spLearner}
#' @param ... optional parameters
#'
#' @return Model summary
#' @method print spLearner
#' @export
"print.spLearner" <- function(x, ...){
  if(any(class(x@spModel)=="subsemble")){
    message("Subsemble model:")
    print(x@spModel$metafit$fit)
    message(paste0("CV R-square:", (1-x@spModel$metafit$fit$x$deviance/x@spModel$metafit$fit$x$null.deviance)))
  }
  if(any(class(x@spModel)=="BaseEnsembleModel")){
    message("BaseEnsembleModel:")
    if(any(class(x@spModel$learner.model$super.model$learner.model)=="lm")){
      print(summary(x@spModel$learner.model$super.model$learner.model))
    }
  }
  message("Variogram model:")
  print(x@vgmModel$vgm)
  message(paste0("Total observations: ", length(x@vgmModel$observations)))
}

#' Predict using spLearner at new locations
#'
#' @param object of type \code{spLearner}.
#' @param predictionLocations \code{SpatialPixelsDataFrame} with values of all features.
#' @param model.error Logical specify if prediction errors should be derived.
#' @param error.type Specify how should be the prediction error be derived.
#' @param t.prob Threshold probability for significant learners; only applyies for meta-learners based on lm model.
#' @param w optional weights vector.
#' @param quantiles Lower and upper quantiles for quantreg forest (0.159 and 0.841 for 1 standard deviation).
#' @param n.cores Number of cores to use (for parallel computation in ranger).
#' @param what A vector of characters indicating what estimates are desired for the quantForestError.
#' @param ... optional parameters.
#'
#' @return Object of class \code{SpatialPixelsDataFrame} with predictions and model error.
#' @method predict spLearner
#' @export
"predict.spLearner" <- function(object, predictionLocations, model.error=TRUE, error.type=c("forestError", "weighted.sd", "quantreg", "interval")[1], t.prob=1/3, w, quantiles = c((1-.682)/2, 1-(1-.682)/2), n.cores = parallel::detectCores(), what = c("mspe", "bias", "interval"), ...){
  if(any(object@spModel$task.desc$type=="classif")){
    error.type <- "weighted.sd"
  }
  if(missing(predictionLocations)){
    predictionLocations <- object@covariates
  }
  out <- predict(object@spModel, newdata=predictionLocations@data[,object@spModel$features])
  if(any(class(object@spModel)=="BaseEnsembleModel")){
    if(model.error==TRUE){
      message("Predicting values using 'getStackedBaseLearnerPredictions'...", immediate. = TRUE)
      out.c <- as.matrix(as.data.frame(mlr::getStackedBaseLearnerPredictions(object@spModel, newdata=predictionLocations@data[,object@spModel$features])))
      ## Classification:
      if(any(object@spModel$task.desc$type=="classif")){
        lvs <- object@spModel$task.desc$class.levels
        ## estimate weights:
        pred.prob <- NULL
        if(any(class(object@spModel$learner.model$super.model$learner.model) %in% "glmnet") & length(lvs)>2){
          w.v = .glmnet.varImp(object@spModel$learner.model$super.model$learner.model)
          ## model error
          message("Deriving model errors using sd of sign. learners...", immediate. = TRUE)
          for(j in lvs){
            out.c0 = out.c[,grep(paste0(".", j), attr(out.c, "dimnames")[[2]])]
            wt = w.v$x[match(attr(out.c0, "dimnames")[[2]], w.v$Group.1)]
            wt = ifelse(is.na(wt), min(wt, na.rm=TRUE), wt)
            pred.prob[[paste0("error.",j)]] <- matrixStats::rowWeightedSds(out.c0, w=wt, na.rm = TRUE)
          }
        } else {
          wt = rep(1, length(lvs))
          for(j in lvs){
            pred.prob[[paste0("error.",j)]] <- matrixStats::rowSds(out.c[,grep(paste0(".", j), attr(out.c, "dimnames")[[2]])], na.rm = TRUE)
          }
        }
        pred <- sp::SpatialPixelsDataFrame(predictionLocations@coords, data=cbind(out$data, data.frame(pred.prob)), grid = predictionLocations@grid, proj4string = predictionLocations@proj4string)
      ## Regression:
      } else {
        pred <- sp::SpatialPixelsDataFrame(predictionLocations@coords, data=out$data, grid=predictionLocations@grid, proj4string=predictionLocations@proj4string)
        if(error.type=="interval" & object@spModel$learner.model$super.model$learner$id=="regr.lm"){
          message("Deriving prediction errors...", immediate. = TRUE)
          ## http://www.sthda.com/english/articles/40-regression-analysis/166-predict-in-r-model-predictions-and-confidence-intervals/
          pred.int = predict(object@spModel$learner.model$super.model$learner.model, newdata = data.frame(out.c), interval = "prediction", level=2/3)
          ## assumes normal distribution
          pred$model.error <- (pred.int[,"upr"]-pred.int[,"lwr"])/2
        } else {
          m.train = object@spModel$learner.model$super.model$learner.model$model
          m.terms = all.vars(object@spModel$learner.model$super.model$learner.model$terms)
          if((error.type=="quantreg" | error.type=="forestError") & object@spModel$learner.model$super.model$learner$id=="regr.lm"){
            if(error.type=="forestError"){
              message("Deriving model errors using forestError package...", immediate. = TRUE)
              ## http://jmlr.org/papers/v22/18-558.html
              pred.q = forestError::quantForestError(object@quantregModel, X.train=m.train[,m.terms[-1]], X.test=as.data.frame(out.c), Y.train=m.train[,m.terms[1]], what = what, alpha = (1-(quantiles[2]-quantiles[1])), n.cores = n.cores)
              if(any(what %in% "mspe")){ pred$model.error <- sqrt(pred.q$mspe) }
              if(any(what %in% "bias")){ pred$model.bias <- pred.q$bias }
              if(any(what %in% "interval")){
                pred@data[,"q.lwr"] <- pred.q[,grep("lower", names(pred.q))]
                pred@data[,"q.upr"] <- pred.q[,grep("upper", names(pred.q))]
              }
            } else {
              message("Deriving model errors using ranger package 'quantreg' option...", immediate. = TRUE)
              pred.q = predict(object@quantregModel, as.data.frame(out.c), type="quantiles", quantiles=quantiles)
              pred$model.error <- (pred.q$predictions[,2]-pred.q$predictions[,1])/2
              pred@data[,"q.lwr"] <- pred.q$predictions[,1]
              pred@data[,"q.upr"] <- pred.q$predictions[,2]
            }
          } else {
            message("Deriving model errors using sd of sign. learners...", immediate. = TRUE)
            ## Linear Models: the absolute value of the t-statistic for each model parameter is used: https://topepo.github.io/caret/variable-importance.html
            ## If coefficient t-value probability >0.05 don't use to derive errors because it could inflate errors
            ## This is an arbitrary decision
            wt <- summary(object@spModel$learner.model$super.model$learner.model)$coefficients[-1,4]
            wt <- which(wt<t.prob)
            if(length(wt)>1){
              ## correction factor:
              eml.MSE0 = matrixStats::rowSds(as.matrix(m.train[,m.terms[-1]]), na.rm=TRUE)^2
              eml.MSE = stats::deviance(object@spModel$learner.model$super.model$learner.model)/stats::df.residual(object@spModel$learner.model$super.model$learner.model)
              eml.cf = eml.MSE/mean(eml.MSE0, na.rm = TRUE)
              pred$model.error <- sqrt(matrixStats::rowSds(out.c[,wt], na.rm=TRUE)^2 * eml.cf)
            } else {
              warning("Number of significant learners <2. Try adding additional learners.")
              pred$model.error <- NA
            }
          }
        }
      }
    return(out <- list(pred=pred, subpred=as.data.frame(out.c), quantiles=quantiles) )
    } else {
      pred <- sp::SpatialPixelsDataFrame(predictionLocations@coords, data=out$data, grid=predictionLocations@grid, proj4string = predictionLocations@proj4string)
      return(out <- list(pred=pred, subpred=NULL, quantiles=NULL))
    }
  }
}

## https://stackoverflow.com/questions/35461839/glmnet-variable-importance
.glmnet.varImp <- function(x, ...) {
  coefList <- glmnet::coef.glmnet(x)
  coefList <- do.call("rbind", lapply(coefList, function(X){data.frame(X@Dimnames[[1]][X@i+1],X@x)}))
  names(coefList) <- c('var','val')
  coefList$aval = abs(coefList$val)
  out = stats::aggregate(coefList$aval, by=list(coefList$var), FUN="mean", na.action=stats::na.omit)
  out = out[!out$Group.1=="",]
  out
}

setMethod("print", signature(x = "spLearner"), print.spLearner)
setMethod("predict", signature(object = "spLearner"), predict.spLearner)

## https://rdrr.io/rforge/meteo/src/R/near.obs.R
.near.obs <- function(
  locations,
  locations.x.y = c(1,2),
  observations,
  observations.x.y = c(1,2),
  zcol = 3,
  n.obs = 10,
  rm.dupl = TRUE
)
{
  if (class(locations) == "SpatialPoints" ||
      class(locations) == "SpatialPointsDataFrame" ||
      class(locations) == "SpatialPixelsDataFrame") {
    locations <- sp::coordinates(locations)
  } else {
    locations <- locations[, locations.x.y]
  }
  if (class(observations) == "SpatialPoints" || class(observations) == "SpatialPointsDataFrame") {
    variable <- observations[[zcol]]
    observations <- sp::coordinates(observations)
  } else {
    variable <- observations[, zcol]
    observations <- observations[, observations.x.y]
  }

  if (nrow(observations) < (n.obs+1)) {
    # return NA
    nl_df <- matrix(NA, nrow = nrow(locations), ncol = (2*n.obs))
  } else {
    # if (identical(locations, observations)){
    if (rm.dupl){
      knn1 <- nabor::knn(observations, locations, k=n.obs+1)
      knn1$nn.idx[round(knn1$nn.dists[, 1]) == 0, 1:n.obs] <- knn1$nn.idx[round(knn1$nn.dists[, 1]) == 0, -1]
      knn1$nn.idx <- knn1$nn.idx[, -(n.obs+1)]
      knn1$nn.dists[round(knn1$nn.dists[, 1]) == 0, 1:n.obs] <- knn1$nn.dists[round(knn1$nn.dists[, 1]) == 0, -1]
      knn1$nn.dists <- knn1$nn.dists[, -(n.obs+1)]
    } else {
      knn1 <- nabor::knn(observations, locations, k=n.obs)
    }

    if(!all(class(knn1$nn.idx)=='integer')) {
      near_o1 <- apply(knn1$nn.idx, 2, function(x) {variable[x]})
      near_o1 <- cbind(near_o1)
      nl_df <- cbind(knn1$nn.dists, near_o1)
    } else {
      near_o1 <- variable[knn1$nn.idx]
      nl_df <- t(c(knn1$nn.dists, near_o1))
    }
  }

  name1 <- c()
  name2 <- c()
  for (i in 1:n.obs) {
    name1 <- c(name1, paste("dist", i, sep = ""))
    name2 <- c(name2, paste("obs", i, sep = ""))
  }
  all_names <- c(name1, name2)
  nl_df <- as.data.frame(nl_df)
  names(nl_df) <- all_names

  return(nl_df)

}
