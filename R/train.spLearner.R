#' Train a spatial prediction and/or interpolation model using Ensemble Machine Learning
#' from a regression/classification matrix
#'
#' @param observations Data frame regression matrix
#' @param formulaString Model formula
#' @param covariates SpatialPixelsDataFrame object
#' @param SL.library List of learners
#' @param family Family e.g. gaussian()
#' @param method Ensemble stacking method (see makeStackedLearner)
#' @param predict.type Prediction type 'prob' or 'response'
#' @param super.learner Ensemble stacking model
#' @param subsets Number of subsets for repeated CV
#' @param cvControl Cross-validation setting (number of folds)
#' @param lambda Target variable transformation (0.5 or 1)
#' @param cov.model Covariance model for variogram fitting
#' @param subsample For large datasets consider random subsetting training data
#' @param parallel Initiate parellel processing
#' @param cell.size Block size for spatial Cross-validation
#' @param id Id column name to control clusters of data
#' @param weights Optional weights (per row) that learners will use to account for variable data quality
#' @param ...
#'
#' @return
#' @export
#'
#' @examples
train.spLearner.matrix <- function(observations, formulaString, covariates, SL.library, family = gaussian(), method = "stack.cv", predict.type, super.learner = "regr.glmnet", subsets = 3, cvControl = 5, lambda = 0.5, cov.model = "exponential", subsample = 5000, parallel = "multicore", cell.size, id = NULL, weights = NULL, ...){
  if(!.Platform$OS.type=="unix") { parallel <- "seq" }
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
      SL.library <- c("regr.xgboost", "regr.ranger", "regr.ksvm", "regr.glmnet")
      if(missing(predict.type)){ predict.type <- "response" }
    }
    if(is.factor(Y) | family$family == "binomial"){
      SL.library <- c("classif.ranger", "classif.multinom", "classif.svm")
      if(missing(predict.type)){ predict.type <- "prob" }
    }
  }
  if(is.numeric(Y)){
    ## variogram fitting:
    if(lambda==1){
      ini.var <- var(log1p(Y), na.rm = TRUE)
    }
    if(lambda==0.5){
      ini.var <- var(Y, na.rm = TRUE)
    }
    ## strip buffer distances from vgm modeling
    covs.vgm <- names(covariates)[-grep(pattern=glob2rx("layer.*$"), names(covariates))]
    formulaString.vgm <- as.formula(paste(tv, "~", paste(covs.vgm, collapse="+")))
    rvgm <- fit.vgmModel(formulaString.vgm, rmatrix = observations, predictionDomain = covariates[covs.vgm], lambda = lambda, ini.var = ini.var, cov.model = cov.model, subsample = subsample)
  } else {
    message("Skipping variogram modeling...", immediate. = TRUE)
    points <- observations
    coordinates(points) <- as.formula(paste("~", paste(xyn, collapse = "+"), sep=""))
    proj4string(points) <- covariates@proj4string
    rvgm <- list(vgm=list(practicalRange=NA, cov.model="nugget", lambda=NA), observations=points)
  }
  if(missing(cell.size)){
    ## automatically determine cell.size using fitted range
    cell.size <- rvgm$vgm$practicalRange/2
    if(is.na(cell.size) | cell.size < (covariates@grid@cellsize[[1]]*2)){
      cell.size <- abs(diff(covariates@bbox[1,]))/20
    }
  }
  ## spatial ID for CV:
  if(is.null(id)){
    grd <- sp::GridTopology(cellcentre.offset=covariates@bbox[,1], cellsize=rep(cell.size,2), cells.dim=c(ceiling(abs(diff(covariates@bbox[1,])/cell.size)), ncols=ceiling(abs(diff(covariates@bbox[2,])/cell.size))))
    r.sp <- sp::SpatialGridDataFrame(grd, proj4string = covariates@proj4string, data=data.frame(gid=1:(grd@cells.dim[1] * grd@cells.dim[2])))
    id <- sp::over(sp::SpatialPoints(observations[,attr(covariates@bbox, "dimnames")[[1]]], proj4string = covariates@proj4string), r.sp)$gid
    id[is.na(id)] = "NULL overlay"
    id <- as.factor(id)
    message("Estimating block size ID for spatial Cross Validation...", immediate. = TRUE)
  }
  ## fit the model:
  r.sel <- complete.cases(observations[, all.vars(formulaString)])
  if(parallel=="multicore"){
    parallelMap::parallelStartSocket(parallel::detectCores())
  }
  if(is.factor(Y) | family$family == "binomial"){
    message("Fitting a spatial learner using 'mlr::makeClassifTask'...", immediate. = TRUE)
    tsk <- mlr::makeClassifTask(data = observations[which(r.sel),all.vars(formulaString)], target = tv, weights = weights[which(r.sel)], coordinates = observations[which(r.sel),xyn])
    lrns <- lapply(SL.library, mlr::makeLearner)
    lrns <- lapply(lrns, setPredictType, "prob")
    init.m <- mlr::makeStackedLearner(base.learners = lrns, predict.type = predict.type, method = "hill.climb", ...)
  } else {
    message("Fitting a spatial learner using 'mlr::makeRegrTask'...", immediate. = TRUE)
    tsk <- mlr::makeRegrTask(data = observations[which(r.sel),all.vars(formulaString)], target = tv, weights = weights[which(r.sel)], coordinates = observations[which(r.sel),xyn], blocking = id[which(r.sel)])
    lrns <- lapply(SL.library, mlr::makeLearner)
    init.m <- mlr::makeStackedLearner(base.learners = lrns, predict.type = predict.type, method = method, super.learner = super.learner, resampling = mlr::makeResampleDesc(method="CV", iters = cvControl), use.feat=TRUE, ...)
  }
  m <- mlr::train(init.m, tsk)
  if(parallel=="multicore"){
    parallelMap::parallelStop()
  }
  out <- new("spLearner", spModel = m, vgmModel = rvgm, covariates = covariates, spID = r.sp)
  return(out)
}

setMethod("train.spLearner", signature(observations = "data.frame", formulaString = "formula", covariates = "SpatialPixelsDataFrame"), train.spLearner.matrix)

#' Train a spatial prediction and/or interpolation model using Ensemble Machine Learning
#'
#' @description Automated spatial predictions and/or interpolation using Ensemble Machine Learning. Extends functionality of the \href{https://github.com/mlr-org/mlr}{mlr} package. Suitable for predicting numeric, binomial and factor-type variables.
#'
#' @param observations SpatialPointsDataFrame.
#' @param formulaString ANY.
#' @param covariates SpatialPixelsDataFrame.
#'
#' @importClassesFrom sp SpatialPixelsDataFrame SpatialPointsDataFrame
#'
#' @return object of class spLearner, which contains fitted model, variogram model and spatial grid
#' used for Cross-validation.
#'
#' @author \href{https://opengeohub.org/people/tom-hengl}{Tom Hengl}
#'
#' @note Incorporates geographical distances when \code{buffer.dist=TRUE}.
#' Effects of adding geographical distances into modeling are explained in detail in \href{https://doi.org/10.7717/peerj.5518}{Hengl et al. (2018)}.
#'
#' @export
#'
#' @examples
#' library(rgdal)
#' library(geoR)
#' library(plotKML)
#' library(raster)
#' library(parallelMap)
#' library(xgboost)
#' library(kernlab)
#' library(mlr)
#' demo(meuse, echo=FALSE)
#' ## Regression:
#' m <- train.spLearner(meuse["lead"], covariates=meuse.grid[,c("dist","ffreq")], lambda = 1)
#' meuse.lead <- predict(m)
#' plot(raster(meuse.lead$pred["model"]), col=R_pal[["rainbow_75"]][4:20], main="spLearner", axes=FALSE, box=FALSE)
#' points(meuse, pch="+")
#' plot(raster(meuse.lead$pred["model.error"]), col=rev(bpy.colors()), main="Model error", axes=FALSE, box=FALSE)
#' points(meuse, pch="+")
#' ## Classification:
#' mC <- train.spLearner(meuse["soil"], covariates=meuse.grid[,c("dist","ffreq")])
#' meuse.soil <- predict(mC)
#' spplot(meuse.soil$pred[grep("prob.", names(meuse.soil$pred))], col=SAGA_pal[["SG_COLORS_YELLOW_RED"]], zlim=c(0,1))
#' spplot(meuse.soil$pred[grep("error.", names(meuse.soil$pred))], col=rev(bpy.colors()))
#' \dontrun{
#' ## SIC1997
#' data("sic1997")
#' X <- sic1997$swiss1km[c("CHELSA_rainfall","DEM")]
#' mR <- train.spLearner(sic1997$daily.rainfall, covariates=X, lambda=1)
#' rainfall1km <- predict(mR)
#' par(mfrow=c(1,2), oma=c(0,0,0,1), mar=c(0,0,4,3))
#' plot(raster(rainfall1km$pred["model"]), col=R_pal[["rainbow_75"]][4:20], main="spLearner", axes=FALSE, box=FALSE)
#' points(sic1997$daily.rainfall, pch="+")
#' plot(raster(rainfall1km$pred["model.error"]), col=rev(bpy.colors()), main="Model error", axes=FALSE, box=FALSE)
#' points(sic1997$daily.rainfall, pch="+")
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
#' mB <- train.spLearner(eberg["Parabraunerde"], covariates=X, family=binomial(), cov.model = "nugget")
#' eberg.Parabraunerde <- predict(mB)
#' plot(raster(eberg.Parabraunerde$pred["prob.1"]), col=SAGA_pal[["SG_COLORS_YELLOW_RED"]], zlim=c(0,1))
#' points(eberg["Parabraunerde"], pch="+")
#'
#' ## Factor variable:
#' data(eberg)
#' coordinates(eberg) <- ~X+Y
#' proj4string(eberg) <- CRS("+init=epsg:31467")
#' X <- eberg_grid[c("PRMGEO6","DEMSRT6","TWISRT6","TIRAST6")]
#' mF <- train.spLearner(eberg["TAXGRSC"], covariates=X, buffer.dist=FALSE)
#' TAXGRSC <- predict(mF)
#' plot(stack(TAXGRSC$pred[grep("prob.", names(TAXGRSC$pred))]), col=SAGA_pal[["SG_COLORS_YELLOW_RED"]], zlim=c(0,1))
#' }
setMethod("train.spLearner", signature(observations = "SpatialPointsDataFrame", formulaString = "ANY", covariates = "SpatialPixelsDataFrame"), function(observations, formulaString, covariates, SL.library, family = gaussian(), method = "stack.cv", predict.type, super.learner = "regr.glmnet", subsets = 3, cvControl = 5, lambda = 0.5, cov.model = "exponential", subsample = 2000, parallel = "multicore", buffer.dist = TRUE, spc = TRUE, id = NULL, weights = NULL, ...){

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
    }
    if(spc==TRUE){
      covariates <- spc(covariates)@predicted
      if(ncol(covariates)>4){
        ## remove last PC as this usually only carries noise?
        covariates <- covariates[-ncol(covariates)]
      }
    }
    if(buffer.dist==TRUE){
      message("Deriving buffer distances to points...", immediate. = TRUE)
      classes <- as.factor(1:nrow(observations))
      covariates <- sp::cbind.Spatial(covariates, buffer.dist(observations[tv], covariates[1], classes))
    }
    formulaString <- as.formula(paste(tv, " ~ ", paste(names(covariates), collapse="+")))
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
  r.sel <- complete.cases(ov[,all.vars(formulaString)])
  n.r.sel <- sum(r.sel)
  if(!n.r.sel==nrow(observations)){
    message(paste0("Subsetting observations to ", round(n.r.sel/nrow(observations), 2)*100, "% complete cases..."), immediate. = TRUE)
  }
  m <- train.spLearner.matrix(observations = ov, formulaString = formulaString, covariates = covariates, SL.library = SL.library, family = family, subsets = subsets, cvControl = cvControl, lambda = lambda, cov.model = cov.model, subsample = subsample, parallel = parallel, id=id, weights=weights, ...)
 return(m)
})

#' Overlay points and grids and prepare regression matrix
#'
#' @param observations
#' @param formulaString
#' @param covariates
#' @param dimensions
#'
#' @return
#' @export
#'
#' @examples
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
#' @param x
#' @param ...
#'
#' @return
#' @export
#'
#' @examples
"print.spLearner" <- function(x, ...){
  message("Ensemble model:")
  if(any(class(x@spModel)=="subsemble")){
    print(x@spModel$metafit$fit)
    print(paste0("CV R-square: ", (1-x@spModel$metafit$fit$object$deviance/x@spModel$metafit$fit$object$null.deviance)))
  }
  if(any(class(x@spModel)=="BaseEnsembleModel")){
    print(x@spModel)
  }
  message("Variogram model:")
  print(x@vgmModel$vgm)
  message("Total observations:")
  nrow(x@vgmModel$observations)
}

#' Predict spLearner
#'
#' @param object
#' @param predictionLocations
#' @param model.error
#' @param ...
#'
#' @return SpatialPixelsDataFrame object with predictions and model error.
#' @export
#'
#' @examples
"predict.spLearner" <- function(object, predictionLocations, model.error=TRUE, ...){
  if(missing(predictionLocations)){
    predictionLocations <- object@covariates
  }
  if(any(class(object@spModel)=="BaseEnsembleModel")){
    if(model.error==TRUE){
      out <- as.matrix(as.data.frame(getStackedBaseLearnerPredictions(object@spModel, newdata=predictionLocations@data)))
      if(any(object@spModel$task.desc$type=="classif")){
        lvs <- object@spModel$task.desc$class.levels
        pred.prob <- NULL
        for(j in lvs){
          pred.prob[[paste0("prob.",j)]] <- matrixStats::rowMeans2(out[,grep(paste0(".", j), attr(out, "dimnames")[[2]])], na.rm = TRUE)
          pred.prob[[paste0("error.",j)]] <- matrixStats::rowSds(out[,grep(paste0(".", j), attr(out, "dimnames")[[2]])], na.rm = TRUE)
        }
        pred <- SpatialPixelsDataFrame(predictionLocations@coords, data=data.frame(pred.prob), grid = predictionLocations@grid, proj4string = predictionLocations@proj4string)
      } else {
        pred <- SpatialPixelsDataFrame(predictionLocations@coords, data=data.frame(model=matrixStats::rowMeans2(out, na.rm = TRUE)), grid = predictionLocations@grid, proj4string = predictionLocations@proj4string)
        ## weighted Sds where weights are the metalearner coefficients
        #wt <- abs(object@spModel$metafit$fit$object$coefficients[-1])
        #pred$model.error <- matrixStats::rowWeightedSds(out, w=wt, na.rm=TRUE)
        pred$model.error <- matrixStats::rowSds(out, na.rm=TRUE)
      }
      return(out <- list(pred=pred, subpred=as.data.frame(out)))
    } else {
      out <- predict(object@spModel, newdata=predictionLocations@data)
      pred <- SpatialPixelsDataFrame(predictionLocations@coords, data=out$data, grid = predictionLocations@grid, proj4string = predictionLocations@proj4string)
      return(out <- list(pred=pred, subpred=NULL))
    }
  }
}

setMethod("predict", signature(object = "spLearner"), predict.spLearner)
