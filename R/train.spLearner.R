#' Train a spatial prediction and/or interpolation model using Ensemble Machine Learning
#' from a regression/classification matrix
#'
#' @param observations Data frame regression matrix,
#' @param formulaString Model formula,
#' @param covariates SpatialPixelsDataFrame object,
#' @param SL.library List of learners,
#' @param family Family e.g. gaussian(),
#' @param method Ensemble stacking method (see makeStackedLearner),
#' @param predict.type Prediction type 'prob' or 'response',
#' @param super.learner Ensemble stacking model usually \code{regr.lm},
#' @param subsets Number of subsets for repeated CV,
#' @param lambda Target variable transformation (0.5 or 1),
#' @param cov.model Covariance model for variogram fitting,
#' @param subsample For large datasets consider random subsetting training data,
#' @param parallel Initiate parellel processing,
#' @param cell.size Block size for spatial Cross-validation,
#' @param id Id column name to control clusters of data,
#' @param weights Optional weights (per row) that learners will use to account for variable data quality,
#' @param ... other arguments that can be passed on to \code{mlr::makeStackedLearner},
#'
#' @return Object of class \code{spLearner}
#' @export
#'
#' @author \href{https://opengeohub.org/people/tom-hengl}{Tom Hengl}
train.spLearner.matrix <- function(observations, formulaString, covariates, SL.library, family = gaussian(), method = "stack.cv", predict.type, super.learner, subsets = 5, lambda = 0.5, cov.model = "exponential", subsample = 10000, parallel = "multicore", cell.size, id = NULL, weights = NULL, ...){
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
        SL.library <- c("regr.ranger", "regr.ksvm", "regr.glmnet", "regr.cubist")
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
    ## variogram fitting:
    if(lambda==1){
      ini.var <- var(log1p(Y), na.rm = TRUE)
    }
    if(lambda==0.5){
      ini.var <- var(Y, na.rm = TRUE)
    }
    ## strip buffer distances from vgm modeling
    rm.n = unlist(sapply(c("rX_*$","rY_*$","layer.*$"), function(i){grep(pattern=glob2rx(i), names(covariates))}))
    if(length(rm.n)>0){
      covs.vgm <- names(covariates)[-rm.n]
    } else {
      covs.vgm <- names(covariates)
    }
    formulaString.vgm <- as.formula(paste(tv, "~", paste(covs.vgm, collapse="+")))
    rvgm <- fit.vgmModel(formulaString.vgm, rmatrix = observations, predictionDomain = covariates[covs.vgm], lambda = lambda, ini.var = ini.var, cov.model = cov.model, subsample = subsample)
  } else {
    message("Skipping variogram modeling...", immediate. = TRUE)
    if(missing(super.learner)){ super.learner <- "classif.glmnet" }
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
    id <- as.factor(make.names(id, unique=TRUE))
    message("Estimating block size ID for spatial Cross Validation...", immediate. = TRUE)
  }
  r.sel <- complete.cases(observations[, all.vars(formulaString)])
  ## fit the mlr model:
  mlr::configureMlr()
  if(parallel=="multicore"){
    parallelMap::parallelStartSocket(parallel::detectCores())
  }
  message(paste0("Using learners: ", paste(SL.library, collapse = ", "), "..."), immediate. = TRUE)
  if(is.factor(Y) | family$family == "binomial"){
    message("Fitting a spatial learner using 'mlr::makeClassifTask'...", immediate. = TRUE)
    if(missing(predict.type)){ predict.type = "prob" }
    if(is.null(weights)){
      tsk <- mlr::makeClassifTask(data = observations[which(r.sel),all.vars(formulaString)], target = tv, coordinates = observations[which(r.sel),xyn])
    } else {
      tsk <- mlr::makeClassifTask(data = observations[which(r.sel),all.vars(formulaString)], target = tv, weights = weights[which(r.sel)], coordinates = observations[which(r.sel),xyn])
    }
    lrns <- lapply(SL.library, mlr::makeLearner)
    lrns <- lapply(lrns, mlr::setPredictType, "prob")
    init.m <- mlr::makeStackedLearner(base.learners = lrns, predict.type = predict.type, method = method, super.learner = super.learner, ...)
  } else {
    message("Fitting a spatial learner using 'mlr::makeRegrTask'...", immediate. = TRUE)
    if(missing(predict.type)){ predict.type = "response" }
    if(is.null(weights)){
      tsk <- mlr::makeRegrTask(data = observations[which(r.sel),all.vars(formulaString)], target = tv, coordinates = observations[which(r.sel),xyn], blocking = id[which(r.sel)])
    } else {
      tsk <- mlr::makeRegrTask(data = observations[which(r.sel),all.vars(formulaString)], target = tv, weights = weights[which(r.sel)], coordinates = observations[which(r.sel),xyn], blocking = id[which(r.sel)])
    }
    lrns <- lapply(SL.library, mlr::makeLearner)
    init.m <- mlr::makeStackedLearner(base.learners = lrns, predict.type = predict.type, method = method, super.learner = super.learner, ...)
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
#' @param SL.library List of learners,
#' @param family Family e.g. gaussian(),
#' @param method Ensemble stacking method (see makeStackedLearner) usually \code{stack.cv},
#' @param predict.type Prediction type 'prob' or 'response',
#' @param super.learner Ensemble stacking model usually \code{regr.lm},
#' @param subsets Number of subsets for repeated CV,
#' @param lambda Target variable transformation (0.5 or 1),
#' @param cov.model Covariance model for variogram fitting,
#' @param subsample For large datasets consider random subsetting training data,
#' @param parallel Initiate parellel processing,
#' @param buffer.dist Specify whether to use buffer distances to points as covariates,
#' @param oblique.coords Specify whether to use oblique coordinates as covariates,
#' @param theta.list List of angles (in radians) used to derive oblique coordinates,
#' @param id Id column name to control clusters of data,
#' @param weights Optional weights (per row) that learners will use to account for variable data quality,
#' @param ... other arguments that can be passed on to \code{mlr::makeStackedLearner},
#'
#' @importClassesFrom sp SpatialPixelsDataFrame SpatialPointsDataFrame
#'
#' @return object of class spLearner, which contains fitted model, variogram model and spatial grid
#' used for Cross-validation.
#'
#' @author \href{https://opengeohub.org/people/tom-hengl}{Tom Hengl}
#'
#' @note By default uses oblique coordinates (rotated coordinates) as described in Moller et al. (2019) "Oblique Coordinates as Covariates for Digital Soil Mapping" to account for geographical distribution of values.
#' Buffer geographical distances can be added by setting \code{buffer.dist=TRUE}.
#' Using either oblique coordinates and/or buffer distances is not recommended for point data set with distinct spatial clustering.
#' Effects of adding geographical distances into modeling are explained in detail in \href{https://doi.org/10.7717/peerj.5518}{Hengl et al. (2018)}.
#' Default learners used for regression are \code{c("regr.ranger", "regr.ksvm", "regr.glmnet", "regr.cubist")}.
#' Default learners used for classification / binomial variables are \code{c("classif.ranger", "classif.svm", "classif.multinom")}, with \code{predict.type="prob"}.
#' When using \code{method = "stack.cv"} each training and prediction round could produce somewhat different results due to randomisation of CV.
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
#' library(ranger)
#' library(deepnet)
#' library(Cubist)
#' demo(meuse, echo=FALSE)
#'
#' ## Regression:
#' m <- train.spLearner(meuse["zinc"], covariates=meuse.grid[,c("dist","ffreq")], lambda = 1)
#' ## Ensemble model (meta-learner):
#' summary(m@spModel$learner.model$super.model$learner.model)
#' meuse.y <- predict(m)
#' plot(raster(meuse.y$pred["response"]), col=R_pal[["rainbow_75"]][4:20],
#'    main="spLearner", axes=FALSE, box=FALSE)
#' points(meuse, pch="+")
#' plot(raster(meuse.y$pred["model.error"]), col=rev(bpy.colors()),
#'    main="Model error", axes=FALSE, box=FALSE)
#' points(meuse, pch="+")
#'
#' ## Classification:
#' SL.library <- c("classif.ranger", "classif.xgboost", "classif.nnTrain")
  #' mC <- train.spLearner(meuse["soil"], covariates=meuse.grid[,c("dist","ffreq")],
#'    SL.library = SL.library, super.learner = "classif.glmnet")
#' meuse.soil <- predict(mC)
#' spplot(meuse.soil$pred[grep("prob.", names(meuse.soil$pred))], col.regions=SAGA_pal[["SG_COLORS_YELLOW_RED"]], zlim=c(0,1))
#' spplot(meuse.soil$pred[grep("error.", names(meuse.soil$pred))], col.regions=rev(bpy.colors()))
#'
#' \dontrun{
#' ## SIC1997
#' data("sic1997")
#' X <- sic1997$swiss1km[c("CHELSA_rainfall","DEM")]
#' mR <- train.spLearner(sic1997$daily.rainfall, covariates=X, lambda=1)
#' rainfall1km <- predict(mR)
#' par(mfrow=c(1,2), oma=c(0,0,0,1), mar=c(0,0,4,3))
#' plot(raster(rainfall1km$pred["response"]), col=R_pal[["rainbow_75"]][4:20], main="spLearner", axes=FALSE, box=FALSE)
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
#' mF <- train.spLearner(eberg["TAXGRSC"], covariates=X)
#' TAXGRSC <- predict(mF)
#' plot(stack(TAXGRSC$pred[grep("prob.", names(TAXGRSC$pred))]), col=SAGA_pal[["SG_COLORS_YELLOW_RED"]], zlim=c(0,1))
#' }
#' @export
setMethod("train.spLearner", signature(observations = "SpatialPointsDataFrame", formulaString = "ANY", covariates = "SpatialPixelsDataFrame"), function(observations, formulaString, covariates, SL.library, family = gaussian(), method = "stack.cv", predict.type, super.learner = "regr.lm", subsets = 5, lambda = 0.5, cov.model = "exponential", subsample = 10000, parallel = "multicore", buffer.dist = FALSE, oblique.coords = TRUE, theta.list=seq(0, 180, length.out = 14)*pi/180, spc = TRUE, id = NULL, weights = NULL, ...){
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
      oblique.xy <- SpatialPixelsDataFrame(covariates@coords, data=as.data.frame(do.call(cbind, oblique.xy)), proj4string = covariates@proj4string)
      R <- expand.grid(c("rX","rY"), round(theta.list, 1))
      names(oblique.xy) <- paste(R$Var1, R$Var2, sep="_")
      covariates <- sp::cbind.Spatial(covariates, oblique.xy)
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
  m <- train.spLearner.matrix(observations = ov, formulaString = formulaString, covariates = covariates, SL.library = SL.library, family = family, subsets = subsets, lambda = lambda, cov.model = cov.model, subsample = subsample, parallel = parallel, id=id, weights=weights, predict.type = predict.type, ...)
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
#' @param object of type \code{spLearner}
#' @param ... optional parameters
#'
#' @return
#' @export
"print.spLearner" <- function(object){
  message("Ensemble model:")
  if(any(class(object@spModel)=="subsemble")){
    print(object@spModel$metafit$fit)
    message(paste0("CV R-square:", (1-object@spModel$metafit$fit$object$deviance/object@spModel$metafit$fit$object$null.deviance)))
  }
  if(any(class(object@spModel)=="BaseEnsembleModel")){
    print(object@spModel$learner.model$super.model$learner.model)
      message(paste0("CV R-square: ", round(1-object@spModel$learner.model$super.model$learner.model$deviance/object@spModel$learner.model$super.model$learner.model$null.deviance, 3)))
  }
  message("Variogram model:")
  print(object@vgmModel$vgm)
  message(paste0("Total observations: ", length(object@vgmModel$observations)))
}

#' Predict spLearner
#'
#' @param object of type \code{spLearner}
#' @param predictionLocations \code{SpatialPixelsDataFrame} with all variables
#' @param model.error Logical speficy if prediction errors should be derived
#' @param error.type Specify how should be the prediction error be derived
#' @param ... optional parameters
#'
#' @return SpatialPixelsDataFrame object with predictions and model error.
#' @export
"predict.spLearner" <- function(object, predictionLocations, model.error=TRUE, error.type=c("weighted.sd", "interval")[1], w, ...){
  if(any(object@spModel$task.desc$type=="classif")){
    error.type <- "weighted.sd"
  }
  if(missing(predictionLocations)){
    predictionLocations <- object@covariates
  }
  out <- predict(object@spModel, newdata=predictionLocations@data)
  if(any(class(object@spModel)=="BaseEnsembleModel")){
    if(model.error==TRUE){
      message("Predicting values using 'getStackedBaseLearnerPredictions'...", immediate. = TRUE)
      out.c <- as.matrix(as.data.frame(mlr::getStackedBaseLearnerPredictions(object@spModel, newdata=predictionLocations@data)))
      if(any(object@spModel$task.desc$type=="classif")){
        lvs <- object@spModel$task.desc$class.levels
        ## estimate weights:
        pred.prob <- NULL
        if(any(class(object@spModel$learner.model$super.model$learner.model) %in% "glmnet") & length(lvs)>2){
          w.v = .glmnet.varImp(object@spModel$learner.model$super.model$learner.model)
          ## model error
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
        pred <- SpatialPixelsDataFrame(predictionLocations@coords, data=cbind(out$data, data.frame(pred.prob)), grid = predictionLocations@grid, proj4string = predictionLocations@proj4string)
      } else {
        pred <- SpatialPixelsDataFrame(predictionLocations@coords, data=out$data, grid=predictionLocations@grid, proj4string=predictionLocations@proj4string)
        if(error.type=="interval"&object@spModel$learner.model$super.model$learner$id=="regr.lm"){
          message("Deriving prediction errors...", immediate. = TRUE)
          pred.int = predict(object@spModel$learner.model$super.model$learner.model, newdata = data.frame(out.c), interval = "prediction", level=2/3)
          ## assumes normal distribution
          pred$model.error <- (pred.int[,"upr"]-pred.int[,"lwr"])/2
        } else {
          message("Deriving model errors using weighted sd...", immediate. = TRUE)
          ## weighted Sds where weights are the metalearner coefficients
          wt <- abs(object@spModel$learner.model$super.model$learner.model$coefficients[-1])
          pred$model.error <- matrixStats::rowWeightedSds(out.c, w=wt, na.rm=TRUE)
        }
      }
    return(out <- list(pred=pred, subpred=as.data.frame(out.c)))
    } else {
      pred <- SpatialPixelsDataFrame(predictionLocations@coords, data=out$data, grid=predictionLocations@grid, proj4string = predictionLocations@proj4string)
      return(out <- list(pred=pred, subpred=NULL))
    }
  }
}

## https://stackoverflow.com/questions/35461839/glmnet-variable-importance
.glmnet.varImp <- function(x, ...) {
  coefList <- glmnet::coef.glmnet(x)
  coefList <- do.call("rbind", lapply(coefList, function(X){data.frame(X@Dimnames[[1]][X@i+1],X@x)}))
  names(coefList) <- c('var','val')
  coefList$aval = abs(coefList$val)
  out = aggregate(coefList$aval, by=list(coefList$var), FUN="mean", na.action=na.rm)
  out = out[!out$Group.1=="",]
  out
}

setMethod("show", signature(object = "spLearner"), print.spLearner)
setMethod("predict", signature(object = "spLearner"), predict.spLearner)
