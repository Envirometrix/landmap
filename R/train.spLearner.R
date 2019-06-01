.train.spLearner.ov <- function(observations, formulaString, covariates, SL.library, family = gaussian(), subsets = 3, cvControl = list(V=5), lambda = 0.5, cov.model = "exponential", subsample = 5000, parallel = "multicore", cell.size, id,  ...){
  if(!.Platform$OS.type=="unix") { parallel <- "seq" }
  tv <- all.vars(formulaString)[1]
  Y <- observations[,tv]
  if(missing(SL.library)){
    if(is.numeric(Y) & family$family == "gaussian"){
      SL.library <- c("SL.xgboost", "SL.ranger", "SL.ksvm")
    }
    if(is.factor(Y) | family$family == "binomial"){
      SL.library <- c("SL.ranger", "SL.ksvm", "SL.glmnet")
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
    points <- observations
    xyn <- attr(covariates@bbox, "dimnames")[[1]]
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
  if(missing(id)){
    grd <- sp::GridTopology(cellcentre.offset=covariates@bbox[,1], cellsize=rep(cell.size,2), cells.dim=c(ceiling(abs(diff(covariates@bbox[1,])/cell.size)), ncols=ceiling(abs(diff(covariates@bbox[2,])/cell.size))))
    r.sp <- sp::SpatialGridDataFrame(grd, proj4string = covariates@proj4string, data=data.frame(gid=1:(grd@cells.dim[1] * grd@cells.dim[2])))
    id <- sp::over(sp::SpatialPoints(observations[,attr(covariates@bbox, "dimnames")[[1]]], proj4string = covariates@proj4string), r.sp)$gid
    message("Using block size ID for spatial Cross Validation...", immediate. = TRUE)
  }
  ## fit the model:
  r.sel <- complete.cases(observations[, all.vars(formulaString)])
  if(is.factor(Y)){
    ## subsemble does not work with Factor vars
    message("Fitting a spatial learner using classification models...", immediate. = TRUE)
    m1 <- ranger::ranger(formulaString, observations[which(r.sel),], importance="impurity", write.forest=TRUE, probability=TRUE)
    m2 <- e1071::svm(formulaString, observations[which(r.sel),], probability=TRUE, cross=5)
    m3 <- nnet::multinom(formulaString, observations[which(r.sel),])
    m <- list(m1, m2, m3)
  } else {
    message("Fitting a spatial learner using 'subsemble'...", immediate. = TRUE)
    m <- subsemble::subsemble(x = observations[which(r.sel), all.vars(formulaString)[-1]], y = Y[which(r.sel)], learner = SL.library, parallel = parallel, id = id[which(r.sel)], cvControl = cvControl, subsets = subsets, family = family, ...)
  }
  out <- new("spLearner", spModel = m, vgmModel = rvgm, covariates = covariates, spID = r.sp)
  return(out)
}

setMethod("train.spLearner", signature(observations = "data.frame", formulaString = "formula", covariates = "SpatialPixelsDataFrame"), .train.spLearner.ov)

#' Train a spatial prediction and/or interpolation model using Ensemble Machine Learning
#'
#' @description Automated spatial predictions and/or interpolation using Ensemble Machine Learning. Extends functionality of the \href{https://github.com/ledell/subsemble}{subsemble} and the \href{https://github.com/ecpolley/SuperLearner}{SuperLearner} packages. Suitable for predicting numeric, binomial and factor-type variables.
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
#' In the case of factor variables, prediction are based on simple average from
#' \code{ranger}, \code{e1071::svm} and \code{nnet::multinom}, which might be suboptimal.
#'
#' @export
#'
#' @examples
#' library(rgdal)
#' library(geoR)
#' library(plotKML)
#' library(raster)
#' library(SuperLearner)
#' library(subsemble)
#' demo(meuse, echo=FALSE)
#' m <- train.spLearner(meuse["lead"], covariates=meuse.grid[,c("dist","ffreq")], lambda = 1, cov.model = "nugget")
#' meuse.lead <- predict(m)
#' plot(raster(meuse.lead$pred["model"]), col=R_pal[["rainbow_75"]][4:20], main="spLearner", axes=FALSE, box=FALSE)
#' points(meuse, pch="+")
#' plot(raster(meuse.lead$pred["model.error"]), col=rev(bpy.colors()), main="Model error", axes=FALSE, box=FALSE)
#' points(meuse, pch="+")
#' \dontrun{
#' ## SIC1997
#' data("sic1997")
#' mR <- train.spLearner(sic1997$daily.rainfall, covariates=sic1997$swiss1km[c("CHELSA_rainfall","DEM")], lambda=1)
#' rainfall1km <- predict(mR)
#' par(mfrow=c(1,2), oma=c(0,0,0,1), mar=c(0,0,4,3))
#' plot(raster(rainfall1km$pred["model"]), col=R_pal[["rainbow_75"]][4:20], main="spLearner", axes=FALSE, box=FALSE)
#' points(sic1997$daily.rainfall, pch="+")
#' plot(raster(rainfall1km$pred["model.error"]), col=rev(bpy.colors()), main="Model error", axes=FALSE, box=FALSE)
#' points(sic1997$daily.rainfall, pch="+")
#'
#' data(eberg)
#' eb.s <- sample.int(nrow(eberg), 1400)
#' eberg <- eberg[eb.s,]
#' coordinates(eberg) <- ~X+Y
#' proj4string(eberg) <- CRS("+init=epsg:31467")
#' ## Binomial variable
#' summary(eberg$TAXGRSC)
#' eberg$Parabraunerde <- ifelse(eberg$TAXGRSC=="Parabraunerde", 1, 0)
#' data(eberg_grid)
#' gridded(eberg_grid) <- ~x+y
#' proj4string(eberg_grid) <- CRS("+init=epsg:31467")
#' mB <- train.spLearner(eberg["Parabraunerde"], covariates=eberg_grid[c("PRMGEO6","DEMSRT6","TWISRT6","TIRAST6")],
#'                         family=binomial(), cov.model = "nugget", SL.library = c("SL.ranger", "SL.ksvm"))
#' ## Factor variable
#' data(eberg)
#' coordinates(eberg) <- ~X+Y
#' proj4string(eberg) <- CRS("+init=epsg:31467")
#' mF <- train.spLearner(eberg["TAXGRSC"], covariates=eberg_grid[c("PRMGEO6","DEMSRT6","TWISRT6","TIRAST6")])
#' TAXGRSC <- predict(mF)
#' plot(stack(TAXGRSC$pred), col=SAGA_pal[["SG_COLORS_YELLOW_RED"]], zlim=c(0,1))
#' }
setMethod("train.spLearner", signature(observations = "SpatialPointsDataFrame", formulaString = "ANY", covariates = "SpatialPixelsDataFrame"), function(observations, formulaString, covariates, SL.library, family = gaussian(), subsets = 3, cvControl = list(V=5), lambda = 0.5, cov.model = "exponential", subsample = 2000, parallel = "multicore", buffer.dist = TRUE, spc = TRUE,  ...){

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
      if(is.factor(observations@data[,1])){
        classes <- observations@data[,1]
      } else {
        classes <- as.factor(1:nrow(observations))
      }
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
  m <- .train.spLearner.ov(observations = ov, formulaString = formulaString, covariates = covariates, SL.library = SL.library, family = family, subsets = subsets, cvControl = cvControl, lambda = lambda, cov.model = cov.model, subsample = subsample, parallel = parallel, ...)
 return(m)
})

## overlay points and grids and prepare regression matrix
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


"print.spLearner" <- function(x, ...){
  message("Ensemble model:")
  if(class(x@spModel)=="subsemble"){
    print(x@spModel$metafit$fit)
    print(paste0("CV R-square: ", (1-x@spModel$metafit$fit$object$deviance/x@spModel$metafit$fit$object$null.deviance)))
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
  if(class(object@spModel)=="subsemble"){
    out <- predict(object@spModel, predictionLocations@data)
    pred <- SpatialPixelsDataFrame(predictionLocations@coords, data=data.frame(model=out$pred), grid = predictionLocations@grid, proj4string = predictionLocations@proj4string)
    if(model.error==TRUE){
      if(class(object@spModel)=="subsemble"){
        ## weighted Sds where weights are the metalearner coefficients
        wt <- abs(object@spModel$metafit$fit$object$coefficients[-1])
        pred$model.error <- matrixStats::rowWeightedSds(as.matrix(out$subpred), w=wt, na.rm=TRUE)
      }
    }
    return(out = list(pred=pred, subpred=out$subpred))
  } else {
    leg <- object@spModel[[1]]$forest$levels
    out1 <- predict(object@spModel[[1]], predictionLocations@data, probability=TRUE, na.action = na.pass)$predictions
    out2 <- attr(predict(object@spModel[[2]], predictionLocations@data, probability=TRUE, na.action = na.pass), "probabilities")
    out3 <- predict(object@spModel[[3]], predictionLocations@data, type="probs", na.action = na.pass)
    lt <- list(out1[,leg], out2[,leg], out3[,leg])
    out <- Reduce("+", lt) / length(lt)
    pred <- SpatialPixelsDataFrame(predictionLocations@coords, data=data.frame(out), grid = predictionLocations@grid, proj4string = predictionLocations@proj4string)
    return(out = list(pred=pred, subpred=lt))
  }
}

setMethod("predict", signature(object = "spLearner"), predict.spLearner)
