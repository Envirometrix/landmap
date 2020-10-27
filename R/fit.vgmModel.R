#' Fit variogram using point data
#'
#' @aliases fit.vgmModel
#'
#' @param formulaString.vgm formula.
#' @param rmatrix data.frame with coordinates and values of covariates.
#' @param predictionDomain SpatialPixelsDataFrame.
#' @param cov.model covariance model type used by the geoR package.
#' @param dimensions optional 2D or 3D dimensions.
#' @param lambda transformation value used by the geoR package.
#' @param psiR range parameter used by the geoR package.
#' @param subsample number of subset of original samples.
#' @param ini.var initial variance (sill) used by the geoR package.
#' @param ini.range initial range parameter used by the geoR package.
#' @param fix.psiA setting used by the geoR package.
#' @param fix.psiR setting used by the geoR package.
#' @param ... optional arguments to pass to the geoR package.
#'
#' @return Fitted variogram model
#' @export
#'
#' @author \href{https://opengeohub.org/people/tom-hengl}{Tom Hengl}
#'
#' @note Extends variogram fitting functionality from the geoR package.
#' Can be used for 2D or 3D point data sets, with and without trend variables.
#' Models need to be in the form \code{zinc ~ dist} and only numeric variables are allowed.
#' Often reports \code{Singular matrix. Covariates may have different orders of magnitude.} if the covariates are perfectly aligned.
#'
#' @examples
#' library(raster)
#' library(rgdal)
#' library(geoR)
#' demo(meuse, echo=FALSE)
#' vgm = fit.vgmModel(zinc~dist, as.data.frame(meuse), meuse.grid["dist"], lambda=1)
#' plot(variog(vgm$geodata))
#' lines(vgm$vgm)
setMethod("fit.vgmModel", signature(formulaString.vgm = "formula", rmatrix = "data.frame", predictionDomain = "SpatialPixelsDataFrame"), function(formulaString.vgm, rmatrix, predictionDomain, cov.model = "exponential", dimensions = list("2D", "3D", "2D+T", "3D+T"), lambda = 0.5, psiR = NULL, subsample = nrow(rmatrix), ini.var, ini.range, fix.psiA = FALSE, fix.psiR = FALSE, ...){

  if(missing(dimensions)){ dimensions <- dimensions[[1]] }
  if(is.na(sp::proj4string(predictionDomain))){ stop("proj4 string required for argument 'predictionDomain'") }
  if(!any(names(rmatrix) %in% all.vars(formulaString.vgm))){
    stop("Variables in the 'formulaString.vgm' not found in the 'rmatrix' object.")
  }
  ## variable names:
  tv <- all.vars(formulaString.vgm)[1]
  if(!is.numeric(rmatrix[,tv])){ stop("Numeric variable expected for 'fit.vgmModel'") }
  if(length(all.vars(formulaString.vgm))>1){
    if(length(all.vars(formulaString.vgm))==2){
      tcovs <- stats::as.formula(paste(" ~ ", all.vars(formulaString.vgm)[2]))
    } else {
      tcovs <- stats::as.formula(paste(" ~ ", paste(all.vars(formulaString.vgm)[-1], collapse = "+")))
    }
  }
  sel.r <-  stats::complete.cases(lapply(all.vars(formulaString.vgm), function(x){rmatrix[,x]}))
  if(!sum(sel.r)==nrow(rmatrix)){ rmatrix <- rmatrix[sel.r,] }
  ## spatial coordinates (column names):
  xyn <- attr(predictionDomain@bbox, "dimnames")[[1]]
  if(!any(names(rmatrix) %in% xyn)){
       stop(paste("Column names:", paste(xyn[which(!(xyn %in% names(rmatrix)))], collapse=", "), "could not be located in the regression matrix"))
  }
  ## add 3D dimension if missing:
  if(dimensions=="3D" & length(xyn)==2){
    xyn <- c(xyn, "altitude")
  }
  ## create spatial points:
  sp::coordinates(rmatrix) <- stats::as.formula(paste("~", paste(xyn, collapse = "+"), sep=""))
  sp::proj4string(rmatrix) <- predictionDomain@proj4string
  points <- rmatrix
  ## subset to speed up the computing:
  if(subsample < nrow(rmatrix)){
    pcnt <- subsample/nrow(rmatrix)
    message(paste0("Subsetting observations to ", signif(pcnt*100, 1), "%..."))
    rmatrix <- rmatrix[stats::runif(nrow(rmatrix))<pcnt,]
  }
  if(length(all.vars(formulaString.vgm))>1){
    x.geo <- geoR::as.geodata(rmatrix[c(tv, all.vars(formulaString.vgm)[-1])], data.col=tv, covar.col=all.vars(formulaString.vgm)[-1])
  } else {
    x.geo <- geoR::as.geodata(rmatrix, data.col=tv)
  }
  if(!cov.model == "nugget"){
    ## guess the dimensions:
    if(missing(dimensions)){
        xyn = attr(rmatrix@bbox, "dimnames")[[1]]
        if(length(xyn)==2) {
           dimensions = "2D"
        } else {
           dimensions = "3D"
        }
    }
    if(dimensions == "3D"){
      ## estimate area extent:
      ini.range = sqrt(sp::areaSpatialGrid(predictionDomain))/3
      ## estimate anisotropy:
      if(is.null(psiR)){
        ## estimate initial range in the vertical direction:
        dr <- abs(diff(range(rmatrix@coords[,3], na.rm=TRUE)))/3
        psiR <- 2*dr/ini.range
      }
    }
    if(dimensions == "2D"){
      diag <- (sqrt((rmatrix@bbox[1,2]-rmatrix@bbox[1,1])**2+(rmatrix@bbox[2,2]-rmatrix@bbox[2,1])**2))/3
      ## check if it is projected object:
      if(!is.na(sp::proj4string(predictionDomain))){
        if(!sp::is.projected(predictionDomain)){
          if(requireNamespace("fossil", quietly = TRUE)){
            ## Haversine Formula for Great Circle distance
            p.1 <- matrix(c(predictionDomain@bbox[1,1], predictionDomain@bbox[1,2]), ncol=2, dimnames=list(1,c("lon","lat")))
            p.2 <- matrix(c(predictionDomain@bbox[2,1], predictionDomain@bbox[2,2]), ncol=2, dimnames=list(1,c("lon","lat")))
            ini.range = fossil::deg.dist(lat1=p.1[,2], long1=p.1[,1], lat2=p.2[,2], long2=p.2[,1])/2
          }
        } else {
          ini.range <- sqrt(sp::areaSpatialGrid(predictionDomain))/3
        }
      } else {
        ini.range = diag/3
      }
    }
    ## initial variogram:
    if(lambda==1){
      ini.var <- stats::var(log1p(x.geo$data), na.rm = TRUE)
    } else {
      ini.var <- stats::var(x.geo$data, na.rm = TRUE)
    }
    ## fit sample variogram:
    rvgm <- list(cov.model="nugget", lambda=lambda, practicalRange=ini.range)
    if(length(all.vars(formulaString.vgm))==1){
      message("Fitting a variogram using 'linkfit'...", immediate. = TRUE)
      try( rvgm <- geoR::likfit(x.geo, lambda = lambda, messages = FALSE, ini = c(ini.stats::var, ini.range), cov.model = cov.model), silent = TRUE )
    } else {
      message("Fitting a variogram using 'linkfit' and trend model...", immediate. = TRUE)
      try( rvgm <- geoR::likfit(x.geo, lambda = lambda, messages = FALSE, trend = tcovs, ini = c(ini.var, ini.range), fix.psiA = FALSE, fix.psiR = FALSE, cov.model = cov.model, lik.met = "REML"), silent = TRUE )
    }
    if(class(.Last.value)[1]=="try-error"){
      warning("Variogram model could not be fitted.")
    }
  } else {
    rvgm <- list(cov.model="nugget", lambda=lambda, practicalRange=NA)
  }
  return(list(vgm=rvgm, observations=points, geodata=x.geo))
})
