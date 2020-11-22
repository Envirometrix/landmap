#' Derive buffer distances for a list of points
#'
#' @rdname buffer.dist-methods
#' @aliases buffer.dist buffer.dist,SpatialPointsDataFrame,SpatialPixelsDataFrame-method
#'
#' @description Derive buffer distances using the \code{raster::distance} function, so that these can be used as predictors for spatial prediction i.e. to account for spatial proximity to low, medium and high values.
#'
#' @param observations SpatialPointsDataFrame.
#' @param predictionDomain SpatialPixelsDataFrame.
#' @param classes vector of selected points as factors.
#' @param width maximum width for buffer distance.
#' @param parallel optional parallelization setting.
#' @param ... optional arguments to pass to \code{raster::distance} function.
#'
#' @return object of class \code{SpatialPixelsDataFrame} with distances to points
#' @export
#'
#' @author \href{https://opengeohub.org/people/tom-hengl}{Tom Hengl}
#'
#' @references
#' \itemize{
#'   \item Hengl, T., Nussbaum, M., Wright, M. N., Heuvelink, G. B., and Gräler, B. (2018) \href{https://doi.org/10.7717/peerj.5518}{Random Forest as a generic framework for predictive modeling of spatial and spatio-temporal variables}. PeerJ 6:e5518.
#' }
#'
#' @examples
#' \donttest{
#' library(raster)
#' library(rgdal)
#' demo(meuse, echo=FALSE)
#' b <- buffer.dist(meuse["zinc"], meuse.grid[1],
#'         classes=as.factor(1:nrow(meuse)), parallel=FALSE)
#' }
setMethod("buffer.dist", signature(observations = "SpatialPointsDataFrame", predictionDomain = "SpatialPixelsDataFrame"), function(observations, predictionDomain, classes, width, parallel=TRUE, ...){

  ## check classes
  if(missing(width)){ width <- sqrt(sp::areaSpatialGrid(predictionDomain)) }
  if(!length(classes)==length(observations)){ stop("Length of 'observations' and 'classes' does not match.") }
  ## remove classes without any points:
  xg <- summary(classes, maxsum=length(levels(classes)))
  selg.levs = attr(xg, "names")[xg > 0]
  if(length(selg.levs)<length(levels(classes))){
    fclasses <- as.factor(classes)
    fclasses[which(!fclasses %in% selg.levs)] <- NA
    classes <- droplevels(fclasses)
  }
  ## subset to points within the area
  r.sel <- stats::complete.cases(sp::over(observations, predictionDomain))
  if(sum(r.sel)<nrow(observations)){
    observations <- observations[r.sel,]
    classes <- classes[r.sel]
    classes <- droplevels(classes)
  }
  ## derive buffer distances
  if(parallel==TRUE & .Platform$OS.type=="unix"){
    s <- parallel::mclapply(1:length(levels(classes)), function(i){ raster::distance(raster::rasterize(observations[which(classes==levels(classes)[i]),1]@coords, y=raster::raster(predictionDomain)), width=width, ...) }, mc.cores = parallel::detectCores())
  } else {
    s <- list(NULL)
    for(i in 1:length(levels(classes))){
      s[[i]] <- raster::distance(raster::rasterize(observations[which(classes==levels(classes)[i]),1]@coords, y=raster::raster(predictionDomain)), width=width, ...)
    }
  }

  ## remove empty slots
  s <- s[sapply(s, function(x){!is.null(x)})]
  s <- raster::brick(s)
  s <- methods::as(s, "SpatialPixelsDataFrame")
  s <- s[predictionDomain@grid.index,]
  return(s)

})
