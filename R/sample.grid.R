#' Sample spatial points by grids
#'
#' @rdname sample.grid
#' @aliases sample.grid sample.grid,SpatialPointsDataFrame-method
#'
#' @description Get a subset of a object of class \code{"SpatialPoints"} or \code{"SpatialPointsDataFrame"} avoiding spatial clustering.
#'
#' @param obj \code{"SpatialPoints*"} object,
#' @param cell.size numeric; the cell size of the overlayed \code{"SpatialGridDataFrame"} in the form of \code{c(x,y)},
#' @param n integer; specifies maximum number points in each grid,
#' @param bbox matrix; the bounding box of output \code{"SpatialPoints"} or \code{"SpatialPointsDataFrame"}; it is set the same as the \code{obj} if missing
#' @param ... other optional arguments that can be passed to \code{over}
#'
#' @return Returns a list of two objects: (1) an object of type \code{"SpatialPoints"} or \code{"SpatialPointsDataFrame"} that contains a subset of the {obj}, and (2) resulting grid.
#'
#' @export
#'
#' @note Spatial points are overlaid with spatial grids with a specified cell size and then get a subset from each grid with a specified number at most. If one grid has less points than the specified number, all the points are taken. If one grid has more points than the specified number, only this number of points are taken by \code{\link{sample}}. This function can be used when there are too much point observations to be handled, especially for spatially clustered observations. The total number of sampled points are determined by \code{cell.size} and \code{n} together. You will get fewer the sampled points when \code{cell.size} is larger, or/and when \code{n} is smaller. Similar sample sizes can be achieved by different combinations of \code{cell.size} and \code{n}.
#'
#' @author Wei Shangguan
#'
#' @references
#' \itemize{
#'   \item Shangguan, W., Hengl, T., de Jesus, J. M., Yuan, H., & Dai, Y. (2017). \href{https://doi.org/10.1002/2016MS000686}{Mapping the global depth to bedrock for land surface modeling}. Journal of Advances in Modeling Earth Systems, 9(1), 65-88.
#' }
#'
#' @examples
#' library(sp)
#' data(edgeroi)
#' profs <- edgeroi[["sites"]]
#' coordinates(profs) <- ~  LONGDA94 + LATGDA94
#' proj4string(profs) <- CRS("+proj=longlat +ellps=GRS80 +towgs84=0,0,0,0,0,0,0 +no_defs")
#' ## sample SpatialPointsDataFrame:
#' prof1 <- sample.grid(profs, cell.size = c(0.02,0.02), n = 1)
#' l0 <- list("sp.points", profs, pch=1, col="red")
#' l1 <- list("sp.points", prof1$subset, pch="+", col="black", cex=1.2)
#' spplot(prof1$grid, scales=list(draw=TRUE),
#'    col.regions="grey", sp.layout=list(l0, l1))
#' ## Subsampling ratio in percent:
#' round(length(prof1$subset)/length(profs)*100, 1)
#' @docType methods
setMethod("sample.grid", signature(obj ="SpatialPointsDataFrame"), function(obj, cell.size, n, bbox, ...){
  #to avoid the gid in obj
  if ("gid" %in% names(obj)){
    gid <- obj$gid
  }else{
    gid <- NULL
  }
  obj$gid <- 1:length(obj)
  if(missing(bbox)) { bbox <- obj@bbox }
  if(missing(cell.size)){
    ## automatically determine width:
    cell.size <- bbox[1,]/400
    message("Assigning 'cell.size'", immediate. = TRUE)
  }
  x <- sp::GridTopology(cellcentre.offset=bbox[,1],
                        cellsize=cell.size,
                        cells.dim=c(floor(abs(diff(bbox[1,])/cell.size[1])),
                                    ncols=floor(abs(diff(bbox[2,])/cell.size[2]))))
  r.sp <- sp::SpatialGrid(x, proj4string = obj@proj4string)
  r <- raster::raster(r.sp)
  grd <- methods::as(raster::rasterize(obj, r, field = "gid"), "SpatialGridDataFrame")
  ov <- sp::over(obj, grd, ...)
  pnts <- tapply(obj$gid, ov, function(x, n){
    if(length(x) <= n){ x
    } else { sample(x,n) }
  }, n=n)
  pnts <- unlist(pnts, use.names = FALSE)
  ret <- list(obj[obj$gid %in% pnts, ], grd)
  names(ret) <- c("subset", "grid")
  if(is.null(gid)){
    ret$subset$gid <- NULL
  }else{
    ret$subset$gid <- gid
  }
  return(ret)
})

