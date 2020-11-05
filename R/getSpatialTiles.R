#' Split a Spatial object into tiles
#'
#' @aliases getSpatialTiles
#' @rdname getSpatialTiles-methods
#'
#' @param obj output of the GDALinfo.
#' @param block.x size of the block in x dimension.
#' @param block.y size of the block in y dimension.
#' @param overlap.percent optional overlap percent between tiles.
#' @param limit.bbox optional bounding box.
#' @param return.SpatialPolygons logical specificies whether to return a data frame or Spatial Polygons.
#'
#' @return List of object result of clipping
#' @export
#'
#' @author \href{https://opengeohub.org/people/tom-hengl}{Tom Hengl}
#'
#' @examples
#' library(sp)
#' data(meuse.grid)
#' gridded(meuse.grid) <- ~x+y
#' tl <- getSpatialTiles(meuse.grid, block.x=1000)
#' image(meuse.grid)
#' lines(as(tl, "SpatialLines"))
#' ## all at once:
#' pix.lst <- tile(meuse.grid, block.x=1000)
#' \donttest{
#'   library(plotKML)
#'   ## raster files via rgdal:
#'   library(rgdal)
#'   fn = system.file("pictures/SP27GTIF.TIF",
#'                    package = "rgdal")
#'   obj <- GDALinfo(fn)
#'   ras.lst <- getSpatialTiles(obj, block.x=1000)
#'   offset <- c(ras.lst$offset.y[1], ras.lst$offset.x[1])
#'   region.dim <- c(ras.lst$region.dim.y[1],
#'                   ras.lst$region.dim.x[1])
#'   ## read the first tile:
#'   SP27GTIF_T1 <- readGDAL(fn, offset=offset,
#'                           region.dim=region.dim)
#'   str(SP27GTIF_T1)
#' }
setMethod("getSpatialTiles", signature(obj = "Spatial"), function(obj, block.x, block.y = block.x, overlap.percent = 0, limit.bbox = TRUE, return.SpatialPolygons = TRUE){

  if(overlap.percent<0){
    stop("'overlap.percent' argument must be a positive number")
  }

  ## check the input bbox:
  if(!(ncol(obj@bbox)==2&nrow(obj@bbox)==2&obj@bbox[1,1]<obj@bbox[1,2]&obj@bbox[2,1]<obj@bbox[2,2])){
    stop("Bounding box with two-column matrix required; the first column should contain the minimum, the second the maximum values;\n rows represent the spatial dimensions required")
  }

  bb <- obj@bbox
  btiles <- makeTiles(bb, block.x, block.y, overlap.percent, limit.bbox)

  if(return.SpatialPolygons == TRUE){
    pol <- .tiles2pol(bb=bb, btiles=btiles, proj4string=obj@proj4string)
  } else {
    pol = btiles
  }

  message(paste("Returning a list of tiles for an object of class", class(obj), "with", signif(overlap.percent, 3), "percent overlap"))
  return(pol)

})

#' Estimate a tiling system
#'
#' @param obj ANY.
#' @param block.x size of the block in x dimension.
#' @param block.y size of the block in y dimension.
#' @param overlap.percent optional overlap percent between tiles.
#' @param limit.bbox optional bounding box.
#' @param return.SpatialPolygons logical specificies whether to return a data frame or Spatial Polygons.
#'
#' @return Tiling system for a spatial object
#' @export
setMethod("getSpatialTiles", signature(obj = "ANY"), function(obj, block.x, block.y = block.x, overlap.percent = 0, limit.bbox = TRUE, return.SpatialPolygons = FALSE){

  if(!class(obj)=="GDALobj"){
    stop("Object of class \"GDALobj\" required.")
  }

  if(overlap.percent<0){
    stop("'overlap.percent' argument must be a positive number")
  }

  ## create bbox:
  bb <- matrix(c(obj[["ll.x"]], obj[["ll.y"]], obj[["ll.x"]]+obj[["columns"]]*obj[["res.x"]], obj[["ll.y"]]+obj[["rows"]]*obj[["res.y"]]), nrow=2)
  attr(bb, "dimnames") <- list(c("x","y"), c("min","max"))
  ## tile using rows and columns:
  btiles <- makeTiles(bb, block.x, block.y, overlap.percent, limit.bbox, rows = obj[["rows"]], columns = obj[["columns"]])

  if(return.SpatialPolygons == TRUE){
    pol <- .tiles2pol(bb=bb, btiles=btiles, proj4string=sp::CRS(attr(obj, "projection")))
  } else {
    pol = btiles
  }

  message(paste("Returning a list of tiles for an object of class", class(obj), "with", signif(overlap.percent, 3), "percent overlap"))
  return(pol)

})
