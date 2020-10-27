#' Tile spatial layers
#'
#' @aliases tile tile,SpatialLinesDataFrame-method tile,SpatialPixelsDataFrame-method tile,SpatialPointsDataFrame-method tile,SpatialPolygonsDataFrame-method
#' @rdname tile-methods
#'
#' @param x RasterLayer.
#' @param y either points, pixels, polygons or lines.
#' @param block.x size of the block in x direction.
#' @param tmp.file temporary file name.
#' @param show.output.on.console shows progress.
#' @param program optional location of the gdalwarp.
#' @param ... optional argument.
#'
#' @usage
#' \S4method{tile}{SpatialPointsDataFrame}(x, y, block.x, \dots)
#' \S4method{tile}{SpatialPixelsDataFrame}(x, y, block.x, \dots)
#' \S4method{tile}{SpatialPolygonsDataFrame}(x, y, block.x, tmp.file = TRUE,
#'                                           program, show.output.on.console = FALSE, \dots)
#' \S4method{tile}{SpatialLinesDataFrame}(x, y, block.x, tmp.file = TRUE,
#'                                        program, show.output.on.console = FALSE, \dots)
#' \S4method{tile}{RasterLayer}(x, y, block.x, tmp.file = TRUE,
#'                              program, show.output.on.console = FALSE, \dots)
#'
#' @importClassesFrom raster RasterLayer
#'
#' @return Regular tiling system
#'
#' @author \href{https://opengeohub.org/people/tom-hengl}{Tom Hengl}
#' @export
#' @docType methods
setMethod("tile", signature(x = "RasterLayer"), function(x, y, block.x, tmp.file = TRUE, program, show.output.on.console = FALSE, ...){

  if(raster::filename(x)==""){
    stop("Function applicable only to 'RasterLayer' objects linking to external raster files")
  }

  if(missing(program)){
    program <- .programPath(utility="gdalwarp")
  }

  if(missing(y)){
    b <- sp::bbox(x)
    pol <- sp::SpatialPolygons(list(sp::Polygons(list(sp::Polygon(matrix(c(b[1,1], b[1,1], b[1,2], b[1,2], b[1,1], b[2,1], b[2,2], b[2,2], b[2,1], b[2,1]), ncol=2))), ID="1")), proj4string=sp::CRS(sp::proj4string(x)))
    y <- getSpatialTiles(pol, block.x = block.x, return.SpatialPolygons = FALSE, ...)
  }
  ## gdalwarp by tiles:
  x.lst <- list()
  message("Clipling raster object using 'gdalwarp'...")
  for(j in 1:nrow(y)){
    if(tmp.file==TRUE){
      outname <- tempfile()
    } else {
      outname <- paste(plotKML::normalizeFilename(deparse(substitute(x, env = parent.frame()))), j, sep="_")
    }
    try(system(paste(program, utils::shortPathName(normalizePath(raster::filename(x))), RSAGA::set.file.extension(outname, ".tif"), '-te',  y[j,1], y[j,2], y[j,3], y[j,4]), show.output.on.console = show.output.on.console))
    try(x.lst[[j]] <- raster::raster(RSAGA::set.file.extension(outname, ".tif")))
  }
  return(x.lst)

})

## tile points:
.subsetTiles <- function(x, y, block.x, ...){
  if(missing(y)){
    y <- getSpatialTiles(x, block.x = block.x, ...)
  }
  ## subset by tiles:
  ov <- sp::over(x, y)
  t.lst <- sapply(methods::slot(y, "polygons"), methods::slot, "ID")
  bbox.lst <- lapply(methods::slot(y, "polygons"), sp::bbox)
  message(paste('Subseting object of class \"', class(x), '\"...', sep=""))
  x.lst <- list()
  for(i in 1:length(y)){
    sel <- ov == t.lst[i]
    if(length(sel)>0&!all(is.na(sel))){
      if(!all(sel==FALSE)){
        x.lst[[i]] <- subset(x, subset=sel)
        x.lst[[i]]@bbox <- bbox.lst[[i]]
      }
    }
  }
  return(x.lst)
}

.tiles2pol <- function(bb, btiles, proj4string){
  ## get coordinates for each tile:
  coords.lst <- lapply(as.list(as.data.frame(t(as.matrix(btiles)))), function(x){matrix(c(x[1], x[1], x[3], x[3], x[1], x[2], x[4], x[4], x[2], x[2]), ncol=2, dimnames=list(paste("p", 1:5, sep=""), attr(bb, "dimnames")[[1]]))})

  ## create an object of class "SpatialPolygons"
  srl = lapply(coords.lst, sp::Polygon)
  Srl = list()
  for(i in 1:length(srl)){ Srl[[i]] <- sp::Polygons(list(srl[[i]]), ID=row.names(btiles)[i]) }
  pol = sp::SpatialPolygons(Srl, proj4string=proj4string)

  return(pol)
}

## tile using the OGR2OGR function:
.clipTiles <- function(x, y, block.x, tmp.file = TRUE, program, show.output.on.console = FALSE, ...){

  if(missing(program)){
    program <- .programPath(utility="ogr2ogr")
  }

  if(missing(y)){
    y <- getSpatialTiles(x, block.x = block.x, return.SpatialPolygons = FALSE, ...)
  }

  ## write to shape file:
  if(tmp.file==TRUE){
    tf <- tempfile()
  } else {
    tf <- plotKML::normalizeFilename(deparse(substitute(x, env = parent.frame())))
  }
  suppressMessages( rgdal::writeOGR(x, RSAGA::set.file.extension(tf, ".shp"), layer=".", driver="ESRI Shapefile", overwrite_layer=TRUE) )

  ## clip by tiles:
  x.lst <- list()
  message("Clipling lines using 'ogr2ogr'...")
  for(j in 1:nrow(y)){
    if(tmp.file==TRUE){
      outname <- tempfile()
    } else {
      outname <- paste(plotKML::normalizeFilename(deparse(substitute(x, env = parent.frame()))), j, sep="_")
    }
    layername <- basename(sub("[.][^.]*$", "", outname, perl=TRUE))

    if(class(x)=="SpatialPolygonsDataFrame"){
      try(system(paste(program, '-where \"OGR_GEOMETRY=\'Polygon\'\" -f \"ESRI Shapefile\"', RSAGA::set.file.extension(outname, ".shp"), RSAGA::set.file.extension(tf, ".shp"), '-clipsrc',  y[j,1], y[j,2], y[j,3], y[j,4], '-skipfailures'), show.output.on.console = show.output.on.console))
      try(x.lst[[j]] <- rgdal::readOGR(normalizePath(RSAGA::set.file.extension(outname, ".shp")), layername, verbose = FALSE))
    }
    if(class(x)=="SpatialLinesDataFrame"){
      try(system(paste(program, '-where \"OGR_GEOMETRY=\'Linestring\'\" -f \"ESRI Shapefile\"', RSAGA::set.file.extension(outname, ".shp"), RSAGA::set.file.extension(tf, ".shp"), '-clipsrc',  y[j,1], y[j,2], y[j,3], y[j,4], '-skipfailures'), show.output.on.console = show.output.on.console))
      try(x.lst[[j]] <- rgdal::readOGR(normalizePath(RSAGA::set.file.extension(outname, ".shp")), layername, verbose = FALSE))
    }
  }
  return(x.lst)
}

.programPath <- function(path, utility){
  if(missing(path)){
    if(!file.exists("C:/PROGRA~1/GDAL/")&.Platform$OS.type == "windows"){
      if(requireNamespace("gdalUtils", quietly = TRUE)){
        path <- getOption("gdalUtils_gdalPath")[[1]]$path
        if(is.null(path)){
          ## force gdal installation:
          gdalUtils::gdal_setInstallation()
          message("Forcing installation of GDAL utilities... this might take time.")
          path <- getOption("gdalUtils_gdalPath")[[1]]$path
        }
      }
    }
    if(file.exists(paste0("C:/PROGRA~1/GDAL/", utility, ".exe"))&.Platform$OS.type == "windows"){
      program = shQuote(utils::shortPathName(normalizePath(file.path("C:/PROGRA~1/GDAL/", paste0(utility, ".exe")))))
    }
  }

  if(.Platform$OS.type == "windows") {
    program = shQuote(utils::shortPathName(normalizePath(file.path(path, paste(utility, ".exe", sep="")))))
  } else {
    program = utility
  }
  return(program)
}

setMethod("tile", signature(x = "SpatialPointsDataFrame"), .subsetTiles)
setMethod("tile", signature(x = "SpatialPixelsDataFrame"), .subsetTiles)
setMethod("tile", signature(x = "SpatialPolygonsDataFrame"), .clipTiles)
setMethod("tile", signature(x = "SpatialLinesDataFrame"), .clipTiles)
