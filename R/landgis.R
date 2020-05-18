#' Access and download layers from OpenLandMap.org (LandGIS data service)
#'
#' @param coverageId Coverage ID.
#' @param filename Download filename.
#' @param scalefactor Scale factor for WCS request.
#' @param subset Subset string for WCS request.
#' @param service URL of the WCS service.
#' @param silent Silent output.
#'
#' @return Locally downloaded GeoTIFF.
#' @export
#'
#' @author \href{https://opengeohub.org/people/tom-hengl}{Tom Hengl}
#'
#' @examples
#' library(plotKML)
#' library(rgdal)
#' data("sic1997")
#' search.landgis(pattern=c("clay", "10..10cm"))
#' search.landgis(pattern=c("clay", "10..10cm"))
#' \dontrun{
#' coverageId = "predicted250m:sol_clay.wfraction_usda.3a1a1a_m_250m_b10..10cm_1950..2017_v0.2"
#' swiss1km.ll <- raster::projectRaster(raster(sic1997$swiss1km[2]), crs = "+init=epsg:4326", res=c(1/120, 1/120))
#' s1 = paste0("Lat(", swiss1km.ll@extent@ymin, ",", swiss1km.ll@extent@ymax,")")
#' s2 = paste0("Long(", swiss1km.ll@extent@xmin, ",", swiss1km.ll@extent@xmax,")")
#' sf = (1/480)/(1/120)
#' download.landgis(coverageId, filename = "clay_ch1km.tif", subset = c(s1,s2), scalefactor = sf)
#' swiss1km.ll1km <- as(swiss1km.ll, "SpatialGridDataFrame")
#' swiss1km.ll1km$clay_10..10cm <- readGDAL("clay_ch1km.tif")$band1
#' swiss1km.ll1km$clay_10..10cm <- ifelse(is.na(swiss1km.ll1km$DEM), NA, swiss1km.ll1km$clay_10..10cm)
#' plot(swiss1km.ll1km["clay_10..10cm"])
#' ## mapview(swiss1km.ll1km["clay_10..10cm"])
#' }
download.landgis <- function(coverageId, filename, scalefactor=NULL, subset=NULL, service="https://geoserver.opengeohub.org/landgisgeoserver/ows?service=WCS&version=2.0.1", silent=TRUE, ...){
  if(!is.null(scalefactor)){ scalefactor <- paste0('&scalefactor=', scalefactor) }
  if(!is.null(subset)){ subset <- paste0('&subset=', subset[1], '&subset=', subset[2]) }
  wcs.url <- paste0(service, '&request=GetCoverage&coverageId=', coverageId, subset, scalefactor)
  x <- download.file(wcs.url, filename, quiet = silent, ...)
  try( obj <- rgdal::GDALinfo(filename, silent = silent) )
  if(!class(obj)=="GDALobj"){
    ## download from zenodo?
    x <- search.landgis(strsplit(coverageId, ":")[[1]][2])
    warning(paste("WCS request exceed memory limitations. Try downloading from zenodo:\n", x[[2]]))
  } else {
    if(silent==FALSE){
      message("Downloaded layer successfully.")
    }
  }
}

#' Search for available landgis layers
#'
#' @param pattern String pattern
#' @param layersURL Default URL with the list of layers
#' @param update Logical specify to update the layer list
#'
#' @return
#' @export
search.landgis <- function(pattern, layersURL="https://landgisapi.opengeohub.org/query/layers", update=FALSE){
  #pattern=c("clay", "10..10cm")
  data("landgis.tables")
  if(update==TRUE){
    x1 <- tempfile(fileext = ".csv")
    t1 <- download.file("https://raw.githubusercontent.com/Envirometrix/LandGISmaps/master/tables/LandGIS_tables_landgis_layers.csv", x1)
    t1 <- read.csv(x1)
    t0.n <- c("dtm_landform_usgs.ecotapestry_c_250m_s0..0cm_2014_v1.0.tif.csv", "dtm_lithology_usgs.ecotapestry_c_250m_s0..0cm_2014_v1.0.tif.csv", "lcv_land.cover_esacci.lc.l4_c.csv", "pnv_biome.type_biome00k_c_1km_s0..0cm_2000..2017_v0.1.tif.csv", "sol_grtgroup_usda.soiltax_c_250m_s0..0cm_1950..2017_v0.1.tif.csv", "sol_texture.class_usda.tt_m_250m_b_1950..2017_v0.1.tif.csv")
    t0 <- lapply(t0.n, function(i){ x <- tempfile(fileext = ".csv"); download.file(paste0("https://raw.githubusercontent.com/Envirometrix/LandGISmaps/master/tables/", i), x); read.csv(x) })
    names(t0) = t0.n
    layers <- rjson::fromJSON(RCurl::getURL(layersURL))
    dep.id = unique(t1$layer_zenodo_deposit)
    dep.id = dep.id[!is.na(dep.id)]
    TOKEN = scan("~/TOKEN_ACCESS", what="character")
    z0 <- list(NULL)
    for(i in 1:length(dep.id)){
      z0[[i]] <- rjson::fromJSON(system(paste0('curl -H \"Accept: application/json\" -H \"Authorization: Bearer ', TOKEN, '\" \"https://www.zenodo.org/api/deposit/depositions/', dep.id[i], '\"'), intern=TRUE))
    }
    landgis.tables <- list(tables=t1, layers=layers, classes=t0, zenodo.files=z0)
    #save(landgis.tables, file = "data/landgis.tables.rda", compress = "xz")
  }
  layers.files <- unlist(landgis.tables$layers)
  layers.files <- layers.files[grep(pattern=glob2rx("*.tif$"), layers.files)]
  out1 <- layers.files[grepl(pattern[1], paste(layers.files))]
  if(length(pattern)==2){ out1 <- out1[grepl(pattern[2], paste(out1))]  }
  zenodo.files <- unlist(lapply(landgis.tables$zenodo.files, function(x){ sapply(x$files, function(i){ i$links$download }) }))
  out2 <- zenodo.files[grepl(pattern[1], paste(zenodo.files))]
  if(length(pattern)==2){ out2 <- out2[grepl(pattern[2], paste(out2))]  }
  list(out1, out2)
}

#over.landgis <- function(coverageId, points){
#}
