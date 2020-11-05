
.onLoad <- function(libname, pkgname) {
  utils::data("soil.classes", "munsell", "USDA.TT.im", "landgis.tables", package=pkgname, envir=parent.env(environment()))
}

.onAttach <- function(lib, pkg)  {
  packageStartupMessage("version: ", utils::packageDescription("landmap", field="Version"), appendLF = TRUE)
}
