\name{landgis.tables}
\docType{data}
\alias{landgis.tables}
\title{OpenLandMap list of GeoTIFFs global mosaics}
\description{Computer names of the layers in the www.OpenLandMap.org.}
\usage{data(landgis.tables)}
\format{
The \code{landgis.tables} data set contains four objects:
  \describe{
  \item{\code{tables}}{data.frame; layer metadata}
  \item{\code{layers}}{list; list of layers at different resolutions}
  \item{\code{classes}}{list; list of classes}
  \item{\code{zenodo.files}}{list; list of URI addressed to download all OpenLandMap.org layers}
}
}
\references{
\itemize{
\item https://openlandmap.org
}
}
\examples{
data(landgis.tables)
str(landgis.tables$tables[1,])
}
\keyword{datasets}
