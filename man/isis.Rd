\name{isis}
\docType{data}
\encoding{latin1}
\alias{isis}
\title{ISRIC Soil Information System}
\description{\href{https://isis.isric.org/}{ISRIC's collection of global soil monoliths} that represent the main soil reference groups of the World Reference Base for Soil Resources (WRB). Includes some 950 monoliths (785 with coordinates) from over 70 countries with detailed soil profile and environmental data.}
\usage{data(isis)}
\format{
The \code{isis} data set contains two data frames --- sites and horizons. Sites table contains the following columns:
  \describe{
  \item{\code{SOURCEID}}{factor; unique ISIS code}
  \item{\code{LONWGS84}}{numeric; longitude in decimal degrees on the WGS84 datum}
  \item{\code{LATWGS84}}{numeric; latitude in decimal degrees on the WGS84 datum}
  \item{\code{TIMESTRR}}{Date; the date on which this particular soil was described or sampled}
  \item{\code{TAXGWRB}}{factor; soil group based on the WRB classification system}
  \item{\code{TAXNUSDA}}{factor; Keys to Soil Taxonomy taxon name e.g. \code{"Natraqualf"}}
  \item{\code{BDRICM}}{numeric; depth to bedrock (R horizon) if observed}
  \item{\code{SOURCEDB}}{factor; source data base}
}
Horizons table contains the following columns:
  \describe{
  \item{\code{SOURCEID}}{factor; unique ISIS code}
  \item{\code{UHDICM}}{numeric; upper horizon depth from the surface in cm}
  \item{\code{LHDICM}}{numeric; lower horizon depth from the surface in cm}
  \item{\code{CRFVOL}}{numeric; volume percentage of coarse fragments (> 2 mm)}
  \item{\code{PHIHOX}}{numeric; pH index measured in water solution}
  \item{\code{PHIKCL}}{numeric; pH index measured in KCl solution}
  \item{\code{ORCDRC}}{numeric; soil organic carbon content in permilles}
  \item{\code{SNDPPT}}{numeric; weight percentage of the sand particles (0.05--2 mm)}
  \item{\code{SLTPPT}}{numeric; weight percentage of the silt particles (0.0002--0.05 mm)}
  \item{\code{CLYPPT}}{numeric; weight percentage of the clay particles (<0.0002 mm)}
  \item{\code{CEC}}{numeric; Cation Exchange Capacity in cmol+/kg}
  \item{\code{BLD}}{bulk density in tonnes per cubic-meter}
}
}
\author{\href{https://isis.isric.org/}{ISRIC --- World Soil Information}}
\examples{
library(rgdal)
library(sp)

data(isis)
sites <- isis$sites
coordinates(sites) <- ~ LONWGS84 + LATWGS84
proj4string(sites) <- "+proj=longlat +datum=WGS84"
\dontrun{
## obtain country borders:
library(maps)
country.m = map('world', plot=FALSE, fill=TRUE)
IDs <- sapply(strsplit(country.m$names, ":"), function(x) x[1])
require(maptools)
country <- as(map2SpatialPolygons(country.m, IDs=IDs), "SpatialLines")
proj4string(country) = "+proj=longlat +datum=WGS84"
## overlay and plot points and maps:
plot(country, col="darkgrey")
points(sites, pch=21, bg="red", cex=.6, col="black")
}
}
\keyword{datasets}
