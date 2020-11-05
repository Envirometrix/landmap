\name{edgeroi}
\docType{data}
\alias{edgeroi}
\alias{edgeroi.grids}
\alias{edgeroi.grids100}
\title{The Edgeroi Data Set}
\description{Soil samples and covariate layers for the Edgeroi area in NSW, Australia (ca 1500 square-km).}
\usage{data(edgeroi)}
\format{
The \code{edgeroi} data set contains two data frames --- sites and horizons. Sites table contains the following columns:
  \describe{
  \item{\code{SOURCEID}}{factor; unique label to help a user identify a particular site (\code{ID} in the \href{http://www.asris.csiro.au/mapping/hyperdocs/NatSoil/}{NatSoil})}
  \item{\code{LONGDA94}}{numeric; longitude in decimal degrees on the GDA94 datum}
  \item{\code{LATGDA94}}{numeric; latitude in decimal degrees on the GDA94 datum}
  \item{\code{TAXGAUC}}{factor; Australian Great Soil Groups (GSG; see details)}
  \item{\code{NOTEOBS}}{character; free-form observation notes}
}
Horizons table contains the following columns:
  \describe{
  \item{\code{SOURCEID}}{factor; unique identifier used in the NatSoil DB}
  \item{\code{LSQINT}}{integer; a layer sequence number 1 to N}
  \item{\code{HZDUSD}}{factor; horizon designation (primary letter)}
  \item{\code{UHDICM}}{numeric; upper horizon depth from the surface in cm}
  \item{\code{LHDICM}}{numeric; lower horizon depth from the surface in cm}
  \item{\code{CLYPPT}}{numeric; weight percentage of the clay particles (<0.0002 mm)}
  \item{\code{SNDPPT}}{numeric; weight percentage of the silt particles (0.0002--0.05 mm)}
  \item{\code{SLTPPT}}{numeric; weight percentage of the sand particles (0.05--2 mm)}
  \item{\code{PHIHO5}}{numeric; pH index measured in water solution(\code{ph_h2o} in the NSCD)}
  \item{\code{ORCDRC}}{numeric; soil organic carbon content in permille}
}
The \code{edgeroi.grids} data frame contains a list of covariates at 250 m resolution:
  \describe{
  \item{\code{DEMSRT5}}{numeric; SRTM DEM}
  \item{\code{TWISRT5}}{numeric; SAGA Topographic Wetness Index based on the SRTM DEM}
  \item{\code{PMTGEO5}}{factor; parent material class based on the National Geological map at scale 1:250,000 --- sand with minor silty sand (\code{"Qd"}), alluvium gravel, sand, silt, clay (\code{"Qrs"}), quartz sandstone obscured by quartenary sands (\code{"Qrt/Jp"}), quartz sandstone obscured by talus material (\code{"Qrt/Rn"}), basalt obscured by talus material (\code{"Qrt/Tv"}), mottled clay, silt, sandstone and gravel (\code{"Ts"}), and basalt, dolerite, trachyte, techenite (\code{"Tv"})}
  \item{\code{EV1MOD5}}{numeric; first principal component of the MODIS EVI (MOD13Q1) time series data (year 2011)}
  \item{\code{EV2MOD5}}{numeric; second principal component of the MODIS EVI (MOD13Q1) time series data (year 2011)}
  \item{\code{EV3MOD5}}{numeric; third principal component of the MODIS EVI (MOD13Q1) time series data (year 2011)}
  \item{\code{x}}{numeric; x-coordinate in the GDA94 / MGA zone 55}
  \item{\code{y}}{numeric; y-coordinate in the GDA94 / MGA zone 55}
}
The \code{edgeroi.grids100} data frame contains a list of covariates at 100 m resolution prepared for the study area:
  \describe{
  \item{\code{LNUABS6}}{factor; Australian National scale land use data}
  \item{\code{MVBSRT6}}{numeric; SAGA GIS Multi-resolution Index of Valley Bottom Flatness based on the SRTM DEM}
  \item{\code{TI1LAN6}}{numeric; principal component 1 for the Landsat band 7 (thermal) based on three periods of the Global Land Survey Landsat images (GLS1990, GLS2000, GLS2005)}
  \item{\code{TI2LAN6}}{numeric; principal component 2 for the Landsat band 7 (thermal) based on three periods of the Global Land Survey Landsat images (GLS1990, GLS2000, GLS2005)}
  \item{\code{PCKGAD6}}{numeric; percentage of Potassium estimated based on the gamma radiometrics radmap09 (GADDS)}
  \item{\code{RUTGAD6}}{numeric; ratio Uranium over Thorium estimated based on the gamma radiometrics radmap09 (GADDS)}
  \item{\code{PCTGAD6}}{numeric; parts per million of Thorium estimated based on the gamma radiometrics radmap09 (GADDS)}
  \item{\code{x}}{numeric; x-coordinate in the GDA94 / MGA zone 55}
  \item{\code{y}}{numeric; y-coordinate in the GDA94 / MGA zone 55}
}
}
\details{The Edgeroi is one of the standard soil data sets used to test soil mapping methods in Australia. Out of 359 profiles, 210 sites were sampled on a systematic, equilateral triangular grid with a spacing of 2.8 km between sites, the other sites are distributed more irregularly or on transects. The data set is described in detail in \href{http://dx.doi.org/10.1016/j.geoderma.2009.10.007}{Malone et al. (2010)} and McGarry et al. (1989). The \code{edgeroi} contains only a subset of the original \href{http://www.asris.csiro.au/mapping/hyperdocs/NatSoil/}{NatSoil} records. Observed soil classes for \code{TAXGAUC} are (alphabetically): Alluvial soil (\code{"A"}), Brown clay (\code{"BC"}), Black earth (\code{"BE"}), Earthy sand (\code{"ES"}), Grey clay (\code{"GC"}), Grey earth (\code{"GE"}), No suitable group (\code{"NSG"}), Prairie soil (\code{"PS"}), Rendzina (\code{"R"}), Red-brown earth (\code{"RBE"}), Red clay (\code{"RC"}), Red earth (\code{"RE"}), Red podzolic soil (\code{"RP"}), Solodic soil (\code{"SC"}), Soloth (\code{"SH"}), Solonchak (\code{"SK"}), Siliceous sand (\code{"SS"}), and Solonetz (\code{"SZ"}).}
\note{The Landsat images and SRTM DEM have been obtained from the \href{https://geog.umd.edu/}{Global Land Cover Facility}. Scanned geology map (paper sheets) has been obtained from the \href{http://www.geoscience.gov.au}{Geoscience Australia}, then georeferenced and rasterized to 250 m resolution. The land use map has been obtained from the Australian Collaborative Land Use and Management program. The Radiometric Map of Australia grids has been downloaded using the Geophysical Archive Data Delivery System (GADDS) on the Australian Government's Geoscience Portal (\href{http://dx.doi.org/10.1071/EG09025}{Mitny et al, 2009}).\cr
Listed gridded layers follow a standard naming convention used by WorlGrids.org (the standard 8.3 filename convention with at most eight characters): first three letter are used for the variable type e.g. \code{DEM} (digital elevation model); the next three letters represent the data source or collection method e.g. \code{SRT} (SRTM mission); the 6th character is the effective scale e.g. \code{5} indicates the 5th standard scale i.e. 1/600 decimal degrees (in this case 250 m).}
\author{The \href{http://www.asris.csiro.au/mapping/hyperdocs/NatSoil/}{original detailed profile description and laboratory analysis} was funded by a Cotton Research and Development Corporation project in the mid-late 1980's by the CSIRO Division of Soils and available via the \href{http://www.asris.csiro.au/mapping/hyperdocs/NatSoil/}{NatSoil} DB. The gamma radiometrics images are property of the NSW Department of Primary Industries --- Mineral Resources.}
\references{
\itemize{
\item Malone, B.P., McBratney, A.B., Minasny, B. (2010) \href{http://dx.doi.org/10.1016/j.geoderma.2009.10.007}{Mapping continuous depth functions of soil carbon storage and available water capacity}. Geoderma 154, 138-152.
\item McGarry, D., Ward, W.T., McBratney, A.B. (1989) Soil Studies in the Lower Namoi Valley: Methods and Data. The Edgeroi Data Set. (2 vols) (CSIRO Division of Soils: Adelaide).
\item Minty, B., Franklin, R., Milligan, P., Richardson, L.M., and Wilford, J., (2009) \href{http://dx.doi.org/10.1071/EG09025}{The Radiometric Map of Australia}. Exploration Geophysics, 40(4), 325-333.
}
}
\examples{
library(rgdal)
library(aqp)
library(sp)

data(edgeroi)
edgeroi$sites[edgeroi$sites$SOURCEID=="399_EDGEROI_ed095_1",]
edgeroi$horizons[edgeroi$horizons$SOURCEID=="399_EDGEROI_ed095_1",]
## spPoints:
sites <- edgeroi$sites
coordinates(sites) <- ~ LONGDA94 + LATGDA94
proj4string(sites) <- CRS("+proj=longlat +ellps=GRS80 +towgs84=0,0,0,0,0,0,0 +no_defs")
sites <- spTransform(sites, CRS("+init=epsg:28355"))

\donttest{
## plot points and grids:
pnts <- list("sp.points", sites, pch="+", col="black")
## load the 250 m grids:
library(RCurl)
rep = "https://raw.github.com/Envirometrix/PredictiveSoilMapping/master/extdata/"
x = tempfile(fileext = ".rda")
con <- download.file(paste0(rep, "edgeroi.grids.rda"), x, method="wget")
load(x)
str(edgeroi.grids)
gridded(edgeroi.grids) <- ~x+y
proj4string(edgeroi.grids) <- CRS("+init=epsg:28355")
spplot(edgeroi.grids[1], sp.layout=pnts)
## load the 100 m grids:
x2 = tempfile(fileext = ".rda")
con2 <- download.file(paste0(rep, "edgeroi.grids100.rda"), x2, method="wget")
load(x2)
str(edgeroi.grids100)
gridded(edgeroi.grids100) <- ~x+y
proj4string(edgeroi.grids100) <- CRS("+init=epsg:28355")
spplot(edgeroi.grids100["TI1LAN6"], sp.layout=pnts)
}
}
\keyword{datasets}
