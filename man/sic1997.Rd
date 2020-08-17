\name{sic1997}
\docType{data}
\alias{sic1997}
\title{The SIC 1997 Data Set}
\description{Daily rainfall dataset from Switzerland used in the Spatial Interpolation Comparison 1997.}
\usage{data(sic1997)}
\format{
The \code{sic1997} data set contains two objects --- stations (points) and covariates (1 km grids). The \code{daily.rainfall} contains the following columns:
  \describe{
  \item{\code{rainfall}}{numeric; daily rainfall measurements at 467 meteo stations on the 8th of May 1986}
  \item{\code{X}}{numeric; easting in the local coordinate system}
  \item{\code{Y}}{numeric; northing in the local coordinate system}
}
Object \code{swiss1km} contains the following columns:
  \describe{
  \item{\code{CHELSA_rainfall}}{numeric; monthly rainfall based on the CHELSA climate grids}
  \item{\code{DEM}}{numeric; digital elevation based on the SRTM DEM}
  \item{\code{border}}{character; country border mask}
}
}
\details{The gathering of the rainfall data, provided by Giovanni Graziani from the Environment Institute of the Joint Research Centre (Ispra, Italy), has been undertaken, under JRC-Ispra funding, by the Air pollution
Group at Imperial College, London. The Digital Elevation Model has been provided by EROS Data Centre from the U.S. Geological Survey (USGS). http://edcwww.cr.usgs.gov/.}
\references{
\itemize{
\item Dubois, G. (1998). Spatial interpolation comparison 97: foreword and introduction. Journal of Geographic Information and Decision Analysis, 2(2), 1-10.
\item https://wiki.52north.org/AI_GEOSTATS/WebHome
}
}
\examples{
data(sic1997)
}
\keyword{datasets}
