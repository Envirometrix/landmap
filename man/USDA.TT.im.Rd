\name{USDA.TT.im}
\docType{data}
\alias{USDA.TT.im}
\title{Probability density for texture triangle}
\description{Probability density for texture triangle (USDA system) based on global soil profile data (see ISRIC WoSIS).}
\usage{data(USDA.TT.im)}
\format{
The \code{USDA.TT.im} data frame contains the following columns:
  \describe{
  \item{\code{v}}{numeric; probability density derived using the \code{soiltexture::TT.kde2d} function and global soil profile data}
  \item{\code{TEXMHT}}{factor; USDA soil texture class estimated by hand (one of the following: \code{"C"}, \code{"SiC"}, \code{"SC"}, \code{"CL"}, \code{"SiCL"}, \code{"SCL"}, \code{"L"}, \code{"SiL"}, \code{"SL"}, \code{"Si"}, \code{"LS"}, \code{"S")}}
  \item{\code{s1}}{numeric; horizontal coordinate (sand content 0--1) in the texture triangle system}
  \item{\code{s2}}{numeric; vertical coordinate (0--0.85) in the texture triangle system}
}
}
\note{Texture by hand class can be converted to sand, silt, clay content fractions by using the \code{TT2tri} function. This function uses the \code{v} column in the \code{USDA.TT.im} (i.e. prior probability densities) to adjust for texture fraction combinations that are more probable.}
\author{ Tomislav Hengl }
\references{
\itemize{
\item Skaggs, T. H., Arya, L. M., Shouse, P. J., Mohanty, B. P., (2001) Estimating Particle-Size Distribution from Limited Soil Texture Data. Soil Science Society of America Journal 65 (4): 1038-1044.
}
}
\examples{
## plot prior probabilities:
library(sp)
data(USDA.TT.im)
gridded(USDA.TT.im) <- ~s1+s2
spplot(USDA.TT.im["v"])
}
\keyword{datasets}
