\name{soil.classes}
\docType{data}
\encoding{latin1}
\alias{soil.classes}
\title{Soil classification tables}
\description{Standard soil classification tables for the United States Department of Agriculture (USDA) and IUSS / FAO World Reference Base (WRB) classification systems, including the tables used to correlate various soil classification systems.}
\usage{data(soil.classes)}
\format{
Contains a list of tables:
\describe{
  \item{\code{Canadian}}{data frame; Canadian soil classification system (Soil Classification Working Group, 1998)}
  \item{\code{FAO1990.WRB}}{data frame; FAO 1990 (FAO-Unesco Soil Classification System) system to WRB (2006)}
  \item{\code{USDA_GreatGroups}}{data frame; list of USDA Great Groups (USDA, 2010)}
  \item{\code{WRB_versions}}{data frame; correlation between various FAO/WRB versions (Krasilnikov et al. 2009)}
  \item{\code{FAO1974.WRB}}{data frame; correlation FAO 1974 system to WRB (IUSS Working Group WRB, 2006)}
  \item{\code{USDA.WRB}}{data frame; correlation USDA system to WRB system}
  \item{\code{Soils_World}}{data frame; referent soil profiles of the world (van Baren and Lof, 1987)}    
}
}
\note{Some of the original tables from the literature have been adjusted / updated by the author. Correlation between various national and international systems often leads to multiple soil classes being possible equivalents. These are separated in tables using \code{"/"} symbol e.g. \code{Dark Gray Chernozem} = \code{Boralfic Boroll / Albolls}. Some national soil classification systems contain classes which are completely unique and hence most likely can not be correlated to any class in the target system.}
\author{Tomislav Hengl}
\references{
\itemize{
\item van Baren, H. and Lof, P. (1987) Soils of the World, published by Elsevier in Association with ISRIC, FAO, and UNESCO ISBN 0444425756
\item Krasilnikov, P., Marti, J., Arnold, R., and Shoba, S., eds. (2009) A Handbook of Soil Terminology, Correlation and Classification, Earthscan LLC, pp. 448 ISBN 9781136546631
\item Soil Classification Working Group, (1998) The Canadian System of Soil Classification. 3rd Ed. Agriculture and Agri-Food Canada Publication 1646, 187 pp. ISBN 0660174049
\item USDA / Soil Survey Staff, (2010) Keys to Soil Taxonomy (Eleventh Edition) U.S. Department Of Agriculture, Natural Resources Conservation Service, 346 pp. ISBN 9781782662112
\item IUSS Working Group WRB (2006) World Reference Base for Soil Resources 2006: A Framework for International Classification, Correlation and Communication, Food and Agriculture Organization of the United Nations 128 pp.
}
}
\examples{
data(soil.classes)
soil.classes$USDA_GreatGroups[1,]
DGC <- which(soil.classes$Canadian$CSSC_Great_Groups=="Dark Gray Chernozem")
soil.classes$Canadian[DGC,]
}
\keyword{datasets}