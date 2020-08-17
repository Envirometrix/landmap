\name{SpatialComponents-class}
\docType{class}
\alias{SpatialComponents-class}
\title{A class for gridded components derived using the \code{spc} method}
\description{A class containing a list of gridded components and results of principal component analysis.}
\section{Slots}{
  \describe{
    \item{\code{predicted}:}{object of class \code{"SpatialPixelsDataFrame"}; predicted values for components}
    \item{\code{pca}:}{object of class \code{"list"}; output objects from the \code{stats::prcomp} process --- contains objects: 'stdev', 'rotation', 'center' and 'scale'}
  }
}
\author{ Tomislav Hengl }
\keyword{classes}
