\name{SpatialMemberships-class}
\docType{class}
\alias{SpatialMemberships-class}
\title{A class for membership maps derived using the \code{fkmeans} classification}
\description{A class containing a list of gridded maps and results of model fitting.}
\section{Slots}{
  \describe{
    \item{\code{predicted}:}{object of class \code{"SpatialPixelsDataFrame"}; predicted values (factor)}
    \item{\code{model}:}{object of class \code{"multinom"}; output object from the \code{nnet::multinom} method}
    \item{\code{mu}:}{object of class \code{"SpatialPixelsDataFrame"}; a list of predicted memberships}
    \item{\code{class.c}:}{object of class \code{"matrix"}; class centres}
    \item{\code{class.sd}:}{object of class \code{"matrix"}; class deviations}
    \item{\code{confusion}:}{object of class \code{"matrix"}; confusion matrix}
  }
}
\author{ Tomislav Hengl }
\seealso{ \code{\link{SpatialComponents-class}}
}
\keyword{classes}
