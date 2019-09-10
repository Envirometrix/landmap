#' Estimate occurrence probabilities of a sampling plan (points)
#'
#' @description Estimates occurrence probabilities as an average between the kernel density estimation (spreading of points in geographical space) and MaxLike analysis (spreading of points in feature space). The output \code{'iprob'} indicates whether the sampling plan has systematically missed some important locations / features, and can be used as an input for geostatistical modelling (e.g. as weights for regression modeling).
#'
#' @param observations SpatialPoints.
#' @param covariates SpatialPixelsDataFrame.
#'
#' @return Returns a list of objects where \code{'iprob'} (\code{"SpatialPixelsDataFrame"}) is the map showing the estimated occurrence probabilities.
#' @export
#'
#' @note Occurrence probabilities for geographical space are derived using kernel density estimator. The sampling intensities are converted to probabilities by deviding the sampling intensity by the maximum sampling intensity for the study area (Baddeley, 2008). The occurrence probabilities for feature space are determined using MaxLike algorithm (Royle et al., 2012). The lower the average occurrence probability for the whole study area, the lower the representation efficiency of a sampling plan. \cr MaxLike function might fail to produce predictions (e.g. if not at least one continuous covariate is provided and if the \code{optim} function is not able to find the global optima) in which case an error message is generated. Running Principal Component analysis i.e. standardizing the covariates prior to running \code{spsample.prob} is, thus, highly recommended.\cr This function can be time consuming for large grids.
#'
#' @author \href{https://opengeohub.org/people/tom-hengl}{Tom Hengl}
#'
#' @examples
#' library(plotKML)
#' library(maxlike)
#' library(spatstat)
#' library(maptools)
#'
#' data(eberg)
#' data(eberg_grid)
#' ## existing sampling plan:
#' sel <- runif(nrow(eberg)) < .2
#' eberg.xy <- eberg[sel,c("X","Y")]
#' coordinates(eberg.xy) <- ~X+Y
#' proj4string(eberg.xy) <- CRS("+init=epsg:31467")
#' ## covariates:
#' gridded(eberg_grid) <- ~x+y
#' proj4string(eberg_grid) <- CRS("+init=epsg:31467")
#' ## convert to continuous independent covariates:
#' formulaString <- ~ PRMGEO6+DEMSRT6+TWISRT6+TIRAST6
#' eberg_spc <- spc(eberg_grid, formulaString)
#'
#' ## derive occurrence probability:
#' covs <- eberg_spc@predicted[1:8]
#' iprob <- spsample.prob(eberg.xy, covs)
#' ## Note: obvious omission areas:
#' hist(iprob[[1]]@data[,1])
#' ## compare with random sampling:
#' rnd <- spsample(eberg_grid, type="random",
#'                 n=length(iprob[["observations"]]))
#' iprob2 <- spsample.prob(rnd, covs)
#' ## compare the two:
#' par(mfrow=c(1,2))
#' plot(raster(iprob[[1]]), zlim=c(0,1), col=SAGA_pal[[1]])
#' points(iprob[["observations"]])
#' plot(raster(iprob2[[1]]), zlim=c(0,1), col=SAGA_pal[[1]])
#' points(iprob2[["observations"]])
#'
#' ## fit a weighted lm:
#' eberg.xy <- eberg[sel,c("SNDMHT_A","X","Y")]
#' coordinates(eberg.xy) <- ~X+Y
#' proj4string(eberg.xy) <- CRS("+init=epsg:31467")
#' eberg.xy$iprob <- over(eberg.xy, iprob[[1]])$iprob
#' eberg.xy@data <- cbind(eberg.xy@data, over(eberg.xy, covs))
#' fs <- as.formula(paste("SNDMHT_A ~ ",
#'                        paste(names(covs), collapse="+")))
#' ## the lower the occurrence probability, the higher the weight:
#' w <- 1/eberg.xy$iprob
#' m <- lm(fs, eberg.xy, weights=w)
#' summary(m)
#' ## compare to standard lm:
#' m0 <- lm(fs, eberg.xy)
#' summary(m)$adj.r.squared
#' summary(m0)$adj.r.squared
setMethod("spsample.prob", signature(observations = "SpatialPoints", covariates = "SpatialPixelsDataFrame"), function(observations, covariates, quant.nndist=.95, n.sigma, ...){

  ## mask out missing combinations:
  covariates <- covariates[stats::complete.cases(covariates@data),]
  ov <- sp::over(observations, covariates)
  observations <- observations[stats::complete.cases(ov),]

  if(requireNamespace("spatstat", quietly = TRUE)&requireNamespace("maxlike", quietly = TRUE)){
    mg_owin <- spatstat::as.owin(data.frame(x = covariates@coords[,1], y = covariates@coords[,2], window = TRUE))
    suppressWarnings( locs.ppp <- spatstat::ppp(x=coordinates(observations)[,1], y=coordinates(observations)[,2], window=mg_owin) )
    dist.locs <- spatstat::nndist(locs.ppp)
    ## inlcusion probabilities geographical space:
    if(missing(n.sigma)){
      n.sigma <- quantile(dist.locs, quant.nndist)
    }
    if(n.sigma < 0.5*sqrt(length(covariates)*covariates@grid@cellsize[1]*covariates@grid@cellsize[2]/length(observations))){
        warning(paste0("'Sigma' set at ", signif(n.sigma, 3), ". Consider increasing the value."))
    }
    message(paste("Deriving kernel density map using sigma", signif(n.sigma, 3), "..."))
    dmap <- maptools::as.SpatialGridDataFrame.im(density(locs.ppp, sigm=n.sigma, ...))
    ## Pixel values are estimated intensity values, expressed in 'points per unit area' (hence multiply by area).
    dmap.max <- max(dmap@data[,1], na.rm=TRUE)
    dmap@data[,1] <- signif(dmap@data[,1]/dmap.max, 3)

    ## occurrence probabilities in feature space:
    message("Deriving inclusion probabilities using MaxLike analysis...")
    fm <- as.formula(paste("~", paste(names(covariates), collapse="+")))
    ml <- maxlike::maxlike(formula=fm, rasters=stack(covariates), points=observations@coords, method="BFGS", savedata=TRUE)
    ## bug in "maxlike" (https://github.com/rbchan/maxlike/issues/1); need to replace this 'by hand':
    ml$call$formula <- fm
    ## TH: this operation can be time consuming and is not recommended for large grids
    ml.p <- predict(ml)
    ml.p <- as(ml.p, "SpatialGridDataFrame")
    ml.max <- max(ml.p@data[,1], na.rm=TRUE)
    ml.p@data[,1] <- signif(ml.p@data[,1]/ml.max, 3)
    ## sum two occurrence probabilities (masks for the two maps need to be exactly the same):
    covariates$iprob <- signif((ml.p@data[covariates@grid.index,1] + dmap@data[covariates@grid.index,1])/2, 3)

    out <- list(prob=covariates["iprob"], observations=as(observations, "SpatialPoints"), density=dmap, maxlike=ml.p, maxlikeFit=ml[-which(names(ml)=="rasters")])
    return(out)
  }  else {
    stop("Missing packages 'maxlike' and/or 'spatstat'")
  }
})
