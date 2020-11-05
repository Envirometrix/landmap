#' Generate Principal Components using SpatialPixelsDataFrame object
#'
#' @aliases spc
#'
#' @description Combines the \code{stats::prcomp} method and predicts a list principal components for an object of type \code{"SpatialPixelsDataFrame"}.
#'
#' @param obj SpatialPixelsDataFrame.
#' @param formulaString optional model definition.
#' @param scale. scale all numbers.
#' @param silent silent output.
#'
#' @return Object of class \code{SpatialComponents}. List of grids with generic names \code{PC1},\dots,\code{PCp}, where \code{p} is the total number of input grids.
#'
#' @note This method assumes that the input covariates are cross-correlated and hence their overlap can be reduced. The input variables are scaled by default and the missing values will be replaced with 0 values to reduce loss of data due to missing pixels.
#'
#' @author \href{https://opengeohub.org/people/tom-hengl}{Tom Hengl}
#'
#' @examples
#' library(plotKML)
#' library(sp)
#' pal = rev(rainbow(65)[1:48])
#' data(eberg_grid)
#' gridded(eberg_grid) <- ~x+y
#' proj4string(eberg_grid) <- CRS("+init=epsg:31467")
#' formulaString <- ~ PRMGEO6+DEMSRT6+TWISRT6+TIRAST6
#' eberg_spc <- spc(eberg_grid, formulaString)
#' names(eberg_spc@predicted) # 11 components on the end;
#' \donttest{
#' ## plot maps:
#' rd = range(eberg_spc@predicted@data[,1], na.rm=TRUE)
#' sq = seq(rd[1], rd[2], length.out=48)
#' spplot(eberg_spc@predicted[1:4], at=sq, col.regions=pal)
#' }
#' @export
#' @docType methods
setMethod("spc", signature(obj = "SpatialPixelsDataFrame"), function(obj, formulaString, scale. = TRUE, silent = FALSE){
  if(missing(formulaString)){
    formulaString <- stats::as.formula(paste("~", paste(names(obj), collapse="+")))
  }
  vars = all.vars(formulaString)
  if(length(vars)< 2){
    stop("At least two covarites required to run Principal Component Analysis")
  }
  obj@data <- obj@data[,vars]

  ## print warning:
  if(silent==FALSE){
  if(nrow(obj)>10e6){
    warning('Operation not recommended for large grids', immediate. = TRUE)
  }}

  ## convert every factor to indicators:
  for(j in 1:length(vars)){
    if(is.factor(obj@data[,vars[j]])){
      # remove classes without pixels:
      obj@data[,vars[j]] <- as.factor(paste(obj@data[,vars[j]]))
      ln <- levels(obj@data[,vars[j]])
      for(k in 1:length(ln)){
        vn <- paste(vars[j], k, sep="_")
        obj@data[,vn] <- ifelse(obj@data[,vars[j]]==ln[k], 1, 0)
      }
    message(paste("Converting", vars[j], "to indicators..."))
    }
  }
  varsn = names(obj)[which(!sapply(obj@data, is.factor))]
  obj@data <- obj@data[,varsn]

  ## filter the missing values:
  if(scale. == TRUE){
    x <- scale(obj@data)
    x[is.na(x)] <- 0
    x <- as.data.frame(x)
    sd.l <- lapply(x, FUN=stats::sd)
    x0 <- sd.l==0
    if(any(x0)){
      message(paste("Columns with zero variance removed:", names(x)[which(x0)]), immediate. = TRUE)
      formulaString.f = stats::as.formula(paste("~", paste(varsn[-which(x0)], collapse="+")))
      ## principal component analysis:
      pcs <- stats::prcomp(formulaString.f, x)
    } else {
      formulaString = stats::as.formula(paste("~", paste(varsn, collapse="+")))
      pcs <- stats::prcomp(formulaString, x)
    }
  } else {
    formulaString = stats::as.formula(paste("~", paste(varsn, collapse="+")))
    pcs <- stats::prcomp(formulaString, obj@data)
  }

  ## copy values:
  obj@data <- as.data.frame(pcs$x)
  sp::proj4string(obj) <- obj@proj4string
  if(silent==FALSE){
    message(paste("Converting covariates to principal components..."))
    summary(pcs)
  }

  pcs <- methods::new("SpatialComponents", predicted = obj, pca = pcs[-which(names(pcs)=="x")])
  return(pcs)

})
