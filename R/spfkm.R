#' Fit a supervised fuzzy kmeans model and predict memberships
#'
#' @aliases spfkm
#' @rdname spfkm
#'
#' @description Runs supervised fuzzy \emph{k}-means (\href{http://dx.doi.org/10.1080/13658810310001620924}{Hengl et al., 2004}) using a list of covariates layers provided as \code{"SpatialPixelsDataFrame-class"} object. If class centres and variances are not provided, it first fits a multinomial logistic regression model (\code{spmultinom}), then predicts the class centres and variances based on the output from the \code{nnet::multinom}.
#'
#' @param formulaString formula.
#' @param observations SpatialPointsDataFrame.
#' @param covariates SpatialPixelsDataFrame.
#' @param class.c class centers (per variable).
#' @param class.sd class standard deviation (per variable).
#' @param fuzzy.e fuzzy coefficient.
#'
#' @return A fuzzy kmeans model
#' @export
#'
#' @author \href{https://opengeohub.org/people/tom-hengl}{Tom Hengl}
#'
#' @examples
#' library(plotKML)
#' library(sp)
#'
#' data(eberg)
#' # subset to 20%:
#' eberg <- eberg[runif(nrow(eberg))<.2,]
#' data(eberg_grid)
#' coordinates(eberg) <- ~X+Y
#' proj4string(eberg) <- CRS("+init=epsg:31467")
#' gridded(eberg_grid) <- ~x+y
#' proj4string(eberg_grid) <- CRS("+init=epsg:31467")
#' # derive soil predictive components:
#' eberg_spc <- spc(eberg_grid, ~PRMGEO6+DEMSRT6+TWISRT6+TIRAST6)
#' # predict memberships:
#' formulaString = soiltype ~ PC1+PC2+PC3+PC4+PC5+PC6+PC7+PC8+PC9+PC10
#' eberg_sm <- spfkm(formulaString, eberg, eberg_spc@predicted)
#' \donttest{
#' # plot memberships:
#'   pal = seq(0, 1, 1/50)
#'   spplot(eberg_sm@mu, col.regions=grey(rev(pal)))
#' }
setMethod("spfkm", signature(formulaString = "formula", observations = "SpatialPointsDataFrame", covariates = "SpatialPixelsDataFrame"), function(formulaString, observations, covariates, class.c = NULL, class.sd = NULL, fuzzy.e = 1.2){

  ## generate formula if missing:
  if(missing(formulaString)) {
    formulaString <- stats::as.formula(paste(names(observations)[1], "~", paste(names(covariates), collapse="+"), sep=""))
  }
  ## check the formula string:
  if(!plyr::is.formula(formulaString)){
      stop("'formulaString' object of class 'formula' required")
  }

  ## selected variables:
  tv = all.vars(formulaString)[1]
  sel = names(covariates) %in% all.vars(formulaString)[-1]
  if(all(sel==FALSE)|length(sel)==0){
      stop("None of the covariates in the 'formulaString' matches the column names in the 'covariates' object")
  }

  ## if available, use class centres:
  check_tc <- !is.null(class.c)&!is.null(class.sd)
  if(check_tc){
    if(!class(class.c)=="matrix"){ stop("Object of type 'matrix' with column names for covariates and row names correspodning to the class names required") }
    if(!class(class.sd)=="matrix"){ stop("Object of type 'matrix' with column names for covariates and row names correspodning to the class names required") }
    mout = list(NULL)
  }
  ## otherwise, estimate class centres using the multinomial logistic regression:
  else {
    message("Trying to estimate the class centres using the 'multinom' method...")
    ## multinomial logistic regression:
    rout <- spmultinom(formulaString=formulaString, observations, covariates, class.stats=TRUE, predict.probs=FALSE)
    mout = rout$model
    if(length(unique(rout$fit))<2){ stop("Predictions resulted in <2 classes. See ?multinom for more info") }
    class.c = rout$class.c
    class.sd = rout$class.sd
  }

  cl <- as.list(row.names(class.c))
  dsf <- NULL
  ## derive distances in feature space:
  for(c in unlist(cl)){
      dsf[[c]] <- data.frame(lapply(names(covariates)[sel], FUN=function(x){rep(NA, length(covariates@data[,1]))}))
      names(dsf[[c]]) <- names(covariates)[sel]
      for(j in names(covariates)[sel]){
         dsf[[c]][,j] <- ((covariates@data[,j]-class.c[c,j])/class.sd[c,j])^2
      }
  }
  ## sum up distances per class:
  ds <- NULL
  ds <- lapply(dsf, FUN=function(x){sqrt(rowSums(x, na.rm=TRUE, dims=1))})
  names(ds) <- unlist(cl)
  ds <- data.frame(ds)
  ## total sum:
  tt <- rowSums(ds^(-2/(fuzzy.e-1)), na.rm=TRUE, dims=1)
  ## derive the fuzzy membership:
  mm <- covariates[1]
  for(c in unlist(cl)){
    mm@data[,c] <- (ds[,c]^(-2/(fuzzy.e-1))/tt)
  }
  mm@data[,names(covariates)[1]] <- NULL

  ## Derive the dominant class:
  maxm <- sapply(data.frame(t(as.matrix(mm@data))), FUN=function(x){max(x, na.rm=TRUE)})
  ## class having the highest membership
  cout <- NULL
  for(c in unlist(cl)){
       cout[which(mm@data[,c] == maxm)] <- c
  }
  cout <- as.factor(cout)

  ## construct a map:
  pm <- covariates[1]
  pm@data[,tv] <- cout
  pm@data[,names(covariates)[1]] <- NULL

  ## overlay observations and covariates:
  ov <- sp::over(observations, pm)
  sel.c <- !is.na(ov[,tv]) & !is.na(observations@data[,tv])

  ## kappa statistics:
  if(requireNamespace("mda", quietly = TRUE)&requireNamespace("psych", quietly = TRUE)){
    cf <- mda::confusion(ov[sel.c,tv], as.character(observations@data[sel.c,tv]))
    ## remove missing classes:
    a <- attr(cf, "dimnames")[[1]] %in% attr(cf, "dimnames")[[2]]
    b <- attr(cf, "dimnames")[[2]] %in% attr(cf, "dimnames")[[1]]
    c.kappa = psych::cohen.kappa(cf[a,b])
    message(paste("Estimated Cohen Kappa (weighted):", signif(c.kappa$weighted.kappa, 4)))
  } else {
    cf <- NULL
  }

  ## create the output object:
  out <- methods::new("SpatialMemberships", predicted = pm, model = mout, mu = mm, class.c = class.c, class.sd = class.sd, confusion = cf)
  return(out)

})
