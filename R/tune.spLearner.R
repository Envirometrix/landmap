#' Optimize spLearner by fine-tuning parameters and running feature selection
#'
#' @rdname tune.spLearner-methods
#' @aliases tune.spLearner tune.spLearner,spLearner-method
#'
#' @param object spLearner object (unoptimized),
#' @param num.trees number of random forest trees,
#' @param blocking blocking columns,
#' @param discrete_ps settings for random forest,
#' @param rdesc resampling method for fine-tuning,
#' @param inner resampling method for feature selection,
#' @param maxit maximum number of iterations for feature selection,
#' @param xg.model_Params xgboost parameter set,
#' @param xg.skip logical, should the tuning of the XGboost should be skipped?
#' @param parallel Initiate parallel processing,
#' @param hzn_depth specify whether horizon depth available in the training dataframe,
#' @param ... other arguments that can be passed on to \code{mlr::makeStackedLearner},
#'
#' @return optimized object of type spLearner
#'
#' @author \href{https://opengeohub.org/people/tom-hengl}{Tom Hengl}
#'
#' @note Currently requires that two base learners are \code{regr.ranger} and
#' \code{regr.xgboost}, and that there are at least 3 base learners in total.
#' Fine-tuning and feature selection can be quite computational
#' and it is highly recommended to start with smaller subsets of data and then measure
#' processing time. The function \code{mlr::makeFeatSelWrapper} can result in
#' errors if the covariates have a low variance or follow a zero-inflated distribution.
#' Reducing the number of features via feature selection and fine-tuning of the
#' Random Forest mtry and XGboost parameters, however, can result in
#' significantly higher prediction speed and accuracy.
#'
#' @export
#'
#' @examples
#' \donttest{
#' library(mlr)
#' library(ParamHelpers)
#' library(geoR)
#' library(xgboost)
#' library(kernlab)
#' library(ranger)
#' library(glmnet)
#' library(boot)
#' library(raster)
#' demo(meuse, echo=FALSE)
#' ## Regression:
#' sl = c("regr.ranger", "regr.xgboost", "regr.ksvm", "regr.cvglmnet")
#' m <- train.spLearner(meuse["lead"], covariates=meuse.grid[,c("dist","ffreq")],
#'       lambda=0, parallel=FALSE, SL.library=sl)
#' summary(m@spModel$learner.model$super.model$learner.model)
#' ## Optimize model:
#' m0 <- tune.spLearner(m, xg.skip = TRUE, parallel=FALSE)
#' summary(m0@spModel$learner.model$super.model$learner.model)
#' }
setMethod("tune.spLearner", signature(object = "spLearner"), function(object, num.trees = 85, blocking, discrete_ps, rdesc = mlr::makeResampleDesc("CV", iters = 2L), inner = mlr::makeResampleDesc("Holdout"), maxit = 20, xg.model_Params, xg.skip = FALSE, parallel = "multicore", hzn_depth = FALSE, ...){
  if(all(c("regr.ranger", "regr.xgboost") %in% object@spModel$learner.model$super.model$features) & object@spModel$learner.model$super.model$learner$id=="regr.lm"){
    rf.no = which(object@spModel$learner.model$super.model$features %in% "regr.ranger")
    xg.no = which(object@spModel$learner.model$super.model$features %in% "regr.xgboost")
    min.mtry = round(sqrt(length(object@spModel$features))/3)
    max.mtry = ifelse(length(object@spModel$features)>10, length(object@spModel$features)-4, length(object@spModel$features)-1)
    if(missing(discrete_ps)){
      discrete_ps <- ParamHelpers::makeParamSet(ParamHelpers::makeDiscreteParam("mtry", values = unique(round(seq(min.mtry, max.mtry, length.out = 10)))))
    }
    if(missing(xg.model_Params)){
      xg.model_Params <- ParamHelpers::makeParamSet(
        ParamHelpers::makeDiscreteParam("nrounds", value=c(10,50)),
        ParamHelpers::makeDiscreteParam("max_depth", value=c(1,2,3,4)),
        ParamHelpers::makeDiscreteParam("eta", value=c(0.3,0.4,0.5,0.6)),
        ParamHelpers::makeDiscreteParam("subsample", value=c(1)),
        ParamHelpers::makeDiscreteParam("min_child_weight", value=c(1)),
        ParamHelpers::makeDiscreteParam("colsample_bytree", value=c(0.8))
      )
    }
    ctrl = mlr::makeTuneControlGrid()
    ctrlF = mlr::makeFeatSelControlRandom(maxit = maxit)
    ## regression matrix
    t.var = all.vars(object@spModel$learner.model$super.model$learner.model$terms)[1]
    rm.df = object@vgmModel$observations@data[,c(t.var, object@spModel$features)]
    if(missing(blocking)){ blocking <- as.factor(sp::over(object@vgmModel$observations, object@spID)[,1]) }
    coordinates <- as.data.frame(object@vgmModel$observations@coords)
    tsk0 <- mlr::makeRegrTask(data = rm.df, target = t.var, blocking = blocking, coordinates = coordinates)
    if(parallel=="multicore"){
      parallelMap::parallelStartSocket(parallel::detectCores())
    }
    ## fine-tune mtry
    message("Running tuneParams for ranger... ", immediate. = TRUE)
    resR.lst = mlr::tuneParams(mlr::makeLearner("regr.ranger", num.threads = round(parallel::detectCores()/length(discrete_ps$pars$mtry$values)), num.trees=num.trees), task = tsk0, resampling = rdesc, par.set = discrete_ps, control = ctrl)
    ## Error: mtry can not be larger than number of variables in data.
    if(resR.lst$x$mtry >= length(object@spModel$features)){
      lrn.rf = mlr::makeLearner("regr.ranger", num.threads = parallel::detectCores(), num.trees=num.trees, importance="impurity")
    } else {
      lrn.rf = mlr::makeLearner("regr.ranger", num.threads = parallel::detectCores(), mtry=resR.lst$x$mtry, num.trees=num.trees, importance="impurity")
    }
    ## feature selection
    lrn1 <- mlr::makeFeatSelWrapper(lrn.rf, resampling = inner, control = ctrlF, show.info=TRUE)
    try( var.mod1 <- mlr::train(lrn1, task = tsk0), silent = TRUE)
    if(!class(.Last.value)[1]=="try-error"){
      var.sfeats1 <- mlr::getFeatSelResult(var.mod1)
    } else {
      var.sfeats1 <- data.frame(x=object@spModel$features)
    }
    ## fine-tune xgboost
    lrn.xg = mlr::makeLearner("regr.xgboost", par.vals = list(objective ='reg:squarederror'))
    if(xg.skip == FALSE){
      message("Running tuneParams for xgboost... ", immediate. = TRUE)
      resX.lst = mlr::tuneParams(mlr::makeLearner("regr.xgboost"), task = tsk0, resampling = rdesc, par.set = xg.model_Params, control = ctrl)
      lrn.xg = mlr::setHyperPars(lrn.xg, par.vals = resX.lst$x)
      lrn2 = mlr::makeFeatSelWrapper(lrn.xg, resampling = inner, control = ctrlF, show.info=TRUE)
      try( var.mod2 <- mlr::train(lrn2, task = tsk0), silent = TRUE)
      if(!class(.Last.value)[1]=="try-error"){
        var.sfeats2 <- mlr::getFeatSelResult(var.mod2)
      } else {
        var.sfeats2 <- data.frame(x=object@spModel$features)
      }
      sel.pars = c(var.sfeats1$x, var.sfeats2$x)
    } else {
      sel.pars = var.sfeats1$x
    }
    ## new shorter formula
    ## we add depth otherwise not a 3D model (for S it was removed)
    if(hzn_depth==TRUE){
      formulaString.y = stats::as.formula(paste(t.var, ' ~ ', paste(unique(c("hzn_depth", sel.pars)), collapse="+")))
    } else {
      formulaString.y = stats::as.formula(paste(t.var, ' ~ ', paste(unique(sel.pars), collapse="+")))
    }
    message("Fiting Stacked Ensemble Model... ", immediate. = TRUE)
    ## final EML model
    tskF <- mlr::makeRegrTask(data = rm.df[,all.vars(formulaString.y)], target = t.var, blocking = blocking, coordinates = coordinates)
    lrnsE <- c(list(lrn.rf, lrn.xg), lapply(object@spModel$learner.model$super.model$features[-c(rf.no, xg.no)], function(i){mlr::makeLearner(i)}))
    init.m <- mlr::makeStackedLearner(base.learners = lrnsE, predict.type = "response", method = "stack.cv", super.learner = "regr.lm", resampling=mlr::makeResampleDesc(method = "CV", blocking.cv=TRUE), ...)
    suppressWarnings( t.m <- mlr::train(init.m, tskF) )
    message("Fitting a quantreg model using 'ranger::ranger'...", immediate. = TRUE)
    quantregModel <- ranger::ranger(t.m$learner.model$super.model$learner.model$terms, t.m$learner.model$super.model$learner.model$model, num.trees=85, importance="impurity", quantreg=TRUE)
    if(parallel=="multicore"){
      parallelMap::parallelStop()
    }
    out <- methods::new("spLearner", spModel = t.m, vgmModel = object@vgmModel, covariates = object@covariates, spID = object@spID, quantregModel = quantregModel)
    return(out)
  } else {
    stop("spLearner model with 'regr.lm' super model and learners 'regr.ranger' and 'regr.xgboost' required")
  }
})
