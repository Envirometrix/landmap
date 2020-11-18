
setClass("spLearner", slots = c(spModel = "ANY", vgmModel = "list", covariates = "SpatialPixelsDataFrame", spID = "SpatialGridDataFrame", quantregModel = "ANY"), validity = function(object) {
    if(!class(object@vgmModel$observations)=="SpatialPointsDataFrame")
      return("Expecting an object of class 'SpatialPointsDataFrame'")
    cn = c("cov.model", "lambda", "practicalRange")
    if(!all(cn %in% names(object@vgmModel$vgm))){
      x <- cn[!(cn %in% names(object@vgmModel$vgm))]
      return(paste("Missing column names:", paste(x, collapse=", ")))
    }
})

setClass("SpatialComponents", representation (predicted = "SpatialPixelsDataFrame", pca = "list"), validity = function(object) {
   cnames <- attr(object@pca$rotation, "dimnames")[[1]]
   pnames <- attr(object@pca$rotation, "dimnames")[[2]]
   if(!length(object@pca$sdev)==length(cnames)|!length(object@pca$sdev)==length(pnames))
      return("Number of components of the 'sdev' and 'rotation' objects do not match")
   # check if column names match:
   if(!all(pnames %in% names(object@predicted)))
      return("Column names in the 'predicted' slot and 'pca' slots do not match")
})

setClass("SpatialMemberships", representation (predicted = "SpatialPixelsDataFrame", model = "list", mu = "SpatialPixelsDataFrame", class.c = "matrix", class.sd = "matrix", confusion = "ANY"), validity = function(object) {
   ## check if column names match:
   if(!all(row.names(object@class.c) %in% levels(object@predicted@data[,1])))
      return("Row names in the 'class.c' slot and 'predicted' slots do not match")
   if(!all(row.names(object@class.sd) %in% levels(object@predicted@data[,1])))
      return("Row names in the 'class.sd' slot and 'predicted' slots do not match")
   if(ncol(object@mu@data)<2)
      return("A minimum of two membership maps required")
   # check if all mu's sum to 1 (plus minus 1%):
   if(!all(rowSums(object@mu@data, na.rm=TRUE)>.99&rowSums(object@mu@data, na.rm=TRUE)<1.01))
      return("Some rows in the 'mu' slot do not sum up to 1")
})

if(!isGeneric("predict")){
  setGeneric("predict", function(object, ...){standardGeneric("predict")})
}

if(!isGeneric("print")){
  setGeneric("print", function(x, ...){standardGeneric("print")})
}

if(!isGeneric("over")){
  setGeneric("over", function(x, y, ...){standardGeneric("over")})
}

if(!isGeneric("getSpatialTiles")){
  setGeneric("getSpatialTiles", function(obj, ...){standardGeneric("getSpatialTiles")})
}

if(!isGeneric("tile")){
  setGeneric("tile", function(x, ...){standardGeneric("tile")})
}

if(!isGeneric("spc")){
  setGeneric("spc", function(obj, ...){standardGeneric("spc")})
}

if(!isGeneric("spsample.prob")){
  setGeneric("spsample.prob", function(observations, covariates, ...){standardGeneric("spsample.prob")})
}

if (!isGeneric("train.spLearner")){
  setGeneric("train.spLearner", function(observations, formulaString, covariates, ...){standardGeneric("train.spLearner")})
}

if (!isGeneric("tune.spLearner")){
  setGeneric("tune.spLearner", function(object, ...){standardGeneric("tune.spLearner")})
}

if (!isGeneric("fit.vgmModel")){
  setGeneric("fit.vgmModel", function(formulaString.vgm, rmatrix, predictionDomain, ...){standardGeneric("fit.vgmModel")})
}

if (!isGeneric("spfkm")){
  setGeneric("spfkm", function(formulaString, observations, covariates, ...){standardGeneric("spfkm")})
}

if (!isGeneric("sample.grid")){
  setGeneric("sample.grid", function(obj, cell.size, n, ...){standardGeneric("sample.grid")})
}

if(!isGeneric("buffer.dist")){
  setGeneric("buffer.dist", function(observations, predictionDomain, ...){standardGeneric("buffer.dist")})
}

if (!isGeneric("spmultinom")){
  setGeneric("spmultinom", function(formulaString, observations, covariates, ...){standardGeneric("spmultinom")})
}
