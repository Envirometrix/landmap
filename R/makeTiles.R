#' Make a tiling system from a bounding box
#'
#' @param bb Bounding Box
#' @param block.x Size of the block in X
#' @param block.y Size of the block in Y
#' @param overlap.percent Percent of overlap; default 0
#' @param limit.bbox Optional limiting bounding box
#' @param columns Optional number of columns
#' @param rows Optional number of rows
#'
#' @return A regular tiling system
#' @export
#'
#' @author \href{https://opengeohub.org/people/tom-hengl}{Tom Hengl}
makeTiles <- function(bb, block.x, block.y, overlap.percent, limit.bbox, columns = NULL, rows = NULL){

  ## number of tiles:
  xn = ceiling(signif(diff(bb[1,]),5)/block.x)
  yn = ceiling(signif(diff(bb[2,]),5)/block.y)

  # number of tiles:
  message(paste("Generating", xn*yn, "tiles..."))
  xminl = bb[1,1]
  yminl = bb[2,1]
  xmaxl = bb[1,1] + (xn-1) * block.x
  ymaxl = bb[2,1] + (yn-1) * block.y
  xminu = bb[1,1] + block.x
  yminu = bb[2,1] + block.y
  xmaxu = bb[1,1] + xn * block.x
  ymaxu = bb[2,1] + yn * block.y

  b.l <- expand.grid(KEEP.OUT.ATTRS=FALSE, xl=seq(xminl, xmaxl, by=block.x), yl=seq(yminl, ymaxl, by=block.y))
  b.u <- expand.grid(KEEP.OUT.ATTRS=FALSE, xu=seq(xminu, xmaxu, by=block.x), yu=seq(yminu, ymaxu, by=block.y))
  btiles <- cbind(b.l, b.u)
  # expand if necessary:
  btiles$xl <- btiles$xl - block.x * overlap.percent/100
  btiles$yl <- btiles$yl - block.y * overlap.percent/100
  btiles$xu <- btiles$xu + block.x * overlap.percent/100
  btiles$yu <- btiles$yu + block.y * overlap.percent/100

  if(limit.bbox == TRUE){
    ## fix min max coordinates:
    btiles$xl <- ifelse(btiles$xl < bb[1,1], bb[1,1], btiles$xl)
    btiles$yl <- ifelse(btiles$yl < bb[2,1], bb[2,1], btiles$yl)
    btiles$xu <- ifelse(btiles$xu > bb[1,2], bb[1,2], btiles$xu)
    btiles$yu <- ifelse(btiles$yu > bb[2,2], bb[2,2], btiles$yu)
  }

  # add offset for rgdal (optional):
  if(!is.null(columns)&!is.null(rows)){
    btiles$offset.y <- round(rows*(bb[2,2]-btiles$yu)/(bb[2,2]-bb[2,1]))
    btiles$offset.x <- columns + round(columns*(btiles$xl-bb[1,2])/(bb[1,2]-bb[1,1]))
    btiles$region.dim.y <- round(rows*(btiles$yu-btiles$yl)/(bb[2,2]-bb[2,1]))
    btiles$region.dim.x <- round(columns*(btiles$xu-btiles$xl)/(bb[1,2]-bb[1,1]))
  }

  return(btiles)

}
