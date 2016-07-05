#' My heatmap
#'
#' This is a cool heatmap function
#'
#' @param eset An ExpressionSet
#' @param n The top most variable genes
#' @param ... passed to heatmap
#'
#' @export
#' @import Biobase
myHm = function(eset,n,...) {
  sds = apply(exprs(eset),1,sd)
  topN = order(sds, decreasing=TRUE)[1:n]
  heatmap(exprs(eset)[topN,],...)
}
