#' @export

GeneFlush = function(genes, genesets){
  gsgenes <- unique(unlist(genesets))
  remgenes <- sapply(setdiff(genes, gsgenes), function(i) { which(i == genes) }, USE.NAMES = FALSE)
  return(remgenes)
}
