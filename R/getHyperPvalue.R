#' @title Hypergeometric pvalue
#'
#' @description calculate hypergeometric pvalue for given genes
#'
#'
#' @param genes genes
#'
#' @export
#'
#'

getHyperPvalue <- function(genes, genesets) {
  A <- length(unique(unlist(genesets)))
  pv <- rep(0, length(genesets))

  for (i in 1:length(genesets)) {
    q <- length(intersect(genesets[[i]], genes))
    m <- length(genesets[[i]])
    n <- A - m
    k <- length(genes)
    pv[i] <- 1 - phyper(q - 1, m, n, k)
  }
  names(pv) <- names(genesets)
  return(pv)
}
