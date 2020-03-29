getbiobj <- function(genes, genesets) {
  gidx <- 1:length(genes)
  names(gidx) <- genes

  res <- matrix(0, length(genes), length(genesets))
  for (i in 1:length(genesets)) {
    res[unname(gidx[genesets[[i]]]), i] <- 1
  }

  rownames(res) <- genes
  colnames(res) <- names(genesets)
  return(res)
}
