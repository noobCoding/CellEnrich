GenesetFlush <- function(genes, genesets){
  cat('GenesetFlush\n')
  for (i in 1:length(genesets)) {
    genesets[[i]] <- intersect(genesets[[i]], genes)
  }
  return(genesets)
}
