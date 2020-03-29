getHyperPvalue <- function(genes, genesets, A, lgs, q0, biobj) {
  lg <- length(genes)
  if(length(genes)==1){
    biobj <- biobj[genes,]
  }
  else{
    biobj <- unname(colSums(biobj[genes, ]))
  }

  pv <- sapply(1:length(genesets), function(i) {
    q <- biobj[i] # selected white ball
    m <- lgs[i] # white ball
    n <- A - m # black ball
    k <- lg # selected ball
    1 - phyper(q - 1, m, n, k)
  })
  # names(pv) <- names(genesets)
  pv <- p.adjust(pv, "fdr")
  return(which(pv < q0))
}
