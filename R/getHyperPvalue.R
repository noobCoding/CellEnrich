getHyperPvalue <- function(genes, genesets, A, lgs, q0) {

  lg <- length(genes)

  if(lg == 0){return(integer(0))}
  gidx <- 1:length(genes)
  names(gidx) <- genes

  # ------ define smaller biobj

  biobj <- matrix(0,length(genes), length(genesets))
  for (i in 1:length(genesets)) {
    biobj[unname(gidx[genesets[[i]]]), i] <- 1
  }

  rownames(biobj) <- genes
  colnames(biobj) <- names(genesets)

  if(length(genes)==1){
    biobj <- biobj[genes,]
  }
  else{
    biobj <- unname(colSums(biobj))
  }

  # ------

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
