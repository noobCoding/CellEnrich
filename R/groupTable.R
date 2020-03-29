groupTable <- function(pres, genesets, dfobj, pres2) {
  # for pres2
  genesetIdx <- sapply(names(pres2), function(i) {
    which(i == names(genesets))
  }, USE.NAMES = FALSE)
  pres2Idx <- pres2
  names(pres2Idx) <- genesetIdx

  groups <- sort(as.character(unique(dfobj$col)))
  res <- data.frame(stringsAsFactors = FALSE)

  tot <- sum(pres2Idx)

  for (i in 1:length(groups)) {
    pathways <- table(unlist(pres[which(dfobj$col == groups[i])]))
    if(length(pathways)<1){next}
    # what genesets are enriched per each group.

    k <- sum(pathways) # selected ball

    gt <- sapply(1:length(pathways), function(j) {
      q <- pathways[j] # selected white ball, 1
      m <- unname(pres2Idx[names(pathways[j])]) # total white ball, 28
      # n <- tot - m # total black ball
      round(1 - phyper(q - 1, m, tot - m, k), 4)
    })
    gt <- gt[which(gt < 0.25)] # pvalue 0.25
    if(length(gt)){
      res <- rbind(res, cbind(groups[i], names(gt), unname(gt)))
    }
  }
  colnames(res) <- c("groups", "genesetidx", "pvalue")

  res$groups <- as.character(res$groups)
  res$genesetidx <- as.numeric(as.character(res$genesetidx))
  res$genesetidx <- sapply(res$genesetidx, function(i) {
    names(genesets)[i]
  })
  res$pvalue <- as.numeric(as.character(res$pvalue))

  return(res)
}
