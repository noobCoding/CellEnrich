#' @export

pathwayPvalue <- function(GroupInfo, pres, pres2) {
  res <- c()
  Cells <- unique(GroupInfo)
  total <- length(GroupInfo)

  for (i in 1:length(Cells)) {
    thisCell <- Cells[i]
    thisCellIdx <- which(GroupInfo == thisCell)

    k <- length(thisCellIdx)

    thisCellPathways <- table(unlist(pres[thisCellIdx]))
    if(nrow(thisCellPathways)<1){next}
    pv <- c()

    for (j in 1:length(thisCellPathways)) {
      thisPathway <- names(thisCellPathways)[j]
      q <- unname(thisCellPathways)[j] # selected white ball
      m <- pres2[names(genesets)[as.numeric(thisPathway)]] # total white ball
      pv[j] <- round(1 - phyper(q - 1, m, total - m, k), 4)
    }
    names(pv) <- names(genesets)[as.numeric(names(thisCellPathways))]

    res <- rbind(res, cbind(thisCell, names(pv), unname(pv)))
  }
  res <- data.frame(res, stringsAsFactors = FALSE)
  colnames(res) <- c("Cell", "Geneset", "Qvalue")

  res$Cell <- as.character(res$Cell)
  res$Geneset <- as.character(res$Geneset)

  res$Qvalue <- round(p.adjust(as.numeric(as.character(res$Qvalue)), "fdr"), 4)

  return(res)
}
