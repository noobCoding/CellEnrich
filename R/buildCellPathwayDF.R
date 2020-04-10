buildCellPathwayDF <- function(GroupInfo, pres, genesets){
  cat('buildCellPathwayDF\n')
  Cells <- unique(GroupInfo)
  CellPathwayDF <- data.frame(stringsAsFactors = FALSE)

  for (i in 1:length(Cells)) {
    thisCell <- Cells[i]
    tt <- table(unlist(pres[which(thisCell == GroupInfo)]))

    if(nrow(tt)) { CellPathwayDF <- rbind(CellPathwayDF, cbind(thisCell, names(tt), unname(tt))) }
  }

  colnames(CellPathwayDF) <- c("Cell", "Geneset", "Count")

  CellPathwayDF$Cell <- as.character(CellPathwayDF$Cell)
  CellPathwayDF$Geneset <- names(genesets)[as.numeric(as.character(CellPathwayDF$Geneset))]
  CellPathwayDF$Count <- as.numeric(as.character(CellPathwayDF$Count))

  # ------ add length column

  # Length <- getlgs(CellPathwayDF$Geneset)
  Length <- getlgs(genesets[as.character(CellPathwayDF$Geneset)])
  CellPathwayDF <- cbind(CellPathwayDF, Length)

  # ------ select genesets with count > 1

  CellPathwayDF <- CellPathwayDF %>%
    dplyr::filter(Count > 1)

  return(CellPathwayDF)
}
