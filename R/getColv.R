getColv <- function(GroupInfo){
  Cells <- unique(sort(GroupInfo))

  UniqueCol <- briterhex(scales::hue_pal()(length(Cells)))
  names(UniqueCol) <- Cells

  x <- c()
  y <- c()
  for (i in 1:length(Cells)) {
    x[i] <- Cells[i]
    y[i] <- length(which(GroupInfo == Cells[i]))
  }

  colV <- unname(UniqueCol[x])
  names(colV) <- Cells
  return(colV)
}
