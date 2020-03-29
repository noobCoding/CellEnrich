getCellHistogram <- function(GroupInfo, colV){
  require(ggplot2)

  Cells <- unique(sort(GroupInfo))

  x <- c()
  y <- c()
  for (i in 1:length(Cells)) {
    x[i] <- Cells[i]
    y[i] <- length(which(GroupInfo == Cells[i]))
  }

  ggobjdf <- data.frame(x, y, stringsAsFactors = FALSE)
  colnames(ggobjdf) <- c("x", "y")

  return(
    ggplot(ggobjdf, aes(x = x, y = y)) +
      geom_bar(stat = "identity", fill = colV)
  )


}
