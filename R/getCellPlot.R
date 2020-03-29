getCellPlot = function(dfobj, Cells){
  colnames(dfobj) <- c("x", "y", "col")
  dfobj <<- dfobj

  source('R/briterhex.R')

  UniqueCol <- briterhex(scales::hue_pal()(length(Cells)))
  names(UniqueCol) <- Cells
  colV <- unname(UniqueCol[dfobj$col])
  return(
    ggplot(dfobj, aes(x = x, y = y)) +
      geom_point(colour = colV)
  )
}
