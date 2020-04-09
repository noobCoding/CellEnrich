getCellPlot = function(dfobj, Cells){

  cat('getCellPlot\n')
  require(highcharter)

  colnames(dfobj) <- c("x", "y", "col")
  dfobj <<- dfobj


  source('R/briterhex.R')

  UniqueCol <- briterhex(scales::hue_pal()(length(Cells)))
  names(UniqueCol) <- Cells

  dfobj$col <- unname(UniqueCol[dfobj$col])

  hc <- hchart(
    dfobj,
    type = 'scatter',
    hcaes(x = x, y = y, color = col)
  ) %>% hc_tooltip(FALSE)

  cat('\n')

  return(hc)

}
