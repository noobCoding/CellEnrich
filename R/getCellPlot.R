getCellPlot = function(dfobj, Cells){

  cat('getCellPlot\n')
  require(ggplot2)
  #require(highcharter)

  colnames(dfobj) <- c("x", "y", "col")
  dfobj <<- dfobj


  source('R/briterhex.R')

  UniqueCol <- briterhex(scales::hue_pal()(length(Cells)))
  names(UniqueCol) <- Cells

  colV <- unname(UniqueCol[dfobj$col])

  cat('\n')
  return(
    ggplot(dfobj, aes(x = x, y = y)) +
      geom_point(colour = colV)
  )

  # highcharter cancel

  #dfobj$col <- unname(UniqueCol[dfobj$col])
  #rownames(dfobj) = NULL

  #hchart(dfobj, type = 'scatter', hcaes(x = x, y = y, color = col)) %>%
    #hc_add_series(data = dfobj[3000,nrow(dfobj)], colorByPoint = TRUE, showInLegend = FALSE) %>%
    #hc_colors(colV)
    #hc_tooltip(FALSE) %>%
    #hc_exporting(enabled = TRUE)

  #return(hc)

}
