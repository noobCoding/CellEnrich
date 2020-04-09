getCellHistogram <- function(GroupInfo, colV){
  cat('getCellHistogram\n')
  #require(ggplot2)
  require(highcharter)
  Cells <- unique(sort(GroupInfo))

  x <- c()
  y <- c()
  for (i in 1:length(Cells)) {
    x[i] <- Cells[i]
    y[i] <- length(which(GroupInfo == Cells[i]))
  }

  hc <- highchart() %>%
    hc_chart(type = 'column', legend = list(enabled  = FALSE)) %>%
    hc_title(text = 'Groups') %>%
    hc_xAxis(categories = x) %>%
    hc_plotOptions(grouping = FALSE) %>%
    hc_add_series(data = y, colorByPoint = TRUE, showInLegend = FALSE, name = 'Count') %>%
    hc_colors(colV) %>%
    hc_exporting(enabled = TRUE)
  return(hc)

}
