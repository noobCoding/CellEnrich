#' @export

emphasize <- function(path = FALSE, inputObj, dfobj, Cells, pres, genesets) {
  cat('emphasize\n')
  buildRlobj <- function(items) {

    rlobj <- data.frame(stringsAsFactors = FALSE)

    for (i in 1:length(items)) {
      kk <- strsplit(items[[i]], " @")[[1]]
      name <- kk[1]
      location <- kk[2]
      rlobj <- rbind(rlobj, cbind(name, location))
    }

    colnames(rlobj) <- c("name", "location")
    rlobj$name <- as.character(rlobj$name)
    rlobj$location <- as.character(rlobj$location)

    return(rlobj)
  }

  getCellValues <- function(rlobj) {

    ret <- list()
    for (i in 1:nrow(rlobj)) {
      thisGeneset <- which(names(genesets) == rlobj[i, 1]) # index
      thisGroup <- rlobj[i, 2]

      thisCellsIdx <- which(dfobj$col == thisGroup)
      if(length(thisCellsIdx)==0){
        ret[[i]] <- c()
        next
      }
      rn <- thisCellsIdx
      res <- c()
      for (j in 1:length(rn)) {
        if (thisGeneset %in% pres[[rn[j]]]) {
          res <- c(res, rn[j])
        }
      }
      if (length(res) == 0) {
        rlobj <- rlobj[-i, ]
        next
      }
      ret[[i]] <- res
      names(ret)[i] <- thisGroup
    }

    #names(ret) <- rlobj$location
    return(ret)
  }

  rlobj <- buildRlobj(inputObj) # split into name, location dataframe

  cellValues <- getCellValues(rlobj) # get cell index for each cell

  dfobj_new <- data.frame(dfobj, stringsAsFactors = FALSE)

  colnames(dfobj_new) <- c("x", "y", "col")

  # define ggobj2 element

  UniqueCol <- briterhex(scales::hue_pal()(length(Cells)))
  names(UniqueCol) <- Cells

  colV <- unname(UniqueCol[dfobj_new$col])

  colV[-unlist(cellValues, use.names = FALSE)] <- "#95A5A6" # gray color

  dfobj_new$col <- colV
  rownames(dfobj_new) = NULL

  hc <- hchart(
      dfobj_new,
      type = 'scatter',
      hcaes(x = x, y = y, color = col)
    ) %>%
    hc_tooltip(FALSE) %>%
    hc_exporting(enabled = TRUE)

  if (path) { # add mean point to path
    dfobj_path = data.frame()
    for (i in 1:length(cellValues)) {
      x <- mean(as.numeric(dfobj_new$x[cellValues[[i]]]))
      y <- mean(as.numeric(dfobj_new$y[cellValues[[i]]]))
      dfobj_path <- rbind(dfobj_path, c(x, y))
    }
    colnames(dfobj_path) = c('x','y')

    hc <- hc %>%
      hc_add_series(dfobj_path, 'line', hcaes(x = x, y = y), showInLegend = FALSE) %>%
      hc_plotOptions(
        line = list(
          marker = FALSE,
          dashStyle = 'ShortDash',
          lineColor = '#5f27cd',
          lineWidth = 2
        ))

    #newIdx <- (nrow(dfobj) + 1):nrow(dfobj_new)
    #cellValues <- c(unname(unlist(cellValues)), newIdx)

    #dfobj_new$x <- round(as.numeric(dfobj_new$x), 4)
    #dfobj_new$y <- round(as.numeric(dfobj_new$y), 4)

    #for (i in 1:(length(newIdx) - 1)) { # add curve
      #newCurve <- paste(
        #" + geom_curve( aes(x = ", "x[newIdx[", i,
        #"]], y = y[newIdx[", i, "]], xend = x[newIdx[", i + 1,
        #"]], yend = y[newIdx[", i + 1, ']]), size = 0.5, linetype = "longdash",',
        #"curvature = 0.1, colour = '#000000', ", 'arrow = arrow(length = unit(0.1,"inches")))'
        # "curvature = 0.1, colour = colV[newIdx[", i, ']], arrow = arrow(length = unit(0.1,"inches")))'
      #)
      #graphString <- paste(graphString, newCurve, sep = "")
    #}
  }

  #  eval(parse(text = graphString))
  return(hc)
}
