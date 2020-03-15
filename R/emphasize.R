emphasize <- function(path = FALSE, inputObj, dfobj, Cells, pres) {
  #print('emphasize : ')
  #print(inputObj)

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
    }

    names(ret) <- rlobj$location
    return(ret)
  }

  rlobj <- buildRlobj(inputObj) # split into name, location dataframe

  cellValues <- getCellValues(rlobj) # get cell index for each cell

  dfobj_new <- data.frame(dfobj, stringsAsFactors = FALSE)

  colnames(dfobj_new) <- c("x", "y", "col")

  # define ggobj2 element

  UniqueCol <- briterhex(scales::hue_pal()(length(Cells)))
  names(UniqueCol) <- Cells

  colV <- unname(UniqueCol[dfobj$col])

  colV[-unlist(cellValues, use.names = FALSE)] <- "#95A5A6" # gray color

  graphString <- "ggobj2 <- ggplot(dfobj_new, aes(x = x, y = y)) + geom_point(colour = colV)"

  if (path) { # add mean point to path

    for (i in 1:length(cellValues)) {
      x <- mean(as.numeric(dfobj_new$x[cellValues[[i]]]))
      y <- mean(as.numeric(dfobj_new$y[cellValues[[i]]]))
      dfobj_new <- rbind(dfobj_new, c(x, y, "meanPoint"))
      colV <- c(colV, "#000000")
    }
    newIdx <- (nrow(dfobj) + 1):nrow(dfobj_new)
    cellValues <- c(unname(unlist(cellValues)), newIdx)

    dfobj_new$x <- round(as.numeric(dfobj_new$x), 4)
    dfobj_new$y <- round(as.numeric(dfobj_new$y), 4)

    for (i in 1:(length(newIdx) - 1)) { # add curve
      newCurve <- paste(
        " + geom_curve( aes(x = ", "x[newIdx[", i,
        "]], y = y[newIdx[", i, "]], xend = x[newIdx[", i + 1,
        "]], yend = y[newIdx[", i + 1, ']]), size = 0.5, linetype = "longdash",',
        "curvature = 0.1, colour = '#000000', ", 'arrow = arrow(length = unit(0.1,"inches")))'
        # "curvature = 0.1, colour = colV[newIdx[", i, ']], arrow = arrow(length = unit(0.1,"inches")))'
      )
      graphString <- paste(graphString, newCurve, sep = "")
    }
  }

  eval(parse(text = graphString))
  return(ggobj2)
}
