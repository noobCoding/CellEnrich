#' @title CellEnrich
#'
#'
#' @importFrom DT dataTableOutput
#'
#' @rawNamespace import(SingleCellExperiment, except = show)
#' @import Rtsne
#' @import shinyCyJS
#' @import dplyr
#' @rawNamespace import(shiny, except = dataTableOutput)
#' @import ggplot2
#' @import uwot
#' @import htmltools
#' @import magrittr
#' @import shinymaterial
#' @import waiter
#' @rawNamespace import(shinyjs, except = runExample)
#' @import scales
#' @import sortable
#' @import scran
#'
#' @export

CellEnrich <- function(CountData, CellInfo) {
  require(shinymaterial)
  require(shiny)
  require(waiter)
  require(Rtsne)
  require(uwot)
  require(ggplot2)
  require(DT)
  require(scales)
  require(sortable)
  require(scran)
  require(dplyr)

  # actionButton with style Attributes, width = NULL
  solvedButton <- function(inputId, label, style = NULL, ...) {
    value <- restoreInput(id = inputId, default = NULL)
    tags$button(
      id = inputId, style = style,
      type = "button", class = "btn btn-default action-button",
      `data-val` = value, list(label),
      ...
    )
  }

  # actionButton with Onclick Attributes
  myButton <- function(inputId, label, width = NULL, onClick = NULL, ...) {
    value <- restoreInput(id = inputId, default = NULL)
    tags$button(
      id = inputId, style = if (!is.null(width)) {
        paste0("width: ", validateCssUnit(width), ";")
      },
      type = "button", class = "btn btn-default action-button",
      onClick = onClick,
      `data-val` = value, list(label),
      ...
    )
  }


  ui <- CellEnrichUI()

  options(useFancyQuotes = FALSE)

  server <- function(input, output, session) {
    A <- ""

    ## FUNCTIONS

    briterhex <- function(colors) {
      res <- c()
      for (i in 1:length(colors)) {
        v <- as.vector(col2rgb(colors[i])) * 1.3
        v <- sapply(v, function(i) {
          min(i, 255)
        })
        res[i] <- rgb(v[1], v[2], v[3], max = 255)
      }
      return(res)
    }




    sortItem <- function(label, tableName) {
      options(useFancyQuotes = FALSE)
      paste0(
        "$('#", tableName, "')",
        ".append(", "`<div class=", "'rank-list-item'", " draggable='true'",
        " style = 'transform: translateZ(0px);'>` + ", label, " + `</div>`)"
      )
    }

    buildDT <- function(pres2) {
      DT::datatable(
        data.frame(
          Geneset = names(pres2),
          Count = as.numeric(pres2)
        ),
        options = list(
          dom = "ltp",
          lengthChange = FALSE
        ),
        rownames = FALSE,
        selection = "single"
      )
    }

    mapColor <- function(v) {
      as.factor(v)
    }

    changeCol <- function(v) {
      cols <- briterhex(scales::hue_pal()(length(unique(v))))
      uv <- unique(v)
      res <- c()
      for (i in 1:length(uv)) {
        res[which(v == uv[i])] <- cols[i]
      }
      return(res)
    }

    pathwayPvalue <- function(q0 = 0.1) {
      res <- c()
      # gs <- sapply(names(pres2), function(i) { which(names(genesets) == i) })
      Cells <- unique(CellInfo)

      total <- length(CellInfo)

      for (i in 1:length(Cells)) {
        thisCell <- Cells[i]
        thisCellIdx <- which(CellInfo == thisCell)

        k <- length(thisCellIdx)

        thisCellPathways <- table(unlist(pres[thisCellIdx]))
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


      res <- res %>% dplyr::filter(Qvalue <= q0) # 1237 * 3

      return(res)
    }

    groupTable <- function() {
      # for pres2
      genesetIdx <- sapply(names(pres2), function(i) {
        which(i == names(genesets))
      }, USE.NAMES = FALSE)
      pres2Idx <- pres2
      names(pres2Idx) <- genesetIdx

      groups <- sort(as.character(unique(dfobj$col)))
      res <- data.frame(stringsAsFactors = FALSE)

      tot <- sum(pres2Idx)

      for (i in 1:length(groups)) {
        pathways <- table(unlist(pres[which(dfobj$col == groups[i])]))
        # what genesets are enriched per each group.

        k <- sum(pathways) # selected ball

        gt <- sapply(1:length(pathways), function(j) {
          q <- pathways[j] # selected white ball, 1
          m <- unname(pres2Idx[names(pathways[j])]) # total white ball, 28
          # n <- tot - m # total black ball
          round(1 - phyper(q - 1, m, tot - m, k), 4)
        })
        gt <- gt[which(gt < 0.25)] # pvalue 0.25
        res <- rbind(res, cbind(groups[i], names(gt), unname(gt)))
      }
      colnames(res) <- c("groups", "genesetidx", "pvalue")

      res$groups <- as.character(res$groups)
      res$genesetidx <- as.numeric(as.character(res$genesetidx))
      res$genesetidx <- sapply(res$genesetidx, function(i) {
        names(genesets)[i]
      })
      res$pvalue <- as.numeric(as.character(res$pvalue))

      return(res)
    }

    getHyperPvalue <- function(genes, genesets, A) {

      lgs <- sapply(1:length(genesets), function(i) {
        length(genesets[[i]])
      })
      lg = length(genes)
      pv <- sapply(1:length(genesets), function(i) {
        q <- length(intersect(genesets[[i]], genes)) # selected white ball
        m <- lgs[i] # white ball
        n <- A - m # black ball
        k <- lg # selected ball
        1 - phyper(q - 1, m, n, k)
      })
      names(pv) <- names(genesets)
      return(pv)
    }


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
        if(length(res)==0){
          rlobj = rlobj[-i,]
          next
        }
        ret[[i]] <- res
      }


      names(ret) <- rlobj$location
      return(ret)
    }

    emphasize <- function(path = FALSE, inputObj) {
      rlobj <- buildRlobj(inputObj) # split into name, location dataframe

      cellValues <- getCellValues(rlobj) # get cell index for each cell

      dfobj_new <- data.frame(dfobj, stringsAsFactors = FALSE)

      colnames(dfobj_new) <- c("x", "y", "col")

      # define ggobj2 element

      # dfobj -> oocyte1, oocyte2

      UniqueCol <- briterhex(scales::hue_pal()(length(unique(CellInfo))))
      names(UniqueCol) <- unique(CellInfo)



      colV <- unname(UniqueCol[dfobj$col])


      colV[-unlist(cellValues, use.names = FALSE)] <- "#95A5A6" # gray color

      # colV <- changeCol(dfobj_new$col)

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

    ### CODES





    shinyjs::runjs(code = '$("#btn5Cell, #btn6Cell, #btn7Cell").prop("disabled",true)')

    # load("c2v7idx.RData")

    myf <- function(genesetnames, pres) {
      idx <- which(names(genesets) == genesetnames)
      res <- c()
      for (i in 1:length(pres)) {
        if (idx %in% pres[[i]]) {
          res <- c(res, i)
        }
      }
      return(res)
    }

    # variable initialize

    ggobj2 <- ggobj <- dtobj <- dfobj <- pres <- CellPathwayDF <- ""
    gt <- pres2 <- dfobj2 <- ggobj3 <- ggobj4 <- ""

    # start

    observeEvent(input$btn, {
      shinyjs::hide("btn")

      if (input$dropdowninput3 == "Curated") load("c2v7.RData")

      if (input$dropdowninput3 == "GeneOntology") load("c5v7.RData")

      q0 <- input$sliderinput3

      w <- Waitress$new(selector = NULL, theme = "overlay")$start()

      # geneset background flushing
      genes <- sort(rownames(CountData))
      for (i in 1:length(genesets)) {
        genesets[[i]] <- intersect(genesets[[i]], genes)
      }
      rm(genes)

      lgs <- sapply(1:length(genesets), function(i) {
        length(genesets[[i]])
      })

      # genesets = genesets[intersect(which(lgs >= 15), which(lgs <= 500))]
      # cat('minimum gene-set size :', input$sliderinput1, '\n')
      # cat('maximum gene-set size :', input$sliderinput2, '\n')

      genesets <- genesets[intersect(which(lgs >= input$sliderinput1), which(lgs <= input$sliderinput2))]
      genesets <<- genesets
      cat("minimum gene-set size :", input$sliderinput1, "\n")
      cat("maximum gene-set size :", input$sliderinput2, "\n")

      cat(length(genesets), "genesets\n")

      A <<- length(unique(unlist(genesets))) # Background genes

      shinyjs::runjs('$("form p label input").attr("disabled",true)') # radio button disable
      shinyjs::runjs("$('.shinymaterial-slider-sliderinput1').attr('disabled',true)")
      shinyjs::runjs("$('.shinymaterial-slider-sliderinput2').attr('disabled',true)")
      shinyjs::runjs("$('.shinymaterial-slider-sliderinput3').attr('disabled',true)")

      v <- CountData

      if (input$dropdowninput != "GSVA") {
        # s <- findSigGenes(v, 'median')
        s <- findSigGenes(v, input$dropdowninput)
        names(s) <- CellInfo # oocyte 1, oocyte 2
      }

      s2 <- findSigGenesGroup(v, CellInfo, q0, TopCutoff = 5)
      s2$FDR = round(as.numeric(s2$FDR),6)

      markerl1 = s2 %>% filter(Top<10)
      markerl1$Group = as.factor(markerl1$Group)

      # marker L1
      output$markerL1 = DT::renderDataTable(
        DT::datatable(markerl1,
          rownames = FALSE,
          filter = 'top',
          options = list(
            autoWidth = TRUE,
            dom = 'ltp',
            lengthChange = FALSE
          ),
          selection = 'none',
        )
      )

      # for test
      res <- list()
      for (i in 1:length(s)) {
        w$inc(1 / length(s))
        pvh <- getHyperPvalue(s[[i]], genesets, A)
        qvh <- p.adjust(pvh, "fdr")
        res[[i]] <- unname(which(qvh < q0))
      }

      w$close()
      pres <<- res
      rm(res)

      # pres : which gene-sets are significant for each cells.

      CellPathwayDF <- data.frame(stringsAsFactors = FALSE)
      Cells <- unique(CellInfo) # 16cell1, 16cell2,
      for (i in 1:length(Cells)) {
        thisCell <- Cells[i]
        tt <- table(unlist(pres[which(thisCell == CellInfo)]))
        CellPathwayDF <- rbind(CellPathwayDF, cbind(thisCell, names(tt), unname(tt)))
      }

      colnames(CellPathwayDF) <- c("Cell", "Geneset", "Count")

      CellPathwayDF$Cell <- as.character(CellPathwayDF$Cell)
      CellPathwayDF$Geneset <- names(genesets)[as.numeric(as.character(CellPathwayDF$Geneset))]
      CellPathwayDF$Count <- as.numeric(as.character(CellPathwayDF$Count))


      # select genesets with count > 1
      CellPathwayDF <- CellPathwayDF %>%
        dplyr::filter(Count > 1)

      Length <- sapply(CellPathwayDF$Geneset, function(i) {
        length(genesets[[i]])
      })
      CellPathwayDF <- cbind(CellPathwayDF, Length)


      # pres2 : for each gene-sets, how many cells are significant that gene-sets.
      pres2 <- sort(table(unlist(pres)), decreasing = T)
      names(pres2) <- names(genesets)[as.numeric(names(pres2))]
      pres2 <<- pres2

      # 2625*4
      PP <- pathwayPvalue(q0)

      CellPathwayDF <- CellPathwayDF %>%
        inner_join(PP) # 1232 * 5

      CellPathwayDF <<- CellPathwayDF

      # l2
      CellMarkers = data.frame()

      for(i in 1:length(Cells)){
        thisCell = Cells[i]
        thisCellPathways = CellPathwayDF %>% filter(Cell==thisCell) %>% select(Geneset)

        thisCellDEs = s2 %>% filter(Group==thisCell) %>% select(genes)

        tcd = thisCellDEs[,1]
        tcp = thisCellPathways[,1]
        tcp = sapply(tcp, function(i){which(names(genesets)==i)}, USE.NAMES = FALSE)
        if(length(tcp)<1){next}
        if(length(tcd)<1){next}
        tcp = table(unlist(genesets[tcp], use.names = FALSE))
        tcd = intersect(names(tcp),tcd)
        tcp = tcp[tcd]

        genes = names(tcp)
        Count = unname(tcp)
        CellMarkers = rbind(CellMarkers, data.frame(cbind(genes, Count, Group = thisCell), stringsAsFactors = FALSE))
      }

      CellMarkers = CellMarkers %>% inner_join(s2) %>% filter( Top < 10 )
      CellMarkers$Group = as.factor(CellMarkers$Group)
      CellMarkers$Count = as.numeric(CellMarkers$Count)
      CellMarkers$Group = as.factor(CellMarkers$Group)

      output$markerL2 = DT::renderDataTable(
        DT::datatable(CellMarkers,
                      rownames = FALSE,
                      filter = 'top',
                      options = list(
                        autoWidth = TRUE,
                        dom = 'ltp',
                        lengthChange = FALSE
                      ),
                      selection = 'none',
        )
      )

      # group 별 significant pathways
      # group 별 DE Genes

      # is counted


      dtobj <<- buildDT(pres2)
      # cellinfo -> oocyte1, oocyte2, ...

      if (input$dropdowninput2 == "t-SNE") {
        tsneE <- Rtsne(t(v), check_duplicates = FALSE, perplexity = 15)
        dfobj <- data.frame(tsneE$Y, col = CellInfo, stringsAsFactors = FALSE)
      }

      if (input$dropdowninput2 == "U-MAP") {
        umapE <- uwot::umap(t(v), fast_sgd = TRUE)
        dfobj <- data.frame(umapE, col = CellInfo, stringsAsFactors = FALSE)
      }

      UniqueCol <- briterhex(scales::hue_pal()(length(unique(CellInfo))))
      names(UniqueCol) <- unique(CellInfo)

      Cells <- unique(CellInfo)

      x <- c()
      y <- c()
      for (i in 1:length(Cells)) {
        x[i] <- Cells[i]
        y[i] <- length(which(CellInfo == Cells[i]))
      }

      colV <- unname(UniqueCol[x])
      names(colV) <- Cells



      ggobjdf <- data.frame(x, y, stringsAsFactors = FALSE)
      colnames(ggobjdf) <- c("x", "y")

      ggobj <<- ggplot(ggobjdf, aes(x = x, y = y)) +
        geom_bar(stat = "identity", fill = colV) # cell histogram

      #if (!is.null(ClustInfo)) {
        #ggobj3 <<- ggplot(data.frame(table(ClustInfo)), aes(x = ClustInfo, y = Freq, fill = ClustInfo)) +
          #geom_bar(stat = "identity") # cell histogram
      #}

      colnames(dfobj) <- c("x", "y", "col")
      dfobj <<- dfobj

      #if (nrow(dfobj2)) {
        #colnames(dfobj2) <- c("x", "y", "col")
        #dfobj2$col <- as.character(dfobj2$col)
        #dfobj2 <<- dfobj2
      #}

      UniqueCol <- briterhex(scales::hue_pal()(length(unique(CellInfo))))
      names(UniqueCol) <- unique(CellInfo)

      colV <- unname(UniqueCol[dfobj$col])

      # scatter plot
      ggobj2 <<- ggplot(dfobj, aes(x = x, y = y)) +
        geom_point(colour = colV) # +
      # scale_color_manual(values = briterhex(scales::hue_pal()(length(unique(dfobj$col)))))

      #if (!is.null(ClustInfo)) {
        #ggobj4 <<- ggplot(dfobj2, aes(x = x, y = y, color = col)) +
          #geom_point() +
          #scale_color_manual(values = briterhex(scales::hue_pal()(length(unique(dfobj2$col)))))
      #}

      output$img1 <- shiny::renderPlot(ggobj2) # CELL SCATTERPLOT
      output$img2 <- shiny::renderPlot(ggobj) # CELL HISTOGRAM
      # output$tab <- DT::renderDataTable(dtobj) # CELL DATATABLE

      output$img3 <- shiny::renderPlot(ggobj4) # GROUP SCATTERPLOT
      output$img4 <- shiny::renderPlot(ggobj3) # GROUP HISOTRAM

      gt <<- groupTable() # 16cell, 16cell

      # create table
      output$dynamicCell <- renderUI({
        Tabs <- unique(CellInfo)
        numTabs <- length(Tabs)
        CardColors <- briterhex(scales::hue_pal()(length(Tabs)))

        tagList(
          material_row(
            lapply(1:numTabs, function(i) {
              material_column(
                # solved material card
                shiny::tags$div(
                  class = paste("card", "z-depth-", "5"), # depth = null, color = null ; color will define in style
                  style = paste0("border : ", "solid 0.5em ", CardColors[i]),
                  shiny::tags$div(
                    class = "card-content",
                    shiny::tags$span(class = "card-title", Tabs[i]), # title
                    shiny::tags$div(class = "divider"), # divider = TRUE
                    DT::dataTableOutput(paste0("dtC", i), width = "100%", height = "500px"),
                    solvedButton(inputId = paste0("tsC", i), label = "Select", style = "position:absolute; top:1em; right:1em;")
                  )
                ),
                width = 4
              )
            })
          )
        )
      })


      # fill table
      g <- unique(CellInfo)

      # ordered by Count , not length;

      for (i in 1:length(g)) {
        t <- paste0(
          "output$dtC", i, " = DT::renderDataTable(datatable(CellPathwayDF[which(CellPathwayDF[,1]==g[", i, "]),-1]", # removed group column
          ", options = list(dom = 'ltp',scroller = TRUE, scrollX = TRUE, autoWidth = TRUE, lengthChange = FALSE, order = list(list(1,'desc'))), rownames = FALSE",
          ", selection = 'single'))"
        )
        eval(parse(text = t))
      }
    })

    observeEvent(input$freqbtn, {
      if (input$freqbtn == 0) {
        return(NULL)
      }

      Cells <- unique(CellInfo)

      res <- c()
      for (i in 1:length(Cells)) {
        thisCell <- Cells[i]

        thisCellData <- CellPathwayDF %>% dplyr::filter(Cell == thisCell)

        if (nrow(thisCellData) > 1) {
          thisCellData <- thisCellData %>% top_n(1, wt = Count)
          if (nrow(thisCellData) > 1) {
            thisCellData <- thisCellData %>% top_n(-1, wt = Length)
            if (nrow(thisCellData) > 1) {
              thisCellData <- thisCellData %>% top_n(-1, wt = Qvalue)
              if (nrow(thisCellData) > 1) {
                thisCellData <- thisCellData %>% top_n(1)
              }
            }
          }

          res <- c(res, paste0(thisCellData$Geneset, " @", thisCellData$Cell))
        }
      }

      ggobj2 <- emphasize(FALSE, res)

      output$img1 <- shiny::renderPlot(ggobj2)
    })

    observeEvent(input$sigbtn, {
      if (input$sigbtn == 0) {
        return(NULL)
      }

      Cells <- sort(unique(CellInfo))
      res <- c()
      for (i in 1:length(Cells)) {
        thisCell <- Cells[i]

        thisCellData <- CellPathwayDF %>% dplyr::filter(Cell == thisCell)

        if (nrow(thisCellData) > 1) {
          thisCellData <- thisCellData %>% top_n(-1, wt = Qvalue)

          if (nrow(thisCellData) > 1) {
            thisCellData <- thisCellData %>% top_n(-1, wt = Length)

            if (nrow(thisCellData) > 1) {
              thisCellData <- thisCellData %>% top_n(1, wt = Count)
              if (nrow(thisCellData) > 1) {
                thisCellData <- thisCellData %>% top_n(1)
              }
            }
          }

          res <- c(res, paste0(thisCellData$Geneset, " @", thisCellData$Cell))
        }
      }


      ggobj2 <- emphasize(FALSE, res)

      output$img1 <- shiny::renderPlot(ggobj2)
    })

    # draw gray color images
    observeEvent(input$graybtn, {
      if (input$graybtn == 0) {
        return(NULL)
      }

      ggobj2 <- ggplot(dfobj, aes(x = x, y = y)) +
        geom_point(colour = "gray")

      output$img1 <- shiny::renderPlot(ggobj2)
    })

    observeEvent(input$colorbtn, {
      if (input$colorbtn == 0) {
        return(NULL)
      }

      UniqueCol <- briterhex(scales::hue_pal()(length(unique(CellInfo))))
      names(UniqueCol) <- unique(CellInfo)

      colV <- unname(UniqueCol[dfobj$col])

      # scatter plot
      ggobj2 <- ggplot(dfobj, aes(x = x, y = y)) +
        geom_point(colour = colV)

      output$img1 <- shiny::renderPlot(ggobj2)
    })

    observeEvent(input$btn3, {
      if (input$btn3 == 0) {
        return(NULL)
      }

      shinyjs::hide("btn3")
    })

    output$dynamic <- renderUI({
      if (input$btn3 == 0) {
        return(NULL)
      }
      input$btn3

      Tabs <- sort(unique(gt[, 1]))
      numTabs <- length(Tabs)
      CardColors <- briterhex(scales::hue_pal()(numTabs))

      tagList(
        material_row(
          lapply(1:numTabs, function(i) {
            material_column(
              # solved material card
              shiny::tags$div(
                class = paste("card", "z-depth-", "5"), # depth = null, color = null ; color will define in style
                style = paste0("border : ", "solid 0.5em ", CardColors[i]),
                shiny::tags$div(
                  class = "card-content",
                  shiny::tags$span(class = "card-title", Tabs[i]), # title
                  shiny::tags$div(class = "divider"), # divider = TRUE
                  DT::dataTableOutput(paste0("dt", i), width = "100%"),
                  actionButton(inputId = paste0("ts", i), label = "Select")
                )
              ),
              width = 4
            )
          })
        )
      )
    })

    # call sort button in cell tab

    observeEvent(input$callSortFunction, {
      if (input$callSortFunction == 0) {
        return(NULL)
      }
      shinyjs::runjs(code = '$("#mysortableCell .rank-list-item").remove()')

      g <- sort(unique(gt[, 1]))

      options(useFancyQuotes = FALSE)

      for (i in 1:length(g)) {
        item <- paste0("$('#dtC", i, " .selected td')[0].innerText")

        shinyjs::runjs(
          code = paste0(
            "$('#tsC", i, "').attr('onClick',", '"', sortItem(paste0(item, " + ' @", g[i], "'"), "mysortableCell"), "; ",
            "$('#tsC", i, "').attr('disabled', true);", '");'
          )
        )
      }
      shinyjs::runjs(code = ' $("#btn5Cell, #btn6Cell, #btn7Cell").prop("disabled",false)')

      shinyjs::runjs(code = '$("#callSortFunction").hide();')
    })

    observeEvent(input$btn4, {
      if (input$btn4 == 0) {
        return(NULL)
      }
      g <- sort(unique(gt[, 1]))

      for (i in 1:length(g)) {
        t <- paste0(
          "output$dt", i, " = DT::renderDataTable(datatable(gt[which(gt[,1]==g[", i, "]),2:3]", # removed group column
          ", options = list(dom = ", "'ltp'", ",scroller = TRUE, scrollX = TRUE, autoWidth = TRUE, lengthChange = FALSE), rownames = FALSE",
          ", selection = ", "'single'", ", colnames =c(", "'Geneset',", "'P-value'", ")))"
        )
        eval(parse(text = t))
      }

      shinyjs::hide("btn3")
      shinyjs::hide("btn4")
    })

    # draw time plot in Cell tab
    observeEvent(input$btn5Cell, {
      if (input$btn5Cell == 0) {
        return(NULL)
      }

      ggobj2 <- emphasize(TRUE, input$rlistCell)

      # modify ggobj2
      output$img1 <- shiny::renderPlot(ggobj2)
    })

    observeEvent(input$btn6Cell, {
      if (input$btn6Cell == 0) {
        return(NULL)
      }
      ggobj2 <- emphasize(FALSE, input$rlistCell)

      output$img1 <- shiny::renderPlot(ggobj2)
    })

    # draw time plot in Group tab
    observeEvent(input$btn5, {
      if (input$btn5 == 0) {
        return(NULL)
      }
      ggobj2 <- emphasize(TRUE, input$rlist)

      output$img3 <- shiny::renderPlot(ggobj2)
    })


    # call sort button in group tab
    observeEvent(input$btn6, {
      shinyjs::runjs(code = '$("#mysortable .rank-list-item").remove()')
      g <- sort(unique(gt[, 1]))
      for (i in 1:length(g)) {
        item <- paste0("$('#dt", i, " .selected td')[0].innerText")
        shinyjs::runjs(
          code = paste0(
            "$('#ts", i, "').attr(", "'onClick'", ",",
            '"', sortItem(paste0(item, " + ' @", g[i], "'"), "mysortable"),
            ";$('#ts", i, "').attr('disabled', true);", '")'
          )
        )
      }
      shinyjs::hide("btn6")
    })

    # clear timelist in Group tab
    observeEvent(input$btn7, {
      shinyjs::runjs(code = '$("#mysortable .rank-list-item").remove(); $("#dynamic button").attr("disabled",false)')
    })

    # clear timelist in Cell tab
    observeEvent(input$btn7Cell, {
      if (input$btn7Cell == 0) {
        return(NULL)
      }
      shinyjs::runjs(code = '$("#mysortableCell .rank-list-item").remove(); $("#dynamicCell button").attr("disabled",false)')
    })
  }
  shiny::shinyApp(ui, server, options = list(launch.browser = TRUE))
}

CellEnrichUI <- function() {

  # actionButton with style Attributes, width = NULL
  solvedButton <- function(inputId, label, style = NULL, ...) {
    value <- restoreInput(id = inputId, default = NULL)
    tags$button(
      id = inputId, style = style,
      type = "button", class = "btn btn-default action-button",
      `data-val` = value, list(label),
      ...
    )
  }

  material_page(
    shinyjs::useShinyjs(),

    # dynamic datatable full width
    tags$head(tags$style(HTML(".display.dataTable.no-footer{width : 100% !important;}"))),
    use_waitress(color = "#697682", percent_color = "#333333"),
    title = "CellEnrich <a href = 'https://github.com/jhk0530/cellenrich' target = '_blank'> <i class='material-icons' style = 'font-size:1.3em;'>info</i> </a>",
    nav_bar_fixed = FALSE,
    nav_bar_color = "light-blue darken-1",
    font_color = "#ffffff",
    include_fonts = FALSE,
    include_nav_bar = TRUE,
    include_icons = FALSE,

    # options in navigator.
    #material_side_nav(
    material_row(
      material_column(
        material_card(
          title = "Options",
          divider = TRUE,
          material_radio_button(
            input_id = "dropdowninput",
            label = "select FC Option",
            choices = c("median", "mean", "zero", "GSVA"),
            selected = "median"
          ),
          material_radio_button(
            input_id = "dropdowninput2",
            label = "select Plot Option",
            choices = c("t-SNE", "U-MAP"),
            selected = "t-SNE"
          ),
          material_radio_button(
            input_id = "dropdowninput3",
            label = "select Gene-set",
            choices = c("Curated", "GeneOntology"),
            selected = "Curated"
          ),
          material_slider(
            input_id = "sliderinput1",
            label = "Minimum gene-set size",
            min_value = 10,
            max_value = 30,
            initial_value = 15,
            step_size = 5
          ),
          material_slider(
            input_id = "sliderinput2",
            label = "Maximum gene-set size",
            min_value = 250,
            max_value = 750,
            initial_value = 500,
            step_size = 5
          ),
          material_slider(
            input_id = "sliderinput3",
            label = "q-value threshold",
            min_value = 0,
            max_value = 0.25,
            initial_value = 0.1,
            step_size = 0.05
          ),
          solvedButton(inputId = "btn", label = "Start CellEnrich", style = 'margin-left:45%;'),
          depth = 3
        ),width = 6, offset = 3
      )
    )
    ,


    #),

    ## Tabs
    #material_tabs(
      #tabs = c(
        #"Cell" = "tab_cell"
      #),
      #color = "blue"
    #),

    #material_tab_content(
      #tab_id = "tab_cell",
      textOutput(outputId = "txtbox"),
      material_row(
        material_column(
          material_card(
            title = "tSNE/UMAP Plot",
            depth = 3,
            material_column(
              plotOutput("img1", height = "700px"),
              width = 6
            ),
            material_column(
              plotOutput("img2", height = "700px"), # cell distribution
              width = 6
            ),
            material_button("colorbtn", "toColor", icon = "color_lens", color = "blue lighten-1"),
            material_button("graybtn", "toGray", icon = "clear", color = "blue lighten-1"),
            material_button("freqbtn", "Frequent", icon = "grain", color = "blue lighten-1"),
            material_button("sigbtn", "Significant", icon = "grade", color = "blue lighten-1")
            # material_button('timeplot', 'timeplot', icon = 'timeline'),
          ),
          width = 12
        ),
        style = "margin : 1em"
      ),
      material_row(
        material_card(
          title = 'MarkerGenes',
          p('DE from each Cell specific'),
          DT::dataTableOutput('markerL1'),
          p('DE - Pathway from each Cell specific'),
          DT::dataTableOutput('markerL2')
        ),
        style = "margin : 1em"
      ),
      material_row(
        material_card(
          title = "",
          material_button("callSortFunction", "activate", icon = "play_arrow", color = "pink"),
          material_card(
            title = "Pathway Emphasize", divider = TRUE,
            tags$p("If list not recognized, please re-move their position"),
            rank_list(text = "Pathways", labels = "Please Activate First", input_id = "rlistCell", css_id = "mysortableCell"),
            material_button("btn5Cell", "Emphasize with Order", icon = "timeline"),
            material_button("btn6Cell", "Emphasize without Order", icon = "bubble_chart"),
            material_button("btn7Cell", "Clear List", icon = "clear_all"),
          ),
          uiOutput("dynamicCell"),
          depth = 3
        ),
        style = "margin : 1em"
      )
    #)
    #,
    #material_tab_content(
      #tab_id = "tab_group",
      #material_card(
        #plotOutput("img3", height = "700px"),
        #actionButton("btn3", "Create Table"),
        #actionButton("btn4", "Fill Table"),
        #actionButton("btn6", "call SortButtons"),
        #material_card(
          #title = "TimePlot", divider = TRUE,
          #tags$p("If list not recognized, please re-move their position"),
          #rank_list(text = "List for TimePlot", labels = "", input_id = "rlist", css_id = "mysortable"),
          #actionButton("btn5", "Generate Time Plot"),
          #actionButton("btn7", "Clear List"),
        #),
        #plotOutput("img4"),
        # DT::dataTableOutput('tab2'),
        #uiOutput("dynamic"),
        #depth = 3
      #)
    #),
  )
}
