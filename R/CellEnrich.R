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

CellEnrich <- function(CountData, CellInfo, ClustInfo = NULL, q0 = 0.1) {
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

  ui <- CellEnrichUI()

  options(useFancyQuotes = FALSE)

  server <- function(input, output, session) {

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

    sortItem <- function(label, tableName) {
      options(useFancyQuotes = FALSE)
      paste0(
        "$('#" ,tableName,  "')",
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
      uv <- sort(unique(v))
      res <- c()
      for (i in 1:length(uv)) {
        res[which(v == uv[i])] <- cols[i]
      }
      return(res)
    }

    pathwayPvalue = function(q0 = 0.1){
      res = c()
      gs = sapply(names(pres2),function(i){ which(names(genesets)== i )})
      Cells = sort(unique(CellInfo))

      for(i in 1:length(Cells)){
        thisCell = Cells[i]
        thisCellIdx = which(CellInfo == thisCell)
        total = length(CellInfo)
        k = length(thisCellIdx)
        thisCellPathways = table(unlist(pres[thisCellIdx]))
        pv = c()

        for(j in 1:length(thisCellPathways)){
          thisPathway = names(thisCellPathways)[j]

          q = unname(thisCellPathways[j]) # selected white ball
          m = pres2[names(genesets)[as.numeric(thisPathway)]] # total white ball
          pv[j] = round(1-phyper(q - 1, m, total - m, k), 4)
        }
        names(pv) = names(genesets)[as.numeric(names(thisCellPathways))]

        res = rbind(res, cbind(thisCell, names(pv), unname(pv)))
      }
      res = data.frame(res, stringsAsFactors = FALSE)
      colnames(res) = c('Cell', 'Geneset', 'Qvalue')
      res$Qvalue = round(p.adjust(as.numeric(res$Qvalue),'fdr'),4)
      res = res %>% dplyr::filter(Qvalue < q0)
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
          round(1-phyper(q - 1, m, tot - m, k), 4)
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

    getHyperPvalue <- function(genes, genesets) {
      pv <- sapply(1:length(genesets), function(i) {
        q <- length(intersect(genesets[[i]], genes)) # selected white ball
        m <- length(genesets[[i]]) # white ball
        n <- A - m # black ball
        k <- length(genes) # selected ball
        1 - phyper(q - 1, m, n, k)
      })
      names(pv) <- names(genesets)
      return(pv)
    }

    load("c2v7.RData")
    # load("c2v7idx.RData")
    A <- length(unique(unlist(genesets))) # Background genes
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
      w <- Waitress$new(selector = NULL, theme = "overlay")$start()

      v <- CountData

      if(input$dropdowninput != 'GSVA'){
        s <- findSigGenes(v, input$dropdowninput)
        names(s) <- CellInfo
      }

      # TODO BUILD GSVA CASE

      if(!is.null(ClustInfo)){
        cat('Clust Info found\n')
        s2 = findSigGenesGroup(v, ClustInfo, q0)
      }

      # for test
      res <- list()
      for (i in 1:length(s)) {
        w$inc(1 / length(s))
        pvh <- getHyperPvalue(s[[i]], genesets)
        qvh <- p.adjust(pvh, "fdr")
        res[[i]] <- unname(which(qvh < q0))
      }

      w$close()
      pres <<- res
      rm(res)

      # pres : which gene-sets are significant for each cells.
      # pres2 : for each gene-sets, how many cells are significant that gene-sets.

      CellPathwayDF = data.frame()
      Cells = unique(CellInfo)
      for(i in 1:length(Cells)){
        thisCell = Cells[i]
        tt = table(unlist(pres[which(thisCell==CellInfo)]))
        CellPathwayDF = rbind(CellPathwayDF, cbind(thisCell, names(tt), unname(tt)))
      }

      colnames(CellPathwayDF) = c('Cell', 'Geneset', 'Count')
      CellPathwayDF$Geneset = names(genesets)[as.numeric(CellPathwayDF$Geneset)]
      CellPathwayDF$Count = as.numeric(CellPathwayDF$Count)
      CellPathwayDF$Cell = as.character(CellPathwayDF$Cell)
      Length = sapply(CellPathwayDF$Geneset, function(i){length(genesets[[i]])})
      CellPathwayDF = cbind(CellPathwayDF, Length)

      # select genesets with count > 1
      CellPathwayDF = CellPathwayDF %>%
        dplyr::filter(Count > 1)

      # select genesets with count == max(Count)
      CellPathwayDF = CellPathwayDF %>%
        right_join( CellPathwayDF %>% group_by(Cell) %>% summarize(Count = max(Count)) )


      pres2 <- sort(table(unlist(pres)), decreasing = T)
      names(pres2) <- names(genesets)[as.numeric(names(pres2))]
      pres2 <<- pres2

      CellPathwayDF = CellPathwayDF %>%
        inner_join( pathwayPvalue() )

      CellPathwayDF <<- CellPathwayDF

      dtobj <<- buildDT(pres2)

      if(input$dropdowninput2 == 't-SNE'){
        tsneE <- Rtsne(t(v), check_duplicates = FALSE, perplexity = 15)
        dfobj <- data.frame(tsneE$Y, col = CellInfo)
        if(!is.null(ClustInfo)){
          dfobj2 <- data.frame(tsneE$Y, col = ClustInfo)
        }
      }

      if(input$dropdowninput2 == 'U-MAP'){
        umapE = uwot::umap(t(v), fast_sgd = TRUE)
        dfobj = data.frame(umapE, col = CellInfo)
        if(!is.null(ClustInfo)){
          dfobj2 <- data.frame(umapE, col = ClustInfo)
        }
      }

      ggobj <<- ggplot(data.frame(table(CellInfo)), aes(x = CellInfo, y = Freq, fill = CellInfo)) +
        geom_bar(stat = "identity") # cell histogram

      if(!is.null(ClustInfo)){
        ggobj3 <<- ggplot(data.frame(table(ClustInfo)), aes(x = ClustInfo, y = Freq, fill = ClustInfo)) +
          geom_bar(stat = "identity") # cell histogram
      }

      colnames(dfobj) <- colnames(dfobj2) <- c("x", "y", "col")
      dfobj <<- dfobj

      dfobj2$col = as.character(dfobj2$col)
      dfobj2 <<- dfobj2

      # scatter plot
      ggobj2 <<- ggplot(dfobj, aes(x = x, y = y, color = col)) +
        geom_point() +
        scale_color_manual(values = briterhex(scales::hue_pal()(length(unique(dfobj$col)))))

      if(!is.null(ClustInfo)){
        ggobj4 <<- ggplot(dfobj2, aes(x = x, y = y, color = col) ) +
          geom_point() +
          scale_color_manual(values = briterhex(scales::hue_pal()(length(unique(dfobj2$col)))))
      }

      output$img1 <- shiny::renderPlot(ggobj2) # CELL SCATTERPLOT
      output$img2 <- shiny::renderPlot(ggobj) # CELL HISTOGRAM
      # output$tab <- DT::renderDataTable(dtobj) # CELL DATATABLE

      output$img3 <- shiny::renderPlot(ggobj4) # GROUP SCATTERPLOT
      output$img4 <- shiny::renderPlot(ggobj3) # GROUP HISOTRAM

      gt <<- groupTable()



      output$dynamicCell <- renderUI({

        Tabs <- sort(unique(CellInfo))
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
                    DT::dataTableOutput(paste0("dtC", i), width = "100%", height = '500px'),
                    actionButton(inputId = paste0("tsC", i), label = "Select")
                  )
                ),
                width = 4
              )
            })
          )
        )
      })

      g <- unique(sort(CellInfo))

      for (i in 1:length(g)) {
        t <- paste0(
          "output$dtC", i, " = DT::renderDataTable(datatable(CellPathwayDF[which(CellPathwayDF[,1]==g[", i, "]),-1]", # removed group column
          ", options = list(dom = 'ltp',scroller = TRUE, scrollX = TRUE, autoWidth = TRUE, lengthChange = FALSE, order = list(list(2,'asc'))), rownames = FALSE",
          ", selection = 'single'))"
        )
        eval(parse(text = t))
      }




    })


    # draw gray color images
    observeEvent(input$graybtn, {
      if (input$graybtn == 0) {
        return(NULL)
      }

      ggobj2 <- ggplot(dfobj, aes(x= x, y = y)) +
        geom_point(colour = 'gray')

      output$img1 <- shiny::renderPlot(ggobj2)
    })

    observeEvent(input$colorbtn, {
      if(input$colorbtn ==0){
        return(NULL)
      }

      ggobj2 <<- ggplot(dfobj, aes(x = x, y = y, color = col)) +
        geom_point() +
        scale_color_manual(values = briterhex(scales::hue_pal()(length(unique(dfobj$col)))))

      output$img1 <- shiny::renderPlot(ggobj2)
    })

    observeEvent(input$btn3, {
      if (input$btn3 == 0) {
        return(NULL)
      }

      shinyjs::hide('btn3')
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
      if(input$callSortFunction ==0 ){return(NULL)}
      shinyjs::runjs(code = '$("#mysortableCell .rank-list-item").remove()')

      g <- sort(unique(gt[, 1]))

      options(useFancyQuotes = FALSE)

      for (i in 1:length(g)) {

        item <- paste0("$('#dtC", i, " .selected td')[0].innerText")

        shinyjs::runjs(
          code = paste0(
            "$('#tsC", i, "').attr('onClick',",'"', sortItem(paste0(item, " + ' @", g[i], "'"), 'mysortableCell'),"; ",
            "$('#tsC", i, "').attr('disabled', true);", '")'
          )
        )
      }

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


    # draw time plot in Group tab
    observeEvent(input$btn5, {
      if (input$btn5 == 0) {
        return(NULL)
      }

      rlobj = data.frame(stringsAsFactors = FALSE)

      items = input$rlist

      for(i in 1:length(items)){
        kk = strsplit(items[[i]], ' @')[[1]]
        name = kk[1]
        location = kk[2]
        rlobj = rbind(rlobj, cbind(name, location))
      }
      colnames(rlobj) = c('name', 'location')
      rlobj$name = as.character(rlobj$name)
      rlobj$location = as.character(rlobj$location)

      getCellValues = function(rlobj){
        ret = list()
        for(i in 1:nrow(rlobj)){

          thisGeneset = which(names(genesets)==rlobj[i,1])
          thisGroup = rlobj[i,2]

          thisCellsIdx = which(dfobj$col==thisGroup)
          thisCells = dfobj[thisCellsIdx,]

          rn = as.numeric(rownames(thisCells))
          res = c()
          for(j in 1:nrow(thisCells)){
            if( thisGeneset %in% pres[[ rn[j] ]] ){
              res = c(res,thisCellsIdx[j])
            }
          }

          ret[[i]] = res
        }
        names(ret) = rlobj$location
        return(ret)
      }
      cellValues <- getCellValues(rlobj)


      dfobj_new <- data.frame(dfobj)

      for(i in 1:length(cellValues)){
        x = mean(as.numeric(dfobj_new$x[cellValues[[i]]]))
        y = mean(as.numeric(dfobj_new$y[cellValues[[i]]]))
        dfobj_new = rbind(dfobj_new, c( x, y, rlobj[i,2]))
      }
      colnames(dfobj_new) = c('x','y','col')

      newIdx = (nrow(dfobj)+1):nrow(dfobj_new)
      cellValues = c(unname(unlist(cellValues)), newIdx)
      dfobj_new$x = round(as.numeric(dfobj_new$x), 4)
      dfobj_new$y = round(as.numeric(dfobj_new$y), 4)
      colV <- changeCol(dfobj_new$col)
      colV[-cellValues] <- "#95A5A6" # gray color

      # define ggobj2 with curve element
      graphString = 'ggobj2 <- ggplot(dfobj_new, aes(x = x, y = y)) + geom_point(colour = colV)'

      for(i in 1:(length(newIdx)-1)){
        newCurve = paste(' + geom_curve( aes(x = ', 'dfobj_new$x[newIdx[',i,
                         ']], y = dfobj_new$y[newIdx[',i,']], xend = dfobj_new$x[newIdx[',i+1,
                         ']], yend = dfobj_new$y[newIdx[',i+1,']]), size = 0.5, linetype = "longdash",',
                         'curvature = 0.1, colour = colV[newIdx[',i,']], arrow = arrow(length = unit(0.1,"inches")))')
        graphString = paste(graphString, newCurve, sep = '')
      }
      eval(parse(text = graphString))

      # modify ggobj2
      output$img3 <- shiny::renderPlot(ggobj2)
    })

    # TODO BUILD Btn5 Cell;
    # btn5 cell : generate time plot

    # TODO BUILD btn6 cell:
    # btn6 cell : generate color plot



    # call sort button in group tab
    observeEvent(input$btn6, {
      shinyjs::runjs(code = '$("#mysortable .rank-list-item").remove()')
      g <- sort(unique(gt[, 1]))
      for (i in 1:length(g)) {
        item <- paste0("$('#dt", i, " .selected td')[0].innerText")
        shinyjs::runjs(
          code = paste0(
            "$('#ts", i, "').attr(", "'onClick'", ",",
            '"', sortItem(paste0(item, " + ' @", g[i], "'"), 'mysortable'),
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
      shinyjs::runjs(code = '$("#mysortableCell .rank-list-item").remove(); $("#dynamicCell button").attr("disabled",false)')
    })


  }
  shiny::shinyApp(ui, server, options = list(launch.browser = TRUE))
}

CellEnrichUI <- function() {
  material_page(
    shinyjs::useShinyjs(),
    # dynamic datatable full width
    tags$head(tags$style(HTML(".display.dataTable.no-footer{width : 100% !important;}"))),
    use_waitress(color = "#697682", percent_color = "#333333"),
    title = "CellEnrich",
    nav_bar_fixed = FALSE,
    nav_bar_color = "light-blue darken-1",
    font_color = "#ffffff",
    include_fonts = FALSE,
    include_nav_bar = TRUE,
    include_icons = FALSE,

    ## Tabs
    material_tabs(
      tabs = c(
        "Start" = "tab_start",
        "Cell" = "tab_cell",
        "Group" = "tab_group"
      ),
      color = "blue"
    ),

    # Define tab content
    material_tab_content(
      tab_id = "tab_start",
      material_row(
        material_column(
          material_card(
            title = "Options",
            material_dropdown(
              input_id = "dropdowninput",
              label = "select FC Option",
              choices = c("median", "mean", "zero", "GSVA"),
              selected = "median"
            ),
            material_dropdown(
              input_id = "dropdowninput2",
              label = "select Plot Option",
              choices = c("t-SNE", "U-MAP"),
              selected = "median"
            ),
            actionButton("btn", "Start CellEnrich")
          ),
          width = 12
        ),
        style = "margin : 1em"
      ),
      textOutput("txt1")
    ),

    material_tab_content(
      tab_id = "tab_cell",
      textOutput(outputId = "txtbox"),
      material_row(
        material_column(
          material_card(
            depth = 3,
            plotOutput("img1", height = "700px"),
            material_card(
              title='',
              material_button('colorbtn', 'toColor', icon = 'brush', color = 'blue lighten-1'),
              material_button('graybtn', 'toGray', color = 'grey darken-1'), # to gray color
              material_button('timeplot', 'timeplot', icon = 'timeline'),
              depth = 2
            )
          ),
          width = 12
        ),
        style = "margin : 1em"
      ),
      material_row(

        material_card(
          # plotOutput("img2"), # cell distribution
          title = '',
          material_button('callSortFunction', 'activate', icon = 'play_arrow', color = 'pink'),
          material_card(
            title = 'Pathway Emphasize', divider = TRUE,
            tags$p('If list not recognized, please re-move their position'),
            rank_list(text = "Pathways", labels = "Please Activate First", input_id = "rlistCell", css_id = "mysortableCell"),
            actionButton("btn5Cell", "Emphasize with Order"),
            actionButton("btn6Cell", "Emphasize without Order"),
            actionButton("btn7Cell", "Clear List"),
          ),
          uiOutput("dynamicCell"),
          depth = 3
        ),
        style = "margin : 1em"
      )
    ),
    material_tab_content(
      tab_id = "tab_group",
      material_card(
        plotOutput("img3", height = '700px'),
        actionButton("btn3", "Create Table"),
        actionButton("btn4", "Fill Table"),
        actionButton("btn6", "call SortButtons"),
        material_card(
          title = 'TimePlot', divider = TRUE,
          tags$p('If list not recognized, please re-move their position'),
          rank_list(text = "List for TimePlot", labels = "", input_id = "rlist", css_id = "mysortable"),
          actionButton("btn5", "Generate Time Plot"),
          actionButton("btn7", "Clear List"),
        ),
        plotOutput("img4"),
        #DT::dataTableOutput('tab2'),
        uiOutput("dynamic"),
        depth = 3
      )
    ),
  )
}



