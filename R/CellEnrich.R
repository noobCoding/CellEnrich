#' @title CellEnrich
#'
#'
#' @importFrom DT dataTableOutput
#'
#' @import SingleCellExperiment
#' @import Rtsne
#' @import shinyCyJS
#' @rawNamespace import(shiny, except = dataTableOutput)
#' @import ggplot2
#' @import Seurat
#' @import htmltools
#' @import magrittr
#' @import shinymaterial
#' @import waiter
#' @import shinyjs
#' @import scales
#' @import sortable
#'
#' @export

CellEnrich <- function(scData) {
  require(shinymaterial)
  require(shiny)
  require(waiter)
  require(Rtsne)
  require(ggplot2)
  require(DT)
  require(scales)
  require(sortable)

  ui <- CellEnrichUI()

  options(useFancyQuotes = FALSE)

  server <- function(input, output, session) {
    buildDT <- function(pres2) {
      DT::datatable(
        data.frame(
          Geneset = names(pres2),
          Count = as.numeric(pres2)
        ),
        options = list(
          dom = "ltp"
        ),
        rownames = FALSE,
        selection = "single"
      )
    }

    mapColor <- function(v) {
      as.factor(v)
    }

    changeCol = function(v){
      cols = briterhex(scales::hue_pal()(length(unique(v))))
      uv = unique(v)
      res = c()
      for(i in 1:length(uv)){
        res[which(v==uv[i])] = cols[i]
      }
      return(res)
    }

    # load sample Data
    # yan = readRDS('yan.rds')
    # load sample geneset Data
    # load("c2v7.RData")

    groupTable <- function() {

      # for pres2
      genesetIdx <- sapply(names(pres2), function(i) {
        which(i == names(genesets))
      }, USE.NAMES = FALSE)
      pres2Idx <- pres2
      names(pres2Idx) <- genesetIdx

      groups <- as.character(unique(dfobj$col))
      res <- data.frame(stringsAsFactors = FALSE)

      for (i in 1:length(groups)) {
        pathways <- rep(0, length(genesets))
        names(pathways) <- 1:length(genesets)
        tt <- table(unlist(pres[which(dfobj$col == groups[i])]))
        pathways[as.numeric(names(tt))] <- unname(tt)

        pathways <- pathways[which(pathways != 0)] # remove zero genesets

        # what genesets are enriched per each group.
        tot <- sum(pres2Idx)
        gt <- sort(sapply(1:length(pathways), function(i) {
          q <- pathways[i] # selected white ball
          m <- unname(pres2Idx[names(pathways[i])]) # total white ball
          n <- tot - m # total black ball
          k <- sum(pathways)
          round(phyper(q - 1, m, n, k), 4)
        }))
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
        q <- length(intersect(genesets[[i]], genes))
        m <- length(genesets[[i]])
        n <- A - m
        k <- length(genes)
        1 - phyper(q - 1, m, n, k)
      })
      names(pv) <- names(genesets)
      return(pv)
    }

    load("c2v7.RData")
    # load("c2v7idx.RData")
    A <- length(unique(unlist(genesets)))
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

    ggobj2 <- ggobj <- dtobj <- dfobj <- pres <- gt <- pres2 <- ""
    observeEvent(input$btn, {
      w <- Waitress$new(selector = NULL, theme = "overlay")$start()
      v <- extractArray(scData)
      n <- cellNames(scData)
      s <- findSigGenes(v)
      names(s) <- n
      # for test
      res <- list()
      for (i in 1:length(s)) {
        w$inc(1 / length(s))
        pvh <- getHyperPvalue(s[[i]], genesets)
        qvh <- p.adjust(pvh, "fdr")
        res[[i]] <- unname(which(qvh < 0.1))
      }

      w$close()
      pres <<- res
      rm(res)

      pres2 <- sort(table(unlist(pres)), decreasing = T)
      names(pres2) <- names(genesets)[as.numeric(names(pres2))]
      pres2 <<- pres2

      dtobj <<- buildDT(pres2)

      tsneE <- Rtsne(t(v), check_duplicates = FALSE, perplexity = 15)

      ggobj <<- ggplot(data.frame(table(n)), aes(x = n, y = Freq, fill = n)) +
        geom_bar(stat = "identity") # cell histogram

      dfobj <- data.frame(tsneE$Y, col = n)
      colnames(dfobj) <- c("x", "y", "col")
      dfobj <<- dfobj
      # scatter plot
      ggobj2 <<- ggplot(dfobj, aes(x = x, y = y, color = col)) +
        geom_point() +
        scale_color_manual( values = briterhex(scales::hue_pal()(length(unique(dfobj$col)))))

      output$img1 <- shiny::renderPlot(ggobj2)
      output$img2 <- shiny::renderPlot(ggobj)
      output$img3 <- shiny::renderPlot(ggobj2)
      output$tab <- DT::renderDataTable(dtobj)
      gt <<- groupTable()
    })

    observeEvent(input$btn2, {
      if (input$btn2 == 0) {
        return(NULL)
      }
      cellValues <- input$tab_cell_clicked

      selectedRow = input$tab_rows_selected # check none selected

      # if not selected : return;
      if(is.null(selectedRow)){
        output$img1 <- shiny::renderPlot(ggobj2)
        return(NULL)
      }
      cellValues <- cellValues$value

      dfobj_new <- data.frame(dfobj[order(dfobj$col),])

      colV = changeCol(dfobj_new$col)
      colV[-myf(cellValues, pres)] = '#95a5a6'

      ggobj2 <- ggplot(dfobj_new, aes(x = x, y = y)) +
        geom_point(colour = colV)

      output$img1 <- shiny::renderPlot(ggobj2)
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
                style = paste0("border : ", 'solid 0.5em ',CardColors[i]),
                shiny::tags$div(
                  class = "card-content",
                  shiny::tags$span(class = "card-title", Tabs[i]), # title
                  shiny::tags$div(class = "divider"), # divider = TRUE
                  DT::dataTableOutput(paste0("dt", i), width = "100%"),
                  actionButton(inputId = paste0('ts',i),label = 'ToSort')
                )
              ),
              width = 4
            )
          })
        )
      )
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
          ", selection = ", "'single'", ", colnames =c(", "'Geneset'",",","'P-value'",")))")
        eval(parse(text = t))
      }

      shinyjs::hide('btn3')
      shinyjs::hide('btn4')
    })

    observeEvent(input$btn5, {
      if (input$btn5 == 0) {
        return(NULL)
      }
      # modify ggobj2
      output$img3 <- shiny::renderPlot(ggobj2)
    })

    observeEvent(input$btn6, {
      g <- sort(unique(gt[, 1]))
      for(i in 1:length(g)){
        shinyjs::runjs(code = paste0("$('#ts",i,"').attr(","'onClick'", ',', sortItem(i), ')'))
      }
      shinyjs::hide('btn6')
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
            div(
              actionButton("btn","Start CellEnrich"),
              style = "margin-left : 45%"
            )
          ),
          width = 6
        ),
        material_column(
          material_card(
            title = "card5",
            material_dropdown(
              input_id = "dropdowninput",
              label = "select Method",
              choices = c("median", "GSVA"),
              selected = "median"
            )
          ),
          width = 6
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
            title = "card 1",
            depth = 3,
            plotOutput("img1", height = "700px")
          ),
          width = 6
        ),
        material_column(
          material_card(
            title = "card 2",
            depth = 3,
            DT::dataTableOutput("tab", height = "600px"),
            material_button("btn2", "btn2"),
          ),
          width = 6
        ),
        style = "margin : 1em"
      ),
      material_row(
        material_card(
          plotOutput("img2"),
          title = "card 3", depth = 3
        ),
        style = "margin : 1em"
      )
    ),
    material_tab_content(
      tab_id = "tab_group",
      material_card(
        plotOutput("img3"),
        actionButton("btn3", "Create Table"),
        actionButton("btn4", "Fill Table"),
        actionButton("btn5", "Colorize"),
        actionButton('btn6', 'call SortButtons'),
        rank_list(text = 'text', labels = '', input_id = 'rlist', css_id = 'mysortable'),
        uiOutput("dynamic"),
        depth = 3
      )
    ),
  )
}

briterhex = function(colors){
  res = c()
  for(i in 1:length(colors)){
    v = as.vector(col2rgb(colors[i])) * 1.3
    v = sapply(v, function(i){min(i,255)})
    res[i] = rgb(v[1],v[2],v[3],max = 255)
  }
  return(res)
}

# actionButton with Onclick Attributes
myButton = function(inputId, label, width = NULL, onClick = NULL, ...){
  value <- restoreInput(id = inputId, default = NULL)
  tags$button(id = inputId, style = if (!is.null(width))
    paste0("width: ", validateCssUnit(width), ";"),
    type = "button", class = "btn btn-default action-button",
    onClick = onClick,
    `data-val` = value, list(label),
    ...)
}

sortItem = function(label){
  paste0(
    "$('#mysortable')", '.append(', "`<div class=", "'rank-list-item'", ' draggable=',"'true'",
    ' style = ', "'transform: translateZ(0px);'", ">", label," </div>`)"
  )
}
