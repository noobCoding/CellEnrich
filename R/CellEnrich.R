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
#'
#' @export

CellEnrich = function(scData){
  require(shinymaterial)
  require(shiny)
  require(waiter)
  require(Rtsne)
  require(ggplot2)

  ui = material_page(
    use_waitress(color = '#697682', percent_color = '#333333'),
    title = "CellEnrich",
    nav_bar_fixed = FALSE,
    nav_bar_color = "light-blue darken-1",
    font_color = '#ffffff',
    include_fonts = FALSE,
    include_nav_bar = TRUE,
    include_icons = FALSE,

    ## Tabs
    material_tabs(
      tabs = c(
        "Start" = "tab_start",
        "Cell" = "tab_cell",
        "Cluster" = "tab_cluster"
      ),
      color = 'blue'
    ),

    ## navigator

    # Define tab content
    material_tab_content(
      tab_id = 'tab_start',
      material_row(
        material_column(
          material_card(title = 'card4',
            fileInput(
              inputId = 'fileinput',
              label = 'upload RDS File',
              accept = c('.RDS','.rds')
            ),
            div(
              actionButton('btn','Submit', style = 'color : #e53935')
              ,style = 'margin-left : 45%'
            )
          ),
          width = 6
        ),
        material_column(
          material_card(title = 'card5',
            material_dropdown(
              input_id = 'dropdowninput',
              label = 'select Method',
              choices = c('mean', 'median', '0'),
              selected = 'mean'
            )
          ),
          width = 6
        ),
        style = 'margin : 1em'
      ),
      textOutput('txt1')
    ),

    material_tab_content(
      tab_id = "tab_cell",
      tags$h1("cell Tab Content"),
      material_button('btn2','btn2'),
      #material_button('btn3','btn3'),
      textOutput(outputId = 'txtbox'),
      material_row(
        material_column(
          material_card(
            title = 'card 1',
            depth = 3,
            plotOutput('img1', height = '600px')
          ),
        width = 6) ,
        material_column(
          material_card(
            title = 'card 2',
            depth = 3,
            DT::dataTableOutput('tab', height = '600px')
            ),
          width =6),
        style = 'margin : 1em'
      ),
      material_row(
        material_card(
          plotOutput('img2'), title = 'card 3', depth = 3
        ), style = 'margin : 1em'
      )
    ),
    material_tab_content(
      tab_id = "tab_cluster",
      tags$h1("clutser Tab Content")
    ),
  )

  server = function(input, output, session){

    load("c2v7.RData")

    myf = function(genesetnames, pres){
      idx = which(names(genesets)== genesetnames)
      res = c()
      for(i in 1:length(pres)){
        if(idx %in% pres[[i]]){
          res = c(res, i)
        }
      }
      return(res)
    }

    ggobj2 = ggobj = dtobj = dfobj = pres = ''
    observeEvent(input$btn, {

      w = Waitress$new(selector = NULL, theme = 'overlay')$start()
      v = extractArray(scData)
      n = cellNames(scData)
      s = findSigGenes(v)
      names(s) = n

      res = list()
      for(i in 1:length(s)){
        w$inc(1/length(s))
        pvh = getHyperPvalue(s[[i]] , genesets)
        qvh = p.adjust(pvh, 'fdr')
        res[[i]] = unname(which(qvh<0.1))
      }

      w$close()
      pres <<- res; rm(res)

      pres2 = sort(table(unlist(pres)), decreasing = T)
      names(pres2) = names(genesets)[as.numeric(names(pres2))]

      dtobj <<- buildDT(pres2)

      tsneE = Rtsne(t(v), check_duplicates = FALSE, perplexity = 15)

      ggobj <<- ggplot(nasDF(n), aes(x = name, fill = name)) +
        geom_histogram(stat = 'count', binwidth = 0.2)

      dfobj = data.frame(tsneE$Y, col = n)
      colnames(dfobj) = c('x','y', 'col')
      dfobj <<- dfobj

      ggobj2 <<- ggplot(dfobj, aes(x = x, y = y, color = col)) +
        geom_point()

      output$img1 = shiny::renderPlot(ggobj2)
      output$img2 = shiny::renderPlot(ggobj)
      output$tab = DT::renderDataTable(dtobj)

    })

    observeEvent(input$btn2,{
      if(input$btn2 == 0){return(NULL)}
      cellValues = input$tab_cell_clicked

      cellValues = cellValues$value

      dfobj_new = data.frame(dfobj, size = 1)
      colnames(dfobj_new) = c('x','y','col','size')

      dfobj_new$size[myf(cellValues, pres)]  = 2
      ggobj2 = ggplot(dfobj_new, aes(x = x, y = y, color = col, size = size)) +
        geom_point()
      output$img1 = shiny::renderPlot( ggobj2 )
    })
  }

  shiny::shinyApp(ui,server, options = list(launch.browser = TRUE))
}

findPathway = function(s, w, genesets){
  res = list()
  for(i in 1:length(s)){
    w$inc(1/length(s))
    pvh = getHyperPvalue(s[[i]] , genesets)
    qvh = p.adjust(pvh, 'fdr')
    res[[i]] = unname(which(qvh<0.1))
  }
  w$hide()
  return(res)
}




# load sample Data
# yan = readRDS('yan.rds')
# load sample geneset Data
# load("c2v7.RData")

