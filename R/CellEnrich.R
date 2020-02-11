#' @import SingleCellExperiment
#' @import Rtsne
#' @import shinyCyJS
#' @import shiny
#' @import DT
#' @import ggplot2
#' @import Seurat
#' @import htmltools
#' @import magrittr
#' @import shinymaterial

CellEnrich = function(yan){
  v = extractArray(yan)
  n = cellNames(yan)
  s = findSigGenes(v)
  names(s) = n

  pres = findPathway(s) ## TAKES LONG TIME

  pres2 = sort(table(unlist(pres)), decreasing = T)
  names(pres2) = names(genesets)[as.numeric(names(pres2))]

  dtobj = buildDT(pres2)

  tsneE = Rtsne(t(v), check_duplicates = FALSE, perplexity = 15)

  ggobj =
    ggplot(nasDF(n), aes(x = name, fill = name)) + geom_histogram(stat = 'count')

  dfobj = data.frame(tsneE$Y, col = n)
  colnames(dfobj) = c('x','y', 'col')

  ggobj2 =
    ggplot(dfobj, aes(x = x, y = y, color = col)) + geom_point()

  ui = material_page(
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
            material_file_input(
              input_id = 'fileinput',
              label = 'upload RDS File'
            ),
            material_button('btn','Submit')
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
            dataTableOutput('tab', height = '600px')
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

    observeEvent(input$btn,{
        in_file = input$fileinput
        if(is.null(in_file)){return(NULL)}
        contents = readLines(in_file$datapath)
        output$txt1 = renderText(contents)
    })

    output$img1 = renderPlot(ggobj2)
    output$img2 = renderPlot(ggobj)
    output$tab = renderDataTable(dtobj)

  }

  shinyApp(ui,server)
}


# load sample Data
# yan = readRDS('yan.rds')
# load sample geneset Data
# load("c2v7.RData")



