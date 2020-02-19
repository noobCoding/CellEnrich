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
#' @import slickR
#'
#' @export

CellEnrich = function(scData){
  require(shinymaterial)
  require(shiny)
  require(waiter)
  require(Rtsne)
  require(ggplot2)
  require(DT)
  require(slickR)

  ui = CellEnrichUI()

  server = function(input, output, session){

    buildDT = function(pres2){
      DT::datatable(
        data.frame(
          Geneset = names(pres2),
          Count = as.numeric(pres2)
        ),
        options = list(
          dom = 'ltp'
        ),
        rownames = FALSE,
        selection = 'single'
      )
    }

    mapColor = function(v){
      as.factor(v)
    }




    # load sample Data
    # yan = readRDS('yan.rds')
    # load sample geneset Data
    # load("c2v7.RData")

    groupTable = function(){

      # for pres2
      genesetIdx = sapply(names(pres2), function(i){which(i==names(genesets))}, USE.NAMES = FALSE)
      pres2Idx = pres2
      names(pres2Idx) = genesetIdx

      groups = as.character(unique(dfobj$col))
      res = data.frame(stringsAsFactors = FALSE)

      for(i in 1:length(groups)){
        pathways = rep(0,length(genesets))
        names(pathways) = 1:length(genesets)
        tt = table(unlist(pres[which(dfobj$col==groups[i])]))
        pathways[as.numeric(names(tt))] = unname(tt)

        pathways = pathways[which(pathways!=0)] # remove zero genesets

        # what genesets are enriched per each group.
        tot = sum(pres2Idx)
        gt = sort(sapply(1:length(pathways), function(i){
          q = pathways[i] # selected white ball
          m = unname(pres2Idx[names(pathways[i])]) # total white ball
          n = tot-m # total black ball
          k = sum(pathways)
          round(phyper(q-1,m,n,k),4)
        }))
        gt = gt[which(gt<0.25)] # pvalue 0.25
        res = rbind(res, cbind(groups[i], names(gt), unname(gt)))
      }
      colnames(res) = c('groups',"genesetidx", 'pvalue')

      res$groups = as.character(res$groups)
      res$genesetidx = as.numeric(as.character(res$genesetidx))
      res$genesetidx = sapply(res$genesetidx, function(i){names(genesets)[i]})
      res$pvalue = as.numeric(as.character(res$pvalue))

      return(res)

    }


    getHyperPvalue <- function(genes, genesets) {
      pv = sapply(1:length(genesets),function(i){
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
    #load("c2v7idx.RData")
    A <- length(unique(unlist(genesets)))
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

    ggobj2 = ggobj = dtobj = dfobj = pres = gt = pres2 = ''
    observeEvent(input$btn, {

      w = Waitress$new(selector = NULL, theme = 'overlay')$start()
      v = extractArray(scData)
      n = cellNames(scData)
      s = findSigGenes(v)
      names(s) = n
      # for test
      res = list()
      for(i in 1:length(s)){
        w$inc(1/length(s))
        pvh = getHyperPvalue(s[[i]], genesets)
        qvh = p.adjust(pvh, 'fdr')
        res[[i]] = unname(which(qvh<0.1))
      }

      w$close()
      pres <<- res; rm(res)

      pres2 = sort(table(unlist(pres)), decreasing = T)
      names(pres2) = names(genesets)[as.numeric(names(pres2))]
      pres2 <<- pres2

      dtobj <<- buildDT(pres2)

      tsneE = Rtsne(t(v), check_duplicates = FALSE, perplexity = 15)

      ggobj <<- ggplot(data.frame(table(n)), aes(x = n,y = Freq, fill = n)) +
        geom_bar(stat = 'identity')

      dfobj = data.frame(tsneE$Y, col = n)
      colnames(dfobj) = c('x','y', 'col')
      dfobj <<- dfobj

      ggobj2 <<- ggplot(dfobj, aes(x = x, y = y, color = col)) +
        geom_point()

      output$img1 = shiny::renderPlot(ggobj2)
      output$img2 = shiny::renderPlot(ggobj)
      output$tab = DT::renderDataTable(dtobj)
      gt <<- groupTable()
      dtobj2 = DT::datatable(gt, options = list(dom = 'ltp'), rownames = FALSE, selection = 'single')

      output$tab2 = DT::renderDataTable(dtobj2)



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

    output$dynamic = renderUI({
      if(input$btn3==0){return(NULL)}
      input$btn3

      Tabs = unique(gt[,1])
      numTabs = length(Tabs)

      tagList(
        material_row(
          lapply(1:numTabs, function(i){
            material_column(
              material_card(
                title = Tabs[i],
                DT::dataTableOutput(paste0('dt',i)),
                #material_button(paste0('btn4card',i),label = paste0('min',i)),
                divider = TRUE
              ),
              width = 4
            )
          })
        )
      )
    })

    observeEvent(input$btn4, {
      if(input$btn4==0){return(NULL)}
      #shinyjs::hide('btn4')
      #shinyjs::hide('btn3')
      g = unique(gt[,1])
      for(i in 1:length(g)){
        t = paste0(
          'output$dt',
          i,
          ' = DT::renderDataTable(datatable(gt[which(gt[,1]==g[',i,']),]',
          ',options = list(dom = ',"'ltp'",', autoWidth = TRUE), rownames = FALSE',
          ', selection = ',"'single'",
          '))')

        eval( parse( text = t ) )
      }
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

CellEnrichUI = function(){
  material_page(
    shinyjs::useShinyjs(),
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
        "Group" = "tab_group"
      ),
      color = 'blue'
    ),

    ## navigator

    # Define tab content
    material_tab_content(
      tab_id = 'tab_start',
      material_row(
        material_column(
          material_card(
            title = 'card4',
            div(
              actionButton(
                'btn',
                'Submit',
                style = 'color : #e53935'
              ),
              style = 'margin-left : 45%'
            )
          ),
          width = 6
        ),
        material_column(
          material_card(
            title = 'card5',
            material_dropdown(
              input_id = 'dropdowninput',
              label = 'select Method',
              choices = c('median', 'GSVA'),
              selected = 'median'
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
      textOutput(outputId = 'txtbox'),
      material_row(
        material_column(
          material_card(
            title = 'card 1',
            depth = 3,
            plotOutput('img1', height = '700px')
          ),
          width = 6) ,
        material_column(
          material_card(
            title = 'card 2',
            depth = 3,
            DT::dataTableOutput('tab', height = '600px'),
            material_button('btn2','btn2'),
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
      tab_id = "tab_group",
      material_card(
        DT::dataTableOutput('tab2'),
        actionButton('btn3', 'btn3'),
        actionButton('btn4', 'btn4'),
        uiOutput('dynamic'),
        depth = 3
      )
    ),
  )
}

