#' @title CellEnrich
#'
#'
#' @importFrom DT dataTableOutput
#' @importFrom Matrix t
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
#' @import highcharter
#'
#' @export

CellEnrich <- function(CountData, GroupInfo, genesets = NULL) {

  require(dplyr)

  options(useFancyQuotes = FALSE)

  server <- function(input, output, session) {

    ### CODES


    # variable initialize

    dtobj <- dfobj <- pres <- pres2 <- ""
    CellPathwayDF <- ""
    gt <- Cells <- A <- ""
    CellScatter <- ""
    CellHistogram <- ""

    observeEvent(input$StartCellEnrich, {

      pt <- proc.time()

      if(is.null(genesets) ){
        if(input$genesetOption == 'User-Geneset'){
          shiny::showNotification('Geneset not given ...', type = 'error', duration = 10)
          return(NULL)
        }
      }

      # ------ Hide Start Button

      shinyjs::hide("StartCellEnrich")

      # ------ Load Genesets

      if(is.null(genesets)){
        if (input$genesetOption == "Human-Curated") load("c2v7.RData")
        if (input$genesetOption == "Human-GO") load("humanGO.RData")
        if (input$genesetOption == "Human-GO-BP") load("humanGOBP.RData")
        if (input$genesetOption == "Human-GO-CC") load("humanGOCC.RData")
        if (input$genesetOption == "Human-GO-MF") load("humanGOMF.RData")
        if (input$genesetOption == "Human-KEGG") load("keggv7.RData")
        if (input$genesetOption == "Mouse-KEGG") load("mouseKEGG.RData")
        if (input$genesetOption == "Mouse-GO") load("mouseGO.RData")
        if (input$genesetOption == "Mouse-GO-BP") load("mouseGOBP.RData")
        if (input$genesetOption == "Mouse-GO-CC") load("mouseGOCC.RData")
        if (input$genesetOption == "Mouse-GO-MF") load("mouseGOMF.RData")
      }

      else{
        shiny::showNotification('User Geneset will be used', type = 'message', duration = 10)
      }

      genesets <<- genesets
      # ------ for test
      # q0 <- 0.05

      q0 <- input$qvalueCutoff

      # ------ Create new Waitress
      w <- Waitress$new(selector = NULL, theme = "overlay-radius")

      source('R/GenesetFlush.R')
      source('R/getlgs.R')

      genes <- rownames(CountData)
      genesets <- GenesetFlush(genes, genesets)
      lgs <- getlgs(genesets)

      # ------ Genesetsize Flush

      source('R/GenesetsizeFlush.R')

      genesets <- GenesetsizeFlush(genesets, lgs, input$minGenesetSize, input$maxGenesetSize)
      # ------ For Tests
      # genesets <- GenesetsizeFlush(genesets, lgs, 15, 500)

      # ------ Gene Flush
      source('R/GeneFlush.R')
      remgenes <- GeneFlush(genes, genesets)
      CountData <- CountData[-remgenes, ]

      genesets <<- genesets

      # ------ Background genes
      source('R/getBackgroundGenes.R')
      A <<- getBackgroundGenes(genesets)

      # ------ Calculate t-SNE / U-MAP First
      # require(Matrix)

      source('R/getTU.R')
      # dfobj <- getTU(CountData, 't-SNE')
      dfobj <- getTU(CountData, input$plotOption)
      dfobj <<- dfobj

      cat('getTU Finished\n')

      # ------ Disable radio button
      shinyjs::runjs('$("form p label input").attr("disabled",true)')
      shinyjs::runjs("$('.shinymaterial-slider-minGenesetSize').attr('disabled',true)")
      shinyjs::runjs("$('.shinymaterial-slider-maxGenesetSize').attr('disabled',true)")
      shinyjs::runjs("$('.shinymaterial-slider-qvalueCutoff').attr('disabled',true)")

      cat('running gc\n')
      gc()

      source('R/findSigGenes.R')
      # ------ Find Significant Genes with Fold Change
      if (input$FCoption != "GSVA") {
        # ------ need to build GSVA CASE

        # s <- findSigGenes(CountData, 'median', GroupInfo)
        s <- findSigGenes(CountData, input$FCoption, GroupInfo)
      }

      cat("s Finished\n")

      source('r/findSigGenesGroup.R')

      # ------ Find Significant Genes with findMarkers
      require(dplyr)
      s2 <- findSigGenesGroup(CountData, GroupInfo, q0, TopCutoff = 5)

      rc <- rownames(CountData)

      # ------ free memory to calculate biobj
      rm(CountData)

      # ------ marker l1
      markerl1 <- s2 %>% filter(Top < 10)
      markerl1$Group <- as.factor(markerl1$Group)

      shinyjs::runjs('$(.markerP).show()')

      output$markerL1 <- DT::renderDataTable(
        DT::datatable(markerl1,
          rownames = FALSE,
          filter = "top",
          options = list(
            autoWidth = TRUE,
            dom = "ltp",
            lengthChange = FALSE
          ),
          selection = "none",
        )
      )

      # ------ Hypergeometric pvalue calculation
      lgs <- getlgs(genesets)
      lens = length(s)
      lens100 = round(lens/100)

      source('R/getbiobj.R')
      biobj <- getbiobj(genes, genesets)

      source('R/getHyperPvalue.R')
      pres <- list()

      w$start()
      for (i in 1:lens) {
        if(i %% lens100 == 0) w$inc(1)
        pres[[i]] <- getHyperPvalue(rc[s[[i]]], genesets, A, lgs, q0, biobj)
      }
      w$close()

      pres <<- pres

      # pres : which gene-sets are significant for each cells.

      # ------ CellPathwayDF

      source('R/buildCellPathwayDF.R')
      CellPathwayDF <- buildCellPathwayDF(GroupInfo, pres, genesets)

      # pres2 : for each gene-sets, how many cells are significant that gene-sets.

      cat('pres2\n')

      pres2 <- sort(table(unlist(pres)), decreasing = T)
      names(pres2) <- names(genesets)[as.numeric(names(pres2))]
      pres2 <<- pres2

      source('R/pathwayPvalue.R')
      # 2625*4
      PP <- pathwayPvalue(GroupInfo, pres, pres2, genesets) # qvalue cutoff removed

      CellPathwayDF <- CellPathwayDF %>%
        inner_join(PP) # 1232 * 5


      CellPathwayDF <<- CellPathwayDF

      # l2
      CellMarkers <- data.frame()

      Cells <- sort(unique(GroupInfo))
      Cells <<- Cells

      for (i in 1:length(Cells)) {
        thisCell <- Cells[i]
        thisCellPathways <- CellPathwayDF %>%
          filter(Cell == thisCell) %>%
          select(Geneset)

        thisCellDEs <- s2 %>%
          filter(Group == thisCell) %>%
          select(genes)

        tcd <- thisCellDEs[, 1]
        tcp <- thisCellPathways[, 1]
        tcp <- sapply(tcp, function(i) {
          which(names(genesets) == i)
        }, USE.NAMES = FALSE)

        # ------ Exception handling with gene-sets name with special character
        tcp <- unlist(tcp)

        if (length(tcp) < 1) {
          next
        }
        if (length(tcd) < 1) {
          next
        }
        tcp <- table(unlist(genesets[tcp], use.names = FALSE))
        tcd <- intersect(names(tcp), tcd)
        tcp <- tcp[tcd]

        genes <- names(tcp)
        Count <- unname(tcp)
        additive = data.frame(cbind(genes, Count, Group = thisCell))

        # ------ first add
        if(ncol(CellMarkers)==0) {
          CellMarkers <- rbind(CellMarkers, data.frame(cbind(genes, Count, Group = thisCell), stringsAsFactors = FALSE))
        }
        else{
          if(ncol(CellMarkers)!= ncol(additive)){
            CellMarkers <- rbind(CellMarkers, data.frame(cbind(genes, Count, Group = thisCell), stringsAsFactors = FALSE))
          }
        }

      }

      if(nrow(CellMarkers)){
        CellMarkers <- CellMarkers %>%
          inner_join(s2) %>%
          filter(Top < 10)

        CellMarkers$Group <- as.factor(CellMarkers$Group)
        CellMarkers$Count <- as.numeric(CellMarkers$Count)
        CellMarkers$Group <- as.factor(CellMarkers$Group)

        output$markerL2 <- DT::renderDataTable(
          DT::datatable(CellMarkers,
                        rownames = FALSE,
                        filter = "top",
                        options = list(
                          autoWidth = TRUE,
                          dom = "ltp",
                          lengthChange = FALSE
                        ),
                        selection = "none",
          )
        )
      }
      else{
        print('CellMarker Not Available')
        output$markerL2 <- DT::renderDataTable(
          DT::datatable(s2,
                        rownames = FALSE,
                        filter = "top",
                        options = list(
                          autoWidth = TRUE,
                          dom = "ltp",
                          lengthChange = FALSE
                        ),
                        selection = "none",
          )
        )
      }

      # group 별 significant pathways
      # group 별 DE Genes

      # is counted
      source('R/buildDT.R')
      dtobj <<- buildDT(pres2)

      source('R/briterhex.R')
      source('R/getColv.R')

      # ------ Color define
      colV <- getColv(GroupInfo)

      source('R/getCellHistogram.R')

      CellHistogram <<- getCellHistogram(GroupInfo, colV)

      output$CellBar <- renderHighchart(CellHistogram) # CELL HISTOGRAM

      source("R/getCellPlot.R")

      CellScatter <<- getCellPlot(dfobj, Cells)

      output$CellPlot <- renderHighchart(CellScatter)

      source('R/groupTable.R')

      gt <<- groupTable(pres, genesets, dfobj, pres2)

      # generate dynamic table

      output$dynamicTable <- renderUI({

        numTabs <- length(Cells)
        CardColors <- briterhex(scales::hue_pal()(length(Cells)))

        options(useFancyQuotes = FALSE)

        tagList(
          material_row(
            lapply(1:numTabs, function(i) {

              item <- paste0("$('#dynamicGroupTable", i, " .selected td')[0].innerText")

              material_column(
                # solved material card
                shiny::tags$div(
                  class = "card z-depth-5", # depth = 5, color = null ; color is white
                  style = paste0("border : solid 0.5em ", CardColors[i]), # border color defined
                  shiny::tags$div(
                    class = "card-content",
                    shiny::tags$span(class = "card-title", Cells[i]), # title
                    shiny::tags$div(class = "divider"), # divider = TRUE
                    DT::dataTableOutput(
                      paste0("dynamicGroupTable", i),
                      width = "100%",
                      height = "500px"
                    ),
                    solvedButton(
                      inputId = paste0("toSortButton", i),
                      label = "Select",
                      onClick = HTML(
                        paste0(
                          sortItem(paste0(item, " + ' @", Cells[i], "'"), "mysortableCell"),
                          "; $('#toSortButton", i, "').attr('disabled', true);"
                          )
                      ),
                      style = "position:absolute; top:1em; right:1em;background-color: #1976d2"
                    )
                  )
                ),
                width = 4 # maximum 3 table in column
              )
            })
          )
        )
      })

      # ------ fill dynamic table
      # ordered by Count , not length;

      for (i in 1:length(Cells)) {
        t <- paste0(
          "output$dynamicGroupTable", i, " = DT::renderDataTable(",
          "DT::datatable(",
            "CellPathwayDF[which(CellPathwayDF[,1]==Cells[", i, "]),-1],", # removed group column
            "options = list(",
              "dom = 'ltp',",
              "scroller = TRUE,",
              "scrollX = TRUE,",
              "autoWidth = TRUE,",
              "lengthChange = FALSE,",
              "order = list(list(1,'desc'))",
              "),",
            "rownames = FALSE,",
            "selection = 'single')",
          ")"
          )
        eval(parse(text = t))
      }
      # StartEnrich Finished

      print(proc.time() - pt)

    })

    # draw frequently colored images
    observeEvent(input$freqbtn, {
      if (input$freqbtn == 0) { # prevent initial click state
        return(NULL)
      }

      res <- c()

      # select each cell's frequent pathway
      for (i in 1:length(Cells)) {
        thisCellData <- CellPathwayDF %>% dplyr::filter(Cell == Cells[i])
        if (nrow(thisCellData) >= 1) {
          thisCellData <- thisCellData %>% top_n(1, wt = Count)
          if (nrow(thisCellData) >= 1) {
            thisCellData <- thisCellData %>% top_n(-1, wt = Length)
            if (nrow(thisCellData) >= 1) {
              thisCellData <- thisCellData %>% top_n(-1, wt = Qvalue)
              if (nrow(thisCellData) >= 1) {
                thisCellData <- thisCellData %>% top_n(1)
              }
            }
          }
          res <- c(res, paste0(thisCellData$Geneset, " @", thisCellData$Cell))
        }
      }
      output$CellPlot <- renderHighchart(emphasize(FALSE, res, dfobj, Cells, pres, genesets))
    })

    # draw significant colored images
    observeEvent(input$sigbtn, {
      if (input$sigbtn == 0) { # prevent default click state
        return(NULL)
      }

      res <- c()

      for (i in 1:length(Cells)) {
        thisCell <- Cells[i]

        thisCellData <- CellPathwayDF %>% dplyr::filter(Cell == thisCell)
        if (nrow(thisCellData) >= 1) {
          thisCellData <- thisCellData %>% top_n(-1, wt = Qvalue)

          if (nrow(thisCellData) >= 1) {
            thisCellData <- thisCellData %>% top_n(-1, wt = Length)

            if (nrow(thisCellData) >= 1) {
              thisCellData <- thisCellData %>% top_n(1, wt = Count)
              if (nrow(thisCellData) >= 1) {
                thisCellData <- thisCellData %>% top_n(1)
              }
            }
          }

          res <- c(res, paste0(thisCellData$Geneset, " @", thisCellData$Cell))
        }
      }

      output$CellPlot <- renderHighchart(emphasize(FALSE, res, dfobj, Cells, pres, genesets))
    })

    # draw gray colored images
    observeEvent(input$graybtn, {
      if (input$graybtn == 0) { # prevent default click state
        return(NULL)
      }

      dfobj_new <- dfobj
      dfobj_new$col <- '#8395a7'

      grayImage <- hchart(
          dfobj_new,
          type = 'scatter',
          hcaes(x = x, y = y, color = col)
        ) %>%
        hc_tooltip(FALSE) %>%
        hc_exporting(enabled = TRUE)

      output$CellPlot <- renderHighchart(grayImage)
    })

    # draw group colored images
    observeEvent(input$colorbtn, {
      if (input$colorbtn == 0) { # prevent default click state
        return(NULL)
      }

      # UniqueCol <- briterhex(scales::hue_pal()(length(Cells)))
      # names(UniqueCol) <- Cells

      # colV <- unname(UniqueCol[dfobj$col])

      # CellScatter <<- getCellPlot(dfobj, Cells)

      output$CellPlot <- renderHighchart(CellScatter)

    })

    # Emphasize with order
    observeEvent(input$OrderEmphasize, {
      if (input$OrderEmphasize == 0) { # prevent default click state
        return(NULL)
      }

      output$CellPlot <- renderHighchart(emphasize(TRUE, input$sortList, dfobj, Cells, pres, genesets))
    })

    # Emphasize without order
    observeEvent(input$Emphasize, {
      if (input$Emphasize == 0) { # prevent default click state
        return(NULL)
      }

      output$CellPlot <- renderHighchart(emphasize(FALSE, input$sortList, dfobj, Cells, pres, genesets))
    })

    # clear timelist in Cell tab
    observeEvent(input$ClearList, {
      if (input$ClearList == 0) { # prevent default click state
        return(NULL)
      }
      shinyjs::runjs(
        code = paste0(
          '$("#mysortableCell .rank-list-item").remove();',
          '$("#dynamicTable button").attr("disabled",false)'
        )
      )
    })
  }

  source('R/CellEnrichUI.R')
  ui <- CellEnrichUI()

  shiny::shinyApp(ui, server, options = list(launch.browser = TRUE))
}
