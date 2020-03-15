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

CellEnrich <- function(CountData, GroupInfo) {

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


    pathwayPvalue <- function(Cells, q0 = 0.1) {
      res <- c()
      Cells <- unique(GroupInfo) # 16cell1, 16cell2,
      total <- length(GroupInfo)

      for (i in 1:length(Cells)) {
        thisCell <- Cells[i]
        thisCellIdx <- which(GroupInfo == thisCell)

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

    A <- ""

    ### CODES

    # variable initialize

    ggobj2 <- ggobj <- dtobj <- dfobj <- pres <- CellPathwayDF <- ""
    gt <- pres2 <- Cells <- ""

    # start

    observeEvent(input$StartCellEnrich, {
      shinyjs::hide("StartCellEnrich")

      if (input$genesetOption == "Curated") load("c2v7.RData")

      if (input$genesetOption == "GeneOntology") load("c5v7.RData")

      if (input$genesetOption == "KEGG") load("keggv7.RData")

      q0 <- input$qvalueCutoff

      w <- Waitress$new(selector = NULL, theme = "overlay")$start()

      # geneset background flushing
      genes <- rownames(CountData)
      for (i in 1:length(genesets)) {
        genesets[[i]] <- intersect(genesets[[i]], genes)
      }

      lgs <- sapply(1:length(genesets), function(i) {
        length(genesets[[i]])
      })

      # gene geneset flushing
      gsgenes <- unique(unlist(genesets))
      remgenes <- sapply(setdiff(genes, gsgenes), function(i) {
        which(i == rownames(CountData))
      }, USE.NAMES = FALSE)
      genes <- rownames(CountData)

      CountData <- CountData[-remgenes, ]

      gidx <- 1:length(genes)
      names(gidx) <- genes

      # calculate background intersection object
      getbiobj <- function(genes, genesets) {
        res <- matrix(0, length(genes), length(genesets))
        for (i in 1:length(genesets)) {
          res[unname(gidx[genesets[[i]]]), i] <- 1
        }
        rownames(res) <- genes
        colnames(res) <- names(genesets)
        return(res)
      }

      biobj <- getbiobj(genes, genesets)

      # genesets = genesets[intersect(which(lgs >= 15), which(lgs <= 500))]

      genesets <- genesets[intersect(which(lgs >= input$minGenesetSize), which(lgs <= input$maxGenesetSize))]
      genesets <<- genesets
      cat("minimum gene-set size :", input$minGenesetSize, "\n")
      cat("maximum gene-set size :", input$maxGenesetSize, "\n")

      cat(length(genesets), "genesets\n")

      A <<- length(unique(unlist(genesets))) # Background genes

      shinyjs::runjs('$("form p label input").attr("disabled",true)') # radio button disable
      shinyjs::runjs("$('.shinymaterial-slider-minGenesetSize').attr('disabled',true)")
      shinyjs::runjs("$('.shinymaterial-slider-maxGenesetSize').attr('disabled',true)")
      shinyjs::runjs("$('.shinymaterial-slider-qvalueCutoff').attr('disabled',true)")



      v <- CountData

      if (input$FCoption != "GSVA") {
        # s <- findSigGenes(v, 'median')
        s <- findSigGenes(v, input$FCoption)
        names(s) <- GroupInfo # oocyte 1, oocyte 2
      }

      s2 <- findSigGenesGroup(v, GroupInfo, q0, TopCutoff = 5)
      s2$FDR <- round(as.numeric(s2$FDR), 6)

      markerl1 <- s2 %>% filter(Top < 10)
      markerl1$Group <- as.factor(markerl1$Group)

      shinyjs::runjs('$(.markerP).show()')

      # marker L1
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

      lgs <- sapply(1:length(genesets), function(i) {
        length(genesets[[i]])
      })

      # for test

      pres <- list()
      for (i in 1:length(s)) {
        w$inc(1 / length(s))
        pres[[i]] <- getHyperPvalue(s[[i]], genesets, A, lgs, q0, biobj)
      }

      w$close()
      pres <<- pres

      # pres : which gene-sets are significant for each cells.

      CellPathwayDF <- data.frame(stringsAsFactors = FALSE)
      Cells <<- unique(GroupInfo) # 16cell1, 16cell2,
      for (i in 1:length(Cells)) {
        thisCell <- Cells[i]
        tt <- table(unlist(pres[which(thisCell == GroupInfo)]))
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
      PP <- pathwayPvalue(Cells, q0)

      CellPathwayDF <- CellPathwayDF %>%
        inner_join(PP) # 1232 * 5

      CellPathwayDF <<- CellPathwayDF

      # l2
      CellMarkers <- data.frame()

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
        CellMarkers <- rbind(CellMarkers, data.frame(cbind(genes, Count, Group = thisCell), stringsAsFactors = FALSE))
      }

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

      # group 별 significant pathways
      # group 별 DE Genes

      # is counted


      dtobj <<- buildDT(pres2)
      # GroupInfo -> oocyte1, oocyte2, ...

      if (input$plotOption == "t-SNE") {
        tsneE <- Rtsne(t(v), check_duplicates = FALSE, perplexity = 15)
        dfobj <- data.frame(tsneE$Y, col = GroupInfo, stringsAsFactors = FALSE)
      }

      if (input$plotOption == "U-MAP") {
        umapE <- uwot::umap(t(v), fast_sgd = TRUE)
        dfobj <- data.frame(umapE, col = GroupInfo, stringsAsFactors = FALSE)
      }

      UniqueCol <- briterhex(scales::hue_pal()(length(Cells)))
      names(UniqueCol) <- Cells

      x <- c()
      y <- c()
      for (i in 1:length(Cells)) {
        x[i] <- Cells[i]
        y[i] <- length(which(GroupInfo == Cells[i]))
      }

      colV <- unname(UniqueCol[x])
      names(colV) <- Cells

      ggobjdf <- data.frame(x, y, stringsAsFactors = FALSE)
      colnames(ggobjdf) <- c("x", "y")

      ggobj <<- ggplot(ggobjdf, aes(x = x, y = y)) +
        geom_bar(stat = "identity", fill = colV) # cell histogram

      colnames(dfobj) <- c("x", "y", "col")
      dfobj <<- dfobj

      UniqueCol <- briterhex(scales::hue_pal()(length(Cells)))
      names(UniqueCol) <- Cells

      colV <- unname(UniqueCol[dfobj$col])

      # scatter plot
      ggobj2 <<- ggplot(dfobj, aes(x = x, y = y)) +
        geom_point(colour = colV) # +
      # scale_color_manual(values = briterhex(scales::hue_pal()(length(unique(dfobj$col)))))

      output$CellPlot <- shiny::renderPlot(ggobj2) # CELL SCATTERPLOT
      output$CellBar <- shiny::renderPlot(ggobj) # CELL HISTOGRAM
      # output$tab <- DT::renderDataTable(dtobj) # CELL DATATABLE

      gt <<- groupTable()

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
                      style = "position:absolute; top:1em; right:1em;"
                    )
                  )
                ),
                width = 4 # maximum 3 table in column
              )
            })
          )
        )
      })

      # fill dynamic table
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

      output$CellPlot <- shiny::renderPlot(emphasize(FALSE, res, dfobj, Cells, pres))
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

      output$CellPlot <- shiny::renderPlot(emphasize(FALSE, res, dfobj, Cells, pres))
    })

    # draw gray colored images
    observeEvent(input$graybtn, {
      if (input$graybtn == 0) { # prevent default click state
        return(NULL)
      }

      gobj <- ggplot(dfobj, aes(x = x, y = y)) +
        geom_point(colour = "gray")

      output$CellPlot <- shiny::renderPlot(gobj)
    })

    # draw group colored images
    observeEvent(input$colorbtn, {
      if (input$colorbtn == 0) { # prevent default click state
        return(NULL)
      }

      UniqueCol <- briterhex(scales::hue_pal()(length(Cells)))
      names(UniqueCol) <- Cells

      colV <- unname(UniqueCol[dfobj$col])

      # scatter plot
      gobj =  ggplot(dfobj, aes(x = x, y = y)) +
        geom_point(colour = colV)

      output$CellPlot <- shiny::renderPlot(gobj)
    })

    # Emphasize with order
    observeEvent(input$OrderEmphasize, {
      if (input$OrderEmphasize == 0) { # prevent default click state
        return(NULL)
      }

      output$CellPlot <- shiny::renderPlot(emphasize(TRUE, input$sortList, dfobj, Cells, pres))
    })

    # Emphasize without order
    observeEvent(input$Emphasize, {
      if (input$Emphasize == 0) { # prevent default click state
        return(NULL)
      }

      output$CellPlot <- shiny::renderPlot(emphasize(FALSE, input$sortList, dfobj, Cells, pres))
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

  shiny::shinyApp(ui, server, options = list(launch.browser = TRUE))
}

CellEnrichUI <- function() {

  material_page(
    shinyjs::useShinyjs(),

    # dynamic datatable full width

    tags$head(tags$style(type = "text/css", ".display.dataTable.no-footer{width : 100% !important;}")),
    tags$head(tags$style(type = "text/css", ".markerP {display : none;}")),

    # waitress declare
    use_waitress(color = "#697682", percent_color = "#333333"),

    title = paste0(
      "CellEnrich ",
      "<a href = 'https://github.com/jhk0530/cellenrich' target = '_blank'> ", # github link
      "<i class='material-icons' style = 'font-size:1.3em;'>info</i> </a>" # icon tag
    ),
    nav_bar_fixed = FALSE,
    nav_bar_color = "light-blue darken-1",
    font_color = "#ffffff",
    include_fonts = FALSE,
    include_nav_bar = TRUE,
    include_icons = FALSE,

    # CellEnrich options
    material_row(
      material_column(
        material_card(
          title = "Options",
          divider = TRUE,
          material_radio_button(
            input_id = "FCoption",
            label = "Select FoldChange Option",
            choices = c("median", "mean", "zero", "GSVA"),
            selected = "median"
          ),
          material_radio_button(
            input_id = "plotOption",
            label = "Select Plot Option",
            choices = c("t-SNE", "U-MAP"),
            selected = "t-SNE"
          ),
          material_radio_button(
            input_id = "genesetOption",
            label = "Select Gene-set",
            choices = c(
              "Curated", # c2
              "GeneOntology", # c5 = GO
              "KEGG"), # KEGG
            selected = "Curated"
          ),
          material_slider(
            input_id = "minGenesetSize",
            label = "Minimum Gene-set Size",
            min_value = 10,
            max_value = 30,
            initial_value = 15,
            step_size = 5
          ),
          material_slider(
            input_id = "maxGenesetSize",
            label = "Maximum Gene-set Size",
            min_value = 250,
            max_value = 750,
            initial_value = 500,
            step_size = 5
          ),
          material_slider(
            input_id = "qvalueCutoff",
            label = "Q-value threshold",
            min_value = 0,
            max_value = 0.25,
            initial_value = 0.1,
            step_size = 0.05
          ),
          solvedButton(
            inputId = "StartCellEnrich",
            label = "Start CellEnrich",
            style = "margin-left:45%;",
            onClick = 'console.log("CellEnrich");'
          ),
          depth = 3
        ),
        width = 6, # Centered Layout
        offset = 3
      )
    ),

    # tSNE/UMAP plot
    material_row(
      material_column(
        material_card(
          title = "Group plot / Distribution",
          depth = 3,
          material_column(
            plotOutput("CellPlot", height = "700px"),
            width = 6
          ),
          material_column(
            plotOutput("CellBar", height = "700px"), # cell distribution
            width = 6
          ),
          material_button("colorbtn", "toColor", icon = "color_lens", color = "blue lighten-1"),
          material_button("graybtn", "toGray", icon = "clear", color = "blue lighten-1"),
          material_button("freqbtn", "Frequent", icon = "grain", color = "blue lighten-1"),
          material_button("sigbtn", "Significant", icon = "grade", color = "blue lighten-1")
        ),
        width = 12
      ),
      style = "margin : 1em"
    ),

    # marker table
    material_row(
      material_card(
        title = "MarkerGenes",
        p("DE from each Cell specific", class = 'markerP'),
        DT::dataTableOutput("markerL1"),
        p("DE - Pathway from each Cell specific", class = 'markerP'),
        DT::dataTableOutput("markerL2")
      ),
      style = "margin : 1em"
    ),

    # emphasize tables
    material_row(
      material_card(
        title = "",
        # material_button("callSortFunction", "Activate", icon = "play_arrow", color = "pink"),
        material_card(
          title = "Pathway Emphasize", divider = TRUE,
          tags$h3("To be recognized by application, Please move element's position"),
          rank_list(text = "Pathways", labels = "Please Clear First", input_id = "sortList", css_id = "mysortableCell"),
          material_button("OrderEmphasize", "Emphasize with Order", icon = "timeline"),
          material_button("Emphasize", "Emphasize without Order", icon = "bubble_chart"),
          material_button("ClearList", "Clear List", icon = "clear_all"),
        ),
        uiOutput("dynamicTable"),
        depth = 3
      ),
      style = "margin : 1em"
    )
  )
}
