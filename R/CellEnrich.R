GenesetFlush <- function(genes, genesets) {
  cat("GenesetFlush\n")
  for (i in 1:length(genesets)) {
    genesets[[i]] <- intersect(genesets[[i]], genes)
  }
  return(genesets)
}

getlgs <- function(genesets) {
  n <- names(genesets)
  cat("getlgs\n")
  sapply(1:length(n), function(i) {
    length(genesets[[n[i]]])
  })
}

GenesetsizeFlush <- function(genesets, lgs, minsize, maxsize) {
  cat("GenesetSizeFlush\n")
  genesets[intersect(which(lgs >= minsize), which(lgs <= maxsize))]
}

GeneFlush <- function(genes, genesets) {
  cat("GeneFlush\n")
  gsgenes <- unique(unlist(genesets))
  remgenes <- sapply(setdiff(genes, gsgenes), function(i) {
    which(i == genes)
  }, USE.NAMES = FALSE)
  return(remgenes)
}

getBackgroundGenes <- function(genesets) {
  cat("getBackgroundGenes\n")
  length(unique(unlist(genesets)))
}

getTU <- function(CountData, GroupInfo, plotOption) {
  cat("transpose in getTU started\n")
  CountData <- as.matrix(Matrix::t(CountData))
  cat("t-SNE / U-MAP started\n")
  if (plotOption == "t-SNE") {
    if (nrow(CountData) < 500) {
      tsneE <- Rtsne::Rtsne(CountData, check_duplicates = FALSE, perplexity = 15)
    } else {
      tsneE <- Rtsne::Rtsne(CountData, check_duplicates = FALSE, perplexity = 15, partial_pca = TRUE) # 15 seconds
    }
    dfobj <- data.frame(tsneE$Y, col = GroupInfo, stringsAsFactors = FALSE)
  }

  if (plotOption == "U-MAP") {
    umapE <- uwot::umap(CountData, fast_sgd = TRUE) # 55 seconds
    dfobj <- data.frame(umapE, col = GroupInfo, stringsAsFactors = FALSE)
  }
  colnames(dfobj) <- c("x", "y", "col")
  return(dfobj)
}

findSigGenes <- function(v, method = "median", Name) {
  if (!method %in% c("median", "mean", "zero", "GSVA")) stop("wrong method")
  # it's already matrix
  # v <- as.matrix(v)
  cat("findSigGenes started\n")

  rownames(v) <- colnames(v) <- NULL

  res <- list()

  # 100 : 0.8 seconds
  # 200 : 1.12
  # 250 : 1.36. find 20000 gene will take less than minute

  cat("scaling\n")

  idx <- floor(nrow(v) / 250)
  v2 <- c()
  for (i in 1:idx) {
    thisIdx <- 1:250 + 250 * (i - 1)
    vv <- as.matrix(v[thisIdx, ]) + 1
    meds <- apply(vv, 1, median)
    vv <- as(log(sweep(vv, 1, meds, "/")), "dgCMatrix")
    v2 <- rbind(v2, vv) # use rbind, not assign
  }

  if (nrow(v) %% 250 != 0) {
    thisIdx <- (idx * 250 + 1):nrow(v)
    vv <- as.matrix(v[thisIdx, ]) + 1
    meds <- apply(vv, 1, median)
    vv <- log(sweep(vv, 1, meds, "/"))
    v2 <- rbind(v2, vv)
  }

  v <- v2
  rm(v2)

  cat("define Lists\n")

  if (method == "median") {
    for (i in 1:ncol(v)) {
      res[[i]] <- which(v[, i] > median(v[, i]))
    }
  }
  if (method == "zero") {
    for (i in 1:ncol(v)) {
      res[[i]] <- which(v[, i] > 0)
    }
  }
  if (method == "mean") {
    for (i in 1:ncol(v)) {
      res[[i]] <- which(v[, i] > mean(v[, i]))
    }
  }

  names(res) <- Name
  return(res)
}

findSigGenesGroup <- function(Count = NULL, ClustInfo = NULL, q0 = 0.1, TopCutoff = 5) {
  require(scran)

  if (is.null(Count)) stop("Count must given")
  if (is.null(ClustInfo)) stop("ClustInfo must given")

  GrpRes <- scran::findMarkers(x = as.matrix(Count), ClustInfo, test = "wilcox", direction = "up")
  Grp <- unique(ClustInfo)

  res <- data.frame(stringsAsFactors = FALSE)

  for (i in 1:length(Grp)) {
    G <- data.frame(
      genes = rownames(GrpRes[[1]]),
      Group = Grp[i],
      GrpRes[[1]],
      row.names = NULL,
      stringsAsFactors = FALSE
    ) %>%
      select(Group, Top, genes, FDR) %>%
      filter(FDR <= q0) %>%
      # filter(Top <= TopCutoff) %>%
      arrange(FDR)

    res <- rbind(res, G)
  }
  res$genes <- as.character(res$genes)
  res$Group <- as.character(res$Group)
  res$FDR <- round(as.numeric(res$FDR), 6)
  return(res)
}

getbiobj <- function(genes, genesets) {
  gidx <- 1:length(genes)
  names(gidx) <- genes

  res <- matrix(0, length(genes), length(genesets))
  for (i in 1:length(genesets)) {
    res[unname(gidx[genesets[[i]]]), i] <- 1
  }

  rownames(res) <- genes
  colnames(res) <- names(genesets)
  return(res)
}

getHyperPvalue <- function(genes, genesets, A, lgs, q0, biobj) {
  lg <- length(genes)

  if (lg == 0) {
    return(integer(0))
  }
  gidx <- 1:length(genes)
  names(gidx) <- genes

  if (length(genes) == 1) {
    biobj <- biobj[genes, ]
  }
  else {
    biobj <- unname(colSums(biobj[genes, ]))
  }

  # ------

  pv <- sapply(1:length(genesets), function(i) {
    q <- biobj[i] # selected white ball
    m <- lgs[i] # white ball
    n <- A - m # black ball
    k <- lg # selected ball
    1 - phyper(q - 1, m, n, k)
  })
  # names(pv) <- names(genesets)
  return(pv)
  # return(which(pv < q0))
}

buildCellPathwayDF <- function(GroupInfo, pres, genesets) {
  cat("buildCellPathwayDF\n")
  Cells <- unique(GroupInfo)
  CellPathwayDF <- data.frame(stringsAsFactors = FALSE)

  for (i in 1:length(Cells)) {
    thisCell <- Cells[i]
    tt <- table(unlist(pres[which(thisCell == GroupInfo)]))

    if (nrow(tt)) {
      CellPathwayDF <- rbind(CellPathwayDF, cbind(thisCell, names(tt), unname(tt)))
    }
  }

  colnames(CellPathwayDF) <- c("Cell", "Geneset", "Count")

  CellPathwayDF$Cell <- as.character(CellPathwayDF$Cell)
  CellPathwayDF$Geneset <- names(genesets)[as.numeric(as.character(CellPathwayDF$Geneset))]
  CellPathwayDF$Count <- as.numeric(as.character(CellPathwayDF$Count))

  # ------ add length column

  # Length <- getlgs(CellPathwayDF$Geneset)
  Length <- getlgs(genesets[as.character(CellPathwayDF$Geneset)])
  CellPathwayDF <- cbind(CellPathwayDF, Length)

  # ------ select genesets with count > 1

  CellPathwayDF <- CellPathwayDF %>%
    dplyr::filter(Count > 1)

  return(CellPathwayDF)
}

pathwayPvalue <- function(GroupInfo, pres, pres2, genesets) {
  cat("pathwayPvalue\n")
  res <- c()
  Cells <- unique(GroupInfo)
  total <- length(GroupInfo)

  for (i in 1:length(Cells)) {
    thisCell <- Cells[i]
    thisCellIdx <- which(GroupInfo == thisCell)

    k <- length(thisCellIdx)

    thisCellPathways <- table(unlist(pres[thisCellIdx]))
    if (nrow(thisCellPathways) < 1) {
      next
    }
    pv <- c()

    for (j in 1:length(thisCellPathways)) {
      thisPathway <- names(thisCellPathways)[j]
      q <- unname(thisCellPathways)[j] # selected white ball
      m <- pres2[names(genesets)[as.numeric(thisPathway)]] # total white ball
      pv[j] <- 1 - phyper(q - 1, m, total - m, k)
    }
    names(pv) <- names(genesets)[as.numeric(names(thisCellPathways))]

    res <- rbind(res, cbind(thisCell, names(pv), unname(pv)))
  }
  res <- data.frame(res, stringsAsFactors = FALSE)
  colnames(res) <- c("Cell", "Geneset", "Qvalue")

  res$Cell <- as.character(res$Cell)
  res$Geneset <- as.character(res$Geneset)
  res$Qvalue[which(res$Qvalue <= 1e-20)] <- 1e-20
  res$Qvalue <- -log10(as.numeric(as.character(res$Qvalue)))

  return(res)
}

# pres : which gene-sets are significant for each cells.
# pres2 : for each gene-sets, how many cells are significant that gene-sets.

# 전체 그룹에서 유의한 회수 20 # pres2[genesets[i]]
# 특정 그룹에서 유의한 회수 6 # pres2[thiscellidx]

# 전체 그룹 Cell 수 : N
# 특정 그룹 Cell 수 : K

# Group_specific_OR = (6/K) / (14/N)

getOddRatio <- function(GroupInfo, pres, pres2, genesets, ratio) {
  cat("getOddRatio\n")

  res <- data.frame(stringsAsFactors = FALSE)
  Cells <- unique(GroupInfo)
  total <- length(GroupInfo)
  for (i in 1:length(Cells)) {
    thisCell <- Cells[i]
    thisCellIdx <- which(GroupInfo == thisCell)
    OR <- unname(sapply(1:length(genesets), function(k) {
      B <- table(unlist(pres[thisCellIdx]))[as.character(k)] # 특정 Cell에서 유의한 회수
      if (is.na(B)) {
        return(0)
      }
      if (B < length(thisCellIdx) * ratio) {
        return(0)
      }
      A <- pres2[names(genesets)[k]] # 전체 Cell에서 유의한 회수
      if (is.na(A)) {
        return(0)
      }
      N <- total # 전체 Cell 수

      K <- length(thisCellIdx)

      return((B / K) / (A / N))
    }))

    OR <- round(OR, 4)
    # Cell, Geneset, OR
    res <- rbind(
      res,
      data.frame(
        Cell = as.character(thisCell),
        Geneset = as.character(names(genesets)),
        OddRatio = as.numeric(OR), stringsAsFactors = FALSE
      )
    )
  }

  colnames(res) <- c("Cell", "Geneset", "OddRatio")

  res <- res %>% filter(OddRatio > 1)

  return(res)
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

briterhex <- function(colors) {
  if (length(colors) == 0) {
    stop("Length of color should be larger than zero --- briterhex")
  }
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

getColv <- function(GroupInfo) {
  Cells <- unique(sort(GroupInfo))

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
  return(colV)
}

getCellHistogram <- function(GroupInfo, colV) {
  cat("getCellHistogram\n")
  # require(ggplot2)
  require(highcharter)
  Cells <- unique(sort(GroupInfo))

  x <- c()
  y <- c()
  for (i in 1:length(Cells)) {
    x[i] <- Cells[i]
    y[i] <- length(which(GroupInfo == Cells[i]))
  }
  colV <- unname(colV)

  hc <- highchart() %>%
    hc_chart(type = "column", legend = list(enabled = FALSE)) %>%
    hc_title(text = "Groups") %>%
    hc_xAxis(categories = x) %>%
    hc_plotOptions(grouping = FALSE) %>%
    hc_add_series(data = y, colorByPoint = TRUE, showInLegend = FALSE, name = "Count") %>%
    hc_colors(colV) %>%
    hc_exporting(enabled = TRUE)
  return(hc)
}

getCellPlot <- function(dfobj, Cells) {
  cat("getCellPlot\n")
  require(ggplot2)
  # require(highcharter)

  colnames(dfobj) <- c("x", "y", "col")
  dfobj <<- dfobj


  UniqueCol <- briterhex(scales::hue_pal()(length(Cells)))
  names(UniqueCol) <- Cells

  colV <- unname(UniqueCol[dfobj$col])

  cat("\n")
  return(
    ggplot(dfobj, aes(x = x, y = y)) +
      geom_point(colour = colV)
  )

  # highcharter cancel

  # dfobj$col <- unname(UniqueCol[dfobj$col])
  # rownames(dfobj) = NULL

  # hchart(dfobj, type = 'scatter', hcaes(x = x, y = y, color = col)) %>%
  # hc_add_series(data = dfobj[3000,nrow(dfobj)], colorByPoint = TRUE, showInLegend = FALSE) %>%
  # hc_colors(colV)
  # hc_tooltip(FALSE) %>%
  # hc_exporting(enabled = TRUE)

  # return(hc)
}

groupTable <- function(pres, genesets, dfobj, pres2) {
  cat("groupTable\n")
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
    if (length(pathways) < 1) {
      next
    }
    # what genesets are enriched per each group.

    k <- sum(pathways) # selected ball

    gt <- sapply(1:length(pathways), function(j) {
      q <- pathways[j] # selected white ball, 1
      m <- unname(pres2Idx[names(pathways[j])]) # total white ball, 28
      # n <- tot - m # total black ball
      round(1 - phyper(q - 1, m, tot - m, k), 4)
    })
    gt <- gt[which(gt < 0.25)] # pvalue 0.25
    if (length(gt)) {
      res <- rbind(res, cbind(groups[i], names(gt), unname(gt)))
    }
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


CellEnrichUI <- function() {
  require(shinymaterial)
  require(highcharter)
  material_page(
    shinyjs::useShinyjs(),

    # dynamic datatable full width

    tags$head(tags$style(type = "text/css", ".display.dataTable.no-footer{width : 100% !important;}")),

    # waitress declare
    use_waitress(color = "#697682", percent_color = "#333333"),
    title = paste0(
      "CellEnrich ",
      "<a href = 'https://github.com/jhk0530/cellenrich' target = '_blank'> ", # github link
      "<i class='material-icons' style = 'font-size:1.3em;'>info</i> </a>" # icon tag
    ),
    nav_bar_fixed = FALSE,
    nav_bar_color = "blue darken-2",
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
          style = "border : solid 0.5em #1976d2",
          material_row(
            material_column(
              material_card(
                material_radio_button(
                  input_id = "FCoption",
                  label = "Select FoldChange Option",
                  choices = c("median", "mean", "zero"),
                  selected = "median",
                  color = "#1976d2"
                ),
                material_radio_button(
                  input_id = "plotOption",
                  label = "Select Plot Option",
                  choices = c("t-SNE", "U-MAP"),
                  selected = "t-SNE",
                  color = "#1976d2"
                )
              ),
              width = 4
            ),
            material_column(
              material_card(
                material_number_box(
                  input_id = "minGenesetSize",
                  label = "Minimum Gene-set Size",
                  min_value = 10,
                  max_value = 30,
                  initial_value = 15,
                  step_size = 5
                ),
                material_number_box(
                  input_id = "maxGenesetSize",
                  label = "Maximum Gene-set Size",
                  min_value = 250,
                  max_value = 750,
                  initial_value = 500,
                  step_size = 5
                ),
                material_number_box(
                  input_id = "ORratio",
                  label = "OddRatio Frequency",
                  min_value = 0,
                  max_value = 0.5,
                  initial_value = 0.1,
                  step_size = 0.05
                ),
                material_number_box(
                  input_id = "qvalueCutoff",
                  label = "Q-value threshold",
                  min_value = 0,
                  max_value = 0.25,
                  initial_value = 0.05,
                  step_size = 0.01
                )
              ),
              width = 4
            ),
            material_column(
              material_card(
                material_radio_button(
                  input_id = "genesetOption",
                  label = "Select Gene-set",
                  color = "#1976d2",
                  choices = c(
                    "User-Geneset",
                    "Human-Curated", # c2
                    "Human-KEGG", # KEGG
                    "Human-GO",
                    "Human-GO-BP",
                    "Human-GO-CC",
                    "Human-GO-MF",
                    "Mouse-KEGG", # Mouse
                    "Mouse-GO",
                    "Mouse-GO-BP",
                    "Mouse-GO-CC",
                    "Mouse-GO-MF"
                  ),
                  selected = NULL
                )
              ),
              width = 4
            )
          ),
          solvedButton(
            inputId = "StartCellEnrich",
            label = "Start CellEnrich",
            style = "margin-left:45%; background-color: #1976d2",
            onClick = 'console.log("CellEnrich");'
          ),
          depth = 3
        ),
        width = 6,
        offset = 3 # center half layout
      )
    ),

    # tSNE/UMAP plot
    material_row(
      material_column(
        material_card(
          title = "Group plot / Distribution",
          depth = 3,
          material_row(
            material_column(
              plotOutput("CellPlot", height = "700px"),
              width = 6
            ),
            material_column(
              highchartOutput("CellBar", height = "700px"), # cell distribution
              width = 6
            )
          ),
          material_row(
            shiny::downloadButton("imgdn", "Save Scatter Plot", icon = "save", style = "background-color : #616161 !important")
          ),
          material_row(
            material_card(
              title = "information",
              DT::dataTableOutput("legendTable"),
              shiny::downloadButton("legenddn", "Save Legend", icon = "save", style = "background-color : #616161 !important; display:none;")
            )
          ),
          material_row(
            material_button("colorbtn", "toColor", icon = "color_lens", color = "blue darken-2"),
            material_button("freqbtn", "Frequent", icon = "grain", color = "blue darken-2"),
            material_button("sigbtn", "Significant", icon = "grade", color = "blue darken-2")
          )
        ),
        width = 12
      ),
      style = "margin : 1em; border : solid 0.5em #1976d2"
    ),

    # marker table
    material_row(
      material_card(
        title = "MarkerGenes",
        material_row(
          material_column(
            material_card(
              title = "DE from each Cell specific",
              DT::dataTableOutput("markerL1")
            ),
            width = 6
          ),
          material_column(
            material_card(
              title = "DE - Pathway from each Cell specific",
              DT::dataTableOutput("markerL2")
            ),
            width = 6
          )
        )
      ),
      style = "margin : 1em; border : solid 0.5em #1976d2"
    ),

    # emphasize tables
    material_row(
      material_card(
        title = "",
        material_card(
          title = "Pathway biplot", divider = TRUE,
          material_row(
            material_column(
              plotOutput("biPlot", height = "700px"),
              width = 9
            ),
            material_column(
              material_row(
                numericInput("biFont", label = "Label Size", value = 3, min = 1, max = 5, step = 1),
                numericInput("biX", label = "X area", value = 5, min = 1, max = 10, step = 1),
                numericInput("biY", label = "Y area", value = 5, min = 1, max = 10, step = 1),
                material_row(
                  material_button("relFreq", "Biplot with Relative Freq.", color = "blue darken-2")
                ),
                material_row(
                  material_button("absFreq", "Biplot with Absolute Freq.", color = "blue darken-2")
                ),
                material_row(
                  material_button("refreshBiplot", "Refresh Biplot Download Image", color = "blue darken-2")
                ),
                shiny::downloadButton("biplotdn", "Save Biplot", icon = "save", style = "background-color : #616161 !important")
              ),
              width = 3
            )
          ),
        ),
        material_card(
          title = "Pathway Emphasize", divider = TRUE,
          tags$h3("To be recognized by application, Please move element's position"),
          rank_list(text = "Pathways", labels = "Please Clear First", input_id = "sortList", css_id = "mysortableCell"),
          material_row(
            material_button("OrderEmphasize", "Emphasize with Order", icon = "timeline", color = "blue darken-2"),
            material_button("Emphasize", "Emphasize without Order", icon = "bubble_chart", color = "blue darken-2"),
            material_button("ClearList", "Clear List", icon = "clear_all", color = "blue darken-2")
          ),
          material_row(
            shiny::downloadButton("tbldn", "Save", icon = "save", style = "background-color : #616161 !important")
          ),
        ),
        uiOutput("dynamicTable"),
        depth = 3
      ),
      style = "margin : 1em; border : solid 0.5em #1976d2"
    )
  )
}

myDnButton <- function(outputId, label = "Download", type = "default", ...) {
  aTag <- tags$a(
    id = outputId, href = "", target = "_blank",
    class = paste0("btn btn-", type, " shiny-download-link"),
    download = NA, icon("download", label), ...
  )
}

emphasize <- function(path = FALSE, inputObj, dfobj, Cells, pres, genesets) {
  cat("emphasize\n")
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
      if (length(thisCellsIdx) == 0) {
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

    # names(ret) <- rlobj$location
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
  rownames(dfobj_new) <- NULL

  # hc <- hchart(
  # dfobj_new,
  # type = 'scatter',
  # hcaes(x = x, y = y, color = col)
  # ) %>%
  # hc_tooltip(FALSE) %>%
  # hc_exporting(enabled = TRUE)

  graphString <- "ggobj2 <- ggplot(dfobj_new, aes(x = x, y = y)) + geom_point(colour = colV)"

  if (path) { # add mean point to path
    dfobj_path <- data.frame()
    for (i in 1:length(cellValues)) {
      x <- mean(as.numeric(dfobj_new$x[cellValues[[i]]]))
      y <- mean(as.numeric(dfobj_new$y[cellValues[[i]]]))

      dfobj_new <- rbind(dfobj_new, c(x, y, "meanPoint"))
      colV <- c(colV, "#000000")

      dfobj_path <- rbind(dfobj_path, c(x, y))
    }
    colnames(dfobj_path) <- c("x", "y")

    # hc <- hc %>%
    # hc_add_series(dfobj_path, 'line', hcaes(x = x, y = y), showInLegend = FALSE) %>%
    # hc_plotOptions(
    # line = list(
    # marker = FALSE,
    # dashStyle = 'ShortDash',
    # lineColor = '#5f27cd',
    # lineWidth = 2
    # ))

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
  # return(hc)
}

sortItem <- function(label, tableName) {
  options(useFancyQuotes = FALSE)
  paste0(
    "$('#", tableName, "')",
    ".append(", "`<div class=", "'rank-list-item'", " draggable='true'",
    " style = 'transform: translateZ(0px);'>` + ", label, " + `</div>`)"
  )
}

solvedButton <- function(inputId, label, style = NULL, onClick = NULL, ...) {
  value <- restoreInput(id = inputId, default = NULL)
  tags$button(
    id = inputId, style = style, onClick = onClick,
    type = "button", class = "btn btn-default action-button",
    `data-val` = value, list(label),
    ...
  )
}


#'
#' @name CellEnrich
#' @title Pathway Enrichment Analysis / Visualize for Single Cell Data
#'
#' @param CountData CountData [dgCMatrix]
#' @param GroupInfo GroupInfo for each samples. [string]
#' @param genesets (optional), user geneset to analysis. [list]
#'
#' @return no return.
#'
#'
#' @importFrom DT dataTableOutput
#' @importFrom Matrix t
#'
#' @rawNamespace import(SingleCellExperiment, except = show)
#' @import Rtsne
#' @import dplyr
#' @rawNamespace import(shiny, except = dataTableOutput)
#' @import shinymaterial
#' @import ggplot2
#' @import uwot
#' @import htmltools
#' @import ggbiplot
#' @import magrittr
#' @import waiter
#' @rawNamespace import(shinyjs, except = runExample)
#' @import scales
#' @import sortable
#' @import scran
#'
#' @export

CellEnrich <- function(CountData, GroupInfo, genesets = NULL) {
  require(dplyr)
  require(plyr)
  require(shiny)
  options(useFancyQuotes = FALSE)

  server <- function(input, output, session) {


    buildbiplot <- function(relative = TRUE, biFont, biX, biY, genesets) {
      require(ggbiplot)

      Cells <- sort(unique(GroupInfo))

      tab <- matrix(0,nrow = length(genesets), ncol = length(Cells))

      for(i in 1:length(Cells)){
        thisCell <- Cells[i]
        thisCellIdx <- which(GroupInfo == thisCell)
        v <- rep(0,length(genesets))

        vs <- table(unlist(pres[thisCellIdx]))
        nvs <- as.numeric(names(vs))
        vs <- unname(vs)
        v[nvs] <- vs
        if(relative){
          tab[,i] = v/length(thisCellIdx)
        }
        else{
          tab[,i] = v
        }
      }
      rownames(tab) <- names(genesets)
      colnames(tab) <- Cells
      tab <- tab[-which(sapply(1:nrow(tab), function(i){sum(tab[i,])==0})),] # remove zero

      # select high in groups

      high <- c()
      for(i in 1:ncol(tab)){
        high <- c(high, names(tab[order(tab[,i], decreasing = TRUE)[1:5],i]))
      }
      high <- unique(high)
      tab <- tab[high,]

      model <- prcomp(tab, scale = TRUE)
      BiPlot <<- ggbiplot::ggbiplot(
        model,
        labels = rownames(tab),
        labels.size = biFont,
        scale = 1,
        var.scale = 1,
        obs.scale = 1) +
        xlim(c(-biX,biX)) +
        ylim(c(-biY,biY))

      return( BiPlot )

    }


    ### CODES


    # variable initialize

    dtobj <- dfobj <- pres <- pres2 <- ""
    CellPathwayDF <- ""
    gt <- Cells <- A <- ""
    CellScatter <- ""
    CellHistogram <- ""
    BiPlot <- ""

    observeEvent(input$StartCellEnrich, {
      pt <- proc.time()

      if (is.null(genesets)) {
        if (input$genesetOption == "User-Geneset") {
          shiny::showNotification("Geneset not given ...", type = "error", duration = 10)
          return(NULL)
        }
      }

      # ------ Hide Start Button

      shinyjs::hide("StartCellEnrich")

      # ------ Load Genesets

      if (is.null(genesets)) {
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

      else {
        shiny::showNotification("User Geneset will be used", type = "message", duration = 10)
      }

      genesets <<- genesets
      # ------ for test
      # q0 <- 0.05

      q0 <- input$qvalueCutoff

      # ------ Create new Waitress
      w <- Waitress$new(selector = NULL, theme = "overlay-radius")

      genes <- rownames(CountData)
      genesets <- GenesetFlush(genes, genesets)
      lgs <- getlgs(genesets)

      # ------ Genesetsize Flush

      genesets <- GenesetsizeFlush(genesets, lgs, input$minGenesetSize, input$maxGenesetSize)
      # ------ For Tests
      # genesets <- GenesetsizeFlush(genesets, lgs, 15, 500)

      # ------ Gene Flush
      remgenes <- GeneFlush(genes, genesets)
      CountData <- CountData[-remgenes, ]

      genesets <<- genesets

      # ------ Background genes
      A <<- getBackgroundGenes(genesets)

      # ------ Calculate t-SNE / U-MAP First
      # require(Matrix)

      # dfobj <- getTU(CountData, GroupInfo, 't-SNE')
      dfobj <- getTU(CountData, GroupInfo, input$plotOption)
      dfobj <<- dfobj

      cat("getTU Finished\n")

      # ------ Disable radio button
      shinyjs::runjs('$("form p label input").attr("disabled",true)')
      shinyjs::runjs("$('.shinymaterial-slider-minGenesetSize').attr('disabled',true)")
      shinyjs::runjs("$('.shinymaterial-slider-maxGenesetSize').attr('disabled',true)")
      shinyjs::runjs("$('.shinymaterial-slider-qvalueCutoff').attr('disabled',true)")

      cat("running gc\n")
      gc()

      # ------ Find Significant Genes with Fold Change
      if (input$FCoption != "GSVA") {
        # ------ need to build GSVA CASE

        # s <- findSigGenes(CountData, 'median', GroupInfo)
        s <- findSigGenes(CountData, input$FCoption, GroupInfo)
      }

      cat("s Finished\n")

      # ------ Find Significant Genes with findMarkers
      require(dplyr)
      s2 <- findSigGenesGroup(CountData, GroupInfo, q0, TopCutoff = 5)

      rc <- rownames(CountData)

      # ------ free memory to calculate biobj
      rm(CountData)

      # ------ marker l1
      markerl1 <- s2 %>% filter(Top < 10)
      markerl1$Group <- as.factor(markerl1$Group)

      shinyjs::runjs("$(.markerP).show()")

      output$markerL1 <- DT::renderDataTable(
        DT::datatable(markerl1,
          rownames = FALSE,
          filter = "top",
          options = list(
            autoWidth = TRUE,
            dom = "ltp",
            lengthChange = FALSE,
            columnDefs = list(list(className = "dt-center", targets = 0:3))
          ),
          selection = "none",
        )
      )

      # ------ Hypergeometric pvalue calculation
      lgs <- getlgs(genesets)
      lens <- length(s)
      lens100 <- round(lens / 100)

      biobj <- getbiobj(genes, genesets)

      pres <- list()

      presTab <- c()

      w$start()
      for (i in 1:lens) {
        if (i %% lens100 == 0) w$inc(1)
        prespv <- getHyperPvalue(rc[s[[i]]], genesets, A, lgs, q0, biobj)

        pres[[i]] <- which(p.adjust(prespv, "fdr") < q0)

        prespv[which(prespv < 1e-20)] <- 1e-20

        presTab <- cbind(presTab, -log10(prespv))
      }
      w$close()

      colnames(presTab) <- colnames(CountData)
      rownames(presTab) <- names(genesets)


      pres <<- pres

      # pres : which gene-sets are significant for each cells.

      # ------ CellPathwayDF

      CellPathwayDF <- buildCellPathwayDF(GroupInfo, pres, genesets)

      # pres2 : for each gene-sets, how many cells are significant that gene-sets.

      cat("pres2\n")

      pres2 <- sort(table(unlist(pres)), decreasing = T)
      names(pres2) <- names(genesets)[as.numeric(names(pres2))]
      pres2 <<- pres2

      # 2625*4
      PP <- pathwayPvalue(GroupInfo, pres, pres2, genesets) # qvalue cutoff removed
      OR <- getOddRatio(GroupInfo, pres, pres2, genesets, input$ORratio)
      # OR <- getOddRatio(GroupInfo, pres, pres2, genesets, 0.1)

      # QVCUTOFF <- 4
      CellPathwayDFP <- CellPathwayDF %>%
        inner_join(PP) %>%
        select(Cell, Geneset, Qvalue) %>%
        filter(Qvalue > 4)

      ggs <- unique(CellPathwayDFP %>% select(Geneset))[, 1]
      ces <- sort(unique(CellPathwayDFP %>% select(Cell))[, 1])

      nr <- length(ggs) # nrow
      nc <- length(ces) # ncol

      output$tbldn <- downloadHandler(
        filename = function() {
          "mytable.csv"
        },
        content = function(file) {
          # pathway - cell ? -log pvalue
          # outputFile = matrix(0,nr,nc)

          # rownames(outputFile) = ggs
          # colnames(outputFile) = ces

          # for(i in 1:length(ces)){
          # tf <- CellPathwayDFP %>% filter(Cell == ces[i]) %>% select(Geneset, Qvalue)
          # outputFile[(tf %>% select(Geneset))[,1],i] = (tf%>% select(Qvalue))[,1]
          # }
          # write.csv(outputFile, file)
          write.csv(presTab, file)
        }
      )


      CellPathwayDF <- CellPathwayDF %>%
        inner_join(OR)

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

        # s2 <- findSigGenesGroup(CountData, GroupInfo, q0, TopCutoff = 5)
        # find markers

        thisCellDEs <- s2 %>%
          filter(Group == thisCell) %>%
          select(genes)

        tcd <- thisCellDEs[, 1] # ThisCellDES
        tcp <- thisCellPathways[, 1] # ThisCellPathways
        tcp <- sapply(tcp, function(i) { # indexed
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
        additive <- data.frame(cbind(genes, Count, Group = thisCell))

        # ------ first add
        if (ncol(CellMarkers) == 0) {
          CellMarkers <- rbind(CellMarkers, data.frame(cbind(genes, Count, Group = thisCell), stringsAsFactors = FALSE))
        }
        else {
          if (ncol(CellMarkers) == ncol(additive)) {
            CellMarkers <- rbind(CellMarkers, data.frame(cbind(genes, Count, Group = thisCell), stringsAsFactors = FALSE))
          }
        }
      }

      if (nrow(CellMarkers)) {
        CellMarkers <- CellMarkers %>%
          inner_join(s2) %>%
          filter(Top < 10)

        CellMarkers$Group <- as.factor(CellMarkers$Group)
        CellMarkers$Count <- as.numeric(CellMarkers$Count)
        CellMarkers$FDR <- as.numeric(CellMarkers$FDR)

        output$markerL2 <- DT::renderDataTable(
          DT::datatable(CellMarkers,
            rownames = FALSE,
            filter = "top",
            options = list(
              autoWidth = TRUE,
              dom = "ltp",
              lengthChange = FALSE,
              columnDefs = list(list(className = "dt-center", targets = 0:4))
            ),
            selection = "none",
          )
        )
      }
      else {
        cat("CellMarker Not Available\n")
        output$markerL2 <- DT::renderDataTable(
          DT::datatable(
            s2,
            rownames = FALSE,
            filter = "top",
            options = list(
              autoWidth = TRUE,
              dom = "ltp",
              lengthChange = FALSE,
              columnDefs = list(list(className = "dt-center", targets = 0:4))
            ),
            selection = "none",
          )
        )
      }

      # group 별 significant pathways
      # group 별 DE Genes

      # is counted
      dtobj <<- buildDT(pres2)

      # ------ Color define
      colV <- getColv(GroupInfo)

      CellHistogram <<- getCellHistogram(GroupInfo, colV)

      output$CellBar <- renderHighchart(CellHistogram) # CELL HISTOGRAM

      CellScatter <<- getCellPlot(dfobj, Cells)

      output$CellPlot <- renderPlot(CellScatter)

      output$legenddn <- downloadHandler(
        filename = 'mylegend.png',
        content = function(file){
          png(file)
          plot(NULL ,xaxt='n',yaxt='n',bty='n',ylab='',xlab='', xlim=0:1, ylim=0:1)
          legend(
            'center',
             legend = c('Sugar maple', 'White ash', 'Black walnut','Red oak', 'Eastern hemlock'),
             pch = 16, pt.cex = 3, cex = 1.5, bty='n',
             col = c('orange', 'red', 'green', 'blue', 'purple')
          )
          dev.off()
        }
      )

      output$imgdn <- downloadHandler(
        filename = function() {
          "myfigure.png"
        },
        content = function(file) {
          ggsave(file, CellScatter, device = "png")
        }
      )

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
          "order = list(list(3,'desc'))", # odd ratio based
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
              thisCellData <- thisCellData %>% top_n(-1, wt = OddRatio)
              if (nrow(thisCellData) >= 1) {
                thisCellData <- thisCellData %>% top_n(1)
              }
            }
          }
          res <- c(res, paste0(thisCellData$Geneset, " @", thisCellData$Cell))
        }
      }

      output$legendTable <- DT::renderDataTable(
        buildLegend(res)
      )

      plotImg <- emphasize(FALSE, res, dfobj, Cells, pres, genesets)

      output$imgdn <- downloadHandler(
        filename = function() {
          "myfigure.png"
        },
        content = function(file) {
          ggsave(file, plotImg, device = "png")
        }
      )
      shinyjs::show('legenddn')

      output$legenddn <- downloadHandler(
        filename = 'mylegend.png',
        content = function(file){
          buildLegend(res, img = TRUE, name = file)
        }
      )

      output$CellPlot <- renderPlot(plotImg)

      # output$CellPlot <- renderPlot(emphasize(FALSE, res, dfobj, Cells, pres, genesets))
    })

    observeEvent(input$refreshBiplot, {
      if(input$refreshBiplot==0){return(NULL)}
      cat('refreshed\n')
      output$biplotdn <- downloadHandler(
        filename = function() {
          "mybiplot.png"
        },
        content = function(file) {
          ggsave(file, BiPlot, device = "png")
        }
      )
    })

    # draw significant colored images
    observeEvent(input$sigbtn, {
      if (input$sigbtn == 0) { # prevent default click state
        return(NULL)
      }

      shinyjs::show('legenddn')
      res <- c()

      for (i in 1:length(Cells)) {
        thisCell <- Cells[i]

        thisCellData <- CellPathwayDF %>% dplyr::filter(Cell == thisCell)
        if (nrow(thisCellData) >= 1) {
          thisCellData <- thisCellData %>% top_n(-1, wt = OddRatio)

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

      output$legendTable <- DT::renderDataTable(
        buildLegend(res)
      )

      plotImg <- emphasize(FALSE, res, dfobj, Cells, pres, genesets)

      output$imgdn <- downloadHandler(
        filename = function() {
          "myfigure.png"
        },
        content = function(file) {
          ggsave(file, plotImg, device = "png")
        }
      )

      output$legenddn <- downloadHandler(
        filename = 'mylegend.png',
        content = function(file){
          buildLegend(res, img = TRUE, name = file)
        }
      )

      output$CellPlot <- renderPlot(plotImg)
    })

    # draw group colored images
    observeEvent(input$colorbtn, {
      if (input$colorbtn == 0) { # prevent default click state
        return(NULL)
      }
      shinyjs::hide('legenddn')
      UniqueCol <- briterhex(scales::hue_pal()(length(Cells)))
      names(UniqueCol) <- Cells

      colV <- unname(UniqueCol[dfobj$col])

      colorImage <- ggplot(dfobj, aes(x = x, y = y)) +
        geom_point(colour = colV)

      # CellScatter <<- getCellPlot(dfobj, Cells)

      output$imgdn <- downloadHandler(
        filename = function() {
          "myfigure.png"
        },
        content = function(file) {
          ggsave(file, colorImage, device = "png")
        }
      )
      o <- data.frame(Pathway = "", Group = "")
      colnames(o) <- c("Pathway", "Group")

      # clear legendtable
      output$legendTable <- DT::renderDataTable(
        DT::datatable(
          o,
          rownames = FALSE,
          options = list(
            autoWidth = TRUE,
            dom = "ltp",
            lengthChange = FALSE,
            columnDefs = list(list(className = "dt-center", targets = 0:1))
          ),
          selection = "none"
        )
      )
      output$CellPlot <- shiny::renderPlot(colorImage)

      # output$CellPlot <- renderHighchart(CellScatter)
    })

    # Emphasize with order
    observeEvent(input$OrderEmphasize, {
      if (input$OrderEmphasize == 0) { # prevent default click state
        return(NULL)
      }

      shinyjs::show('legenddn')
      plotImg <- emphasize(TRUE, input$sortList, dfobj, Cells, pres, genesets)
      output$CellPlot <- renderPlot(plotImg)

      output$legendTable <- DT::renderDataTable(
        buildLegend(input$sortList)
      )

      output$imgdn <- downloadHandler(
        filename = function() {
          "myfigure.png"
        },
        content = function(file) {
          ggsave(file, plotImg, device = "png")
        }
      )

      output$legenddn <- downloadHandler(
        filename = 'mylegend.png',
        content = function(file){
          buildLegend(res, img = TRUE, name = file)
        }
      )
    })

    # Emphasize without order
    observeEvent(input$Emphasize, {
      if (input$Emphasize == 0) { # prevent default click state
        return(NULL)
      }
      shinyjs::show('legenddn')
      plotImg <- emphasize(FALSE, input$sortList, dfobj, Cells, pres, genesets)
      output$CellPlot <- renderPlot(plotImg)

      output$legendTable <- DT::renderDataTable(
        buildLegend(input$sortList)
      )

      output$imgdn <- downloadHandler(
        filename = function() {
          "myfigure.png"
        },
        content = function(file) {
          ggsave(file, plotImg, device = "png")
        }
      )

      output$legenddn <- downloadHandler(
        filename = 'mylegend.png',
        content = function(file){
          buildLegend(res, img = TRUE, name = file)
        }
      )
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

    observeEvent(input$relFreq,{
      if(input$relFreq == 0){
        return(NULL)
      }


      output$biPlot <- renderPlot(buildbiplot(relative = TRUE, input$biFont, input$biX, input$biY, genesets))

      output$biplotdn <- downloadHandler(
        filename = function() {
          "mybiplot.png"
        },
        content = function(file) {
          ggsave(file, device = "png")
        }
      )
    })

    observeEvent(input$absFreq,{
      if(input$absFreq == 0){
        return(NULL)
      }

      output$biPlot <- renderPlot(buildbiplot(relative = FALSE, input$biFont, input$biX, input$biY, genesets))

      output$biplotdn <- downloadHandler(
        filename = function() {
          "mybiplot.png"
        },
        content = function(file) {
          ggsave(file, device = "png")
        }
      )
    })
  }

  ui <- CellEnrichUI()

  shiny::shinyApp(ui, server, options = list(launch.browser = TRUE))
}


buildLegend <- function(sortList, img = FALSE, name = NULL) {

  colV <- getColv(GroupInfo)

  rlobj <- data.frame(stringsAsFactors = FALSE)

  for (i in 1:length(sortList)) {
    kk <- strsplit(sortList[[i]], " @")[[1]]
    Pathway <- kk[1]
    Group <- kk[2]
    rlobj <- rbind(rlobj, cbind(Pathway, Group))
  }

  colnames(rlobj) <- c("Pathway", "Group")

  if(img){
    png(name)
    plot(NULL ,xaxt='n',yaxt='n',bty='n',ylab='',xlab='', xlim = 0:1, ylim = 0:1)
    legend(
      'center',
      legend = rlobj$Pathway,
      pch = 15, pt.cex = 2, cex = 1.2, bty='n',
      col = colV[rlobj$Group]
    )
    dev.off()
    return()
  }


  rlobj$Pathway <- paste0(
    sapply(
      colV[rlobj$Group],
      function(i) {
        paste0('<div style="background: ', i, '; display: inline-block; width: 1em;height: 1em;"></div>')
      }
    ),
    " ", as.character(rlobj$Pathway)
  )
  rlobj$Group <- as.character(rlobj$Group)

  return(
    DT::datatable(
      rlobj,
      escape = FALSE,
      rownames = FALSE,
      options = list(
        autoWidth = TRUE,
        dom = "ltp",
        lengthChange = FALSE,
        columnDefs = list(list(className = "dt-center", targets = 0:1))
      ),
      selection = "none"
    )
  )
}
