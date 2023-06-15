## 23.06.08

if(!require(waiter)){
  install.packages('waiter') # install 'waiter' if not installed.
} 

if(!require(farver)){
  install.packages('farver') # install 'farver' if not installed.
  
}
library(farver)

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
  }, USE.NAMES = TRUE)  # changed to TRUE - Hai
  return(remgenes)
}

getBackgroundGenes <- function(genesets) {
  cat("getBackgroundGenes\n")
  length(unique(unlist(genesets)))
}

getTU <- function(CountData, GroupInfo, plotOption='UMAP') {
  cat("Mapping is started\n")
  library(Seurat)

  seu <- CreateSeuratObject(CountData)
  # store mitochondrial percentage in object meta data
  seu <- PercentageFeatureSet(seu, pattern = "^MT-", col.name = "percent.mt")
  # run sctransform
  seu <- SCTransform(seu, vars.to.regress = "percent.mt", verbose = FALSE)

  # Add cell type annotation to metadata
  seu <- AddMetaData(seu, GroupInfo, col.name = "cell_type")
  # celltype <- unique(GroupInfo)

  # Dimension reduction
  # These are now standard steps in the Seurat workflow for visualization and clustering
  seu <- RunPCA(seu, verbose = FALSE)

  # TSNE
  if (plotOption == "TSNE") {
    seu <- RunTSNE(seu, dims = 1:30)
  }

  # UMAP
  if (plotOption == "UMAP") {
    seu <- RunUMAP(seu, dims = 1:30, uwot.sgd = TRUE)
  }
  seu <- FindNeighbors(seu, verbose = FALSE, dims = 1:30)
  seu <- FindClusters(seu, algorithm = 3, random.seed = 7968, resolution = 0.5)
  return (seu)
}

# BiocManager::install('scMerge', force=TRUE)
library(scMerge)
gnm <- function(v) {
  out <- scMerge:::gammaNormMix(as.matrix(v), plot = FALSE )
  mat_prob <- matrix(out$probExpressed, nrow(v), ncol(v))
  mat_discretised <- 1 * (mat_prob > 0.5)
  return(mat_discretised)
}

findSigGenes <- function(v, method = "CellEnrich - median", Name) {
  if (!method %in% c("CellEnrich - median", "CellEnrich - mixture", "Fisher")) stop("wrong method")

  cat("findSigGenes started\n")
  rownames(v) <- colnames(v) <- NULL

  res <- list()

  if (method == "Fisher") {
    return(res)
  }
  if (method == "CellEnrich - mixture") {
    v <- gnm(v)
    for (i in 1:ncol(v)) {
      res[[i]] <- which(v[, i] > 0)
    }
  }
  else { # median
    cat("scaling\n")

    cat("define Lists\n")
    med2 <- function(v) {
      v <- v[which(v > 0)]
      return(median(v) / 2)
    }

    if (method == "CellEnrich - median") {
      for (i in 1:ncol(v)) {
        res[[i]] <- which(v[, i] > med2(v[, i]))
      }
    }
  }

  names(res) <- Name
  return(res)
}

findSigGenesGroup <- function(Count = NULL, ClustInfo = NULL, q0 = 0.1, TopCutoff = 5) {
  library(scran)
  if (is.null(Count)) stop("Count must given")
  if (is.null(ClustInfo)) stop("ClustInfo must given")

  GrpRes <- scran::findMarkers(x = as.matrix(Count), ClustInfo, test.type = "wilcox", direction = "up")
  Grp <- unique(ClustInfo)

  res <- data.frame(stringsAsFactors = FALSE)

  for (i in 1:length(Grp)) {
    G <- data.frame(
      genes = rownames(GrpRes[[i]]),
      Group = Grp[i],
      GrpRes[[i]],
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

  if(length(pres) == length(Cells)){ # FISHER
    for (i in 1:length(Cells)) {
      thisCell <- Cells[i]
      tt <- table(pres[[i]])

      if (length(tt)) {
        CellPathwayDF <- rbind(CellPathwayDF, cbind(thisCell, names(tt), unname(tt)))
      }
    }
  }
  else{
    for (i in 1:length(Cells)) {
      thisCell <- Cells[i]
      tt <- table(unlist(pres[which(thisCell == GroupInfo)]))

      if (nrow(tt)) {
        CellPathwayDF <- rbind(CellPathwayDF, cbind(thisCell, names(tt), unname(tt)))
      }
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

  if(length(pres) != length(Cells)){
    CellPathwayDF <- CellPathwayDF %>%
      dplyr::filter(Count > 1)
  }

  return(CellPathwayDF)
}

pathwayPvalue <- function(GroupInfo, pres, pres2, genesets) {
  cat("pathwayPvalue\n")
  res <- c()
  Cells <- unique(GroupInfo)
  total <- length(GroupInfo)

  if (length(pres) == length(Cells)) { # FISHER
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
        m <- pres2[thisPathway] # total white ball
        pv[j] <- 1 - phyper(q - 1, m, total - m, k)
      }
      names(pv) <- names(thisCellPathways)

      res <- rbind(res, cbind(thisCell, names(pv), unname(pv)))
    }
  }
  else {
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
    v <- as.vector(col2rgb(colors[i])) #* 1.3
    v <- sapply(v, function(i) {
      min(i, 255)
    })
    res[i] <- rgb(v[1], v[2], v[3], max = 255)
  }
  return(res)
}

getColv <- function(GroupInfo) {
  Cells <- unique(sort(GroupInfo))

  UniqueCol <- briterhex(scales::hue_pal(h = c(20, 350),
                                         c = 100, l = 65, h.start = 0,
                                         direction = 1)(length(Cells)))
  # print(UniqueCol)
  # "#FF998D" "#FFBD00" "#BFDD00" "#00F148" "#00FACE" "#00F0FF" "#7ECAFF" "#FF94FF" "#FF7EFD"
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
  # library(ggplot2)
  library(highcharter)
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
    hc_title(text = "Cell Group distribution") %>%
    hc_xAxis(categories = x) %>%
    hc_plotOptions(grouping = FALSE) %>%
    hc_add_series(data = y, colorByPoint = TRUE, showInLegend = FALSE, name = "Count") %>%
    hc_colors(colV) %>%
    hc_exporting(enabled = TRUE)
  return(hc)
}

getCellPlot <- function(dfobj, Cells) {
  cat("getCellPlot\n")
  library(ggplot2)
  # library(highcharter)

  colnames(dfobj) <- c("x", "y", "col")
  dfobj <<- dfobj


  UniqueCol <- briterhex(scales::hue_pal(h = c(20, 350), c = 100, l = 65, h.start = 0,
                                         direction = 1)(length(Cells)))
  names(UniqueCol) <- Cells

  colV <- unname(UniqueCol[dfobj$col])

  cat("\n")
  ap <- colV
  ap <- ifelse(ap=='#E5E5E5', 0.2, 1)

  p <- ggplot(dfobj, aes(x = x, y = y)) +    geom_point(colour = colV, alpha=ap)
  # p <- p + scale_color_manual(values = method_col)
  # p <- p + labs(x= NULL, y= "Percentage")
  # p <- p + theme(legend.position=0)+ coord_flip()
  # p <- p + theme(axis.title.y = element_text(size = 15, vjust= 0.5))
  # p <- p + theme(axis.text = element_text(size = 12))
  # p <- p + geom_hline(data=df, aes(yintercept=median),linetype="dashed", color='black')
  # p <- p + theme(axis.text.y = element_text(size=14, color = a))

  p <- p + theme(
    # axis.line = element_line(color = "black", size = 0.5, linetype = "solid"),
    # axis.title.x = element_text(size=16),
    # axis.title.y = element_blank(),
    # axis.text.x = element_text(size=14, colour = 'black', hjust = 1),
    # panel.grid.major = element_line(colour = "grey86"),
    # panel.grid.major.x = element_blank(),
    # panel.grid.minor.x = element_blank(),
    # panel.grid.minor.y = element_blank(),
    # panel.background = element_blank(),
    panel.background = element_rect(fill = 'white', colour = 'white'),
    panel.border = element_rect(fill = NA, colour = 'black', size=0.25 ),
    # panel.spacing.x = unit(0.5, "lines"),
    # panel.spacing.y = unit(1, "lines"),
    # strip.text = element_text(size=17, color="black"),
    # strip.background.x = element_rect(fill="#CDE8DF"),
    # strip.background.x = element_blank(),
    # strip.background.y = element_blank(),
    legend.position="none")

  return(p)
}

groupTable <- function(pres, genesets, dfobj, pres2) {
  cat("groupTable\n")
  # for pres2
  genesetIdx <- sapply(names(pres2), function(i) {
    v <- which(i == names(genesets))
    return(v[1])
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

# values <- c("Carbs" = "carbs", "Proteins" = "prots", "BMI" = "bmi")
default_genesets <- c(
  "Human-WikiPathway", # 2021
  "Human-KEGG", # KEGG 2021
  "Human-GO",
  "Human-GO-BP",
  "Human-GO-CC",
  "Human-GO-MF",
  "Mouse-WikiPathway", # 2019
  "Mouse-KEGG", # 2019
  "Mouse-GO",
  "Mouse-GO-BP",
  "Mouse-GO-CC",
  "Mouse-GO-MF")

CellEnrichUI <- function() {
  library(shinymaterial)
  library(highcharter)
  library(sortable)

  material_page(
    shinyjs::useShinyjs(),
    shinyFeedback::useShinyFeedback(feedback = TRUE, toastr = TRUE),
    # dynamic datatable full width

    tags$head(tags$style(type = "text/css", ".display.dataTable.no-footer{width : 100% !important;}
                                            ")),
    # waitress declare
    use_waitress(color = "#1976d2", percent_color = "#333333"),
    title = paste0(
      "CellEnrich ",
      "<a href = 'https://github.com/noobCoding/cellenrich' target = 'https://github.com/noobCoding/CellEnrich'> ", # github link
      "<i class='material-icons' style = 'font-size:1.3em;'>info</i> </a>" # icon tag
    ),
    nav_bar_fixed = FALSE,
    nav_bar_color = "blue darken-2",
    font_color = "#1976d2",
    include_fonts = TRUE,
    include_nav_bar = TRUE,
    include_icons = FALSE,

    # CellEnrich options
    material_row(
      material_column(
        material_card(
          title = shiny::tags$h4("Options"),
          divider = TRUE,
          style = "border : solid 0.5em #1976d2",
          material_row(
            material_column(
              material_card(
                radioButtons(
                  "FCoption",
                  label = HTML("<font color='black' size='5'>Methods</font>"),#"Methods",
                  choiceNames = list(
                    HTML("<font color='black'>CellEnrich - median</font>"),
                    HTML("<font color='black'>CellEnrich - mixture</font>"),
                    HTML("<font color='black'>Fisher</font>")
                  ),
                  choiceValues = c("CellEnrich - median", "CellEnrich - mixture", "Fisher"),
                  selected = "CellEnrich - median",
                  # color = "#1976d2"
                )
              ),
              material_card(
                radioButtons(
                  "plotOption",
                  label = HTML("<font color='black' size='5'>Scatter Plot</font>"),#"Scatter Plot",
                  choiceNames = list(
                    tags$span(style = "color:black", "UMAP"),
                    tags$span(style = "color:black", "TSNE")
                  ),
                  choiceValues = c("UMAP", "TSNE"),
                  selected = "TSNE",
                  # color = "#1976d2"
                )
              ),
              width = 4
            ),
            material_column(
              material_card(
                material_number_box(
                  input_id = "minGenesetSize",
                  label = HTML("<font color='black' size='4'>Minimum Gene-set Size</font>"), #"Minimum Gene-set Size",
                  min_value = 10,
                  max_value = 30,
                  initial_value = 15,
                  step_size = 5
                ),
                material_number_box(
                  input_id = "maxGenesetSize",
                  label = HTML("<font color='black' size='4'>Maximum Gene-set Size</font>"), #"Maximum Gene-set Size",
                  min_value = 250,
                  max_value = 750,
                  initial_value = 500,
                  step_size = 5
                ),
                material_number_box(
                  input_id = "ORratio",
                  label = HTML("<font color='black' size='4'>Pathway Frequency</font>"),
                  min_value = 0,
                  max_value = 0.5,
                  initial_value = 0.1,
                  step_size = 0.05
                ),
                material_number_box(
                  input_id = "qvalueCutoff",
                  label = HTML("<font color='black' size='4'>Q-value threshold</font>"),
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
                radioButtons(
                  "genesetOption",
                  label=HTML("<font color='black' size='5'>Genesets</font>"),#"Genesets",
                  choiceNames = list(
                    HTML("<font color='black'>Human-WikiPathway</font>"),
                    HTML("<font color='black'>Human-KEGG</font>"),
                    HTML("<font color='black'>Human-GO</font>"),
                    tags$span(style = "color:black", "Human-GO-BP"),
                    tags$span(style = "color:black", "Human-GO-CC"),
                    tags$span(style = "color:black", "Human-GO-MF"),
                    tags$span(style = "color:black", "Mouse-WikiPathway"),
                    tags$span(style = "color:black", "Mouse-KEGG"),
                    tags$span(style = "color:black", "Mouse-GO"),
                    tags$span(style = "color:black", "Mouse-GO-BP"),
                    tags$span(style = "color:black", "Mouse-GO-CC"),
                    tags$span(style = "color:black", "Mouse-GO-MF")
                  ),
                  choiceValues = default_genesets,
                  # color= "#1976d2",
                  selected= default_genesets[1] #"Human-WikiPathway"
                ),
                sidebarPanel(
                  tags$style("
                               .btn-file {
                                  background-color:#1976d2;
                                  border-color: none;
                               }
                               "),
                  HTML("<font color='black' size='3'>User's geneset input</font>"),
                  fileInput("other", "", placeholder = "RData format is required!"),
                  # fileInput("other", "User's geneset input", accept=".Rdata, .RData"),
                  # checkboxInput("header", "Header", TRUE),
                  actionButton("addgeneset", "Add geneset", style="color: #ffffff; background-color: #1976d2;")
                )
              ),
              width = 4
            )
          ),
          solvedButton(
            inputId = "StartCellEnrich",
            label = "Start",
            style = "margin-left:45%; background-color: #1976d2",
            onClick = 'console.log("CellEnrich");'
          ),
          depth = 3
        ),
        width = 6,
        offset = 3 # center half layout
      )
    ),

    # tSNE/UMAP plot - SlingShot - Monocle
    material_row(
      material_column(

        material_card(
          title = shiny::tags$h4("Scatter & Bar"), depth = 3,

          # Comparison - SlingShot
          material_row(
            material_column(width = 6,
                            plotOutput("Comparison", height = "480px") # cell distribution
            ),
            material_column(width = 6,
                            plotOutput("SlingShot", height = "480px")
            )
          ),
          # Cell - Bar
          material_row(
            material_column(width = 6,
                            plotOutput("CellPlot", height = "480px")
            ),
            material_column(width = 6,
                            highchartOutput("CellBar", height = "480px")
            )
          ),
          material_row(
            shinyjs::hidden(
              material_button("slingDraw", "SlingShot",  color = "blue darken-2", ),
              material_button("slingRedraw", "Redraw",  color = "blue darken-2", )
            ),
            # p('Coloring by'),
            material_button("colorbtn", "Cell groups", icon = shiny::icon("color_lens"), color = "blue darken-2"),
            material_button("freqbtn", "Frequency", icon = shiny::icon("grain"), color = "blue darken-2"),
            material_button("sigbtn", "Odds Ratio", icon = shiny::icon("grade"), color = "blue darken-2")

          ),
          material_row(
            shiny::downloadButton("imgdn", "Save Plot", icon = shiny::icon("save"), style = "background-color : #616161 !important"),
            shiny::downloadButton("sppcdn", "Significant pathways ", icon = shiny::icon("save"), style = "background-color : #616161 !important")

          ),
          material_row(
            material_card(
              title = "",
              DT::dataTableOutput("legendTable"),
              shiny::downloadButton("legenddn", "Save Legend", icon = shiny::icon("save"), style = "background-color : #616161 !important; display:none;")
            )
          ),
          material_row(
            material_card(
              title = shiny::tags$h4("Highlighted  pathways"), divider = TRUE,
              #tags$h4("To be recognized by application, Please move element's position"),
              my_ranklist<-rank_list(text = "", labels = "Please drag&drop to change pathway positions at least once!",
                                     input_id = "sortList", css_id = "mysortableCell"),
              material_row(
                # material_button("OrderEmphasize", "Emphasize with Order", icon = "timeline", color = "blue darken-2"),
                material_button("Emphasize", "Emphasize", icon = shiny::icon("bubble_chart"), color = "blue darken-2"),
                material_button("ClearList", "Clear List", icon = shiny::icon("clear_all"), color = "blue darken-2")
              ),
              material_row(
                shiny::downloadButton("tbldn", "Save Highlighted Pathways", icon = shiny::icon("save"), style = "background-color : #616161 !important")
              ),
            ),
            shiny::uiOutput("dynamicTable"), depth = 3
          )
        ), width = 12
      ), style = "margin : 1em; border : solid 0.5em #1976d2"
    ),

    # Biplot
    material_row(
      material_card(
        title = "",
        material_card(
          title = shiny::tags$h4("Biplot between pathways and cell groups"), divider = TRUE,
          material_row(
            material_column(
              plotOutput("biPlot", height = "700px"),
              width = 10
            ),
            # material_column(
            # DT::dataTableOutput('bitable'),
            # width = 2
            # ),
            material_column(
              material_row(
                numericInput("biCount", label = "Pathways uses in each group", value = 5, min = 1, max = 10, step = 1),
                numericInput("biFont", label = "Label Size", value = 5, min = 1, max = 10, step = 1),
                numericInput("gsFont", label = "Pathway Size", value = 5, min = 1, max = 10, step = 1),
                numericInput("biX", label = "Range of X-axis", value = 5, min = 1, max = 10, step = 1),
                numericInput("biY", label = "Range of Y-axis", value = 2, min = 1, max = 10, step = 1),
                numericInput("axlab", label = "Axes Title", value = 15, min = 10, max = 30, step = 1),
                numericInput("axtxt", label = "Axes Index", value = 13, min = 10, max = 30, step = 1),
                material_row(
                  material_button("freqbp", "Biplot with Frequency", color = "blue darken-2")
                ),
                material_row(
                  material_button("orbp", "Biplot with Odds Ratio", color = "blue darken-2")
                )
              ),
              material_row(
                shiny::downloadButton("biplotdn", "Save Biplot", icon = shiny::icon("save"), style = "background-color : #616161 !important")
              ),
              width = 2
            )
          ),
        ),
        depth = 3
      ),
      style = "margin : 1em; border : solid 0.5em #1976d2"
    ),

    # marker table
    material_row(
      material_card(
        title = shiny::tags$h4("Marker Genes"),
        material_row(
          material_column(
            material_card(
              title = "FindMarker function (scran)",
              DT::dataTableOutput("markerL1")
            ),
            width = 6
          ),
          material_column(
            material_card(
              title = "Frequently up-regulated in each group",
              DT::dataTableOutput("markerL2")
            ),
            width = 6
          )
        )
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

emphasize <- function(path = FALSE, inputObj, dfobj, Cells, pres, genesets, seu, presTab) {
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

      if (length(thisGeneset) > 1) {
        thisGeneset <- thisGeneset[1]
      }

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
  UniqueCol <- briterhex(scales::hue_pal(h = c(20, 350), c = 100, l = 65, h.start = 0,
                                         direction = 1)(length(Cells)))
  names(UniqueCol) <- Cells
  colV <- unname(UniqueCol[dfobj_new$col])

  colV[-unlist(cellValues, use.names = FALSE)] <- "#E5E5E5" # gray color

  dfobj_new$col <- colV
  ## Get -log10(p_value) of a cell based on which genesets are chosen

  cell_pval <- c()
  maxall <- max(presTab)
  minall <- min(presTab)
  # print(maxall)
  cellidx <- unlist(cellValues, use.names = FALSE)
  for (cell in cellidx){
    sgs <- names(genesets)[pres[[cell]]]
    for (i in 1:nrow(rlobj)) {
      thisGeneset <- which(sgs == rlobj[i, 1]) # index
      if (length(thisGeneset) > 1) {
        thisGeneset <- thisGeneset[1]
      }
      if (length(thisGeneset) > 0){
        gid <- which(names(genesets) == sgs[thisGeneset])
        tmp<- presTab[gid, cell]
        cell_pval<- c(cell_pval, tmp)
        break
      }
    }
  }

  # # scale to [0, 1]
  # pval_scaling <- function(x){(x-min(x))/(max(x)-min(x))}
  # cell_pval <- c(cell_pval, maxall)
  # cell_pval <- pval_scaling(cell_pval)
  # cell_pval <-ifelse(is.nan(cell_pval),0. ,cell_pval)
  # cell_pval <- cell_pval[1:(length(cell_pval)-1)]

  # absolute scale 0~0.001; 0.001~0.01; 0.01~0.05; 0.05~0.1; 0.1~1
  pval_scaling <- function(x){(x-min(x))/(max(x)-min(x))}
  cell_pval <- c(cell_pval, c(maxall, minall))
  cell_pval <- pval_scaling(cell_pval)
  cell_pval <-ifelse(is.nan(cell_pval),0. ,cell_pval)
  cell_pval <- cell_pval[1:(length(cell_pval)-2)]

  for (i in 1:length(cellidx)){
    # colV[cellidx[i]]<- col2hcl(colV[cellidx[i]],
    #                            c=min(100, 100*cell_pval[i]),
    #                            l=min(100, 100*cell_pval[i]))

    tmp <- decode_colour(colV[cellidx[i]], to="hcl")
    # tmp[2] <- min (100, tmp[2] * cell_pval[i])
    # tmp[3] <- min (100, tmp[3] + 100 * (1 - cell_pval[i]))
    tmp[3] <- min (100, 100 * cell_pval[i])
    colV[cellidx[i]] <- encode_colour(tmp, from = 'hcl')
  }
  ap <- dfobj_new$col
  ap <- ifelse(ap=='#E5E5E5', 0.2, 1)

  rownames(dfobj_new) <- NULL
  #, alpha=ap
  graphString <- "ggobj2 <-  ggplot(dfobj_new, aes(x = x, y = y)) + geom_point(colour = colV, alpha=ap) +
                            theme(panel.background = element_rect(fill = 'white', colour = 'white'),
                            panel.border = element_rect(colour = 'black', fill=NA, size=0.25)
                                )"

  if (path) { # add mean point to path
    dfobj_path <- data.frame()
    for (i in 1:length(cellValues)) {
      x <- mean(as.numeric(dfobj_new$x[cellValues[[i]]]))
      y <- mean(as.numeric(dfobj_new$y[cellValues[[i]]]))

      dfobj_new <- rbind(dfobj_new, c(x, y, "meanPoint"))
      # colV <- c(colV, "#aaaa00")
      colV <- c(colV, UniqueCol[i])

      dfobj_path <- rbind(dfobj_path, c(x, y))
    }
    colnames(dfobj_path) <- c("x", "y")

    newIdx <- (nrow(dfobj) + 1):nrow(dfobj_new)
    cellValues <- c(unname(unlist(cellValues)), newIdx)

    dfobj_new$x <- round(as.numeric(dfobj_new$x), 4)
    dfobj_new$y <- round(as.numeric(dfobj_new$y), 4)

    for (i in 1:(length(newIdx) - 1)) { # add curve
      newCurve <- paste(
        " + geom_curve( aes(x = ", "x[newIdx[", i,
        "]], y = y[newIdx[", i, "]], xend = x[newIdx[", i + 1,
        "]], yend = y[newIdx[", i + 1, ']]), size = 0.5, linetype = "longdash",',
        # "curvature = 0.1, colour = '#ffff00', ", 'arrow = arrow(length = unit(0.1,"inches")))'
        "curvature = 0.1, colour = colV[newIdx[", i, ']], arrow = arrow(length = unit(0.1,"inches")))'
      )
      graphString <- paste(graphString, newCurve, sep = "")
    }
  }
  eval(parse(text = graphString))

  return(ggobj2)
}

emphasizePathway <-
  function(inputObj, dfobj, Cells, pres, genesets, seu, presTab) {
    cat("emphasize pathways\n")
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

        if (length(thisGeneset) > 1) {
          thisGeneset <- thisGeneset[1]
        }

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
    UniqueCol <- briterhex(scales::hue_pal(h = c(20, 350), c = 100, l = 65, h.start = 0,
                                           direction = 1)(length(Cells)))
    names(UniqueCol) <- Cells
    colV <- unname(UniqueCol[dfobj_new$col])
    colV[-unlist(cellValues, use.names = FALSE)] <- "#E5E5E5" # gray color
    dfobj_new$col <- colV
    ## Get -log10(p_value) of a cell based on which genesets are chosen

    cell_pval <- c()
    cellidx <- unlist(cellValues, use.names = FALSE)

    for (cell in cellidx){
      sgs <- names(genesets)[pres[[cell]]]
      theseGenesets <- which(sgs %in% rlobj[,1]) # index

      for (i in 1:nrow(rlobj)) {
        thisGeneset <- which(sgs == rlobj[i, 1]) # index
        if (length(thisGeneset) > 1) thisGeneset <- thisGeneset[1]
        if (length(thisGeneset) > 0){
          gid <- which(names(genesets) == sgs[thisGeneset])
          tmp<- presTab[gid, cell]
          cell_pval<- c(cell_pval, tmp)
          break
        }
      }

      if (length(theseGenesets) > 0){
        gsid <- which(names(genesets) %in% sgs[theseGenesets])

        # find the best value
        tmp <- presTab[gsid, cell]
        # cell_pval <- c(cell_pval, max(tmp))

        # find the best geneset
        bestgsid <- gsid[which(tmp == max(tmp))]
        bestgs <- which(rlobj[,1] %in% names(genesets)[bestgsid])
        bestgs <- bestgs[1]
        # obtain the corresponding color
        colV[cell] <- unname(UniqueCol[rlobj[bestgs, 2]])
      }
    }

    # absolute scale 0~0.001; 0.001~0.01; 0.01~0.05; 0.05~0.1; 0.1~1
    pval_scaling <- function(x){(x-min(x))/(max(x)-min(x))}
    cell_pval <- c(cell_pval, c(max(presTab), min(presTab)))
    cell_pval <- pval_scaling(cell_pval)
    cell_pval <-ifelse(is.nan(cell_pval),0. ,cell_pval)
    cell_pval <- cell_pval[1:(length(cell_pval)-2)]

    for (i in 1:length(cellidx)){
      # colV[cellidx[i]]<- col2hcl(colV[cellidx[i]],
      #                            c=min(100, 100*cell_pval[i]),
      #                            l=min(100, 100*cell_pval[i]))

      tmp <- decode_colour(colV[cellidx[i]], to="hcl")
      # tmp[2] <- min (100, tmp[2] * cell_pval[i])
      tmp[3] <- min (100, 100 * cell_pval[i])
      colV[cellidx[i]] <- encode_colour(tmp, from = 'hcl')
    }
    ap <- dfobj_new$col
    ap <- ifelse(ap=='#E5E5E5', 0.2, 1)

    rownames(dfobj_new) <- NULL
    graphString <- "ggobj2 <-  ggplot(dfobj_new, aes(x = x, y = y)) + geom_point(color = colV, alpha=ap ) +
                             theme(panel.background = element_rect(fill = 'white', colour = 'white'),
                                  panel.border = element_rect(colour = 'black', fill=NA, size=0.25)
                                )"


    eval(parse(text = graphString))
    return(ggobj2)
  }

emphasizeSlingShot <- function(inputObj, dfobj, Cells, pres, genesets, seu, presTab) {
  cat("emphasize SlingShot \n")
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

      if (length(thisGeneset) > 1) {
        thisGeneset <- thisGeneset[1]
      }

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
    return(ret)
  }

  rlobj <- buildRlobj(inputObj) # split into name, location dataframe

  cellValues <- getCellValues(rlobj) # get cell index for each cell

  dfobj_new <- data.frame(dfobj, stringsAsFactors = FALSE)
  colnames(dfobj_new) <- c("x", "y", "col")

  # print(dfobj_new)
  dfobj_new_coord <- cbind(dfobj_new$x, dfobj_new$y)
  dfobj_new_coord <- dfobj_new_coord[unlist(cellValues, use.names = FALSE), 1:2]
  # dfobj_new_label <- dfobj_new$col[unlist(cellValues, use.names = FALSE)]
  dfobj_new_label <- seu$seurat_clusters[unlist(cellValues, use.names = FALSE)]

  # sds <- slingshot(dfobj_new_coord, clusterLabels = dfobj_new_label)
  sds <- getLineages(dfobj_new_coord, clusterLabels = dfobj_new_label)
  # print (sds@lineages)
  # print (length(sds@lineages))

  # define ggobj2 element
  UniqueCol <- briterhex(scales::hue_pal(h = c(20, 350), c = 100, l = 65, h.start = 0,
                                         direction = 1)(length(Cells)))
  names(UniqueCol) <- Cells
  colV <- unname(UniqueCol[dfobj_new$col])
  colV[-unlist(cellValues, use.names = FALSE)] <- "#E5E5EEE5" # gray color

  dfobj_new$col <- colV
  ## Get -log10(p_value) of a cell based on which genesets are chosen

  rownames(dfobj_new) <- NULL
  # graphString <- "ggobj2 <- ggplot(dfobj_new, aes(x = x, y = y)) + geom_point(colour = colV)"

  now_obj <- cbind(dfobj_new$x, dfobj_new$y)
  colnames(now_obj) <- c("x", "y")

  myslingshot <-
    list (plot(now_obj, col = colV, pch = 16, cex = 0.8),
          lines(SlingshotDataSet(sds), lwd = 1.5, type = 'lineages', col='red4'))

  return(myslingshot)
}

sortItem <- function(label, tableName) {
  options(useFancyQuotes = FALSE)
  paste0(
    "$('#", tableName, "')",
    ".append(", "`<div class=", "'rank-list-item'", " draggable='false'",
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
#' @importFrom scMerge gammaNormMix
#'
#' @rawNamespace import(SingleCellExperiment, except = show)
#' @import plyr
#' @import dplyr
#' @rawNamespace import(shiny, except = dataTableOutput)
#' @import shiny
#' @import shinymaterial
#' @import ggplot2
#' @import uwot
#' @import htmltools
#' @import ggbiplot
#' @import magrittr
#' @import waiter
#' @import faver
#' @rawNamespace import(shinyjs, except = runExample)
#' @import scales
#' @import sortable
#' @import scran
#' @import ggrepel
#' @import shinyFeedback
#'
#' @export

CellEnrich <- function(CountData, GroupInfo, genesets = NULL) {
  library(dplyr)
  library(shiny)

  if(!require(ggbiplot)){
    remotes::install_github('vqv/ggbiplot')
  }

  library(ggbiplot)
  library(ggrepel)
  options(useFancyQuotes = FALSE)

  server <- function(input, output, session) {

    buildbiplot <- function(biFont, biX, biY, genesets, TOPN = 5, oddratio = FALSE, gsFont=5, axtxt=13, axlab=15) {
      Cells <- sort(unique(GroupInfo))
      # pres : which gene-sets are significant for each cells.
      # pres2 : for each gene-sets, how many cells are significant that gene-sets.

      if (oddratio) { ## ODDRATIO
        total <- length(GroupInfo)

        dat <- OR %>%
          group_by(Cell) %>%
          arrange(Cell) %>%
          top_n(TOPN)

        gs <- unique(dat$Geneset)

        tab <- matrix(0, nrow = length(gs), ncol = length(Cells))
        rownames(tab) <- gs
        colnames(tab) <- Cells

        gs <- sapply(gs, function(i) {
          which(names(genesets) == i)
        })

        for (i in 1:length(Cells)) {
          thisCell <- Cells[i]
          thisCellIdx <- which(GroupInfo == thisCell)

          tab[, i] <- round(unname(
            sapply(1:length(gs), function(k) {
              k <- gs[k]
              B <- table(unlist(pres[thisCellIdx]))[as.character(unname(k))] # 특정 Cell에서 유의한 회수
              if (is.na(B)) {
                return(0)
              }
              A <- pres2[names(k)]

              if (is.na(A)) {
                return(0)
              }
              N <- total
              K <- length(thisCellIdx)

              return((B / K) / (A / N))
            })
          ), 4)
          # Cell, Geneset, OR
        }
      }
      else { # FREQUENCY

        tab <- matrix(0, nrow = length(genesets), ncol = length(Cells))

        for (i in 1:length(Cells)) {
          thisCell <- Cells[i]
          thisCellIdx <- which(GroupInfo == thisCell)
          v <- rep(0, length(genesets))

          vs <- table(unlist(pres[thisCellIdx]))
          nvs <- as.numeric(names(vs))
          vs <- unname(vs)
          v[nvs] <- vs

          tab[, i] <- v / length(thisCellIdx)
        }
        rownames(tab) <- names(genesets)
        colnames(tab) <- Cells
        tab <- tab[-which(sapply(1:nrow(tab), function(i) {
          sum(tab[i, ]) == 0
        })), ] # remove zero

        # select high in groups
        high <- c()
        if (TOPN == 1){
          for (i in 1:ncol(tab)) {
            high <- c(high, rownames(tab[order(-tab[,i]),])[1])
          }
        }
        else{
          for (i in 1:ncol(tab)) {
            high <- c(high, rownames(tab[order(-tab[, i]),])[1:TOPN])
          }
        }

        high <- unique(high)
        tab <- tab[high, ]
      }
      labels <- rownames(tab)

      model <- prcomp(tab, scale = TRUE)
      library(ggplot2)
      BiPlot <<-
        ggbiplot(
          model,
          labels = NULL,
          labels.size = 0,
          varname.size = biFont,
          ellipse = TRUE,
          scale = 1
          # var.scale = 1
          # obs.scale = 1

        ) +
        xlim(c(-biX, biX)) +
        ylim(c(-biY, biY)) +
        geom_text_repel(
          label = labels,
          size = gsFont,
          box.padding = 1,
          point.padding = 1,
          max.overlaps = 30
          )+
        labs(y = 'PC2', x = 'PC1') +
        theme(axis.line = element_blank(), #element_line(color = "black", size = 0.5, linetype = "solid"),
                 axis.title.x = element_text(size=axlab),
                 axis.title.y = element_text(size=axlab), #element_blank(),
                 axis.text.x = element_text(size=axtxt, colour = 'black', hjust = 1),
                 axis.text.y = element_text(size=axtxt, colour = 'black', hjust = 1),
                 panel.grid.major = element_line(colour = "grey86"),
                 # panel.grid.major.x = element_blank(),
                 # panel.grid.minor.x = element_blank(),
                 # panel.grid.minor.y = element_blank(),
                 panel.background = element_blank(),
                 panel.border = element_blank(),
                 panel.spacing.x = unit(0.5, "lines"),
                 panel.spacing.y = unit(1, "lines"),
                 strip.text = element_text(size=17, color="black"),
                 strip.background.x = element_rect(fill="#CDE8DF"),
                 # strip.background.x = element_blank(),
                 # strip.background.y = element_blank(),
                 legend.position="none")

      return(BiPlot)
    }

    buildSlingShot <- function(seu, plotOption) {
      # BiocManager::install('slingshot')
      library(slingshot)
      # TSNE
      if (plotOption == "TSNE") {
        sds <- getLineages(Embeddings(seu, "tsne"), clusterLabels = seu$seurat_clusters)
      }

      # UMAP
      if (plotOption == "UMAP") {
        sds <- getLineages(Embeddings(seu, "umap"), clusterLabels = seu$seurat_clusters)
      }

      cell_pal <- function(cell_vars, pal_fun,...) {
        if (is.numeric(cell_vars)) {
          pal <- pal_fun(100, ...)
          return(pal[cut(cell_vars, breaks = 100)])
        } else {
          categories <- sort(unique(cell_vars))
          pal <- setNames(pal_fun(length(categories), ...), categories)
          return(pal[cell_vars])
        }
      }
      cell_colors <- cell_pal(seu$cell_type, hue_pal())
      # cell_colors_clust <- cell_pal(seu$seurat_clusters,  brewer_pal("qual", "Set1"))

      myslingshot <<-
        list (plot(sds@elementMetadata@listData$reducedDim, col = cell_colors, pch = 16, cex = .8),
              lines(SlingshotDataSet(sds), lwd = 1.5, type = 'lineages', col='red4'))
      return(myslingshot)
    }

    usergs <- ""

    observeEvent(input$addgeneset, {
      other_geneset <- input$other
      # print(other_geneset)
      ext <- tools::file_ext(other_geneset$datapath)
      # print(ext)

      req(other_geneset)
      # validate(need(!ext %in% c("csv", 'txt', 'RData', 'Rdata'), "Please upload a CSV, TXT, Rdata or RData file!"))
      inFile <- other_geneset$datapath

      otherVal <- other_geneset$name
      usergs <<- inFile
      updatedValues <- c(default_genesets, otherVal)

      # for(i in 1:(length(updatedValues)-1) ){
      #   # list(
      #   #   HTML("<font color='red'>Normal</font>"),
      #   #   tags$span(style = "color:red", "Uniform"),
      #   #   "Log-normal", "Exponential"
      #   # )
      #   updatedValues[i] <- paste0("tags$span(style = 'color:green'", updatedValues[i] , ")")
      # }
      # updatedValues[length(updatedValues)] <- paste0("tags$span(style = 'color:red'", updatedValues[length(updatedValues)] , ")")
      updateRadioButtons(session, "genesetOption", choices = updatedValues, selected= otherVal)
      print(usergs)
      cat("add genesets done!")
    })

    ### CODES

    # variable initialize

    dtobj <- dfobj <- pres <- pres2 <- presTab <- ""
    CellPathwayDF <- ""
    gt <- Cells <- A <- ""
    CellScatter <- ""
    CellHistogram <- ""
    BiPlot <- OR <- ""

    observeEvent(input$StartCellEnrich, {
      pt <- proc.time()

      if (input$FCoption == "Fisher") {
        shinyjs::runjs('$("#colorbtn").attr("disabled",true)')
        shinyjs::runjs('$("#freqbtn").attr("disabled",true)')
        shinyjs::runjs('$("#sigbtn").attr("disabled",true)')
        shinyjs::runjs('$("#Emphasize").attr("disabled",true)')
        shinyjs::runjs('$("#freqbp").attr("disabled",true)')
        shinyjs::runjs('$("#orbp").attr("disabled",true)')

        shinyFeedback::showToast(
          type = "error",
          message = "Emphasize / Biplot will not be available with Fisher",
          .options = list(timeOut = 20000)
        )
      }

      # ------ Hide Start Button
      shinyjs::hide("StartCellEnrich")

      #Disable Emphasize 1st
      shinyjs::runjs('$("#Emphasize").attr("disabled",true)')

      # ------ Load Genesets
      # Check that data object exists and is data frame.
      if (is.null(genesets)) {
        if (input$genesetOption == "Human-WikiPathway") load("WikiPathways_2021_Human.RData")
        if (input$genesetOption == "Human-KEGG") load("KEGG_2021_Human.RData")
        if (input$genesetOption == "Human-GO") load("humanGO.RData")
        if (input$genesetOption == "Human-GO-BP") load("humanGOBP.RData")
        if (input$genesetOption == "Human-GO-CC") load("humanGOCC.RData")
        if (input$genesetOption == "Human-GO-MF") load("humanGOMF.RData")
        if (input$genesetOption == "Mouse-WikiPathway") load("WikiPathways_2019_Mouse.RData")
        if (input$genesetOption == "Mouse-KEGG") load("KEGG_2019_Mouse.RData")
        if (input$genesetOption == "Mouse-GO") load("mouseGO.RData")
        if (input$genesetOption == "Mouse-GO-BP") load("mouseGOBP.RData")
        if (input$genesetOption == "Mouse-GO-CC") load("mouseGOCC.RData")
        if (input$genesetOption == "Mouse-GO-MF") load("mouseGOMF.RData")
      }

      read_Input_geneset<-function(gs){
        genesets <- ''
        if (endsWith(gs, '.xlsx')){
          genesets <- read_excel(gs)
        } else if (endsWith(gs, '.csv')){
          genesets <- read_csv(gs)
        } else
          if (endsWith(gs, '.RData') || endsWith(gs, 'Rdata')){
            load(gs)
          }
        return(genesets)
      }

      if (is.null(genesets)) {
        # load(usergs)
        # load("hmgobp-usergs.RData")
        genesets <- read_Input_geneset(usergs)

        if (is.null(genesets)){
          shiny::showNotification("Geneset file is invalid!", type = "error", duration = 10)
          return(NULL)
        }
      }

      if (is.null(genesets)) {
        shiny::showNotification("Geneset file is missing!", type = "error", duration = 10)
        return(NULL)
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
      #CountData <- CountData[-remgenes, ]
      CountData <- CountData[!(rownames(CountData) %in% names(remgenes)),]
      genes <- genes[!genes %in% names(remgenes)]
      rm(remgenes)

      genesets <<- genesets

      # ------ Background genes
      A <<- getBackgroundGenes(genesets)

      # ------ Calculate TSNE / UMAP First
      # library(Matrix)

      seu <- getTU(CountData, GroupInfo, input$plotOption)
      seu <<- seu

      # TSNE
      if (input$plotOption == "TSNE") {
        dfobj <- data.frame(Embeddings(seu, 'tsne'), col = GroupInfo, stringsAsFactors = FALSE)
      }
      # UMAP
      if (input$plotOption == "UMAP") {
        dfobj <- data.frame(Embeddings(seu, 'umap'), col = GroupInfo, stringsAsFactors = FALSE)
      }
      colnames(dfobj) <- c("x", "y", "col")
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
      CountData <- NormalizeData(CountData, normalization.method = 'LogNormalize', scale.factor = 1e6)

      if (input$FCoption != "GSVA") {
        # ------ need to build GSVA CASE

        # s <- findSigGenes(CountData, 'median', GroupInfo)
        s <- findSigGenes(CountData, input$FCoption, GroupInfo)
      }
      cat("s Finished\n")

      # ------ Find Significant Genes with findMarkers
      s2 <- findSigGenesGroup(CountData, GroupInfo, q0, TopCutoff = 5)
      rc <- rownames(CountData)

      # ------ free memory to calculate biobj
      rm(CountData)

      # ------ marker l1
      markerl1 <- s2 %>% filter(Top <= 10)
      markerl1$Group <- as.factor(markerl1$Group)
      colnames(markerl1)[4] <- "FDR < 0.01"
      shinyjs::runjs("$(.markerP).show()")

      output$markerL1 <- DT::renderDataTable(
        DT::datatable(markerl1,
                      rownames = FALSE,
                      filter = "top",
                      options = list(
                        autoWidth = TRUE,
                        dom = "ltp",
                        # lengthChange = FALSE,
                        columnDefs = list(list(className = "dt-center", targets = 0:3))
                      ),
                      selection = "none",
        )
      )
      tmp_df <- data.frame()

      cat("biobj \n")
      biobj <- getbiobj(genes, genesets)

      if (length(s) == 0) {
        tmp_cells <- unique(s2$Group)
        tmp_pres <- matrix(0, length(genesets), length(tmp_cells))
        tG <- table(GroupInfo)
        # define s
        for (i in 1:length(tmp_cells)) {
          tmp_genes <- s2 %>%
            filter(Group == tmp_cells[i]) %>%
            dplyr::select(genes) %>%
            unlist() %>%
            unname()

          tmp_pv <- getHyperPvalue(tmp_genes, genesets, A, lgs, q0, biobj)
          sigidx <- which(p.adjust(tmp_pv, "fdr") <= q0)
          tmp_pres[sigidx, i] <- unname(tG[tmp_cells[i]])
          tmp_pv[which(tmp_pv < 1e-12)] <- 1e-12

          tmp_genes <- sapply(tmp_genes, function(i) {
            which(rc == i)
          }, USE.NAMES = FALSE)

          s[[i]] <- tmp_genes

          names(s)[i] <- tmp_cells[i]
        }

        colnames(tmp_pres) <- tmp_cells
        rownames(tmp_pres) <- names(genesets)

        # fisher ODD RATIO
        ors <- rowSums(tmp_pres) / sum(tG)
        ors[which(ors != 0)] <- 1 / ors[which(ors != 0)]

        for (i in 1:ncol(tmp_pres)) {
          pathways <- names(which(tmp_pres[, i] != 0))
          tmp_df <- rbind(
            tmp_df,
            data.frame(
              cell = colnames(tmp_pres)[i],
              pathway = pathways,
              oddratio = unname(ors[pathways])
            )
          )
        }
      }

      # ------ Hypergeometric pvalue calculation
      lgs <- getlgs(genesets)
      lens <- length(s)
      lens100 <- round(lens / 100)

      pres <- list()

      cat("pres declare\n")
      presTab <- c()

      if (length(s) >= 100) {
        w$start()
        for (i in 1:lens) {
          if (i %% lens100 == 0) w$inc(1)
          prespv <- getHyperPvalue(rc[s[[i]]], genesets, A, lgs, q0, biobj)

          pres[[i]] <- which(p.adjust(prespv, "fdr") <= q0)

          prespv[which(prespv < 1e-12)] <- 1e-12

          presTab <- cbind(presTab, -log10(prespv))
        }
        w$close()
        colnames(presTab) <- colnames(CountData)
      }
      else {
        for (i in 1:lens) {
          prespv <- getHyperPvalue(rc[s[[i]]], genesets, A, lgs, q0, biobj)

          pres[[i]] <- which(p.adjust(prespv, "fdr") <= q0)

          prespv[which(prespv < 1e-12)] <- 1e-12

          presTab <- cbind(presTab, -log10(prespv))
        }
        colnames(presTab) <- names(s)
      }

      cat("pres defined\n")
      rownames(presTab) <- names(genesets)
      presTab <<- presTab
      # write.csv(presTab, 'presTab.csv')

      pres <<- pres
      # [[686]]
      # [1]  11  62 102 118 161 197 215 216 218 225 228 229 230
      #
      # [[687]]
      # [1]  62  69 140 216 229
      # pres : which gene-sets are significant for each cells.

      # ------ CellPathwayDF

      CellPathwayDF <- buildCellPathwayDF(GroupInfo, pres, genesets)

      # group / pathway count

      # pres2 : for each gene-sets, how many cells are significant that gene-sets.

      cat("pres2\n")

      pres2 <- sort(table(unlist(pres)), decreasing = T)
      if (length(s) != 0) {
        names(pres2) <- names(genesets)[as.numeric(names(pres2))]
      }
      pres2 <<- pres2
      # print(pres2)

      # 2625*4
      PP <- pathwayPvalue(GroupInfo, pres, pres2, genesets) # qvalue cutoff removed

      if (nrow(tmp_df) > 0) {
        OR <- tmp_df
        colnames(OR) <- c("Cell", "Geneset", "OddRatio")
        OR <- OR %>% filter(OddRatio > 1)
        OR <<- OR
      }
      else {
        # OR <- getOddRatio(GroupInfo, pres, pres2, genesets, 0.1)
        OR <<- getOddRatio(GroupInfo, pres, pres2, genesets, input$ORratio)
      }

      # OR -> # Group / PATHWAY / ODDRATIO

      # QVCUTOFF <- 4

      # group / pathway /
      CellPathwayDFP <- CellPathwayDF %>%
        inner_join(PP) %>%
        dplyr::select(Cell, Geneset, Qvalue) %>%
        filter(Qvalue > 4)
      write.csv(CellPathwayDFP, 'cellpathwaydf_p.csv')

      ggs <- unique(CellPathwayDFP %>% dplyr::select(Geneset))[, 1]
      ces <- sort(unique(CellPathwayDFP %>% dplyr::select(Cell))[, 1])

      nr <- length(ggs) # nrow
      nc <- length(ces) # ncol

      output$tbldn <- downloadHandler(
        filename = "myhighlight.csv",
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
          # write.csv(presTab, file)

          # name <- location <- c()
          # for (i in 1:length(input$sortList)) {
          #   kk <- strsplit(input$sortList[[i]], " @")[[1]]
          #   name <- c(name, kk[1])
          #   location <- c(location, kk[2])
          # }
          # newdf <- do.call(rbind, Map(data.frame, Pathway=name, Cell=location))
          # write.csv(newdf, file, row.names = FALSE)
          write.csv(CellPathwayDF, file, row.names = FALSE)
        }
      )

      # presTab
      output$sppcdn <- downloadHandler(
        filename = "sig_pw_cell.csv",
        content = function(file) {
          write.csv(presTab, file, row.names = TRUE)
        }
      )

      CellPathwayDF <- CellPathwayDF %>%
        inner_join(OR)

      CellPathwayDF <<- CellPathwayDF
      write.csv(CellPathwayDF, "cellpathwaydf_innerjoin.csv")


      # l2
      CellMarkers <- data.frame()

      Cells <- sort(unique(GroupInfo))
      Cells <<- Cells

      for (i in 1:length(Cells)) {
        thisCell <- Cells[i]
        thisCellPathways <- CellPathwayDF %>%
          filter(Cell == thisCell) %>%
          dplyr::select(Geneset)

        thisCellDEs <- s2 %>%
          filter(Group == thisCell) %>%
          dplyr::select(genes)

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
        Count <- as.numeric(unname(tcp))

        additive <- data.frame(cbind(genes, Count, Group = thisCell))
        additive$Count <- as.numeric(additive$Count)
        additive <- additive %>% arrange(dplyr::desc(Count))
        additive <- additive[1:min(nrow(additive), 20), ]

        # ------ first add
        if (ncol(CellMarkers) == 0) {
          CellMarkers <- additive
        }
        else {
          if (ncol(CellMarkers) == ncol(additive)) {
            CellMarkers <- rbind(CellMarkers, additive)
          }
        }
      }

      if (nrow(CellMarkers)) {

        # CellMarkers <- Genes Count Group

        # CellMarkers <- CellMarkers %>%
        # inner_join(s2) %>%
        # filter(Top < 10)

        CellMarkers$Group <- as.factor(CellMarkers$Group)
        CellMarkers$Count <- as.numeric(CellMarkers$Count)
        # CellMarkers$FDR <- as.numeric(CellMarkers$FDR)

        output$markerL2 <- DT::renderDataTable(
          DT::datatable(CellMarkers,
                        rownames = FALSE,
                        filter = "top",
                        options = list(
                          autoWidth = TRUE,
                          dom = "ltp",
                          # lengthChange = FALSE,
                          columnDefs = list(list(className = "dt-center", targets = 0:2))
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
      # group 별DE Genes

      # is counted
      dtobj <<- buildDT(pres2)

      # ------ Color define
      colV <- getColv(GroupInfo)

      CellHistogram <<- getCellHistogram(GroupInfo, colV)

      output$CellBar <- renderHighchart(CellHistogram) # CELL HISTOGRAM

      CellScatter <<- getCellPlot(dfobj, Cells)

      output$CellPlot <- renderPlot(CellScatter)

      output$legenddn <- downloadHandler(
        filename = "mylegend.png",
        content = function(file) {
          png(file)
          plot(NULL, xaxt = "n", yaxt = "n", bty = "n", ylab = "", xlab = "", xlim = 0:1, ylim = 0:1)
          legend(
            "center",
            legend = c("Sugar maple", "White ash", "Black walnut", "Red oak", "Eastern hemlock"),
            pch = 16, pt.cex = 3, cex = 1.5, bty = "n",
            col = c("orange", "red", "green", "blue", "purple")
          )
          dev.off()
        }
      )

      output$imgdn <- downloadHandler(
        filename = "myfigure.png",
        content = function(file) {
          ggsave(file, CellScatter, device = "png")
        }
      )

      gt <<- groupTable(pres, genesets, dfobj, pres2)

      # generate dynamic table

      cat("dynT\n")

      output$dynamicTable <- renderUI({
        numTabs <- length(Cells)
        CardColors <- briterhex(scales::hue_pal(h = c(20, 350), c = 100, l = 65, h.start = 0,
                                                direction = 1)(length(Cells)))
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
                          "; $('#toSortButton", i, "').attr('disabled', true); ",
                          'Shiny.onInputChange("input", this.input); '
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
      shinyjs::click('ClearList')
      shinyjs::click('freqbp')
      shinyjs::click('slingDraw')
      shiny::showNotification("Please   DRAG&DROP   the pathway positions ONCE to activate 'EMPHASIZE' button!",
                              type = "error", duration = NULL )

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
              thisCellData <- thisCellData %>% top_n(1, wt = OddRatio)
              if (nrow(thisCellData) >= 1) {
                thisCellData <- thisCellData %>% top_n(1)
              }
            }
          }
          res <- c(res, paste0(thisCellData$Geneset, " @", thisCellData$Cell))
        }
      }

      output$legendTable <- DT::renderDataTable(
        buildGradientLegend(res, Cells = Cells)
      )

      plotImg <- emphasize(FALSE, res, dfobj, Cells, pres, genesets, seu, presTab)
      output$SlingShot <- renderPlot(emphasizeSlingShot(res, dfobj, Cells, pres, genesets, seu, presTab))
      output$Comparison <- renderPlot(emphasizePathway(res, dfobj, Cells, pres, genesets, seu, presTab))

      output$imgdn <- downloadHandler(
        filename = function() {
          "myfigure.png"
        },
        content = function(file) {
          ggsave(file, plotImg, device = "png")
        }
      )
      shinyjs::show("legenddn")

      output$legenddn <- downloadHandler(
        filename = "mylegend.png",
        content = function(file) {
          buildGradientLegend(res, img = TRUE, name = file, Cells = Cells)
        }
      )

      output$CellPlot <- renderPlot(plotImg)
    })

    # draw significant colored images
    observeEvent(input$sigbtn, {
      if (input$sigbtn == 0) { # prevent default click state
        return(NULL)
      }

      shinyjs::show("legenddn")
      res <- c()
      # write.csv(CellPathwayDF, "cellpathwaydf_innerjoin_after.csv")
      for (i in 1:length(Cells)) {
        thisCell <- Cells[i]

        thisCellData <- CellPathwayDF %>% dplyr::filter(Cell == thisCell)
        # print (thisCellData)
        if (nrow(thisCellData) >= 1) {
          thisCellData <- thisCellData %>% top_n(1, wt = OddRatio)
          # print (thisCellData)
          if (nrow(thisCellData) >= 1) {
            thisCellData <- thisCellData %>% top_n(-1, wt = Length)
            if (nrow(thisCellData) >= 1) {
              thisCellData <- thisCellData %>% top_n(1, wt = Count)
              if (nrow(thisCellData) >= 1) {
                thisCellData <- thisCellData %>% top_n(1)
              }
            }
          }
          # print (thisCellData)
          res <- c(res, paste0(thisCellData$Geneset, " @", thisCellData$Cell))
          # print (res)
        }
      }

      output$legendTable <- DT::renderDataTable(
        buildGradientLegend(res, Cells = Cells)
      )

      plotImg <- emphasize(FALSE, res, dfobj, Cells, pres, genesets, seu, presTab)
      output$SlingShot <- renderPlot(emphasizeSlingShot(res, dfobj, Cells, pres, genesets, seu, presTab))
      output$Comparison <- renderPlot(emphasizePathway(res, dfobj, Cells, pres, genesets, seu, presTab))

      output$imgdn <- downloadHandler(
        filename = function() {
          "myfigure.png"
        },
        content = function(file) {
          ggsave(file, plotImg, device = "png")
        }
      )

      output$legenddn <- downloadHandler(
        filename = "mylegend.png",
        content = function(file) {
          buildLegend(res, img = TRUE, name = file, GroupInfo = GroupInfo)
        }
      )

      output$CellPlot <- renderPlot(plotImg)
    })

    # draw group colored images
    observeEvent(input$colorbtn, {
      if (input$colorbtn == 0) { # prevent default click state
        return(NULL)
      }
      shinyjs::hide("legenddn")
      UniqueCol <- briterhex(scales::hue_pal(h = c(20, 350), c = 100, l = 65, h.start = 0,
                                             direction = 1)(length(Cells)))
      names(UniqueCol) <- Cells

      colV <- unname(UniqueCol[dfobj$col])
      ap <- colV
      ap <- ifelse(ap=='#E5E5E5', 0.2, 1)

      colorImage <- ggplot(dfobj, aes(x = x, y = y)) +
        theme(panel.background = element_rect(fill = 'white', colour = 'white'),
              panel.border = element_rect(colour = 'black', fill=NA, size=0.25))+
        geom_point(colour = colV, alpha=ap)


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
      output$SlingShot <- renderPlot(buildSlingShot(seu, input$plotOption))
    })

    observeEvent(input$sortList, {
      shinyjs::runjs('$("#Emphasize").attr("disabled",false)')
    })

    observeEvent(input$Emphasize, {
      if (input$Emphasize == 0) { # prevent default click state
        return(NULL)
      }
      shinyjs::runjs('$("#Emphasize").attr("disabled",true)')

      # updateSelectInput(session, inputId = "sortList")

      res<-input$sortList

      shinyjs::show("legenddn")
      plotImg <- emphasize(FALSE, res, dfobj, Cells, pres, genesets, seu, presTab)
      output$CellPlot <- renderPlot(plotImg)
      output$SlingShot <- renderPlot(emphasizeSlingShot(res, dfobj, Cells, pres, genesets, seu, presTab))
      output$Comparison <- renderPlot(emphasizePathway(res, dfobj, Cells, pres, genesets, seu, presTab))

      output$legendTable <- DT::renderDataTable(
        buildGradientLegend(res, Cells = Cells)
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
        filename = "mylegend.png",
        content = function(file) {
          buildLegend(res, img = TRUE, name = file, GroupInfo = GroupInfo)
        }
      )
    })

    # clear list in Cell tab
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
    observeEvent(input$freqbp, {
      if (input$freqbp == 0) {
        return(NULL)
      }

      output$biPlot <- renderPlot(buildbiplot(input$biFont, input$biX, input$biY, genesets, TOPN = input$biCount,
                                              gsFont = input$gsFont, axtxt = input$axtxt, axlab = input$axlab,
                                              oddratio = FALSE))
      output$biplotdn <- downloadHandler(
        filename = function() {
          "mybiplot.png"
        },
        content = function(file) {
          ggsave(file, device = "png")
        }
      )
    })

    observeEvent(input$orbp, {
      if (input$orbp == 0) {
        return(NULL)
      }

      output$biPlot <- renderPlot(buildbiplot(input$biFont, input$biX, input$biY, genesets, TOPN = input$biCount,
                                              gsFont = input$gsFont, axtxt = input$axtxt, axlab = input$axlab,
                                              oddratio = TRUE))
      output$biplotdn <- downloadHandler(
        filename = function() {
          "mybiplot.png"
        },
        content = function(file) {
          ggsave(file, device = "png")
        }
      )
    })

    observeEvent(input$slingDraw, {
      if (input$slingDraw == 0) {
        return(NULL)
      }

      output$SlingShot <- renderPlot(buildSlingShot(seu, input$plotOption))

      # output$slingDownload <- downloadHandler(
      #   filename = function() {
      #     "myslingshot.png"
      #   },
      #   content = function(file) {
      #     ggsave(file, device = "png")
      #   }
      # )
    })

    observeEvent(input$slingRedraw, {
      if (input$slingRedraw == 0) {
        return(NULL)
      }

      res <- input$sortList
      output$biPlot <- renderPlot(buildbiplot(input$biFont, input$biX, input$biY, genesets, TOPN = input$biCount,
                                              gsFont = input$gsFont, axtxt = input$axtxt, axlab = input$axlab,
                                              oddratio = FALSE))
      output$SlingShot <- renderPlot(emphasizeSlingShot(res, dfobj, Cells, pres, genesets, seu, presTab))
      output$Comparison <- renderPlot(emphasizePathway(res, dfobj, Cells, pres, genesets, seu, presTab))
      # output$slingDownload <- downloadHandler(
      #   filename = function() {
      #     "myslingshot.png"
      #   },
      #   content = function(file) {
      #     ggsave(file, device = "png")
      #   }
      # )
    })

  }

  ui <- CellEnrichUI()

  shiny::shinyApp(ui, server, options = list(launch.browser = TRUE))
}

buildGradientLegend <- function(sortList, img = FALSE, name = NULL, Cells) {
  colV <- getColv(Cells)

  rlobj <- data.frame(stringsAsFactors = FALSE)
  for (i in 1:length(sortList)) {
    kk <- strsplit(sortList[[i]], " @")[[1]]
    Scale <- ""
    Pathway <- kk[1]
    Group <- kk[2]
    rlobj <- rbind(rlobj, cbind(Scale, Pathway, Group))
  }
  colnames(rlobj) <- c("Scale", "Pathway", "Group")
  currentGroup <- as.vector(rlobj$Group)
  # currentGroup <- rlobj$Group

  if (img) {
    png(name)
    plot(NULL, xaxt = "n", yaxt = "n", bty = "n", ylab = "", xlab = "", xlim = 0:1, ylim = 0:1)
    legend(
      "center",
      legend = rlobj$Pathway,
      pch = 15, pt.cex = 2, cex = 1.2, bty = "n",
      col = colV[currentGroup]
    )
    dev.off()
    return()
  }

  rlobj$Pathway <- paste0(
    sapply(
      colV[currentGroup],
      function(i) {
        paste0('<div style="background: ', i, '; display: inline-block; width: 1em;height: 1em;"></div>')
      }
    ),
    " ", as.character(rlobj$Pathway)
  )

  rlobj$Scale <- paste0(
    sapply(
      colV[currentGroup],
      function(i) {
        paste0(
          '<a style="color:black">low </a>',
          '<div style="background: linear-gradient(to right, ', col2hcl(i, l=0),' -40%, ', col2hcl(i),' 100%); display: inline-block; width: 8em;height: 1em;"></div>',
          '<div style="background: linear-gradient(to right, ', col2hcl(i),' 0%, ', col2hcl(i, l=100),' 95%); display: inline-block; width: 4em;height: 1em;"></div>',
          '<a style="color:black"> high</a>'
        )
      }
    ),
    " ", as.character(rlobj$Scale)
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

buildLegend <- function(sortList, img = FALSE, name = NULL, GroupInfo) {
  colV <- getColv(GroupInfo)

  rlobj <- data.frame(stringsAsFactors = FALSE)
  for (i in 1:length(sortList)) {
    kk <- strsplit(sortList[[i]], " @")[[1]]
    Pathway <- kk[1]
    Group <- kk[2]
    rlobj <- rbind(rlobj, cbind(Pathway, Group))
  }

  colnames(rlobj) <- c("Pathway", "Group")
  currentGroup <- rlobj$Group

  if (img) {
    png(name)
    plot(NULL, xaxt = "n", yaxt = "n", bty = "n", ylab = "", xlab = "", xlim = 0:1, ylim = 0:1)
    legend(
      "center",
      legend = rlobj$Pathway,
      pch = 15, pt.cex = 2, cex = 1.2, bty = "n",
      col = colV[currentGroup]
    )
    dev.off()
    return()
  }

  rlobj$Pathway <- paste0(
    sapply(
      colV[currentGroup],
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
