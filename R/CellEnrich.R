## 23.11.16
if(!require(farver)){
  install.packages('waiter') # install 'waiter' if not installed.
}
library(waiter)

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

getTU <- function(CountData, GroupInfo, plotOption='UMAP', topdims= 50) {
  cat("Mapping is started\n")
  library(Seurat)

  # CountData is normalized
  seu <- CreateSeuratObject(CountData)
  
  nfeat <- nrow(CountData) 
  nfeat <- min(100*round(nfeat/100,0), 3000)
  
  seu <- NormalizeData(seu)
  seu@assays$RNA@layers$data <- seu@assays$RNA@layers$counts # restoration
  seu <- ScaleData(seu)
  seu <- FindVariableFeatures(seu, nfeatures = nfeat)

  # Add cell type annotation to metadata
  seu <- AddMetaData(seu, GroupInfo, col.name = "cell_type")
  
  # Dimension reduction
  # These are now standard steps in the Seurat workflow for visualization and clustering
  seu <- RunPCA(seu, npcs=topdims, verbose = FALSE)

  # TSNE
  seu <- RunTSNE(seu, dims = 1:topdims)

  # UMAP
  seu <- RunUMAP(seu, dims = 1:topdims, uwot.sgd = TRUE)
  seu <- FindNeighbors(seu, verbose = FALSE, dims = 1:30)
  seu <- FindClusters(seu, algorithm = 3, random.seed = 7968, resolution = 0.5)
  
  # DimPlot(seu)
  return (seu)
}

library(scMerge)
gnm <- function(v) {
  out <- scMerge:::gammaNormMix(as.matrix(v), plot = FALSE )
  mat_prob <- matrix(out$probExpressed, nrow(v), ncol(v))
  mat_discretised <- 1 * (mat_prob > 0.5)
  return(mat_discretised)
}

findSigGenes <- function(v, method = "CellEnrich - HALF median", Name) {
  if (!method %in% c("CellEnrich - HALF median", "CellEnrich - median", "CellEnrich - mixture")) stop("wrong method")

  cat("findSigGenes started\n")
  rownames(v) <- colnames(v) <- NULL

  res <- list()
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
    
    medfull <- function(v) {
      v <- v[which(v > 0)]
      return(median(v))
    }

    if (method == "CellEnrich - HALF median") {
      for (i in 1:ncol(v)) {
        res[[i]] <- which(v[, i] > med2(v[, i]))
      }
    } else if (method == "CellEnrich - median") {
      for (i in 1:ncol(v)) {
        res[[i]] <- which(v[, i] > medfull(v[, i]))
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
    ) %>% select(Group, Top, genes, FDR) %>%
      filter(FDR <= q0) %>%
      # filter(Top <= TopCutoff) %>%
      arrange(FDR)
    
    res <- rbind(res, G)
  }
  res$genes <- as.character(res$genes)
  res$Group <- as.character(res$Group)
  res$FDR <- round(as.numeric(res$FDR), 4)
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

  colnames(CellPathwayDF) <- c("Cell", "Geneset", "Frequency")
  CellPathwayDF$Cell <- as.character(CellPathwayDF$Cell)
  CellPathwayDF$Geneset <- names(genesets)[as.numeric(as.character(CellPathwayDF$Geneset))]

  CellPathwayDF$Frequency <- as.numeric(as.character(CellPathwayDF$Frequency))

  # ------ add length column

  # Length <- getlgs(CellPathwayDF$Geneset)
  Size <- getlgs(genesets[as.character(CellPathwayDF$Geneset)])
  CellPathwayDF <- cbind(CellPathwayDF, Size)

  # ------ select genesets with count > 1

  if(length(pres) != length(Cells)){
    CellPathwayDF <- CellPathwayDF %>%
      dplyr::filter(Frequency > 1)
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
      thisCell <- Cells[i]  # this cell means this cell type
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
  colnames(res) <- c("Cell", "Geneset", "Pvalue") 

  res$Cell <- as.character(res$Cell)
  res$Geneset <- as.character(res$Geneset)
  res$Pvalue <- as.numeric(as.character(res$Pvalue))  
  
  res$Qvalue <- res$Pvalue
  res$Qvalue[which(res$Qvalue < 1e-20)] <- 1e-20
  res$Qvalue <- round(-log10(res$Qvalue), 4)  ## -log10(p-value) ->> Q-value

  return(res)
}

# pres : which gene-sets are significant for each cells.
# pres2 : for each gene-sets, how many cells are significant that gene-sets.

# 전체 그룹에서 유의한 회수 20 # pres2[genesets[i]]
# 특정 그룹에서 유의한 회수 6 # pres2[thiscell[idx]

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

  # res <- res %>% filter(OddRatio > 1)

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
  p <- p + theme(axis.title.y = element_text(size = 15, vjust= 0.5))
  p <- p + theme(axis.text = element_text(size = 12))

  p <- p + theme(
    axis.title.x = element_text(size=15),
    panel.background = element_rect(fill = 'white', colour = 'white'),
    panel.border = element_rect(fill = NA, colour = 'black', size=0.25 ),
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

    k <- sum(pathways) # selected ball

    gt <- sapply(1:length(pathways), function(j) {
      q <- pathways[j] # selected white ball, 1
      m <- unname(pres2Idx[names(pathways[j])]) # total white ball, 28
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

default_genesets <- c(
  "Reactome", # 2022
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
  
  if(!require(farver)){
    install.packages('farver') # install 'farver' if not installed.
    
  }
  library(farver)

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
                    HTML("<font color='black'>CellEnrich - HALF median (coef = 0.5)</font>"),
                    HTML("<font color='black'>CellEnrich - Median (coef = 1)</font>"),
                    HTML("<font color='black'>CellEnrich - Mixture</font>")
                    # HTML("<font color='black'>Fisher</font>")
                  ),
                  # choiceValues = c("CellEnrich - median", "CellEnrich - mixture", "Fisher"),
                  choiceValues = c("CellEnrich - HALF median", "CellEnrich - median", "CellEnrich - mixture"),
                  selected = "CellEnrich - HALF median",
                  # color = "#1976d2"
                )
              ),
              material_card(
                radioButtons(
                  "plotOption",
                  label = HTML("<font color='black' size='5'>Scatter Plot</font>"),#"Scatter Plot",
                  choiceNames = list(
                    tags$span(style = "color:black", "PCA"),
                    tags$span(style = "color:black", "TSNE"),
                    tags$span(style = "color:black", "UMAP")
                  ),
                  choiceValues = c("PCA", "TSNE", "UMAP"),
                  selected = "UMAP",
                  # color = "#1976d2"
                ),
                
                material_number_box(
                  input_id = "topdims",
                  label = HTML("<font color='black' size='5'>Top-N dims</font>"), # top dims
                  min_value = 30,
                  max_value = 100,
                  initial_value = 50,
                  step_size = 10
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
                    HTML("<font color='black'>Reactome</font>"),
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
                  selected= default_genesets[1] #"Reactome_2022"
                ),
                sidebarPanel(
                  tags$style("
                               .btn-file {
                                  background-color:#1976d2;
                                  border-color: none;
                               }
                               "),
                  HTML("<font color='black' size='3'>User defined geneset</font>"),
                  fileInput("other", "", placeholder = "RData format is required!"),
                )
              ),
              width = 4
            )
          ),
          solvedButton(
            inputId = "StartCellEnrich",
            label = "RUN",
            style = "margin-left:40%; background-color: #1976d2; height:60px; width:360px; font-size : 32px;",
            onClick = 'console.log("CellEnrich");'
          ),
          depth = 3
        ),
        width = 12, # 6
      )
    ),

    # tSNE/UMAP plot
    material_row(
      material_column(

        material_card(
          title = shiny::tags$h4("Scatter & Bar"), depth = 3,

          # Bar
          material_row(
            material_column(width = 6,
                            plotOutput("CellScatter", height = "480px")
            ),
            material_column(width = 6,
                          highchartOutput("CellBar", height = "480px")
            )
          ),
          # Cell - Comparison
          material_row(
            material_column(width = 6,
                            plotOutput("CellPlot", height = "480px")
            ),
           
            material_column(width = 6,
                            plotOutput("Comparison", height = "480px") # cell distribution
            )
          ),
          material_row(
            material_column(
              material_button("sigbtn", "Odds Ratio", icon = "sort", color = "blue darken-2"), 
              material_button("freqbtn", "Frequency", icon = "menu", color = "blue darken-2"),
              width = 3
            ), 
            material_column(
              shiny::downloadButton("sctdn", "Save Scatter", style = "background-color : #616161 !important"),
              shiny::downloadButton("imgdn", "Pathway Significance in each Group", style = "background-color : #616161 !important"),
              shiny::downloadButton("cmpdn", "Pathway Significance in whole data", style = "background-color : #616161 !important")
              , width = 6
            )
          ),
          material_row(
            material_column( 
                            width = 3
            ),
            material_column(
              shiny::downloadButton("sppcdn", "Cell-Pathway (P-values)", style = "background-color : #616161 !important"),
              shiny::downloadButton("tbldn", " All Significant Pathways", style = "background-color : #616161 !important")
              , width = 6
            )

          ),
          material_row(
            material_card(
              title = "",
              DT::dataTableOutput("legendTable"),
              shiny::downloadButton("legenddn", "Save Legend", style = "background-color : #616161 !important; display:none;")
            )
          ),
          material_row(
            material_card(
              title = shiny::tags$h4("User chosen Pathways"), divider = TRUE,
              my_ranklist<-rank_list(text = "", labels = "Switch any pathway position ONCE to activate the plot button",
                                     input_id = "sortList", css_id = "mysortableCell"),
              material_row(
                material_button("Emphasize", "Plot the Selected Pathways", color = "blue darken-2"),
                material_button("ClearList", "Clear List", color = "blue darken-2"),
                material_button("nonrun", "Please switch any element's position to activate the plot button", color = "orange")
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
          title = shiny::tags$h4("Biplot between top pathways and groups"), divider = TRUE,
          material_row(
            material_column(
              plotOutput("biPlot", height = "900px"),
              width = 10
            ),
            material_column(
              material_row(
                numericInput("biCount", label = "Top Pathways in each Group", value = 5, min = 1, max = 10, step = 1),
                numericInput("biFont", label = "Cell Type Label", value = 5, min = 1, max = 10, step = 1),
                numericInput("gsFont", label = "Pathway Label", value = 5, min = 1, max = 10, step = 1),
                numericInput("biX", label = "Range of X-axis", value = 5, min = 1, max = 10, step = 1),
                numericInput("biY", label = "Range of Y-axis", value = 2, min = 1, max = 10, step = 1),
                numericInput("axlab", label = "Axes Title Font", value = 15, min = 10, max = 30, step = 1),
                numericInput("axtxt", label = "Axes Index Font", value = 13, min = 10, max = 30, step = 1),
                material_row(
                  material_button("orbp", "Odds Ratio based Plot", color = "blue darken-2")
                ),
                material_row(
                  material_button("freqbp", "Frequency based Plot", color = "blue darken-2")
                )
              ),
              material_row(
                shiny::downloadButton("biplotdn", "Save Biplot", style = "background-color : #616161 !important")
              ),
              width = 2
            )
          ),
        ),
        material_card(
          title = shiny::tags$h4("Heatmap between top pathways and groups"), divider = TRUE,
          material_row(
            material_column(
              plotOutput("heatPlot", height = "900px"),
              width = 10
            ),
            material_column(
              material_row(
                shiny::downloadButton("heatplotdn", "Save Heatmap", style = "background-color : #616161 !important")
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

emphasize <- function(path = FALSE, inputObj, dfobj, Cells, pres, genesets, seu, presTab, maptitle='Odds Ratio') {
  
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

  # absolute scale 0~0.001; 0.001~0.01; 0.01~0.05; 0.05~0.1; 0.1~1
  pval_scaling <- function(x){(x-min(x))/(max(x)-min(x))}
  cell_pval <- c(cell_pval, c(maxall, minall))
  cell_pval <- pval_scaling(cell_pval)
  cell_pval <-ifelse(is.nan(cell_pval),0. ,cell_pval)
  cell_pval <- cell_pval[1:(length(cell_pval)-2)]

  for (i in 1:length(cellidx)){
    tmp <- farver::decode_colour(colV[cellidx[i]], to="hcl")
    tmp[3] <- min (100, 100 * cell_pval[i])
    colV[cellidx[i]] <- encode_colour(tmp, from = 'hcl')
  }
  ap <- dfobj_new$col
  ap <- ifelse(ap=='#E5E5E5', 0.2, 1)

  rownames(dfobj_new) <- NULL
  graphString <- "ggobj2 <- ggplot(dfobj_new, aes(x = x, y = y)) + geom_point(colour = colV, alpha=ap) +
                            theme(panel.background = element_rect(fill = 'white', colour = 'white'),
                            panel.border = element_rect(colour = 'black', fill=NA, size=0.25),
                            plot.title = element_text(size = 18, face = 'bold'))"
  
  graphString <-  paste0(graphString, '+ggtitle(maptitle)')
  
  
  eval(parse(text = graphString))

  return(ggobj2)
}

emphasizePathway <-
  function(inputObj, dfobj, Cells, pres, genesets, seu, presTab, maptitle=' ') {
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
      tmp <- farver::decode_colour(colV[cellidx[i]], to="hcl")
      tmp[3] <- min (100, 100 * cell_pval[i])
      colV[cellidx[i]] <- encode_colour(tmp, from = 'hcl')
    }
    ap <- dfobj_new$col
    ap <- ifelse(ap=='#E5E5E5', 0.2, 1)

    rownames(dfobj_new) <- NULL
    graphString <- "ggobj2 <-  ggplot(dfobj_new, aes(x = x, y = y)) + geom_point(color = colV, alpha=ap ) +
                             theme(panel.background = element_rect(fill = 'white', colour = 'white'),
                              panel.border = element_rect(colour = 'black', fill=NA, size=0.25),
                            plot.title = element_text(size = 18, face = 'bold'))"
    
    graphString <-  paste0(graphString, '+ggtitle(maptitle)')


    eval(parse(text = graphString))
    return(ggobj2)
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
#' @import farver
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
    
    toptab <- function (genesets, TOPN = 5, oddratio = TRUE, myplot='biplot'){
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
          thisCellIdx <- which(GroupInfo == Cells[i])
    
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
          
          # adjusted to original OR values
          for (i in 1:length(Cells)){
            corGeneset <- dat$Geneset[dat$Cell == Cells[i]]
            corOddRatio <- dat$OddRatio[dat$Cell == Cells[i]]
            for (g in corGeneset){
              tab[g, Cells[i]] <- corOddRatio[corGeneset==g]
            }
          }
          
      }
      else { # FREQUENCY
        
        tab <- matrix(0, nrow = length(genesets), ncol = length(Cells))
        
        for (i in 1:length(Cells)) {
          thisCellIdx <- which(GroupInfo == Cells[i])
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
      
      saveRDS(tab, "ortab.rds")
      return(tab)
    }
    

    buildbiplot <- function(biFont, biX, biY, genesets, TOPN = 5, oddratio = TRUE, gsFont=5, axtxt=13, axlab=15,
                            myplot='biplot', tab=matrix(0, 0, 0)) {
    
      hmtype <- 'Odds Ratio based Biplot'
      if (!oddratio){
        hmtype <- 'Frequency based Biplot'
      }
      ##########  BiPlot
      labels <- rownames(tab)
      
      if (length(which(apply(tab, 2, var)==0)) > 0){
        model <- prcomp(tab)  # non scaling because of constant/zero column(s)
      } else {
        model <- prcomp(tab, scale. = T)
      }
      
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
        labs(y = 'PC2', x = 'PC1', title = hmtype) +
        theme(   plot.title = element_text(size=20, hjust = 0.5),
                 axis.line = element_blank(), #element_line(color = "black", size = 0.5, linetype = "solid"),
                 axis.title.x = element_text(size=axlab),
                 axis.title.y = element_text(size=axlab), #element_blank(),
                 axis.text.x = element_text(size=axtxt, colour = 'black', hjust = 1),
                 axis.text.y = element_text(size=axtxt, colour = 'black', hjust = 1),
                 panel.grid.major = element_line(colour = "grey86"),
                 panel.background = element_blank(),
                 panel.border = element_blank(),
                 panel.spacing.x = unit(0.5, "lines"),
                 panel.spacing.y = unit(1, "lines"),
                 strip.text = element_text(size=17, color="black"),
                 strip.background.x = element_rect(fill="#CDE8DF"),
                 legend.position="none")

     
      return(BiPlot)
    }
    
    buildheatplot <- function(biFont, biX, biY, genesets, TOPN = 5, oddratio = TRUE, gsFont=5, axtxt=13, axlab=15, 
                              myplot='biplot', tab=maxtrix(0, 0, 0)) {
      
      ###########   Heatmap 
      hmtype <- 'Odds Ratio based Heatmap'
      if (!oddratio){
        hmtype <- 'Frequency based Heatmap'
      }
      
      # First define your breaks
      col_breaks <- c(seq(0, ceiling(max(tab)), length.out=101))
      my_palette <- c(colorRampPalette(c("white", "red"))(length(col_breaks)-1))
      
      mat_data <- round(tab, 2) 
      
      library(gplots)
      HeatPlot <<- heatmap.2(mat_data,
                main = hmtype, # heat map title
                density.info="none",  # turns off density plot inside color legend
                trace="none",         # turns off trace lines inside the heat map
                margins =c(15,70),      # widens margins around plot
                col=my_palette,       # use on color palette defined earlier
                scale = "none",
                breaks=col_breaks,    # enable color transition at specified limits
                dendrogram="row",    # only draw a row dendrogram
                Colv=NA,            # turn off column clustering
                
                # additional control of the presentation
                lhei = c(2, 13),       # adapt the relative areas devoted to the matrix
                lwid = c(2, 10),
                cexRow = 1.5,
                cexCol = 2,
                key.title = NA,
                key.xlab = NA,
                key.ylab = NA,
                key.par = list(mar=c(4, 1, 2, 1),
                               mgp=c(1.5, 0.5, 0),
                               cex=1)
      )  
      # saveRDS(HeatPlot, 'hm.rds')
      return(HeatPlot)
    }
    
    usergs <- ""
    
    observeEvent(input$other,{
      other_geneset <- input$other
      ext <- tools::file_ext(other_geneset$datapath)
      
      req(other_geneset)
      inFile <- other_geneset$datapath
      
      otherVal <- other_geneset$name
      usergs <<- inFile
      updatedValues <- c(default_genesets, otherVal)
   
      updateRadioButtons(session, "genesetOption", choices = updatedValues, selected= otherVal)
      print(usergs)
      cat("Geneset added!")
    })
   
    ### CODES
    # variable initialize

    dtobj <- dfobj <- pres <- pres2 <- presTab <- ""
    CellPathwayDF <- ""
    gt <- Cells <- A <- ""
    CellScatter <- ""
    CellHistogram <- ""
    BiPlot <- HeatPlot <- OR <- ""
    

    observeEvent(input$StartCellEnrich, {
      pt <- proc.time()

      # ------ Hide Start Button
      shinyjs::disable("StartCellEnrich")

      #Disable Emphasize 1st
      shinyjs::runjs('$("#Emphasize").attr("disabled",true)')

      # ------ Load Genesets
      # Check that data object exists and is data frame.
      if (is.null(genesets)) {
        if (input$genesetOption == "Reactome") load("Reactome_2022.RData")
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
        genesets <- read_Input_geneset(usergs)
        
        if (is.null(genesets)){
          shiny::showNotification("Geneset file is invalid!", type = "error", duration = 20)
          return(NULL)
        }
      }

      if (is.null(genesets)) {
        shiny::showNotification("Geneset file is missing!", type = "error", duration = 20)
        return(NULL)
      }
      
      if (length(rownames(CountData))==0) {
        shiny::showNotification("This dataset has not enough significant genes!", type = "error", duration = 60)
        return(NULL)
        stop("This dataset has not enough significant genes!")
      }
      
      if (length(genesets)==0) {
        shiny::showNotification("This dataset has no significant pathway detected based on current genesets!", type = "error", duration = 60)
        return(NULL)
        stop("This dataset has no significant pathway detected!")
      }

      genesets <<- genesets
      
      # ------ for test
      q0 <- input$qvalueCutoff

      # ------ Create new Waitress
      w <- Waitress$new(selector = NULL, theme = "overlay-radius")

      genes <- rownames(CountData)
      genesets <- GenesetFlush(genes, genesets)
      lgs <- getlgs(genesets)

      # ------ Genesetsize Flush
      genesets <- GenesetsizeFlush(genesets, lgs, input$minGenesetSize, input$maxGenesetSize)
      
      # ------ Gene Flush
      remgenes <- GeneFlush(genes, genesets)
      CountData <- CountData[!(rownames(CountData) %in% names(remgenes)),]
      genes <- genes[!genes %in% names(remgenes)]
      rm(remgenes)

      genesets <<- genesets
      
      if (length(rownames(CountData))==0) {
        shiny::showNotification("This dataset has not enough significant genes!", type = "error", duration = 60)
        return(NULL)
        stop("This dataset has not enough significant genes!")
      }
      
      if (length(genesets)==0) {
        shiny::showNotification("This dataset has no significant pathway detected based on current genesets!", type = "error", duration = 60)
        return(NULL)
        stop("This dataset has no significant pathway detected!")
      }

      # ------ Background genes
      A <<- getBackgroundGenes(genesets)

      # ------ Calculate TSNE / UMAP First
      # library(Matrix)

      seu <- getTU(CountData, GroupInfo, input$plotOption, topdims=input$topdims)
      seu <<- seu
      
      #PCA
      if (input$plotOption == "PCA") {
        dfobj <- data.frame(Embeddings(seu, 'pca')[,1:2], col = GroupInfo, stringsAsFactors = FALSE)
      }
      
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

      cat("running gc\n")
      gc()

      # ------ Find Significant Genes with Fold Change
      if (input$FCoption != "GSVA") {
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
                        columnDefs = list(list(className = "dt-center", targets = 0:3))
                      ),
                      selection = "none",
        )
      )
      tmp_df <- data.frame()

      cat("biobj \n")
      biobj <- getbiobj(genes, genesets)

      if (length(s) == 0) {
        tmp_cells <- unique(s2$Group) ## cells containing significant genes
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
          sigidx <- which(p.adjust(tmp_pv, "fdr") <= q0)  # q-values cutoff 0.05
          tmp_pres[sigidx, i] <- unname(tG[tmp_cells[i]])
          
          tmp_pv[which(tmp_pv < 1e-20)] <- 1e-20  # -12

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
      presTab_pval <- presTab <- c()
      

      if (length(s) >= 100) {
        w$start()
        for (i in 1:lens) {
          if (i %% lens100 == 0) w$inc(1)
          prespv <- getHyperPvalue(rc[s[[i]]], genesets, A, lgs, q0, biobj)

          pres[[i]] <- which(p.adjust(prespv, "fdr") <= q0)

          presTab_pval <- cbind(presTab_pval, prespv)
        }
        w$close()
        colnames(presTab_pval) <- colnames(CountData)
      }
      else {
        for (i in 1:lens) {
          prespv <- getHyperPvalue(rc[s[[i]]], genesets, A, lgs, q0, biobj)

          pres[[i]] <- which(p.adjust(prespv, "fdr") <= q0)
          
          presTab_pval <- cbind(presTab_pval, prespv) ### presTab_pval-> hyper P-value
        }
        colnames(presTab_pval) <- names(s)
      }
      rownames(presTab_pval) <- names(genesets)
      
      # presTab_pval
      output$sppcdn <- downloadHandler(
        filename = "cell_pw_pval.csv",
        content = function(file) {
          write.csv(presTab_pval, file, row.names = TRUE)
        }
      )
      
      cat("pres defined\n")
      presTab <- presTab_pval ### presTab -log10(p-value) -> Q-value
      presTab[which(presTab < 1e-20)] <- 1e-20
      presTab <- round(-log10(presTab), 4) ### presTab -log10(p-value) -> Q-value
      presTab <<- presTab

      # pres : which gene-sets are significant for each cells.
      pres <<- pres
      
      # pres2 : for each gene-sets, how many cells are significant that gene-sets.
      cat("pres2\n")
      pres2 <- sort(table(unlist(pres)), decreasing = T) 
      if (length(s) != 0) {
        names(pres2) <- names(genesets)[as.numeric(names(pres2))]
      }
      pres2 <<- pres2
      
      # ------ CellPathwayDF
      CellPathwayDF <- buildCellPathwayDF(GroupInfo, pres, genesets)
      
      PP <- pathwayPvalue(GroupInfo, pres, pres2, genesets) 

      # OR -> # Group / PATHWAY / ODDRATIO        
      if (nrow(tmp_df) > 0) {
        OR <- tmp_df
        colnames(OR) <- c("Cell", "Geneset", "OddRatio")
        OR <<- OR
      }
      else {
        OR <<- getOddRatio(GroupInfo, pres, pres2, genesets, input$ORratio)
      }
      
      #### Ceiling the max OR 
      non_zero_or <- OR$OddRatio[OR$OddRatio > 0]
      outlier_or <- round(mean(non_zero_or) + 2*sd(non_zero_or), 4)
      
      # retrieval
      extreme_values <- OR$OddRatio[OR$OddRatio > outlier_or]
      extreme_min <- max(OR$OddRatio[OR$OddRatio < min(extreme_values)])
      # rescale extreme values respecting to the ranks
      variation <- 1:length(extreme_values) /10 + sapply(1:length(extreme_values), function(i){ return (rbeta(1, 2, 8))})
      converted_extreme <- round(extreme_min + 0.1*sd(non_zero_or) + variation, 4)
      ordx <- order(extreme_values)
      OR$OddRatio[OR$OddRatio > outlier_or] <- sort(converted_extreme)[order(ordx)]
      
      
      CellPathwayDFP <- CellPathwayDF %>%
        inner_join(PP) %>% inner_join(OR)  %>% select(-Pvalue)  # %>% filter(Qvalue > 1) 
      
      CellPathwayDFP <- CellPathwayDFP %>% filter(Qvalue > -log10(q0)) 
      output$tbldn <- downloadHandler(
        filename = "all_sig_pathways.csv",
        content = function(file) {
          write.csv(CellPathwayDFP, file, row.names = FALSE)
        }
      )

      ## significant OR > 1
      OR <<- OR %>% filter(OddRatio > 1)
      CellPathwayDF <- CellPathwayDFP %>% filter(OddRatio > 1)
      
      CellPathwayDF <<- CellPathwayDF

      
      CellMarkers <- data.frame()
      Cells <- sort(unique(GroupInfo))
      Cells <<- Cells

      for (i in 1:length(Cells)) {
        thisCellPathways <- CellPathwayDF %>%
          filter(Cell == Cells[i]) %>%
          dplyr::select(Geneset)

        thisCellDEs <- s2 %>%
          filter(Group == Cells[i]) %>%
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

        additive <- data.frame(cbind(genes, Count, Group = Cells[i]))
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
        CellMarkers$Group <- as.factor(CellMarkers$Group)
        CellMarkers$Count <- as.numeric(CellMarkers$Count)

        output$markerL2 <- DT::renderDataTable(
          DT::datatable(CellMarkers,
                        rownames = FALSE,
                        filter = "top",
                        options = list(
                          autoWidth = TRUE,
                          dom = "ltp",
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

      output$CellScatter <- renderPlot(CellScatter)

      output$sctdn <- downloadHandler(
        filename = "scatter.pdf",
        content = function(file) {
          ggsave(file, CellScatter, device = "pdf")
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
                      height = "100%"
                      # height = "500px"
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
                width = "100%" # maximum 3 table in column
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
          "order = list(list(4,'desc'))", # odds ratio based
          "),",
          "rownames = FALSE,",
          "selection = 'single')",
          ")"
        )
        eval(parse(text = t))
      }

      # StartEnrich Finished
      shinyjs::click('ClearList')
      shinyjs::click('sigbtn')
      shinyjs::click('orbp')
      shinyjs::click('heat_orbp')
      
      print(proc.time() - pt)
      
      # ------ Enable Start Button Again
      shinyjs::enable("StartCellEnrich")
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
          thisCellData <- thisCellData %>% top_n(1, wt = Frequency)
          if (nrow(thisCellData) >= 1) {
            thisCellData <- thisCellData %>% top_n(-1, wt = Size)
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

      plotImg <- emphasize(FALSE, res, dfobj, Cells, pres, genesets, seu, presTab, maptitle='Frequency')
      compareImg <- emphasizePathway(res, dfobj, Cells, pres, genesets, seu, presTab)
      
      output$cmpdn <- downloadHandler(
        filename = function() {
          "frequency_global.pdf"
        },
        content = function(file) {
          ggsave(file, compareImg, device = "pdf")
        }
      )

      output$imgdn <- downloadHandler(
        filename = function() {
          "frequency.pdf"
        },
        content = function(file) {
          ggsave(file, plotImg, device = "pdf")
        }
      )
      
      shinyjs::show("legenddn")
      
      output$legenddn <- downloadHandler(
        filename = "mylegend.pdf",
        content = function(file) {
          buildLegend(res, img = TRUE, name = file, GroupInfo = GroupInfo)
        }
      )

      output$CellPlot <- renderPlot(plotImg)
      output$Comparison <- renderPlot(compareImg)
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

        thisCellData <- CellPathwayDF %>% dplyr::filter(Cell == Cells[i])
        
        if (nrow(thisCellData) >= 1) {
          thisCellData <- thisCellData %>% top_n(1, wt = OddRatio)
        
          if (nrow(thisCellData) >= 1) {
            thisCellData <- thisCellData %>% top_n(-1, wt = Size)
            if (nrow(thisCellData) >= 1) {
              thisCellData <- thisCellData %>% top_n(1, wt = Frequency)
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
      compareImg <- emphasizePathway(res, dfobj, Cells, pres, genesets, seu, presTab)
      
      output$cmpdn <- downloadHandler(
        filename = function() {
          "oddsratio_global.pdf"
        },
        content = function(file) {
          ggsave(file, compareImg, device = "pdf")
        }
      )

      output$imgdn <- downloadHandler(
        filename = function() {
          "oddsratio.pdf"
        },
        content = function(file) {
          ggsave(file, plotImg, device = "pdf")
        }
      )

      output$legenddn <- downloadHandler(
        filename = "mylegend.pdf",
        content = function(file) {
          buildLegend(res, img = TRUE, name = file, GroupInfo = GroupInfo)
        }
      )

      output$CellPlot <- renderPlot(plotImg)
      output$Comparison <- renderPlot(compareImg)
    })
    
    
    # ##### trigger value
    
    for(i in 1:length(Cells)){
      tmp <- paste0("toSortButton", i)
      observeEvent(input[[tmp]],{
        rv$select_clicked <- rv$select_clicked + 1
      })
    }
    
    observeEvent(input$sortList, {
      shinyjs::runjs('$("#Emphasize").attr("disabled",false)')
    })

    observeEvent(input$Emphasize, {
      if (input$Emphasize == 0) { # prevent default click state
        return(NULL)
      }
      shinyjs::runjs('$("#Emphasize").attr("disabled",true)')
      res<-input$sortList
      
      plotImg <- emphasize(FALSE, res, dfobj, Cells, pres, genesets, seu, presTab)
      compareImg <- emphasizePathway(res, dfobj, Cells, pres, genesets, seu, presTab)
      
      output$CellPlot <- renderPlot(plotImg)
      output$Comparison <- renderPlot(compareImg)

      output$legendTable <- DT::renderDataTable(
        buildGradientLegend(res, Cells = Cells)
      )

      output$imgdn <- downloadHandler(
        filename = function() {
          "my_pw_group.pdf"
        },
        content = function(file) {
          ggsave(file, plotImg, device = "pdf")
        }
      )
      
      output$cmpdn <- downloadHandler(
        filename = function() {
          "my_pw_global.pdf"
        },
        content = function(file) {
          ggsave(file, compareImg, device = "pdf")
        }
      )

      output$legenddn <- downloadHandler(
        filename = "mylegend.pdf",
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
                                              oddratio = FALSE, 
                                              tab=toptab(genesets, TOPN = input$biCount, oddratio = FALSE)))
      
      output$heatPlot <- renderPlot(buildheatplot(input$biFont, input$biX, input$biY, genesets, TOPN = input$biCount,
                                              gsFont = input$gsFont, axtxt = input$axtxt, axlab = input$axlab,
                                              oddratio = FALSE,
                                              tab=toptab(genesets, TOPN = input$biCount, oddratio = FALSE)))
   
      output$biplotdn <- downloadHandler(
        filename = function() {
          "mybiplot.pdf"
        },
        content = function(file) {
          ggsave(file, device = "pdf")
        }
      )
      
      output$heatplotdn <- downloadHandler(
        filename = function() {
          "myheatmap.pdf"
        },
        content = function(file) {
          
          revRowInd <- match(c(1:length(HeatPlot$rowInd)), HeatPlot$rowInd)
          revColInd <- match(c(1:length(HeatPlot$colInd)), HeatPlot$colInd)
          
          pdf(file, width = 24, height = 16)
          heatmap.2(t(HeatPlot$carpet)[revRowInd, revColInd],
                    main = "Frequency based Heatmap", # heat map title
                    density.info="none",  # turns off density plot inside color legend
                    trace="none",         # turns off trace lines inside the heat map
                    margins =c(15,70),      # widens margins around plot
                    col=HeatPlot$col,       # use on color palette defined earlier
                    scale = "none",
                    breaks=HeatPlot$breaks,    # enable color transition at specified limits
                    dendrogram = "row",
                    Rowv = HeatPlot$rowDendrogram,
                    Colv=NA,            # turn off column clustering
                    
                    # # additional control of the presentation
                    lhei = c(2, 13),       # adapt the relative areas devoted to the matrix
                    lwid = c(2, 10),
                    cexRow = 1.5,
                    cexCol = 2,
                    key.title = NA,
                    key.xlab = NA,
                    key.ylab = NA,
                    key.par = list(mar=c(5, 1, 2, 1),
                                   mgp=c(1.5, 0.5, 0),
                                   cex=1)
          )  
          
          dev.off()
        }
      )
      
    })### input$freqbb

    observeEvent(input$orbp, {
      if (input$orbp == 0) {
        return(NULL)
      }
      
      output$biPlot <- renderPlot(buildbiplot(input$biFont, input$biX, input$biY, genesets, TOPN = input$biCount,
                                              gsFont = input$gsFont, axtxt = input$axtxt, axlab = input$axlab,
                                              oddratio = TRUE,
                                              tab = toptab(genesets, TOPN = input$biCount, oddratio = TRUE)))

      output$heatPlot <- renderPlot(buildheatplot(input$biFont, input$biX, input$biY, genesets, TOPN = input$biCount,
                                                gsFont = input$gsFont, axtxt = input$axtxt, axlab = input$axlab,
                                                oddratio = TRUE,
                                                tab=toptab(genesets, TOPN = input$biCount, oddratio = TRUE)))
      
      output$biplotdn <- downloadHandler(
        filename = function() {
          "mybiplot.pdf"
        },
        content = function(file) {
          ggsave(file, device = "pdf")
        }
      )
      
      output$heatplotdn <- downloadHandler(
        filename = function() {
          "myheatmap.pdf"
        },
        content = function(file) {
          
          revRowInd <- match(c(1:length(HeatPlot$rowInd)), HeatPlot$rowInd)
          revColInd <- match(c(1:length(HeatPlot$colInd)), HeatPlot$colInd)
          
          pdf(file, width = 24, height = 16)
          heatmap.2(t(HeatPlot$carpet)[revRowInd, revColInd],
                    main = "Odds Ratio based Heatmap", # heat map title
                    density.info="none",  # turns off density plot inside color legend
                    trace="none",         # turns off trace lines inside the heat map
                    margins =c(15,70),      # widens margins around plot
                    col=HeatPlot$col,       # use on color palette defined earlier
                    scale = "none",
                    breaks=HeatPlot$breaks,    # enable color transition at specified limits
                    dendrogram = "row",
                    Rowv = HeatPlot$rowDendrogram,
                    Colv=NA,            # turn off column clustering
                    
                    # # additional control of the presentation
                    lhei = c(2, 13),       # adapt the relative areas devoted to the matrix
                    lwid = c(2, 10),
                    cexRow = 1.5,
                    cexCol = 2,
                    key.title = NA,
                    key.xlab = NA,
                    key.ylab = NA,
                    key.par = list(mar=c(5, 1, 2, 1),
                                   mgp=c(1.5, 0.5, 0),
                                   cex=1)
          )  
          
          dev.off()
        }
      )
      
    })  ### input$orbp
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

  if (img) {
    pdf(name, width=16, height = 9)
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
          '<a style="color:black"> High </a>',
          '<div style="background: linear-gradient(to right, ', col2hcl(i, l=0),' -40%, ', col2hcl(i),' 100%); display: inline-block; width: 8em;height: 1em;"></div>',
          '<div style="background: linear-gradient(to right, ', col2hcl(i),' 0%, ', col2hcl(i, l=100),' 95%); display: inline-block; width: 4em;height: 1em;"></div>',
          '<a style="color:black"> Low </a>'
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
    pdf(name, width=16, height = 9)
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
