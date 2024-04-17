## 24.04.04
if(!require(waiter)){
  install.packages('waiter') # install 'waiter' if not installed.
}
if(!require(Seurat)){
  install.packages('Seurat')
}
library(waiter)
library(Seurat)

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
  # CountData is normalized
  # seu <- CreateSeuratObject(pbmc)
  seu <- CreateSeuratObject(CountData)


  nfeat <- nrow(CountData)
  nfeat <- min(100*round(nfeat/100,0), 3000)

  seu <- NormalizeData(seu)
  seu@assays$RNA@layers$data <- seu@assays$RNA@layers$counts
  seu <- ScaleData(seu, do.center=FALSE)
  seu <- FindVariableFeatures(seu, nfeatures = nfeat)

  # Add cell type annotation to metadata
  seu <- AddMetaData(seu, GroupInfo, col.name = "cell_type")
  seu <- RunPCA(seu, npcs=topdims, verbose = FALSE)

  # TSNE
  seu <- RunTSNE(seu, dims = 1:topdims)
  # UMAP
  seu <- RunUMAP(seu, dims = 1:topdims, uwot.sgd = TRUE)

  seu <- FindNeighbors(seu, verbose = FALSE, dims = 1:30)
  seu <- FindClusters(seu, algorithm = 3, random.seed = 7968, resolution = 0.5)
  return (seu)
}

library(scMerge)
gnm <- function(v) {
  out <- scMerge:::gammaNormMix(as.matrix(v), plot = FALSE )
  mat_prob <- matrix(out$probExpressed, nrow(v), ncol(v))
  mat_discretised <- 1 * (mat_prob > 0.5)
  return(mat_discretised)
}

findSigGenes <- function(v, method = "CellEnrich - median", Name, coef=1) {

  if (!method %in% c("CellEnrich - median",
                     "CellEnrich - FGSEA")) stop("wrong method")

  cat("findSigGenes started\n")
  rownames(v) <- colnames(v) <- NULL

  res <- list()
  cat("scaling\n")

  cat("define Lists\n")
  med <- function(v, coef=1) {
    v <- v[which(v > 0)]
    return(median(v) * coef)
  }

  if (method == "CellEnrich - median") {
    for (i in 1:ncol(v)) {
      res[[i]] <- which(v[, i] > med(v[, i], coef=coef))
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

  pv <- sapply(1:length(genesets), function(i) {
    q <- biobj[i] # selected white ball
    m <- lgs[i] # white ball
    n <- A - m # black ball
    k <- lg # selected ball
    1 - phyper(q - 1, m, n, k)
  })
  return(pv)
}

buildCellPathwayDF <- function(GroupInfo, pres, genesets, pwFrequency=0.1) {
  # pwFrequency = 0.1
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
      # i = 1
      thisCell <- Cells[i]
      tt <- table(unlist(pres[which(thisCell == GroupInfo)]))

      if (nrow(tt)) {
        CellPathwayDF <- rbind(CellPathwayDF, cbind(thisCell, names(tt), unname(tt)))
      }
    }
  }

  colnames(CellPathwayDF) <- c("Cell", "Pathway", "Frequency")
  CellPathwayDF$Cell <- as.character(CellPathwayDF$Cell)
  CellPathwayDF$Pathway <- names(genesets)[as.numeric(as.character(CellPathwayDF$Pathway))]
  CellPathwayDF$Frequency <- as.numeric(as.character(CellPathwayDF$Frequency))

  nSample_per_group <- CellPathwayDF$Frequency
  for (i in 1:length(nSample_per_group)) {
    nid <- which(GroupInfo == CellPathwayDF$Cell[i])
    nSample_per_group[i] <- length(nid)
  }
  CellPathwayDF$N_cell <- paste0("/", as.character(nSample_per_group))

  # ------ add length column
  Size <- getlgs(genesets[as.character(CellPathwayDF$Pathway)])
  CellPathwayDF <- cbind(CellPathwayDF, Size)
  # ------ select genesets with count > 1

  CellPathwayDF <- CellPathwayDF %>% dplyr::filter(Frequency > 1)

  return(CellPathwayDF)
}

pathwayPvalue <- function(GroupInfo, pres, pres2, genesets) {
  # GroupInfo <- pbmcClustInfo
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
  colnames(res) <- c("Cell", "Pathway", "Pvalue")

  res$Cell <- as.character(res$Cell)
  res$Pathway <- as.character(res$Pathway)
  res$Pvalue <- as.numeric(as.character(res$Pvalue))
  res$Pvalue[is.na(res$Pvalue)] <- 1

  res$Qvalue <- res$Pvalue

  #### Extreme Qvalue
  common_qv <- res$Qvalue[res$Qvalue >= 1e-20]
  extreme_id <- which(res$Qvalue < 1e-20)
  extreme_qv <- res$Qvalue[extreme_id]
  extreme_min <- min(common_qv)

  res$Qvalue[extreme_id] <- 1e-20
  res$Qvalue <- round(-log10(res$Pvalue), 4)

  # rescale extreme values respecting to the ranks
  variation <- 1:length(extreme_qv) /10 + sapply(1:length(extreme_qv), function(i){ return (rbeta(1, 2, 8))})
  converted_extreme <- round(-log10(extreme_min) + 0.1*sd(-log10(common_qv)) + variation, 4)
  ordx <- order(extreme_qv)
  res$Qvalue[extreme_id] <- sort(converted_extreme)[order(ordx)]

  return(res)
}

# pres : which gene-sets are significant for each cells.
# pres2 : for each gene-sets, how many cells are significant that gene-sets.

# 전체 그룹에서 유의한 회수 20 # pres2[genesets[i]]
# 특정 그룹에서 유의한 회수 6 # pres2[thiscell[idx]

# 전체 그룹 Cell 수 : N
# 특정 그룹 Cell 수 : K

# Group_specific_OR = (6/(K-6)) / (14/(N-14))

getOddRatio <- function(GroupInfo, pres, pres2, genesets, pwFrequency) {
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
      if (B < length(thisCellIdx) * pwFrequency) {
        return(0)
      }
      A <- pres2[names(genesets)[k]] # 전체 Cell에서 유의한 회수
      if (is.na(A)) {
        return(0)
      }

      N <- total # 전체 Cell 수
      K <- length(thisCellIdx)

      if (is.na(N - A) || (N==A) ){
        return (0)
      }

      if (is.na(K - B) || (K==B) ){
        return (0)
      }

      return( (B/(K - B)) / (A/(N - A)) )
    }))

    OR <- round(OR, 4)
    # Cell, Pathway, OR
    res <- rbind(
      res,
      data.frame(
        Cell = as.character(thisCell),
        Pathway = as.character(names(genesets)),
        OddsRatio = as.numeric(OR), stringsAsFactors = FALSE
      )
    )
  }
  colnames(res) <- c("Cell", "Pathway", "OddsRatio")
  return(res)
}

buildDT <- function(pres2) {
  DT::datatable(
    data.frame(
      Pathway = names(pres2),
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

briterhex <- function(colors, alpha = 1.0) {
  if (length(colors) == 0) {
    stop("Length of color should be larger than zero --- briterhex")
  }
  res <- c()
  for (i in 1:length(colors)) {
    v <- as.vector(col2rgb(colors[i])) * alpha #* 1.3
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
  p <- p + theme(axis.title.y = element_text(size = 15, vjust= 0.5))
  p <- p + theme(axis.text = element_text(size = 12))

  p <- p + theme(
    axis.title.x = element_text(size=15),
    panel.background = element_rect(fill = 'white', colour = 'white'),
    panel.border = element_rect(fill = NA, colour = 'black', linewidth =0.25 ),
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

default_genesets <- c(
  "Human-Reactome", # 2022
  "Human-WikiPathway", # 2021
  "Human-KEGG", # KEGG 2021
  "Human-GOBP",
  "Human-GOCC",
  "Human-GOMF",
  "Mouse-Reactome",
  "Mouse-WikiPathway", # 2019
  "Mouse-KEGG", # 2019
  "Mouse-GOBP",
  "Mouse-GOCC",
  "Mouse-GOMF")

CellEnrichUI <- function(GroupInfo) {
  library(shinymaterial)
  library(highcharter)
  library(sortable)
  if(!require(farver)){
    install.packages('farver') # install 'farver' if not installed.

  }
  library(farver)

  Cells <- sort(unique(GroupInfo))
  CardColors <- briterhex(scales::hue_pal(h = c(20, 350), c = 100, l = 65, h.start = 0, direction = 1)(length(Cells)), alpha = 1.1)
  tab_ids <- lapply(1:length(Cells), function(i){
    paste0("tab_id_", i)
  })

  custom_css <- sprintf("
  .runbutton {
    text-align: center;
  }

  .tabs .tab a {
    font-weight: 800;
    font-size: 20px !important;
    text-align: center;
    color: #fff;
    border: 2px solid;
  }

  .tabs .tab a:hover {
    color: #00ffff;
    cursor: pointer;
  }

  .tabs .tab a[class*='active'] {
    background-color: #000; !important;
    color: #223344;
    border-color: #223344;
  }

    %s
  ", paste(sprintf(".tabs .tab a[href*=%s] { background-color: %s }", tab_ids, CardColors), collapse = "\n"))

  material_page(
    shinyjs::useShinyjs(),
    shinyFeedback::useShinyFeedback(feedback = TRUE, toastr = TRUE),
    # dynamic datatable full width

    tags$head(tags$style(type = "text/css", ".display.dataTable.no-footer{width : 100% !important;}
                                            ")),
    # waitress declare
    use_waitress(color = "#1976d2", percent_color = "#223344"),
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
                    HTML("<font color='black' size='4'>CellEnrich - Median</font>"),
                    HTML("<font color='black' size='4'>CellEnrich - Fast GSEA</font>")
                  ),
                  choiceValues = c("CellEnrich - median",
                                   "CellEnrich - FGSEA"),

                  selected = "CellEnrich - median",
                )
              ),
              material_card(
                radioButtons(
                  "plotOption",
                  label = HTML("<font color='black' size='5'>Scatter Plot</font>"),#"Scatter Plot",
                  choiceNames = list(
                    HTML("<font color='black' size='4'>PCA</font>"),
                    HTML("<font color='black' size='4'>TSNE</font>"),
                    HTML("<font color='black' size='4'>UMAP</font>")
                  ),
                  choiceValues = c("PCA", "TSNE", "UMAP"),
                  selected = "UMAP"
                ),

                material_number_box(
                  input_id = "topdims",
                  label = HTML("<font color='black' size='4'>Top-N dims</font>"), # top dims
                  min_value = 30,
                  max_value = 100,
                  initial_value = 50,
                  step_size = 10
                )
              )
              , width = 2 ## column width
            ),

            material_column(
              material_card(  title = HTML("<font color='black' size='4'> </font>"),
                material_number_box(
                  input_id = "medianCoefficient",
                  label = HTML("<font color='black' size='4'>Median coefficient</font>"),
                  min_value = 0,
                  max_value = 2,
                  initial_value = 1,
                  step_size = 0.05
                ),
                material_number_box(
                  input_id = "fgseaNsample",
                  label = HTML("<font color='black' size='4'>FGSEA - N(%) of Top-depth samples</font>"),
                  min_value = 0,
                  max_value = 100,
                  initial_value = 10,
                  step_size = 1
                )
              )
              , width = 3 ## column width
            ),
            material_column(
              material_card(
                radioButtons(
                  "genesetOption",
                  label=HTML("<font color='black' size='5'>Genesets</font>"),
                  choiceNames = list(
                    HTML("<font color='black'>Human-Reactome</font>"),
                    HTML("<font color='black'>Human-WikiPathway</font>"),
                    HTML("<font color='black'>Human-KEGG</font>"),
                    tags$span(style = "color:black", "Human-GOBP"),
                    tags$span(style = "color:black", "Human-GOCC"),
                    tags$span(style = "color:black", "Human-GOMF"),
                    tags$span(style = "color:black", "Mouse-Reactome"),
                    tags$span(style = "color:black", "Mouse-WikiPathway"),
                    tags$span(style = "color:black", "Mouse-KEGG"),
                    tags$span(style = "color:black", "Mouse-GOBP"),
                    tags$span(style = "color:black", "Mouse-GOCC"),
                    tags$span(style = "color:black", "Mouse-GOMF")
                  ),
                  choiceValues = default_genesets,
                  selected= default_genesets[1]
                ),
                sidebarPanel(
                  tags$style("
                               .btn-file {
                                  background-color:#1976d2;
                                  border-color: none;
                               }
                               "),
                  HTML("<font color='black' size='3'>User defined geneset</font>"),
                  fileInput("user_gs", "", placeholder = "RData format is required!"),
                )
              ),
              width = 3
            ),
            material_column(
              material_card( title = HTML("<font color='black' size='5'> </font>"),
                material_number_box(
                  input_id = "minGenesetSize",
                  label = HTML("<font color='black' size='4.5'>Minimum Geneset Size</font>"), #"Minimum Gene-set Size",
                  min_value = 10,
                  max_value = 30,
                  initial_value = 15,
                  step_size = 5
                ),
                material_number_box(
                  input_id = "maxGenesetSize",
                  label = HTML("<font color='black' size='4.5'>Maximum Geneset Size</font>"), #"Maximum Gene-set Size",
                  min_value = 250,
                  max_value = 750,
                  initial_value = 500,
                  step_size = 5
                ),
                material_number_box(
                  input_id = "pwFrequency",
                  label = HTML("<font color='black' size='4.5'>Pathway Frequency</font>"),
                  min_value = 0,
                  max_value = 0.5,
                  initial_value = 0.1,
                  step_size = 0.05
                ),
                material_number_box(
                  input_id = "qvalueCutoff",
                  label = HTML("<font color='black' size='4.5'>Q-value threshold</font>"),
                  min_value = 0,
                  max_value = 1,
                  initial_value = 0.05,
                  step_size = 0.01
                )
              ),
              width = 2
            )
          ),
          shiny::tags$div(class = "runbutton",
                          solvedButton(
                            inputId = "StartCellEnrich",
                            label = "RUN",
                            style = "background-color: #1976d2; height:60px; width:360px; font-size : 32px;",
                            onClick = 'console.log("CellEnrich");'
                          )
          )
          ,
          depth = 3
        ),
        width = 12,
      )
    ),

    # tSNE/UMAP plot
    material_row(
      material_column(

        material_card(
          title = shiny::tags$h4("Scatter & Bar"), depth = 3,
          material_row(
            shiny::downloadButton("sctdn", "Scatter", style = "background-color : #616161 !important")
          ),
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
              material_button("sigbtn", "Odds Ratio", icon = "sort", color = "blue darken-2"),
              material_button("freqbtn", "Frequency", icon = "menu", color = "blue darken-2")
          ),

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
              width = 2
            ),
            material_column(
              shiny::downloadButton("imgdn", "Pathway Significance in each Group", style = "background-color : #616161 !important")
              , width = 4
            ),
            material_column(
              width = 2
            ),
            material_column(
              shiny::downloadButton("cmpdn", "Pathway Significance in whole data", style = "background-color : #616161 !important")
              , width = 4
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
              DT::DTOutput("userpw"),
              material_row(
                material_button("Emphasize", "Plot the Selected Pathways", color = "blue darken-2"),
                material_button("ClearList", "Clear List", color = "blue darken-2")
              ),
              material_row(
                  shiny::downloadButton("sppcdn", "Cell-Pathway (P-values)", style = "background-color : #616161 !important"),
                  shiny::downloadButton("tbldn", " All Significant Pathways", style = "background-color : #616161 !important")
              )
            ),
            tags$head(
              tags$style(HTML(custom_css))
            ),
            title = NULL,
            material_tabs(tabs = setNames(tab_ids, Cells), color="#223344"),
            lapply(1:length(Cells), function(i) {
              material_tab_content(
                tab_id = tab_ids[i],
                shiny::tags$div(
                  class = paste0("card z-depth-5 color: null"),
                  style = paste0("border : solid 0.5em ", CardColors[i]), # border color defined
                  shiny::tags$div(
                    class = "card-content",
                    shiny::tags$span(class = "card-title", Cells[i]), # title
                    shiny::tags$div(class = "divider"), # divider = TRUE
                    DT::dataTableOutput(
                      paste0("dynamicGroupTable", i),
                      width = "100%",
                      height = "100%"
                    )
                  )
                )
              )
            })
          ),

          material_card(
            title = shiny::tags$h4("Gene Activity (GA) Map"), divider = TRUE,
            material_row(
              material_column(
                plotOutput("geneAct", height = "1000px")
                , width = 12
              )
            ),
            material_row(
              material_column(
                width = 3
              ),
              material_column(
                shiny::downloadButton("actScore", "Save GA Map", style = "background-color : #616161 !important"),
                shiny::downloadButton("actScoredata", "Save GA Data", style = "background-color : #616161 !important")
                , width = 6
              ),
              material_column(
                 width = 3
              )
            )
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
                numericInput("biCount", label = "Top Pathways in each Group", value = 3, min = 1, max = 10, step = 1),
                numericInput("biFont", label = "Cell Type Label", value = 7, min = 1, max = 10, step = 1),
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
        ),
        material_row(
          shiny::downloadButton("markerdownload", "Markers of Groups", style = "background-color : #616161 !important")
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
                            panel.border = element_rect(colour = 'black', fill=NA, linewidth=0.25),
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
        # cell_pval <- c(cell_pval, max(tmp))

        # find the best Pathway
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
                              panel.border = element_rect(colour = 'black', fill=NA, linewidth=0.25),
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
#' @param genesets (optional), user Pathway to analysis. [list]
#'
#' @return no return.
#'
#'
#' @importFrom DT dataTableOutput
#' @importFrom Matrix t
#'
#' @rawNamespace import(SingleCellExperiment, except = show)
#' @import tidyr
#' @import purrr
#' @import magrittr
#' @import scales
#' @import plyr
#' @import dplyr
#' @rawNamespace import(shiny, except = dataTableOutput)
#' @import shiny
#' @import shinymaterial
#' @import ggplot2
#' @import uwot
#' @import htmltools
#' @import ggbiplot
#' @import waiter
#' @import Seurat
#' @import farver
#' @rawNamespace import(shinyjs, except = runExample)
#' @import sortable
#' @import scran
#' @import ggrepel
#' @import fgsea
#' @import shinyFeedback
#'
#' @export

CellEnrich <- function(CountData, GroupInfo, genesets = NULL, use.browser=TRUE) {
  library(tidyr)
  library(dplyr)
  library(shiny)

  if(!require(ggbiplot)){
    remotes::install_github('vqv/ggbiplot')
  }

  library(ggbiplot)
  library(ggrepel)
  options(useFancyQuotes = FALSE)

  if (is.factor(GroupInfo)){
    GroupInfo <- unfactor(GroupInfo)
  }

  server <- function(input, output, session) {

    toptab <- function (genesets, TOPN = 3, OddsRatio = TRUE, myplot='biplot'){
      Cells <- sort(unique(GroupInfo))
      # pres : which gene-sets are significant for each cells.
      # pres2 : for each gene-sets, how many cells are significant that gene-sets.

      if (OddsRatio) { ## OddsRatio
        total <- length(GroupInfo)

        dat <- OR %>%
          group_by(Cell) %>%
          arrange(Cell) %>%
          top_n(TOPN)

        gs <- unique(dat$Pathway)

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

              if (is.na(N - A) || (N==A) ){
                return (0)
              }

              if (is.na(K - B) || (K==B) ){
                return (0)
              }

              return( (B/(K - B)) / (A/(N - A)) )
            })
          ), 4)
          # Cell, Pathway, OR
        }

        # adjusted to original OR values
        for (i in 1:length(Cells)){
          corGeneset <- dat$Pathway[dat$Cell == Cells[i]]
          corOddRatio <- dat$OddsRatio[dat$Cell == Cells[i]]
          for (g in corGeneset){
            tab[g, Cells[i]] <- corOddRatio[corGeneset==g]
          }
        }

        ## log scale for OR
        tab <- log1p(tab)
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

      # saveRDS(tab, "ortab.rds")
      return(tab)
    }

    buildbiplot <- function(biFont, biX, biY, genesets, TOPN = 5, OddsRatio = TRUE, gsFont=5, axtxt=13, axlab=15,
                            myplot='biplot', tab=matrix(0, 0, 0)) {

      hmtype <- 'Odds Ratio based Biplot'
      if (!OddsRatio){
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

    buildheatplot <- function(biFont, biX, biY, genesets, TOPN = 5, OddsRatio = TRUE, gsFont=5, axtxt=13, axlab=15,
                              myplot='biplot', tab=maxtrix(0, 0, 0)) {

      ###########   Heatmap
      hmtype <- 'Odds Ratio based Heatmap'
      if (!OddsRatio){
        hmtype <- 'Frequency based Heatmap'
      }

      # Define breaks
      mat_data <- round(tab, 2)
      col_breaks <- c(seq(0, max(mat_data), length.out=101))
      my_palette <- c(colorRampPalette(c("white", "red"))(length(col_breaks)-1))


      library(gplots)
      HeatPlot <<- heatmap.2(mat_data,
                             main = hmtype,
                             density.info="none",
                             trace="none",
                             margins =c(15,70),
                             col=my_palette,
                             scale = "none",
                             breaks=col_breaks,
                             dendrogram="row",
                             Colv=NA,
                             symkey=F,

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
      return(HeatPlot)
    }

    buildGeneActivity <- function(selected_pathway, group, genesets, logtab) {

      if (!is.null(selected_pathway) && !is.null(group)){

        memGenes <- genesets[selected_pathway]
        memGenes <- unlist(memGenes)
        memGenes <- memGenes[memGenes %in% rownames(logtab)]

        ### fgsea??
        # if (!is.null(rdf$fgsea_sample_idx)){
        #   fullpw_data <- logtab[memGenes, rdf$fgsea_sample_idx]
        #   fullpw_data <- fullpw_data[ , which(colnames(fullpw_data)==group)]
        #
        # } else { # NULL
          fullpw_data <- logtab[memGenes, which(colnames(logtab)==group)]
        # }

          saveRDS(fullpw_data, "pw_data.rds")

        ### sorting genes
        avga <- rowMeans(fullpw_data)
        memGenes <- memGenes[order(avga, decreasing = T)]
        fullpw_data <- fullpw_data[memGenes, ] ## sorted decreasingly
        fullpw_data <<- fullpw_data

        if (length(memGenes) > 60){  ## trimming
          memGenes <- memGenes[1:60]
        }
        mat_data <- fullpw_data[memGenes, ]

        saveRDS(mat_data, "mat_data.rds")

        mat_data <- round(mat_data, 2)
        col_breaks <- c(seq(min(mat_data), max(mat_data), length.out=51))
        my_palette <- c(colorRampPalette(c("blue2", "white", "red2"))(length(col_breaks)-1))

        library(gplots)
        actMap <- heatmap.2(mat_data,
                               main = paste0(selected_pathway, " @", group),
                               density.info="none",
                               trace="none",
                               margins =c(1,10),
                               col=my_palette,
                               scale = "none",
                               breaks=col_breaks,
                               dendrogram="col",
                               Colv=TRUE,
                               Rowv = FALSE,
                               labRow = NULL,
                               labCol = NA,
                               symkey=T,

                               # additional control of the presentation
                               lhei = c(2, 15),       # adapt the relative areas devoted to the matrix
                               lwid = c(2, 20),
                               cexRow = 1.3,
                               cexCol = 1,
                               key.title = NA,
                               key.xlab = NA,
                               key.ylab = NA,
                               key.par = list(mar=c(4, 1, 2, 1),
                                              mgp=c(1.5, 0.5, 0),
                                              cex=1)
        )
        return(actMap)
      } else {
        return(NULL)
      }
    }

    glob_options <- reactiveValues(usergs = NULL,
                                   reRun = FALSE
    )
    observeEvent(input$user_gs,{
      other_geneset <- input$user_gs
      ext <- tools::file_ext(other_geneset$datapath)

      req(other_geneset)
      glob_options$usergs  <- other_geneset$datapath

      otherVal <- other_geneset$name
      updatedValues <- c(default_genesets, otherVal)

      updateRadioButtons(session, "genesetOption", choices = updatedValues, selected= otherVal)
      print(glob_options$usergs )
      cat("Geneset added!")
    })

    ### CODES
    # variable initialize

    dtobj <- dfobj <- pres <- pres2 <- presTab <- logtab <- ""
    CellPathwayDF <- ""
    gt <- Cells <- A <- ""
    CellScatter <- ""
    CellHistogram <- ""
    BiPlot <- HeatPlot <- OR <- ""


    # ##### trigger value
    observeEvent(input$StartCellEnrich, {
      pt <- proc.time()

      # ------ Hide Start Button
      shinyjs::disable("StartCellEnrich")

      #Disable Emphasize 1st
      shinyjs::runjs('$("#Emphasize").attr("disabled",true)')
      shinyjs::runjs('$("#Emphasize_freq").attr("disabled",true)')


      # ------ Load Genesets
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

      if (glob_options$reRun){
        if (!is.null(glob_options$user_gs)) {
          genesets <- read_Input_geneset(glob_options$user_gs)
          if (is.null(genesets)){
            shiny::showNotification("Geneset file is invalid!", type = "error", duration = 20)
            return(NULL)
            # ------ Enable Start Button Again
            shinyjs::enable("StartCellEnrich")
          }
        }
      }
      # Check that data object exists and is data frame.
      if (is.null(genesets)) {
        gs_file <- switch( input$genesetOption,
                           "Human-Reactome" = "Human_Reactome.RData",
                           "Human-WikiPathway" ="Human_WikiPathways.RData",
                           "Human-KEGG" ="Human_KEGG.RData",
                           "Human-GOBP" ="Human_GOBP.RData",
                           "Human-GOCC" ="Human_GOCC.RData",
                           "Human-GOMF" ="Human_GOMF.RData",
                           "Mouse-Reactome" ="Mouse_Reactome.RData",
                           "Mouse-WikiPathway" ="Mouse_WikiPathways.RData",
                           "Mouse-KEGG" ="Mouse_KEGG.RData",
                           "Mouse-GOBP" ="Mouse_GOBP.RData",
                           "Mouse-GOCC" ="Mouse_GOCC.RData",
                           "Mouse-GOMF" ="Mouse_GOMF.RData", "Human_Reactome.RData")
        load(gs_file)
      }

      if (is.null(genesets)) {
        genesets <- read_Input_geneset(glob_options$usergs)

        if (is.null(genesets)){
          shiny::showNotification("Geneset file is invalid!", type = "error", duration = 20)
          return(NULL)
          # ------ Enable Start Button Again
          shinyjs::enable("StartCellEnrich")
        }
      }
      if (is.null(genesets)) {
        shiny::showNotification("Geneset file is missing!", type = "error", duration = 20)
        return(NULL)
        # ------ Enable Start Button Again
        shinyjs::enable("StartCellEnrich")
      }

      ## -- remove non-expressed genes ####
      rs <- rowSums(CountData)
      non_exped_genes <- rownames(CountData)[rs==0]
      CountData <- CountData[!(rownames(CountData) %in% names(non_exped_genes)),]

      # CountData <- pbmc
      genes <- rownames(CountData)
      genesets <- GenesetFlush(genes, genesets)
      lgs <- getlgs(genesets)

      # ------ Genesetsize Flush
      # genesets <- GenesetsizeFlush(genesets, lgs, 15, 500)
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
        # ------ Enable Start Button Again
        shinyjs::enable("StartCellEnrich")
        # stop("This dataset has not enough significant genes!")
      }

      if (length(genesets)==0) {
        shiny::showNotification("This dataset has no significant pathway detected based on current genesets!", type = "error", duration = 60)
        return(NULL)
        # ------ Enable Start Button Again
        shinyjs::enable("StartCellEnrich")
        # stop("This dataset has no significant pathway detected!")
      }
      glob_options$reRun <- TRUE

      # ------ for test
      q0 <- input$qvalueCutoff

      # ------ Create new Waitress
      w <- Waitress$new(selector = NULL, theme = "overlay-radius")

      # ------ Background genes
      A <<- getBackgroundGenes(genesets)

      # ------ Calculate TSNE / UMAP First
      # seu <- getTU(CountData, GroupInfo, "CellEnrich - FGSEA", topdims=50)
      seu <- getTU(CountData, GroupInfo, input$plotOption, topdims=input$topdims)

      #PCA
      if (input$plotOption == "PCA") {
        ori_dfobj <- data.frame(Embeddings(seu, 'pca')[,1:2], col = GroupInfo, stringsAsFactors = FALSE)
      }
      # TSNE
      if (input$plotOption == "TSNE") {
        ori_dfobj <- data.frame(Embeddings(seu, 'tsne'), col = GroupInfo, stringsAsFactors = FALSE)
      }
      # UMAP
      if (input$plotOption == "UMAP") {
        ori_dfobj <- data.frame(Embeddings(seu, 'umap'), col = GroupInfo, stringsAsFactors = FALSE)
      }

      # ----- FGSEA ####
      if (input$FCoption == "CellEnrich - FGSEA") {
        library(fgsea)
        library(purrr)
        library(tidyr)

        ###  genes & cells  matching
        scaleCount <- as.matrix(seu@assays$RNA$scale.data)
        rownames(scaleCount) <- rownames(seu)
        scaleCount <- scaleCount[rownames(scaleCount) %in% unique(unlist(genesets)),]
        colnames(scaleCount) <- colnames(seu)

        ### sampling cells if N-cell > Nmax
        Nmax <- round(ncol(CountData) / 100 * input$fgseaNsample)
        # Nmax <- 50

        if (is.null(Nmax)){
          shiny::showNotification("At least 50 samples are required to estimate via FGSEA. 50 samples are used.", type = "error", duration = 30)
          Nmax = 50
        }

        if (Nmax < 50){
          shiny::showNotification("At least 50 samples are required to estimate via FGSEA. 50 samples are used.", type = "error", duration = 30)
          Nmax = 50
        }

        if (Nmax > ncol(scaleCount)){
          shiny::showNotification("Cannot sample more than available cells. All cells are used.", type = "error", duration = 30)
          Nmax <- ncol(scaleCount)
        }
        sample_idx <- NULL

        if (ncol(scaleCount) > Nmax){
          n_group <- unique(GroupInfo)
          sample_per_group <- Nmax %/% length(n_group)

          sample_idx <- lapply(n_group, function(g){
            id_g <- which(GroupInfo==g)

            if(length(id_g) > sample_per_group){
              # id_g <- sample(id_g, sample_per_group, replace = F)
              group_lib <- colSums(scaleCount[, id_g])
              group_lib <- sort(group_lib, decreasing = T)
              top <- names(group_lib[1:sample_per_group])
              return(which(colnames(scaleCount) %in% top))
            }
            return(id_g)
          })
          sample_idx <- unlist(sample_idx)

          ## add missing slots
          if (length(sample_idx) < Nmax){
            x <- 1:ncol(scaleCount)
            x <- setdiff(x, sample_idx)

            group_lib <- colSums(scaleCount[, x])
            group_lib <- sort(group_lib, decreasing = T)
            top <- names(group_lib[1:(Nmax-length(sample_idx))])
            more_idx <- which(colnames(scaleCount) %in% top)

            sample_idx <- c(sample_idx, more_idx)
          }
          sample_idx <- sort(sample_idx)
          scaleCount <- scaleCount[, sample_idx]

          rdf$fgsea_sample_idx <- sample_idx
        }

        conditions <- colnames(scaleCount)
        names(conditions) <- conditions
        permtimes = 1000

        #
        library(tictoc)
        tic("sleeping")

        finalres <- map_dfr(.x = conditions, .f = ~ {
          stats <- scaleCount[, .x]
          options <- list(
            pathways = genesets,
            stats = stats,
            nPermSimple = permtimes,
            nproc = 8,
            scoreType='std',
            eps=1e-16
          )
          withr::with_seed(seed <- sample.int(.Machine$integer.max, 1L), {
            result <- suppressWarnings(do.call(what = fgsea::fgsea, args = options))

            result$leadingEdge = sapply(seq_len(nrow(result)), function(x)paste0(result$leadingEdge[x][[1]], collapse = ', '))

            tres <- result %>% filter(padj < q0) ### CRITICAL!
            if (nrow(tres) == 0){
              tres <- result %>% arrange(pval) %>% top_n(-1, pval)
            }

            tmp <- unlist(tres$leadingEdge)
            tmp <- lapply(tmp, function(x){
              strsplit(x, ", ")[[1]]
            })
            tmp <- unlist(tmp)
            ct <- data.frame(table(tmp))

            res <- ct[ct$Freq >= 1,] #quantile(ct$Freq)["0%"],]
            idx <- paste0(sapply(res$tmp, function(x){which(rownames(CountData)==x)}), collapse = ', ')
          })
        }, .id = "condition")

        saveRDS(finalres, "finalres3.rds")

        s <- lapply(finalres, function(x){
          tmp <- unlist(x)
          tmp <- lapply(tmp, function(x){
            strsplit(x, ", ")[[1]]
          })
          tmp <- as.numeric(unlist(tmp))
        })

        saveRDS(s, "s.rds")

        mytoc <- toc()
        print(mytoc)
        # saveRDS(mytoc, 'mytoc.rds')

        # subseting samples ####
        CountData <- CountData[, sample_idx]
        GroupInfo <- GroupInfo[sample_idx]
        seu <- seu[, sample_idx]
        seu <<- seu

        # if (is.null(sample_idx)){
        names(s) <- GroupInfo

        # } else {
        #   mylist <- sapply(GroupInfo,function(x) {})
        #   mylist[sample_idx] <- s
        #
        #   for (i in 1:length(mylist)){
        #     # i = 1
        #     if (!i %in% sample_idx){    ## length == 0, empty
        #       ti <- which(GroupInfo[sample_idx] == GroupInfo[i])
        #
        #       best_next <- which(sample_idx[ti] > i)
        #
        #       if (length(best_next) > 0){
        #         best_next <- min(best_next)
        #
        #       } else {
        #         best_next <- length(ti) ## 240325
        #       }
        #       mylist[i] <- s[best_next]
        #     }
        #   }
        #   s <- mylist
        #   names(s) <- GroupInfo
        # }

      }
      else { ### other cases
        s <- findSigGenes(CountData, input$FCoption, GroupInfo, coef = input$medianCoefficient)
      }
      cat("s Finished\n")

      # apply logtab ####
      logtab <- apply(seu@assays$RNA$counts, 1, function(v){
        if (is.na(median(v[v > 0]))){
          v
        } else {
          v <- v - median(v[v > 0])
        }
      })
      logtab <- t(logtab)
      colnames(logtab) <- GroupInfo
      logtab <<- logtab

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

      # ------ Find Significant Genes with findMarkers ####

      s2 <- findSigGenesGroup(CountData, GroupInfo, q0, TopCutoff = 5)
      rc <- rownames(CountData)
      # ----- Free memory to calculate biobj ####
      rm(CountData)

      # ------ marker l1
      markerl1 <- s2 %>% filter(Top <= 10)
      markerl1$Group <- as.factor(markerl1$Group)
      colnames(markerl1)[4] <- "FDR"
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

      # saveRDS(s2, 'sigGeneGroup.rds')
      output$markerdownload <- downloadHandler(
        filename = "markers_of_groups.csv",
        content = function(file) {
          write.csv(s2, file, row.names = FALSE)
        }
      )

      ############################################################################################3
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

          tmp_pv[which(tmp_pv < 1e-20)] <- 1e-20

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
              OddsRatio = unname(ors[pathways])
            )
          )
        }
      }

      # ------ Hypergeometric pvalue calculation
      lgs <- getlgs(genesets)
      lens <- length(s)
      lens100 <- round(lens / 100)

      cat("pres declare\n")
      pres <- list()
      presTab_pval <- presTab <- c()

      if (lens >= 100) {
        w$start()
        for (i in 1:lens) {
          if (i %% lens100 == 0) w$inc(1)
          prespv <- getHyperPvalue(rc[s[[i]]], genesets, A, lgs, q0, biobj)
          pres[[i]] <- which(p.adjust(prespv, "fdr") < q0)
          presTab_pval <- cbind(presTab_pval, p.adjust(prespv, "fdr"))
        }
        w$close()
        colnames(presTab_pval) <- colnames(seu)
      }
      else {
        for (i in 1:lens) {
          prespv <- getHyperPvalue(rc[s[[i]]], genesets, A, lgs, q0, biobj)
          pres[[i]] <- which(p.adjust(prespv, "fdr") < q0)
          presTab_pval <- cbind(presTab_pval, p.adjust(prespv, "fdr"))
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
      presTab <- presTab_pval
      presTab[which(presTab < 1e-20)] <- 1e-20
      presTab <- round(-log10(presTab), 4) ### presTab -log10(p-value) -> Q-value
      presTab <<- presTab

      pres <<- pres
      # saveRDS(pres, "pres.rds")

      # pres : which gene-sets are significant for each cells.
      # pres2 : for each gene-sets, how many cells are significant that gene-sets.
      cat("pres2\n")
      pres2 <- sort(table(unlist(pres)), decreasing = T)
      if (length(s) != 0) {
        names(pres2) <- names(genesets)[as.numeric(names(pres2))]
      }
      pres2 <<- pres2
      # saveRDS(pres2, "pres2.rds")

      # ------ CellPathwayDF
      CellPathwayDF <- buildCellPathwayDF(GroupInfo, pres, genesets, input$pwFrequency)
      # saveRDS(genesets, "filtered_gs.rds")

      PP <- pathwayPvalue(GroupInfo, pres, pres2, genesets)

      # OR -> # Group / PATHWAY / OddsRatio
      if (nrow(tmp_df) > 0) {
        OR <- tmp_df
        colnames(OR) <- c("Cell", "Pathway", "OddsRatio")
        OR <<- OR
      }
      else {
        OR <<- getOddRatio(GroupInfo, pres, pres2, genesets, input$pwFrequency)
      }

      CellPathwayDFP <- CellPathwayDF %>%
        inner_join(PP) %>% inner_join(OR)  %>% select(-Pvalue)  # %>% filter(Qvalue > 1)

      CellPathwayDFP <- CellPathwayDFP %>% filter(Qvalue > -log10(q0))
      output$tbldn <- downloadHandler(
        filename = "all_sig_pathways.csv",
        content = function(file) {
          write.csv(CellPathwayDFP, file, row.names = FALSE)
        }
      )

      OR <<- OR %>% filter(OddsRatio > 1)
      CellPathwayDF <- CellPathwayDFP %>% filter(OddsRatio > 1)
      CellPathwayDF <<- CellPathwayDF

      # l2
      CellMarkers <- data.frame()

      Cells <- sort(unique(GroupInfo))
      Cells <<- Cells

      for (i in 1:length(Cells)) {
        # thisCell <- Cells[i]
        thisCellPathways <- CellPathwayDF %>%
          filter(Cell == Cells[i]) %>%
          dplyr::select(Pathway)

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

      CellScatter <<- getCellPlot(ori_dfobj, Cells)  ## updated
      # CellScatter <<- getCellPlot(dfobj, Cells)

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
        CardColors <- briterhex(scales::hue_pal(h = c(20, 350), c = 100, l = 65, h.start = 0,
                                                direction = 1)(length(Cells)))
        options(useFancyQuotes = FALSE)
        tagList(
          material_row(
            lapply(1:length(Cells), function(i) {
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
                    )

                  )
                ),
                width = "100%"
              )
            })
          )
        )
      })


      # ------ fill dynamic table
      # ordered by Count , not length;
      lapply(1:length(Cells), function(i) {
        output[[paste0("dynamicGroupTable", i)]] <- DT::renderDataTable(
          DT::datatable(CellPathwayDF[which(CellPathwayDF$Cell==Cells[i]), -1],
                        colnames = c("-log10 Qvalue" = "Qvalue", " "="N_cell"),
                        rownames=FALSE,
                        selection = 'single',
                        options=list(dom='ltp',
                                     scroller = TRUE,
                                     scrollX = TRUE,
                                     autoWidth = TRUE,
                                     lengthChange = FALSE,
                                     order = list(list(5,'desc'))) # odds ratio based
          ))
      })

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
              thisCellData <- thisCellData %>% top_n(1, wt = OddsRatio)
              if (nrow(thisCellData) >= 1) {
                thisCellData <- thisCellData %>% top_n(1)
              }
            }
          }
          res <- c(res, paste0(thisCellData$Pathway, " @", thisCellData$Cell))
        }
      }

      output$legendTable <- DT::renderDataTable(
        buildGradientLegend(res, Cells = Cells)
      )

      plotImg <- emphasize(FALSE, res, dfobj, Cells, pres, genesets, seu, presTab, maptitle = 'Frequency')
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
          # buildGradientLegend(res, Cells = Cells)
          buildGradientLegend(res, img = TRUE, name = file, Cells = Cells)

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

      for (i in 1:length(Cells)) {

        thisCellData <- CellPathwayDF %>% dplyr::filter(Cell == Cells[i])

        if (nrow(thisCellData) >= 1) {
          thisCellData <- thisCellData %>% top_n(1, wt = OddsRatio)

          if (nrow(thisCellData) >= 1) {
            thisCellData <- thisCellData %>% top_n(-1, wt = Size)
            if (nrow(thisCellData) >= 1) {
              thisCellData <- thisCellData %>% top_n(1, wt = Frequency)
              if (nrow(thisCellData) >= 1) {
                thisCellData <- thisCellData %>% top_n(1)
              }
            }
          }
          res <- c(res, paste0(thisCellData$Pathway, " @", thisCellData$Cell))
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


    Cells <- sort(unique(GroupInfo))

    # rdf ####
    rdf <- reactiveValues(df = data.frame(Pathway=character(), Group=character()),
                          current_pw=NULL,
                          current_group=NULL,
                          fgsea_sample_idx=NULL)

    lapply(1:length(Cells), function(i) {
      observeEvent(input[[paste0("dynamicGroupTable", i, "_cell_clicked")]],{
        info = input[[paste0("dynamicGroupTable", i, "_cell_clicked")]]

        if (is.null(info$value) || info$col != 0) return()

        rdf$current_pw = info$value
        rdf$current_group = Cells[i]

        # replace available pw with the latest selection
        exist_ct <- which(rdf$df$Group == Cells[i])

        if (length(exist_ct) == 1){
          if (rdf$df$Pathway[exist_ct] == info$value){ # de-select a pathway
            rdf$df <- rdf$df[-exist_ct, ]
          } else {
            rdf$df$Pathway[exist_ct] <- info$value
          }

        } else {
          rdf$df <- rbind(rdf$df, data.frame(Pathway=info$value, Group=Cells[i]))
        }

        if (nrow(rdf$df) > 0){
          shinyjs::runjs('$("#Emphasize").attr("disabled",false)')
          shinyjs::runjs('$("#Emphasize_freq").attr("disabled",false)')
        } else {
          shinyjs::runjs('$("#Emphasize").attr("disabled",true)')
          shinyjs::runjs('$("#Emphasize_freq").attr("disabled",true)')
        }

        #---- Gene Activity triggered ####
        output$geneAct <- renderPlot(buildGeneActivity(rdf$current_pw, rdf$current_group, genesets, logtab))

        output$actScoredata <-  downloadHandler(
          filename = "GA_data.csv",
          content = function(file) {
            # write.csv(actMap$carpet, file, row.names = FALSE)
            write.csv(fullpw_data, file, row.names = FALSE)
          }
        )

        output$actScore <- downloadHandler(
          filename = function() {
            "GA_map.pdf"
          },
          content = function(file) {

            revRowInd <- match(c(1:length(actMap$rowInd)), actMap$rowInd)
            revColInd <- match(c(1:length(actMap$colInd)), actMap$colInd)

            pdf(file, width = 24, height = 16)
            heatmap.2(t(actMap$carpet)[revRowInd, revColInd],
                      main = paste0(rdf$current_pw, " @", rdf$current_group),
                      density.info="none",
                      trace="none",
                      margins =c(2,20),
                      col=actMap$col,
                      scale = "none",
                      breaks=actMap$breaks,
                      dendrogram = "col",
                      Colv=TRUE,
                      Rowv = FALSE,
                      labRow = NULL,
                      labCol = NA,
                      symkey=T,

                      # additional control of the presentation
                      lhei = c(2, 15),       # adapt the relative areas devoted to the matrix
                      lwid = c(2, 20),
                      cexRow = 1.5,
                      cexCol = 2,
                      key.title = NA,
                      key.xlab = NA,
                      key.ylab = NA,
                      key.par = list(mar=c(4, 1, 2, 1),
                                     mgp=c(1.5, 0.5, 0),
                                     cex=1)
            )

            dev.off()
          }
        )


      })
    })

    output$userpw <- DT::renderDT(
      rdf$df[nrow(rdf$df):1,], options = list(dom = 'tp',
                                              scroller = TRUE,
                                              scrollX = TRUE,
                                              autoWidth = TRUE,
                                              pageLength = 10
      ),
      rownames=FALSE,
      selection='none'
    )

    observeEvent(input$Emphasize, {
      if (input$Emphasize == 0) { # prevent default click state
        return(NULL)
      }

      res <- c()
      for (i in 1:nrow(rdf$df)){
        res <- c(res, paste0(rdf$df[i, 1], " @", rdf$df[i, 2]))
      }

      shinyjs::show("legenddn")
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

    observeEvent(input$Emphasize_freq, {
      if (input$Emphasize_freq == 0) { # prevent default click state
        return(NULL)
      }

      res <- c()
      for (i in 1:nrow(rdf$df)){
        res <- c(res, paste0(rdf$df[i, 1], " @", rdf$df[i, 2]))
      }

      shinyjs::show("legenddn")
      plotImg <- emphasize(FALSE, res, dfobj, Cells, pres, genesets, seu, presTab, maptitle = 'Frequency')
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
          buildLegend(res, img = TRUE, name = file)
        }
      )
    })

    # clear list in Cell tab
    observeEvent(input$ClearList, {
      if (input$ClearList == 0) { # prevent default click state
        return(NULL)
      }
      rdf$df <- rdf$df[0,]
      shinyjs::runjs('$("#Emphasize").attr("disabled",true)')
      shinyjs::runjs('$("#Emphasize_freq").attr("disabled",true)')
    })

    observeEvent(input$freqbp, {
      if (input$freqbp == 0) {
        return(NULL)
      }

      output$biPlot <- renderPlot(buildbiplot(input$biFont, input$biX, input$biY, genesets, TOPN = input$biCount,
                                              gsFont = input$gsFont, axtxt = input$axtxt, axlab = input$axlab,
                                              OddsRatio = FALSE,
                                              tab=toptab(genesets, TOPN = input$biCount, OddsRatio = FALSE)))

      output$heatPlot <- renderPlot(buildheatplot(input$biFont, input$biX, input$biY, genesets, TOPN = input$biCount,
                                                  gsFont = input$gsFont, axtxt = input$axtxt, axlab = input$axlab,
                                                  OddsRatio = FALSE,
                                                  tab=toptab(genesets, TOPN = input$biCount, OddsRatio = FALSE)))

      output$biplotdn <- downloadHandler(
        filename = function() {
          "mybiplot.pdf"
        },
        content = function(file) {
          ggsave(file, device = "pdf", width = 24, height = 16)
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
                    density.info="none",
                    trace="none",
                    margins =c(15,70),
                    col=HeatPlot$col,
                    scale = "none",
                    breaks=HeatPlot$breaks,
                    dendrogram = "row",
                    Rowv = HeatPlot$rowDendrogram,
                    Colv=NA,
                    symkey=F,

                    # # additional control of the presentation
                    lhei = c(2, 13),
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
                                              OddsRatio = TRUE,
                                              tab = toptab(genesets, TOPN = input$biCount, OddsRatio = TRUE)))

      output$heatPlot <- renderPlot(buildheatplot(input$biFont, input$biX, input$biY, genesets, TOPN = input$biCount,
                                                  gsFont = input$gsFont, axtxt = input$axtxt, axlab = input$axlab,
                                                  OddsRatio = TRUE,
                                                  tab=toptab(genesets, TOPN = input$biCount, OddsRatio = TRUE)))

      output$biplotdn <- downloadHandler(
        filename = function() {
          "mybiplot.pdf"
        },
        content = function(file) {
          ggsave(file, device = "pdf", width = 24, height = 16)
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
                    density.info="none",
                    trace="none",
                    margins =c(15,70),
                    col=HeatPlot$col,
                    scale = "none",
                    breaks=HeatPlot$breaks,
                    dendrogram = "row",
                    Rowv = HeatPlot$rowDendrogram,
                    Colv=NA,
                    symkey=F,

                    # # additional control of the presentation
                    lhei = c(2, 13),
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

  ui <- CellEnrichUI(GroupInfo)


  if (use.browser){
    shiny::shinyApp(ui, server, options = list(launch.browser = TRUE))
  } else {
    shiny::shinyApp(ui, server)
  }


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
