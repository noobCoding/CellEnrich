# JH

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

gnm <- function(v) {
  out <- scMerge:::gammaNormMix(as.matrix(v), plot = FALSE)
  mat_prob <- matrix(out$probExpressed, nrow(v), ncol(v))
  mat_discretised <- 1 * (mat_prob > 0.5)
  return(mat_discretised)
}

findSigGenes <- function(v, method = "CellEnrich - median", Name) {
  if (!method %in% c("CellEnrich - median", "CellEnrich - mixture", "Fisher")) stop("wrong method")
  # it's already matrix
  # v <- as.matrix(v)
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
  require(scran)
  
  if (is.null(Count)) stop("Count must given")
  if (is.null(ClustInfo)) stop("ClustInfo must given")
  
  GrpRes <- scran::findMarkers(x = as.matrix(Count), ClustInfo, test = "wilcox", direction = "up")
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
#' @import shiny
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
#' @import ggrepel
#' @import shinyFeedback
#'
#' @export

CellEnrich <- function(CountData, GroupInfo, genesets = NULL) {
  require(dplyr)
  require(shiny)
  
  if(!require(ggbiplot)){
    remotes::install_github('vqv/ggbiplot')
  }
  
  require(ggbiplot)
  require(ggrepel)
  options(useFancyQuotes = FALSE)
  
  server <- function(input, output, session) {
     ### CODES
    
    # variable initialize
    
    dtobj <- dfobj <- pres <- pres2 <- ""
    CellPathwayDF <- ""
    gt <- Cells <- A <- ""
    CellScatter <- ""
    CellHistogram <- ""
    BiPlot <- OR <- ""
    
    observeEvent(input$StartCellEnrich, {
      pt <- proc.time()
    
      genesets <<- genesets
      # ------ for test
      # q0 <- 0.05
      
      q0 <- input$qvalueCutoff
      
      genes <- rownames(CountData)
      genesets <- GenesetFlush(genes, genesets)
      lgs <- getlgs(genesets)
      
      genesets <- GenesetsizeFlush(genesets, lgs, input$minGenesetSize, input$maxGenesetSize)
      
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
      
      s2 <- findSigGenesGroup(CountData, GroupInfo, q0, TopCutoff = 5)
      
      rc <- rownames(CountData)
      
      # ------ free memory to calculate biobj
      rm(CountData)
      
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
          sigidx <- which(p.adjust(tmp_pv, "fdr") < q0)
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
          
          pres[[i]] <- which(p.adjust(prespv, "fdr") < q0)
          
          prespv[which(prespv < 1e-20)] <- 1e-20
          
          presTab <- cbind(presTab, -log10(prespv))
        }
        w$close()
        colnames(presTab) <- colnames(CountData)
      }
      else {
        for (i in 1:lens) {
          prespv <- getHyperPvalue(rc[s[[i]]], genesets, A, lgs, q0, biobj)
          
          pres[[i]] <- which(p.adjust(prespv, "fdr") < q0)
          
          prespv[which(prespv < 1e-20)] <- 1e-20
          
          presTab <- cbind(presTab, -log10(prespv))
        }
        colnames(presTab) <- names(s)
      }
      
      cat("pres defined\n")
      rownames(presTab) <- names(genesets)
      
      pres <<- pres
      
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
      
      ggs <- unique(CellPathwayDFP %>% dplyr::select(Geneset))[, 1]
      ces <- sort(unique(CellPathwayDFP %>% dplyr::select(Cell))[, 1])
      
      nr <- length(ggs) # nrow
      nc <- length(ces) # ncol
      
      output$tbldn <- downloadHandler(
        filename = "mytable.csv",
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
          write.csv(CellPathwayDF, file, row.names = FALSE)
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
          dplyr::select(Geneset)
        
        # s2 <- findSigGenesGroup(CountData, GroupInfo, q0, TopCutoff = 5)
        # find markers
        
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
        additive <- additive %>% arrange(desc(Count))
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
      
      
      # group 별 significant pathways
      # group 별 DE Genes
      
      # is counted
      dtobj <<- buildDT(pres2)
      
    }
  }
    
      
    