getOddRatio = function(GroupInfo, pres, pres2, genesets){
  cat('getOddRatio\n')

  # pres : which gene-sets are significant for each cells.
  # pres2 : for each gene-sets, how many cells are significant that gene-sets.

  # 전체 그룹에서 유의한 회수 20 # pres2[genesets[i]]
  # 특정 그룹에서 유의한 회수 6 # pres2[thiscellidx]

  # 전체 그룹 Cell 수 : N
  # 특정 그룹 Cell 수 : K

  # Group_specific_OR = (6/K) / (14/N)

  res <- data.frame(stringsAsFactors = FALSE)
  Cells <- unique(GroupInfo)
  total <- length(GroupInfo)

  for(i in 1:length(Cells)){

    thisCell = Cells[i]
    thisCellIdx = which(GroupInfo == thisCell)

    OR <- unname(sapply(1:length(genesets), function(k){
      A <- pres2[names(genesets)[k]] # 전체 Cell에서 유의한 회수
      if(is.na(A)) return(0)
      N <- total # 전체 Cell 수

      B <- table(unlist(pres[thisCellIdx]))[as.character(k)] # 특정 Cell에서 유의한 회수
      if(is.na(B)) return(0)
      K <- length(thisCellIdx)

      return( (B/K) / (A/N) )
    }))
    OR <- round(OR, 4)
    # Cell, Geneset, OR
    res <- rbind(res,
                 data.frame(
                   Cell = as.character(thisCell),
                  Geneset =as.character(names(genesets)),
                  OddRatio = as.numeric(OR), stringsAsFactors = FALSE)
                 )
  }

  colnames(res) <- c('Cell', "Geneset", "OddRatio")

  res <- res %>% filter(OddRatio>0)

  return(res)

}
