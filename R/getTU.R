getTU <- function(CountData, plotOption){
  cat('transpose in getTU started\n')
  CountData <- as.matrix(Matrix::t(CountData))
  cat('t-SNE / U-MAP started\n')
  if (plotOption == "t-SNE") {
    if(nrow(CountData)<500) tsneE <- Rtsne::Rtsne(CountData, check_duplicates = FALSE, perplexity = 15)
    else {
      tsneE <- Rtsne::Rtsne(CountData, check_duplicates = FALSE, perplexity = 15, partial_pca = TRUE) # 15 seconds
    }
    dfobj <- data.frame(tsneE$Y, col = GroupInfo, stringsAsFactors = FALSE)
  }

  if (plotOption == "U-MAP") {
    umapE <- uwot::umap(CountData, fast_sgd = TRUE) # 55 seconds
    dfobj <- data.frame(umapE, col = GroupInfo, stringsAsFactors = FALSE)
  }
  return(dfobj)
}
