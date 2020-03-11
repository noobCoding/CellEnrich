#' @export
#'
transpose = function(CountData, genesets){
  cat('CountData :', nrow(CountData), '\n')
  cat('Genesets :', length(genesets), '\n')

  genes = rownames(CountData)
  AllGenes = intersect(unique(unlist(genesets)), genes)

  CountData = CountData[which(rownames(CountData) %in% AllGenes),]
  for(i in 1:length(genesets)){
    genesets[[i]] = intersect(genesets[[i]], AllGenes)
  }

  cat("intersect :\n")
  cat('CountData :', nrow(CountData), '\n')
  cat('Genesets :', length(genesets), '\n')

  return( list(CountData = CountData, genesets = genesets) )
}
