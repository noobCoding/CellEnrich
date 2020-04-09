GenesetsizeFlush = function(genesets, lgs, minsize, maxsize){
  cat('GenesetSizeFlush\n')
  genesets[intersect(which(lgs >= minsize), which(lgs <= maxsize))]
}
