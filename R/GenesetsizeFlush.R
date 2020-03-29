GenesetsizeFlush = function(genesets, lgs, minsize, maxsize){
  genesets[intersect(which(lgs >= minsize), which(lgs <= maxsize))]
}
