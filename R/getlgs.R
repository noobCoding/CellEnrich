getlgs <- function(genesets){
  sapply(1:length(genesets), function(i) { length(genesets[[i]]) })
}
