getlgs <- function(genesets){
  n <- names(genesets)
  cat('getlgs\n')
  sapply(1:length(n), function(i) { length(genesets[[n[i]]]) })
}
