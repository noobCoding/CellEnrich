getlgs <- function(n){
  sapply(1:length(n), function(i) { length(genesets[[n[i]]]) })
}
