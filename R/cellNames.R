#' @export

cellNames = function(v){
  v = v@colData
  v = sapply(1:nrow(v), function(i){paste0(v[i,1], '-', v[i,2]) }) # n * 2 ?
  return(v)
}
