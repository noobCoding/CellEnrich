#' @export
extractArray = function(v){
  return(v@assays$data$normcounts)
}
