#' @export
#'
#'
findPathway = function(s){
  res = list()
  for(i in 1:length(s)){
    pvh = getHyperPvalue(s[[i]] , genesets)
    qvh = p.adjust(pvh, 'fdr')
    res[[i]] = unname(which(qvh<0.1))
  }

  return(res)
}
