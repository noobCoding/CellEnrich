#' @title findSigGenes
#'
#' @description define significant genes along cell.
#'
#' @param v count data
#' @param method mean or pos
#'
#' @export
#'
findSigGenes = function(v, method = 'median', Name){
  if(!method %in% c('median', 'pos')) stop('wrong method')

  v = as.matrix(v)
  rownames(v) = colnames(v) = NULL
  Additive = 1
  v = v + Additive
  meds = apply(v, 1, median)

  v = log(sweep(v, 1, meds, '/'))

  res = list()

  meds = apply(v, 1, median)

  if(method=='median'){
    for(i in 1:nrow(v)){
      res[[i]] = which(v[i,] > median(v[i,]))
    }
  }

  if(method == 'zero'){
    for(i in 1:nrow(v)){
      idx = which(v[i,]<=0)
      v[i,] = 1
      v[i,idx] = 0
    }
  }
  if(method == 'mean'){
    for(i in 1:nrow(v)){
      idx = which(v[i,]<= mean(v[i,]))
      v[i,] = 1
      v[i,idx] = 0
    }
  }
  names(res) <- Name
  return(res)
}
