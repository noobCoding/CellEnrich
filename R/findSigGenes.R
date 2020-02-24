#' @title findSigGenes
#'
#' @description define significant genes along cell.
#'
#' @param v count data
#' @param method mean or pos
#'
#' @export
#'
findSigGenes = function(v, method = 'median'){
  if(!method %in% c('median', 'pos')) stop('wrong method')

  Additive = 1

  for(i in 1:nrow(v)){
    v[i,] = log( ( v[i,] + Additive ) / median( v[i,] ) )
  }

  res = list()

  if(method=='median'){
    for(i in 1:nrow(v)){
      idx = which(v[i,]<=median(v[i,]))
      v[i,] = 1
      v[i,idx] = 0
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

  for(i in 1:ncol(v)){
    res[[i]] = sort(names(which(v[,i] > 0)))
  }

  return(res)
}
