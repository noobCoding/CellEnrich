#' @title findSigGenes
#'
#' @description define significant genes along cell.
#'
#' @param v v
#' @param method mean or median or abs
#'
#' @export
#'
findSigGenes = function(v, method = 'mean'){
  if(!method %in% c('mean','median', 'abs')) stop('wrong method')

  Additive = 1

  for(i in 1:nrow(v)){
    M = mean(v[i,])
    v[i,] = log((v[i,]+Additive)/M)
  }

  res = list()
  for(i in 1:ncol(v)){
    qv90 = quantile(v[,i], probs = 0.95)
    qv10 = quantile(v[,i], probs = 0.05)

    res[[i]] = sort(c(names(which(v[,i] < qv10)) ,names(which(v[,i] > qv90))))
  }
  return(res)
}
