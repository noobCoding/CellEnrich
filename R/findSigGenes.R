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
  if(!method %in% c('median', 'mean', 'zero')) stop('wrong method')
  # it's already matrix
  # v <- as.matrix(v)
  cat('findSigGenes started\n')

  rownames(v) <- colnames(v) <- NULL

  res = list()

  # 100 : 0.8 seconds
  # 200 : 1.12
  # 250 : 1.36. find 20000 gene will take less than minute

  cat('scaling\n')

  idx <- floor(nrow(v) / 250)
  v2 = c()
  for(i in 1:idx){
    thisIdx <- 1:250 + 250*(i-1)
    vv <- as.matrix(v[thisIdx,])+1
    meds <- apply(vv, 1, median)
    vv <- as(log(sweep(vv, 1, meds, '/')), 'dgCMatrix')
    v2 <- rbind(v2, vv) # use rbind, not assign
  }

  if(nrow(v)%%250!=0){
    thisIdx <- (idx*250+1):nrow(v)
    vv <- as.matrix(v[thisIdx,])+1
    meds <- apply(vv, 1, median)
    vv <- log(sweep(vv, 1, meds, '/'))
    v2 <- rbind(v2, vv)
  }

  v <- v2
  rm(v2)

  cat('define Lists\n')

  if(method=='median'){
    for(i in 1:ncol(v)){
      res[[i]] = which(v[,i] > median(v[,i]))
    }
  }
  if(method == 'zero'){
    for(i in 1:ncol(v)){
      res[[i]] = which(v[,i] > 0)
    }
  }
  if(method == 'mean'){
    for(i in 1:ncol(v)){
      res[[i]] = which(v[,i]> mean(v[,i]))
    }
  }

  names(res) <- Name
  return(res)

}
