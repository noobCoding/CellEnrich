#' @title findSigGenesGroup
#'
#' @description define significant genes along Group.
#'
#' @param Count Count data
#' @param q0 qvalue cutoff, default is 0.1
#' @param ClustInfo
#'
#' @export
#' @import scran
#'
#'
#'
findSigGenesGroup = function(Count = NULL, ClustInfo = NULL, q0 = 0.1){
  require(scran)

  if(is.null(Count)) stop('Count must given')
  if(is.null(ClustInfo)) stop('ClustInfo must given')

  GrpRes = scran::findMarkers(x = Count, ClustInfo, test = 'wilcox', direction = 'up')
  Grp = unique(ClustInfo)
  res = list()

  for(i in 1:length(Grp)){
    Genes = rownames(GrpRes[[i]])
    res[[i]] = Genes[which(p.adjust(GrpRes[[i]]$p.value,'fdr')<q0)]
  }
  return(res)
}
