#' @title findSigGenesGroup
#'
#' @description define significant genes along Group.
#'
#' @param Count Count data
#' #' @param ClustInfo Group information
#' @param q0 qvalue cutoff, default is 0.1
#' @param TopCutoff cutoff for Top in findMarkers, default is 5.
#'
#' @export
#' @import scran
#'
#'
#'
findSigGenesGroup = function(Count = NULL, ClustInfo = NULL, q0 = 0.1, TopCutoff = 5){
  require(scran)

  if(is.null(Count)) stop('Count must given')
  if(is.null(ClustInfo)) stop('ClustInfo must given')

  GrpRes = scran::findMarkers(x = Count, ClustInfo, test = 'wilcox', direction = 'up')
  Grp = unique(ClustInfo)

  res = data.frame(stringsAsFactors = FALSE)

  for(i in 1:length(Grp)){
    G = data.frame(
          genes = rownames(GrpRes[[1]]),
          Group = Grp[i],
          GrpRes[[1]],
          row.names = NULL,
          stringsAsFactors = FALSE
        ) %>%
        select(Group, Top, genes, FDR) %>%
        filter(FDR <= q0) %>%
        # filter(Top <= TopCutoff) %>%
        arrange(FDR)

    res = rbind(res, G)
  }
  res$genes = as.character(res$genes)
  res$Group = as.character(res$Group)
  res$FDR <- round(as.numeric(res$FDR), 6)
  return(res)
}
