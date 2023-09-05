#'
#' @title Download additional CellEnrich Data
#'
#' @usage DownloadData()
#'
#'
#' @export

DownloadData <- function(type='test') {

  if (type=='all') {

    filelist <- c(
      # Gene-sets
      "Reactome_2022.RData",
      "WikiPathways_2021_Human.RData",
      "WikiPathways_2019_Mouse.RData",
      "KEGG_2021_Human.RData",
      "KEGG_2019_Mouse.RData",
      "humanGO.RData",
      "humanGOBP.RData",
      "humanGOCC.RData",
      "humanGOMF.RData",
      "humanKEGG.RData",
      "mouseGO.RData",
      "mouseGOBP.RData",
      "mouseGOCC.RData",
      "mouseGOMF.RData",
      "mouseKEGG.RData",
      # Datasets
      'Alzheimer_Counts_sampled.RDS',
      'Alzheimer_CellType_sampled.RDS',
      "koh.RData",
      "kohInfo.RData",
      "mouseClustInfo.RData",
      "mouseData.RData",
      "pbmcClustInfo.RData",
      "pbmcData.RData"
    )
  } else if (type=='test') {
    filelist <- c(
      # Gene-sets
      "Reactome_2022.RData",
      "WikiPathways_2021_Human.RData",
      "KEGG_2021_Human.RData",
      "humanGO.RData",
      "humanGOBP.RData",
      "humanGOCC.RData",
      "humanGOMF.RData",
      "humanKEGG.RData",

      # Datasets
      'Alzheimer_Counts_sampled.RDS',
      'Alzheimer_CellType_sampled.RDS',
      "pbmcClustInfo.RData",
      "pbmcData.RData"
    )
  }

  filelist <- setdiff(filelist, dir())
  urls <- paste0("https://github.com/noobcoding/CellEnrich/blob/master/data/",filelist)
  for (i in 1:length(urls)) {
    download.file(urls[i], filelist[i])
  }

}

