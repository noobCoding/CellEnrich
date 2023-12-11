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
      "Human-Reactome", # 2022
      "Human-WikiPathway", # 2021
      "Human-KEGG", # KEGG 2021
      "Human-GOBP",
      "Human-GOCC",
      "Human-GOMF",
      "Mouse-Reactome",
      "Mouse-WikiPathway", # 2019
      "Mouse-KEGG", # 2019
      "Mouse-GOBP",
      "Mouse-GOCC",
      "Mouse-GOMF",
      # Datasets
      'Alzheimer_Counts_sampled.RDS',
      'Alzheimer_CellType_sampled.RDS',
      "LPS_exp.rds",
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
      "Human-Reactome",
      "Mouse-WikiPathway",
      
      # Datasets
      'Alzheimer_Counts_sampled.RDS',
      'Alzheimer_CellType_sampled.RDS',
      "LPS_exp.rds",
      "pbmcClustInfo.RData",
      "pbmcData.RData"
    )
  }

  filelist <- setdiff(filelist, dir())
  urls <- paste0("https://github.com/noobcoding/CellEnrich/raw/master/data/",filelist)
  for (i in 1:length(urls)) {
    download.file(urls[i], filelist[i])
  }

}

