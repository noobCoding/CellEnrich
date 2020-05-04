#'
#' @title Download additional CellEnrich Data
#'
#' @usage DownloadData()
#'
#'
#' @export

#' @export

DownloadData <- function() {

  filelist <- c(
    "c2v7.RData",
    "c5v7.RData",
    "humanGO.RData",
    "humanGOBP.RData",
    "humanGOCC.RData",
    "humanGOMF.RData",
    "keggv7.RData",
    "mouseGO.RData",
    "mouseGOBP.RData",
    "mouseGOCC.RData",
    "mouseGOMF.RData",
    "mouseKEGG.RData"
  )

  filelist <- setdiff(filelist, dir())
  urls <- paste0("https://github.com/unistbig/CellEnrich/raw/master/",filelist)
  for (i in 1:length(urls)) {
    download.file(urls[i], filelist[i])
  }

}
