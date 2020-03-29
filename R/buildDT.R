#' @export

buildDT <- function(pres2) {
  DT::datatable(
    data.frame(
      Geneset = names(pres2),
      Count = as.numeric(pres2)
    ),
    options = list(
      dom = "ltp",
      lengthChange = FALSE
    ),
    rownames = FALSE,
    selection = "single"
  )
}
