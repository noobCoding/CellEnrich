#' @export

buildDT = function(pres2){
  datatable(
    data.frame(
      Name = names(pres2),
      Count = as.numeric(pres2)
    )
  )
}
