#' @export

sortItem <- function(label, tableName) {
  options(useFancyQuotes = FALSE)
  paste0(
    "$('#", tableName, "')",
    ".append(", "`<div class=", "'rank-list-item'", " draggable='true'",
    " style = 'transform: translateZ(0px);'>` + ", label, " + `</div>`)"
  )
}
