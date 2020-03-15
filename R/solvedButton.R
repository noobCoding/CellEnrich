# actionButton with style and onClick Attributes, width = NULL
solvedButton <- function(inputId, label, style = NULL, onClick = NULL, ...) {
  value <- restoreInput(id = inputId, default = NULL)
  tags$button(
    id = inputId, style = style, onClick = onClick,
    type = "button", class = "btn btn-default action-button",
    `data-val` = value, list(label),
    ...
  )
}
