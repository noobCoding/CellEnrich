briterhex <- function(colors) {
  res <- c()
  for (i in 1:length(colors)) {
    v <- as.vector(col2rgb(colors[i])) * 1.3
    v <- sapply(v, function(i) {
      min(i, 255)
    })
    res[i] <- rgb(v[1], v[2], v[3], max = 255)
  }
  return(res)
}
