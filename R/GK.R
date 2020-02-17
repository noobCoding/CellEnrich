
#' @title gaussian kernal transform (GSVA)
#'
#' @description transform count data with gaussian kernal using pre calculated cdf
#' this takes ~ 3 min with 20000 * 90 size datas
#'
#' @seealso https://bmcbioinformatics.biomedcentral.com/articles/10.1186/1471-2105-14-7
#' @export
#'
#' @param v count data
#'
#'
#'

GK <- function(v) {
  # pre computed pnorm with 1/10000 resolution
  ppnorms <- sapply(1:100000, function(i) {
    pnorm(i / 10000)
  })

  n <- ncol(v)
  res <- matrix(NA, nrow(v), n)
  for (i in 1:nrow(v)) {
    hi <- sd(v[i, ]) / 4 # gene's standard deviation
    if (hi == 0) hi <- 0.0001
    ei <- unname(v[i, ])
    res[i, ] <- sapply(1:n, function(j) {
      mean(sapply(round((-ei + ei[j]) / hi, 4), ppnorm))
    })
  }
  rownames(res) <- rownames(v)
  colnames(res) <- colnames(v)
  return(res)
}


ppnorm <- function(v) { # precomputed pnorm
  if (v == 0) {
    return(0.5)
  }
  if (v < (-10)) {
    return(0)
  }
  if (v > 10) {
    return(1)
  }
  cdf <- ppnorms[abs(v) * 10000]
  if (v < 0) {
    return(1 - cdf)
  }
  return(cdf)
}


