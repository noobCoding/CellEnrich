#'
#' @name GK
#'
#' @title Gaussian Kernal transform (GSVA)
#'
#' @param v CountData [dgCMatrix]
#'
#' @description Transform count data with gaussian kernal using pre-calculated cdf
#' This takes ~ 3 min with 20000 * 90 size datas, (Really Really Slow)
#'
#' @return transformed CountData
#'
#' @seealso https://bmcbioinformatics.biomedcentral.com/articles/10.1186/1471-2105-14-7
#'
#' @export
#'

GK <- function(v) {
  # pre computed pnorm with 1/10000 resolution

  ppnorms <- sapply(1:100000, function(i) {
    pnorm(i / 10000)
  })

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

  n <- ncol(v)

  v2 = c()

  for(i in 1:nrow(v)){
    vi <- unname(v[i,])
    hi <- sd(vi)
    if(hi == 0){ # all values are same,
      v2 <- rbind(v2, rep(0.5, n))
    }
    else{

      zeroIdx <- unname(which(vi==0))
      Idx <- unname(which(vi!=0))

      vs <- vi[Idx]
      L <- n - length(Idx)

      # zero values are same
      vi[zeroIdx] <- ( sum(sapply(-vs/hi, ppnorm)) + 0.5 * L) / n

      # non - zero values
      vi[Idx] <- unname(sapply(vs, function(j){
        ( sum(sapply((j-vs)/hi, ppnorm)) + ppnorm(j/hi) * L ) / n
      }))

      v2 <- rbind(v2, vi)
    }

  }

  return(v2)
}


