#' Nadaraya-Watson Kernel Quanitle Regression
#'
#' Performs  Nadaraya-Watson quantile regression via C code.
#'
#' @param x Vector of predictor values (x-coordinates).
#' @param y Vector of response values (y-coordinates).
#' @param xp Vector of points at which to evaluate the smoother.
#' @param kernel Integer indicating the kernel (1=Box, 2=Gaussian).
#' @param bandwidth The smoothing bandwidth.
#' @return A list with components 'x' (xp) and 'y' (the smoothed values).
#' @export
nw_q <- function(x, y, quantile_level, kernel = c("box", "normal"), bandwidth = 0.5,
                 range.x = range(x), n.points = max(100L, length(x)), x.points)
{
  ## box is [-0.5, 0.5]. normal is sd = 1.4826/4
  if(missing(y) || is.null(y))
    stop("numeric y must be supplied.\nFor density estimation use density()")
  kernel <- match.arg(kernel)
  krn <- switch(kernel, "box" = 1L, "normal" = 2L)
  x.points <-
    if(missing(x.points))
      seq.int(range.x[1L], range.x[2L], length.out = n.points)
  else { n.points <- length(x.points); sort(x.points) }
  ord <- order(x)
  .Call("nw_q", x[ord], y[ord], x.points, krn, bandwidth, quantile_level)
}

