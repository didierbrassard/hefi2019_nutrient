#' @title Calculate weighted statistics and percentiles
#'
#' @description This function calculates the weighted mean, standard deviation, and percentiles of a numeric vector.
#'
#' @param x A numeric vector.
#' @param w A numeric vector of weights. Default is NULL, which assigns equal weights.
#' @param min_p The minimum percentile to calculate. Default is 1.
#' @param max_p The maximum percentile to calculate. Default is 99.
#' @return A data frame containing the weighted mean, standard deviation, and percentiles.
#'
#' @importFrom stats weighted.mean
#' @importFrom Hmisc wtd.var
#' @importFrom DescTools Quantile
#'
#' @export
#'
#' @examples
#' x <- rnorm(100)
#' w <- rnorm(100)
#' statDistrib(x)
#' statDistrib(x, w)
#' statDistrib(x, w, 5, 95)


statDistrib <- function(x, w=NULL, min_p=1, max_p=99) {

  # Validate input
  if (!is.numeric(x)) {
    stop("statDistrib: Input must be numeric")
  }
  if (!is.null(w) && !is.numeric(w)) {
    stop("statDistrib: Weight vector must be numeric")
  }
  if (!is.null(w) && length(w) != length(x)) {
    stop("statDistrib: Weight vector must have the same length as input")
  }

  # Assign dummy weights if needed
  if (is.null(w)) {
    w <- rep(1, length(x))
    message("statDistrib: Rows were assigned a dummy weight value of 1")
  }

  # Calculate weighted mean, standard deviation, and percentiles
  mean <- stats::weighted.mean(x, w)
  stddev <- sqrt(Hmisc::wtd.var(x, w))
  pctiles <- DescTools::Quantile(x, w, seq(min_p, max_p)/100, type = 7)

  # Note: DescTools::Quantile approx. consistent with SAS estimates

  # Format output
  distrib <- data.frame(
    Mean = mean,
    StdDev = stddev,
    setNames(as.list(pctiles), paste0("Pctile", min_p:max_p))
  )

  return(distrib)
}

