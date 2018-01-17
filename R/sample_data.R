#' Mock Time Series Data Set
#'
#' A dataset with a simple data structure O = (A, Y, W), where exposure (A) and outcome (Y) are
#' binary. This is a simple dataset designed specifically to illustrate the TMLE estimation procedure.
#'
#' @format A \code{data.frame} with 1 column.
#' \describe{
#'   \item{Y}{A binary variable representing an outcome of interest.}
#'   \item{A}{A binary variable representing an intervention of interest.}
#'   \item{W1}{A binary variable representing a covariate of interest.}
#'   \item{W2}{A continuous variable representing a covariate of interest.}
#'   \item{W3}{A continuous variable representing a covariate of interest.}
#' }

"sim_ts_s1"

#' Smaller mock Time Series Data Set
#'
#' A dataset with a simple data structure O = (A, Y, W), where exposure (A) and outcome (Y) are
#' binary. This is a simple, smaller, dataset designed specifically to illustrate the
#' TMLE estimation procedure.
#'
#' @format A \code{data.frame} with 1 column.
#' \describe{
#'   \item{Y}{A binary variable representing an outcome of interest.}
#'   \item{A}{A binary variable representing an intervention of interest.}
#'   \item{W1}{A binary variable representing a covariate of interest.}
#'   \item{W2}{A continuous variable representing a covariate of interest.}
#'   \item{W3}{A continuous variable representing a covariate of interest.}
#' }

"sim_ts_s1_n50"
