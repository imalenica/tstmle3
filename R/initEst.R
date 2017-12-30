#' Initial Estimation with sl3 and specified C_{o(t)} dimension
#'
#' This function relies on the stacked ensemble learner in order to estimate relevant
#' parts of the likelihood as guided by the efficient influence curve.
#' (best used for short-term dependence).
#'
#' @param Y data.frame object containing the outcome variable.
#' @param X data.frame object containing the relevant covariates.
#' @param folds user-specified list of folds- it should correspond to an element of \code{origami}.
#' In case it is not specified, it will defined internally.
#' @param fold_fn cross-validation scheme, as defined by \code{origami}. See \code{origami::fold_funs}
#' for detailed explanations. For time-series, implemented cross-validation schemes are
#' \code{folds_rolling_origin} and \code{folds_rolling_window}.
#' @param SL.library list of \code{sl3} algorithms to be used for estimation. For the list of available
#' learners for time-series, see \code{sl3::sl3_list_learners(c("timeseries"))}.
#'
#'
#' @return An object of class \code{tstmle}.
#' \describe{
#' \item{fitW}{Fit object for W part of the likelihood.}
#' }
#'
#' @importFrom origami make.folds
#'
#' @export
#

initEst <- function(Y, X ) {










}

#' Make initial dataframe
#'
#' Function to create the initial data.frame \code{sl3} would expect with i.i.d data.
#'
#' @param data data.frame object containing the time series with relevant time ordering.
#' @param Cy numeric specifying possible Markov order for Y nodes.
#' @param Ca numeric specifying possible Markov order for A nodes.
#'
#' @return An object of class \code{tstmle}.
#' \describe{
#' \item{Y}{Observed outcome with relevant lagged time-points as specified by Cy.}
#' \item{A}{FObserved exposure with relevant lagged time-points as specified by Ca.}
#' \item{Ca}{Input Markov order for A.}
#' \item{Cy}{Input Markov order for Y.}
#' \item{data}{Input data.}
#' }
#'
#' @importFrom Hmisc Lag
#' @importFrom stats complete.cases
#'
#' @export
#

initFrame <- function(data, Cy, Ca){

  # How many in a batch:
  step <- length(grep("_1$", row.names(data), value = TRUE))

  names<-row.names(data)

  # Lag past
  res_y <- lapply(seq_len(Cy), function(x) {
    Hmisc::Lag(data[, 1], x)
  })
  res_a <- lapply(seq_len(Ca), function(x) {
    Hmisc::Lag(data[, 1], x)
  })

  res_y <- data.frame(matrix(unlist(res_y), nrow = length(res_y[[1]])),
                    stringsAsFactors = FALSE)
  res_a <- data.frame(matrix(unlist(res_a), nrow = length(res_a[[1]])),
                      stringsAsFactors = FALSE)

  data_full_y <- cbind.data.frame(data = data, res_y)
  data_full_a <- cbind.data.frame(data = data, res_a)

  #Start will depend on values of Cy and Ca.
  #Keep only samples that have full history
  cc <- stats::complete.cases(data_full_y)
  data_full_y <- data_full_y[cc, ]

  cc <- stats::complete.cases(data_full_a)
  data_full_a <- data_full_a[cc, ]

  #Separate by process:
  Y <- data_full_y[grep("Y", row.names(data_full_y), value = TRUE), ]
  A <- data_full_a[grep("A", row.names(data_full_a), value = TRUE), ]

  #Shorten the longer one if not equal:
  if(nrow(Y) != nrow(A)){

    warning("Note that having different dimensions of Cy and Ca will result in a smaller n here.")

    if(nrow(Y) < nrow(A)){
      A<-A[(nrow(A)-nrow(Y)+1):nrow(A),]
    }else if(nrow(Y) > nrow(A)){
      Y<-Y[(nrow(Y)-nrow(A)+1):nrow(Y),]
    }

  }

  out <- list(Y = Y, A = A, Cy=Cy, Ca=Ca, data=data)
  return(out)

}

