#' Context-specific causal effect of single-time point intervention
#'
#' This function estimates the causal effect of a single time-point intervention on the outcome
#' within the same time point. Here, time is defined as observing a single intervention and outcome,
#' and possibly multiple covariates. Intervention is impossed within each unit of time.
#'
#' @param data data.frame object containing the time series with relevant time ordering.
#' @param folds user-specified list of folds- it should correspond to an element of \code{origami}.
#' In case it is not specified, it will defined internally.
#' @param fold_fn cross-validation scheme, as defined by \code{origami}. See \code{origami::fold_funs}
#' for detailed explanations. For time-series, implemented cross-validation schemes are
#' \code{folds_rolling_origin} and \code{folds_rolling_window}.
#' @param window in case \code{fold_fn} was set to \code{folds_rolling_window}, specify the
#' number of observations in each training sample.
#' @param skip in case the time-series considered is very long, it is possible there will be many
#' folds to consider. This parameter allows for few nodes to be skipped. Default is 0, which
#' corresonds to no nodes skipped.
#' @param SL.library list of \code{sl3} algorithms to be used for estimation. For the list of available
#' learners for time-series, see \code{sl3::sl3_list_learners(c("timeseries"))}.
#'
#'
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

tstmleSI <- function(data, Cy=NULL, Ca=NULL, SL.library, folds=NULL, fold_fn="folds_rolling_origin", window=NULL,
                     skip=0){

  ##################################################
  # Estimation with sl3 and specified Cy and Ca.
  ##################################################

  data<-initFrame(data=data, Cy=Cy, Ca=Ca)

  Y<-data.frame(Y=data$Y[,1])
  X<-data$Y[,-1]


  #Fit Q:
  message("Fitting Q")
  estQ<-initEst(data=data,fitQ=TRUE,j=1,folds=folds,fold_fn=fold_fn,window=window,skip=skip,
                SL.library=SL.library)

  #Fit g:
  message("Fitting g")
  estg<-initEst(data=data,fitQ=FALSE,j=1,folds=folds,fold_fn=fold_fn,window=window,skip=skip,
                SL.library=SL.library)

  #We also need Q(C_0,1) and Q(C_0,0):





}

