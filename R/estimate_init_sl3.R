#' Initial Estimation with sl3
#'
#' This function relies on the stacked ensemble learner in order to estimate relevant
#' parts of the likelihood as guided by the efficient influence curve. In particular, it
#' utilizes \code{sl3} package and implemented time-series algorithms in order to
#' estimate Q and g, without having to specify the fixed dimensional summary measure of the past
#' or the past window to condition on.
#'
#' @param data data.frame object containing the time series with relevant time ordering.
#' @param fitQ corresponds to the logical TRUE if we are  estimating Q
#' (conditional density of the outcome). Set to logical FALSE if we are estimating g
#' (conditional density of the exposure) part of the likelihood.
#' @param j size of the artificial batch. This is used in cases where we want to consider
#' multiple time-point interventions, or when we want to define time as multiple occurences
#' of A,Y,W nodes. For single time-point interventions default is 1.
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
#' @return An object of class \code{tstmle}.
#' \describe{
#' \item{fitW}{Fit object for W part of the likelihood.}
#' }
#'
#'
#' @export
#

initEst_sl3 <- function(data, fitQ=TRUE, j=1, folds=NULL, fold_fn="folds_rolling_origin", window=NULL,
                    skip=0, SL.library) {

  # How many in a batch:
  step <- length(grep("_1$", row.names(data), value = TRUE))

  if(j<1){
    stop("j<0 is essentially a non-existent time-series! j should be a positive number.")
  }

  #TO DO: allow for this to happen. Essentially, see where in the time-series you can start.
  #This will cause a shorter time-series, but will allow for a bigger window size.
  if(!is.null(window)){
    if(window>step){
      stop("The window size is larger than the number of available time points.")
    }
  }

  #Need the size of t.
  #If j=1, then t we can get batch size from just counting the different elements with t.
  #Otherwise multiply by j.
  step<-step*j

  #Skip is set to 0 so it is intuitive...
  skip=skip+1

  if (is.null(folds)) {

    #Need another function for blip estimation. Or, make all Ys now the D, for example? Should work.

    #Option to skip some
    batch=step*skip

    if(fitQ==TRUE){
      #Skip the first t: not enough datal start with Y(1).
      first<-step+1
    }else{
      #Skip the first t: not enough data; but start with A(1).
      #Very much assumes we have chronological order A,Y,W
      first<-step
    }

    if(is.null(window)){
      window=first
    }

    if(fold_fn=="folds_rolling_origin"){
      folds <- make_folds(ts(data), fold_fun = origami::folds_rolling_origin, first_window = first,
                                   validation_size = 1, batch=batch)

    }else if(fold_fn=="folds_rolling_window"){
      folds <- make_folds(ts(data), fold_fun = origami::folds_rolling_window, window_size = window,
                                   validation_size = 1, batch=batch)
    }

  }

  covars<-names(data)
  outcome<-names(data)

  #Create sl3 task:
  #TO DO: add option for weights
  task <- sl3::make_sl3_Task(data, covariates = covars, outcome = outcome,
                           outcome_type=NULL, folds=folds)

  #Return sl3 trained object:
  sl3_fit<-sl3.fit(task=task,SL.library=SL.library)

  sl3_fit$sl.fit$predict(task)

  out <- list(cvRisk = sl3_fit$risk, valY = sl3_fit$pred, SL.library = SL.library, folds = folds,
              fullFit = sl3_fit$sl.fit, step = step)
  class(out) <- c("initEst_sl3", "sl3")
  return(out)

}

#' Fit Q with specified exposure
#'
#' Function to make estimation of Q with specified exposure easier.
#'
#' @param data data.frame object containing the time series with relevant time ordering.
#' @param fit fit object from \code{initEst}.
#' @param setA specify the value each A node should be set to.
#'

initEst_sl3_setA<-function(data, fit, setA){

  data_setA<-data
  set<-seq(1, nrow(data), fit$step)

  #Set A node to specified intervention:
  data_setA[set,] = setA

  covars<-names(data_setA)
  outcome<-names(data_setA)

  #Create sl3 task:
  #TO DO: add option for weights
  task_setA <- sl3::make_sl3_Task(data_setA, covariates = covars, outcome = outcome,
                                  outcome_type=NULL, folds=fit$folds)

  #How to predict on edited data in time series?
  sl_preds <- fit$fullFit$predict(task_setA)


}

