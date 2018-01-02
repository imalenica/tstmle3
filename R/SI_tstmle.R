#' Context-specific causal effect of single-time point intervention
#'
#' This function estimates the causal effect of a single time-point intervention on the outcome
#' within the same time point. Here, time is defined as observing a single intervention and outcome,
#' and possibly multiple covariates. Intervention is impossed within each unit of time.
#'
#' @param data data.frame object containing the time series with relevant time ordering.
#' @param Cy numeric specifying possible Markov order for Y nodes.
#' @param Ca numeric specifying possible Markov order for A nodes.
#' @param folds user-specified list of folds- it should correspond to an element of \code{origami}.
#' In case it is not specified, it will defined internally.
#' @param V number of cross-validation folds used.
#' @param stratifyAY logical: should we stratify the cross-validation based on (A,Y) pairs
#' @param Q_library list of \code{sl3} algorithms for the fit of E(Y|A,Cy) (Q)
#' @param g_library list of \code{sl3} algorithms for the fit of P(A|Ca) (g)
#' @param Q_SLlibrary list of \code{sl3} algorithms to be used for estimation of E(Y|A,Cy) (Q).
#' For the list of available learners for time-series, see \code{sl3::sl3_list_learners()}.
#' @param q_SLlibrary list of \code{sl3} algorithms to be used for estimation of P(A|Ca) (g).
#' For the list of available learners for time-series, see \code{sl3::sl3_list_learners()}.
#'
#'
#'
#' @param fold_fn cross-validation scheme, as defined by \code{origami}. See \code{?origami::fold_funs}
#' for detailed explanations. For time-series, implemented cross-validation schemes are
#' \code{folds_rolling_origin} and \code{folds_rolling_window}.
#' @param window in case \code{fold_fn} was set to \code{folds_rolling_window}, specify the
#' number of observations in each training sample.
#' @param skip in case the time-series considered is very long, it is possible there will be many
#' folds to consider. This parameter allows for few nodes to be skipped. Default is 0, which
#' corresonds to no nodes skipped.
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

tstmleSI <- function(data,Cy=NULL,Ca=NULL,folds=NULL,V=5,stratifyAY = TRUE,
                     Q_library=list("Lrnr_mean", "Lrnr_glm_fast", "Lrnr_glmnet"),
                     g_library=list("Lrnr_mean", "Lrnr_glm_fast", "Lrnr_glmnet"),
                     Q_SLlibrary=list("Lrnr_mean", "Lrnr_arima", "Lrnr_expSmooth"),
                     g_SLlibrary=list("Lrnr_mean", "Lrnr_arima", "Lrnr_expSmooth")
                     ){

  ##################################################
  # Estimation with sl3 and specified Cy and Ca.
  ##################################################

  #Create relevant data:
  data<-initFrame(data=data, Cy=Cy, Ca=Ca)

  Q<-data$Y
  g<-data$A

  QY<-data.frame(Y=data$Y[,1])
  QX<-data$Y[,-1]

  gY<-data.frame(A=data$A[,1])
  gX<-data$A[,-1]

  #Make folds
  if (is.null(folds)) {
    if(stratifyAY){
      AYstrata <- sprintf("%s %s", QX[, 1], QY[, 1])
      #Stratified folds:
      folds <- origami::make_folds(strata_ids = AYstrata, V = V)
    }else{
      folds <- origami::make_folds(QY, V)
    }
  }

  #Fit Q:
  message("Fitting Q")
  estQ<-initEst(Y=QY,X=QX,folds=folds,SL.library=Q_library)

  #Fit g:
  message("Fitting g")
  estg<-initEst(Y=gY,X=gX,folds=folds,SL.library=g_library)

  #Fit Q(C_0,0):
  message("Fitting Q(C_0,0)")
  estQ0<-initEst_setA(data=Q, fit=estQ, setA=0)

  #Fit Q(C_0,1):
  message("Fitting Q(C_0,1)")
  estQ1<-initEst_setA(data=Q, fit=estQ, setA=1)

  #For now, use gentmle
  tmledata <- data.frame(A=g$A, Y=Q$Y, gk=estg$valY, Qk=estQ$valY, Q1k=estQ1$valY, Q0k=estQ0$valY)























  ##################################################
  # Estimation with just sl3
  ##################################################

  #Fit Q:
  message("Fitting Q")
  estQ<-initEst_sl3(data=data,fitQ=TRUE,j=1,folds=folds,fold_fn=fold_fn,window=window,skip=skip,
                SL.library=SL.library)

  #Fit g:
  message("Fitting g")
  estg<-initEst_sl3(data=data,fitQ=FALSE,j=1,folds=folds,fold_fn=fold_fn,window=window,skip=skip,
                SL.library=SL.library)

  #We also need Q(C_0,1) and Q(C_0,0):





}

