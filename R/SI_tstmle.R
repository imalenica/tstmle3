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
#' @param gbounds bounds for the q estimates.
#' @param Qbounds bounds for the Q estimates.
#' @param tol tolerance parameter for the iterative TMLE.
#' @param maxIter maximum number of iterations for the iterative TMLE.
#'
#'
#' @return An object of class \code{tstmle}.
#' \describe{
#' \item{tmlePsi}{Average treatment effect estimated using TMLE.}
#' \item{iptwPsi}{Average treatment effect estimated using IPTW.}
#' \item{tmleSD}{Standard deviation for the TMLE estimated parameter.}
#' \item{tmleCI}{Confidence Interval for the TMLE estimated parameter.}
#' \item{IC}{Influence function.}
#' \item{steps}{Number of steps until convergence of the iterative TMLE.}
#' \item{initialData}{Initial estimates of g, Q.}
#' \item{tmleData}{Final updates estimates of g, Q and clever covariates.}
#' }
#'
#' @importFrom stats qnorm
#'
#' @export
#

tstmleSI <- function(data,Cy=NULL,Ca=NULL,folds=NULL,V=5,stratifyAY = TRUE,
                     Q_library=list("Lrnr_mean", "Lrnr_glm_fast", "Lrnr_glmnet"),
                     g_library=list("Lrnr_mean", "Lrnr_glm_fast", "Lrnr_glmnet"),
                     Q_SLlibrary=list("Lrnr_mean", "Lrnr_arima", "Lrnr_expSmooth"),
                     g_SLlibrary=list("Lrnr_mean", "Lrnr_arima", "Lrnr_expSmooth"),
                     gbounds=c(1e-4,1-1e-4), Qbounds=c(1e-4,1-1e-4),
                     tol=10^-3, maxIter=1000){

  ##################################################
  # Estimation with sl3 and specified Cy and Ca.
  ##################################################

  #Create relevant data:
  data<-initFrame(data=data, Cy=Cy, Ca=Ca)

  #Get full Q and g data-sets:
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
      folds <- make_folds(strata_ids = AYstrata, V = V)
    }else{
      folds <- make_folds(QY, V)
    }
  }

  #Using tmle3 package:
  #tmle_nodes <- list(tmle3::define_node("A", names(g[1]), names(g[,-1])),
  #                   tmle3::define_node("Y", names(Q[1]), names(Q[,-1])))

  #task <- tmle3::tmle_Task$new(Q, tmle_nodes=tmle_nodes, folds=folds)

  #Fit Q:
  message("Fitting Q")
  estQ<-initEst(Y=QY,X=QX,folds=folds,SL.library=Q_library)
  estQ$valY<-bound(estQ$valY, Qbounds)

  #Fit g:
  message("Fitting g")
  estg<-initEst(Y=gY,X=gX,folds=folds,SL.library=g_library)
  estg$valY<-bound(estg$valY, gbounds)

  #Fit Q(C_0,0):
  message("Fitting Q(C_0,0)")
  estQ0<-initEst_setA(data=Q, fit=estQ, setA=0)
  estQ0$valY<-bound(estQ0$valY, Qbounds)

  #Fit Q(C_0,1):
  message("Fitting Q(C_0,1)")
  estQ1<-initEst_setA(data=Q, fit=estQ, setA=1)
  estQ1$valY<-bound(estQ1$valY, Qbounds)

  g_all<-ifelse(g$A==1,estg$valY,1-estg$valY)

  #For now, use gentmle setup
  tmledata <- data.frame(A=g$A, Y=Q$Y, gk=estg$valY, Qk=estQ$valY,
                         Q1k=estQ1$valY, Q0k=estQ0$valY, g=g_all)

  #IPTW:
  iptwEst <- with(tmledata, mean(as.numeric(A==1)*Y/g) - mean(as.numeric(A==0)*Y/g))

  #TMLE:
  tmleEst<-tmleSI(tmledata, tol=tol, maxIter=maxIter,gbounds=gbounds, Qbounds=Qbounds)

  tmlePsi <- tmleEst$psi
  tmleSD <- sqrt(tmleEst$ED2)/sqrt(nrow(tmleEst$tmledata))
  lower <- tmlePsi - stats::qnorm(0.975) * tmleSD
  upper <- tmlePsi + stats::qnorm(0.975) * tmleSD

  return(list(tmlePsi=tmlePsi, iptwPsi=iptwEst, tmleSD=tmleSD, tmleCI=c(lower,upper),
              IC=tmleEst$IC, steps=tmleEst$steps, initialData=tmledata, tmleData=tmleEst$tmledata))

  ##################################################
  # Estimation with just sl3
  ##################################################

  #Additional parameters:

  # @param fold_fn cross-validation scheme, as defined by \code{origami}. See \code{?origami::fold_funs}
  # for detailed explanations. For time-series, implemented cross-validation schemes are
  # \code{folds_rolling_origin} and \code{folds_rolling_window}.
  # @param window in case \code{fold_fn} was set to \code{folds_rolling_window}, specify the
  # number of observations in each training sample.
  # @param skip in case the time-series considered is very long, it is possible there will be many
  # folds to consider. This parameter allows for few nodes to be skipped. Default is 0, which
  # corresonds to no nodes skipped.

  #Fit Q:
  #message("Fitting Q")
  #estQ<-initEst_sl3(data=data,fitQ=TRUE,j=1,folds=folds,fold_fn=fold_fn,window=window,skip=skip,
  #              SL.library=SL.library)

  #Fit g:
  #message("Fitting g")
  #estg<-initEst_sl3(data=data,fitQ=FALSE,j=1,folds=folds,fold_fn=fold_fn,window=window,skip=skip,
  #              SL.library=SL.library)

}

