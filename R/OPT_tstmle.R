#' Context-specific optimal treatment effect for a single time series
#'
#' This function estimates the optimal individualized treatment effect for a single time series.
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
#' @param gbounds bounds for the q estimates.
#' @param Qbounds bounds for the Q estimates.
#' @param maxIter maximum number of iterations for the iterative TMLE.
#'
#' @return An object of class \code{tstmleOPT}.
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

tstmleOPT <- function(data,Co=TRUE,Cy=NULL,Ca=NULL,folds=NULL,V=5,stratifyAY = TRUE,
                     Q_library=list("Lrnr_mean", "Lrnr_glm_fast", "Lrnr_glmnet"),
                     g_library=list("Lrnr_mean", "Lrnr_glm_fast", "Lrnr_glmnet"),
                     gbounds=c(1e-4,1-1e-4), Qbounds=c(1e-4,1-1e-4), maxIter=1000){

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

  #Fit Q:
  message("Fitting Q")
  estQ<-initEst(Y=QY,X=QX,folds=folds,SL.library=Q_library)
  estQ$valY<-bound(estQ$valY, Qbounds)

  #Fit g:
  message("Fitting g")
  estg<-initEst(Y=gY,X=gX,folds=folds,SL.library=g_library)
  estg$valY<-bound(estg$valY, gbounds)

  #Split-specific predictions:
  message("Generating split-specific predictions")
  estSplit <- cross_validate(cv_split, folds, Q, estQ, estg)
















}
