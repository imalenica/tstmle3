#' Context-specific causal effect of multiple time point interventions
#'
#' This function estimates the causal effect of multiple time-point interventions in an ordered
#' longitudinal data structure within time unit t. Here, time is defined as observing multiple ordered
#' logitudinal data structures of user-specified dimension, with possibly multiple intervention and
#' outcome/covariate nodes within it. Time-block defined in this manner could be artificially created
#' for the purpose of estimation of a particular causal effect, of they could have a unique meaning such
#' as representing measurments over a day/week/cycle.
#'
#' @param data data.frame object containing the time series with relevant time ordering.
#' @param Co numeric specifying possible Markov order for the O(t) nodes.
#' @param block size of the new O(t) block: multiple of the original O(t). Default of 1 will give the
#' original time-series data-frame.
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
#' @return An object of class \code{tstmleSI}.
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

tstmleMI <- function(data,Co=NULL,block=1,folds=NULL,V=5,stratifyAY = TRUE,
                     Q_library=list("Lrnr_mean", "Lrnr_glm_fast", "Lrnr_glmnet"),
                     g_library=list("Lrnr_mean", "Lrnr_glm_fast", "Lrnr_glmnet"),
                     gbounds=c(1e-4,1-1e-4), Qbounds=c(1e-4,1-1e-4), maxIter=1000,
                     fold_fn="folds_rolling_origin", window=1, skip=0){

    ###########################################################
    # Estimation with sl3 and specified Co dimension.
    ###########################################################

    #Create relevant data:
    data<-initFrame_mi(data=data, Co=Co, block=block)

    #Number of samples per process
    num<-nrow(data$final[[1]])

    #Get full Q and g data-sets:
    Q<-data$Y
    g<-data$A

    QY<-data.frame(Y=Q$X1)
    QX<-Q[,-1,drop=FALSE]

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
    tmleEst<-tmleSI(tmledata, maxIter=maxIter, Qbounds=Qbounds)

    tmlePsi <- tmleEst$psi
    tmleSD <- sqrt(tmleEst$ED2)/sqrt(nrow(tmleEst$tmledata))
    lower <- tmlePsi - stats::qnorm(0.975) * tmleSD
    upper <- tmlePsi + stats::qnorm(0.975) * tmleSD

  return(list(tmlePsi=tmlePsi, iptwPsi=iptwEst, tmleSD=tmleSD, tmleCI=c(lower,upper),
              IC=tmleEst$IC, steps=tmleEst$steps, initialData=tmledata, tmleData=tmleEst$tmledata))

}
