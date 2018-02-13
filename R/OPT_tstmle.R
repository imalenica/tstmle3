#' Context-specific optimal treatment effect for a single time series
#'
#' This function estimates the optimal individualized treatment effect for a single time series,
#' under a single intervention on the exposure A. In addition, it estimates the
#' context-specific mean of the outcome under other user-specified rules,
#' not necessarily optimal regimes.
#'
#' @param data data.frame object containing the time series with relevant time ordering.
#' @param Cy numeric specifying possible Markov order for Y nodes.
#' @param Ca numeric specifying possible Markov order for A nodes.
#' @param V number of cross-validation folds used.
#' @param stratifyAY logical: should we stratify the cross-validation based on (A,Y) pairs
#' @param Q_library list of \code{sl3} algorithms for the fit of E(Y|A,Cy) (Q)
#' @param g_library list of \code{sl3} algorithms for the fit of P(A|Ca) (g)
#' @param blip_library list of \code{sl3} algorithms for the fit of the blip function.
#' @param gbounds bounds for the q estimates.
#' @param Qbounds bounds for the Q estimates.
#' @param maxIter maximum number of iterations for the iterative TMLE.
#'
#' @return An object of class \code{tstmle}.
#' \describe{
#' \item{tmlePsi}{Context-specific mean of the outcome under user-specified \code{ruleA} estimated
#' using TMLE methodology. In particular, it returns psi under the optimal regime, observed exposure,
#' A=1 and A=0.}
#' \item{tmleSD}{Standard deviation of the context-specific mean of the outcome under
#' user-specified \code{ruleA} estimated using TMLE methodology. In particular, it returns sd under the
#' optimal regime, observed exposure, A=1 and A=0.}
#' \item{tmleCI}{Confidence interval of the context-specific mean of the outcome under
#' user-specified \code{ruleA} estimated using TMLE methodology. It returns CI under the
#' optimal regime, observed exposure, A=1 and A=0.}
#' \item{IC}{Influence curve for the context-specific parameter under user-specified \code{ruleA}.
#' It returns IC under the optimal regime, observed exposure, A=1 and A=0.}
#' \item{rule}{Used rule for the exposure.}
#' \item{steps}{Number of steps until convergence of the iterative TMLE for each rule.}
#' \item{initialData}{Initial estimates of g and Q, and observed A and Y.}
#' \item{tmleData}{Final updates estimates of g, Q and clever covariates.}
#' \item{all}{Full results of \code{tmleOPT}.}
#' }
#'
#' @export
#

tstmleOPT <- function(data,Cy=1,Ca=1,V=5,stratifyAY = TRUE,
                     Q_library=list("Lrnr_mean", "Lrnr_glm_fast", "Lrnr_glmnet"),
                     g_library=list("Lrnr_mean", "Lrnr_glm_fast", "Lrnr_glmnet"),
                     blip_library=list("Lrnr_mean", "Lrnr_glm_fast", "Lrnr_glmnet","Lrnr_randomForest","Lrnr_xgboost"),
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
  if(stratifyAY){
    AYstrata <- sprintf("%s %s", QX[, 1], QY[, 1])
    #Stratified folds:
    folds <- origami::make_folds(strata_ids = AYstrata, V = V)
    }else{
      folds <- origami::make_folds(QY, V)
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
  estSplt<-estSplit(folds, Q, g, estQ, estg)

  #Fit blip:
  message("Fitting the blip function")
  estBlp<-estBlip(folds=folds, estSplt=estSplt, Q=Q, blip_library=blip_library)

  #Run TMLE for all rules:
  rule0 <- ruletmle(estSplt$valSplit$A, estSplt$valSplit$Y, estSplt$valSplit$pA1,
                      estSplt$valSplit$Q0W, estSplt$valSplit$Q1W, ruleA=0)
  rule1 <- ruletmle(estSplt$valSplit$A, estSplt$valSplit$Y, estSplt$valSplit$pA1,
                      estSplt$valSplit$Q0W, estSplt$valSplit$Q1W, ruleA=1)
  ruleA <- ruletmle(estSplt$valSplit$A, estSplt$valSplit$Y, estSplt$valSplit$pA1,
                      estSplt$valSplit$Q0W, estSplt$valSplit$Q1W, ruleA=unlist(estSplt$valSplit$A))
  ruleOpt <- ruletmle(estSplt$valSplit$A, estSplt$valSplit$Y, estSplt$valSplit$pA1,
                        estSplt$valSplit$Q0W, estSplt$valSplit$Q1W, ruleA=estBlp$optA)
  res <- list(rule0=rule0,rule1=rule1,ruleA=ruleA,ruleOpt=ruleOpt)

  #Extract psi for each rule:
  res_fin<-extract_res(res)

  return(list(tmlePsi=res_fin$tmlePsi,tmleSD=res_fin$tmleSD,tmleCI=res_fin$tmleCI,
              IC=res_fin$IC,rule=res_fin$rule,steps=res_fin$steps,initialData=res_fin$initialData,
              tmleData=res_fin$tmleData,all=res_fin$all))

}
