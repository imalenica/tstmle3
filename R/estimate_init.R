#' Initial Estimation with sl3 and specified C_{o(t)} dimension
#'
#' This function relies on the stacked ensemble learner in order to estimate relevant
#' parts of the likelihood as guided by the efficient influence curve.
#' (best used for short-term dependence).
#'
#' @param Y data.frame object containing the outcome variable.
#' @param X data.frame object containing the relevant covariates.
#' @param folds user-specified list of folds- it should correspond to an element of \code{origami}.
#' @param SL.library list of \code{sl3} algorithms to be used for estimation.
#'
#'
#' @return An object of class \code{tstmle}.
#' \describe{
#' \item{cvRisk}{CV Risk for each algorithm and Super Learner.}
#' \item{valY}{Prediction based on the estimated SL fit.}
#' \item{SL.library}{list of \code{sl3} algorithms to be used for estimation.}
#' \item{folds}{user-specified list of folds- it should correspond to an element of \code{origami}.}
#' \item{fullFit}{Fit on all samples.}
#' \item{cvFit}{CV fit.}
#' }
#'
#' @importFrom sl3 make_sl3_Task
#'
#' @export
#

initEst <- function(Y, X, folds=NULL,SL.library, outcome_type="binary") {

  covars<-names(X)
  outcome<-names(Y)
  data<-cbind.data.frame(Y,X)

  #Create sl3 task:
  #TO DO: add option for weights
  task <- sl3::make_sl3_Task(data, covariates = covars, outcome = outcome,
                             outcome_type=outcome_type, folds=folds)

  #Return sl3 trained object:
  sl3_fit<-sl3.fit(task=task,SL.library=SL.library)

  out <- list(cvRisk = sl3_fit$risk, valY = sl3_fit$pred, SL.library = SL.library, folds = folds,
              fullFit = sl3_fit$sl.fit, cvFit=sl3_fit$cv.fit)

  class(out) <- c("initEst", "sl3")

  return(out)

}

#' Make initial dataframe
#'
#' Function to create the initial data.frame \code{sl3} would expect with i.i.d data.
#' The newly created data will rely on the specified Cy and Ca.
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
  name<-row.names(data)

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

  mY<-match(row.names(Y)[1],name)
  mA<-match(row.names(A)[1],name)

  names(Y)[2:ncol(Y)]<-name[(mY-1):(mY-Cy)]
  names(Y)[1:2]<-c("Y","A")

  names(A)[2:ncol(A)]<-name[(mA-1):(mA-Ca)]
  names(A)[1]<-c("A")

  row.names(Y)<-NULL
  row.names(A)<-NULL

  out <- list(Y = Y, A = A, Cy=Cy, Ca=Ca, data=data)
  return(out)

}

#' Make initial dataframe for multiple intervention setting
#'
#' Function to create the initial data.frame \code{sl3} would expect with i.i.d data.
#' The newly created data will rely on the specified Cl and Ca, as well as on the size of the O(t)
#' block.
#'
#' @param data data.frame object containing the time series with relevant time ordering.
#' @param Co numeric specifying possible Markov order for the O(t) nodes.
#' @param block size of the new O(t) block: multiple of the original O(t). Default of 1 will give the
#' original time-series data-frame.
#'
#' @return An object of class \code{tstmle}.
#' \describe{
#' \item{Y}{data set containing the \code{sl3} ready input for the output and relevant covariates.}
#' \item{A}{data set containing the \code{sl3} ready input for all the exposures and relevant covariates.}
#' \item{final}{Final dataset with each list component representing a data set with targeted node
#' list of covariates as specified by Co and batch t.}
#' \item{original}{Final dataset with each list component representing a data set with targeted node
#' list of covariates as specified by Co.}
#' \item{step}{Number of nodes in an initial batch.}
#' \item{card}{Cardinality of the newly created block.}
#' }
#'
#' @importFrom Hmisc Lag
#' @importFrom stats complete.cases
#' @importFrom plyr ldply
#'
#' @export
#

initFrame_mi <- function(data, Co, block=1){

  # How many in a single, initial batch:
  step <- length(grep("_1$", row.names(data), value = TRUE))
  name<-row.names(data)

  #Get initial batch:
  name_ib<-row.names(data)[1:step]
  nY<-step-which(name_ib=="Y_0")
  nA<-which(name_ib=="A_0")

  #Cardinality of the newly created block:
  card<-block*step

  #Final Y in a block:
  fY<-card-nY

  #Indices for A processes:
  iA<-seq(from=nA,to=card,by = step)

  #Determine block-specific Co:
  res_Co <- lapply(seq_len(Co)+card, function(x) {
    Hmisc::Lag(data[, 1], x)
  })

  res_Co <- data.frame(matrix(unlist(res_Co), nrow = length(res_Co[[1]])),
                      stringsAsFactors = FALSE)
  data_full_Co <- cbind.data.frame(data = data, res_Co)

  #Start will depend on the dimension of Co
  #Keep only samples that have full history
  cc <- stats::complete.cases(data_full_Co)
  data_full_Co <- data_full_Co[cc, ]

  #Define Co:
  #it will be defined by the lags of A, since they are the first ones in the new batch.
  #it also represent the closest history.
  Co_data<-data_full_Co[grep("A", row.names(data_full_Co), value = TRUE), ]
  first_batch<-row.names(Co_data)[1]

  #Time-series data used for estimation will start at first_batch:
  data_new<-data.frame(data=data[which(row.names(data) == first_batch):nrow(data),1,drop=FALSE])

  r<-NULL
  for(i in 1:(nrow(data_new)/card)){
    rn<-data.frame(rep(i,card))
    r<-rbind.data.frame(r,rn)
  }
  names(r)<-"batch"

  data_new<-data.frame(data=data_new[1:nrow(r),,drop=FALSE])
  data_new<-cbind.data.frame(data_new,batch=r)

  #Separate by the size of a batch
  batches<-split(data_new, data_new$batch, drop = TRUE)

  #Lag prior values within the batch
  for(i in 1:length(batches)){
    res_batches <- lapply(seq_len(card)-1, function(x) {
      Hmisc::Lag(batches[[i]]$data, x)
    })

    res_batches <- data.frame(matrix(unlist(res_batches), nrow = length(res_batches)),
                         stringsAsFactors = FALSE)
    batches[[i]] <- cbind.data.frame(batches[[i]], res_batches)
  }

  data_comb<-plyr::ldply(batches, data.frame)
  row.names(data_comb)<-row.names(data_new)

  #Separate by process and batch:
  data_comb$.id<-rep(seq_len(card),nrow(data_comb)/card)
  data_orig<-split(data_comb, data_comb$.id, drop = TRUE)

  #Create Co based on the first A of each batch:
  Co_names<-row.names(data_orig[[1]])
  Co_fin<-Co_data[match(Co_names, row.names(Co_data)),]
  names(Co_fin)<-c("original", paste("Co",1:Co, sep ="_"))

  #Make output prettier:
  data_final<-lapply(1:length(data_orig), function(x){
    #Remove unnecessary columns
    data_orig[[x]]<-data_orig[[x]][,-c(1,2,3)]
    data_orig[[x]][sapply(data_orig[[x]], function(t) all(is.na(t)))] <- NULL

    #Add Co columns:
    data_orig[[x]]<-cbind.data.frame(data_orig[[x]],Co_fin[,-1])
  })

  #Simple data transformation for ltmle



  out <- list(Y=data_final[[fY]],A=data_final[iA],final=data_final,original=data_orig,
              step=step,card=card)
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
#' @return An object of class \code{tstmle}.
#' \describe{
#' \item{valY}{Prediction based on the input \code{fit} object.}
#' }
#'
#' @importFrom sl3 make_sl3_Task
#'
#' @export
#'

initEst_setA<-function(data, fit, setA){

  data_setA<-data

  #Set A node to specified intervention:
  data_setA[,"A"] = setA

  covars<-names(data_setA[,-1])
  outcome<-names(data_setA[,1])

  #Create sl3 task:
  #TO DO: add option for weights
  task_setA <- sl3::make_sl3_Task(data_setA, covariates = covars, outcome = outcome,
                                  outcome_type=NULL, folds=fit$folds)

  sl_preds <- fit$fullFit$predict(task_setA)

  out <- list(valY = sl_preds)

  return(out)

}

#' Split-specific predictions
#'
#' Function to estimate split-specific predictions for Q and g, and relevant parts for the
#' optimal treatment regime.
#'
#' @param fold one of the folds from the folds object.
#' @param Q data.frame object for Q containing the time series with C_o as the relevant covariates.
#' @param g data.frame object for g containing the time series with C_o as the relevant covariates.
#' @param estQ object of class \code{initEst} with \code{sl3} results for Q.
#' @param estg object of class \code{initEst} with \code{sl3} results for g.
#'
#' @return An object of class \code{tstmle}.
#' \describe{
#' \item{Y}{Observed outcomes (Y).}
#' \item{A}{Observed exposure (A)}
#' \item{QAW}{Split-specific prediction of Q.}
#' \item{Q1W}{Split-specific prediction of Q(1,W).}
#' \item{Q0W}{Split-specific prediction of Q(0,W).}
#' \item{pA1}{Split-specific prediction of g.}
#' \item{B}{Split-specific prediction of the blip.}
#' \item{Rule}{Split-specific prediction of the maximizing rule.}
#' \item{K}{Split-specific absolute value of the blip.}
#' }
#'
#' @importFrom sl3 make_sl3_Task
#'
#' @export
#'

cv_split <- function(fold, Q, g, estQ, estg){

  #Observed data:
  Y<-Q[,1]
  A<-Q[,2]

  #Get fold index
  v <- origami::fold_index()

  #Get the coefficients:
  coef<-estQ$fullFit$coefficients
  coefg<-estg$fullFit$coefficients

  #Get split-specific fits:
  splitQ <- estQ$cvFit[[v]]
  splitg <- estg$cvFit[[v]]

  #Predict on all the data:
  new_folds<-origami::make_folds(Q, fold_fun = origami::folds_resubstitution)[[1]]

  #QAW
  covars<-names(Q[,-1])
  outcome<-names(Q[,1])

  task<-sl3::make_sl3_Task(Q, covariates = covars, outcome = outcome,
                           outcome_type=NULL, folds=new_folds)
  QAW <- as.matrix(splitQ$predict(task)) %*% as.matrix(coef)

  #Q1W:
  new_data<-Q
  new_data[,"A"] = 1

  task<-sl3::make_sl3_Task(new_data, covariates = covars, outcome = outcome,
                           outcome_type=NULL, folds=new_folds)
  Q1W <- as.matrix(splitQ$predict(task)) %*% as.matrix(coef)

  #Q0W:
  new_data<-Q
  new_data[,"A"] = 0

  task<-sl3::make_sl3_Task(new_data, covariates = covars, outcome = outcome,
                           outcome_type=NULL, folds=new_folds)
  Q0W <- as.matrix(splitQ$predict(task)) %*% as.matrix(coef)

  #p1A:
  covars<-names(g[,-1])
  outcome<-names(g[,1])

  task<-sl3::make_sl3_Task(g, covariates = covars, outcome = outcome,
                           outcome_type=NULL, folds=new_folds)
  pA1 <- as.matrix(splitg$predict(task)) %*% as.matrix(coefg)

  #Blip
  A<-Q$A
  Y<-Q$Y

  B <- (A/pA1 - (1 - A)/(1 - pA1)) * (Y - QAW) + Q1W - Q0W

  #Maximize:
  Rule <- as.numeric(B > 0)

  K <- as.vector(abs(B))

  return(list(Y=Y, A=A, QAW = QAW, Q1W = Q1W, Q0W = Q0W, pA1 = pA1, B = B, Rule = Rule, K = K))

}

#' Split-specific blip predictions
#'
#' Function to estimate split-specific predictions for the blip function.
#'
#' @param fold one of the folds from the folds object.
#' @param estSplit result of the cross-validated \code{cv_split} function.
#' @param data data.frame object containing the time series with C_o as the relevant covariates.
#' @param SL.library list of \code{sl3} algorithms for the fit of the blip function.
#'
#' @return An object of class \code{tstmle}.
#' \describe{
#' \item{fullPred}{Split-specific blip with full V-fold cross-validation and Super Learner results.}
#' \item{cvPred}{Split-specific blip with prediction for only a specific fold.}
#' \item{valSet}{Validation samples used for each fold.}
#' \item{B}{Estimated blip for the validation samples.}
#' }
#'
#' @importFrom sl3 make_sl3_Task
#' @importFrom origami fold_index
#'
#' @export
#'

cv_split_blip<-function(fold,estSplit,data,SL.library,outcome_type){

  v <- origami::fold_index()

  #Get the corresponding Blip for the fold v:
  B<-data.frame(B=estSplit$B[[v]])
  folds<-estSplit$folds

  res<-initEst(Y = B, X =  data[,-1], folds = folds, SL.library = SL.library, outcome_type=outcome_type)

  #Grab just the v fit prediction:
  valset<-data[fold$validation_set,]

  task<-sl3::make_sl3_Task(valset, covariates = names(valset[,-1]), outcome = names(valset[,1]),
                           outcome_type=outcome_type, folds=fold)

  valset_res<-res$cvFit[[v]]$predict(task)
  row.names(valset_res)<-fold$validation_set

  return(list(fullPred=res,cvPred=valset_res,valSet=fold$validation_set,B=B[fold$validation_set,]))
}
