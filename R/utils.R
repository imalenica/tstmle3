#' Make a sl3 stack
#'
#' Function to make passing algorithms to sl3 easier.
#'
#' @param learner_lists list of algorithms, possibly with additional parameters.
#'
#' @importFrom sl3 make_learner
#'
#

make_stack <- function(learner_lists) {
  #learner_lists <- list(...)

  learners <- lapply(learner_lists, function(learner_list) {
    if (!is.list(learner_list)) {
      learner_list <- list(learner_list)
    }

    learner <- do.call(sl3::make_learner, learner_list)
    return(learner)
  })

  # generate output Stack object and return model stack
  out <- make_learner(Stack, learners)
  return(out)
}

#' Fit sl3 Super Learner
#'
#' Function to make Super Learning with \code{sl3} easier.
#'
#' @param task task object of \code{sl3::make_sl3_Task} form.
#' @param SL.library list of \code{sl3} algorithms to be used for estimation.
#'
#' @importFrom sl3 make_learner Lrnr_sl
#'
#' @export
#

sl3.fit<-function(task,SL.library){

  #Create a stack:
  sl_stack <- make_stack(learner_lists=SL.library)

  #Train and generate predictions on the validation set:
  #cv_stack <- sl3::Lrnr_cv$new(sl_stack)
  #cv_fit <- cv_stack$train(task)
  #cv_preds <- cv_fit$predict(task_pred)
  #risks <- cv_fit$cv_risk(loss_squared_error)

  #cv_fit<-cv_fit$fit_object$fold_fits

  #Super Learner:
  metalearner <- sl3::make_learner(sl3::Lrnr_nnls)
  sl <- sl3::Lrnr_sl$new(learners = sl_stack, metalearner = metalearner)
  sl_fit <- sl$train(task)

  #Save risk:
  risk<-sl_fit$cv_risk(sl3::loss_squared_error)

  #Save predictions:
  sl_preds <- sl_fit$predict()

  #Save CV fits:
  cv_fit<-sl_fit$fit_object$cv_fit$fit_object$fold_fits

  return(list(pred=sl_preds, risk=risk, coefs=sl_fit$coefficients, sl=sl,
              sl.fit=sl_fit, cv.fit=cv_fit))
}

#' Bounds for Q and g
#'
#' Function to bound Q and g estimates with user specified bounds.
#'
#' @param x observed values of either Q or g.
#' @param bds list containing the upper and lower bound for x.
#'
#' @export
#'

bound <- function(x, bds){
  x[x > max(bds)] <- max(bds)
  x[x < min(bds)] <- min(bds)
  x
}

#' Cross-validate and extract validation samples
#'
#' Function to extract and cross-validate validation samples from prediction on all
#' samples for use in CV-TMLE.
#'
#' @param folds user-specified list of folds. It should correspond to an element of \code{origami}.
#' @param split_preds Cross-validated result of \code{cv_split}.
#'
#'

extract_vals <- function(folds, split_preds) {

  val_preds <- origami::cross_validate(extract_val, folds, split_preds)$preds
  val_preds <- val_preds[order(val_preds$index), ]

}

#' Extract validation samples
#'
#' Function to extract the validation sets from split based predictions for use in CV-TMLE.
#'
#' @param fold one fold from a list of folds.
#' @param split_preds Cross-validated result of \code{cv_split}.
#'

extract_val <- function(fold, split_preds) {

  #Extract fold and corresponding validation samples.
  v <- origami::fold_index()
  valid_idx <- origami::validation()

  val_preds <- sapply(split_preds, function(split_pred) {
    split_pred[[v]][valid_idx]
  })

  #Add index to it (aka, which sample is in question)
  val_preds <- as.data.frame(val_preds)
  val_preds$index <- valid_idx
  result <- list(preds = val_preds)
  return(result)
}

#' Wrapper for split-specific predictions
#'
#' Function to make split-specific predictions more streamlined.
#'
#' @param folds user-specified list of folds. It should correspond to an element of \code{origami}.
#' @param Q data.frame containg the relevant outcome and covariates for estimating Q.
#' @param g data.frame containg the relevant outcome and covariates for estimating g.
#' @param estQ Q result of \code{initEst} format.
#' @param estg g result of \code{initEst} format.
#'
#' @export
#

estSplit<-function(folds, Q, g, estQ, estg){

  estSplit <- origami::cross_validate(cv_split, folds, Q, g, estQ, estg, .combine = F)
  estSplit$errors<-NULL
  estSplit$folds<-folds

  estSplit_val<-extract_vals(folds, estSplit)

  return(list(estSplit=estSplit,valSplit=estSplit_val))
}

#' Wrapper for split-specific Blip
#'
#' Function to make split-specific blip estimation more streamlined.
#'
#' @param folds user-specified list of folds- it should correspond to an element of \code{origami}.
#' @param Q data.frame containg the relevant outcome and covariates for estimating Q.
#' @param estSplt result object of \code{estSplit}.
#' @param blip_library list of \code{sl3} algorithms for the fit of the blip function.
#'
#' @importFrom stats coef
#' @importFrom nnls nnls
#'
#' @export
#'

estBlip<-function(folds, estSplt, Q, blip_library){

  blipSplit<-origami::cross_validate(cv_split_blip, folds, estSplt$estSplit, Q, blip_library,
                                     outcome_type = "continuous", .combine = F)

  #Construct SL prediction of the blip:
  #For now, just use non-negative linear least squares.
  x<-do.call(rbind, blipSplit$cvPred)
  y<-unlist(blipSplit$B)

  fit_coef <- stats::coef(nnls::nnls(as.matrix(x), as.matrix(y)))
  fit_coef <- fit_coef/sum(fit_coef)

  pred<-as.matrix(x) %*% fit_coef

  #Treatment from Blip based on MSE:
  #TO DO: add other options.
  optA<-as.numeric(pred > 0)

  return(list(blip=pred,optA=optA,coef=fit_coef,blipSplit=blipSplit,x=x,y=y))

}

#' Wrapper for final TMLE results
#'
#' Function to make extraction of final TMLE results for the optimal individualized
#' treatment regime parameter more streamlined. In particular, it should provide inference and
#' results for 4 different scenarious, where the exposure is set to both binary possibilities,
#' observed exposure, and learner optimal rule.
#'
#' @param res results from multiple calls to \code{ruletmle}.
#'
#' @export
#'

extract_res<-function(res){

  #Extract Psi for each rule:
  psi <- data.frame(unlist(lapply(seq_len(4), function(x) {res[[x]]$tmlePsi})))
  row.names(psi)<-c("A=0","A=1","A=A","A=optA")
  names(psi)<-"Psi"

  #Extract SD for each rule:
  sd <- data.frame(unlist(lapply(seq_len(4), function(x) {res[[x]]$tmleSD})))
  row.names(sd)<-c("A=0","A=1","A=A","A=optA")
  names(sd)<-"SD"

  #Extract CI for each rule:
  lower <- data.frame(unlist(lapply(seq_len(4), function(x) {res[[x]]$tmleCI[1]})))
  upper <- data.frame(unlist(lapply(seq_len(4), function(x) {res[[x]]$tmleCI[2]})))
  CI<-cbind.data.frame(lower=lower,upper=upper)
  names(CI)<-c("lower","upper")
  row.names(CI)<-c("A=0","A=1","A=A","A=optA")

  #Extract IC for each rule:
  IC<-lapply(seq_len(4), function(x) {res[[x]]$IC$Dstar_psi})
  IC<-data.frame(do.call("cbind", IC))
  names(IC)<-c("A=0","A=1","A=A","A=optA")

  #Extract rule:
  rule<-lapply(seq_len(4), function(x) {res[[x]]$rule})
  rule<-data.frame(do.call("cbind", rule))
  names(rule)<-c("A=0","A=1","A=A","A=optA")

  #Extract number of steps until convergence:
  steps <- data.frame(unlist(lapply(seq_len(4), function(x) {res[[x]]$steps})))
  row.names(steps)<-c("A=0","A=1","A=A","A=optA")
  names(steps)<-"steps"

  #Extract initial data:
  initialData<-lapply(seq_len(4), function(x) {res[[x]]$initialData})
  names(initialData)<-c("rule0","rule1","ruleA","ruleOpt")

  #Extract final data:
  tmleData<-lapply(seq_len(4), function(x) {res[[x]]$tmleData})
  names(tmleData)<-c("rule0","rule1","ruleA","ruleOpt")

  #Extract all results:
  all<-lapply(seq_len(4), function(x) {res[[x]]$all})
  names(all)<-c("rule0","rule1","ruleA","ruleOpt")

  return(list(tmlePsi=psi,tmleSD=sd,tmleCI=CI,IC=IC,rule=rule,steps=steps,initialData=initialData,
              tmleData=tmleData,all=all))

}
