#' Context-specific causal effect of single-time point intervention
#'
#' This function estimates the causal effect of a single time-point intervention on the outcome
#' within the same time point. Here, time is defined as observing a single intervention and outcome,
#' and possibly multiple covariates. Intervention is imposed within each unit of time.
#'
#' @param data data.frame object containing the data.
#' @param node_list node list reflecting the relationships between variables.
#' @param parameter target parameter for targeting. Default is the average over time
#' Context-Specific Average Treatment Effect.
#' @param learner_list learner list containing the learners used for the conditional expectation
#' of outcome and propensity score.
#' @param Co user-specified Markov order for the fixed dimensional summary measure.
#' @param cvtmle default is \code{FALSE}.
#' @param fold_fn cross-validation scheme, as defined by \code{origami}. See \code{?origami::fold_funs}
#'  for detailed explanations. For time-series, implemented cross-validation schemes are
#' \code{folds_rolling_origin} and \code{folds_rolling_window}.
#' @param first_window first window size used for training. Only relevant if the set cross-validation is
#' \code{folds_rolling_origin}.
#' @param validation_size number of time points used for validation.
#' @param gap number of time points between training and validation set.
#' @param batch number of time points added in the next next cross-validation fold.
#' @param window_size window size used for training. Only relevant if the set cross-validation is
#' \code{folds_rolling_window}.
#'
#'
#' @return An object of class \code{tstmle}.
#' \describe{
#' \item{tmlePsi}{Average treatment effect estimated using TMLE.}
#' \item{iptwPsi}{Average treatment effect estimated using IPTW.}
#' \item{tmleSE}{Standard error for the TMLE estimated parameter.}
#' \item{tmleSD}{Standard deviation for the TMLE estimated parameter.}
#' \item{tmleCI}{Confidence Interval for the TMLE estimated parameter.}
#' \item{IC}{Influence function.}
#' \item{steps}{Number of steps until convergence of the iterative TMLE.}
#' \item{initialData}{Initial estimates of g, Q.}
#' \item{tmleData}{Final updates estimates of g, Q and clever covariates.}
#' \item{tmle_fit}{The full tmle3 fit object.}
#' }
#'
#' @importFrom stats qnorm
#'
#' @export
#

tstmle3 <- function(data, node_list, parameter="ATE",
                    learner_list, Co=5, cvtmle=FALSE,
                    fold_fn="folds_rolling_origin",
                    #Params for folds_rolling_origin
                    first_window=100,validation_size=50,gap=0,batch=50,
                    #Params for folds_rolling_window
                    window_size=100){

  #Additional params
  n <- nrow(data)
  data <- data.frame(data)

  A <- node_list$A
  Y <- node_list$Y
  W <- node_list$W
  X <- node_list$X

  Y_original <- data[,Y]

  #Make folds
  if(fold_fn=="folds_rolling_origin"){
    folds <- origami::make_folds(folds_rolling_origin,
                                 n=n,
                                 first_window=first_window,
                                 validation_size=validation_size,
                                 gap=gap,
                                 batch=batch)
  }else if(fold_fn=="folds_rolling_window"){
    folds <- origami::make_folds(folds_rolling_window,
                                 n=n,
                                 window_size=window_size,
                                 validation_size=validation_size,
                                 gap=gap,
                                 batch=batch)
  }else{
    folds <- origami::make_folds(folds_vfold, n=n,
                                 V=10)
  }

  if(cvtmle){
    folds <- origami::make_folds(folds_vfold, n=n,
                                 V=10)
  }

  if(parameter=="ATE"){
    #Initialized the ATE Spec
    ts_spec <- tmle3_Spec_ts_ATE$new(treatment_level = 1,
                                     control_level = 0,
                                     Co=Co)

    #Initialize the tasks
    tmle_task <- ts_spec$make_tmle_task(data, node_list, folds=folds)

    A_task <- tmle_task$get_regression_task("A")
    Y_task <- tmle_task$get_regression_task("Y")

    A <- A_task$Y
    Y <- Y_task$Y

    initial_likelihood <- ts_spec$make_initial_likelihood(
      tmle_task,
      learner_list
    )

    if(cvtmle){
      updater <- ts_spec$make_updater(cvtmle=TRUE)
    }else{
      updater <- ts_spec$make_updater(cvtmle=FALSE)
    }

    targeted_likelihood <- ts_spec$make_targeted_likelihood(initial_likelihood, updater)
    tmle_params <- ts_spec$make_params(tmle_task, likelihood=targeted_likelihood)

    #TMLE:
    tmle_fit <- fit_tmle3(tmle_task, targeted_likelihood,
                          tmle_params, updater)


    g <- initial_likelihood$get_likelihood(tmle_task, node = "A")
    g_all<-ifelse(A==1,g,1-g)

    #IPTW:
    iptwEst <- mean(as.numeric(A==1)*Y_original/g_all) -
      mean(as.numeric(A==0)* Y_original/g_all)

    #One-step:
    EIC         <- tmle_fit$estimates[[1]]$IC
    initial_psi <- tmle_fit$initial_psi
    oneStepEst  <- (initial_psi + mean(EIC))
  }

  tmlePsi <- tmle_fit$estimates[[1]]$psi
  tmleSE  <- tmle_fit$summary$se
  lower   <- tmle_fit$summary$lower
  upper   <- tmle_fit$summary$upper

  return(list(tmlePsi=tmlePsi, iptwPsi=iptwEst, onestepPsi=oneStepEst,
              tmleSE=tmleSE, tmleSD=tmleSE*sqrt(n),tmleCI=c(lower,upper),
              IC=EIC, steps=tmle_fit$steps,
              initialData=tmle_fit$tmle_task$data,
              tmleData=tmle_fit$tmle_task,
              tmle_fit=tmle_fit))
}

