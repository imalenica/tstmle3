#' Make a sl3 stack
#'
#' Function to make passing algorithms to sl3 easier.
#'
#' @param learner_list list of algorithms, possibly with additional parameters.
#'
#' @importFrom sl3 make_learner
#

make_stack <- function(learner_lists) {
  #learner_lists <- list(...)

  learners <- lapply(learner_lists, function(learner_list) {
    if (!is.list(learner_list)) {
      learner_list <- list(learner_list)
    }

    learner <- do.call(make_learner, learner_list)
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
#' @param SL.library list of \code{sl3} algorithms to be used for estimation. For the list of available
#' learners for time-series, see \code{sl3::sl3_list_learners(c("timeseries"))}.
#'
#' @importFrom sl3 make_learner Lrnr_sl
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
  metalearner <- sl3::make_learner(Lrnr_nnls)
  sl <- sl3::Lrnr_sl$new(learners = sl_stack, metalearner = metalearner)
  sl_fit <- sl$train(task)

  #Save risk:
  risk<-sl_fit$cv_risk(loss_squared_error)

  #Save predictions:
  sl_preds <- sl_fit$predict()

  #Save CV fits:
  cv_fit<-sl_fit$fit_object$cv_fit$fit_object$fold_fits

  return(list(pred=sl_preds, risk=risk, coefs=sl_fit$coefficients, sl=sl,
              sl.fit=sl_fit, cv.fit=cv_fit))
}
