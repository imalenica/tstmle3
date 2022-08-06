#' Defines a TML Estimator (except for the data)
#'
#'
#' @importFrom R6 R6Class
#'
#' @export
#
tmle3_Spec_ts_ATE <- R6Class(
  classname = "tmle3_Spec_ts_ATE",
  portable = TRUE,
  class = TRUE,
  inherit = tmle3_Spec,
  public = list(
    initialize = function(treatment_level, control_level, Co, ...) {
      super$initialize(
        treatment_level = treatment_level,
        control_level = control_level, 
        Co = Co, ...
      )
    },
    make_tmle_task = function(data, node_list, ...) {
      
      #Include lagged values as specified by Co
      Co <- self$options$Co
      
      if(length(node_list$W)==1){
        data_Co <- self$add_summary(Co = Co, 
                                    W  = data.frame(W=data[,node_list$W]),
                                    A  = data.frame(A=data[,node_list$A]),
                                    Y  = data.frame(Y=data[,node_list$Y]),
                                    X  = data.frame(data[,node_list$X]),
                                    combine = TRUE)
      }else if(length(node_list$X)==1){
        data_Co <- self$add_summary(Co = Co, 
                                    W  = data.frame(data[,node_list$W]),
                                    A  = data.frame(A=data[,node_list$A]),
                                    Y  = data.frame(Y=data[,node_list$Y]),
                                    X  = data.frame(X=data[,node_list$X]),
                                    combine = TRUE)
      }else if(length(node_list$W)==1 & length(node_list$X)==1){
        data_Co <- self$add_summary(Co = Co, 
                                    W  = data.frame(W=data[,node_list$W]),
                                    A  = data.frame(A=data[,node_list$A]),
                                    Y  = data.frame(Y=data[,node_list$Y]),
                                    X  = data.frame(X=data[,node_list$X]),
                                    combine = TRUE)
      }else{
        data_Co <- self$add_summary(Co = Co, 
                                    W  = data.frame(data[,node_list$W]),
                                    A  = data.frame(A=data[,node_list$A]),
                                    Y  = data.frame(Y=data[,node_list$Y]),
                                    X  = data.frame(data[,node_list$X]),
                                    combine = TRUE)
      }
      
      
      variable_types <- self$options$variable_types
      setDT(data_Co)
      
      #Create a new node list:
      #node_list_Co <- list(
      #  W = colnames(data_Co)[-(which(names(data_Co)==node_list$A | names(data_Co)==node_list$Y))],
      #  A = node_list$A,
      #  Y = node_list$Y
      #)
      
      node_list_Co <- list(
        W = colnames(data_Co)[-(which(names(data_Co)=="A" | names(data_Co)=="Y"))],
        A = node_list$A,
        Y = node_list$Y
      )
      
      npsem <- point_tx_npsem(node_list_Co, variable_types)
      tmle_task <- tmle3_Task$new(data_Co, npsem = npsem, ...)
      
      return(tmle_task)
    },
    #Generate lagged columns
    add_summary = function(Co,W,A,Y,X,combine=FALSE){
      
      #Get A lags:
      lag_A <- lapply(seq_len(Co), function(x) {
        Hmisc::Lag(A[, 1], x)
      })
      lag_As <- data.frame(matrix(unlist(lag_A), nrow = length(lag_A[[1]])),
                           stringsAsFactors = FALSE)
      names(lag_As)  <- paste0("A_lag_",seq(Co))
      
      #Get X lags:
      lag_X <- names_lag <- list()
      for(i in 1:ncol(X)){
        lag_X[[i]] <- lapply(seq_len(Co), function(x) {
          Hmisc::Lag(X[,i], x)
        })
        #names(lag_X[[i]]) <- paste0(names(X)[i],"_lag_",seq(Co))
        names_lag[[i]] <- paste0(names(X)[i],"_lag_",seq(Co))
      }
      inter <- do.call(cbind,lag_X)
      lag_Xs <- do.call(cbind,inter)
      lag_Xs <- data.frame(lag_Xs)
      
      names_lag <- do.call(cbind, names_lag)
      names(lag_Xs) <- names_lag
      
      #Get Y lags:
      lag_Y <- lapply(seq_len(Co), function(x) {
        Hmisc::Lag(Y[, 1], x)
      })
      lag_Ys <- data.frame(matrix(unlist(lag_Y), nrow = length(lag_Y[[1]])),
                           stringsAsFactors = FALSE)
      names(lag_Ys)  <- paste0("Y_lag_",seq(Co))
      
      #Assemble dataset:
      data_Co <- cbind.data.frame(lag_As,lag_Xs,lag_Ys)
      
      #Instead of NA, have 0
      data_Co[is.na(data_Co)] <- 0
      
      if(combine){
        data_Co <- cbind.data.frame(W,A,Y,X,data_Co)
      }
      
      return(data_Co)
    },
    make_params = function(tmle_task, likelihood) {
      treatment_value <- self$options$treatment_level
      control_value <- self$options$control_level
      A_levels <- tmle_task$npsem[["A"]]$variable_type$levels
      if (!is.null(A_levels)) {
        treatment_value <- factor(treatment_value, levels = A_levels)
        control_value <- factor(control_value, levels = A_levels)
      }
      treatment <- define_lf(LF_static, "A", value = treatment_value)
      control <- define_lf(LF_static, "A", value = control_value)
      ate <- Param_ATE$new(likelihood, treatment, control)
      tmle_params <- list(ate)
      return(tmle_params)
    }
  ),
  active = list(),
  private = list()
)

#' All Treatment Specific Means
#'
#' O=(W,A,Y)
#' W=Covariates
#' A=Treatment (binary or categorical)
#' Y=Outcome (binary or bounded continuous)
#' @importFrom sl3 make_learner Lrnr_mean
#' @param treatment_level the level of A that corresponds to treatment
#' @param control_level the level of A that corresponds to a control or reference level
#' @export
tmle3_ts_ATE <- function(treatment_level, control_level, Co) {
  # TODO: unclear why this has to be in a factory function
  tmle3_Spec_ts_ATE$new(treatment_level, control_level, Co)
}