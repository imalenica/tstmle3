#library(tstmle)
library(tensorflow)
library(keras)
library(kerasR)
library(reticulate)
library(sl3)
library(origami)

#Can we use python modules?

reticulate::py_module_available("keras.models")
reticulate::import("keras.models")

#Load data:
load("~/Dropbox/Berkeley_Projects/Software/tstmle/data/sim_ts_s1.rda")
load("~/Dropbox/Berkeley_Projects/Software/tstmle/data/sim_ts_s1_n50.rda")

data<-sim_ts_s1_n50
fitQ=TRUE
j=1
folds=NULL
skip=0
fold_fn="folds_rolling_origin"
window=NULL
SL.library <- list("Lrnr_mean", list("Lrnr_arima", order=c(3,0,4)), "Lrnr_expSmooth")

SL.library <- list("Lrnr_mean", "Lrnr_expSmooth", list("Lrnr_arima", order=c(2,0,2)),
                   list("Lrnr_lstm", epochs = 50))
fold_fn="folds_rolling_window"
