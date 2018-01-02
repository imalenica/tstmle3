#library(tstmle)
library(sl3)
library(origami)

#Load data:
#load("~/Dropbox/Berkeley_Projects/Software/tstmle/data/sim_ts_s1.rda")
load("~/Dropbox/Berkeley_Projects/Software/tstmle/data/sim_ts_s1_n50.rda")

data<-sim_ts_s1_n50
Cy=6
Ca=5
folds=NULL
fold_fn="folds_vfold"
V=3
stratifyAY = TRUE

#To create better library, see:
sl3_list_learners(c("binomial"))
Q_library=list("Lrnr_mean", "Lrnr_glm_fast", "Lrnr_glmnet","Lrnr_randomForest","Lrnr_xgboost")
g_library=list("Lrnr_mean", "Lrnr_glm_fast", "Lrnr_glmnet","Lrnr_randomForest","Lrnr_xgboost")














