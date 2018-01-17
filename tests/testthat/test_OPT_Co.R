library(sl3)
library(origami)

load("~/Dropbox/Berkeley_Projects/Software/tstmle/data/sim_ts_s1_n50.rda")

#Set library:
Q_library=list("Lrnr_mean", "Lrnr_glm_fast", "Lrnr_glmnet","Lrnr_randomForest","Lrnr_xgboost")
g_library=list("Lrnr_mean", "Lrnr_glm_fast", "Lrnr_glmnet","Lrnr_randomForest","Lrnr_xgboost")
blip_library=list("Lrnr_glm_fast", "Lrnr_glmnet","Lrnr_randomForest","Lrnr_xgboost", "Lrnr_nnls")

test_that("Single Intervention Context-Specific ATE with default parameters works", {

  res<-tstmleOPT(sim_ts_s1_n50, Cy=6, Ca=5, V=3, Q_library, g_library, blip_library, folds=NULL,
                stratifyAY = TRUE, gbounds=c(1e-4,1-1e-4), Qbounds=c(1e-4,1-1e-4), maxIter=1000)

  expect_equal(res$tmlePsi[1,], 0.5427367, tolerance = 0.01)
  expect_equal(res$tmleSD[1,], 0.1068668, tolerance = 0.01)
  expect_equal(res$steps, 1, tolerance = 0.01)

})