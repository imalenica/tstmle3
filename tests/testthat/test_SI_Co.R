context("Test Single Intervention Context-Specific ATE")

if (FALSE) {
  setwd("..")
  setwd("..")
  getwd()
  library("devtools")
  document()
  load_all("./")  # load all R files in /R and datasets in /data. Ignores NAMESPACE:
  # devtools::check() # runs full check
  setwd("..")
  install("tstmle", build_vignettes = FALSE, dependencies = FALSE)  # INSTALL W/ devtools:
}

library(testthat)
library(sl3)
library(tstmle)
library(origami)

set.seed(1234)

#load the data
data("sim_ts_s1_n50.rda")

#Set library:
Q_library=list("Lrnr_mean", "Lrnr_glm_fast", "Lrnr_glmnet","Lrnr_randomForest","Lrnr_xgboost")
g_library=list("Lrnr_mean", "Lrnr_glm_fast", "Lrnr_glmnet","Lrnr_randomForest","Lrnr_xgboost")

test_that("Single Intervention Context-Specific ATE with default parameters works", {

  res<-tstmleSI(sim_ts_s1_n50, Co=TRUE, Cy=6, Ca=5, V=3, Q_library, g_library, folds=NULL,
                stratifyAY = TRUE, gbounds=c(1e-4,1-1e-4), Qbounds=c(1e-4,1-1e-4), maxIter=1000)

  expect_equal(res$tmlePsi, 0.3166318, tolerance = 0.01)
  expect_equal(res$iptwPsi, 0.2819559, tolerance = 0.01)
  expect_equal(res$steps, 1, tolerance = 0.01)

})
