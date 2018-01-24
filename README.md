
<!-- README.md is generated from README.Rmd. Please edit that file -->
R/`tstmle`
==========

[![Travis-CI Build Status](https://travis-ci.org/podTockom/tstmle.svg?branch=master)](https://travis-ci.org/podTockom/tstmle) [![AppVeyor Build Status](https://ci.appveyor.com/api/projects/status/github/podTockom/tstmle?branch=master&svg=true)](https://ci.appveyor.com/project/podTockom/tstmle) [![Coverage Status](https://img.shields.io/codecov/c/github/podTockom/tstmle/master.svg)](https://codecov.io/github/podTockom/tstmle?branch=master) [![Project Status: WIP - Initial development is in progress, but there has not yet been a stable, usable release suitable for the public.](http://www.repostatus.org/badges/latest/wip.svg)](http://www.repostatus.org/#wip) [![License: MIT](https://img.shields.io/badge/License-MIT-yellow.svg)](https://opensource.org/licenses/MIT)

> Data-adaptive Estimation and Inference for Causal Effects with a Single Time Series

**Authors:** [Ivana Malenica](https://github.com/podTockom)

What's `tstmle`?
----------------

The `tstmle` package implements robust estimation and provides inference for data-dependent causal effects based observing a single time series.

Consider the case where one observes a single time-series, denoted as a single sequence of dependent random variables *O*(1),…*O*(*N*) where each *O*(*t*) with *t* ∈ {1, …, *N*} takes values in **R**<sup>*p*</sup>. Further, we assume that at each time *t*, we have a chronological order of the treatment or exposure *A*(*t*), outcome of interest *Y*(*t*), and possibly other covariates *W*(*t*). While studying time-series data, one might be interested in what the conditional mean of the outcome would have been had we intervened on one or more of the treatment nodes in the observed time-series. Additionally, one might also want to learn the optimal treatment rule for the single unit over time.

The `tstmle` package focuses on a class of statistical target parameters defined as the average over time *t* of context-specific pathwise differentiable target parameters of the conditional distribution of the time-series (van der Laan and Malenica 2018). In particular, it implements several context-specific causal parameters that can be estimated in a double robust manner and therefore fully utilize the sequential randomization.

In particular, `tstmle` implements 3 different context-specific parameters:

1.  Average over time of context-specific causal effect of a single time point intervention.

2.  Average over time of context-specific causal effect of multiple time point interventions.

3.  Adaptive design learning the optimal individualized rule within a single time-series.

Here, initial estimation is based on the [sl3](https://github.com/jeremyrcoyle/sl3) package, which constructs ensemble models with proven optimality properties for time-series data (Malenica and van der Laan 2018).

------------------------------------------------------------------------

Installation
------------

You can install a stable release of `tstmle` from GitHub via [`devtools`](https://www.rstudio.com/products/rpackages/devtools/) with:

``` r
devtools::install_github("podTockom/tstmle")
```

Note that in order to run `tstmle` you will also need `sl3`:

``` r
devtools::install_github("jeremyrcoyle/sl3")
```

<!--

In the future, the package will be available from
[CRAN](https://cran.r-project.org/) and can be installed via


```r
install.packages("tstmle")
```

-->

------------------------------------------------------------------------

Examples
--------

To illustrate how `tstmle` may be used to ascertain the effect of an intervention on a single time series, consider the following examples.

#### Context-specific causal effect of single-time point intervention

In this section we utilize a simple, short data-set in order to estimate the causal effect of a single time-point intervention on the next outcome.

``` r
#Load relevant packages:
library(tstmle)
#> Warning in as.POSIXlt.POSIXct(Sys.time()): unknown timezone 'zone/tz/2017c.
#> 1.0/zoneinfo/America/Los_Angeles'
library(sl3)
library(origami)
#> origami: Generalized Cross-Validation Framework
#> Version: 0.8.2

#set seed:
set.seed(12)

#Load the data:
data("sim_ts_s1")

#Set library:
Q_library=list("Lrnr_mean", "Lrnr_glm_fast", "Lrnr_glmnet","Lrnr_randomForest","Lrnr_xgboost")
g_library=list("Lrnr_mean", "Lrnr_glm_fast", "Lrnr_glmnet","Lrnr_randomForest","Lrnr_xgboost")

#Obtain estimates:
res<-tstmleSI(sim_ts_s1, Co=TRUE, stratifyAY = TRUE, folds=NULL, 
              Cy=6, Ca=5, V=10, Q_library, g_library)
#> Fitting Q
#> Fitting g
#> Fitting Q(C_0,0)
#> Fitting Q(C_0,1)

#TMLE:
res$tmlePsi
#>       psi 
#> 0.3140865

#IPTW:
res$iptwPsi
#> [1] 0.314914
```

#### Adaptive design learning the optimal individualized treatment rule

Similarly to the last example, we again use the same short time-series data-set. However, in this example we are interested in adaptive design learning the optimal individualized treatment rule within a single time-series.

``` r
#Load relevant packages:
library(sl3)
library(origami)

#Load the data:
data("sim_ts_s1")

#Set libraries:
Q_library=list("Lrnr_mean", "Lrnr_glm_fast", "Lrnr_glmnet","Lrnr_randomForest","Lrnr_xgboost")
g_library=list("Lrnr_mean", "Lrnr_glm_fast", "Lrnr_glmnet","Lrnr_randomForest","Lrnr_xgboost")
blip_library=list("Lrnr_glm_fast", "Lrnr_glmnet","Lrnr_randomForest","Lrnr_xgboost", "Lrnr_nnls")

#Obtain estimates:
res<-tstmleOPT(sim_ts_s1_n50, Cy=6, Ca=5, V=10, Q_library, g_library, blip_library)

#TMLE:
res$tmlePsi
```

------------------------------------------------------------------------

Issues
------

If you encounter any bugs or have any specific feature requests, please [file an issue](https://github.com/podTockom/tstmle/issues).

------------------------------------------------------------------------

License
-------

© 2018 [Ivana Malenica](https://github.com/podTockom)

The contents of this repository are distributed under the MIT license. See below for details:

    The MIT License (MIT)

    Copyright (c) 2017-2018

    Permission is hereby granted, free of charge, to any person obtaining a copy
    of this software and associated documentation files (the "Software"), to deal
    in the Software without restriction, including without limitation the rights
    to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
    copies of the Software, and to permit persons to whom the Software is
    furnished to do so, subject to the following conditions:

    The above copyright notice and this permission notice shall be included in all
    copies or substantial portions of the Software.

    THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
    IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
    FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
    AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
    LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
    OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE
    SOFTWARE.

------------------------------------------------------------------------

References
----------

Malenica, Ivana, and Mark J van der Laan. 2018. “Oracle Inequality for Cross-Validation Estimator Selector for Dependent Time-Ordered Experiments.”

van der Laan, Mark J, and Ivana Malenica. 2018. “Robust Estimation of Data-Dependent Causal Effects Based on Observing a Single Time-Series.”
