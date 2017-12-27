
<!-- README.md is generated from README.Rmd. Please edit that file -->
R/`tstmle`
==========

> Data-adaptive Estimation and Inference for Causal Effects with Single Time Series

**Authors:** Ivana Malenica

What's `tstmle`?
----------------

The `tstmle` package implements targeted maximum likelihood estimation (TMLE) of different causal effects based on the observation of a single time series. We consider the case where we observe a single sequence of dependent random variables *O*(1),…*O*(*n*), where each *O*(*t*) with *t* ∈ {1, …*n*} takes values in **R**<sup>*p*</sup>. Further, we assume that at each time *t*, we have a chronological order of the covariate vector *W*(*t*), treatment or exposure *A*(*t*), and outcome *Y*(*t*).

The `tstmle` package focuses on estimation of target parameters of the conditional distribution of *O*(*t*) given *C*<sub>*o*(*t*)</sub>, where *C*<sub>*o*(*t*)</sub> is some fixed-dimensional summary function of the past *O*(*t* − 1),…*O*(1) (van der Laan and Malenica 2018) In particular, `tstmle` provides estimation and inference for the following:

-   data-dependent, *C*<sub>*o*(*t*)</sub>− specific, causal effect within a single time series

-   adaptive design for learning the optimal treatment rule within a single time series

Here, initial estimation is based on the [sl3](https://github.com/jeremyrcoyle/sl3) package, which constructs ensemble models with proven optimality properties for time-series data (Malenica and van der Laan 2018).

------------------------------------------------------------------------

Installation
------------

You can install a stable release of `tstmle` from GitHub via [`devtools`](https://www.rstudio.com/products/rpackages/devtools/) with:

``` r
devtools::install_github("podTockom/tstmle")
```

<!--

In the future, the package will be available from
[CRAN](https://cran.r-project.org/) and can be installed via


```r
install.packages("tstmle")
```

-->

------------------------------------------------------------------------

Issues
------

If you encounter any bugs or have any specific feature requests, please [file an issue](https://github.com/podTockom/tstmle/issues).

------------------------------------------------------------------------

Example
-------

To illustrate how `tstmle` may be used to ascertain the effect of an intervention on a single time series, consider the following example:

``` r
# Forthcoming...check back later.
```

------------------------------------------------------------------------

License
-------

© 2017 [Ivana Malenica](https://github.com/podTockom)

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
