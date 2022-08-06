
<!-- README.md is generated from README.Rmd. Please edit that file -->

# R/`tstmle3`

[![License:
MIT](https://img.shields.io/badge/License-MIT-yellow.svg)](https://opensource.org/licenses/MIT)

> Data-adaptive Estimation and Inference for Causal Effects with a
> Single Time Series

**Authors:** [Ivana Malenica](https://github.com/podTockom)

## What’s `tstmle3`?

The `tstmle3` implements robust estimation and provides inference for
data-dependent causal effects based observing a single time series. It’s
an adapter/extension R package in the `tlverse` ecosystem.

Consider the case where one observes a single time-series, denoted as a
single sequence of dependent random variables
![O(1), \\dots O(N)](https://latex.codecogs.com/png.image?%5Cdpi%7B110%7D&space;%5Cbg_white&space;O%281%29%2C%20%5Cdots%20O%28N%29 "O(1), \dots O(N)")
where each
![O(t)](https://latex.codecogs.com/png.image?%5Cdpi%7B110%7D&space;%5Cbg_white&space;O%28t%29 "O(t)")
with
![t \\in \\{1, \\dots ,N\\}](https://latex.codecogs.com/png.image?%5Cdpi%7B110%7D&space;%5Cbg_white&space;t%20%5Cin%20%5C%7B1%2C%20%5Cdots%20%2CN%5C%7D "t \in \{1, \dots ,N\}")
takes values in
![\\mathbf{R}^p](https://latex.codecogs.com/png.image?%5Cdpi%7B110%7D&space;%5Cbg_white&space;%5Cmathbf%7BR%7D%5Ep "\mathbf{R}^p").
Further, we assume that at each time
![t](https://latex.codecogs.com/png.image?%5Cdpi%7B110%7D&space;%5Cbg_white&space;t "t"),
we have a chronological order of the treatment or exposure
![A(t)](https://latex.codecogs.com/png.image?%5Cdpi%7B110%7D&space;%5Cbg_white&space;A%28t%29 "A(t)"),
outcome of interest
![Y(t)](https://latex.codecogs.com/png.image?%5Cdpi%7B110%7D&space;%5Cbg_white&space;Y%28t%29 "Y(t)"),
and possibly other covariates
![W(t)](https://latex.codecogs.com/png.image?%5Cdpi%7B110%7D&space;%5Cbg_white&space;W%28t%29 "W(t)").
While studying time-series data, one might be interested in what the
conditional mean of the outcome would have been had we intervened on one
or more of the treatment nodes in the observed time-series.

The `tstmle3` package focuses on a class of statistical target
parameters defined as the average over time
![t](https://latex.codecogs.com/png.image?%5Cdpi%7B110%7D&space;%5Cbg_white&space;t "t")
of context-specific pathwise differentiable target parameters of the
conditional distribution of the time-series (Malenica and van der Laan
2018b). In particular, it implements several context-specific causal
parameters that can be estimated in a double robust manner and therefore
fully utilize the sequential randomization.

In particular, `tstmle3` implements few different context-specific
parameters:

1.  Average over time of context-specific ATE of a single time point
    intervention.

2.  Average over time of context-specific TSM of a single time point
    intervention.

Here, initial estimation is based on the
[sl3](https://github.com/tlverse/sl3) package, which constructs ensemble
models with proven optimality properties for time-series data (Malenica
and van der Laan 2018a).

------------------------------------------------------------------------

## Installation

You can install a stable release of `tstmle` from GitHub via
[`devtools`](https://www.rstudio.com/products/rpackages/devtools/) with:

``` r
devtools::install_github("imalenica/tstmle3")
```

Note that in order to run `tstmle` you will also need `sl3` and `tmle3`:

``` r
devtools::install_github("tlverse/sl3")
devtools::install_github("tlverse/tmle3")
```

<!--

In the future, the package will be available from
[CRAN](https://cran.r-project.org/) and can be installed via


```r
install.packages("tstmle")
```

-->

------------------------------------------------------------------------

## Issues

If you encounter any bugs or have any specific feature requests, please
[file an issue](https://github.com/imalenica/tstmle3/issues).

------------------------------------------------------------------------

## Citation

After using the tstmle3 R package, please cite the following:

``` r
@software{malenica2022tstmle3,
      author = {Malenica, Ivana and {van der Laan}, Mark J},
      title = {{tstmle3}: Context-Specific Targeted Learning for time-series},
      year  = {2022},
      doi = {},
      url = {https://github.com/imalenica/tstmle3},
      note = {R package version 1.0.0}
    }
```

## License

© 2022 [Ivana Malenica](https://github.com/imalenica)

The contents of this repository are distributed under the MIT license.
See below for details:

    The MIT License (MIT)

    Copyright (c) 2022-2023

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

## References

<div id="refs" class="references csl-bib-body hanging-indent">

<div id="ref-c3" class="csl-entry">

Malenica, Ivana, and Mark J van der Laan. 2018a. “Oracle Inequality for
Cross-Validation Estimator Selector for Dependent Time-Ordered
Experiments.”

</div>

<div id="ref-c2" class="csl-entry">

———. 2018b. “Robust Estimation of Data-Dependent Causal Effects Based on
Observing a Single Time-Series.”

</div>

</div>
