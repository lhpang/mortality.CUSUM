# Mortality-CUSUM

`mortality-CUSUM` is an R package for drawing CUSUM plots to monitor a medical center's survival outcomes. CUSUM (CUmulative SUMmation) is a real-time monitor for CUmulative SUMmation of risk-adjusted survival outcomes. We provide O-E CUSUM, one-sided CUSUM, V-mask and reversed CUSUM plots. The resulting plot is designed to simultaneously monitor for failure time outcomes that are worse than expected or better than expected.

# Introduction
(motivation-current problem-what we did)
Large-scale time-to-event data derived from national disease registries arise rapidly in medical studies. Detecting and accounting for time-varying effects is particularly important, as time-varying effects have already been reported in the clinical literature. However, there are currently no formal R packages for estimating the time-varying effects without pre-assuming the time-dependent function. Inaccurate pre-assumptions can greatly influence the estimation, leading to unreliable results. To address this issue, we developed a time-varying model using spline terms with penalization that does not require pre-assumption of the true time-dependent function, and implemented it in R.

(benefits, computational burden)
Our package offers several benefits over traditional methods. Firstly, traditional methods for modeling time-varying survival models often rely on expanding the original data into a repeated measurement format. However, even with moderate sample sizes, this leads to a large and computationally burdensome working dataset. Our package addresses this issue by proposing a computationally efficient Kronecker product-based proximal algorithm, which allows for the evaluation of time-varying effects in large-scale studies. Additionally, our package allows for parallel computing and can handle moderate to large sample sizes more efficiently than current methods.

(what problem of current algorithm we modified)
In our statistical software tutorial, we address a common issue encountered when analyzing data with binary covariates with near-zero variation. For example, in the SEER prostate cancer data, only 0.6% of the 716,553 patients had their tumors regional to the lymph nodes. In such cases, the associated observed information matrix of a Newton-type method may have a minimum eigenvalue close to zero and a large condition number. Inverting this nearly singular matrix can lead to numerical instability and the corresponding Newton updates may be confined within a small neighborhood of the initial value, resulting in estimates that are far from the optimal solutions. To address this problem, our proposed Proximal-Newtown method utilizes a modified Hessian matrix, which allows for accurate estimation in these scenarios.
# Methods

O-E CUSUM
one-sided CUSUM
V-mask
reversed CUSUM plots

# Detailed Options
risk-adjusted
grouped
dynamic control limit

# Usage

# Datasets
## Simulated Datasets
## Real Datasets
## Installation

You can install the development version of mypackage from [GitHub](https://github.com/) with:

``` r
# install.packages("devtools")
devtools::install_github("lhpang/mortality-CUSUM")
```

``` r
library(mortality-CUSUM)
## basic example code
```
# Detailed tutorial
# Getting Help
# References
