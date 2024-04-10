
<!-- README.md is generated from README.Rmd. Please edit that file -->

# crrf

<!-- badges: start -->
<!-- badges: end -->

The goal of crrf is to enable fitting of the Firth-penalized Fine-Gray
model for competing risks regression. The Firth penalty allows for
parameter estimation in the case of monotone likelihood. This work was
inspired by a SAS macro by Kohl et al., 2016
(<https://doi.org/10.1016/j.cmpb.2014.11.009>).

## Installation

You can install the development version of crrf from
[GitHub](https://github.com/) with:

``` r
# install.packages("devtools")
devtools::install_github("annaSeffernick/crrf")
```

## Example

This is a basic example which shows you how to solve a common problem:

``` r
library(crrf)
# Generate Data
set.seed(12345)
n=100
tm=rexp(n)
ev=factor(rbinom(n,2,0.3))

x=rnorm(n)
y=rbinom(n,1,0.5)
u=runif(n)

dset=cbind.data.frame(tm=tm,ev=ev,x=x,y=y,u=u)
# fit the model
crrf.fit=crrf(Surv(tm,ev)~x+y+u,etype="1",dset,firth=TRUE,CI=TRUE)
crrf.fit$CI.tbl
#>     ml.beta       lower     upper
#> x 0.1656148 -0.16418727 0.4948082
#> y 0.6365446  0.02979408 1.2726049
#> u 0.1264342 -0.94058592 1.2238442
```
