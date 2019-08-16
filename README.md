# safeBart

<!-- badges: start -->
<!-- badges: end -->

The goal of safeBart is to to provide implementations of Bayesian Additive Regression Trees and Bayesian Causal Forest using importance sampling. This can be viewed as an importance sampling approach to BART-BMA (Hernandez et al. 2018). The data-independent sampling from the prior and use of the power likelihood follows the safe-Bayesian Random Forest method described by Quadrianto and Ghahramani (2015).

Hernández, B., Raftery, A. E., Pennington, S. R., & Parnell, A. C. (2018). Bayesian additive regression trees using Bayesian model averaging. Statistics and computing, 28(4), 869-890.

Quadrianto, N., & Ghahramani, Z. (2014). A very simple safe-Bayesian random forest. IEEE transactions on pattern analysis and machine intelligence, 37(6), 1297-1303.

## Installation

``` r
library(devtools)
install_github("EoghanONeill/safeBart")
```

## Example


``` r
library(safeBart)
beta_par <- 0.5

N <- 100
p<- 5
set.seed(100)

epsilon <- rnorm(N)

xcov <- matrix(runif(N*p), nrow=N)

y <- sin(pi*xcov[,1]*xcov[,2]) + 20*(xcov[,3]-0.5)^2+10*xcov[,4]+5*xcov[,5]+epsilon
# <- rep(1,N) + epsilon

epsilontest <- rnorm(N)

xcovtest <- matrix(runif(N*p), nrow=N)
ytest <- sin(pi*xcovtest[,1]*xcovtest[,2]) + 20*(xcovtest[,3]-0.5)^2+10*xcovtest[,4]+5*xcovtest[,5]+epsilontest
#ytest <- rep(1,N) + epsilontest




Num_split_vars <- 10

lambda <- 0.45
Num_models <- 10000
num_trees1 <- 5

seed1 <- 42
ncores <- 7



examplepreds1 <- safeBart_parallel(seed1,
  y, xcov,xcovtest,
  lambda=0.45,
  num_models=Num_models,
  num_trees=num_trees1,
  beta_par=beta_par,
  ncores=ncores,
  outsamppreds=1,
  nu=3,
  a=3,
  sigquant=0.9)

cbind(examplepreds1,ytest )
```

