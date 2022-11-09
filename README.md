
<!-- README.md is generated from README.Rmd. Please edit that file -->

# MOTRL

<!-- badges: start -->
<!-- badges: end -->

The goal of MOTRL is to directly estimate optimal dynamic treatment
regime (DTR) or tolerant DTR (tDTR) in a multi-stage multi-treatment
setting. The user can self-define their preference on different
priorities and the tolerant rate at each stages.

## Installation

You can install the development version of MOTRL from
[GitHub](https://github.com/) with:

``` r
# install.packages("devtools")
devtools::install_github("SelinaSong0412/MOTRL")
```

## Example

See below for two examples of using MOTRL to estimating tolerant DTR.

``` r
library(MOTRL)
```

**(a). Bi-objective scenario: one stage, three treatments**

Here we simulate some data of 1000 observeations for 6 variables. The
weights on the two objective are set as (0.7, 0.3).

``` r
# Simulate 6 covariates
set.seed(300)
N = 1000
x1<-rnorm(N) # each covariates follows N(0,1)
x2<-rnorm(N)
x3<-rnorm(N)
x4<-rnorm(N)
x5<-rnorm(N)
x6<-answer(N, x = c("No", "Yes"), name = "Smoke")
X0<-cbind(x1,x2,x3,x4,x5,x6) # All of the covariates

# Simulate observed treatment and optimal treatment  
# observed treatment (A1) with total possible number of treatment = 3
pi10 <- rep(1, N)
pi11 <- exp(0.5*x4 + 0.5*x1)
pi12 <- exp(0.5*x5 - 0.5*x1)
# weights matrix
matrix.pi1 <- cbind(pi10, pi11, pi12)
A1 <- A.sim(matrix.pi1)
class.A1 <- sort(unique(A1))
# propensity stage 1
pis1.hat <- M.propen(A1, cbind(x1,x4,x5))
g1.opt <- (x2 <= 0.5)*(x1 > 0.5) + 2*(x2 > 0.5)*(x1 > -1)

# outcome 1
Y11 <- exp(1.68 + 0.2*x4 - abs(0.5*x1 - 1)*((A1 - g1.opt)^2)) - 3*(A1 == 1) + rnorm(N,0,1)
# outcome 2 
Y12 <- 2.37 + x4 + 2*x5 + 2*(A1 == 0)*(2*(g1.opt == 0) - 1) + 1.5*(A1 == 2)*(2*(g1.opt == 2) -1) + 1.5*(exp((A1 == 1)) -1) + rnorm(N,0,1)
Ys1 = cbind(Y11, Y12)
```

Then, we start to grow the DTR tree and the 60% tolerant DTR tree:

``` r
w0 = 0.7         # set the weight on the first objective as 0.7
wt = c(w0, 1-w0) # the weight vector is (0.7, 0.3)
MOTRL0 = MO.tol.DTRtree(Ys1, w = wt, A1, H=X0, delta = 0, pis.hat=pis1.hat, lambda.pct=0.05, minsplit=max(0.05*N,20),depth = 3)
MOTRL0$tree      # return the MOTRL tree with zero tolerance (only the optimal treatment provided) 
#>   node  X     cutoff      mEy opt.trt  avg.mEy tol.trt avg.mE(Y11) avg.mE(Y12)
#> 1    1  2  0.5169215 4.183206      NA       NA      NA          NA          NA
#> 2    2  1  0.5098781 4.440945      NA       NA      NA          NA          NA
#> 3    3  1 -1.0087522 4.827241      NA       NA      NA          NA          NA
#> 4    4 NA         NA 5.082690       0 5.082690       0    5.542246    4.010394
#> 5    5 NA         NA 3.015530       1 3.015530       1    2.322286    4.633101
#> 6    6 NA         NA 4.832783       0 4.832783       0    4.810765    4.884157
#> 7    7 NA         NA 4.826251       2 4.826251       2    5.201142    3.951505

MOTRL1 = MO.tol.DTRtree(Ys1, w = wt, A1, H=X0, delta = 0.4, pis.hat=pis1.hat, lambda.pct=0.05, minsplit=max(0.05*N,20),depth = 3)
MOTRL1$tree      # return the 60% tolerant DTR tree
#>   node  X     cutoff      mEy opt.trt  avg.mEy tol.trt avg.mE(Y11) avg.mE(Y12)
#> 1    1  2  0.5169215 4.183206      NA       NA      NA          NA          NA
#> 2    2  1  0.5098781 4.440945      NA       NA      NA          NA          NA
#> 3    3  1 -1.0087522 4.827241      NA       NA      NA          NA          NA
#> 4    4 NA         NA 5.082690       0 5.082690       0    5.542246    4.010394
#> 5    5 NA         NA 3.015530       1 3.015530       1    2.322286    4.633101
#> 6    6 NA         NA 4.832783       0 4.832783       0    4.810765    4.884157
#> 7    7 NA         NA 4.826251       2 4.826251       2    5.201142    3.951505

# Run this to get the pseudo outcomes:
# MOTRL1$POs
# Run this to get the loss on pseudo outcome, which is used when generate the DTR tree for the previous stage:
# MOTRL1$PO.loss
```

**(b). Bi-objective scenario: two stage, three treatments**

Using the above simulated data as stage 1 data, we continue simulate a
stage two data:

``` r
# observed stage 2 treatment, A2, with K2 = 3
pi20 <- rep(1,N)
pi21 <- exp(0.2*Y12 - 0.5)
pi22 <- exp(0.2*Y11 - 0.5*x2)
matrix.pi2 <- cbind(pi20, pi21, pi22)
A2 <- A.sim(matrix.pi2)
class.A2 <- sort(unique(A2))
# propensity stage 2
pis2.hat<-M.propen(A2,cbind(Y11,Y12,x2))
# optimal stage 2 treatment
g2.opt <- (x1 <= 0)*((x5 > -0.5) + (x5 > 1.5)) + 2*(x1 > 0)*(x3 > - 1)
  
### stage 2 outcome 1 
Y21 <- 1 + exp(1.5 + 0.2*x2 - abs(1.5*x3 + 2)*(A2 - g2.opt)^2) - 2*(A2 == 1) + rnorm(N,0,1)
### stage 2 outcome 2 
Y22 <- 1.7 + x2 + 2*(A2 == 0)*(2*(g2.opt == 0) - 1) + (A2 == 1)*(2*(g2.opt == 1) - 1) +
  3*(A2 == 2)*(2*(g2.opt == 2) - 1) + 1.2*(exp((A2 == 1)) -1) + rnorm(N,0,1)
Ys2 = cbind(Y21, Y22)
```

We start the estimation from the second stage. Same weight vector (0.7,
0.3) is applied.

``` r
# Stage 2 estimation by backward induction
w21 = c(w0, 1-w0) # same weights 
MOTRL20 = MO.tol.DTRtree(Ys2, w = w21, A2, H=cbind(X0,A1,Ys1), delta = 0, pis.hat=pis1.hat, 
                         lambda.pct=0.05, minsplit=max(0.05*N,20),depth = 3)
MOTRLtree20 = MOTRL20$tree
MOTRLtree20
#>   node  X      cutoff      mEy opt.trt  avg.mEy tol.trt avg.mE(Y21) avg.mE(Y22)
#> 1    1  1 -0.03085405 3.851566      NA       NA      NA          NA          NA
#> 2    2  5 -0.50321659 4.125072      NA       NA      NA          NA          NA
#> 3    3  3 -0.99077712 5.403090      NA       NA      NA          NA          NA
#> 4    4 NA          NA 5.076621       0 5.076621       0    5.598989    3.857762
#> 5    5 NA          NA 3.696732       1 3.696732       1    3.249158    4.741071
#> 6    6 NA          NA 5.115887       0 5.115887       0    5.683258    3.792024
#> 7    7 NA          NA 5.462397       2 5.462397       2    5.801681    4.670734

MOTRL21 <- MO.tol.DTRtree(Ys2, w = w21, A2, H=cbind(X0,A1,Ys1), delta = 0.4, pis.hat=pis1.hat, 
                          lambda.pct=0.05, minsplit=max(0.05*N,20),depth = 3)
MOTRLtree21 = MOTRL21$tree
MOTRLtree21
#>   node  X      cutoff      mEy opt.trt  avg.mEy tol.trt avg.mE(Y21) avg.mE(Y22)
#> 1    1  1 -0.03085405 3.851566      NA       NA      NA          NA          NA
#> 2    2  5 -0.50321659 4.125072      NA       NA      NA          NA          NA
#> 3    3  3 -0.99077712 5.403090      NA       NA      NA          NA          NA
#> 4    4 NA          NA 5.076621       0 5.076621       0    5.598989    3.857762
#> 5    5 NA          NA 3.696732       1 3.696732       1    3.249158    4.741071
#> 6    6 NA          NA 5.115887       0 5.115887       0    5.683258    3.792024
#> 7    7 NA          NA 5.462397       2 5.462397       2    5.801681    4.670734
```

Then, calculate the accumulated mean pseudo outcomes, and use it as the
outcomes (`Ys`).

``` r
# calculate average pseudo outcome 
g2.MOTRLtree20 = predict_tol.DTR(MOTRLtree20, newdata = cbind(X0,A1,Ys1))
g2.MOTRLtree21 = predict_tol.DTR(MOTRLtree21, newdata = cbind(X0,A1,Ys1))
PO.MOTRL20 = Ys1 + MOTRL20$PO.loss
PO.MOTRL21 = Ys1 + MOTRL21$PO.loss

w11 = w21
MOTRL10 = MO.tol.DTRtree(PO.MOTRL20, w = w11, A1, H=X0, delta = 0, pis.hat=pis1.hat, 
                         lambda.pct=0.05, minsplit=max(0.05*N,20),depth = 3)
MOTRL10$tree
#>   node  X     cutoff      mEy opt.trt  avg.mEy tol.trt avg.mE(Y11) avg.mE(Y12)
#> 1    1  2  0.5169215 6.349854      NA       NA      NA          NA          NA
#> 2    2  1  0.5098781 6.657724      NA       NA      NA          NA          NA
#> 3    3  1 -1.0087522 7.986483      NA       NA      NA          NA          NA
#> 4    4 NA         NA 6.849552       0 6.849552       0    6.826217    6.904000
#> 5    5 NA         NA 6.231644       1 6.231644       1    5.911325    6.979053
#> 6    6 NA         NA 7.496037       0 7.496037       0    7.301645    7.949619
#> 7    7 NA         NA 8.074062       2 8.074062       2    8.674228    6.673677
MOTRL11 <- MO.tol.DTRtree(PO.MOTRL21, w = w11, A1, H=X0, delta = 0.4,pis.hat=pis1.hat, 
                            lambda.pct=0.05, minsplit=max(0.05*N,20),depth = 3)
MOTRL11$tree
#>   node  X     cutoff      mEy opt.trt  avg.mEy tol.trt avg.mE(Y11) avg.mE(Y12)
#> 1    1  2  0.5169215 6.349854      NA       NA      NA          NA          NA
#> 2    2  1  0.5098781 6.657724      NA       NA      NA          NA          NA
#> 3    3  1 -1.0087522 7.986483      NA       NA      NA          NA          NA
#> 4    4 NA         NA 6.849552       0 6.849552       0    6.826217    6.904000
#> 5    5 NA         NA 6.231644       1 6.231644       1    5.911325    6.979053
#> 6    6 NA         NA 7.496037       0 7.496037       0    7.301645    7.949619
#> 7    7 NA         NA 8.074062       2 8.074062       2    8.674228    6.673677
```
