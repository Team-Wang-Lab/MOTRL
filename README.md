
<!-- README.md is generated from README.Rmd. Please edit that file -->

# MOT-RL

<!-- badges: start -->
<!-- badges: end -->

## Introduction

The goal of MOT-RL is to directly estimate optimal dynamic treatment
regime (DTR) or tolerant DTR (tDTR) in a multi-stage multi-treatment
setting. The users can self-define their preference on different
priorities and the tolerant rate at each stage.

In this README file, we provide a detailed instruction on the
installation of the MOTRL package, its usage, and the interpretation of
the result.

More example code about the simulation in the manuscript can be found in
the folder ‘vignettes’ at `Simulation_2D.Rmd` (2 objective cases) and
`Simulation_3D.Rmd` (3 objective case).

## Installation

You can install the development version of MOTRL from
[GitHub](https://github.com/) with:

``` r
install.packages("devtools")
devtools::install_github("Team-Wang-Lab/MOTRL")
```

## Examples

See below for two examples of using MOT-RL to estimating tolerant DTR.

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
#> Warning in if (m.method == "AIPW") {: the condition has length > 1 and only the
#> first element will be used
MOTRL0$tree      # return the MOT-RL tree with zero tolerance (only the optimal treatment provided) 
#>   node  X     cutoff      mEy opt.trt  avg.mEy tol.trt avg.mE(Y11) avg.mE(Y12)
#> 1    1  2  0.5169215 4.192301      NA       NA      NA          NA          NA
#> 2    2  1  0.5098781 4.446549      NA       NA      NA          NA          NA
#> 3    3  1 -1.0087522 4.832788      NA       NA      NA          NA          NA
#> 4    4 NA         NA 5.087414       0 5.087414       0    5.542246    4.010394
#> 5    5 NA         NA 3.023088       1 3.023088       1    2.322286    4.633101
#> 6    6 NA         NA 4.855931       0 4.855931       0    4.810765    4.884157
#> 7    7 NA         NA 4.828655       2 4.828655       2    5.201142    3.951505

MOTRL1 = MO.tol.DTRtree(Ys1, w = wt, A1, H=X0, delta = 0.4, pis.hat=pis1.hat, lambda.pct=0.05, minsplit=max(0.05*N,20),depth = 3)
#> Warning in if (m.method == "AIPW") {: the condition has length > 1 and only the
#> first element will be used
MOTRL1$tree      # return the 60% tolerant DTR tree
#>   node  X     cutoff      mEy opt.trt  avg.mEy tol.trt avg.mE(Y11) avg.mE(Y12)
#> 1    1  2  0.5169215 4.192301      NA       NA      NA          NA          NA
#> 2    2  1  0.5098781 4.446549      NA       NA      NA          NA          NA
#> 3    3  1 -1.0087522 4.832788      NA       NA      NA          NA          NA
#> 4    4 NA         NA 5.087414       0 5.087414       0    5.542246    4.010394
#> 5    5 NA         NA 3.023088       1 3.023088       1    2.322286    4.633101
#> 6    6 NA         NA 4.855931       0 4.855931       0    4.810765    4.884157
#> 7    7 NA         NA 4.828655       2 4.828655       2    5.201142    3.951505

# Run this to get the pseudo outcomes:
# MOTRL1$POs
# Run this to get the loss on pseudo outcome, which is used when generate the DTR tree for the previous stage:
# MOTRL1$PO.loss
```

To interpret the result, we take `MOTRL1$tree` as an example.

-   The `node` column of the output represents the node number, the
    number is marked from top to bottom and from left to right in a
    binary split tree graph. Examples of how to number nodes can be
    found in Figure 4 in Section 5 of the manuscript.

-   the `X` column represents the tailoring variable (e.g. the 2 in the
    first row means $X_2$ is the tailoring variable at node 1),

-   the `cutoff` column represents the cutoff value in this tailoring
    variable for node splitting (e.g. for node 1, subjects with
    $X_2 < 0.5169215$ goes to node 2, otherwise, goes to node 3),

-   the `mEy` column represents the weighted counterfactual mean outcome
    for all subjects in this node,

-   the `opt.trt` column represents the optimal treatment found by
    MOT-RL in the node (only for leaves),

-   the `ave.mEy` column represents the weighted counterfactual mean
    outcome for all subjects in this node (same as `mEy` but only for
    leaves),

-   the `tol.trt` column represents the tolerant treatment found by
    MOT-RL in the node (only for leaves),

-   the `ave.mEy(*)` column represents the counterfactual mean outcome of
    objective (\*) for all subjects in this node.

**(b). Bi-objective scenario: two stages, three treatments**

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
#> Warning in if (m.method == "AIPW") {: the condition has length > 1 and only the
#> first element will be used
MOTRLtree20 = MOTRL20$tree
MOTRLtree20
#>   node  X      cutoff      mEy opt.trt  avg.mEy tol.trt avg.mE(Y21) avg.mE(Y22)
#> 1    1  1 -0.03085405 5.129994      NA       NA      NA          NA          NA
#> 2    2 NA          NA 4.089238       1 4.089238       1    2.301012    4.032908
#> 3    3  3 -0.99077712 7.019420      NA       NA      NA          NA          NA
#> 6    6 NA          NA 4.705674       0 4.705674       0    5.683258    3.792024
#> 7    7 NA          NA 7.497201       2 7.497201       2    5.801681    4.670734

MOTRL21 <- MO.tol.DTRtree(Ys2, w = w21, A2, H=cbind(X0,A1,Ys1), delta = 0.4, pis.hat=pis1.hat, 
                          lambda.pct=0.05, minsplit=max(0.05*N,20),depth = 3)
#> Warning in if (m.method == "AIPW") {: the condition has length > 1 and only the
#> first element will be used
MOTRLtree21 = MOTRL21$tree
MOTRLtree21
#>   node  X      cutoff      mEy opt.trt  avg.mEy tol.trt avg.mE(Y21) avg.mE(Y22)
#> 1    1  1 -0.03085405 5.129994      NA       NA      NA          NA          NA
#> 2    2 NA          NA 4.089238       1 4.089238       1    2.301012    4.032908
#> 3    3  3 -0.99077712 7.019420      NA       NA      NA          NA          NA
#> 6    6 NA          NA 4.705674       0 4.705674       0    5.683258    3.792024
#> 7    7 NA          NA 7.497201       2 7.497201       2    5.801681    4.670734
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
#> Warning in if (m.method == "AIPW") {: the condition has length > 1 and only the
#> first element will be used
MOTRL10$tree
#>    node  X     cutoff      mEy opt.trt  avg.mEy tol.trt avg.mE(Y11) avg.mE(Y12)
#> 1     1  2  0.5169215 5.793848      NA       NA      NA          NA          NA
#> 2     2  1  0.5098781 6.056146      NA       NA      NA          NA          NA
#> 3     3  1 -1.0087522 7.318712      NA       NA      NA          NA          NA
#> 4     4 NA         NA 5.998375       0 5.998375       0    5.753969    6.529313
#> 5     5 NA         NA 6.184466       1 6.184466       1    5.816791    6.966871
#> 6     6 NA         NA 6.042648       0 6.042648       0    5.557057    7.287680
#> 7     7  2  1.6684238 8.212634      NA       NA      NA          NA          NA
#> 14   14 NA         NA 8.264803       2 8.264803       2    8.633874    7.370360
#> 15   15 NA         NA 7.978385       0 7.978385       0    9.647820    3.933660
MOTRL11 <- MO.tol.DTRtree(PO.MOTRL21, w = w11, A1, H=X0, delta = 0.4,pis.hat=pis1.hat, 
                            lambda.pct=0.05, minsplit=max(0.05*N,20),depth = 3)
#> Warning in if (m.method == "AIPW") {: the condition has length > 1 and only the
#> first element will be used
MOTRL11$tree
#>    node  X     cutoff      mEy opt.trt  avg.mEy tol.trt avg.mE(Y11) avg.mE(Y12)
#> 1     1  2  0.5169215 5.793848      NA       NA      NA          NA          NA
#> 2     2  1  0.5098781 6.056146      NA       NA      NA          NA          NA
#> 3     3  1 -1.0087522 7.318712      NA       NA      NA          NA          NA
#> 4     4 NA         NA 5.998375       0 5.998375       0    5.753969    6.529313
#> 5     5 NA         NA 6.184466       1 6.184466       1    5.816791    6.966871
#> 6     6 NA         NA 6.042648       0 6.042648       0    5.557057    7.287680
#> 7     7  2  1.6684238 8.212634      NA       NA      NA          NA          NA
#> 14   14 NA         NA 8.264803       2 8.264803       2    8.633874    7.370360
#> 15   15 NA         NA 7.978385       0 7.978385       0    9.647820    3.933660
```

**(c). Tri-objective scenario: one stage, three treatments**

``` r
N<-1000 # sample size of training data
N2<-1000 # sample size of test data
iter <- 5 # replication
w1 = 0.35
w2 = 0.35
w3 = 1 - w1 - w2

# Simulation begin
set.seed(300)
x1<-rnorm(N)              # each covariates follows N(0,1)
x2<-rnorm(N)
x3<-rnorm(N)
x4<-rnorm(N)
x5<-rnorm(N)
x6<-rnorm(N)              # each covariates follows N(0,1)
x7<-rnorm(N)
x8<-rnorm(N)
x9<-rnorm(N)
x10<-answer(N, x = c("No", "Yes"), name = "Smoke")
  
X0<-cbind(x1,x2,x3,x4,x5,x6,x7,x8,x9,x10) # All of the covariates

############### stage 1 data simulation ##############
# simulate A1, true stage 1 treatment with K1=3
pi10 <- rep(1, N)
pi11 <- exp(0.5*x4 + 0.5*x1 + 0.05*x3)
pi12 <- exp(0.5*x5 - 0.5*x1 + 0.5*x2)
  
# weights matrix
matrix.pi1 <- cbind(pi10, pi11, pi12)
A1 <- A.sim(matrix.pi1)
class.A1 <- sort(unique(A1))
# propensity stage 1
pis1.hat <- M.propen(A1, cbind(x1,x2, x3,x4,x5))
  
# g1.opt <- (x2 <= 0.5)*(x1 > 0.5) + 2*(x2 > 0.5)*(x1 > -1)
# change !
g1.opt <- (x1 <= 0.5)*(x2 > -0.2) + (x1 > 0.5)*2*(1 - (x3 > -1)*(x3 < -0.5))
  
# outcome 1
Y11 <- 0.57 + exp(1.67 + 0.2*x6 - abs(1.5*x7 + x4 - 1)*((A1 - g1.opt)^2)) - 
    3*(A1 == 1) + rnorm(N,0,1) # noise on A = 1 and 2
# outcome 2
Y12 <- 1.9 + x5 + 0.5*x6 + 2*(A1 == 0)*(2*(g1.opt == 0) - 1) + 1.5*(A1 == 2)*(2*(g1.opt == 2) -1) + 0.5*(A1 == 1)*(2*(g1.opt == 1) -1) +
    1.8*(exp((A1 == 1)) -1) + rnorm(N,0,1) # noise on A = 1 and 2
# outcome 3
Y13 <- 5.32 + x8 - exp(0.1 + 2*(x10 == "No"))*(1*(A1 == 1) + 0.5*(A1 == 0) + 0.1*(A1 == 2)) + rnorm(N,0,1) # noise on A = 1 and 2
  
# stage 1 outcome
Ys1 = cbind(Y11, Y12, Y13)


# MODTR tree and average potential outcome for each of the three objectives (with tolerant rate 100%, 60%)
wt = c(w1, w2, w3)
MOTRL0 = MO.tol.DTRtree(Ys1, w = wt, A1, H=X0, delta = 0, pis.hat=pis1.hat, lambda.pct=0.02, minsplit=20,depth = 4)
#> Warning in if (m.method == "AIPW") {: the condition has length > 1 and only the
#> first element will be used
MOTRL0$tree
#>   node  X     cutoff      mEy opt.trt  avg.mEy tol.trt avg.mE(Y11) avg.mE(Y12)
#> 1    1  2 -0.2473256 3.182414      NA       NA      NA          NA          NA
#> 2    2  1  0.4898486 4.339314      NA       NA      NA          NA          NA
#> 3    3  1  0.4738403 3.412394      NA       NA      NA          NA          NA
#> 4    4 NA         NA 4.187836       0 4.187836       0    5.761938    3.607119
#> 5    5 NA         NA 4.658819       2 4.658819       2    6.380001    2.660269
#> 6    6 NA         NA 2.939191       1 2.939191       1    2.692072    5.324977
#> 7    7 NA         NA 4.288142       2 4.288142       2    4.779649    3.177743
#>   avg.mE(Y13)
#> 1          NA
#> 2          NA
#> 3          NA
#> 4   3.0543405
#> 5   4.9355880
#> 6   0.4460766
#> 7   4.8472891
MOTRL1 <- MO.tol.DTRtree(Ys1, w = wt, A1, H=X0, delta = 0.4,pis.hat=pis1.hat, lambda.pct=0.02, minsplit=20,depth = 4)
#> Warning in if (m.method == "AIPW") {: the condition has length > 1 and only the
#> first element will be used
MOTRL1$tree
#>   node  X     cutoff      mEy opt.trt  avg.mEy tol.trt avg.mE(Y11) avg.mE(Y12)
#> 1    1  2 -0.2473256 3.182414      NA       NA      NA          NA          NA
#> 2    2  1  0.4898486 4.339314      NA       NA      NA          NA          NA
#> 3    3  1  0.4738403 3.412394      NA       NA      NA          NA          NA
#> 4    4 NA         NA 4.187836       0 4.187836       0    5.761938    3.607119
#> 5    5 NA         NA 4.658819       2 4.658819       2    6.380001    2.660269
#> 6    6 NA         NA 2.939191       1 2.698213    1, 2    2.560962    2.839759
#> 7    7 NA         NA 4.288142       2 4.288142       2    4.779649    3.177743
#>   avg.mE(Y13)
#> 1          NA
#> 2          NA
#> 3          NA
#> 4    3.054341
#> 5    4.935588
#> 6    2.693742
#> 7    4.847289
```
