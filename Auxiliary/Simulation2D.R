## ---- include = FALSE---------------------------------------------------------
knitr::opts_chunk$set(
  collapse = TRUE,
  comment = "#>"
)

## ----setup, include=FALSE-----------------------------------------------------
knitr::opts_chunk$set(echo = TRUE)
library(MOTRL)
library(rpart)
library(randomForest)
library(wakefield)
library(dplyr)
library(psych)
library(psychTools)
library(readr)
library(scales)
library(ggplot2)

## -----------------------------------------------------------------------------
# 1. function to summarize simulation results
summary2 <- function(x) { # x is a numerical vector or matrix
  s1 <- summary(x)
  sd2 <- function(y){
    return(sd(y,na.rm = T))
  }
  if(is.matrix(x)) {
    SD<-apply(x,2,sd2)
  } else {
    SD <-sd2(x)
  }
  s2<-list(s1,SD)
  names(s2)<-c("summary","SD") # return summary and sd to each column
  return(s2)
}

print.summary <- function(x) { # x is a numerical vector or matrix
  avg.x = round(mean(x), 2)
  sd2 <- function(y){
    return(sd(y,na.rm = T))
  }
  if(is.matrix(x)) {
    SD <-round(apply(x,2,sd2), 2)
  } else {
    SD <- round(sd2(x), 2)
  }
  out = paste(avg.x, " (", SD, ")", sep = "")
  return(out)
}

# 2. function to sample treatment A
# input matrix.pi as a matrix of sampling probabilities, which could be non-normalized
A.sim <- function(matrix.pi) {
  N <- nrow(matrix.pi) # sample size
  K <- ncol(matrix.pi) # treatment options
  if (N<=1 | K<=1) stop("Sample size or treatment options are insufficient!")
  if (min(matrix.pi)<0) stop("Treatment probabilities should not be negative!")
  
  # normalize probabilities to add up to 1 and simulate treatment A for each row
  probs <- t(apply(matrix.pi,1,function(x){x/sum(x,na.rm = TRUE)}))
  A <- apply(probs, 1, function(x) sample(0:(K-1), 1, prob = x))
  return(A)
}

## -----------------------------------------------------------------------------
N<-1000 # sample size of training data
N2<-1000 # sample size of test data
iter <- 10 # replication
w0 = 0.8

perc.TRL11 = perc.TRL12 = perc.MOTRL10 = perc.MOTRL11 = perc.MOTRL12 = perc.MOTRL13 = rep(NA,iter)
EYs.TRL11.a = EYs.TRL12.a = EYs.TRL11.b = EYs.TRL12.b = 
  EYs.MOTRL10.a = EYs.MOTRL10.b = 
  EYs.MOTRL11.a = EYs.MOTRL11.b = 
  EYs.MOTRL12.a = EYs.MOTRL12.b = 
  EYs.MOTRL13.a = EYs.MOTRL13.b = rep(NA,iter) # estimated mean counterfactual outcome

for (i in 1:iter) {
  # Simulation begin
  set.seed(i+300)
  x1<-rnorm(N)              # each covariates follows N(0,1)
  x2<-rnorm(N)
  x3<-rnorm(N)
  x4<-rnorm(N)
  x5<-rnorm(N)
  x6<-answer(N, x = c("No", "Yes"), name = "Smoke")
  X0<-cbind(x1,x2,x3,x4,x5,x6) # All of the covariates

  ############### stage 1 data simulation ##############
  # simulate A1, true stage 1 treatment with K1=3
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
  
  # 3 models #
  ################# Objective 1 ####################
  # simulate stage 1 optimal g11.opt for reward1
  Y11 <- exp(1.68 + 0.2*x4 - abs(0.5*x1 - 1)*((A1 - g1.opt)^2)) - 3*(A1 == 1) + rnorm(N,0,1)
 
  ################# Objective 2 ####################
  # simulate stage 1 optimal g12.opt for reward2
  Y12 <- 2.37 + x4 + 2*x5 + 2*(A1 == 0)*(2*(g1.opt == 0) - 1) + 1.5*(A1 == 2)*(2*(g1.opt == 2) -1) + 1.5*(exp((A1 == 1)) -1) + rnorm(N,0,1)
  
  # Y11 = rescale(Y11, to = c(0,5))
  # Y12 = rescale(Y12, to = c(0,5))

  ############ Multi-Objective weighted-sum reward at stage 1 ##########
  # stage 1 outcome
  Ys1 = cbind(Y11, Y12)
  
  ################ Grow the trees ######################
  # DTRtree on reward1
  TRLtree11 <- TRL::DTRtree(Y11, A1, H=X0, pis.hat=pis1.hat, lambda.pct=0.05, minsplit=max(0.05*N,20))
  # DTRtree on reward2
  TRLtree12 <- TRL::DTRtree(Y12, A1, H=X0, pis.hat=pis1.hat, lambda.pct=0.05, minsplit=max(0.05*N,20))
  
  # MODTRtree on overall reward (with tolerant rate 100%, 90%, 70%, 50%,
  
  w11 = c(w0, 1-w0)
  MOTRL10 = MO.tol.DTRtree(Ys1, w = w11, A1, H=X0, delta = 0, pis.hat=pis1.hat, lambda.pct=0.05, minsplit=max(0.05*N,20),depth = 3)
  MOTRLtree10 = MOTRL10$tree
  MOTRL11 <- MO.tol.DTRtree(Ys1, w = w11, A1, H=X0, delta = 0.1,pis.hat=pis1.hat, lambda.pct=0.05, minsplit=max(0.05*N,20),depth = 3)
  MOTRLtree11 = MOTRL11$tree
  MOTRL12 <- MO.tol.DTRtree(Ys1, w = w11, A1, H=X0, delta = 0.3,pis.hat=pis1.hat, lambda.pct=0.05, minsplit=max(0.05*N,20),depth = 3)
  MOTRLtree12 = MOTRL12$tree
  MOTRL13 <- MO.tol.DTRtree(Ys1, w = w11, A1, H=X0, delta = 0.5,pis.hat=pis1.hat, lambda.pct=0.05, minsplit=max(0.05*N,20),depth = 3)
  MOTRLtree13 = MOTRL13$tree
  
  ############################################
  # prediction using new data
  ############################################
  set.seed(i+10000)
  x1.test = rnorm(N2)
  x2.test = rnorm(N2)
  x3.test = rnorm(N2)
  x4.test = rnorm(N2)
  x5.test = rnorm(N2)
  x6.test = answer(N2, x = c("No", "Yes"), name = "Smoke")
  X0.test = cbind(x1.test,x2.test,x3.test,x4.test,x5.test,x6.test) # All of the covariates
  
  #new:
  g1.opt.test <- (x2.test <= 0.5)*(x1.test > 0.5) + 2*(x2.test > 0.5)*(x1.test > -1)
  
  ####### stage 1 prediction #######
  # predict selection %
  g1.TRL11 = TRL::predict_DTR(TRLtree11,newdata=data.frame(X0.test))
  g1.TRL12 = TRL::predict_DTR(TRLtree12,newdata=data.frame(X0.test))
  g1.MOTRL10 = predict_tol.DTR(MOTRLtree10, newdata=data.frame(X0.test)) 
  g1.MOTRL11 = predict_tol.DTR(MOTRLtree11, newdata=data.frame(X0.test)) 
  g1.MOTRL12 = predict_tol.DTR(MOTRLtree12, newdata=data.frame(X0.test)) 
  g1.MOTRL13 = predict_tol.DTR(MOTRLtree13, newdata=data.frame(X0.test)) 
  
  # percentage of true prediction
  # TRL
  perc.TRL11[i] = 100 * mean(g1.opt.test == g1.TRL11)
  perc.TRL12[i] = 100 * mean(g1.opt.test == g1.TRL12)
  # MOTRL
  perc.MOTRL10[i] = 100 * mean(mapply(function(x, y) x %in% y, g1.opt.test, g1.MOTRL10))
  perc.MOTRL11[i] = 100 * mean(mapply(function(x, y) x %in% y, g1.opt.test, g1.MOTRL11))
  perc.MOTRL12[i] = 100 * mean(mapply(function(x, y) x %in% y, g1.opt.test, g1.MOTRL12))
  perc.MOTRL13[i] = 100 * mean(mapply(function(x, y) x %in% y, g1.opt.test, g1.MOTRL13))

  
                                        # counterfactual mean outcome for TRL
  ##### reward 1
  # 对应Y11.1  # 4.48
  EYs.TRL11.a[i] = mean(exp(1.68 + 0.2*x4.test - abs(0.5*x1.test - 1)*((g1.TRL11 - g1.opt.test)^2)) - 3*(g1.TRL11 == 1) + rnorm(N2,0,1))
  EYs.TRL12.a[i] = mean(exp(1.68 + 0.2*x4.test - abs(0.5*x1.test - 1)*((g1.TRL12 - g1.opt.test)^2)) - 3*(g1.TRL12 == 1) + rnorm(N2,0,1))
  
  ##### reward 2  # 2.5
  EYs.TRL11.b[i] = mean(2.36 + x4.test + 2*x5.test + 2*(g1.TRL11 == 0)*(2*(g1.opt.test == 0) - 1) +
    1.5*(g1.TRL11 == 2)*(2*(g1.opt.test ==2) - 1) + 1.5*(exp((g1.TRL11 == 1)) -1) + rnorm(N2,0,1))
  EYs.TRL12.b[i] = mean(2.36 + x4.test + 2*x5.test + 2*(g1.TRL12 == 0)*(2*(g1.opt.test == 0) - 1) +
    1.5*(g1.TRL12 == 2)*(2*(g1.opt.test ==2) - 1) + 1.5*(exp((g1.TRL12 == 1)) -1) + rnorm(N2,0,1))
  
   
                                             # counterfactual mean outcome for MOTRL
  ##### reward 1
  EYs.MOTRL10.a[i] = mean(mapply(function(x, y, a, b) {exp(1.68 + 0.2*a - abs(0.5*b - 1)*mean((unlist(x) - y)^2)) - mean(3*(unlist(x) == 1)) + rnorm(N2,0,1)},
                                 g1.MOTRL10, g1.opt.test, x4.test, x1.test))
  EYs.MOTRL11.a[i] = mean(mapply(function(x, y, a, b) {exp(1.68 + 0.2*a - abs(0.5*b - 1)*mean((unlist(x) - y)^2)) - mean(3*(unlist(x) == 1)) + rnorm(N2,0,1)},
                                 g1.MOTRL11, g1.opt.test, x4.test, x1.test))
  EYs.MOTRL12.a[i] = mean(mapply(function(x, y, a, b) {exp(1.68 + 0.2*a - abs(0.5*b - 1)*mean((unlist(x) - y)^2)) - mean(3*(unlist(x) == 1)) + rnorm(N2,0,1)},
                                 g1.MOTRL12, g1.opt.test, x4.test, x1.test))
  EYs.MOTRL13.a[i] = mean(mapply(function(x, y, a, b) {exp(1.68 + 0.2*a - abs(0.5*b - 1)*mean((unlist(x) - y)^2)) - mean(3*(unlist(x) == 1)) + rnorm(N2,0,1)},
                                 g1.MOTRL13, g1.opt.test, x4.test, x1.test))


  ##### reward 2
  EYs.MOTRL10.b[i] = mean(mapply(function(x, y, a, b, c) 
    {2.37 + a + 2*b + mean(2*(unlist(x) == 0)*(2*(y == 0) - 1)) + mean(1.5*(unlist(x) == 2)*(2*(y ==2) -1)) + mean(1.5*(exp((unlist(x) == 1)) -1)) + rnorm(N2,0,1)}, 
    g1.MOTRL10, g1.opt.test, x4.test, x5.test, x6.test))
  EYs.MOTRL11.b[i] = mean(mapply(function(x, y, a, b, c) 
    {2.37 + a + 2*b + mean(2*(unlist(x) == 0)*(2*(y == 0) - 1)) + mean(1.5*(unlist(x) == 2)*(2*(y ==2) -1)) + mean(1.5*(exp((unlist(x) == 1)) -1)) + rnorm(N2,0,1)},   
    g1.MOTRL11, g1.opt.test, x4.test, x5.test, x6.test))
  EYs.MOTRL12.b[i] = mean(mapply(function(x, y, a, b, c) 
    {2.37 + a + 2*b + mean(2*(unlist(x) == 0)*(2*(y == 0) - 1)) + mean(1.5*(unlist(x) == 2)*(2*(y ==2) -1)) + mean(1.5*(exp((unlist(x) == 1)) -1)) + rnorm(N2,0,1)}, 
    g1.MOTRL12, g1.opt.test, x4.test, x5.test, x6.test))
  EYs.MOTRL13.b[i] = mean(mapply(function(x, y, a, b, c) 
    {2.37 + a + 2*b + mean(2*(unlist(x) == 0)*(2*(y == 0) - 1)) + mean(1.5*(unlist(x) == 2)*(2*(y ==2) -1)) + mean(1.5*(exp((unlist(x) == 1)) -1)) + rnorm(N2,0,1)}, 
    g1.MOTRL13, g1.opt.test, x4.test, x5.test, x6.test))
}

Result = data.frame(matrix(NA, ncol = 3, nrow = 6))
names(Result) = c("perc.opt", "mE{Y1*(g_hat)}", "mE{Y2*(g_hat)}")

Result$perc.opt = c(print.summary(perc.TRL11), print.summary(perc.TRL12), print.summary(perc.MOTRL10), print.summary(perc.MOTRL11), print.summary(perc.MOTRL12), print.summary(perc.MOTRL13))

Result$`mE{Y1*(g_hat)}` = c(print.summary(EYs.TRL11.a), print.summary(EYs.TRL12.a), print.summary(EYs.MOTRL10.a), print.summary(EYs.MOTRL11.a), print.summary(EYs.MOTRL12.a), print.summary(EYs.MOTRL13.a))

Result$`mE{Y2*(g_hat)}` = c(print.summary(EYs.TRL11.b), print.summary(EYs.TRL12.b), print.summary(EYs.MOTRL10.b), print.summary(EYs.MOTRL11.b), print.summary(EYs.MOTRL12.b), print.summary(EYs.MOTRL13.b))

Result

