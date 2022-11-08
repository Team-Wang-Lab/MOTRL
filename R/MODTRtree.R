#' Multi-Objective Tree-based Reinforcement Learning for estimating optimal DTR.
#'
#' a tree-based reinforcement learning (T-RL) method to directly
#' estimate multi-objective optimal DTRs in a multi-stage multi-treatment setting.
#' The tree can be gorwn as a tolerant tree, where the tolerant rate can be defined by user.
#' At each stage, T-RL builds an unsupervised decision tree that directly handles
#' the problem of optimization with multiple treatment comparisons, through a
#' purity measure constructed with augmented inverse probability weighted estimators.
#'
#' @param Ys A matrix of outcomes of interest. The number of column is the number of objectives. Each columns contains the value of outcome of that objective for all subjects.
#' @param w A vector of weights for each objectives.
#' @param A A vector of observed treatment options.
#' @param H A matrix of covariates before assigning final treatment, excluding previous treatment variables.
#' @param delta A scalar indicates the minimum loss in pseudo-outcome that can be tolerated. The default value is 0. (non-tolerant)
#' @param tolerate A Boolean indicates whether the DTR tree is a tolerant tree or not.
#' @param pis.hat Estimated propensity score matrix.
#' @param m.method Method for calculating estimated conditional mean.
#' @param mus.reg Regression-based conditional mean outcome.
#' @param depth Maximum tree depth.
#' @param lambda.pct Minimal percent change in purity measure for split.
#' @param minsplit Minimal node size.
#' @param lookahead Whether or not to look into a further step of splitting to find the best split.
#'
#' @return Multi-Objective DTR Tree
#'
#' @name MODTRtree
#'
#' @export

MODTRtree<-function(Ys,w,A,H, delta = 0, tolerate = TRUE,
                    pis.hat=NULL,m.method=c("AIPW","randomForest"),
                    mus.reg=NULL,depth=5,lambda.pct=0.05,minsplit=20,lookahead=F) {

  # Checking the input, if not valid, return error
  if (length(w) == 1 || !is.null(dim(w))) {
    stop("weight (w) should be a vector with more than 1 element!")
  }
  if (is.null(Ys) || dim(Ys)[2] == 1) {
    stop("Outcome matrix (Ys) should be a matirx with more than 1 columns!")
  }
  if (dim(Ys)[2] != length(w)) {
    stop("Dimension of outcome matrix (Ys) and weight (w) does not match!")
  }

  # Calculate the overall reward
  Obj.dim = length(w)
  Overall.Y = tcrossprod(Ys, t(w))

  # Two cases for DTR tree
  # The tolerant tree case:
  if (tolerate) {
    tol.rate = 1-delta
    tree = tol.DTRtree(Overall.Y,A,H,tol.rate,
                       pis.hat=NULL,m.method=c("AIPW","randomForest"),
                       mus.reg=NULL,depth=5,lambda.pct=0.05,minsplit=20,lookahead=F)
    # outcomes = tree[,6] %o% w

  # The non-tolerant tree case:
  } else {
    tree = DTRtree(Overall.Y,A,H,
                   pis.hat=NULL,m.method=c("AIPW","randomForest"),
                   mus.reg=NULL,depth=5,lambda.pct=0.05,minsplit=20,lookahead=F)
    # outcomes = tree[,4] %o% w
  }

  # Get outcomes matrix
  outcomes = matrix(NA, nrow(tree), Obj.dim)
  for (i in 1:Obj.dim) {
    for (j in 1:nrow(tree)) {
      if (!is.na(tree$tol.trt[j]) || length(tree$tol.trt[j]) == 1) {
        opt.trt = tree$tol.trt[j]
        mE.Obj.j =
      } else {

      }
    }



    # Tolerant outcome at jth node
    node.j.tol <- Tol.A(A = A[I.node==j], mus.hat = mus.hat[I.node==j,], tol.rate) # c("Ey.opt1", "trt.opt1", "trt.tol1", "Ey.tol.avg1")
    Eyj.tol <- node.j.tol$Ey.tol.avg1                 # the AVERAGE expected counterfactual outcome for all tolerant treatments at jth node
    output[j, 6] = Eyj.tol                            # return at 6th column
    output$tol.trt[j] = list(node.j.tol$trt.tol1)     # and the output of 7th column is all of the tolerant treatments


  }


  colnames(outcomes) = apply(matrix(1:ncol(outcomes)), 1, function(x) paste("mE(Objective.", as.character(x), ")", sep = ""))


  # return tree and outcomes at each non-ending node
  return(cbind(tree, outcomes))
}








