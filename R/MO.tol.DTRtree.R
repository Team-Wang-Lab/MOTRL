#' Multi-Objective Tree-based Reinforcement Learning for estimating tolerant DTR.
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
#' @param m.method Method for calculating estimated conditional mean. Avaliabel options are "AIPW" and "randomForest".
#' @param mus.reg Regression-based conditional mean outcome.
#' @param depth Maximum tree depth.
#' @param lambda.pct Minimal percent change in purity measure for split.
#' @param minsplit Minimal node size.
#' @param lookahead Whether or not to look into a further step of splitting to find the best split.
#' @param PO.loss ...
#'
#' @return Multi-Objective DTR Tree
#' @export

MO.tol.DTRtree<-function(Ys,w,A,H, delta = 0, tolerate = TRUE,
                         pis.hat=NULL,m.method=c("AIPW","randomForest"),
                         mus.reg=NULL,depth=5,lambda.pct=0.05,minsplit=20,
                         lookahead=F, PO.loss = NULL) {

  `%notin%` <- Negate(`%in%`)
  # Checking the input, if not valid, return error
  if (!is.null(dim(w))) {
    stop("weight (w) should be a vector with 1 or more than 1 element!")
  }
  if (is.null(Ys)) {
    stop("Outcome matrix (Ys) should be a matirx with 1 or more than 1 columns!")
  }
  if (dim(Ys)[2] != length(w)) {
    stop("Dimension of outcome matrix (Ys) and weight (w) does not match!")
  }

  # Calculate the overall reward
  Obj.dim = length(w)
  if (Obj.dim == 1) {
    Overall.Y = c(Ys)[[1]]
  } else {
    Overall.Y = tcrossprod(as.matrix(Ys), t(w))
  }

  if (!is.null(PO.loss)) {
    Ys = Ys + PO.loss
  }

  # - Two cases for DTR tree -
  # 1. The tolerant tree case:
  if (tolerate) {
    tol.rate = 1 - delta

    # Initialization
    n <- nrow(Ys) # number of people
    I.node <- rep(1,n) # indicator of nodes
    class.A <- sort(unique(A))
    tree <- as.data.frame(matrix(NA,2^(depth+1)-1,7+Obj.dim))

    if (is.null(colnames(Ys))) {
      outcome.names = apply(matrix(1:Obj.dim), 1, function(x) paste("avg.mE(Objective.", as.character(x), ")", sep = ""))
    } else {
      outcome.names = apply(matrix(colnames(Ys)), 1, function(x) paste("avg.mE(", x, ")", sep = ""))
    }
    names(tree)<-c("node","X","cutoff","mEy","opt.trt", "avg.mEy", "tol.trt", outcome.names)
    tree[,1] <- 1L:nrow(tree)  # this does not equal to the number of rows(2^k+1)

    # Estimate mus.hat if not given
    list.mus.hat <- vector(mode='list', length = Obj.dim)
    if (m.method[1] == "AIPW") {
      # estimate propensity matrix if not given, using all data
      # same propensity for all subset data
      if (is.null(pis.hat)) {
        pis.hat <- M.propen(A = A,Xs = H)
      }
      # estimate mus.hat AIPW for each objectives
      for (i in 1:Obj.dim) {
        mus.reg.i <- Reg.mu(Y = Ys[,i],As = A, H = H)$mus.reg
        list.mus.hat[[i]] <- mus.AIPW(Y = Ys[,i], A = A, pis.hat = pis.hat, mus.reg = mus.reg.i)
      }
      # the compound AIPW estimates
      mus.hat = matrix(0, nrow = n, ncol = length(class.A))
      for (i in 1:Obj.dim) {
        mus.hat = mus.hat + as.matrix(list.mus.hat[[i]]*w[i])
      }

      # Random Forest
    } else if (m.method[1] == "randomForest") {
      require(randomForest)
      # estimate mus.hat with RF for each objectives
      for (i in 1:Obj.dim) {
        RF.i <- randomForest(Ys[,i] ~., data=data.frame(A,H))
        mus.hat <- matrix(0, nrow = n, ncol = length(class.A))
        for (a in 1L:length(class.A)) {
          mus.hat[,a] <- predict(RF.i, newdata=data.frame(A=rep(class.A[a],n),H))
        }
        list.mus.hat[[i]] <- mus.hat
      }
      # the compound RF estimates
      mus.hat = matrix(0, nrow = n, ncol = length(class.A))
      for (i in 1:Obj.dim) {
        mus.hat = mus.hat + as.matrix(list.mus.hat[[i]]*w[i])
      }
    } else {
      stop("The method for estimating conditional means is not available!")
    }

    # estimate mus.hat for each objectives
    list.mus.hat <- vector(mode='list', length = Obj.dim)
    for (i in 1:Obj.dim) {
      pis.hat.i <- M.propen(A = A,Xs = H)
      mus.reg.i <- Reg.mu(Y = Ys[,i],As = A, H = H)$mus.reg
      list.mus.hat[[i]] <- mus.AIPW(Y = Ys[,i], A = A, pis.hat = pis.hat.i, mus.reg = mus.reg.i)
    }

    # Expected outcome at root
    root <- Opt.A(A, mus.hat)
    root.opt.trt = root$trt.opt1
    Ey0 <- root$Ey.opt1 # counterfactual outcome for A.opt
    # Tolerant outcome at root
    root.tol <- Tol.A(A, mus.hat, tol.rate) # c("Ey.opt1", "trt.opt1", "trt.tol1", "Ey.tol.avg1")
    root.tol.trt = list(root.tol$trt.tol1)
    Ey0.tol <- root.tol$Ey.tol.avg1 # the AVERAGE counterfactual outcome for all tolerant treatments at root node
    # split if improved at least lambda, as a percent of Ey0
    lambda <- abs(Ey0) * lambda.pct

    for (k in 1L:(depth+1)) { # depth is the most number of split to reach one terminal node

      # The 1st node
      if (k==1L) {
        # apply look-ahead to the first split, the most important split
        # only to first split so as to save computation time
        # use a larger minsplit for the first split
        if (lookahead) {
          best.H.1<-best.H.lh(H=H,A=A,mus.hat=mus.hat,minsplit=0.15*n)
        } else {
          best.H.1<-best.H(H=H,A=A,mus.hat=mus.hat,minsplit=minsplit) # choose the best covariate for splitting
        }
        if (is.null(best.H.1) == FALSE && best.H.1$mEy.opt1 > Ey0 + lambda) { # If meet the split criteria

          # if numerical split X
          if (is.numeric(best.H.1$X.subset[[1]])) {
            tree[k,2] <- best.H.1$X
            tree[k,3] <- best.H.1$X.subset
            tree[k,4] <- best.H.1$mEy.opt1
            tree[k,5] <- NA
            I.node[I.node == k & H[,best.H.1$X] <= best.H.1$X.subset] <- 2*k
            tree[2*k,2:5] <- c(NA,NA,NA,best.H.1$trt.L)
            I.node[I.node == k & H[,best.H.1$X] > best.H.1$X.subset] <- 2*k+1
            tree[2*k+1,2:5] <- c(NA,NA,NA,best.H.1$trt.R)

            # else if non-numerical X (a list)
          } else {
            tree[k,2] = best.H.1$X
            tree$cutoff[k] = list(best.H.1$X.subset)
            tree[k,4] = best.H.1$mEy.opt1
            I.node[I.node == k & H[,best.H.1$X] %in% unlist(best.H.1$X.subset)] <- 2*k
            tree[2*k,2:5] <- c(NA,NA,NA,best.H.1$trt.L)
            I.node[I.node == k & H[,best.H.1$X] %notin% unlist(best.H.1$X.subset)] <- 2*k+1 ### notin! **********
            tree[2*k+1,2:5] <- c(NA,NA,NA,best.H.1$trt.R)
          }

        } else {                                # If NOT meet the split criteria, then terminate at this node
          tree[k,4] <-Ey0                       # previous method: if end at the first node, directly return the Ey0 and opt.trt for everybody
          tree[k,5] <-root.opt.trt
          tree[k, 6] = Ey0.tol
          tree$tol.trt[k] = list(root.tol.trt)  # and output for tol.trt is all of the tolerant treatments for everyone in node1
          tree[k, 8:(7+Obj.dim)] = avg.Obj(class.A, root.tol.trt, list.mus.hat)
          break                                 # Finish the splitting at this first node
        }

        # The j^{th} split
      } else {
        for (j in (2^(k-1)):(2^k-1)) {
          if (!is.na(tree[trunc(j/2),2])) {     # if the parent node is not terminated
            best.H.j <- best.H(H = H[I.node==j,], A = A[I.node==j], mus.hat = mus.hat[I.node==j,], minsplit = minsplit)
            if (is.null(best.H.j) == FALSE && as.numeric(best.H.j$mEy.opt1) > as.numeric(tree[trunc(j/2),4]) + lambda
                && j<2^depth && length(unique(A[I.node == j])) == length(class.A)) { # meet the split criteria

              # if non-numerical X
              if (is.numeric(best.H.j$X.subset[[1]])) {
                tree[j,2] <- best.H.j$X
                tree[j,3] <- best.H.j$X.subset
                tree[j,4] <- best.H.j$mEy.opt1
                tree[j,5] <- NA
                I.node[I.node == j & H[,best.H.j$X] <= best.H.j$X.subset] <- 2*j
                tree[2*j,2:5]<-c(NA,NA,NA,best.H.j$trt.L)
                I.node[I.node == j & H[,best.H.j$X] > best.H.j$X.subset] <- 2*j+1
                tree[2*j+1,2:5]<-c(NA,NA,NA,best.H.j$trt.R)
              }
              # else if numerical split X
              else {
                tree[j,2] = best.H.j$X
                tree$cutoff[j] = list(best.H.j$X.subset)
                tree[j,4] = best.H.j$mEy.opt1
                tree[j,5] = NA
                I.node[I.node == j & H[,best.H.j$X] %in% unlist(best.H.j$X.subset)] <- 2*j
                tree[2*j,2:5] <- c(NA,NA,NA,best.H.j$trt.L)
                I.node[I.node == j & H[,best.H.j$X] %notin% unlist(best.H.j$X.subset)] <- 2*j+1
                tree[2*j+1,2:5] <- c(NA,NA,NA,best.H.j$trt.R)
              }

            } else { # If not meet spliting criterion, we end at this jth node
              node.j = Opt.A(A = A[I.node==j], mus.hat = mus.hat[I.node==j,])
              tree[j,4]<-node.j$Ey.opt1          # return the mEy for everybody at node j
              tree[j,5]<-node.j$trt.opt1         # return the opt.trt for everybody at node j

              # Tolerant outcome at jth node
              node.j.tol <- Tol.A(A = A[I.node==j], mus.hat = mus.hat[I.node==j,], tol.rate)
              Eyj.tol <- node.j.tol$Ey.tol.avg1               # the AVERAGE expected counterfactual outcome for all tolerant treatments at jth node
              tree[j, 6] = Eyj.tol                            # return at 6th column
              tree$tol.trt[j] = list(node.j.tol$trt.tol1)     # and the output of 7th column is all of the tolerant treatments

              # Return mE(Obj) values
              list.mus.hat.j = lapply(list.mus.hat, function(x) x[I.node == j,])
              tree[j, 8:(7+Obj.dim)] = avg.Obj(class.A, node.j.tol$trt.tol1, list.mus.hat.j)
            }
          }
        }
      }
    }
    tree = tree[!is.na(tree[,2]) | !is.na(tree[,5]),]
    tree$mEy = as.numeric(tree$mEy)

    #  Return a matrix for PO.loss
    tol.regimes = predict_tol.DTR(tree, newdata = H)
    avg.PO = matrix(0, nrow = n, ncol = Obj.dim)

    if (nrow(tree) == 1) {
      for (i in 1L:n) {
        for (d in 1L:Obj.dim) {
          avg.PO[i,d] = mean(list.mus.hat[[d]][i,which(class.A %in% tol.regimes[[i]][[1]])])
        }
      }
    } else {
      for (i in 1L:n) {
        for (d in 1L:Obj.dim) {
          avg.PO[i,d] = mean(list.mus.hat[[d]][i,which(class.A %in% tol.regimes[[i]])])
        }
      }
    }
    PO.loss = as.matrix(avg.PO - Ys)


    # 2. The non-tolerant tree case:
  } else {
    tree = DTRtree(Overall.Y,A,H,
                   pis.hat,m.method,
                   mus.reg,depth,lambda.pct,minsplit,lookahead)
  }
  return(list(tree = tree,
              PO.loss = PO.loss,
              POs = avg.PO))
}












