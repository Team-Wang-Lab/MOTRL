#' Tree-based Reinforcement Learning for estimating optimal tolerant DTR.
#'
#' a tree-based reinforcement learning (T-RL) method to directly
#' estimate optimal DTRs in a multi-stage multi-treatment setting. At
#' each stage, T-RL builds an unsupervised decision tree that directly handles
#' the problem of optimization with multiple treatment comparisons, through a
#' purity measure constructed with augmented inverse probability weighted estimators.
#'
#' @param Y A vector of outcome of interest.
#' @param A A vector of observed treatment options.
#' @param H A matrix of covariates before assigning final treatment, excluding previous treatment variables.
#' @param tol.rate A scalar indicates the tolerant rate
#' @param pis.hat Estimated propensity score matrix.
#' @param m.method Method for calculating estimated conditional mean.
#' @param mus.reg Regression-based conditional mean outcome.
#' @param depth Maximum tree depth.
#' @param lambda.pct Minimal percent change in purity measure for split.
#' @param minsplit Minimal node size.
#' @param lookahead Whether or not to look into a further step of splitting to find the best split.
#'
#' @importFrom randomForest randomForest
#' @importFrom stats predict
#'
#' @export
tol.DTRtree <- function(Y,A,H,tol.rate,
                      pis.hat=NULL,m.method=c("AIPW","randomForest"),
                      mus.reg=NULL,depth=5,lambda.pct=0.05,minsplit=20,lookahead=F) {
  # initialization
  # indicator for subset data
  n <- length(Y) # number of people
  I.node <- rep(1,n) # indicator of nodes
  class.A <- sort(unique(A))

  # output <- matrix(NA,1,5)
  # *Modification: matrix -> data.frame
  # *Modification: add two column for avg.mEy and tol.trt
  output <- as.data.frame(matrix(NA,1,7))

  # colnames(output)<-c("node","X","cutoff","mEy","opt.trt")
  # *Modification: change colnames to names
  # *Modification: add two column for avg.mEy and tol.trt
  names(output)<-c("node","X","cutoff","mEy","opt.trt", "avg.mEy", "tol.trt")

  # estimate mus.hat if not given
  if(m.method[1]=="AIPW"){
    # estimate propensity matrix if not given, using all data
    # same propensity for all subset data
    if(is.null(pis.hat)) pis.hat<-M.propen(A=A,Xs=H)
    if(is.null(mus.reg)) mus.reg<-Reg.mu(Y=Y,As=A,H=H)$mus.reg
    mus.hat<-mus.AIPW(Y=Y,A=A,pis.hat=pis.hat,mus.reg=mus.reg)
  } else if(m.method[1]=="randomForest"){
    # require(randomForest)
    RF<-randomForest(Y~., data=data.frame(A,H))
    mus.hat<-matrix(NA,n,length(class.A))
    for(i in 1L:length(class.A)) mus.hat[,i]<-predict(RF,newdata=data.frame(A=rep(class.A[i],n),H))
  } else{
    stop("The method for estimating conditional means is not available!")
  }

  # Expected outcome at root
  root <- Opt.A(A,mus.hat)
  root.opt.trt = root$trt.opt1
  Ey0 <- root$Ey.opt1 # counterfactual outcome for A.opt
  # Tolerant outcome at root
  root.tol <- Tol.A(A, mus.hat, tol.rate) # c("Ey.opt1", "trt.opt1", "trt.tol1", "Ey.tol.avg1")
  Ey0.tol <- root.tol$Ey.tol.avg1 # the AVERAGE counterfactual outcome for all tolerant treatments at root node

  # split if improved at least lambda, as a percent of Ey0
  lambda<-abs(Ey0)*lambda.pct

  for (k in 1L:depth) { # depth is the most number of split to reach one terminal node

    # output <- rbind(output,matrix(NA,2^k,5) ) # 2^k??????? originally output=1*5
    # *Modification: since we have data.frame, we add 2^k rows of NA to output
    output[nrow(output) + 1:2^k, ] <- NA

    output[,1] <- 1L:(2^(k+1)-1) # this does not equal to the number of rows(2^k+1)

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
        output[k,2:5] <- c(best.H.1$X, best.H.1$X.subset, best.H.1$mEy.opt1, NA)
        I.node[I.node == k & H[,best.H.1$X] <= best.H.1$X.subset] <- 2*k
        output[2*k,2:5] <- c(NA,NA,NA,best.H.1$trt.L)
        I.node[I.node == k & H[,best.H.1$X] > best.H.1$X.subset] <- 2*k+1
        output[2*k+1,2:5] <- c(NA,NA,NA,best.H.1$trt.R)
      } else { # If NOT meet the split criteria (terminate node)
        output[k,4:5]<-c(Ey0,root.opt.trt)           # previous method: if end at the first node, directly return the Ey0 and opt.trt for everybody
        output[k, 6] = Ey0.tol                       # *Modification: We changed the expected outcome at root Ey0 to
                                                     # the average of expected counterfactual outcome for all tol.A;
        output$tol.trt[k] = list(root.tol$trt.tol1)  # and output for tol.trt is all of the tolerant treatments for everyone in node1
        break # Finish the splitting at this first node
      }

    # The jth split
    } else {
      for (j in (2^(k-1)):(2^k-1)) {
        if (!is.na(output[trunc(j/2),2])) { # if the parent node is not terminated
          best.H.j <- best.H(H = H[I.node==j,], A = A[I.node==j], mus.hat = mus.hat[I.node==j,], minsplit = minsplit)
          if (is.null(best.H.j) == FALSE && best.H.j$mEy.opt1 > output[trunc(j/2),4] + lambda) { # meet the split criteria
            output[j,2:5]<-c(best.H.j$X, best.H.j$X.subset, best.H.j$mEy.opt1, NA)
            I.node[I.node == j & H[,best.H.j$X] <= best.H.j$X.subset] <- 2*j
            output[2*j,2:5]<-c(NA,NA,NA,best.H.j$trt.L)
            I.node[I.node == j & H[,best.H.j$X] > best.H.j$X.subset] <- 2*j+1
            output[2*j+1,2:5]<-c(NA,NA,NA,best.H.j$trt.R)
          } else { # *Modification: if not meet spliting criterion, we end at this jth node,
                   # the previous method do nothing here, since the jth row already have the output as c(j, NA, NA, NA, opt.trt)
                   # however, to get the tolerant regimes, we need to modify output here (to add the mEy, opt.trt(should be same), avg.mEy, and tol.trt)

            # First, change the mEy from NA to the expected optimal outcome for all subject in node j
            node.j = Opt.A(A = A[I.node==j], mus.hat = mus.hat[I.node==j,])
            output[j,4:5]<-c(node.j$Ey.opt1, node.j$trt.opt1) # return the mEy and opt.trt for everybody at node j
            # Tolerant outcome at jth node
            node.j.tol <- Tol.A(A = A[I.node==j], mus.hat = mus.hat[I.node==j,], tol.rate) # c("Ey.opt1", "trt.opt1", "trt.tol1", "Ey.tol.avg1")
            Eyj.tol <- node.j.tol$Ey.tol.avg1                 # the AVERAGE expected counterfactual outcome for all tolerant treatments at jth node
            output[j, 6] = Eyj.tol                            # return at 6th column
            output$tol.trt[j] = list(node.j.tol$trt.tol1)     # and the output of 7th column is all of the tolerant treatments
          }
        }
      }
      if (sum(is.na(output[(2^(k-1)):(2^k-1),2])) == 2^(k-1)) break
    }
  }
  output <- output[!is.na(output[,2]) | !is.na(output[,5]),]
  return(output)
}



