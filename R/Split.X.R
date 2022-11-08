#' Split a parent node into two child nodes by a covariate.
#' Split the parent node in order to maximize the outcome of interest by signing different treatment in child nodes. Calculate new outcome means for each child node.
#' @param X A vector of a covariate for each patients.
#' @param A A vector of observed treatment indicators.
#' @param mus.hat Estimated conditional mean outcome.
#' @param minsplit Minimal node size.
#'
#' @importFrom utils combn
#' @importFrom stats quantile
#'
Split.X<-function(X,A,mus.hat,minsplit=20) { # X is a covariate for patients; minsplit is the minimum cases in each node
  n <- length(X)               # the number of patients
  X.val <- unique(X)           # all possible value of X
  n.X.val <- length(X.val)     # the number of unique possible value of X
  class.A <- sort(unique(A))   # the treatment indicators

  # Check the stopping rules
  if (n < 2*minsplit || n.X.val<2L || length(class.A) < 2L) return(NULL)

  # Case 1: X is numerical, ordered categorical
  if (is.numeric(X) == TRUE || is.ordered(X)==TRUE) { # is.ordered = ordered factorial feature
    X.val <- sort(X.val)
    # reduce computation by using quantiles
    if (n.X.val > 100L) {
      X.val <- quantile(X, 1:100/100) # change the unique x to quantiles of x (only test 100 possible x as candidates)
      n.X.val <- 100L
    }
    # initialize E(Y|optX)
    Ey.opt1.X <- trt.left.X <- trt.right.X <- rep(NA,n.X.val-1)

    for (i in 1L:(n.X.val-1)) {
      left <- which(X <= X.val[i]) # left <- index of X's that is less than the evaluated optX candidate
      if (length(left) >= minsplit && length(left) <= n - minsplit) { # Make sure after split the resulting two nodes has cases more than minsplit
        left.est <- Opt.A(A[left],mus.hat[left,which(class.A %in% unique(A[left]))])
        right.est <- Opt.A(A[-left],mus.hat[-left,which(class.A %in% unique(A[-left]))])
        trt.left.X[i] <- left.est$trt.opt1
        trt.right.X[i] <- right.est$trt.opt1
        Ey.opt1.X[i] <- length(left)/n*left.est$Ey.opt1+(1-length(left)/n)*right.est$Ey.opt1 # Average the optimum Ey for each cutoff point
      }
    }
    # pick the best split of X
    if (sum(!is.na(Ey.opt1.X)) > 0L) { # check if one of the EY for candidate cutoffs is non NA
      mEy.opt1 <- max(Ey.opt1.X, na.rm = T) # max among Ey.opt1.X
      cutoff1 <- which(Ey.opt1.X == mEy.opt1)[1] # take the minimum cutoff location to find
      X.cutoff <- X.val[cutoff1]                 # the cutoff value of X
      trt.L <- trt.left.X[cutoff1]  # find the opt treatment for the left group
      trt.R <- trt.right.X[cutoff1] # find the opt treatment for the right group

      output <- data.frame(X.cutoff, mEy.opt1, trt.L, trt.R)
      names(output) <- c("X.subset","mEy.opt1","trt.L","trt.R")
    } else {
      return(NULL)
    }
  }

  # # Case 2: X is categorical with 2 classes
  # if (is.numeric(X) == F && is.ordered(X) == F && n.X.val == 2L) { # Modification (for 2 class categorical)
  #
  #
  #   levels(X) = 1:n.X.val
  #   X = as.numeric(X)
  #
  #   X.combo <- combn(X.val,1)
  #
  #   # pick the best split of X
  #   if (sum(!is.na(Ey.opt1.X)) > 0L) {
  #     mEy.opt1 <- max(Ey.opt1.X, na.rm=T)
  #     cutoff1 <- which(Ey.opt1.X==mEy.opt1)[1]
  #     X.subset <- X.combo[,cutoff1]
  #     # change a vector into a single string while removing NA's
  #     X.subset<-paste(X.subset[!is.na(X.subset)], collapse=",")
  #     trt.L<-trt.left.X[cutoff1]
  #     trt.R<-trt.right.X[cutoff1]
  #
  #     output<-data.frame(X.subset, mEy.opt1, trt.L, trt.R) # data.frame only has one row
  #     names(output)<-c("X.subset","mEy.opt1","trt.L","trt.R")
  #
  # } else {
  #   return(NULL)
  # }

  # Case 3: X is NOT numerical or ordered categorical or categorical with 2 classes
  # which means that X is categorical with more than 2 classes
  if (is.numeric(X) == F && is.ordered(X) == F && n.X.val >= 2L) {
    # n.X.combo <- 2^(n.X.val-1)-1 # Assume there are c classes for X, each class we can either pick or not pick. -1 for not pick any class. -1 for power for the ????

    # NOT sure what n.X.combo means
    X.combo <- as.matrix(combn(X.val,1), nrow = 1) # ???????unique(X.val)?
    # if (n.X.val > 3L && n.X.val%%2 == 1L) { # ODD CASES
    #   for (k in 2L:((n.X.val-1)/2)) {
    #     X.combo <- combine.mat(X.combo,combn(X.val,k))
    #   }
    # }
    # if (n.X.val > 3L && n.X.val%%2 == 0L) { # EVEN CASES
    #   for (k in 2L:(n.X.val/2)) {
    #     if (k < (n.X.val/2)) {
    #       X.combo<-combine.mat(X.combo,combn(X.val,k))
    #     }
    #     if (k==(n.X.val/2)) {
    #       temp.mat<-combn(X.val[-1],k-1)
    #       first.row<-rep(X.val[1],ncol(temp.mat))
    #       X.combo<-combine.mat(X.combo,rbind(first.row,temp.mat))
    #     }
    #   }
    # }

    # Modification
    if (n.X.val > 2) {
      for (k in 2L:(n.X.val-1)) {
        X.combo <- combine.mat(X.combo,combn(X.val,k))
      }
    }

    # Ey.opt1.X <- trt.left.X <- trt.right.X <- rep(NA, n.X.combo)
    Ey.opt1.X <- trt.left.X <- trt.right.X <- rep(NA, ncol(X.combo))

    # for (i in 1L:n.X.combo) { # original
    for (i in 1L:ncol(X.combo)) {
      left <- which(X %in% X.combo[,i])
      if (length(left) >= minsplit && length(left) <= n-minsplit) { # this combo can be split as one left child node

        left.est <- Opt.A(A[left],mus.hat[left,which(class.A %in% unique(A[left]))])
        right.est <- Opt.A(A[-left],mus.hat[-left,which(class.A %in% unique(A[-left]))])

        trt.left.X[i]<-left.est$trt.opt1
        trt.right.X[i]<-right.est$trt.opt1
        Ey.opt1.X[i]<-length(left)/n*left.est$Ey.opt1+(1-length(left)/n)*right.est$Ey.opt1
      }
    }

    # pick the best split of X
    if (sum(!is.na(Ey.opt1.X)) > 0L) {
      mEy.opt1 <- max(Ey.opt1.X, na.rm=T)
      cutoff1 <- which(Ey.opt1.X==mEy.opt1)[1]
      X.subset <- X.combo[,cutoff1]
      # change a vector into a single string while removing NA's
      # X.subset <- paste(X.subset[!is.na(X.subset)], collapse=",")

      trt.L<-trt.left.X[cutoff1]
      trt.R<-trt.right.X[cutoff1]

      output <- as.data.frame(matrix(NA,1,4)) # data.frame only has one row
      names(output) <- c("X.subset","mEy.opt1","trt.L","trt.R")
      output$X.subset[1] = list(X.subset[!is.na(X.subset)])
      output[1,2] = mEy.opt1
      output[1,3:4] = c(trt.L, trt.R)
    } else {
      return(NULL)
    }
  }
  return(output)
  #RETURN all the avg Y*|cutoff=each X
}
