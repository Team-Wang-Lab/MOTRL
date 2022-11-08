#' Choose the optimal treatment with given pseudo outcome matrix in one given stage.
#'
#' @param A Treatment options given in this stage
#' @param mus.hat The final counterfactual outcome according to the treatment.
Opt.A<-function(A, mus.hat){
  class.A<-sort(unique(A))
  if (length(class.A) == 1) {
    trt.opt1<-class.A
    Ey.opt1<-mean(mus.hat)
  } else {
    if (length(A)!= nrow(mus.hat)) {
      stop("Treatment options and mean matrix dimension do not match!")
    }
    if (length(class.A)!= ncol(mus.hat)) {
      stop("Treatment options in one node do not inclued all treatment options! Try a bigger minsplit")
      # mus.hat = mus.hat[,class.A]
      # mus.hat = mus.hat[,which(all.A %in% class.A)]
    }

    # pick a single best treatment for all patients
    c.means<-apply(mus.hat,2,mean)
    Ey.opt1<-max(c.means)
    trt.opt1<-class.A[which(c.means == Ey.opt1)]###Problem: if the order of class.A is different from mus.hat because of the sort() this will cause problems
  }
  outs<-list(Ey.opt1, trt.opt1)
  names(outs)<-c("Ey.opt1", "trt.opt1")
  # outs = as.data.frame(matrix(NA,1,2))
  # names(outs) <- c("Ey.opt1","trt.opt1")
  # output[1,1] = Ey.opt1
  # output[1,2] = trt.opt1
  return(outs)
}


# Opt.A<-function(A,mus.hat){
#   class.A<-sort(unique(A))
#   if (length(class.A) == 1) {
#     trt.opt1<-class.A
#     Ey.opt1<-mean(mus.hat)
#   } else {
#     if (length(A)!= nrow(mus.hat)) {
#       stop("Treatment options and mean matrix dimension do not match!")
#     }
#     if (length(class.A)!= ncol(mus.hat)) {
#       mus.hat = mus.hat[,class.A]
#     }
#
#     # pick a single best treatment for all patients
#     c.means<-apply(mus.hat,2,mean)
#     Ey.opt1<-max(c.means)
#     trt.opt1<-class.A[which(c.means==Ey.opt1)]###Problem: if the order of class.A is different from mus.hat because of the sort() this will cause problems
#   }
#   outs<-list(Ey.opt1, trt.opt1)
#   names(outs)<-c("Ey.opt1", "trt.opt1")
#   # outs = as.data.frame(matrix(NA,1,2))
#   # names(outs) <- c("Ey.opt1","trt.opt1")
#   # output[1,1] = Ey.opt1
#   # output[1,2] = trt.opt1
#   return(outs)
# }

# id = which(all.A %in% tol.A)
