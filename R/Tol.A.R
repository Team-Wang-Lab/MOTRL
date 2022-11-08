#' The set of tolerant treatments with given tolerant rate and pseudo outcome matrix in one given stage.
#'
#' @param A Treatment options given in this stage
#' @param mus.hat A matrix indicates the final counterfactual outcome according to the treatments.
#' @param tol.rate A scalar indicates the tolerant rate
Tol.A <- function(A, mus.hat, tol.rate) {
  class.A <- sort(unique(A))
  # If there are only one treatment option
  if (length(class.A) == 1) {
    trt.opt1 <- class.A
    Ey.opt1 <- mean(mus.hat)
  } else {
    # if (length(A)!= nrow(mus.hat) || length(class.A)!= ncol(mus.hat)) {
    #   stop("Treatment options and mean matrix dimension do not match!")
    # }
    # Modification:
    if (length(A)!= nrow(mus.hat)) {
      stop("Treatment options and mean matrix dimension do not match!")
    }
    if (length(class.A)!= ncol(mus.hat)) {
      stop("Treatment options in one node do not inclued all treatment options! Try a bigger minsplit")
      # mus.hat = mus.hat[,class.A]
      # mus.hat = mus.hat[,which(all.A %in% class.A)]
    }

    # pick a single best treatment for all patients
    c.means <- apply(mus.hat, 2, mean)
    Ey.opt1 <- max(c.means)
    trt.opt1<-class.A[which(c.means == Ey.opt1)] ### Problem: if the order of class.A is different from mus.hat because of the sort() this will cause problems

    # pick other treatments that in the tolerant range
    if (tol.rate == 1) { # The most strict case
      Ey.tol.avg1 = Ey.opt1
      trt.tol1 = trt.opt1
    } else {
      Ey.min1 = min(c.means)
      Ey.range = Ey.opt1 - Ey.min1
      Ey.tol.cut = Ey.min1 + tol.rate*Ey.range
      Ey.tol.avg1 = mean(c.means[which(c.means >= Ey.tol.cut)])
      trt.tol1 = class.A[which(c.means >= Ey.tol.cut)]
    }
  }
  outs<-list(Ey.opt1, trt.opt1, trt.tol1, Ey.tol.avg1)
  names(outs)<-c("Ey.opt1", "trt.opt1", "trt.tol1", "Ey.tol.avg1")
  # outs = as.data.frame(matrix(NA,1,4))
  # names(outs) <- c("Ey.opt1", "trt.opt1", "trt.tol1", "Ey.tol.avg1")
  # output[1,1] = Ey.opt1
  # output[1,2] = trt.opt1
  # output[1,3] = trt.tol1
  # output[1,4] = Ey.tol.avg1
  return(outs)
}








