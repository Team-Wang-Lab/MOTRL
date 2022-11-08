#' Predict tolerant treatment using the output from MO.tol.DTRtree
#'
#' @param treeout An object (tree) outputted from MO.tol.DTRtree function.
#' @param newdata New data containing H (history).
#' @export

predict_tol.DTR<-function(treeout,newdata){
  n <- nrow(newdata)
  # predicts <- data.frame(NA, n, 2)
  # predicts = as.data.frame(matrix(NA,n,1))
  predicts <- vector(mode='list', length = n)
  # names(predicts) = c("avg.PO")

  # treeout is supposed to be a matrix
  # if there is no split
  if (nrow(treeout) == 1) {
    for (i in 1:n) {
      predicts[i] = treeout$tol.trt
    }
  } else { # if there are splits
    newdata <- as.data.frame(newdata)

    for (i in 1:n) {
      nd <- 1
      while (is.na(treeout$opt.trt[treeout$node == nd])) { # if this node have NO opt.trt assigned
        # Modification : if split is numerical
        # if (is.numeric(treeout$cutoff[nd][[1]])) {
        if (is.numeric(newdata[,as.numeric(treeout$X[treeout$node == nd])])) {
          if (newdata[i,as.numeric(treeout$X[treeout$node==nd])] <= unlist(treeout$cutoff[treeout$node==nd])) {#if the node<= cutoff
            nd = 2 * nd # yes proceed first
          } else {
            nd = 2 * nd + 1 # then no
          }
        } else {
          if (newdata[i,as.numeric(treeout$X[treeout$node==nd])] %in% unlist(treeout$cutoff[treeout$node==nd][[1]])) {#if the node<= cutoff
            nd = 2 * nd # yes proceed first
          } else {
            nd = 2 * nd + 1 # then no
          }
        }
      }
      predicts[i] = treeout$tol.trt[treeout$node==nd]
    }

  }
  return(predicts)
}


