#' Predict optimal treatment using the output from DTRtree.
#'
#' @param treeout An object outputted from DTRtree function.
#' @param newdata New data containing H (history).
#' @export

predict_opt.DTR <- function(treeout, newdata){
  # n<-nrow(newdata)
  # predicts<-rep(NA,n)
  n <- nrow(newdata)
  # predicts <- data.frame(NA, n, 2)
  # predicts = as.data.frame(matrix(NA,n,1))
  predicts <- vector(mode='list', length = n)
  # names(predicts) = c("avg.PO")


  # treeout is supposed to be a matrix
  # if there is no split
  if (nrow(treeout) == 1) {
    predicts <- rep(treeout$opt[1], n)
  } else { # if there are splits
    newdata <- as.data.frame(newdata)
  # treeout is supposed to be a matrix
  # if there is no split
  # if(length(treeout)==5){
  #   predicts<-rep(treeout[5],n)
  # } else{ # if there are splits
  #   treeout<-as.data.frame(treeout)
  #   newdata<-as.data.frame(newdata)

    # for(i in 1:n){
    #   nd<-1
    #   while(is.na(treeout$trt[treeout$node==nd])){
    #     if(newdata[i,treeout$X[treeout$node==nd]] <= treeout$cutoff[treeout$node==nd]){#if the node<= cutoff
    #       nd=2*nd #yes proceed first
    #     } else{
    #       nd=2*nd+1#then no
    #     }
    #   }
    #   predicts[i]<-treeout$trt[treeout$node==nd]
    # }

    for (i in 1:n) {
      nd <- 1
      while (is.na(treeout$opt.trt[treeout$node == nd])) { # if this node have NO opt.trt assigned
        # Modification : if split is numerical
        if (is.numeric(newdata[,as.numeric(treeout$X[nd])])) {
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
      predicts[i] = treeout$opt.trt[treeout$node==nd]
    }
  }
  return(predicts)
}

