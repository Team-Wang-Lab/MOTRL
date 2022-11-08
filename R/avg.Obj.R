#' ??
#'
#' ??
#'
#' @param all.A A vector of all treatments that is sorted???????
#' @param tol.A A list of treatments that are tolerated
#' @param list.mus.hat A list of matrix of estimated conditional mean outcomes. Each element matrix is for one objectives,
#'
#' @import dplyr
#'
#' @return a vector indicates the average pseudo outcome of each objectives
#'
avg.Obj <- function(all.A, tol.A, list.mus.hat) {
  # class.A <- sort(unique(A))
  # Obj.dim = length(list.mus.hat)
  #
  # if (length(class.A) == 1) {
  #   E.Objs <- unlist(lapply(list.mus.hat, function(x) mean(x[,class.A])))
  # }

  id = which(all.A %in% tol.A)
  # id = which(tol.A %in% all.A)
  E.Objs <- unlist(lapply(list.mus.hat, function(x) mean(x[,id])))
  return(E.Objs)
}



# avg.Obj <- function(all.A, tol.A, list.mus.hat) {
#
#   n = length(tol.A)
#   Obj.dim = length(list.mus.hat)
#   avg.PO = matrix(0, nrow = n, ncol = Obj.dim)
#
#   for (i in 1L:n) {
#     for (d in 1L:Obj.dim) {
#       avg.PO[i,d] = mean(list.mus.hat[[d]][i,which(class.A %in% unlist(tol.regimes[i]))])
#     }
#   }
#
#
#   id = which(all.A %in% tol.A)
#   # id = which(tol.A %in% all.A)
#   E.Objs <- unlist(lapply(list.mus.hat, function(x) mean(x[,id])))
#   return(E.Objs)
# }







