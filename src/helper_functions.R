nb2mat_to_nodes <- function(adj_mat) {
  N_A <- nrow(adj_mat)
  N_edges <- sum(adj_mat != 0) / 2
  n1 <- vector(mode="numeric", length = N_edges)
  n2 <- vector(mode="numeric", length = N_edges)
  k <- 1
  for (i in 1:N_A) {
    for (j in i:N_A) {
      if (adj_mat[i, j] != 0) {
        n1[k] <- i
        n2[k] <- j
        k <- k + 1
      }
    }
  }
  return(list(n1 = n1, n2 = n2))
}

prepare_bym2 <- function(adj_mat) {
  nodes <- nb2mat_to_nodes(adj_mat)
  inla_adj <- sparseMatrix(i = nodes$n1, j = nodes$n2,
                           x = 1, symmetric = T)
  # ICAR precision matrix
  Q <- Diagonal(nrow(adj_mat), Matrix::rowSums(inla_adj)) - inla_adj
  Q_jit = Q + Diagonal(nrow(adj_mat)) * max(diag(Q)) * sqrt(.Machine$double.eps)
  
  Q_inv = inla.qinv(Q_jit, constr=list(A = matrix(1, 1, nrow(adj_mat)), e=0))
  
  #Compute the geometric mean of the variances, which are on the diagonal of Q.inv
  scl = exp(mean(log(diag(Q_inv))))
  return(list(n1 = nodes$n1, n2 = nodes$n2, scaling_factor = scl))
}