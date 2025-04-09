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

wilson_lower <- function(phat, n, alpha = 0.05, z = NULL){
  # compute lower bound of wilson score interval 
  # phat : weighted sample mean (proportion)
  # n: (effective) sample size
  # alpha : confidence level 
  if (is.null(z))   z <-  qnorm(1-alpha/2)
  
  if (is.na(phat)){
    return(NA)
  } else if (phat==0) {
    return(0)
  } else if (phat==1){
    return(n/(n+z^2))
  } else {
    return(
      (1/(1+(1/n)*(z^2)))*(phat + (z^2)/(2*n) - ((z)/(2*n))*sqrt(4*n*phat*(1-phat) +z^2))
    ) 
  }
}
wilson_lower <- Vectorize(wilson_lower)
wilson_upper <- function(phat, n,  alpha = 0.05, z = NULL){
  # compute upper bound of wilson score interval 
  # phat : (weighted) sample mean (proportion)
  # n: (effective) sample size 
  # alpha : confidence level 
  if (is.null(z))   z <-  qnorm(1-alpha/2)
  if (is.na(phat)){
    return(NA)
  } else {
    if (phat==0) {
      return(z^2 / (n+z^2))
    } else if (phat==1) {
      return(1)
    } else{
      return(
        (1/(1+(1/n)*(z^2)))*(phat + (z^2)/(2*n) + ((z)/(2*n))*sqrt(4*n*phat*(1-phat) + z^2))
      ) 
    }
  }
}
wilson_upper <- Vectorize(wilson_upper)

compute_wald <- function(est,se, z){
  return(est + z*se)
}
compute_wilson <- function(est,n,z){
  if (z<0) wilson_lower(phat=est,n=n,z=-z) else wilson_upper(phat=est, n=n, z=z)
}

check_overlap<- function(xmin, xmax, ymin, ymax) {
  # Ensure inputs are vectors of the same length
  if (length(xmin) != length(xmax) || length(ymin) != length(ymax)) {
    stop("Input vectors must have the same length.")
  }
  
  # Check for NA values
  overlap <- ifelse(
    is.na(xmin) | is.na(xmax) | is.na(ymin) | is.na(ymax),
    NA, # Return NA if any input in the pair is NA
    !(xmax < ymin | ymax < xmin) # Otherwise, compute overlap
  )
  
  return(overlap)
}

