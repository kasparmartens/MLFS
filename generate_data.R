rsparsematrix = function(m, n, binary_vector){
  A = matrix(rnorm(m*n), m, n)
  return(A * binary_vector)
}


generate_data = function(H, C = 2, N = 100, d = c(5, 5, 5, 5)){
  R = nrow(H)
  
  V = matrix(rnorm(N*R), N, R)
  
  W_list = lapply(1:ncol(H), function(j)rsparsematrix(R, d[j], H[, j]))
  
  
  U_list = lapply(W_list, function(W){
    V %*% W
  })
  
  beta = matrix(rnorm(R*C), R, C)
  z = V %*% beta
  y = apply(z, 1, which.max)
  
  return(list(U_list = U_list, y = y, W_list = W_list))
}


# W = do.call("cbind", W_list)


# Generate X
generate_X = function(U_list, type, g_list = NULL){
  X_list = list()
  for(j in 1:length(type)){
    if(type[j] == "gaussian"){
      X_list[[j]] = U_list[[j]]
    }
    else if(type[j] == "ordinal"){
      u = U_list[[j]]
      g = g_list[[j]]
      temp = matrix(NA, nrow(u), ncol(u))
      for(k in 1:(length(g)-1)){
        condition = (u > g[k]) & (u < g[k+1])
        temp[condition] = k
      }
      if(sum(is.na(temp)) > 0) stop("some u values outside of range (-G, G)")
      X_list[[j]] = temp
    }
  }
  return(X_list)
}

