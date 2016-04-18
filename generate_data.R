rsparsematrix = function(m, n, binary_vector){
  A = matrix(rnorm(m*n), m, n)
  return(A * binary_vector)
}


generate_latent_subspace = function(H, N = 100, d = c(5, 5, 5, 5)){
  R = nrow(H)
  
  V = matrix(rnorm(N*R), N, R)
  
  W_list = lapply(1:ncol(H), function(j)rsparsematrix(R, d[j], H[, j]))
  
  U_list = lapply(W_list, function(W) V %*% W)
  
  return(list(U_list = U_list, W_list = W_list, V = V, beta = beta))
}

generate_y = function(V, C = 2){
  R = ncol(V)
  beta = matrix(rnorm(R*C), R, C)
  z = V %*% beta
  y = apply(z, 1, which.max)
  return(y)
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
        condition = (u > g[k]) & (u <= g[k+1])
        temp[condition] = k
      }
      if(sum(is.na(temp)) > 0) stop("some u values outside of range (-G, G)")
      X_list[[j]] = temp
    }
  }
  return(X_list)
}

split_into_train_and_test = function(X_list, y, prop = 0.8){
  N = nrow(X_list[[1]])
  trainind = sample(1:N, prop*N)
  trainX = lapply(X_list, function(mat)mat[trainind, ])
  testX = lapply(X_list, function(mat)mat[-trainind, ])
  trainy = y[trainind]
  testy = y[-trainind]
  return(list(trainX = trainX, trainy = trainy, testX = testX, testy = testy))
}

add_irrelevant_features = function(X_list, type, n_features){
  N = nrow(X_list[[1]])
  p = n_features
  for(j in 1:length(type)){
    U = matrix(rnorm(N*p), N, p)
    if(type[j] == "gaussian"){
      X_list[[j]] = cbind(X_list[[j]], U)
    }
    else if(type[j] == "ordinal"){
      n_levels = max(X_list[[j]])
      g = c(-Inf, quantile(U, 1:(n_levels-1) / n_levels), Inf)
      temp = matrix(NA, N, p)
      for(k in 1:(length(g)-1)){
        condition = (U > g[k]) & (U <= g[k+1])
        temp[condition] = k
      }
      if(sum(is.na(temp)) > 0) stop("some u values outside of range (-G, G)")
      X_list[[j]] = cbind(X_list[[j]], temp)
    }
  }
  return(X_list)
}
