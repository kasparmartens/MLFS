library(mvtnorm)
library(truncnorm)

MLFS_mcmc = function(y, X_list, type, R, max_iter=10, rotate=TRUE, d_sim = 5, verbose = TRUE){
  N = length(y)
  M = length(X_list)
  d = ifelse(type != "similarity", sapply(X_list, ncol), d_sim)
  C = max(y)
  if(sum(!(y %in% 1:max(y))) > 0) stop("y must have labels 1, 2, ..., C")
  if(sum(sapply(X_list, class) != "matrix") > 0) stop("X_list must contain matrices only!")
  
  a0 = 0.01
  b0 = 0.01
  a_gamma = 0.01
  b_gamma = 0.01
  rho = 0.01
  
  # alpha, gamma
  alpha = matrix(1, M, R)
  gamma = rep(1, M)
  
  # initialise pi, H, W
  pi = rbeta(M, 1, 1)
  H = matrix(runif(M*R) < pi, M, R)
  W = list()
  for(j in 1:M){
    fill = rep(H[j, ], d[j])*rnorm(d[j]*R)
    W[[j]] = matrix(fill, R, d[j])
  }
  
  # initialise V
  V = matrix(rnorm(N*R), N, R)
  
  # initalise U
  U = list()
  for(j in 1:M){
    if(type[j] == "gaussian") U[[j]] = X_list[[j]]
  }
  
  # initialise z, beta
  z = rep(0, N)
  beta = rep(0, R)
  beta_trace = matrix(NA, max_iter, R)
  pred_trace = matrix(NA, max_iter, N)
  
  for(iter in 1:max_iter){
    temp = update_H_and_W(U, V, W, alpha, gamma, pi, d, M, R)
    H = temp$H
    W = temp$W
    
    # update pi
    H_rowsums = rowSums(H)
    pi = rbeta(M, 0.01 + H_rowsums, 0.01 + R - H_rowsums)
    
    ### update alpha
    for(j in 1:M){
      temp0 = rgamma(d[j], a0, b0)
      temp1 = rgamma(d[j], a0 + 0.5*d[j], b0 + 0.5*sum(W[[j]]**2))
      alpha[j, ] = ifelse(H[j, ] == 1, temp1, temp0)
    }
    
    ### update gamma
    for(j in 1:M){
      temp = norm(U[[j]] - V %*% W[[j]], "F")
      gamma[j] = rgamma(1, a_gamma + 0.5*d[j]*N, b_gamma + 0.5*temp**2)
    }
    
    ### sample V
    sigmainv_V = diag(R) + matrix_list_sum(lapply(1:M, function(j){
      gamma[j] * W[[j]] %*% t(W[[j]])
    }))
    sigma_V = solve(sigmainv_V)
    
    
    for(i in 1:N){
      mu_tmp = update_V_mu_individual_i(i, U, W, gamma, M)
      mu_i = (z[i]*beta + mu_tmp) %*% sigma_V
      V[i, ] = rmvnorm(1, mu_i, sigma_V)
    }
    
    ### For classification follow Albert & Chib 1993
    # sample beta
    sigma_beta = solve(rho * diag(R) + t(V) %*% V)
    mu_beta = sigma_beta %*% t(V) %*% z
    beta = as.numeric(rmvnorm(1, mu_beta, sigma_beta))
    beta_trace[iter, ] = beta
      
    # sample z
    z_mu = V %*% beta
    subset = (y == 2)
    if(sum(subset)>0) z[subset] = rtruncnorm(sum(subset), a=0, mean=z_mu[subset], sd=1)
    subset = (y == 1)
    if(sum(subset)>0) z[subset] = rtruncnorm(sum(subset), b=0, mean=z_mu[subset], sd=1)
    
    pred = (z_mu > 0) + 1
    pred_trace[iter, ] = pred
    
    if(iter %% 10 == 0){
      cat("pred acc:", mean(pred == y), "\n")
    }
  }
  return(pred_trace)
}

### update H
# update_H_and_W = function(U, V, W, alpha, gamma, pi, d, M, R){
#   H = matrix(0, M, R)
#   for(j in 1:M){
#     s2 = 1/alpha[j, ] + gamma[j] * diag(t(V) %*% V) # sapply(1:ncol(V), function(r)t(V[, r]) %*% V[, r])
#     Ures = U[[j]] - V %*% W[[j]]
#     mu = gamma[j] * t(Ures) %*% V %*% diag(1 / s2)
#     z = logit(pi[j]) + 0.5*d[j]*log(s2 * alpha[j, ]) + 0.5*diag(t(mu) %*% mu) / s2
#     acceptance_probs = 1 / (1 + exp(-z))
#     print(acceptance_probs)
#     u = runif(R)
#     H[j, ] = ifelse(u < acceptance_probs, 1, 0)
#     W_fill = mu + rnorm(R*d[j])*sqrt(rep(s2, each=d[j]))
#     W[[j]] = matrix(W_fill * rep(H[j, ], d[j]), R, d[j])
#   }
#   return(list(H = H, W = W))
# }

update_H_and_W = function(U, V, W, alpha, gamma, pi, d, M, R){
  H = matrix(0, M, R)
  for(j in 1:M){
    for(r in 1:R){
      s2 = 1/(alpha[j, r] + gamma[j] * sum(V[, r] * V[, r]))
      Ures = U[[j]] - V %*% W[[j]]
      mu = gamma[j] * s2 * t(Ures) %*% V[, r]
      z = logit(pi[j]) + 0.5*d[j]*log(s2 * alpha[j, r]) + 0.5*diag(t(mu) %*% mu) / s2
      acceptance_prob = 1 / (1 + exp(-z))
      u = runif(1)
      H[j, r] = ifelse(u < acceptance_prob, 1, 0)
      W_fill = mu + rnorm(d[j])*sqrt(s2)
      W[[j]][r, ] = W_fill * H[j, r]
    }
  }
  return(list(H = H, W = W))
}




