library(ggplot2)
library(dplyr)
library(reshape2)

source("helpers.R")
source("par_updates.R")
source("rotation.R")
source("ordinal_cutpoints.R")
source("MLFS.R")
source("generate_data.R")

# H = rbind(c(1, 1), 
#           c(1, 1), 
#           c(1, 0), 
#           c(1, 0), 
#           c(0, 1), 
#           c(0, 1))
# 
# set.seed(0)
# R = nrow(H)
# d = c(5, 5)
# type = c("gaussian", "ordinal")

# Perfect scenario
experiment0 = function(H, d, type, R, n_rep = 3, n_levels = 3, rotate=TRUE){
  pred_acc_train = rep(NA, n_rep)
  pred_acc = rep(NA, n_rep)
  for(i in 1:n_rep){
    # generate latent variables
    latent_data = generate_latent_subspace(H, N = 200, d = d)
    y = generate_y(latent_data$V, C = 2)
    g_list = lapply(latent_data$U_list, function(x)c(-Inf, quantile(x, 1:(n_levels-1) / n_levels), Inf))

    X0 = generate_X(latent_data$U_list, type, g_list)
    data = split_into_train_and_test(X0, y, prop=0.5)
    MLFSobj = MLFS(data$trainy, data$trainX, type, R, max_iter=20, rotate=rotate, verbose=FALSE)
    pred = pred_out_of_sample(data$testX, MLFSobj)
    pred_acc[i] = mean(pred == data$testy)
    pred_acc_train[i] = MLFSobj$pred_acc_train
  }
  df = data.frame(pred_acc, pred_acc_train)
  names(df) = c("test", "train")
  return(df)
}

# repeat experiment
three_experiments = function(experiment, H, d, R, n_rep = 5, n_views = 2){
  type = rep(c("gaussian", "gaussian"), n_views/2)
  res1 = experiment(H, d, type, R, n_rep = n_rep)
  type = rep(c("gaussian", "ordinal"), n_views/2)
  res2 = experiment(H, d, type, R, n_rep = n_rep)
  type = rep(c("ordinal", "ordinal"), n_views/2)
  res3 = experiment(H, d, type, R, n_rep = n_rep)
  df = data.frame(rbind(res1, res2, res3), 
                   type = rep(c("gaussian + gaussian", "gaussian + ordinal", "ordinal + ordinal"), each=nrow(res1)))
  return(df)
}

# Convergence of lower bound
experiment_lowerbound = function(H, d, type, n_levels = 3){
  # generate latent variables
  latent_data = generate_latent_subspace(H, N = 200, d = d)
  y = generate_y(latent_data$V, C = 2)
  g_list = lapply(latent_data$U_list, function(x)c(-Inf, quantile(x, 1:(n_levels-1) / n_levels), Inf))
  
  X0 = generate_X(latent_data$U_list, type, g_list)
  data = split_into_train_and_test(X0, y, prop=0.5)
  MLFSobj1 = MLFS(data$trainy, data$trainX, type, R, max_iter=25, rotate=TRUE, verbose=FALSE)
  MLFSobj2 = MLFS(data$trainy, data$trainX, type, R, max_iter=25, rotate=FALSE, verbose=FALSE)
  bound1 = MLFSobj1$lowerbound
  bound2 = MLFSobj2$lowerbound
  
  df = data.frame(lowerbound = c(bound1, bound2), iterations = rep(1:length(bound1), 2), rotate = rep(c("TRUE", "FALSE"), each=length(bound1)))
  return(df)
}

# Experiment with gamma
experiment_gamma = function(gammas, H, d, type, R, n_rep = 3, n_levels = 3, rotate=TRUE){
  df = data.frame()
  for(k in 1:length(gammas)){
    pred_acc_train = rep(NA, n_rep)
    pred_acc = rep(NA, n_rep)
    for(i in 1:n_rep){
      # generate latent variables
      latent_data = generate_latent_subspace(H, N = 200, d = d, gamma = gammas[k])
      y = generate_y(latent_data$V, C = 2)
      g_list = lapply(latent_data$U_list, function(x)c(-Inf, quantile(x, 1:(n_levels-1) / n_levels), Inf))
      
      X0 = generate_X(latent_data$U_list, type, g_list)
      data = split_into_train_and_test(X0, y, prop=0.5)
      MLFSobj = MLFS(data$trainy, data$trainX, type, R, max_iter=20, rotate=rotate, verbose=FALSE)
      pred = pred_out_of_sample(data$testX, MLFSobj)
      pred_acc[i] = mean(pred == data$testy)
      pred_acc_train[i] = MLFSobj$pred_acc_train
    }
    temp = data.frame(pred_acc, pred_acc_train, gamma = gammas[k])
    df = rbind(df, temp)
  }
  names(df)[1:2] = c("test", "train")
  return(df)
}


# Experiment with the number of irrelevant features
experiment1 = function(n_irrelevant_features, H, d, type, n_rep = 3, n_levels = 3){
  K = length(n_irrelevant_features)
  pred_acc_train = rep(NA, K)
  pred_acc = rep(NA, K)
  df = data.frame()
  for(i in 1:n_rep){
    # generate latent variables
    latent_data = generate_latent_subspace(H, N = 200, d = d)
    y = generate_y(latent_data$V, C = 2)
    g_list = lapply(latent_data$U_list, function(x)c(-Inf, quantile(x, 1:(n_levels-1) / n_levels), Inf))
    
    for(k in 1:K){
      X0 = generate_X(latent_data$U_list, type, g_list)
      X1 = add_irrelevant_features(X0, type, n_irrelevant_features[k])
      data = split_into_train_and_test(X1, y, prop=0.5)
      MLFSobj = MLFS(data$trainy, data$trainX, type, R, max_iter=20, rotate=TRUE, verbose=FALSE)
      pred = pred_out_of_sample(data$testX, MLFSobj)
      pred_acc[k] = mean(pred == data$testy)
      pred_acc_train[k] = MLFSobj$pred_acc_train
    }
    temp = data.frame(n_irrelevant_features, pred_acc, pred_acc_train)
    df = rbind(df, temp)
  }
  names(df)[2:3] = c("test", "train")
  return(df)
}

# Experiment with the dimensionality of latent space
experiment2 = function(R_values, H, d, type, n_rep = 3, n_levels = 3){
  K = length(R_values)
  pred_acc_train = rep(NA, K)
  pred_acc = rep(NA, K)
  df = data.frame()
  for(i in 1:n_rep){
    # generate latent variables
    latent_data = generate_latent_subspace(H, N = 200, d = d)
    y = generate_y(latent_data$V, C = 2)
    g_list = lapply(latent_data$U_list, function(x)c(-Inf, quantile(x, 1:(n_levels-1) / n_levels), Inf))
    
    for(k in 1:K){
      X0 = generate_X(latent_data$U_list, type, g_list)
      data = split_into_train_and_test(X0, y, prop=0.5)
      MLFSobj = MLFS(data$trainy, data$trainX, type, R_values[k], max_iter=30, rotate=FALSE, verbose=FALSE)
      pred = pred_out_of_sample(data$testX, MLFSobj)
      pred_acc[k] = mean(pred == data$testy)
      pred_acc_train[k] = MLFSobj$pred_acc_train
    }
    temp = data.frame(R_values, pred_acc, pred_acc_train)
    df = rbind(df, temp)
  }
  names(df) = c("latent_dim", "test", "train")
  return(df)
}

# Multiple correlated features

rsparsematrix_mod = function(n_added, m, n, binary_vector, rho=0.9){
  A = matrix(rnorm(m*n), m, n)
  B = matrix(0, m, n_added*n)
  for(k in 1:n_added){
    for(j in 1:n){
      B[, (k-1)*n+j] = rho*A[, j] + rnorm(m)
      # print(cor(B[, (k-1)*n+j], A[, j]))
    }
  }
  out = cbind(A, B)
  return(out*binary_vector)
}


generate_latent_subspace_with_correlated_features = function(n_added, H, N = 100, d = c(5, 5, 5, 5)){
  R = nrow(H)
  V = matrix(rnorm(N*R), N, R)
  W_list = lapply(1:length(d), function(j)rsparsematrix_mod(n_added, R, d[j], H[, j]))
  U_list = lapply(W_list, function(W) V %*% W)
  return(list(U_list = U_list, W_list = W_list, V = V, beta = beta))
}

experiment_correlated_features = function(n_added, H, d, type, n_rep = 5, n_levels = 3){
  K = length(n_added)
  pred_acc_train = rep(NA, K)
  pred_acc = rep(NA, K)
  df = data.frame()
  for(i in 1:n_rep){
    for(k in 1:K){
      # generate latent variables
      latent_data = generate_latent_subspace_with_correlated_features(n_added[k], H, N = 200, d = d)
      y = generate_y(latent_data$V, C = 2)
      g_list = lapply(latent_data$U_list, function(x)c(-Inf, quantile(x, 1:(n_levels-1) / n_levels), Inf))
      X0 = generate_X(latent_data$U_list, type, g_list)
      data = split_into_train_and_test(X0, y, prop=0.5)
      MLFSobj = MLFS(data$trainy, data$trainX, type, R, max_iter=20, rotate=TRUE, verbose=FALSE)
      pred = pred_out_of_sample(data$testX, MLFSobj)
      pred_acc[k] = mean(pred == data$testy)
      pred_acc_train[k] = MLFSobj$pred_acc_train
    }
    temp = data.frame(n_added, pred_acc, pred_acc_train)
    df = rbind(df, temp)
  }
  names(df)[2:3] = c("test", "train")
  return(df)
}


generate_correlated_features = function(n_added, U_list, rho=0.9){
  out = lapply(U_list, function(U){
    B = matrix(0, nrow(U), ncol(U)*n_added)
    for(k in 1:n_added){
      for(j in 1:ncol(U)){
        B[, (k-1)*ncol(U)+j] = rho*U[, j] + rnorm(nrow(U))
        print(cor(U[, j], B[, (k-1)*ncol(U)+j]))
      }
    }
    return(cbind(U, B))
  })
  return(out)
}

experiment_correlated_features2 = function(n_added, H, d, type, n_rep = 5, n_levels = 3){
  K = length(n_added)
  pred_acc_train = rep(NA, K)
  pred_acc = rep(NA, K)
  df = data.frame()
  for(i in 1:n_rep){
    for(k in 1:K){
      # generate latent variables
      latent_data = generate_latent_subspace(H, N = 200, d = d)
      y = generate_y(latent_data$V, C = 2)
      g_list = lapply(latent_data$U_list, function(x)c(-Inf, quantile(x, 1:(n_levels-1) / n_levels), Inf))
      latent_data$U_list = generate_correlated_features(n_added[k], latent_data$U_list)
      X0 = generate_X(latent_data$U_list, type, g_list)
      data = split_into_train_and_test(X0, y, prop=0.5)
      MLFSobj = MLFS(data$trainy, data$trainX, type, R, max_iter=20, rotate=TRUE, verbose=FALSE)
      pred = pred_out_of_sample(data$testX, MLFSobj)
      pred_acc[k] = mean(pred == data$testy)
      pred_acc_train[k] = MLFSobj$pred_acc_train
    }
    temp = data.frame(n_added, pred_acc, pred_acc_train)
    df = rbind(df, temp)
  }
  names(df)[2:3] = c("test", "train")
  return(df)
}

generate_latent_subspace_t = function(H, N = 100, d = c(5, 5, 5, 5), gamma=1000){
  R = nrow(H)
  
  V = matrix(rnorm(N*R), N, R)
  W_list = lapply(1:length(d), function(j)rsparsematrix(R, d[j], H[, j]))
  
  U_list = lapply(W_list, function(W) V %*% W + 1/sqrt(gamma)*rt(N*ncol(W), 1))
  
  return(list(U_list = U_list, W_list = W_list, V = V, beta = beta))
}

# Experiment with gaussianity
experiment_gaussianity = function(H, d, type, R, n_rep = 5, n_levels = 3, rotate=TRUE){
  pred_acc_train = rep(NA, n_rep)
  pred_acc = rep(NA, n_rep)
  for(i in 1:n_rep){
    # generate latent variables
    latent_data = generate_latent_subspace_t(H, N = 200, d = d, gamma = 1)
    y = generate_y(latent_data$V, C = 2)
    g_list = lapply(latent_data$U_list, function(x)c(-Inf, quantile(x, 1:(n_levels-1) / n_levels), Inf))
    
    X0 = generate_X(latent_data$U_list, type, g_list)
    data = split_into_train_and_test(X0, y, prop=0.5)
    MLFSobj = MLFS(data$trainy, data$trainX, type, R, max_iter=20, rotate=rotate, verbose=FALSE)
    pred = pred_out_of_sample(data$testX, MLFSobj)
    pred_acc[i] = mean(pred == data$testy)
    pred_acc_train[i] = MLFSobj$pred_acc_train
  }
  df = data.frame(pred_acc, pred_acc_train)
  names(df)[1:2] = c("test", "train")
  return(df)
}