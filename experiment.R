library(ggplot2)
library(dplyr)
library(reshape2)

source("helpers.R")
source("par_updates.R")
source("rotation.R")
source("ordinal_cutpoints.R")
source("MLFS.R")
source("MLFS_regression.R")
source("MLFS_mcmc.R")
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
experiment0 = function(H, d, type, R, n_rep = 3, n_levels = 3, rotate=TRUE, continuous = FALSE){
  pred_acc_train = rep(NA, n_rep)
  pred_acc = rep(NA, n_rep)
  for(i in 1:n_rep){
    # generate latent variables
    latent_data = generate_latent_subspace(H, N = 200, d = d)
    y = generate_y(latent_data$V, C = 2, continuous = continuous)
    g_list = lapply(latent_data$U_list, function(x)c(-Inf, quantile(x, 1:(n_levels-1) / n_levels), Inf))

    X0 = generate_X(latent_data$U_list, type, g_list)
    data = split_into_train_and_test(X0, y, prop=0.5)
    if(continuous){
      MLFSobj = MLFS_regression(data$trainy, data$trainX, type, R, max_iter=20, rotate=rotate, verbose=FALSE)
      pred = pred_out_of_sample_regression(data$testX, MLFSobj)
      pred_acc[i] = cor(pred, data$testy)**2
      pred_acc_train[i] = MLFSobj$pred_acc_train      
    } else{
      MLFSobj = MLFS(data$trainy, data$trainX, type, R, max_iter=20, rotate=rotate, verbose=FALSE)
      pred = pred_out_of_sample(data$testX, MLFSobj)
      pred_acc[i] = mean(pred == data$testy)
      pred_acc_train[i] = MLFSobj$pred_acc_train
    }
  }
  df = data.frame(pred_acc, pred_acc_train)
  names(df) = c("test", "train")
  return(df)
}

# repeat experiment
three_experiments = function(experiment, H, d, R, n_views = 2, ...){
  type = rep(c("gaussian", "gaussian"), n_views/2)
  res1 = experiment(H, d, type, R, ...)
  type = rep(c("gaussian", "ordinal"), n_views/2)
  res2 = experiment(H, d, type, R, ...)
  type = rep(c("ordinal", "ordinal"), n_views/2)
  res3 = experiment(H, d, type, R, ...)
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
experiment_gamma = function(gammas, H, d, type, R, n_rep = 3, n_levels = 3, rotate=TRUE, continuous = FALSE){
  df = data.frame()
  for(k in 1:length(gammas)){
    pred_acc_train = rep(NA, n_rep)
    pred_acc = rep(NA, n_rep)
    for(i in 1:n_rep){
      # generate latent variables
      latent_data = generate_latent_subspace(H, N = 200, d = d, gamma = gammas[k])
      y = generate_y(latent_data$V, C = 2, continuous = continuous)
      g_list = lapply(latent_data$U_list, function(x)c(-Inf, quantile(x, 1:(n_levels-1) / n_levels), Inf))
      
      X0 = generate_X(latent_data$U_list, type, g_list)
      data = split_into_train_and_test(X0, y, prop=0.5)
      if(continuous){
        MLFSobj = MLFS_regression(data$trainy, data$trainX, type, R, max_iter=20, rotate=rotate, verbose=FALSE)
        pred = pred_out_of_sample_regression(data$testX, MLFSobj)
        pred_acc[i] = cor(pred, data$testy)**2
        pred_acc_train[i] = MLFSobj$pred_acc_train      
      } else{
        MLFSobj = MLFS(data$trainy, data$trainX, type, R, max_iter=20, rotate=rotate, verbose=FALSE)
        pred = pred_out_of_sample(data$testX, MLFSobj)
        pred_acc[i] = mean(pred == data$testy)
        pred_acc_train[i] = MLFSobj$pred_acc_train
      }
    }
    temp = data.frame(pred_acc, pred_acc_train, gamma = gammas[k])
    df = rbind(df, temp)
  }
  names(df)[1:2] = c("test", "train")
  return(df)
}


# Experiment with the number of irrelevant features
library(doParallel)
library(foreach)
registerDoParallel(cores = 8)
experiment1 = function(n_irrelevant_features, H, d, type, n_rep = 3, n_levels = 3, continuous = FALSE){
  K = length(n_irrelevant_features)
  pred_acc_train = rep(NA, K)
  pred_acc = rep(NA, K)
  df = data.frame()
  df = foreach(i = 1:n_rep, .combine="rbind") %dopar% {
    # generate latent variables
    latent_data = generate_latent_subspace(H, N = 200, d = d)
    y = generate_y(latent_data$V, C = 2, continuous = continuous)
    g_list = lapply(latent_data$U_list, function(x)c(-Inf, quantile(x, 1:(n_levels-1) / n_levels), Inf))
    for(k in 1:K){
      X0 = generate_X(latent_data$U_list, type, g_list)
      X1 = add_irrelevant_features(X0, type, n_irrelevant_features[k])
      data = split_into_train_and_test(X1, y, prop=0.5)
      if(continuous){
        MLFSobj = MLFS_regression(data$trainy, data$trainX, type, R, max_iter=20, rotate=TRUE, verbose=FALSE)
        pred = pred_out_of_sample_regression(data$testX, MLFSobj)
        pred_acc[k] = cor(pred, data$testy)**2
        pred_acc_train[k] = MLFSobj$pred_acc_train      
      } else{
        MLFSobj = MLFS(data$trainy, data$trainX, type, R, max_iter=20, rotate=TRUE, verbose=FALSE)
        pred = pred_out_of_sample(data$testX, MLFSobj)
        pred_acc[k] = mean(pred == data$testy)
        pred_acc_train[k] = MLFSobj$pred_acc_train
      }
    }
    data.frame(n_irrelevant_features, pred_acc, pred_acc_train)
  }
  names(df)[2:3] = c("test", "train")
  return(df)
}

# Experiment with the dimensionality of latent space
experiment2 = function(R_values, H, d, type, n_rep = 3, n_levels = 3, continuous = FALSE){
  K = length(R_values)
  pred_acc_train = rep(NA, K)
  pred_acc = rep(NA, K)
  df = foreach(i = 1:n_rep, .combine="rbind") %dopar% {
    # generate latent variables
    latent_data = generate_latent_subspace(H, N = 200, d = d)
    y = generate_y(latent_data$V, C = 2, continuous = continuous)
    g_list = lapply(latent_data$U_list, function(x)c(-Inf, quantile(x, 1:(n_levels-1) / n_levels), Inf))
    
    for(k in 1:K){
      X0 = generate_X(latent_data$U_list, type, g_list)
      data = split_into_train_and_test(X0, y, prop=0.5)
      if(continuous){
        MLFSobj = MLFS_regression(data$trainy, data$trainX, type, R, max_iter=20, rotate=TRUE, verbose=FALSE)
        pred = pred_out_of_sample_regression(data$testX, MLFSobj)
        pred_acc[k] = cor(pred, data$testy)**2
        pred_acc_train[k] = MLFSobj$pred_acc_train      
      } else{
        MLFSobj = MLFS(data$trainy, data$trainX, type, R, max_iter=20, rotate=TRUE, verbose=FALSE)
        pred = pred_out_of_sample(data$testX, MLFSobj)
        pred_acc[k] = mean(pred == data$testy)
        pred_acc_train[k] = MLFSobj$pred_acc_train
      }
    }
    data.frame(R_values, pred_acc, pred_acc_train)
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


generate_latent_subspace_with_correlated_features = function(n_added, H, N = 100, d = c(5, 5, 5, 5), gamma=10){
  R = nrow(H)
  V = matrix(rnorm(N*R), N, R)
  W_list = lapply(1:length(d), function(j)rsparsematrix_mod(n_added, R, d[j], H[, j]))
  U_list = lapply(W_list, function(W) V %*% W + rnorm(N*ncol(W), 0, 1/sqrt(gamma)))
  return(list(U_list = U_list, W_list = W_list, V = V, beta = beta))
}

experiment_correlated_features = function(n_added, H, d, type, n_rep = 5, n_levels = 3, continuous = FALSE){
  K = length(n_added)
  pred_acc_train = rep(NA, K)
  pred_acc = rep(NA, K)
  df = foreach(i = 1:n_rep, .combine="rbind") %dopar% {
    for(k in 1:K){
      # generate latent variables
      latent_data = generate_latent_subspace_with_correlated_features(n_added[k], H, N = 200, d = d)
      y = generate_y(latent_data$V, C = 2, continuous = continuous)
      g_list = lapply(latent_data$U_list, function(x)c(-Inf, quantile(x, 1:(n_levels-1) / n_levels), Inf))
      X0 = generate_X(latent_data$U_list, type, g_list)
      data = split_into_train_and_test(X0, y, prop=0.5)
      if(continuous){
        MLFSobj = MLFS_regression(data$trainy, data$trainX, type, R, max_iter=20, rotate=TRUE, verbose=FALSE)
        pred = pred_out_of_sample_regression(data$testX, MLFSobj)
        pred_acc[k] = cor(pred, data$testy)**2
        pred_acc_train[k] = MLFSobj$pred_acc_train      
      } else{
        MLFSobj = MLFS(data$trainy, data$trainX, type, R, max_iter=20, rotate=TRUE, verbose=FALSE)
        pred = pred_out_of_sample(data$testX, MLFSobj)
        pred_acc[k] = mean(pred == data$testy)
        pred_acc_train[k] = MLFSobj$pred_acc_train
      }
    }
    data.frame(n_added, pred_acc, pred_acc_train)
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

generate_latent_subspace_t = function(H, N = 100, d = c(5, 5, 5, 5), gamma=1000, df=1){
  R = nrow(H)
  
  V = matrix(rnorm(N*R), N, R)
  W_list = lapply(1:length(d), function(j)rsparsematrix(R, d[j], H[, j]))
  
  U_list = lapply(W_list, function(W) V %*% W + 1/sqrt(gamma)*rt(N*ncol(W), df))
  
  return(list(U_list = U_list, W_list = W_list, V = V, beta = beta))
}

# Experiment with gaussianity
experiment_gaussianity = function(H, d, type, R, n_rep = 5, n_levels = 3, rotate = TRUE, use_t_distr = TRUE, transformation = FALSE, continuous = FALSE, df = 1){
  pred_acc_train = rep(NA, n_rep)
  pred_acc = rep(NA, n_rep)
  for(i in 1:n_rep){
    # generate latent variables
    if(use_t_distr){
      latent_data = generate_latent_subspace_t(H, N = 200, d = d, gamma = 1, df = df)
    } else{
      latent_data = generate_latent_subspace(H, N = 200, d = d, gamma = 1)
    }
    y = generate_y(latent_data$V, C = 2, continuous = continuous)
    g_list = lapply(latent_data$U_list, function(x)c(-Inf, quantile(x, 1:(n_levels-1) / n_levels), Inf))
    X0 = generate_X(latent_data$U_list, type, g_list)
    if(transformation){
      for(j in 1:length(X0)){
        if(type[j] == "gaussian") X0[[j]] = apply(X0[[j]], 2, inversenormal)
      }
    }
    data = split_into_train_and_test(X0, y, prop=0.5)
    if(continuous){
      MLFSobj = MLFS_regression(data$trainy, data$trainX, type, R, max_iter=20, rotate=TRUE, verbose=FALSE)
      pred = pred_out_of_sample_regression(data$testX, MLFSobj)
      pred_acc[i] = cor(pred, data$testy)**2
      pred_acc_train[i] = MLFSobj$pred_acc_train      
    } else{
      MLFSobj = MLFS(data$trainy, data$trainX, type, R, max_iter=20, rotate=TRUE, verbose=FALSE)
      pred = pred_out_of_sample(data$testX, MLFSobj)
      pred_acc[i] = mean(pred == data$testy)
      pred_acc_train[i] = MLFSobj$pred_acc_train
    }
  }
  df = data.frame(pred_acc, pred_acc_train)
  names(df)[1:2] = c("test", "train")
  return(df)
}

# Experiment with label switching
generate_label_switching = function(U_list, prop_individuals = 0.05, n_views = 1){
  if(prop_individuals == 0) return(U_list)
  out = U_list
  N = nrow(U_list[[1]])
  M = length(U_list)
  selected_individuals = sample(1:N, floor(prop_individuals*N))
  reordering_selected_individuals = c(selected_individuals[-1], selected_individuals[1])
  for(i in 1:length(selected_individuals)){
    for(j in sample(1:M, n_views)){
      switch_indexes = c(selected_individuals[i], reordering_selected_individuals[i])
      out[[j]][rev(switch_indexes), ] = U_list[[j]][switch_indexes, ]
    }
  }
  return(out)
}

experiment_label_switching = function(prop_individuals = 0.1, H, d, type, R, n_views = 1, t_distribution = FALSE, continuous = FALSE){
  K = length(prop_individuals)
  df = data.frame()
  latent_data = generate_latent_subspace(H, N = 200, d = d, gamma = 10)
  if(t_distribution) latent_data = generate_latent_subspace_t(H, N = 200, d = d, gamma = 10)
  y = generate_y(latent_data$V, C = 2, continuous = continuous)
  X0 = generate_X(latent_data$U_list, type)
  data = split_into_train_and_test(X0, y, prop=0.5)
  df = foreach(j = 1:length(n_views), .combine="rbind") %dopar% {
    pred_acc_train = rep(NA, K)
    pred_acc = rep(NA, K)
    for(i in 1:K){
      X_train_mod = generate_label_switching(data$trainX, prop_individuals[i], n_views[j])
      if(continuous){
        MLFSobj = MLFS_regression(data$trainy, X_train_mod, type, R, max_iter=20, rotate=TRUE, verbose=FALSE)
        pred = pred_out_of_sample_regression(data$testX, MLFSobj)
        pred_acc[i] = cor(pred, data$testy)**2
        pred_acc_train[i] = MLFSobj$pred_acc_train      
      } else{
        MLFSobj = MLFS(data$trainy, X_train_mod, type, R, max_iter=20, rotate=TRUE, verbose=FALSE)
        pred = pred_out_of_sample(data$testX, MLFSobj)
        pred_acc[i] = mean(pred == data$testy)
        pred_acc_train[i] = MLFSobj$pred_acc_train
      }
    }
    data.frame(test = pred_acc, train = pred_acc_train, prop_switched = prop_individuals, n_views = n_views[j])
  }
  return(df)
}

generate_missing_data = function(X, prop, n_views){
  if(prop == 0) return(X)
  N = nrow(X[[1]])
  M = length(X)
  for(j in sample(1:M, n_views)){
    selected_individuals = sample(1:N, floor(prop*N))
    X[[j]][selected_individuals, ] = NA
  }
  return(X)
}

experiment_missing_data = function(prop_individuals = 0.1, H, d, type, R, n_views = 1, n_rep = 5){
  K = length(prop_individuals)
  df = foreach(k = 1:n_rep, .combine="rbind") %dopar% {
    latent_data = generate_latent_subspace(H, N = 200, d = d, gamma = 10)
    y = generate_y(latent_data$V, C = 2, continuous = FALSE)
    X0 = generate_X(latent_data$U_list, type)
    data = split_into_train_and_test(X0, y, prop=0.5)
    foreach(j = 1:length(n_views), .combine = "rbind") %dopar% {
      pred_acc_train = rep(NA, K)
      pred_acc = rep(NA, K)
      for(i in 1:K){
        Xmis = generate_missing_data(data$trainX, prop_individuals[i], n_views[j])
        MLFSobj = MLFS_mcmc(data$trainy, Xmis, data$testy, data$testX, type, R, max_iter=200, verbose=FALSE)
        pred_acc_train[i] = MLFSobj$pred_acc_train
        pred_acc[i] = MLFSobj$pred_acc_test
      }
      data.frame(test = pred_acc, train = pred_acc_train, prop_missing = prop_individuals, n_views = n_views[j])
    }
  }
  return(df)
}
