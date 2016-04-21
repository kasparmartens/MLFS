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
experiment0 = function(H, d, type, n_rep = 3, n_levels = 3){
  pred_acc_train = rep(NA, n_rep)
  pred_acc = rep(NA, n_rep)
  for(i in 1:n_rep){
    # generate latent variables
    latent_data = generate_latent_subspace(H, N = 200, d = d)
    y = generate_y(latent_data$V, C = 2)
    g_list = lapply(latent_data$U_list, function(x)c(-Inf, quantile(x, 1:(n_levels-1) / n_levels), Inf))

    X0 = generate_X(latent_data$U_list, type, g_list)
    data = split_into_train_and_test(X0, y, prop=0.5)
    MLFSobj = MLFS(data$trainy, data$trainX, type, R, max_iter=20, rotate=TRUE, verbose=FALSE)
    pred = pred_out_of_sample(data$testX, MLFSobj)
    pred_acc[i] = mean(pred == data$testy)
    pred_acc_train[i] = MLFSobj$pred_acc_train
  }
  df = data.frame(pred_acc, pred_acc_train)
  names(df) = c("test", "train")
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
