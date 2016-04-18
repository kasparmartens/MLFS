library(ggplot2)
library(dplyr)
library(reshape2)

source("helpers.R")
source("par_updates.R")
source("rotation.R")
source("ordinal_cutpoints.R")
source("MLFS.R")
source("generate_data.R")

H = rbind(c(1, 1), 
          c(1, 1), 
          c(1, 0), 
          c(0, 1), 
          c(0, 1), 
          c(0, 1))

set.seed(0)
R = nrow(H)
d = c(5, 5)
type = c("gaussian", "ordinal")

# Experiment with the number of irrelevant features
experiment = function(n_irrelevant_features, n_rep, U, type, n_levels = 3){
  K = length(n_irrelevant_features)
  pred_acc_train = rep(NA, K)
  pred_acc = rep(NA, K)
  g_list = lapply(latent_data$U_list, function(x)c(-Inf, quantile(x, 1:(n_levels-1) / n_levels), Inf))
  df = data.frame()
  for(i in 1:n_rep){
    for(k in 1:K){
      X0 = generate_X(U, type, g_list)
      X1 = add_irrelevant_features(X0, type, n_irrelevant_features[k])
      data = split_into_train_and_test(X1, y, prop=0.5)
      MLFSobj = MLFS(data$trainy, data$trainX, type, R, max_iter=10, rotate=TRUE)
      pred = pred_out_of_sample(data$testX, MLFSobj)
      pred_acc[k] = mean(pred == data$testy)
      pred_acc_train[k] = MLFSobj$pred_acc_train
    }
    temp = data.frame(n_irrelevant_features, pred_acc, pred_acc_train)
    df = rbind(df, temp)
  }
  return(df)
}

n_irrelevant_features = c(0, 10, 100, 1000, 10000)
res = experiment(n_irrelevant_features, n_rep = 5, U = latent_data$U_list, type = type)
df.m = res %>%
  melt(id.vars = c("n_irrelevant_features")) 
df = df.m %>%
  group_by(n_irrelevant_features, variable) %>%
  summarise(mean = mean(value), min = min(value), max = max(value))

ggplot(df.m, aes(factor(n_irrelevant_features), value, col=variable)) + 
  geom_boxplot() + 
  theme_bw() + ylim(0, 1)
