source("helpers.R")
source("par_updates.R")
source("rotation.R")
source("ordinal_cutpoints.R")
source("MLFS.R")
source("generate_data.R")

H = rbind(c(1, 1, 1, 1), 
          c(1, 1, 1, 0), 
          c(1, 0, 1, 1), 
          c(1, 1, 0, 1), 
          c(1, 0, 0, 1), 
          c(0, 1, 1, 0))

set.seed(0)
R = nrow(H)
d = c(5, 5, 5, 5)
latent_data = generate_latent_subspace(H, N = 200, d = d)
y = generate_y(latent_data$V, C = 2)

# generate gaussian data
type = rep("gaussian", length(d))
X0 = generate_X(latent_data$U_list, type)
X1 = add_irrelevant_features(X0, type, 1)
data = split_into_train_and_test(X1, y, prop=0.5)
MLFSobj = MLFS(data$trainy, data$trainX, type, R, max_iter=10, rotate=TRUE)
pred = pred_out_of_sample(data$testX, MLFSobj)
mean(pred == data$testy)

# generate ordinal data
type = rep("ordinal", length(d))
n_levels = 3
g_list = lapply(latent_data$U_list, function(x)c(-Inf, quantile(x, 1:(n_levels-1) / n_levels), Inf))
X0 = generate_X(latent_data$U_list, type, g_list)
X1 = add_irrelevant_features(X0, type, 1)
data = split_into_train_and_test(X1, y, prop=0.5)
MLFSobj = MLFS(data$trainy, data$trainX, type, R, max_iter=10, rotate=TRUE)
pred = pred_out_of_sample(data$testX, MLFSobj)
mean(pred == data$testy)

