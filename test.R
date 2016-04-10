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
data = generate_data(H, C = 2, N = 100, d = d)

# generate gaussian data
type = rep("gaussian", length(d))
X_list = generate_X(data$U_list, type)

MLFS(data$y, X_list, type, R, max_iter=10, rotate=TRUE)

# generate ordinal data
type = rep("ordinal", length(d))
n_levels = 3
g_list = lapply(data$U_list, function(x)c(-Inf, quantile(x, 1:(n_levels-1) / n_levels), Inf))
X_list = generate_X(data$U_list, type, g_list)

MLFS(data$y, X_list, type, R, max_iter=20, rotate=TRUE)
