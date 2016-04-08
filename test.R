source("helpers.R")
source("par_updates.R")
source("MLFS.R")
source("generate_data.R")

# generate data
H = rbind(c(1, 1, 1, 1), 
          c(1, 1, 1, 0), 
          c(1, 0, 1, 1), 
          c(1, 1, 0, 1), 
          c(1, 0, 0, 1), 
          c(0, 1, 1, 0))
R = nrow(H)
d = c(5, 5, 5, 5)
type = rep("gaussian", length(d))
data = generate_gaussian_views(H, C = 2, N = 100, d = d)
X_list = generate_X(data$U_list, type)

MLFS(data$y, X_list, data$type, R, max_iter=10, rotate=TRUE)
