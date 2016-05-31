library(reshape2)
library(ggplot2)
library(gridExtra)
library(doParallel)
library(foreach)
registerDoParallel(cores = 10)

source("MLFS_all.R")
source("helpers_label_switching.R")

data = split_into_train_and_test(X_list, y, prop=0.7)

n_rep = 10

plot_list = foreach(i = 1:n_rep) %dopar% {
  N = length(data$trainy)
  trainXmod = data$trainX
  two_indexes = sample(1:N, 2)
  trainXmod[[2]][two_indexes, ] = trainXmod[[2]][rev(two_indexes), ]
  MLFSobj = MLFS_mcmc(data$trainy, trainXmod, data$testy, data$testX, rep("gaussian", length(data$trainX)), 
                      R = 20, max_iter=10000, verbose=TRUE, label_switching = TRUE)
  mat = MLFSobj$label_state_mat[[2]]
  create_plot_label_switching(mat, two_indexes)
}
args_list = c(plot_list, list("ncol" = 1))
g = do.call("grid.arrange", args_list)
ggsave("fig_label_switching.pdf", g, height=28, width=6)
