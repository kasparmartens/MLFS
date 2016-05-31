library(reshape2)
library(ggplot2)
library(doParallel)
library(foreach)
registerDoParallel(cores = 10)

source("helpers_baseline_classifiers.R")
source("MLFS_all.R")

n_rep = 8
set.seed(0)
# Make sure the data_splitted is the same as beforehand. The rest of the code should (hopefully) work fine
data_splitted = lapply(1:n_rep, function(x)split_into_train_and_test(X0, y, prop=0.7))


# get results for baseline classifiers
df.baseline = foreach(j = 1:n_rep, .combine="rbind") %dopar% {
  data = data_splitted[[j]]
  baseline_classifiers(data$trainX, data$trainy-1, data$testX, data$testy-1)
}

ggplot(df.baseline, aes(factor(1), value, col=variable)) + 
  geom_boxplot() + facet_wrap(~ type, nrow=1) + 
  theme_bw() + scale_x_discrete(breaks=NULL) +
  scale_color_brewer("Method", palette="Set1") + 
  xlab("") + ylim(0, 1) + ylab("Prediction accuracy")

# now results for multiview
df.multiview = foreach(j = 1:n_rep, .combine="rbind") %dopar% {
  data = data_splitted[[j]]
  MLFSobj = MLFS_mcmc(data$trainy, data$trainX, data$testy, data$testX, rep("gaussian", length(data$trainX)), 
                      R = 6, max_iter=1000, verbose=TRUE)
  data.frame(type = "multi-view", variable = "multi-view", value = MLFSobj$pred_acc_test)
}

df = rbind(df.baseline, df.multiview)


p = ggplot(df, aes(factor(1), value, col=variable)) + 
  geom_boxplot() + facet_wrap(~ type, nrow=1) + 
  theme_bw() + scale_x_discrete(breaks=NULL) +
  scale_color_brewer("Method", palette="Set1") + 
  xlab("") + ylim(0, 1) + ylab("Prediction accuracy")

ggsave("fig_baseline_comparison.pdf", p, width=12, height=4)


df.multiview = foreach(j = 1:n_rep, .combine="rbind") %dopar% {
  foreach(R = c(3, 4, 6, 10), .combine="rbind") %do%{
    data = data_splitted[[j]]
    MLFSobj = MLFS_mcmc(data$trainy, data$trainX, data$testy, data$testX, rep("gaussian", length(data$trainX)), 
                        R = R, max_iter=1000, verbose=TRUE)
    data.frame(R = R, value = MLFSobj$pred_acc_test)
  }
}

ggplot(df, aes(factor(R), value, col="x")) + 
  geom_boxplot() + 
  theme_bw() + scale_color_brewer("Method", palette="Set1") + theme(legend.position="none") +
  xlab("Dimensionality of latent space") + ylim(0, 1) + ylab("Prediction accuracy")
ggsave("fig_dimensionality_latent.pdf", width=6, heigh=3)
