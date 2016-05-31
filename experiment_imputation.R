library(reshape2)
library(ggplot2)
library(doParallel)
library(foreach)
registerDoParallel(cores = 10)

source("helpers_imputation.R")
source("MLFS_all.R")
source("MLFS_imputation.R")

n_rep = 10
set.seed(0)
data_splitted = lapply(1:n_rep, function(x)split_into_train_and_test(X_list, y, prop=0.7))

missingness_train = list(c(0.25, 0, 0, 0), 
                      c(0.50, 0, 0, 0), 
                      c(0.75, 0, 0, 0), 
                      c(0, 0.25, 0, 0), 
                      c(0, 0.50, 0, 0), 
                      c(0, 0.75, 0, 0), 
                      c(0, 0, 0.25, 0), 
                      c(0, 0, 0.50, 0), 
                      c(0, 0, 0.75, 0))

missingness_test = list(c(0, 1, 0, 0), 
                     c(0, 1, 0, 0),
                     c(0, 1, 0, 0), 
                     c(1, 0, 1, 0),
                     c(1, 0, 1, 0),
                     c(1, 0, 1, 0),
                     c(0, 1, 0, 0),
                     c(0, 1, 0, 0),
                     c(0, 1, 0, 0))

df = foreach(j = 1:n_rep, .combine="rbind") %dopar% {
  data = data_splitted[[j]]
  imputation_experiment(missingness_train, missingness_test, 
                        data$trainy, data$trainX, data$testy, data$testX, 
                        R=20, max_iter = 3000, burnin=500)
}


names(df)[1:5] = c("multi-view impute", "multi-view available views", "multi-view impute train", "random forest", "random forest impute")
df.m = melt(df, id.vars = c("missing_pattern"))

p = ggplot(df.m) + 
  geom_boxplot(aes(factor(1), value, col=variable)) + 
  theme_bw() + scale_color_brewer("Method", palette = "Set1") + 
  facet_wrap(~ missing_pattern) + 
  ylim(0, 1) + scale_x_discrete(breaks=NULL) + 
  xlab("") + ylab("Prediction accuracy") 

ggsave("fig_imputation.pdf", p, width=9, height=6)
