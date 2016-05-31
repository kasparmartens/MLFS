library(randomForest)
library(reshape2)

baseline_classifiers = function(trainX, trainy, testX, testy){
  # each of the views separately
  df1 = foreach(i = 1:length(trainX), .combine="rbind") %do% {
    acc1 = logreg_predacc(trainX[[i]], trainy, testX[[i]], testy)
    acc2 = rf_predacc(trainX[[i]], trainy, testX[[i]], testy)
    data.frame("logreg" = acc1, "randomforest" = acc2, type = paste0("view ", i))
  }
  # combine separate classifiers
  acc1 = logreg_avg_predacc(trainX, trainy, testX, testy)
  acc2 = rf_avg_predacc(trainX, trainy, testX, testy)
  df2 = data.frame("logreg" = acc1, "randomforest" = acc2, type = "averaged prediction")
  # joint classification
  acc1 = logreg_predacc(do.call("cbind", trainX), trainy, do.call("cbind", testX), testy)
  acc2 = rf_predacc(do.call("cbind", trainX), trainy, do.call("cbind", testX), testy)
  df3 = data.frame("logreg" = acc1, "randomforest" = acc2, type = "concatenated views")
  df = rbind(df1, df2, df3)
  df.m = melt(df, id.vars = "type")
  df.m$variable = factor(df.m$variable, levels = c("logreg", "randomforest"), labels = c("logistic regression", "random forest"))
  return(df.m)
}

logreg_predacc = function(trainX, trainy, testX, testy){
  m = glm(trainy ~ ., data=data.frame(trainX), family=binomial())
  pred = 1*(predict(m, newdata = data.frame(testX)) > 0)
  return(mean(pred == testy))
}

rf_predacc = function(trainX, trainy, testX, testy){
  m = randomForest(trainX, factor(trainy), ntree = 1000)
  pred = predict(m, newdata = testX)
  return(mean(pred == testy))
}

logreg_avg_predacc = function(trainX_list, trainy, testX_list, testy){
  preds = foreach(i = 1:length(trainX_list), .combine="cbind") %do% {
    m = glm(trainy ~ ., data=data.frame(trainX_list[[i]]), family=binomial())
    predict(m, newdata = data.frame(testX_list[[i]]))
  }
  pred = 1*(apply(preds, 1, mean) > 0)
  return(mean(pred == testy))
}

rf_avg_predacc = function(trainX_list, trainy, testX_list, testy){
  preds = foreach(i = 1:length(trainX_list), .combine="cbind") %do% {
    m = randomForest(trainX_list[[i]], factor(trainy), ntree = 1000)
    predict(m, newdata = testX_list[[i]], type="prob")[,2]
  }
  logitprobs = logit(preds)
  pred = 1*(apply(logitprobs, 1, mean) > 0)
  return(mean(pred == testy))
}
